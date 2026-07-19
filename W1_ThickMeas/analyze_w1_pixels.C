#include <TFile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TStyle.h>

#include "strip_map.h"

#include <array>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

struct Hit {
  Int_t channel = -1;
  Int_t board = -1;
  ULong64_t timestamp = 0;
  UShort_t energy = 0;
  UShort_t energyShort = 0;
  UInt_t flags = 0;
};

struct PixelSpectrum {
  std::unique_ptr<TH1D> hist;
  double peak = 0.0;
  double sigma = 0.0;
  double fwhm = 0.0;
  bool fitted = false;
};

static constexpr double kTimeDiffHistogramMaxPs = 400000.0;
static constexpr double kFallbackCoincidenceGatePs = 5000.0;
static constexpr double kDefaultAdcThreshold = 100.0;
static constexpr double kDefaultEnergyMin = 0.0;
static constexpr double kDefaultEnergyMax = 4096.0;
static constexpr double kDefaultEnergyBinWidth = 1.0;

static void SortHitsByTimestamp(std::vector<Hit> &hits) {
  std::sort(hits.begin(), hits.end(), [](const Hit &lhs, const Hit &rhs) {
    if (lhs.timestamp != rhs.timestamp) {
      return lhs.timestamp < rhs.timestamp;
    }
    if (lhs.board != rhs.board) {
      return lhs.board < rhs.board;
    }
    return lhs.channel < rhs.channel;
  });
}

static void PruneOldHits(std::vector<Hit> &hits, ULong64_t currentTimestamp, double gatePs) {
  const ULong64_t minTimestamp = currentTimestamp > static_cast<ULong64_t>(gatePs) ? currentTimestamp - static_cast<ULong64_t>(gatePs) : 0;
  while (!hits.empty() && hits.front().timestamp < minTimestamp) {
    hits.erase(hits.begin());
  }
}

static int FindClosestHitIndex(const std::vector<Hit> &hits, ULong64_t timestamp, double gatePs) {
  int bestIndex = -1;
  ULong64_t bestDt = 0;
  for (int i = static_cast<int>(hits.size()) - 1; i >= 0; --i) {
    const ULong64_t dt = hits[i].timestamp > timestamp ? hits[i].timestamp - timestamp : timestamp - hits[i].timestamp;
    if (dt > static_cast<ULong64_t>(gatePs)) {
      continue;
    }
    if (bestIndex < 0 || dt < bestDt) {
      bestIndex = i;
      bestDt = dt;
    }
  }
  return bestIndex;
}

static bool FindNearestTimestampDifference(const std::vector<Hit> &hits, ULong64_t timestamp, ULong64_t &bestDt) {
  if (hits.empty()) {
    return false;
  }

  const auto lower = std::lower_bound(hits.begin(), hits.end(), timestamp, [](const Hit &hit, ULong64_t value) {
    return hit.timestamp < value;
  });

  bool found = false;
  bestDt = 0;

  auto consider = [&](const Hit &hit) {
    const ULong64_t dt = hit.timestamp > timestamp ? hit.timestamp - timestamp : timestamp - hit.timestamp;
    if (!found || dt < bestDt) {
      bestDt = dt;
      found = true;
    }
  };

  if (lower != hits.end()) {
    consider(*lower);
  }
  if (lower != hits.begin()) {
    consider(*(lower - 1));
  }

  return found;
}

static std::filesystem::path ResolveRunDirectory(const char *inputFile) {
  const std::filesystem::path inputPath = std::filesystem::absolute(std::filesystem::path(inputFile)).lexically_normal();
  const std::filesystem::path rawDir = inputPath.parent_path();
  if (rawDir.filename() == "RAW" && !rawDir.parent_path().empty()) {
    return rawDir.parent_path();
  }
  return rawDir;
}

static std::string ExtractTrailingDigits(const std::string &text) {
  const auto pos = text.find_last_not_of("0123456789");
  if (pos == std::string::npos) {
    return text;
  }
  if (pos + 1 >= text.size()) {
    return std::string();
  }
  return text.substr(pos + 1);
}

static std::string MakeOutputName(const char *inputFile, const char *suffix, const char *fallbackName) {
  const std::filesystem::path inputPath(inputFile);
  const std::filesystem::path runDir = ResolveRunDirectory(inputFile);
  const std::string runFolder = runDir.filename().string();
  const std::string runNumber = ExtractTrailingDigits(runFolder);
  const std::string stem = inputPath.stem().string();

  if (stem.empty()) {
    return (runDir / fallbackName).string();
  }

  if (!runNumber.empty()) {
    return (runDir / (stem + "_" + runNumber + suffix)).string();
  }

  return (runDir / (stem + suffix)).string();
}

static std::string MakePeaksName(const char *summaryPdf) {
  const std::filesystem::path pdfPath(summaryPdf);
  std::string stem = pdfPath.stem().string();
  constexpr const char *kSummarySuffix = "_summary";
  constexpr std::size_t kSummarySuffixLen = 8;
  if (stem.size() >= kSummarySuffixLen &&
      stem.compare(stem.size() - kSummarySuffixLen, kSummarySuffixLen, kSummarySuffix) == 0) {
    stem.replace(stem.size() - kSummarySuffixLen, kSummarySuffixLen, "_peaks");
  } else {
    stem += "_peaks";
  }
  return (pdfPath.parent_path() / (stem + ".txt")).string();
}

static bool WritePixelPeaksFile(const char *peaksFile,
                                const std::array<std::array<PixelSpectrum, kNY>, kNX> &pixelSummary) {
  std::ofstream out(peaksFile);
  if (!out) {
    std::cerr << "Cannot write peaks file: " << peaksFile << std::endl;
    return false;
  }

  out << "# x y peak_adc fwhm resolution_frac entries\n";
  out << "# resolution_frac = fwhm / peak (unitless; multiply by 100 for %)\n";
  out << std::setprecision(8);
  for (int ix = 0; ix < kNX; ++ix) {
    for (int iy = 0; iy < kNY; ++iy) {
      if (!pixelSummary[ix][iy].fitted || !(pixelSummary[ix][iy].peak > 0.0)) {
        continue;
      }
      const double peak = pixelSummary[ix][iy].peak;
      const double fwhm = pixelSummary[ix][iy].fwhm;
      const double entries = pixelSummary[ix][iy].hist ? pixelSummary[ix][iy].hist->GetEntries() : 0.0;
      out << ix << ' ' << iy << ' ' << peak << ' ' << fwhm << ' ' << (fwhm / peak) << ' ' << entries << '\n';
    }
  }
  return true;
}

static double DetermineCoincidenceGatePs(TH1D &hTimeDiff) {
  if (hTimeDiff.GetEntries() < 20) {
    return kFallbackCoincidenceGatePs;
  }

  TF1 fit("time_fit", "gaus", 0.0, kTimeDiffHistogramMaxPs);
  fit.SetParameters(hTimeDiff.GetMaximum(), hTimeDiff.GetMean(), std::max(1.0, hTimeDiff.GetRMS()));
  hTimeDiff.Fit(&fit, "QNR");

  const double meanPs = fit.GetParameter(1);
  const double sigmaPs = std::abs(fit.GetParameter(2));
  if (!(sigmaPs > 0.0) || !std::isfinite(sigmaPs)) {
    return kFallbackCoincidenceGatePs;
  }

  const double gatePs = meanPs + 5.0 * sigmaPs;
  return std::min(kTimeDiffHistogramMaxPs, std::max(kFallbackCoincidenceGatePs, gatePs));
}

static std::unique_ptr<TFile> OpenInputFile(const char *inputFile) {
  auto tryOpen = [](const std::filesystem::path &path) {
    return std::unique_ptr<TFile>(TFile::Open(path.string().c_str(), "READ"));
  };

  std::filesystem::path requested(inputFile);
  auto file = tryOpen(requested);
  if (file && !file->IsZombie()) {
    return file;
  }

  const std::filesystem::path fileName = requested.filename();
  const std::array<std::filesystem::path, 3> candidates = {
      std::filesystem::path(fileName),
      std::filesystem::path("RAW") / fileName,
      std::filesystem::path(".") / fileName};

  for (const auto &candidate : candidates) {
    file = tryOpen(candidate);
    if (file && !file->IsZombie()) {
      return file;
    }
  }

  return nullptr;
}

static void WriteSummaryPdf(const char *pdfName,
                            const std::array<std::unique_ptr<TH1D>, kNX> &xSpectra,
                            const std::array<std::unique_ptr<TH1D>, kNY> &ySpectra,
                            TH2D &hHitMap,
                            TH1D &hTimeDiff,
                            TH2D &hResolution,
                            TH2D &hMean) {
  TCanvas c("c_summary", "W1 summary", 1400, 1000);

  c.Print(Form("%s[", pdfName));

  gStyle->SetOptStat(1);

  c.Clear();
  c.Divide(4, 4);
  for (int i = 0; i < kNX; ++i) {
    c.cd(i + 1);
    xSpectra[i]->Draw();
  }
  c.Print(pdfName);

  c.Clear();
  c.Divide(4, 4);
  for (int i = 0; i < kNY; ++i) {
    c.cd(i + 1);
    ySpectra[i]->Draw();
  }
  c.Print(pdfName);

  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("g");

  c.Clear();
  hHitMap.Draw("COLZ TEXT");
  c.Print(pdfName);

  c.Clear();
  hTimeDiff.Draw();
  c.Print(pdfName);

  c.Clear();
  const std::unique_ptr<TH2D> hResolutionPercent(static_cast<TH2D *>(hResolution.Clone("hResolutionPercent")));
  for (int x = 1; x <= hResolutionPercent->GetNbinsX(); ++x) {
    for (int y = 1; y <= hResolutionPercent->GetNbinsY(); ++y) {
      hResolutionPercent->SetBinContent(x, y, 100.0 * hResolutionPercent->GetBinContent(x, y));
    }
  }
  gStyle->SetPaintTextFormat("4.2f%%");
  hResolutionPercent->SetTitle("Pixel resolution (FWHM/peak) [%];X;Y");
  hResolutionPercent->Draw("COLZ TEXT");
  c.Print(pdfName);

  c.Clear();
  gStyle->SetPaintTextFormat("g");
  hMean.Draw("COLZ TEXT");
  c.Print(pdfName);

  c.Print(Form("%s]", pdfName));
}

int analyze_w1_pixels(const char *inputFile = "RAW/DataR_thick_meas_calib_3.root",
                      double coincidenceGatePs = 0.0,
                      double adcThreshold = kDefaultAdcThreshold,
                      double energyMin = kDefaultEnergyMin,
                      double energyMax = kDefaultEnergyMax,
                      double energyBinWidth = kDefaultEnergyBinWidth,
                      const char *channelMapFile = "channel_map.txt",
                      const char *summaryPdf = "") {
  TH1::AddDirectory(kFALSE);

  StripMap stripMap;
  if (!LoadStripMap(channelMapFile, stripMap)) {
    return 5;
  }

  if (!(energyMax > energyMin) || !(energyBinWidth > 0.0)) {
    std::cerr << "Invalid energy spectrum axis: min=" << energyMin
              << " max=" << energyMax << " binWidth=" << energyBinWidth << std::endl;
    return 3;
  }
  const int nEnergyBins = std::max(1, static_cast<int>(std::llround((energyMax - energyMin) / energyBinWidth)));
  const double energyAxisMax = energyMin + nEnergyBins * energyBinWidth;

  auto in = OpenInputFile(inputFile);
  if (!in || in->IsZombie()) {
    std::cerr << "Cannot open input file: " << inputFile << std::endl;
    return 1;
  }

  const std::string derivedPdfOutput = MakeOutputName(inputFile, "_w1_pixel_summary.pdf", "w1_pixel_summary.pdf");
  const char *actualSummaryPdf = (summaryPdf && summaryPdf[0] != '\0') ? summaryPdf : derivedPdfOutput.c_str();
  std::cout << "Using ADC threshold: " << adcThreshold << std::endl;
  std::cout << "Energy spectra: [" << energyMin << ", " << energyAxisMax << "), "
            << nEnergyBins << " bins, width=" << energyBinWidth << std::endl;

  auto *tree = dynamic_cast<TTree *>(in->Get("Data_R"));
  if (!tree) {
    std::cerr << "Cannot find TTree Data_R in " << inputFile << std::endl;
    return 2;
  }

  UShort_t channel = 0;
  UShort_t board = 0;
  ULong64_t timestamp = 0;
  UShort_t energy = 0;
  UShort_t energyShort = 0;
  UInt_t flags = 0;

  tree->SetBranchAddress("Channel", &channel);
  tree->SetBranchAddress("Board", &board);
  tree->SetBranchAddress("Timestamp", &timestamp);
  tree->SetBranchAddress("Energy", &energy);
  tree->SetBranchAddress("EnergyShort", &energyShort);
  tree->SetBranchAddress("Flags", &flags);

  std::array<std::unique_ptr<TH1D>, kNX> hXSpectra;
  std::array<std::unique_ptr<TH1D>, kNY> hYSpectra;
  for (int i = 0; i < kNX; ++i) {
    hXSpectra[i] = std::make_unique<TH1D>(Form("hX_%02d", i), Form("X strip %02d", i),
                                          nEnergyBins, energyMin, energyAxisMax);
  }
  for (int i = 0; i < kNY; ++i) {
    hYSpectra[i] = std::make_unique<TH1D>(Form("hY_%02d", i), Form("Y strip %02d", i),
                                          nEnergyBins, energyMin, energyAxisMax);
  }

  std::vector<Hit> hits;
  hits.reserve(tree->GetEntries());

  const Long64_t nEntries = tree->GetEntries();
  for (Long64_t entry = 0; entry < nEntries; ++entry) {
    tree->GetEntry(entry);
    if (static_cast<double>(energy) < adcThreshold) {
      continue;
    }

    MappedStrip mapped;
    if (!stripMap.Lookup(board, channel, mapped)) {
      continue;
    }

    Hit hit;
    hit.channel = channel;
    hit.board = board;
    hit.timestamp = timestamp;
    hit.energy = energy;
    hit.energyShort = energyShort;
    hit.flags = flags;
    hits.push_back(hit);

    if (mapped.axis == StripAxis::kX) {
      hXSpectra[mapped.strip]->Fill(energy);
    } else if (mapped.axis == StripAxis::kY) {
      hYSpectra[mapped.strip]->Fill(energy);
    }
  }

  SortHitsByTimestamp(hits);

  std::vector<Hit> xHits;
  std::vector<Hit> yHits;
  xHits.reserve(hits.size());
  yHits.reserve(hits.size());
  for (const auto &hit : hits) {
    MappedStrip mapped;
    if (!stripMap.Lookup(hit.board, hit.channel, mapped)) {
      continue;
    }
    if (mapped.axis == StripAxis::kX) {
      xHits.push_back(hit);
    } else if (mapped.axis == StripAxis::kY) {
      yHits.push_back(hit);
    }
  }

  std::array<std::array<std::unique_ptr<TH1D>, kNY>, kNX> pixelEnergy;
  std::array<std::array<std::unique_ptr<TH1D>, kNY>, kNX> pixelDiff;
  std::array<std::array<PixelSpectrum, kNY>, kNX> pixelSummary;

  for (int ix = 0; ix < kNX; ++ix) {
    for (int iy = 0; iy < kNY; ++iy) {
      pixelEnergy[ix][iy] = std::make_unique<TH1D>(Form("hPixE_%02d_%02d", ix, iy),
                                                   Form("Pixel X=%02d Y=%02d calibrated energy", ix, iy),
                                                   nEnergyBins, energyMin, energyAxisMax);
      pixelDiff[ix][iy] = std::make_unique<TH1D>(Form("hPixD_%02d_%02d", ix, iy),
                                                 Form("Pixel X=%02d Y=%02d X-Y difference", ix, iy),
                                                 1000, -500, 500);
      pixelSummary[ix][iy].hist = std::make_unique<TH1D>(Form("hPixSum_%02d_%02d", ix, iy),
                                                         Form("Pixel X=%02d Y=%02d X-only calibrated energy", ix, iy),
                                                         nEnergyBins, energyMin, energyAxisMax);
    }
  }

  TH2D hHitMap("hHitMap", "Pixel hit map;X;Y", kNX, -0.5, kNX - 0.5, kNY, -0.5, kNY - 0.5);
  TH1D hTimeDiff("hTimeDiff", "Nearest X-Y timestamp difference;|#Delta t| [ps];Counts", 400, 0, kTimeDiffHistogramMaxPs);

  for (const auto &xHit : xHits) {
    ULong64_t nearestDt = 0;
    if (FindNearestTimestampDifference(yHits, xHit.timestamp, nearestDt)) {
      hTimeDiff.Fill(static_cast<double>(nearestDt));
    }
  }

  const double fittedCoincidenceGatePs = DetermineCoincidenceGatePs(hTimeDiff);
  const double effectiveCoincidenceGatePs = coincidenceGatePs > 0.0 ? coincidenceGatePs : fittedCoincidenceGatePs;
  std::cout << "Coincidence gate from hTimeDiff fit: " << fittedCoincidenceGatePs << " ps" << std::endl;
  std::cout << "Using coincidence gate: " << effectiveCoincidenceGatePs << " ps" << std::endl;

  std::vector<Hit> xBuffer;
  std::vector<Hit> yBuffer;
  xBuffer.reserve(128);
  yBuffer.reserve(128);

  std::vector<std::pair<Hit, Hit>> coincidences;
  coincidences.reserve(hits.size() / 2);

  for (const auto &hit : hits) {
    MappedStrip mapped;
    if (!stripMap.Lookup(hit.board, hit.channel, mapped)) {
      continue;
    }
    const bool isX = (mapped.axis == StripAxis::kX);
    const bool isY = (mapped.axis == StripAxis::kY);
    if (!isX && !isY) {
      continue;
    }

    auto &thisBuffer = isX ? xBuffer : yBuffer;
    auto &otherBuffer = isX ? yBuffer : xBuffer;

    PruneOldHits(xBuffer, hit.timestamp, effectiveCoincidenceGatePs);
    PruneOldHits(yBuffer, hit.timestamp, effectiveCoincidenceGatePs);

    const int matchIndex = FindClosestHitIndex(otherBuffer, hit.timestamp, effectiveCoincidenceGatePs);
    if (matchIndex >= 0) {
      const Hit &other = otherBuffer[matchIndex];
      if (isX) {
        coincidences.emplace_back(hit, other);
      } else {
        coincidences.emplace_back(other, hit);
      }
      otherBuffer.erase(otherBuffer.begin() + matchIndex);
      continue;
    }

    thisBuffer.push_back(hit);
  }

  for (const auto &pair : coincidences) {
    const Hit &xHit = pair.first;
    const Hit &yHit = pair.second;

    MappedStrip xMapped;
    MappedStrip yMapped;
    if (!stripMap.Lookup(xHit.board, xHit.channel, xMapped) || xMapped.axis != StripAxis::kX) {
      continue;
    }
    if (!stripMap.Lookup(yHit.board, yHit.channel, yMapped) || yMapped.axis != StripAxis::kY) {
      continue;
    }

    const int xStrip = xMapped.strip;
    const int yStrip = yMapped.strip;
    const double xE = xHit.energy;
    const double yE = yHit.energy;

    const double avgEnergy = 0.5 * (xE + yE);
    const double diffEnergy = xE - yE;

    hHitMap.Fill(xStrip, yStrip);
    pixelEnergy[xStrip][yStrip]->Fill(avgEnergy);
    pixelDiff[xStrip][yStrip]->Fill(diffEnergy);
    pixelSummary[xStrip][yStrip].hist->Fill(xE);
  }

  if (coincidences.empty()) {
    std::cout << "No X/Y coincidences were found within the time gate." << std::endl;
  }

  for (int ix = 0; ix < kNX; ++ix) {
    for (int iy = 0; iy < kNY; ++iy) {
      auto *hist = pixelSummary[ix][iy].hist.get();
      if (hist->GetEntries() < 1) {
        continue;
      }

      if (hist->GetEntries() < 20) {
        pixelSummary[ix][iy].peak = hist->GetMean();
        const double minSigma = 0.5 * hist->GetXaxis()->GetBinWidth(1);
        pixelSummary[ix][iy].sigma = std::max(hist->GetRMS(), minSigma);
        pixelSummary[ix][iy].fwhm = 2.354820045 * pixelSummary[ix][iy].sigma;
      } else {
        const int maxBin = hist->GetMaximumBin();
        const double peak = hist->GetBinCenter(maxBin);
        const double fitMin = std::max(0.0, peak - 80.0);
        const double fitMax = peak + 80.0;
        TF1 fit("pixel_gaus", "gaus", fitMin, fitMax);
        hist->Fit(&fit, "QNR");

        pixelSummary[ix][iy].peak = fit.GetParameter(1);
        pixelSummary[ix][iy].sigma = std::abs(fit.GetParameter(2));
        pixelSummary[ix][iy].fwhm = 2.354820045 * pixelSummary[ix][iy].sigma;
      }
      pixelSummary[ix][iy].fitted = true;
    }
  }

  TH2D hResolution("hResolution", "Pixel resolution (FWHM/peak);X;Y",
                   kNX, -0.5, kNX - 0.5, kNY, -0.5, kNY - 0.5);
  TH2D hMean("hPixelPeak", "Pixel peak energy;X;Y", kNX, -0.5, kNX - 0.5, kNY, -0.5, kNY - 0.5);
  for (int ix = 0; ix < kNX; ++ix) {
    for (int iy = 0; iy < kNY; ++iy) {
      if (!pixelSummary[ix][iy].fitted) {
        continue;
      }
      const double peak = pixelSummary[ix][iy].peak;
      if (!(peak > 0.0)) {
        continue;
      }
      hResolution.SetBinContent(ix + 1, iy + 1, pixelSummary[ix][iy].fwhm / peak);
      hMean.SetBinContent(ix + 1, iy + 1, peak);
    }
  }

  WriteSummaryPdf(actualSummaryPdf, hXSpectra, hYSpectra, hHitMap, hTimeDiff, hResolution, hMean);

  const std::string peaksFile = MakePeaksName(actualSummaryPdf);
  if (!WritePixelPeaksFile(peaksFile.c_str(), pixelSummary)) {
    return 4;
  }

  in->Close();

  std::cout << "Wrote " << actualSummaryPdf << std::endl;
  std::cout << "Wrote " << peaksFile << std::endl;
  std::cout << "Board/Channel -> X/Y strip mapping is defined in " << channelMapFile
            << "; edit that file to correct cabling order or axis flips." << std::endl;
  std::cout << "Hits below the ADC threshold are rejected before spectra and coincidence matching." << std::endl;
  std::cout << "The timestamp gate is fit from hTimeDiff in ps, then used as a one-to-one coincidence window." << std::endl;
  std::cout << "hTimeDiff stores the nearest Y-side timestamp difference for every X hit, in ps." << std::endl;
  return 0;
}
