#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TTree.h>

#include "strip_map.h"

#include <array>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

struct Hit {
  Int_t channel = -1;
  Int_t board = -1;
  ULong64_t timestamp = 0;
  UShort_t energy = 0;
};

struct PixelCalib {
  double slope = 0.0;
  double offset = 0.0;
  bool valid = false;
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
static constexpr double kDefaultAdcThreshold = 50.0;
static constexpr double kDefaultEnergyMinKeV = 2000.0;
static constexpr double kDefaultEnergyMaxKeV = 7000.0;
static constexpr double kDefaultEnergyBinWidthKeV = 1.0;
static constexpr Long64_t kDefaultMinPixelEntries = 500;

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
  const ULong64_t minTimestamp =
      currentTimestamp > static_cast<ULong64_t>(gatePs) ? currentTimestamp - static_cast<ULong64_t>(gatePs) : 0;
  while (!hits.empty() && hits.front().timestamp < minTimestamp) {
    hits.erase(hits.begin());
  }
}

static int FindClosestHitIndex(const std::vector<Hit> &hits, ULong64_t timestamp, double gatePs) {
  int bestIndex = -1;
  ULong64_t bestDt = 0;
  for (int i = static_cast<int>(hits.size()) - 1; i >= 0; --i) {
    const ULong64_t dt =
        hits[i].timestamp > timestamp ? hits[i].timestamp - timestamp : timestamp - hits[i].timestamp;
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

static std::vector<std::string> SplitCommaSeparated(const char *text) {
  std::vector<std::string> parts;
  if (!text) {
    return parts;
  }
  std::stringstream ss(text);
  std::string item;
  while (std::getline(ss, item, ',')) {
    const auto start = item.find_first_not_of(" \t");
    const auto end = item.find_last_not_of(" \t");
    if (start == std::string::npos) {
      continue;
    }
    parts.push_back(item.substr(start, end - start + 1));
  }
  return parts;
}

static bool LoadPixelCalib(const char *path, std::array<std::array<PixelCalib, kNY>, kNX> &calib) {
  std::ifstream in(path);
  if (!in) {
    std::cerr << "Cannot open calibration file: " << path << std::endl;
    return false;
  }

  std::string line;
  int lineNumber = 0;
  int loaded = 0;
  while (std::getline(in, line)) {
    ++lineNumber;
    if (line.empty() || line[0] == '#') {
      continue;
    }
    std::istringstream iss(line);
    int x = -1;
    int y = -1;
    double slope = 0.0;
    double offset = 0.0;
    if (!(iss >> x >> y >> slope >> offset)) {
      std::cerr << path << ":" << lineNumber << ": failed to parse calibration line" << std::endl;
      return false;
    }
    if (x < 0 || x >= kNX || y < 0 || y >= kNY) {
      std::cerr << path << ":" << lineNumber << ": pixel index out of range" << std::endl;
      return false;
    }
    calib[x][y].slope = slope;
    calib[x][y].offset = offset;
    calib[x][y].valid = true;
    ++loaded;
  }

  std::cout << "Loaded " << loaded << " pixel calibration entries from " << path << std::endl;
  return loaded > 0;
}

static bool AppendHitsFromFile(const char *inputFile,
                               double adcThreshold,
                               const StripMap &stripMap,
                               std::vector<Hit> &hits) {
  std::unique_ptr<TFile> in(TFile::Open(inputFile, "READ"));
  if (!in || in->IsZombie()) {
    std::cerr << "Cannot open input file: " << inputFile << std::endl;
    return false;
  }

  auto *tree = dynamic_cast<TTree *>(in->Get("Data_R"));
  if (!tree) {
    std::cerr << "Cannot find TTree Data_R in " << inputFile << std::endl;
    return false;
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

  const Long64_t nEntries = tree->GetEntries();
  hits.reserve(hits.size() + static_cast<size_t>(nEntries));
  Long64_t kept = 0;
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
    hits.push_back(hit);
    ++kept;
  }

  std::cout << "Read " << inputFile << ": kept " << kept << " / " << nEntries << " hits" << std::endl;
  in->Close();
  return true;
}

static std::pair<double, double> FitPeakTwoStage(TH1D &hist,
                                                 double axisMin,
                                                 double axisMax,
                                                 const char *fitName,
                                                 bool storeFineFit = false) {
  if (hist.GetEntries() < 20) {
    const double peak = hist.GetMean();
    const double sigma = std::max(hist.GetRMS(), 0.5 * hist.GetXaxis()->GetBinWidth(1));
    return {peak, sigma};
  }

  const int maxBin = hist.GetMaximumBin();
  const double peakSeed = hist.GetBinCenter(maxBin);
  const double histMax = hist.GetMaximum();
  const double roughHalf = std::max(50.0, 0.05 * std::max(std::abs(peakSeed), 1.0));
  const double roughMin = std::max(axisMin, peakSeed - roughHalf);
  const double roughMax = std::min(axisMax, peakSeed + roughHalf);

  TF1 rough(Form("%s_rough", fitName), "gaus", roughMin, roughMax);
  rough.SetParameters(histMax, peakSeed, roughHalf / 4.0);
  hist.Fit(&rough, "QNR");

  double peak = rough.GetParameter(1);
  double sigma = std::abs(rough.GetParameter(2));
  if (!(sigma > 0.0) || !std::isfinite(sigma) || !std::isfinite(peak)) {
    peak = peakSeed;
    sigma = roughHalf / 4.0;
  }

  const double fineHalf = std::max(3.0 * sigma, 2.0 * hist.GetXaxis()->GetBinWidth(1));
  const double fineMin = std::max(axisMin, peak - fineHalf);
  const double fineMax = std::min(axisMax, peak + fineHalf);

  TF1 fine(fitName, "gaus", fineMin, fineMax);
  fine.SetParameters(histMax, peak, sigma);
  hist.Fit(&fine, storeFineFit ? "QR" : "QNR");
  peak = fine.GetParameter(1);
  sigma = std::abs(fine.GetParameter(2));
  if (!(sigma > 0.0) || !std::isfinite(sigma)) {
    sigma = fineHalf / 3.0;
  }
  return {peak, sigma};
}

static void SetEnergyDisplayRange(TH1D &hist, double peakKeV, double sigmaKeV) {
  if (hist.GetEntries() < 1 || !(peakKeV > 0.0)) {
    return;
  }
  const double binWidth = hist.GetXaxis()->GetBinWidth(1);
  const double sigma = std::max(sigmaKeV, binWidth);
  const double fwhm = 2.354820045 * sigma;
  const double halfWidth = std::max(150.0, 4.0 * fwhm);

  int binLo = hist.FindBin(peakKeV - halfWidth);
  int binHi = hist.FindBin(peakKeV + halfWidth);
  binLo = std::max(1, binLo);
  binHi = std::min(hist.GetNbinsX(), binHi);
  if (binHi <= binLo) {
    binLo = std::max(1, hist.FindBin(peakKeV) - 75);
    binHi = std::min(hist.GetNbinsX(), hist.FindBin(peakKeV) + 75);
  }
  // Use bin-index range so displayed edges stay on true bin boundaries.
  hist.GetXaxis()->SetRange(binLo, binHi);
}

static void WriteSummaryPdf(const char *pdfName,
                            TH2D &hHitMap,
                            TH1D &hTimeDiff,
                            TH2D &hPeak,
                            TH2D &hResolution,
                            TH2D &hRelDev,
                            TH1D &hPeakDist,
                            TH1D &hAllEnergy,
                            double allEnergyPeakKeV,
                            double allEnergySigmaKeV,
                            std::array<std::array<PixelSpectrum, kNY>, kNX> &pixelSummary) {
  TCanvas c("c_uniformity", "Pixel energy uniformity", 1400, 1000);
  gStyle->SetOptStat(1);
  gStyle->SetPaintTextFormat("g");

  c.Print(Form("%s[", pdfName));

  c.Clear();
  gStyle->SetOptStat(0);
  hHitMap.Draw("COLZ TEXT");
  c.Print(pdfName);

  c.Clear();
  gStyle->SetOptStat(1);
  hTimeDiff.Draw();
  c.Print(pdfName);

  // Per-pixel calibrated energy spectra [keV], one X-strip page with 4x4 Y pads.
  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(111);
  for (int ix = 0; ix < kNX; ++ix) {
    c.Clear();
    c.Divide(4, 4);
    for (int iy = 0; iy < kNY; ++iy) {
      c.cd(iy + 1);
      if (!pixelSummary[ix][iy].hist || pixelSummary[ix][iy].hist->GetEntries() < 1) {
        continue;
      }
      // Draw a clone so display zoom does not mutate the fill histogram.
      auto *drawHist = static_cast<TH1D *>(pixelSummary[ix][iy].hist->Clone(Form("hPixDraw_%02d_%02d", ix, iy)));
      drawHist->SetDirectory(nullptr);
      drawHist->SetTitle(Form("Pixel X=%02d Y=%02d;Energy [keV];Counts", ix, iy));
      drawHist->GetXaxis()->SetTitle("Energy [keV]");
      drawHist->GetXaxis()->SetNoExponent();
      if (pixelSummary[ix][iy].fitted && pixelSummary[ix][iy].peak > 0.0) {
        SetEnergyDisplayRange(*drawHist, pixelSummary[ix][iy].peak, pixelSummary[ix][iy].sigma);
      }
      drawHist->Draw();
    }
    c.Print(pdfName);
  }

  c.Clear();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  hPeak.Draw("COLZ TEXT");
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
  gStyle->SetPaintTextFormat("4.2f%%");
  hRelDev.Draw("COLZ TEXT");
  c.Print(pdfName);

  c.Clear();
  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(111);
  gStyle->SetPaintTextFormat("g");
  hPeakDist.Draw();
  c.Print(pdfName);

  c.Clear();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(1111);
  auto *allDraw = static_cast<TH1D *>(hAllEnergy.Clone("hAllEnergyDraw"));
  allDraw->SetDirectory(nullptr);
  allDraw->GetXaxis()->SetTitle("Energy [keV]");
  allDraw->GetXaxis()->SetNoExponent();
  SetEnergyDisplayRange(*allDraw, allEnergyPeakKeV, allEnergySigmaKeV);
  allDraw->Draw();
  c.Print(pdfName);

  c.Print(Form("%s]", pdfName));
}

int check_pixel_uniformity(const char *calibFile = "CalPixE/pixel_calib.txt",
                           const char *inputFiles = "current.root",
                           double coincidenceGatePs = 0.0,
                           double adcThreshold = kDefaultAdcThreshold,
                           double energyMinKeV = kDefaultEnergyMinKeV,
                           double energyMaxKeV = kDefaultEnergyMaxKeV,
                           double energyBinWidthKeV = kDefaultEnergyBinWidthKeV,
                           Long64_t minPixelEntries = kDefaultMinPixelEntries,
                           const char *channelMapFile = "channel_map.txt",
                           const char *summaryPdf = "UnifPix/pixel_uniformity.pdf") {
  TH1::AddDirectory(kFALSE);

  if (!(energyMaxKeV > energyMinKeV) || !(energyBinWidthKeV > 0.0)) {
    std::cerr << "Invalid energy axis" << std::endl;
    return 1;
  }
  if (minPixelEntries < 1) {
    std::cerr << "minPixelEntries must be >= 1" << std::endl;
    return 1;
  }
  std::cout << "Minimum entries per pixel: " << minPixelEntries << std::endl;

  StripMap stripMap;
  if (!LoadStripMap(channelMapFile, stripMap)) {
    return 6;
  }

  std::array<std::array<PixelCalib, kNY>, kNX> calib{};
  if (!LoadPixelCalib(calibFile, calib)) {
    return 2;
  }

  const auto files = SplitCommaSeparated(inputFiles);
  if (files.empty()) {
    std::cerr << "No input ROOT files given" << std::endl;
    return 3;
  }

  std::vector<Hit> hits;
  for (const auto &file : files) {
    if (!AppendHitsFromFile(file.c_str(), adcThreshold, stripMap, hits)) {
      return 4;
    }
  }
  if (hits.empty()) {
    std::cerr << "No hits kept after ADC threshold cut" << std::endl;
    return 5;
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

  const int nEnergyBins =
      std::max(1, static_cast<int>(std::llround((energyMaxKeV - energyMinKeV) / energyBinWidthKeV)));
  const double energyAxisMax = energyMinKeV + nEnergyBins * energyBinWidthKeV;

  std::array<std::array<PixelSpectrum, kNY>, kNX> pixelSummary;
  for (int ix = 0; ix < kNX; ++ix) {
    for (int iy = 0; iy < kNY; ++iy) {
      pixelSummary[ix][iy].hist = std::make_unique<TH1D>(
          Form("hPixE_%02d_%02d", ix, iy), Form("Pixel X=%02d Y=%02d energy [keV]", ix, iy), nEnergyBins,
          energyMinKeV, energyAxisMax);
    }
  }

  TH2D hHitMap("hHitMap", "Pixel hit map;X;Y", kNX, -0.5, kNX - 0.5, kNY, -0.5, kNY - 0.5);
  TH1D hTimeDiff("hTimeDiff", "Nearest X-Y timestamp difference;|#Delta t| [ps];Counts", 400, 0,
                 kTimeDiffHistogramMaxPs);
  TH1D hAllEnergy("hAllEnergy", "All-pixel calibrated X energy;Energy [keV];Counts", nEnergyBins, energyMinKeV,
                  energyAxisMax);

  for (const auto &xHit : xHits) {
    ULong64_t nearestDt = 0;
    if (FindNearestTimestampDifference(yHits, xHit.timestamp, nearestDt)) {
      hTimeDiff.Fill(static_cast<double>(nearestDt));
    }
  }

  const double fittedGate = DetermineCoincidenceGatePs(hTimeDiff);
  const double gatePs = coincidenceGatePs > 0.0 ? coincidenceGatePs : fittedGate;
  std::cout << "Using coincidence gate: " << gatePs << " ps" << std::endl;

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
    PruneOldHits(xBuffer, hit.timestamp, gatePs);
    PruneOldHits(yBuffer, hit.timestamp, gatePs);

    const int matchIndex = FindClosestHitIndex(otherBuffer, hit.timestamp, gatePs);
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

  Long64_t used = 0;
  Long64_t skippedNoCalib = 0;
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

    const int ix = xMapped.strip;
    const int iy = yMapped.strip;
    if (!calib[ix][iy].valid) {
      ++skippedNoCalib;
      continue;
    }
    const double energyKeV = calib[ix][iy].offset + calib[ix][iy].slope * xHit.energy;
    hHitMap.Fill(ix, iy);
    pixelSummary[ix][iy].hist->Fill(energyKeV);
    hAllEnergy.Fill(energyKeV);
    ++used;
  }

  std::cout << "Coincidences used with calibration: " << used << std::endl;
  if (skippedNoCalib > 0) {
    std::cout << "Skipped coincidences without calib: " << skippedNoCalib << std::endl;
  }

  Long64_t skippedLowStats = 0;
  for (int ix = 0; ix < kNX; ++ix) {
    for (int iy = 0; iy < kNY; ++iy) {
      auto *hist = pixelSummary[ix][iy].hist.get();
      if (hist->GetEntries() < minPixelEntries) {
        if (hist->GetEntries() > 0) {
          ++skippedLowStats;
        }
        continue;
      }

      if (hist->GetEntries() < 20) {
        pixelSummary[ix][iy].peak = hist->GetMean();
        const double minSigma = 0.5 * hist->GetXaxis()->GetBinWidth(1);
        pixelSummary[ix][iy].sigma = std::max(hist->GetRMS(), minSigma);
      } else {
        const auto fitted =
            FitPeakTwoStage(*hist, energyMinKeV, energyAxisMax, Form("pixel_gaus_%02d_%02d", ix, iy), true);
        pixelSummary[ix][iy].peak = fitted.first;
        pixelSummary[ix][iy].sigma = fitted.second;
      }
      pixelSummary[ix][iy].fwhm = 2.354820045 * pixelSummary[ix][iy].sigma;
      pixelSummary[ix][iy].fitted = true;
    }
  }
  if (skippedLowStats > 0) {
    std::cout << "Skipped pixels below entry threshold: " << skippedLowStats << std::endl;
  }

  TH2D hPeak("hPixelPeak", "Calibrated pixel peak energy [keV];X;Y", kNX, -0.5, kNX - 0.5, kNY, -0.5,
             kNY - 0.5);
  TH2D hResolution("hResolution", "Pixel resolution (FWHM/peak);X;Y", kNX, -0.5, kNX - 0.5, kNY, -0.5,
                   kNY - 0.5);
  TH2D hRelDev("hRelDev", "Peak relative deviation from mean [%];X;Y", kNX, -0.5, kNX - 0.5, kNY, -0.5,
               kNY - 0.5);

  std::vector<double> peaks;
  peaks.reserve(kNX * kNY);
  for (int ix = 0; ix < kNX; ++ix) {
    for (int iy = 0; iy < kNY; ++iy) {
      if (!pixelSummary[ix][iy].fitted || !(pixelSummary[ix][iy].peak > 0.0)) {
        continue;
      }
      const double peak = pixelSummary[ix][iy].peak;
      peaks.push_back(peak);
      hPeak.SetBinContent(ix + 1, iy + 1, peak);
      hResolution.SetBinContent(ix + 1, iy + 1, pixelSummary[ix][iy].fwhm / peak);
    }
  }

  double peakMean = 0.0;
  double peakRms = 0.0;
  double allEnergyPeakKeV = 0.0;
  double allEnergySigmaKeV = 50.0;
  std::unique_ptr<TH1D> hPeakDist;
  if (!peaks.empty()) {
    for (double peak : peaks) {
      peakMean += peak;
    }
    peakMean /= static_cast<double>(peaks.size());
    for (double peak : peaks) {
      const double d = peak - peakMean;
      peakRms += d * d;
    }
    peakRms = std::sqrt(peakRms / static_cast<double>(peaks.size()));

    double sigmaSum = 0.0;
    int sigmaCount = 0;
    for (int ix = 0; ix < kNX; ++ix) {
      for (int iy = 0; iy < kNY; ++iy) {
        if (!pixelSummary[ix][iy].fitted || !(pixelSummary[ix][iy].sigma > 0.0)) {
          continue;
        }
        sigmaSum += pixelSummary[ix][iy].sigma;
        ++sigmaCount;
      }
    }
    if (sigmaCount > 0) {
      allEnergySigmaKeV = sigmaSum / static_cast<double>(sigmaCount);
    }
    allEnergyPeakKeV = peakMean;

    const auto minmax = std::minmax_element(peaks.begin(), peaks.end());
    const double peakMin = *minmax.first;
    const double peakMax = *minmax.second;
    const double span = std::max(peakMax - peakMin, 1.0);
    const double halfWidth = std::max({2.0, 5.0 * peakRms, 0.75 * span, 0.002 * peakMean});
    const double distMin = peakMean - halfWidth;
    const double distMax = peakMean + halfWidth;
    const double binWidth = std::max(0.05, halfWidth / 40.0);
    const int nDistBins = std::max(20, static_cast<int>(std::llround((distMax - distMin) / binWidth)));

    hPeakDist = std::make_unique<TH1D>("hPeakDist", "Distribution of pixel peak energies;Peak energy [keV];Pixels",
                                       nDistBins, distMin, distMax);
    for (double peak : peaks) {
      hPeakDist->Fill(peak);
    }

    // Rough then fine Gaussian fit of the pixel-peak distribution.
    FitPeakTwoStage(*hPeakDist, distMin, distMax, "hPeakDist_gaus", true);

    for (int ix = 0; ix < kNX; ++ix) {
      for (int iy = 0; iy < kNY; ++iy) {
        if (!pixelSummary[ix][iy].fitted || !(pixelSummary[ix][iy].peak > 0.0) || !(peakMean > 0.0)) {
          continue;
        }
        hRelDev.SetBinContent(ix + 1, iy + 1, 100.0 * (pixelSummary[ix][iy].peak - peakMean) / peakMean);
      }
    }
  } else {
    hPeakDist = std::make_unique<TH1D>("hPeakDist", "Distribution of pixel peak energies;Peak energy [keV];Pixels",
                                       100, energyMinKeV, energyAxisMax);
    if (hAllEnergy.GetEntries() > 0) {
      allEnergyPeakKeV = hAllEnergy.GetBinCenter(hAllEnergy.GetMaximumBin());
    }
  }

  WriteSummaryPdf(summaryPdf, hHitMap, hTimeDiff, hPeak, hResolution, hRelDev, *hPeakDist, hAllEnergy,
                  allEnergyPeakKeV, allEnergySigmaKeV, pixelSummary);

  std::cout << "Wrote " << summaryPdf << std::endl;
  std::cout << "Pixels with fitted peaks: " << peaks.size() << std::endl;
  if (!peaks.empty()) {
    const auto minmax = std::minmax_element(peaks.begin(), peaks.end());
    std::cout << "Peak mean = " << peakMean << " keV" << std::endl;
    std::cout << "Peak RMS  = " << peakRms << " keV (" << (100.0 * peakRms / peakMean) << " %)" << std::endl;
    std::cout << "Peak min/max = " << *minmax.first << " / " << *minmax.second << " keV" << std::endl;
    if (auto *fit = hPeakDist->GetFunction("hPeakDist_gaus")) {
      std::cout << "hPeakDist fit mean/sigma = " << fit->GetParameter(1) << " / " << std::abs(fit->GetParameter(2))
                << " keV" << std::endl;
    }
  }
  return 0;
}
