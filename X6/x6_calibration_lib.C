#include "x6_calibration.h"

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace X6Calib {

static bool ParseMappingLine(const char *line,
                             int &board,
                             int &channel,
                             int &detector,
                             char &side,
                             int &detAna,
                             int &strip,
                             int &pin) {
  while (*line != '\0' && std::isspace(static_cast<unsigned char>(*line))) {
    ++line;
  }
  if (*line == '\0' || *line == '#') {
    return false;
  }

  char sideToken[8] = {};
  const int nRead = std::sscanf(line, "%d %d %d %7s %d %d %d", &board, &channel,
                                &detector, sideToken, &detAna, &strip, &pin);
  if (nRead != 7) {
    return false;
  }

  side = sideToken[0];
  return side == 'J' || side == 'O';
}

static long long ChannelKey(int board, int channel) {
  return static_cast<long long>(board) * 1000 + channel;
}

static std::vector<DetectorMap> LoadDetectorMapFromFile(const TString &path,
                                                        std::map<long long, ChannelInfo> &channelLookup) {
  FILE *fp = fopen(path.Data(), "r");
  if (fp == nullptr) {
    std::cerr << "Could not open mapping file: " << path << "\n";
    return {};
  }

  using StripKey = std::tuple<int, int>;
  std::map<int, DetectorMap> detectorGroups;
  std::map<StripKey, std::map<int, int>> stripPins;
  channelLookup.clear();

  char buffer[512];
  int lineNumber = 0;
  while (fgets(buffer, sizeof(buffer), fp) != nullptr) {
    ++lineNumber;

    int board = 0;
    int channel = 0;
    int detector = 0;
    char side = '\0';
    int detAna = 0;
    int strip = 0;
    int pin = 0;
    if (!ParseMappingLine(buffer, board, channel, detector, side, detAna, strip,
                          pin)) {
      continue;
    }

    (void)detAna;

    if (detector < 0 || detector >= kNumDetectors) {
      std::cerr << "Invalid detector index " << detector << " at line "
                << lineNumber << " in mapping file: " << path << "\n";
      fclose(fp);
      return {};
    }

    if (detectorGroups.find(detector) == detectorGroups.end()) {
      DetectorMap dm;
      dm.board = board;
      dm.detector = detector;
      for (int i = 0; i < kNumOhmic; ++i) {
        dm.ohmic[i].ohmicIndex = i;
        dm.ohmic[i].digiChannel = -1;
      }
      detectorGroups[detector] = dm;
    } else if (detectorGroups[detector].board != board) {
      std::cerr << "Detector " << detector << " mapped to multiple boards at line "
                << lineNumber << " in mapping file: " << path << "\n";
      fclose(fp);
      return {};
    }

    ChannelInfo info;
    info.board = board;
    info.digiChannel = channel;
    info.detector = detector;
    channelLookup[ChannelKey(board, channel)] = info;

    if (side == 'O') {
      if (strip < 0 || strip >= kNumOhmic || pin != -1) {
        std::cerr << "Invalid ohmic strip/pin at line " << lineNumber
                  << " in mapping file: " << path << "\n";
        fclose(fp);
        return {};
      }

      detectorGroups[detector].ohmic[strip].ohmicIndex = strip;
      detectorGroups[detector].ohmic[strip].digiChannel = channel;
      info.isOhmic = true;
      info.ohmicIndex = strip;
      channelLookup[ChannelKey(board, channel)] = info;
      continue;
    }

    if (strip < 0 || strip >= kNumStrips || (pin != 0 && pin != 1)) {
      std::cerr << "Invalid junction strip/pin at line " << lineNumber
                << " in mapping file: " << path << "\n";
      fclose(fp);
      return {};
    }

    info.isOhmic = false;
    info.strip = strip;
    info.pin = pin;
    channelLookup[ChannelKey(board, channel)] = info;
    stripPins[StripKey(detector, strip)][pin] = channel;
  }
  fclose(fp);

  for (const auto &entry : stripPins) {
    const int detector = std::get<0>(entry.first);
    const int strip = std::get<1>(entry.first);
    const auto &pins = entry.second;

    const auto pin0 = pins.find(0);
    const auto pin1 = pins.find(1);
    if (pin0 == pins.end() || pin1 == pins.end()) {
      std::cerr << "Incomplete junction strip mapping for detector=" << detector
                << " strip=" << strip << "\n";
      continue;
    }

    JunctionPair jp;
    jp.board = detectorGroups[detector].board;
    jp.detector = detector;
    jp.stripIndex = strip;
    jp.chA = pin0->second;
    jp.chB = pin1->second;
    jp.pinA = 0;
    jp.pinB = 1;
    detectorGroups[detector].pairs.push_back(jp);
  }

  std::vector<DetectorMap> maps;
  maps.reserve(detectorGroups.size());
  for (auto &entry : detectorGroups) {
    auto &dm = entry.second;
    std::sort(dm.pairs.begin(), dm.pairs.end(),
              [](const JunctionPair &a, const JunctionPair &b) {
                return a.stripIndex < b.stripIndex;
              });
    maps.push_back(dm);
  }

  std::sort(maps.begin(), maps.end(),
            [](const DetectorMap &a, const DetectorMap &b) {
              return a.detector < b.detector;
            });

  std::cout << "Loaded " << maps.size() << " detector group(s) from mapping file: "
            << path << "\n";
  return maps;
}

static const DetectorMap *FindDetectorMap(
    const std::vector<DetectorMap> &detectorMaps, int detector) {
  for (const auto &dm : detectorMaps) {
    if (dm.detector == detector) {
      return &dm;
    }
  }
  return nullptr;
}

static const ChannelInfo *FindChannelInfo(
    const std::map<long long, ChannelInfo> &channelLookup, int board, int channel) {
  const auto it = channelLookup.find(ChannelKey(board, channel));
  if (it == channelLookup.end()) {
    return nullptr;
  }
  return &it->second;
}

static void LoadChain(TChain &chain, const TString &pattern) {
  const int nAdded = chain.Add(pattern);
  std::cout << "Added " << nAdded << " file(s) from pattern: " << pattern << "\n";
  if (nAdded <= 0) {
    std::cerr << "No ROOT files matched the requested pattern.\n";
  }
}

static std::vector<Event> ReadEvents(TChain &chain, bool useEnergyShort = false) {
  std::vector<Event> events;
  events.reserve(chain.GetEntries());

  TTreeReader reader(&chain);
  TTreeReaderValue<UShort_t> board(reader, "Board");
  TTreeReaderValue<UShort_t> channel(reader, "Channel");
  TTreeReaderValue<ULong64_t> timestamp(reader, "Timestamp");
  TTreeReaderValue<UShort_t> energy(reader, "Energy");
  TTreeReaderValue<UShort_t> energyShort(reader, "EnergyShort");
  TTreeReaderValue<UInt_t> flags(reader, "Flags");

  while (reader.Next()) {
    Event ev;
    ev.board = *board;
    ev.channel = *channel;
    ev.timestamp = *timestamp;
    ev.energy = *energy;
    ev.energyShort = *energyShort;
    ev.flags = *flags;
    if (useEnergyShort) {
      ev.energy = ev.energyShort;
    }
    events.push_back(ev);
  }

  std::sort(events.begin(), events.end(),
            [](const Event &a, const Event &b) { return a.timestamp < b.timestamp; });

  std::cout << "Loaded " << events.size() << " entries and time-sorted them.\n";
  return events;
}

static bool IsValidPulse(const Event &ev, double adcThreshold) {
  // Require ADC above threshold; skip saturated full-scale entries.
  const double thr = std::max(0.0, adcThreshold);
  return (static_cast<double>(ev.energy) >= thr && ev.energy < 4095);
}

static void CollectTimingDiagnostics(
    const std::vector<Event> &events,
    const std::vector<DetectorMap> &detectorMaps,
    const std::map<long long, ChannelInfo> &channelLookup,
    TimingDiagnostics &timing,
    double adcThreshold = 50.0,
    ULong64_t probeWindow = 400000,
    std::size_t maxSamples = 100000) {
  timing.dtJunctionPair.clear();
  timing.dtOhmicPin0.clear();
  timing.dtOhmicPin1.clear();
  timing.dtOhmicJun.clear();
  timing.dtPin0Pin1.clear();
  timing.maxSpan.clear();
  timing.dtJunctionPair.reserve(maxSamples);
  timing.dtOhmicPin0.reserve(maxSamples);
  timing.dtOhmicPin1.reserve(maxSamples);
  timing.dtOhmicJun.reserve(2 * maxSamples);
  timing.dtPin0Pin1.reserve(maxSamples);
  timing.maxSpan.reserve(maxSamples);

  if (events.empty()) {
    return;
  }

  std::size_t left = 0;

  for (std::size_t i = 0; i < events.size(); ++i) {
    const Event &ev = events[i];
    if (!IsValidPulse(ev, adcThreshold)) {
      continue;
    }

    const ChannelInfo *info =
        FindChannelInfo(channelLookup, ev.board, ev.channel);
    if (info == nullptr) {
      continue;
    }

    const DetectorMap *dm = FindDetectorMap(detectorMaps, info->detector);
    if (dm == nullptr || dm->board != ev.board) {
      continue;
    }

    while (left < i && ev.timestamp > events[left].timestamp &&
           ev.timestamp - events[left].timestamp > probeWindow) {
      ++left;
    }

    std::size_t right = i;
    while (right + 1 < events.size() &&
           events[right + 1].timestamp >= ev.timestamp &&
           events[right + 1].timestamp - ev.timestamp <= probeWindow) {
      ++right;
    }

    if (!info->isOhmic) {
      // pin0-pin1 pair on one strip. Either pin may arrive first.
      for (const auto &pair : dm->pairs) {
        if (ev.channel != pair.chA && ev.channel != pair.chB) {
          continue;
        }

        const Event *partner = nullptr;
        for (std::size_t j = left; j < i; ++j) {
          const Event &cand = events[j];
          if (cand.board != ev.board || !IsValidPulse(cand, adcThreshold)) {
            continue;
          }
          const ChannelInfo *candInfo =
              FindChannelInfo(channelLookup, cand.board, cand.channel);
          if (candInfo == nullptr || candInfo->detector != info->detector ||
              candInfo->isOhmic) {
            continue;
          }

          const bool matchedPair =
              (cand.channel == pair.chA && ev.channel == pair.chB) ||
              (cand.channel == pair.chB && ev.channel == pair.chA);
          if (!matchedPair) {
            continue;
          }

          if (partner == nullptr ||
              std::abs(static_cast<double>(ev.timestamp) -
                       static_cast<double>(cand.timestamp)) <
                  std::abs(static_cast<double>(ev.timestamp) -
                           static_cast<double>(partner->timestamp))) {
            partner = &cand;
          }
        }

        if (partner != nullptr) {
          timing.dtJunctionPair.push_back(
              std::abs(static_cast<double>(ev.timestamp) -
                       static_cast<double>(partner->timestamp)));
        }
      }
    } else {
      // Fill only when both pins of one strip are found near this ohmic.
      // Arrival order of ohmic/pin0/pin1 does not matter.
      for (const auto &pair : dm->pairs) {
        const Event *hitPin0 = nullptr;
        const Event *hitPin1 = nullptr;

        for (std::size_t j = left; j <= right; ++j) {
          if (j == i) {
            continue;
          }

          const Event &cand = events[j];
          if (cand.board != ev.board || !IsValidPulse(cand, adcThreshold)) {
            continue;
          }
          const ChannelInfo *candInfo =
              FindChannelInfo(channelLookup, cand.board, cand.channel);
          if (candInfo == nullptr || candInfo->detector != info->detector ||
              candInfo->isOhmic) {
            continue;
          }

          if (cand.channel == pair.chA) {
            if (hitPin0 == nullptr ||
                std::abs(static_cast<double>(ev.timestamp) -
                         static_cast<double>(cand.timestamp)) <
                    std::abs(static_cast<double>(ev.timestamp) -
                             static_cast<double>(hitPin0->timestamp))) {
              hitPin0 = &cand;
            }
          } else if (cand.channel == pair.chB) {
            if (hitPin1 == nullptr ||
                std::abs(static_cast<double>(ev.timestamp) -
                         static_cast<double>(cand.timestamp)) <
                    std::abs(static_cast<double>(ev.timestamp) -
                             static_cast<double>(hitPin1->timestamp))) {
              hitPin1 = &cand;
            }
          }
        }

        if (hitPin0 == nullptr || hitPin1 == nullptr) {
          continue;
        }

        const double dt0 =
            std::abs(static_cast<double>(ev.timestamp) -
                     static_cast<double>(hitPin0->timestamp));
        const double dt1 =
            std::abs(static_cast<double>(ev.timestamp) -
                     static_cast<double>(hitPin1->timestamp));
        const double dt01 =
            std::abs(static_cast<double>(hitPin0->timestamp) -
                     static_cast<double>(hitPin1->timestamp));

        timing.dtOhmicPin0.push_back(dt0);
        timing.dtOhmicPin1.push_back(dt1);
        timing.dtOhmicJun.push_back(dt0);
        timing.dtOhmicJun.push_back(dt1);
        timing.dtPin0Pin1.push_back(dt01);
        timing.maxSpan.push_back(std::max(dt01, std::max(dt0, dt1)));
      }
    }

    if (timing.dtOhmicJun.size() >= 2 * maxSamples &&
        timing.dtJunctionPair.size() >= maxSamples) {
      break;
    }
  }

  std::cout << "Timing diagnostics collected:\n"
            << "  |t_pin0 - t_pin1|           : " << timing.dtJunctionPair.size() << "\n"
            << "  |t_ohmic - t_pin0|          : " << timing.dtOhmicPin0.size() << "\n"
            << "  |t_ohmic - t_pin1|          : " << timing.dtOhmicPin1.size() << "\n"
            << "  |t_ohmic - t_pin| (both)    : " << timing.dtOhmicJun.size() << "\n"
            << "  3-hit strip candidates     : " << timing.maxSpan.size() << "\n";
}

static TF1 *FitMainPeak(TH1D *hist, const char *fitName) {
  if (hist == nullptr || hist->GetEntries() < 10) {
    return nullptr;
  }

  const int peakBin = hist->GetMaximumBin();
  const double peakCenter = hist->GetBinCenter(peakBin);
  const double xMin = hist->GetXaxis()->GetXmin();
  const double xMax = hist->GetXaxis()->GetXmax();
  const double rms = std::max(hist->GetRMS(), hist->GetBinWidth(1) * 3.0);
  const double fitHalfWidth = std::max(20.0, 1.5 * rms);
  const double fitMin = std::max(xMin, peakCenter - fitHalfWidth);
  const double fitMax = std::min(xMax, peakCenter + fitHalfWidth);

  TF1 *fit = new TF1(fitName, "gaus", fitMin, fitMax);
  fit->SetParameters(hist->GetMaximum(), peakCenter, std::max(5.0, 0.5 * rms));
  // |dt| histograms are non-negative; keep the mean in range.
  fit->SetParLimits(1, xMin, xMax);
  fit->SetParLimits(2, 0.0, xMax);
  hist->Fit(fit, "RQ0");
  return fit;
}

// Find up to nPeaks local maxima by iterative masking, sorted low→high in ADC.
static std::vector<double> FindPeakCentersAdc(TH1D *hist,
                                              int nPeaks,
                                              double minSeparationAdc) {
  std::vector<double> centers;
  if (hist == nullptr || nPeaks <= 0 || hist->GetEntries() < 10) {
    return centers;
  }

  TH1D *work = static_cast<TH1D *>(hist->Clone(Form("%s_peakfind", hist->GetName())));
  work->SetDirectory(nullptr);

  const double xMin = hist->GetXaxis()->GetXmin();
  const double xMax = hist->GetXaxis()->GetXmax();
  const double binW = hist->GetBinWidth(1);
  const double sep = std::max(minSeparationAdc, 3.0 * binW);
  const double minHeight =
      std::max(3.0, 0.02 * hist->GetMaximum());  // ignore tiny noise bumps

  for (int ip = 0; ip < nPeaks; ++ip) {
    const int peakBin = work->GetMaximumBin();
    const double height = work->GetBinContent(peakBin);
    if (height < minHeight) {
      break;
    }
    const double center = work->GetBinCenter(peakBin);
    centers.push_back(center);

    const double lo = std::max(xMin, center - sep);
    const double hi = std::min(xMax, center + sep);
    const int binLo = work->FindBin(lo);
    const int binHi = work->FindBin(hi);
    for (int b = binLo; b <= binHi; ++b) {
      work->SetBinContent(b, 0.0);
    }
  }

  delete work;
  std::sort(centers.begin(), centers.end());
  return centers;
}

static bool FitLinearAdcVsEnergy(const std::vector<double> &meansAdc,
                                 const std::vector<double> &energiesMeV,
                                 double &calibA,
                                 double &calibB) {
  calibA = 0.0;
  calibB = 0.0;
  if (meansAdc.size() != energiesMeV.size() || meansAdc.empty()) {
    return false;
  }

  if (meansAdc.size() == 1) {
    // Single point: force intercept 0 → ADC = B * E
    if (energiesMeV[0] == 0.0) {
      return false;
    }
    calibA = 0.0;
    calibB = meansAdc[0] / energiesMeV[0];
    return calibB != 0.0;
  }

  TGraph gr(static_cast<int>(meansAdc.size()));
  for (std::size_t i = 0; i < meansAdc.size(); ++i) {
    gr.SetPoint(static_cast<int>(i), energiesMeV[i], meansAdc[i]);
  }

  TF1 line("calib_line", "pol1", energiesMeV.front(), energiesMeV.back());
  const int status = gr.Fit(&line, "Q0");
  if (status != 0) {
    return false;
  }
  calibA = line.GetParameter(0);
  calibB = line.GetParameter(1);
  return true;
}

// Fit nPeaks Gaussians (ordered low→high ADC), assign sorted MeV energies,
// and build ADC = A + B*E calibration.
static TF1 *FitCalibrationPeaks(TH1D *hist,
                                const char *fitName,
                                int nPeaks,
                                const std::vector<double> &energiesMeVIn,
                                PeakFitResult &result) {
  result.meansAdc.clear();
  result.sigmasAdc.clear();
  result.energiesMeV.clear();
  result.nPeaksFitted = 0;
  result.fitted = false;
  result.calibrated = false;
  result.calibA = 0.0;
  result.calibB = 0.0;

  if (hist == nullptr || nPeaks <= 0 || hist->GetEntries() < 10) {
    return nullptr;
  }

  std::vector<double> energies = energiesMeVIn;
  if (static_cast<int>(energies.size()) < nPeaks) {
    energies.resize(static_cast<std::size_t>(nPeaks), 0.0);
  } else if (static_cast<int>(energies.size()) > nPeaks) {
    energies.resize(static_cast<std::size_t>(nPeaks));
  }
  std::sort(energies.begin(), energies.end());

  const double rms = std::max(hist->GetRMS(), hist->GetBinWidth(1) * 3.0);
  const double minSep = std::max(30.0, 0.5 * rms);
  const std::vector<double> centers = FindPeakCentersAdc(hist, nPeaks, minSep);
  if (centers.empty()) {
    return nullptr;
  }

  const int nUse = static_cast<int>(centers.size());
  const double xMin = hist->GetXaxis()->GetXmin();
  const double xMax = hist->GetXaxis()->GetXmax();

  // Build sum-of-Gaussians formula: gaus(0)+gaus(3)+...
  TString formula = "gaus(0)";
  for (int i = 1; i < nUse; ++i) {
    formula += Form("+gaus(%d)", 3 * i);
  }

  // Fit range covering outermost peaks with some margin.
  const double halfW = std::max(25.0, 1.2 * rms);
  const double fitMin = std::max(xMin, centers.front() - halfW);
  const double fitMax = std::min(xMax, centers.back() + halfW);

  TF1 *fit = new TF1(fitName, formula.Data(), fitMin, fitMax);
  for (int i = 0; i < nUse; ++i) {
    const int peakBin = hist->FindBin(centers[static_cast<std::size_t>(i)]);
    const double amp = std::max(1.0, hist->GetBinContent(peakBin));
    const double sigma0 = std::max(5.0, 0.15 * rms);
    fit->SetParameter(3 * i + 0, amp);
    fit->SetParameter(3 * i + 1, centers[static_cast<std::size_t>(i)]);
    fit->SetParameter(3 * i + 2, sigma0);
    fit->SetParLimits(3 * i + 0, 0.0, 10.0 * hist->GetMaximum());
    fit->SetParLimits(3 * i + 1, xMin, xMax);
    fit->SetParLimits(3 * i + 2, hist->GetBinWidth(1) * 0.5, 0.5 * (xMax - xMin));
  }

  hist->Fit(fit, "RQ0");

  result.meansAdc.reserve(static_cast<std::size_t>(nUse));
  result.sigmasAdc.reserve(static_cast<std::size_t>(nUse));
  result.energiesMeV.reserve(static_cast<std::size_t>(nUse));
  for (int i = 0; i < nUse; ++i) {
    result.meansAdc.push_back(fit->GetParameter(3 * i + 1));
    result.sigmasAdc.push_back(std::abs(fit->GetParameter(3 * i + 2)));
    result.energiesMeV.push_back(energies[static_cast<std::size_t>(i)]);
  }
  // Keep ADC/energy pairs ordered by ADC after the fit (means can swap slightly).
  std::vector<std::size_t> order(static_cast<std::size_t>(nUse));
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
    return result.meansAdc[a] < result.meansAdc[b];
  });
  std::vector<double> meansSorted;
  std::vector<double> sigmasSorted;
  meansSorted.reserve(order.size());
  sigmasSorted.reserve(order.size());
  for (std::size_t idx : order) {
    meansSorted.push_back(result.meansAdc[idx]);
    sigmasSorted.push_back(result.sigmasAdc[idx]);
  }
  result.meansAdc.swap(meansSorted);
  result.sigmasAdc.swap(sigmasSorted);
  // Energies were sorted low→high and assigned to low→high ADC peaks.
  result.energiesMeV.assign(energies.begin(),
                            energies.begin() + static_cast<std::ptrdiff_t>(nUse));

  result.nPeaksFitted = nUse;
  result.fitted = true;
  result.calibrated =
      FitLinearAdcVsEnergy(result.meansAdc, result.energiesMeV, result.calibA,
                           result.calibB);
  return fit;
}

static double EstimateCoincidenceWindowFromFit(TH1D *hist,
                                               TF1 *&fitOut,
                                               double nSigma = 5.0,
                                               ULong64_t probeWindow = 400000) {
  fitOut = nullptr;
  if (hist == nullptr || hist->GetEntries() < 10) {
    std::cout << "Not enough timing entries to fit. Falling back to 100000.\n";
    return 100000.0;
  }

  fitOut = FitMainPeak(hist, "f_dt_max_span_gaus");
  if (fitOut == nullptr) {
    std::cout << "Timing fit failed. Falling back to 100000.\n";
    return 100000.0;
  }

  const double mean = std::max(0.0, fitOut->GetParameter(1));
  const double sigma = std::abs(fitOut->GetParameter(2));
  if (nSigma <= 0.0) {
    nSigma = 5.0;
  }

  double window = mean + nSigma * sigma;

  if (window < 1000.0) {
    window = 1000.0;
  }
  if (window > static_cast<double>(probeWindow)) {
    window = static_cast<double>(probeWindow);
  }

  std::cout << "Timing fit (max |dt| of ohmic+pin0+pin1):\n"
            << "  mean  = " << mean << "\n"
            << "  sigma = " << sigma << "\n"
            << "  nSigma = " << nSigma << "\n"
            << "  window = mean + nSigma*sigma = " << window << "\n";
  return window;
}

static std::vector<CoincidenceEvent> FindCoincidences(
    const std::vector<Event> &events,
    const std::vector<DetectorMap> &detectorMaps,
    const std::map<long long, ChannelInfo> &channelLookup,
    double coincidenceWindow,
    double adcThreshold = 50.0) {
  std::vector<CoincidenceEvent> out;
  if (events.empty()) {
    return out;
  }

  const ULong64_t window = static_cast<ULong64_t>(std::ceil(coincidenceWindow));
  std::size_t left = 0;

  // Ohmic is only the trigger label, not required to arrive first.
  // For each ohmic, search BOTH earlier and later hits in [t-window, t+window].
  for (std::size_t i = 0; i < events.size(); ++i) {
    const Event &ev = events[i];
    if (!IsValidPulse(ev, adcThreshold)) {
      continue;
    }

    const ChannelInfo *ohmicInfo =
        FindChannelInfo(channelLookup, ev.board, ev.channel);
    if (ohmicInfo == nullptr || !ohmicInfo->isOhmic) {
      continue;
    }

    const DetectorMap *currentMap =
        FindDetectorMap(detectorMaps, ohmicInfo->detector);
    if (currentMap == nullptr || currentMap->board != ev.board) {
      continue;
    }

    while (left < i && ev.timestamp > events[left].timestamp &&
           ev.timestamp - events[left].timestamp > window) {
      ++left;
    }

    std::size_t right = i;
    while (right + 1 < events.size() &&
           events[right + 1].timestamp >= ev.timestamp &&
           events[right + 1].timestamp - ev.timestamp <= window) {
      ++right;
    }

    CoincidenceEvent best;
    bool found = false;
    double bestSpan = std::numeric_limits<double>::max();

    for (const auto &pair : currentMap->pairs) {
      // pin0/pin1 are channel identities from mapping, not arrival order.
      const Event *hitPin0 = nullptr;
      const Event *hitPin1 = nullptr;

      for (std::size_t j = left; j <= right; ++j) {
        if (j == i) {
          continue;
        }

        const Event &cand = events[j];
        if (cand.board != ev.board || !IsValidPulse(cand, adcThreshold)) {
          continue;
        }

        const ChannelInfo *candInfo =
            FindChannelInfo(channelLookup, cand.board, cand.channel);
        if (candInfo == nullptr || candInfo->detector != ohmicInfo->detector ||
            candInfo->isOhmic) {
          continue;
        }

        if (cand.channel == pair.chA) {
          if (hitPin0 == nullptr ||
              std::abs(static_cast<double>(ev.timestamp) -
                       static_cast<double>(cand.timestamp)) <
                  std::abs(static_cast<double>(ev.timestamp) -
                           static_cast<double>(hitPin0->timestamp))) {
            hitPin0 = &cand;
          }
        } else if (cand.channel == pair.chB) {
          if (hitPin1 == nullptr ||
              std::abs(static_cast<double>(ev.timestamp) -
                       static_cast<double>(cand.timestamp)) <
                  std::abs(static_cast<double>(ev.timestamp) -
                           static_cast<double>(hitPin1->timestamp))) {
            hitPin1 = &cand;
          }
        }
      }

      if (hitPin0 == nullptr || hitPin1 == nullptr) {
        continue;
      }

      const double dtA =
          std::abs(static_cast<double>(ev.timestamp) -
                   static_cast<double>(hitPin0->timestamp));
      const double dtB =
          std::abs(static_cast<double>(ev.timestamp) -
                   static_cast<double>(hitPin1->timestamp));
      const double dtAB =
          std::abs(static_cast<double>(hitPin0->timestamp) -
                   static_cast<double>(hitPin1->timestamp));
      const double maxSpan = std::max(dtAB, std::max(dtA, dtB));

      if (dtA <= coincidenceWindow && dtB <= coincidenceWindow &&
          dtAB <= coincidenceWindow && maxSpan < bestSpan) {
        best.board = currentMap->board;
        best.detector = currentMap->detector;
        best.stripIndex = pair.stripIndex;
        best.ohmicIndex = ohmicInfo->ohmicIndex;
        best.energyPin0 = hitPin0->energy;
        best.energyPin1 = hitPin1->energy;
        best.ohmicEnergy = ev.energy;
        best.junctionSum = static_cast<double>(hitPin0->energy) +
                           static_cast<double>(hitPin1->energy);
        if (best.junctionSum > 0.0) {
          best.junctionPos =
              (best.energyPin0 - best.energyPin1) / best.junctionSum;
        } else {
          best.junctionPos = 0.0;
        }
        best.anchorTs = ev.timestamp;
        best.dtOhmicPin0 = dtA;
        best.dtOhmicPin1 = dtB;
        best.dtPin0Pin1 = dtAB;
        best.maxSpan = maxSpan;
        bestSpan = maxSpan;
        found = true;
      }
    }

    if (found) {
      out.push_back(best);
    }
  }

  std::cout << "Found " << out.size()
            << " 2-junction + 1-ohmic coincidence candidates "
            << "with window=" << coincidenceWindow << "\n";
  return out;
}

static TH1D *MakeHistogram(const char *name,
                           const char *title,
                           const std::vector<double> &values,
                           int nBins = 400) {
  double minVal = std::numeric_limits<double>::max();
  double maxVal = std::numeric_limits<double>::lowest();

  for (double v : values) {
    if (v < minVal) minVal = v;
    if (v > maxVal) maxVal = v;
  }

  if (!std::isfinite(minVal) || !std::isfinite(maxVal) || values.empty()) {
    minVal = 0.0;
    maxVal = 1.0;
  }

  if (maxVal <= minVal) {
    maxVal = minVal + 1.0;
  }

  const double padding = 0.05 * (maxVal - minVal);
  TH1D *hist = new TH1D(name, title, nBins, minVal - padding, maxVal + padding);
  for (double v : values) {
    hist->Fill(v);
  }
  return hist;
}

static TH1D *MakeTimingHistogram(const char *name,
                                 const char *title,
                                 const std::vector<double> &values,
                                 double xMax = 1.0e7,
                                 int nBins = 200) {
  if (xMax <= 0.0) {
    xMax = 1.0e7;
  }

  TH1D *hist = new TH1D(name, title, nBins, 0.0, xMax);
  for (double v : values) {
    hist->Fill(v);
  }
  return hist;
}

static void ChoosePadGrid(int nPads, int &nCols, int &nRows) {
  if (nPads <= 1) {
    nCols = 1;
    nRows = 1;
  } else if (nPads <= 4) {
    nCols = 2;
    nRows = 2;
  } else if (nPads <= 6) {
    nCols = 3;
    nRows = 2;
  } else if (nPads <= 8) {
    nCols = 4;
    nRows = 2;
  } else if (nPads <= 9) {
    nCols = 3;
    nRows = 3;
  } else if (nPads <= 12) {
    nCols = 4;
    nRows = 3;
  } else if (nPads <= 16) {
    nCols = 4;
    nRows = 4;
  } else {
    nCols = 4;
    nRows = (nPads + nCols - 1) / nCols;
  }
}

// One histogram pad for DrawDetectorPagesPdf (1D with optional fit, or 2D COLZ).
struct DetectorPad {
  TH1 *hist = nullptr;
  TF1 *fit = nullptr;
  int color = kBlue + 1;
  bool is2D = false;
};

// Flexible multi-page PDF: one detector (one vector of pads) per page.
// Grid is chosen automatically from the number of pads.
static void DrawDetectorPagesPdf(const std::vector<std::vector<DetectorPad>> &detPages,
                                 const TString &outDir,
                                 const TString &fileName) {
  if (detPages.empty()) {
    return;
  }

  std::size_t maxPads = 0;
  for (const auto &pads : detPages) {
    maxPads = std::max(maxPads, pads.size());
  }
  if (maxPads == 0) {
    return;
  }

  int nCols = 1;
  int nRows = 1;
  ChoosePadGrid(static_cast<int>(maxPads), nCols, nRows);

  const TString pdfPath = outDir + "/" + fileName;
  TCanvas *c = new TCanvas(Form("c_%s", fileName.Data()), fileName, 1600, 1200);
  c->Print(pdfPath + "[", "pdf");

  for (std::size_t idet = 0; idet < detPages.size(); ++idet) {
    const auto &pads = detPages[idet];
    c->Clear();
    c->Divide(nCols, nRows);

    for (std::size_t ip = 0; ip < pads.size(); ++ip) {
      c->cd(static_cast<int>(ip) + 1);
      gPad->SetLeftMargin(pads[ip].is2D ? 0.14 : 0.15);
      gPad->SetRightMargin(pads[ip].is2D ? 0.14 : 0.05);
      gPad->SetBottomMargin(0.14);
      gPad->SetTopMargin(0.10);

      TH1 *hist = pads[ip].hist;
      if (hist == nullptr) {
        continue;
      }

      if (pads[ip].is2D) {
        hist->SetStats(0);
        hist->GetXaxis()->SetTitleSize(0.05);
        hist->GetYaxis()->SetTitleSize(0.05);
        hist->GetXaxis()->SetLabelSize(0.045);
        hist->GetYaxis()->SetLabelSize(0.045);
        hist->Draw("COLZ");
      } else {
        hist->SetLineColor(pads[ip].color);
        hist->SetLineWidth(2);
        hist->SetTitleSize(0.06, "");
        hist->GetXaxis()->SetLabelSize(0.05);
        hist->GetYaxis()->SetLabelSize(0.05);
        hist->GetXaxis()->SetTitleSize(0.05);
        hist->GetYaxis()->SetTitleSize(0.05);
        hist->Draw();
        if (pads[ip].fit != nullptr) {
          pads[ip].fit->SetLineColor(kRed + 1);
          pads[ip].fit->SetLineWidth(2);
          pads[ip].fit->Draw("same");
        }
      }
    }
    c->Print(pdfPath, "pdf");
  }

  c->Print(pdfPath + "]", "pdf");
  std::cout << "Wrote PDF (one detector/page, " << nCols << "x" << nRows
            << " pads): " << pdfPath << "\n";
}

static void DrawAndSaveTimingPdf(const std::vector<TH1D *> &hists,
                                 const std::vector<int> &colors,
                                 const TString &outDir,
                                 const TString &fileName = "dt_timing.pdf",
                                 TF1 *fit = nullptr,
                                 int fitHistIndex = -1) {
  if (hists.empty()) {
    return;
  }

  const TString pdfPath = outDir + "/" + fileName;
  TCanvas *c = new TCanvas("c_dt_timing", "Timing diagnostics", 900, 650);

  // Timing diagnostics: one histogram per page.
  c->Print(pdfPath + "[", "pdf");
  for (std::size_t i = 0; i < hists.size(); ++i) {
    c->cd();
    c->Clear();
    gPad->SetLogy(0);
    hists[i]->SetLineColor(colors[i]);
    hists[i]->SetLineWidth(2);
    hists[i]->Draw();
    if (fit != nullptr && static_cast<int>(i) == fitHistIndex) {
      fit->SetLineColor(kRed + 1);
      fit->SetLineWidth(2);
      fit->Draw("same");
    }
    c->Print(pdfPath, "pdf");
  }
  c->Print(pdfPath + "]", "pdf");

  std::cout << "Wrote timing PDF (one histogram per page): " << pdfPath << "\n";
}

static void CreateEnergyHistograms(TH1D *hOhmic[kNumDetectors][kNumOhmic],
                                   TH1D *hJunction[kNumDetectors][kNumStrips],
                                   TH2D *hPin0Pin1[kNumDetectors][kNumStrips],
                                   TH2D *hPosEnergy[kNumDetectors][kNumStrips]) {
  for (int det = 0; det < kNumDetectors; ++det) {
    for (int ohm = 0; ohm < kNumOhmic; ++ohm) {
      hOhmic[det][ohm] = new TH1D(
          Form("h_ohmic_d%d_o%d", det, ohm),
          Form("Ohmic energy det=%d ohmic=%d;ADC;Counts", det, ohm),
          kEnergyBins, kEnergyMin, kEnergyMax);
    }
    for (int strip = 0; strip < kNumStrips; ++strip) {
      hJunction[det][strip] = new TH1D(
          Form("h_junction_d%d_s%d", det, strip),
          Form("Junction sum det=%d strip=%d;ADC sum;Counts", det, strip),
          kEnergyBins, kEnergyMin, kEnergyMax);
      hPin0Pin1[det][strip] = new TH2D(
          Form("h_pin0_pin1_d%d_s%d", det, strip),
          Form("E_{pin0} vs E_{pin1} det=%d strip=%d;E_{pin0} [ADC];E_{pin1} [ADC]",
               det, strip),
          kPinShareBins, kEnergyMin, kEnergyMax, kPinShareBins, kEnergyMin,
          kEnergyMax);
      hPosEnergy[det][strip] = new TH2D(
          Form("h_pos_energy_d%d_s%d", det, strip),
          Form("Position vs energy det=%d strip=%d;"
               "(E_{0}-E_{1})/(E_{0}+E_{1});E_{0}+E_{1} [ADC]",
               det, strip),
          kPositionBins, kPositionMin, kPositionMax, kEnergyBins, kEnergyMin,
          kEnergyMax);
    }
  }
}

static void FillEnergyHistograms(
    const std::vector<CoincidenceEvent> &coincidences,
    TH1D *hOhmic[kNumDetectors][kNumOhmic],
    TH1D *hJunction[kNumDetectors][kNumStrips],
    TH2D *hPin0Pin1[kNumDetectors][kNumStrips],
    TH2D *hPosEnergy[kNumDetectors][kNumStrips]) {
  for (const auto &co : coincidences) {
    if (co.detector < 0 || co.detector >= kNumDetectors) {
      continue;
    }
    if (co.ohmicIndex >= 0 && co.ohmicIndex < kNumOhmic) {
      hOhmic[co.detector][co.ohmicIndex]->Fill(co.ohmicEnergy);
    }
    if (co.stripIndex >= 0 && co.stripIndex < kNumStrips) {
      hJunction[co.detector][co.stripIndex]->Fill(co.junctionSum);
      hPin0Pin1[co.detector][co.stripIndex]->Fill(co.energyPin0, co.energyPin1);
      hPosEnergy[co.detector][co.stripIndex]->Fill(co.junctionPos, co.junctionSum);
    }
  }
}

static void FitEnergyHistograms(
    TH1D *hOhmic[kNumDetectors][kNumOhmic],
    TH1D *hJunction[kNumDetectors][kNumStrips],
    TF1 *fOhmic[kNumDetectors][kNumOhmic],
    TF1 *fJunction[kNumDetectors][kNumStrips],
    int nPeaks,
    const std::vector<double> &peakEnergiesMeV,
    std::vector<PeakFitResult> &ohmicFits,
    std::vector<PeakFitResult> &junctionFits) {
  ohmicFits.clear();
  junctionFits.clear();

  for (int det = 0; det < kNumDetectors; ++det) {
    for (int ohm = 0; ohm < kNumOhmic; ++ohm) {
      PeakFitResult o;
      o.detector = det;
      o.index = ohm;
      o.entries = hOhmic[det][ohm]->GetEntries();
      fOhmic[det][ohm] = FitCalibrationPeaks(
          hOhmic[det][ohm], Form("f_ohmic_d%d_o%d", det, ohm), nPeaks,
          peakEnergiesMeV, o);
      ohmicFits.push_back(o);
    }

    for (int strip = 0; strip < kNumStrips; ++strip) {
      PeakFitResult j;
      j.detector = det;
      j.index = strip;
      j.entries = hJunction[det][strip]->GetEntries();
      fJunction[det][strip] = FitCalibrationPeaks(
          hJunction[det][strip], Form("f_junction_d%d_s%d", det, strip), nPeaks,
          peakEnergiesMeV, j);
      junctionFits.push_back(j);
    }
  }
}

static void WritePeakFitSummaryBlock(FILE *fp,
                                     const char *title,
                                     const char *indexName,
                                     const std::vector<PeakFitResult> &fits) {
  fprintf(fp, "\n%s\n", title);
  fprintf(fp, "  det %s entries nFit\n", indexName);
  for (const auto &r : fits) {
    fprintf(fp, "  %3d %5d %7.0f %4d", r.detector, r.index, r.entries,
            r.nPeaksFitted);
    if (!r.fitted || r.meansAdc.empty()) {
      fprintf(fp, "  (no peaks)\n");
      continue;
    }
    fprintf(fp, "\n");
    for (std::size_t ip = 0; ip < r.meansAdc.size(); ++ip) {
      const double e = (ip < r.energiesMeV.size()) ? r.energiesMeV[ip] : 0.0;
      const double sig = (ip < r.sigmasAdc.size()) ? r.sigmasAdc[ip] : 0.0;
      fprintf(fp, "           peak %zu: E=%.6f MeV  mean=%.3f ADC  sigma=%.3f\n",
              ip, e, r.meansAdc[ip], sig);
    }
  }
}

static void SaveCalibrationCoeffs(const std::vector<PeakFitResult> &ohmicFits,
                                  const std::vector<PeakFitResult> &junctionFits,
                                  int nPeaks,
                                  const std::vector<double> &peakEnergiesMeV,
                                  const TString &outDir) {
  const TString txtPath = outDir + "/x6_calibration_coeffs.txt";
  FILE *fp = fopen(txtPath.Data(), "w");
  if (fp == nullptr) {
    std::cerr << "Could not write calibration coeffs file: " << txtPath << "\n";
    return;
  }

  fprintf(fp, "X6 calibration coefficients\n");
  fprintf(fp, "===========================\n");
  fprintf(fp, "Form: ADC = A + B * E[MeV]\n");
  fprintf(fp, "      E[MeV] = (ADC - A) / B\n");
  fprintf(fp, "Number of peaks: %d\n", nPeaks);
  fprintf(fp, "Peak energies [MeV]:");
  for (double e : peakEnergiesMeV) {
    fprintf(fp, " %.6f", e);
  }
  fprintf(fp, "\n\n");

  fprintf(fp, "Ohmic (det, ohmic 0-3)\n");
  fprintf(fp, "  det ohmic calibrated            A            B\n");
  for (const auto &r : ohmicFits) {
    fprintf(fp, "  %3d %5d %10d %12.6f %12.6f\n", r.detector, r.index,
            r.calibrated ? 1 : 0, r.calibA, r.calibB);
  }

  fprintf(fp, "\nJunction sum (det, strip 0-7)\n");
  fprintf(fp, "  det strip calibrated            A            B\n");
  for (const auto &r : junctionFits) {
    fprintf(fp, "  %3d %5d %10d %12.6f %12.6f\n", r.detector, r.index,
            r.calibrated ? 1 : 0, r.calibA, r.calibB);
  }

  fclose(fp);
  std::cout << "Wrote calibration coeffs: " << txtPath << "\n";
}

static const PeakFitResult *FindPeakFitResult(const std::vector<PeakFitResult> &fits,
                                              int detector,
                                              int index) {
  for (const auto &r : fits) {
    if (r.detector == detector && r.index == index) {
      return &r;
    }
  }
  return nullptr;
}

// Draw ADC vs E linear calibration: one detector per page.
// Pad order: ohmic 0-3, then junction strips 0-7.
static void DrawCalibrationLinearPdf(const std::vector<PeakFitResult> &ohmicFits,
                                     const std::vector<PeakFitResult> &junctionFits,
                                     const TString &outDir) {
  const int nPads = kNumOhmic + kNumStrips;
  int nCols = 1;
  int nRows = 1;
  ChoosePadGrid(nPads, nCols, nRows);

  const TString pdfPath = outDir + "/calibration_linear_fits.pdf";
  TCanvas *c = new TCanvas("c_calib_linear", "calibration linear fits", 1600, 1200);
  c->Print(pdfPath + "[", "pdf");

  for (int det = 0; det < kNumDetectors; ++det) {
    c->Clear();
    c->Divide(nCols, nRows);

    for (int pad = 0; pad < nPads; ++pad) {
      c->cd(pad + 1);
      gPad->SetLeftMargin(0.16);
      gPad->SetRightMargin(0.04);
      gPad->SetBottomMargin(0.14);
      gPad->SetTopMargin(0.12);

      const bool isOhmic = (pad < kNumOhmic);
      const int index = isOhmic ? pad : (pad - kNumOhmic);
      const PeakFitResult *r =
          isOhmic ? FindPeakFitResult(ohmicFits, det, index)
                  : FindPeakFitResult(junctionFits, det, index);

      const TString title =
          isOhmic ? Form("Det %d Ohmic %d", det, index)
                  : Form("Det %d Junction strip %d", det, index);

      if (r == nullptr || !r->calibrated || r->meansAdc.empty()) {
        TH1D frame(Form("frame_empty_d%d_p%d", det, pad), title + ";E [MeV];ADC",
                   10, 0.0, 1.0);
        frame.SetStats(0);
        frame.SetMinimum(0.0);
        frame.SetMaximum(1.0);
        frame.DrawCopy();
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.07);
        latex.DrawLatex(0.25, 0.5, "no calibration");
        continue;
      }

      double eMin = r->energiesMeV.front();
      double eMax = r->energiesMeV.front();
      double adcMin = r->meansAdc.front();
      double adcMax = r->meansAdc.front();
      for (std::size_t i = 0; i < r->meansAdc.size(); ++i) {
        const double e =
            (i < r->energiesMeV.size()) ? r->energiesMeV[i] : 0.0;
        eMin = std::min(eMin, e);
        eMax = std::max(eMax, e);
        adcMin = std::min(adcMin, r->meansAdc[i]);
        adcMax = std::max(adcMax, r->meansAdc[i]);
      }

      const double ePad = std::max(0.15 * std::max(eMax - eMin, 0.2), 0.05);
      const double xLo = std::max(0.0, eMin - ePad);
      const double xHi = eMax + ePad;
      const double yLineLo = r->calibA + r->calibB * xLo;
      const double yLineHi = r->calibA + r->calibB * xHi;
      adcMin = std::min(adcMin, std::min(yLineLo, yLineHi));
      adcMax = std::max(adcMax, std::max(yLineLo, yLineHi));
      const double yPad = std::max(0.15 * (adcMax - adcMin), 5.0);
      const double yLo = std::max(0.0, adcMin - yPad);
      const double yHi = adcMax + yPad;

      TH1D frame(Form("frame_d%d_p%d", det, pad), title + ";E [MeV];ADC", 10, xLo,
                 xHi);
      frame.SetStats(0);
      frame.SetMinimum(yLo);
      frame.SetMaximum(yHi);
      frame.GetXaxis()->SetTitleSize(0.05);
      frame.GetYaxis()->SetTitleSize(0.05);
      frame.GetXaxis()->SetLabelSize(0.045);
      frame.GetYaxis()->SetLabelSize(0.045);
      frame.DrawCopy();

      TGraphErrors *gr = new TGraphErrors(static_cast<int>(r->meansAdc.size()));
      gr->SetName(Form("gr_calib_d%d_%s%d", det, isOhmic ? "o" : "s", index));
      for (std::size_t i = 0; i < r->meansAdc.size(); ++i) {
        const double e =
            (i < r->energiesMeV.size()) ? r->energiesMeV[i] : 0.0;
        const double sig =
            (i < r->sigmasAdc.size()) ? r->sigmasAdc[i] : 0.0;
        gr->SetPoint(static_cast<int>(i), e, r->meansAdc[i]);
        gr->SetPointError(static_cast<int>(i), 0.0, sig);
      }
      gr->SetMarkerStyle(20);
      gr->SetMarkerSize(1.0);
      gr->SetMarkerColor(isOhmic ? kBlue + 1 : kGreen + 2);
      gr->SetLineColor(isOhmic ? kBlue + 1 : kGreen + 2);
      gr->Draw("P SAME");

      TF1 *line = new TF1(Form("f_line_d%d_%s%d", det, isOhmic ? "o" : "s", index),
                          "pol1", xLo, xHi);
      line->SetParameters(r->calibA, r->calibB);
      line->SetLineColor(kRed + 1);
      line->SetLineWidth(2);
      line->Draw("SAME");

      TLatex latex;
      latex.SetNDC();
      latex.SetTextSize(0.055);
      latex.DrawLatex(0.20, 0.82, Form("A=%.3f", r->calibA));
      latex.DrawLatex(0.20, 0.74, Form("B=%.4f", r->calibB));
    }

    c->Print(pdfPath, "pdf");
  }

  c->Print(pdfPath + "]", "pdf");
  std::cout << "Wrote PDF (one detector/page, " << nCols << "x" << nRows
            << " pads): " << pdfPath << "\n";
}

static void SaveSummary(const std::vector<CoincidenceEvent> &coincidences,
                        double coincidenceWindow,
                        double nSigma,
                        double adcThreshold,
                        int nPeaks,
                        const std::vector<double> &peakEnergiesMeV,
                        TF1 *fDt,
                        const TimingDiagnostics &timing,
                        const std::vector<PeakFitResult> &ohmicFits,
                        const std::vector<PeakFitResult> &junctionFits,
                        const TString &outDir) {
  const TString txtPath = outDir + "/x6_calibration_summary.txt";
  FILE *fp = fopen(txtPath.Data(), "w");
  if (fp == nullptr) {
    std::cerr << "Could not write summary file: " << txtPath << "\n";
    return;
  }

  fprintf(fp, "X6 calibration summary\n");
  fprintf(fp, "======================\n");
  fprintf(fp, "ADC threshold: %.3f\n", adcThreshold);
  fprintf(fp, "Number of calibration peaks: %d\n", nPeaks);
  fprintf(fp, "Peak energies [MeV] (assigned low→high ADC):");
  for (double e : peakEnergiesMeV) {
    fprintf(fp, " %.6f", e);
  }
  fprintf(fp, "\n");
  fprintf(fp, "Coincidence window used (timestamp units): %.3f\n", coincidenceWindow);
  fprintf(fp, "  Applied equally to |t_ohmic-t_pin0|, |t_ohmic-t_pin1|, |t_pin0-t_pin1|\n");
  if (fDt != nullptr) {
    fprintf(fp, "Timing fit (max |dt| of ohmic+pin0+pin1)\n");
    fprintf(fp, "  mean   = %.6f\n", fDt->GetParameter(1));
    fprintf(fp, "  sigma  = %.6f\n", fDt->GetParameter(2));
    fprintf(fp, "  nSigma = %.3f\n", nSigma);
    fprintf(fp, "  window = mean + nSigma*sigma = %.6f\n", coincidenceWindow);
  }
  fprintf(fp, "Coincidence count: %zu\n", coincidences.size());
  fprintf(fp, "\nTiming diagnostics sample counts\n");
  fprintf(fp, "  |t_pin0 - t_pin1| : %zu\n", timing.dtJunctionPair.size());
  fprintf(fp, "  |t_ohmic - t_pin0|: %zu\n", timing.dtOhmicPin0.size());
  fprintf(fp, "  |t_ohmic - t_pin1|: %zu\n", timing.dtOhmicPin1.size());
  fprintf(fp, "  |t_ohmic - t_pin| : %zu\n", timing.dtOhmicJun.size());
  fprintf(fp, "  pin0-pin1 w/ ohmic: %zu\n", timing.dtPin0Pin1.size());
  fprintf(fp, "  max span          : %zu\n", timing.maxSpan.size());

  WritePeakFitSummaryBlock(fp, "Ohmic peak fits (per det, ohmic 0-3)", "ohmic",
                           ohmicFits);
  WritePeakFitSummaryBlock(fp, "Junction-sum peak fits (per det, strip 0-7)",
                           "strip", junctionFits);

  fclose(fp);
  std::cout << "Wrote summary: " << txtPath << "\n";

  SaveCalibrationCoeffs(ohmicFits, junctionFits, nPeaks, peakEnergiesMeV, outDir);
  DrawCalibrationLinearPdf(ohmicFits, junctionFits, outDir);
}


// =============================================================================
// Config I/O
// =============================================================================

static std::string Trim(const std::string &s) {
  std::size_t b = 0;
  while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b]))) {
    ++b;
  }
  std::size_t e = s.size();
  while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1]))) {
    --e;
  }
  return s.substr(b, e - b);
}

static std::string StripComment(const std::string &s) {
  bool inQuote = false;
  for (std::size_t i = 0; i < s.size(); ++i) {
    if (s[i] == '"') {
      inQuote = !inQuote;
    }
    if (!inQuote && s[i] == '#') {
      return Trim(s.substr(0, i));
    }
  }
  return Trim(s);
}

static bool ParseBoolToken(const std::string &v, bool &out) {
  const std::string t = Trim(v);
  if (t == "1" || t == "true" || t == "TRUE" || t == "yes" || t == "YES" ||
      t == "on" || t == "ON") {
    out = true;
    return true;
  }
  if (t == "0" || t == "false" || t == "FALSE" || t == "no" || t == "NO" ||
      t == "off" || t == "OFF") {
    out = false;
    return true;
  }
  return false;
}

static std::vector<double> ParseDoubleList(const std::string &v) {
  std::vector<double> out;
  std::string token;
  for (std::size_t i = 0; i <= v.size(); ++i) {
    const char c = (i < v.size()) ? v[i] : ',';
    if (c == ',' || i == v.size()) {
      const std::string t = Trim(token);
      if (!t.empty()) {
        out.push_back(std::atof(t.c_str()));
      }
      token.clear();
    } else {
      token.push_back(c);
    }
  }
  return out;
}

static std::string Unquote(const std::string &v) {
  std::string t = Trim(v);
  if (t.size() >= 2 && ((t.front() == '"' && t.back() == '"') ||
                        (t.front() == '\'' && t.back() == '\''))) {
    return t.substr(1, t.size() - 2);
  }
  return t;
}

bool LoadConfig(const std::string &path, RunConfig &cfg) {
  FILE *fp = fopen(path.c_str(), "r");
  if (fp == nullptr) {
    std::cerr << "Could not open config file: " << path << "\n";
    return false;
  }

  char buffer[1024];
  int lineNumber = 0;
  while (fgets(buffer, sizeof(buffer), fp) != nullptr) {
    ++lineNumber;
    std::string line = StripComment(buffer);
    if (line.empty()) {
      continue;
    }

    const std::size_t eq = line.find('=');
    if (eq == std::string::npos) {
      std::cerr << "Config line " << lineNumber << ": missing '=' in: " << line
                << "\n";
      fclose(fp);
      return false;
    }

    const std::string key = Trim(line.substr(0, eq));
    const std::string value = Unquote(Trim(line.substr(eq + 1)));

    if (key == "input_pattern") {
      cfg.inputPattern = value;
    } else if (key == "output_dir") {
      cfg.outputDir = value;
    } else if (key == "mapping_file") {
      cfg.mappingFile = value;
    } else if (key == "use_energy_short") {
      if (!ParseBoolToken(value, cfg.useEnergyShort)) {
        std::cerr << "Config line " << lineNumber
                  << ": bad bool for use_energy_short\n";
        fclose(fp);
        return false;
      }
    } else if (key == "coincidence_window") {
      cfg.coincidenceWindow = std::atof(value.c_str());
    } else if (key == "n_sigma") {
      cfg.nSigma = std::atof(value.c_str());
    } else if (key == "adc_threshold") {
      cfg.adcThreshold = std::atof(value.c_str());
    } else if (key == "n_peaks") {
      cfg.nPeaks = std::atoi(value.c_str());
    } else if (key == "peak_energies_mev") {
      cfg.peakEnergiesMeV = ParseDoubleList(value);
    } else if (key == "timing_probe_window") {
      cfg.timingProbeWindow = static_cast<ULong64_t>(std::atoll(value.c_str()));
    } else {
      std::cerr << "Config line " << lineNumber << ": unknown key '" << key
                << "'\n";
      fclose(fp);
      return false;
    }
  }
  fclose(fp);

  if (cfg.nPeaks < 1) {
    cfg.nPeaks = 1;
  }
  if (cfg.nPeaks > 10) {
    cfg.nPeaks = 10;
  }
  if (cfg.peakEnergiesMeV.empty()) {
    cfg.peakEnergiesMeV.push_back(5.486);
  }
  if (static_cast<int>(cfg.peakEnergiesMeV.size()) < cfg.nPeaks) {
    cfg.peakEnergiesMeV.resize(static_cast<std::size_t>(cfg.nPeaks),
                               cfg.peakEnergiesMeV.back());
  } else if (static_cast<int>(cfg.peakEnergiesMeV.size()) > cfg.nPeaks) {
    cfg.peakEnergiesMeV.resize(static_cast<std::size_t>(cfg.nPeaks));
  }
  std::sort(cfg.peakEnergiesMeV.begin(), cfg.peakEnergiesMeV.end());
  return true;
}

void PrintConfig(const RunConfig &cfg) {
  std::cout << "\nX6 calibration config\n";
  std::cout << "=====================\n";
  std::cout << "  input_pattern        : " << cfg.inputPattern << "\n";
  std::cout << "  output_dir           : " << cfg.outputDir << "\n";
  std::cout << "  mapping_file         : " << cfg.mappingFile << "\n";
  std::cout << "  use_energy_short     : " << (cfg.useEnergyShort ? "yes" : "no")
            << "\n";
  std::cout << "  coincidence_window   : " << cfg.coincidenceWindow
            << "  (-1 = auto)\n";
  std::cout << "  n_sigma              : " << cfg.nSigma << "\n";
  std::cout << "  adc_threshold        : " << cfg.adcThreshold << "\n";
  std::cout << "  n_peaks              : " << cfg.nPeaks << "\n";
  std::cout << "  peak_energies_mev    :";
  for (double e : cfg.peakEnergiesMeV) {
    std::cout << " " << e;
  }
  std::cout << "\n";
  std::cout << "  timing_probe_window  : " << cfg.timingProbeWindow << "\n\n";
}

// Resolve relative paths with string-based '..' using $PWD (logical path).
// Needed when the package is a symlink inside a CoMPASS project: Linux follows
// the real directory for '..', so "../DAQ" would miss the project's DAQ/.
static std::string ResolvePath(const std::string &path) {
  if (path.empty() || path[0] == '/') {
    return path;
  }

  const char *pwdEnv = std::getenv("PWD");
  std::string base = (pwdEnv != nullptr) ? std::string(pwdEnv) : std::string();
  if (base.empty()) {
    return path;
  }

  std::string joined = base + "/" + path;
  std::vector<std::string> parts;
  std::string token;
  for (std::size_t i = 0; i <= joined.size(); ++i) {
    const char c = (i < joined.size()) ? joined[i] : '/';
    if (c == '/') {
      if (token.empty() || token == ".") {
        // skip
      } else if (token == "..") {
        if (!parts.empty()) {
          parts.pop_back();
        }
      } else {
        parts.push_back(token);
      }
      token.clear();
    } else {
      token.push_back(c);
    }
  }

  std::string out = "/";
  for (std::size_t i = 0; i < parts.size(); ++i) {
    if (i > 0) {
      out += "/";
    }
    out += parts[i];
  }
  return out;
}

void RunCalibration(const RunConfig &cfg) {
  const std::string resolvedInput = ResolvePath(cfg.inputPattern);
  const std::string resolvedOutput = ResolvePath(cfg.outputDir);
  const std::string resolvedMapping = ResolvePath(cfg.mappingFile);

  if (resolvedInput != cfg.inputPattern || resolvedMapping != cfg.mappingFile ||
      resolvedOutput != cfg.outputDir) {
    std::cout << "Resolved paths (via $PWD):\n"
              << "  input_pattern : " << resolvedInput << "\n"
              << "  output_dir    : " << resolvedOutput << "\n"
              << "  mapping_file  : " << resolvedMapping << "\n";
  }

  const TString tInputPattern(resolvedInput.c_str());
  const TString tOutputDir(resolvedOutput.c_str());
  const TString tMappingFile(resolvedMapping.c_str());

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(1111);

  gSystem->mkdir(tOutputDir, kTRUE);

  TChain chain("Data_R");
  LoadChain(chain, tInputPattern);
  if (chain.GetEntries() <= 0) {
    return;
  }

  std::map<long long, ChannelInfo> channelLookup;
  const auto detectorMaps = LoadDetectorMapFromFile(tMappingFile, channelLookup);
  if (detectorMaps.empty()) {
    std::cerr << "No detector mapping loaded. Check mapping file: "
              << tMappingFile << "\n";
    return;
  }
  const auto events = ReadEvents(chain, cfg.useEnergyShort);

  TimingDiagnostics timing;
  CollectTimingDiagnostics(events, detectorMaps, channelLookup, timing,
                           cfg.adcThreshold, cfg.timingProbeWindow);

  const double timingXMax = 4.0e5;
  TH1D *hDtMaxSpan = MakeTimingHistogram(
      "h_dt_max_span",
      "max |#Delta t| of ohmic+pin0+pin1;max |#Delta t| [ps];Counts",
      timing.maxSpan, timingXMax);

  TF1 *fDtMaxSpan = nullptr;
  double coincidenceWindow = cfg.coincidenceWindow;
  if (coincidenceWindow <= 0.0) {
    coincidenceWindow = EstimateCoincidenceWindowFromFit(
        hDtMaxSpan, fDtMaxSpan, cfg.nSigma, cfg.timingProbeWindow);
  } else {
    fDtMaxSpan = FitMainPeak(hDtMaxSpan, "f_dt_max_span_gaus");
  }

  std::cout << "Actual coincidence window used: " << coincidenceWindow
            << " timestamp units\n"
            << "  (applied equally to |t_ohmic-t_pin0|, |t_ohmic-t_pin1|, "
            << "|t_pin0-t_pin1|)\n";

  const auto coincidences = FindCoincidences(
      events, detectorMaps, channelLookup, coincidenceWindow, cfg.adcThreshold);

  std::vector<double> acceptedMaxSpan;
  acceptedMaxSpan.reserve(coincidences.size());
  for (const auto &co : coincidences) {
    acceptedMaxSpan.push_back(co.maxSpan);
  }

  TH1D *hOhmic[kNumDetectors][kNumOhmic] = {};
  TH1D *hJunction[kNumDetectors][kNumStrips] = {};
  TH2D *hPin0Pin1[kNumDetectors][kNumStrips] = {};
  TH2D *hPosEnergy[kNumDetectors][kNumStrips] = {};
  TF1 *fOhmic[kNumDetectors][kNumOhmic] = {};
  TF1 *fJunction[kNumDetectors][kNumStrips] = {};
  CreateEnergyHistograms(hOhmic, hJunction, hPin0Pin1, hPosEnergy);
  FillEnergyHistograms(coincidences, hOhmic, hJunction, hPin0Pin1, hPosEnergy);

  std::vector<PeakFitResult> ohmicFits;
  std::vector<PeakFitResult> junctionFits;
  FitEnergyHistograms(hOhmic, hJunction, fOhmic, fJunction, cfg.nPeaks,
                      cfg.peakEnergiesMeV, ohmicFits, junctionFits);

  TH1D *hDtOhmicJun = MakeTimingHistogram(
      "h_dt_ohmic_jun", "|t_{ohmic} - t_{pin}|;|#Delta t| [ps];Counts",
      timing.dtOhmicJun, timingXMax);
  TH1D *hDtOhmicPin0 = MakeTimingHistogram(
      "h_dt_ohmic_pin0", "|t_{ohmic} - t_{pin0}|;|#Delta t| [ps];Counts",
      timing.dtOhmicPin0, timingXMax);
  TH1D *hDtOhmicPin1 = MakeTimingHistogram(
      "h_dt_ohmic_pin1", "|t_{ohmic} - t_{pin1}|;|#Delta t| [ps];Counts",
      timing.dtOhmicPin1, timingXMax);
  TH1D *hDtJunctionPair = MakeTimingHistogram(
      "h_dt_junction_pair", "|t_{pin0} - t_{pin1}|;|#Delta t| [ps];Counts",
      timing.dtJunctionPair, timingXMax);
  TH1D *hDtPin0Pin1 = MakeTimingHistogram(
      "h_dt_pin0_pin1_in_triple",
      "|t_{pin0} - t_{pin1}| (with ohmic);|#Delta t| [ps];Counts",
      timing.dtPin0Pin1, timingXMax);
  TH1D *hDtAcceptedSpan = MakeTimingHistogram(
      "h_dt_accepted_max_span",
      "max |#Delta t| of accepted coincidences;max |#Delta t| [ps];Counts",
      acceptedMaxSpan, timingXMax);

  DrawAndSaveTimingPdf(
      {hDtOhmicJun, hDtOhmicPin0, hDtOhmicPin1, hDtJunctionPair, hDtPin0Pin1,
       hDtMaxSpan, hDtAcceptedSpan},
      {kBlue + 1, kOrange + 7, kOrange + 1, kMagenta + 1, kCyan + 2, kViolet + 1,
       kPink + 2},
      tOutputDir, "dt_timing.pdf", fDtMaxSpan, 5);

  std::vector<std::vector<DetectorPad>> ohmicPages(kNumDetectors);
  std::vector<std::vector<DetectorPad>> junctionPages(kNumDetectors);
  std::vector<std::vector<DetectorPad>> pin0Pin1Pages(kNumDetectors);
  std::vector<std::vector<DetectorPad>> posEnergyPages(kNumDetectors);
  for (int det = 0; det < kNumDetectors; ++det) {
    for (int ohm = 0; ohm < kNumOhmic; ++ohm) {
      DetectorPad p;
      p.hist = hOhmic[det][ohm];
      p.fit = fOhmic[det][ohm];
      p.color = kBlue + 1;
      ohmicPages[det].push_back(p);
    }
    for (int strip = 0; strip < kNumStrips; ++strip) {
      DetectorPad j;
      j.hist = hJunction[det][strip];
      j.fit = fJunction[det][strip];
      j.color = kGreen + 2;
      junctionPages[det].push_back(j);

      DetectorPad p01;
      p01.hist = hPin0Pin1[det][strip];
      p01.is2D = true;
      pin0Pin1Pages[det].push_back(p01);

      DetectorPad pe;
      pe.hist = hPosEnergy[det][strip];
      pe.is2D = true;
      posEnergyPages[det].push_back(pe);
    }
  }
  DrawDetectorPagesPdf(ohmicPages, tOutputDir, "ohmic_peaks.pdf");
  DrawDetectorPagesPdf(junctionPages, tOutputDir, "junction_peaks.pdf");
  DrawDetectorPagesPdf(pin0Pin1Pages, tOutputDir, "junction_epin0_epin1.pdf");
  DrawDetectorPagesPdf(posEnergyPages, tOutputDir, "junction_position_energy.pdf");

  TFile outFile(tOutputDir + "/x6_calibration_histograms.root", "RECREATE");
  hDtOhmicJun->Write();
  hDtOhmicPin0->Write();
  hDtOhmicPin1->Write();
  hDtJunctionPair->Write();
  hDtPin0Pin1->Write();
  hDtMaxSpan->Write();
  hDtAcceptedSpan->Write();
  if (fDtMaxSpan != nullptr) {
    fDtMaxSpan->Write();
  }
  for (int det = 0; det < kNumDetectors; ++det) {
    for (int ohm = 0; ohm < kNumOhmic; ++ohm) {
      hOhmic[det][ohm]->Write();
      if (fOhmic[det][ohm] != nullptr) {
        fOhmic[det][ohm]->Write();
      }
    }
    for (int strip = 0; strip < kNumStrips; ++strip) {
      hJunction[det][strip]->Write();
      if (fJunction[det][strip] != nullptr) {
        fJunction[det][strip]->Write();
      }
      hPin0Pin1[det][strip]->Write();
      hPosEnergy[det][strip]->Write();
    }
  }
  outFile.Close();

  SaveSummary(coincidences, coincidenceWindow, cfg.nSigma, cfg.adcThreshold,
              cfg.nPeaks, cfg.peakEnergiesMeV, fDtMaxSpan, timing, ohmicFits,
              junctionFits, tOutputDir);

  std::cout << "\nDone. Results written to: " << tOutputDir << "\n";
}

}  // namespace X6Calib
