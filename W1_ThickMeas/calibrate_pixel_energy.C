#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TStyle.h>

#include "strip_map.h"

#include <array>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

struct PixelPeak {
  double peak = 0.0;
  double fwhm = 0.0;
  double resolution = 0.0;
  double entries = 0.0;
  bool valid = false;
};

struct PixelCalib {
  double slope = 0.0;
  double offset = 0.0;
  double adc1 = 0.0;
  double adc2 = 0.0;
  bool valid = false;
};

static bool LoadPeaksFile(const char *path, std::array<std::array<PixelPeak, kNY>, kNX> &peaks) {
  std::ifstream in(path);
  if (!in) {
    std::cerr << "Cannot open peaks file: " << path << std::endl;
    return false;
  }

  std::string line;
  int lineNumber = 0;
  while (std::getline(in, line)) {
    ++lineNumber;
    if (line.empty() || line[0] == '#') {
      continue;
    }

    std::istringstream iss(line);
    int x = -1;
    int y = -1;
    PixelPeak pixel;
    if (!(iss >> x >> y >> pixel.peak >> pixel.fwhm >> pixel.resolution)) {
      std::cerr << path << ":" << lineNumber << ": failed to parse line: " << line << std::endl;
      return false;
    }
    iss >> pixel.entries;  // optional

    if (x < 0 || x >= kNX || y < 0 || y >= kNY) {
      std::cerr << path << ":" << lineNumber << ": pixel index out of range: " << x << " " << y
                << std::endl;
      return false;
    }
    if (!(pixel.peak > 0.0)) {
      continue;
    }

    pixel.valid = true;
    peaks[x][y] = pixel;
  }

  return true;
}

static std::string ReplaceExtension(const std::string &path, const char *newExt) {
  const auto dot = path.find_last_of('.');
  if (dot == std::string::npos) {
    return path + newExt;
  }
  return path.substr(0, dot) + newExt;
}

static void WriteSummaryPdf(const char *pdfName,
                            TH2D &hSlope,
                            TH2D &hOffset,
                            TH2D &hAdc1,
                            TH2D &hAdc2,
                            const std::array<std::array<PixelCalib, kNY>, kNX> &calib,
                            double energy1KeV,
                            double energy2KeV) {
  TCanvas c("c_calib_summary", "Pixel calibration summary", 1400, 1000);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111);
  gStyle->SetPaintTextFormat("g");

  c.Print(Form("%s[", pdfName));

  c.Clear();
  hAdc1.Draw("COLZ TEXT");
  c.Print(pdfName);

  c.Clear();
  hAdc2.Draw("COLZ TEXT");
  c.Print(pdfName);

  // One page per X strip: 4x4 pads of Energy vs ADC linear fits for each Y strip.
  std::vector<std::unique_ptr<TGraph>> graphs;
  std::vector<std::unique_ptr<TF1>> fits;
  graphs.reserve(kNX * kNY);
  fits.reserve(kNX * kNY);

  for (int ix = 0; ix < kNX; ++ix) {
    c.Clear();
    c.Divide(4, 4);
    for (int iy = 0; iy < kNY; ++iy) {
      c.cd(iy + 1);
      if (!calib[ix][iy].valid) {
        continue;
      }

      const double adc1 = calib[ix][iy].adc1;
      const double adc2 = calib[ix][iy].adc2;
      const double adcMin = std::min(adc1, adc2);
      const double adcMax = std::max(adc1, adc2);
      const double margin = std::max(20.0, 0.1 * (adcMax - adcMin));

      auto graph = std::make_unique<TGraph>(2);
      graph->SetName(Form("gCalib_%02d_%02d", ix, iy));
      graph->SetTitle(Form("Pixel X=%02d Y=%02d;ADC;Energy [keV]", ix, iy));
      graph->SetPoint(0, adc1, energy1KeV);
      graph->SetPoint(1, adc2, energy2KeV);
      graph->SetMarkerStyle(20);
      graph->SetMarkerSize(1.2);
      graph->GetXaxis()->SetLimits(adcMin - margin, adcMax + margin);
      graph->SetMinimum(std::min(energy1KeV, energy2KeV) * 0.9);
      graph->SetMaximum(std::max(energy1KeV, energy2KeV) * 1.1);
      graph->Draw("AP");

      auto fit = std::make_unique<TF1>(Form("fCalib_%02d_%02d", ix, iy), "pol1", adcMin - margin, adcMax + margin);
      fit->SetParameters(calib[ix][iy].offset, calib[ix][iy].slope);
      fit->SetLineColor(kRed);
      fit->Draw("SAME");

      graphs.push_back(std::move(graph));
      fits.push_back(std::move(fit));
    }
    c.Print(pdfName);
  }

  gStyle->SetOptFit(0);

  c.Clear();
  hSlope.Draw("COLZ TEXT");
  c.Print(pdfName);

  c.Clear();
  hOffset.Draw("COLZ TEXT");
  c.Print(pdfName);

  c.Print(Form("%s]", pdfName));
}

int calibrate_pixel_energy(const char *peaksFile1 = "AnaW1Pix/thick_meas_1_w1_pixel_peaks.txt",
                           double energy1KeV = 5486.0,
                           const char *peaksFile2 = "AnaW1Pix/thick_meas_4_w1_pixel_peaks.txt",
                           double energy2KeV = 3180.0,
                           const char *outputFile = "CalPixE/pixel_calib.txt") {
  if (energy1KeV == energy2KeV) {
    std::cerr << "energy1 and energy2 must differ" << std::endl;
    return 1;
  }

  std::array<std::array<PixelPeak, kNY>, kNX> peaks1{};
  std::array<std::array<PixelPeak, kNY>, kNX> peaks2{};
  if (!LoadPeaksFile(peaksFile1, peaks1) || !LoadPeaksFile(peaksFile2, peaks2)) {
    return 2;
  }

  std::ofstream out(outputFile);
  if (!out) {
    std::cerr << "Cannot write calibration file: " << outputFile << std::endl;
    return 3;
  }

  const std::string rootOutput = ReplaceExtension(outputFile, ".root");
  const std::string pdfOutput = ReplaceExtension(outputFile, ".pdf");

  TH2D hSlope("hSlope", "Pixel calibration slope [keV/ADC];X;Y",
              kNX, -0.5, kNX - 0.5, kNY, -0.5, kNY - 0.5);
  TH2D hOffset("hOffset", "Pixel calibration offset [keV];X;Y",
               kNX, -0.5, kNX - 0.5, kNY, -0.5, kNY - 0.5);
  TH2D hAdc1("hAdc1", Form("Peak ADC (E=%.3g keV);X;Y", energy1KeV),
             kNX, -0.5, kNX - 0.5, kNY, -0.5, kNY - 0.5);
  TH2D hAdc2("hAdc2", Form("Peak ADC (E=%.3g keV);X;Y", energy2KeV),
             kNX, -0.5, kNX - 0.5, kNY, -0.5, kNY - 0.5);

  std::array<std::array<PixelCalib, kNY>, kNX> calib{};

  out << "# x y slope_keV_per_adc offset_keV adc1 adc2 energy1_kev energy2_kev\n";
  out << "# peaks1=" << peaksFile1 << " E1=" << energy1KeV << " keV\n";
  out << "# peaks2=" << peaksFile2 << " E2=" << energy2KeV << " keV\n";
  out << "# Energy_keV = offset_keV + slope_keV_per_adc * ADC_x\n";

  int calibrated = 0;
  int skippedIdentical = 0;
  int only1 = 0;
  int only2 = 0;

  for (int ix = 0; ix < kNX; ++ix) {
    for (int iy = 0; iy < kNY; ++iy) {
      const bool has1 = peaks1[ix][iy].valid;
      const bool has2 = peaks2[ix][iy].valid;
      if (has1 && !has2) {
        ++only1;
        continue;
      }
      if (!has1 && has2) {
        ++only2;
        continue;
      }
      if (!has1 && !has2) {
        continue;
      }

      const double adc1 = peaks1[ix][iy].peak;
      const double adc2 = peaks2[ix][iy].peak;
      if (adc1 == adc2) {
        ++skippedIdentical;
        continue;
      }

      const double slope = (energy1KeV - energy2KeV) / (adc1 - adc2);
      const double offset = energy1KeV - slope * adc1;

      out << ix << ' ' << iy << ' ' << slope << ' ' << offset << ' ' << adc1 << ' ' << adc2 << ' '
          << energy1KeV << ' ' << energy2KeV << '\n';

      calib[ix][iy].slope = slope;
      calib[ix][iy].offset = offset;
      calib[ix][iy].adc1 = adc1;
      calib[ix][iy].adc2 = adc2;
      calib[ix][iy].valid = true;

      hSlope.SetBinContent(ix + 1, iy + 1, slope);
      hOffset.SetBinContent(ix + 1, iy + 1, offset);
      hAdc1.SetBinContent(ix + 1, iy + 1, adc1);
      hAdc2.SetBinContent(ix + 1, iy + 1, adc2);
      ++calibrated;
    }
  }

  WriteSummaryPdf(pdfOutput.c_str(), hSlope, hOffset, hAdc1, hAdc2, calib, energy1KeV, energy2KeV);

  std::cout << "Wrote " << outputFile << std::endl;
  std::cout << "Wrote " << pdfOutput << std::endl;
  std::cout << "Calibrated pixels: " << calibrated << std::endl;
  if (skippedIdentical > 0) {
    std::cout << "Skipped pixels (identical ADC peaks): " << skippedIdentical << std::endl;
  }
  if (only1 > 0 || only2 > 0) {
    std::cout << "Pixels only in peaks1: " << only1 << "; only in peaks2: " << only2 << std::endl;
  }
  if (calibrated == 0) {
    std::cerr << "No common pixels with valid peaks in both files" << std::endl;
    return 4;
  }
  return 0;
}
