#ifndef X6_CALIBRATION_H
#define X6_CALIBRATION_H

#include <RtypesCore.h>

#include <string>
#include <vector>

// X6 silicon detector calibration helpers (CAEN V274X + CoMPASS RAW trees).
namespace X6Calib {

// -----------------------------------------------------------------------------
// Geometry / histogram constants
// -----------------------------------------------------------------------------
constexpr int kNumDetectors = 9;
constexpr int kNumStrips = 8;
constexpr int kNumPins = 2;
constexpr int kNumOhmic = 4;
constexpr int kEnergyBins = 512;
constexpr double kEnergyMin = 0.0;
constexpr double kEnergyMax = 2048.0;

// -----------------------------------------------------------------------------
// Run configuration (loaded from text config file)
// -----------------------------------------------------------------------------
struct RunConfig {
  std::string inputPattern =
      "../DAQ/X6decay_test_2/RAW/DataR_X6decay_test_0002*.root";
  std::string outputDir = "output";
  std::string mappingFile = "mapping.txt";
  bool useEnergyShort = false;
  double coincidenceWindow = -1.0;  // -1 => auto from timing fit
  double nSigma = 5.0;
  double adcThreshold = 50.0;
  int nPeaks = 1;
  std::vector<double> peakEnergiesMeV = {5.486};
  ULong64_t timingProbeWindow = 400000;  // ~400 ns if timestamp unit is ps
};

// -----------------------------------------------------------------------------
// Data structures
// -----------------------------------------------------------------------------
struct Event {
  UShort_t board = 0;
  UShort_t channel = 0;
  ULong64_t timestamp = 0;
  UShort_t energy = 0;
  UShort_t energyShort = 0;
  UInt_t flags = 0;
};

struct JunctionPair {
  int board = -1;
  int detector = -1;
  int stripIndex = -1;
  int chA = -1;
  int chB = -1;
  int pinA = 0;
  int pinB = 1;
};

struct OhmicChannel {
  int ohmicIndex = -1;
  int digiChannel = -1;
};

struct DetectorMap {
  int board = -1;
  int detector = -1;
  std::vector<JunctionPair> pairs;
  OhmicChannel ohmic[4];
};

struct ChannelInfo {
  int board = -1;
  int digiChannel = -1;
  int detector = -1;
  bool isOhmic = false;
  int strip = -1;
  int pin = -1;
  int ohmicIndex = -1;
};

struct CoincidenceEvent {
  int board = -1;
  int detector = -1;
  int stripIndex = -1;
  int ohmicIndex = -1;
  double energyPin0 = 0.0;
  double energyPin1 = 0.0;
  double ohmicEnergy = 0.0;
  double junctionSum = 0.0;
  ULong64_t anchorTs = 0;
  double dtOhmicPin0 = 0.0;
  double dtOhmicPin1 = 0.0;
  double dtPin0Pin1 = 0.0;
  double maxSpan = 0.0;
};

struct TimingDiagnostics {
  std::vector<double> dtJunctionPair;
  std::vector<double> dtOhmicPin0;
  std::vector<double> dtOhmicPin1;
  std::vector<double> dtOhmicJun;
  std::vector<double> dtPin0Pin1;
  std::vector<double> maxSpan;
};

struct PeakFitResult {
  int detector = -1;
  int index = -1;  // ohmic 0-3 or junction strip 0-7
  double entries = 0.0;
  std::vector<double> meansAdc;
  std::vector<double> sigmasAdc;
  std::vector<double> energiesMeV;
  int nPeaksFitted = 0;
  bool fitted = false;
  // ADC = calibA + calibB * E[MeV]
  double calibA = 0.0;
  double calibB = 0.0;
  bool calibrated = false;
};

// -----------------------------------------------------------------------------
// Public API
// -----------------------------------------------------------------------------
bool LoadConfig(const std::string &path, RunConfig &cfg);
void PrintConfig(const RunConfig &cfg);
void RunCalibration(const RunConfig &cfg);

}  // namespace X6Calib

#endif  // X6_CALIBRATION_H
