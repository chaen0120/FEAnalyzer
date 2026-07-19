// X6 calibration entry point.
// Edit settings in x6_calibration_config.txt, then run:
//   root -l x6_calibration_macro.C
// Optional custom config:
//   root -l 'x6_calibration_macro.C("my_config.txt")'

#include "x6_calibration.h"
#include "x6_calibration_lib.C"

void run_calibration(const char *configFile = "config/config.txt") {
  X6Calib::RunConfig cfg;
  if (!X6Calib::LoadConfig(configFile, cfg)) {
    return;
  }
  X6Calib::PrintConfig(cfg);
  X6Calib::RunCalibration(cfg);
}
