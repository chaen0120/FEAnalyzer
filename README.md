# FEAnalyzer

ROOT-based analysis tools for front-end silicon detectors read out with
**CAEN digitizers + CoMPASS** (`Data_R` trees).

This repository currently includes:

| Package | Purpose |
|---------|---------|
| `X6/` | X6 detector energy calibration (timing, coincidences, multi-peak fits, ADC→MeV) |
| `W1_ThickMeas/` | W1 DSSD pixel analysis for thickness measurement (peaks, energy calib, uniformity) |

## Usage with CoMPASS projects

Macros expect `../DAQ` relative to the package directory (`X6/`, `W1_ThickMeas/`).
Symlink the package you need into a CoMPASS project (next to that project’s `DAQ/`):

```bash
cd /path/to/<CoMPASS_project>
ln -s /path/to/FEAnalyzer/X6 X6
# ln -s /path/to/FEAnalyzer/W1_ThickMeas W1_ThickMeas
```

Then run from the symlink, e.g. `cd X6 && root -l -b -q run_calibration.C`.

Configs live in this repository and are shared across projects that link the same
package. For project-specific settings, pass a custom config path to the entry
macro.

## Requirements

- [ROOT](https://root.cern/) with Cling (interpreted macros)
- Bash (for the W1 helper scripts)
- CoMPASS RAW ROOT files with tree `Data_R` (branches such as `Board`, `Channel`, `Timestamp`, `Energy`, …)

## Repository layout

```
FEAnalyzer/
├── README.md
├── .gitignore
├── X6/
│   ├── run_calibration.C      # entry point
│   ├── x6_calibration.h
│   ├── x6_calibration_lib.C
│   └── config/
│       ├── config.txt         # run settings
│       └── mapping.txt        # board/channel → detector map
└── W1_ThickMeas/
    ├── analyze_w1_pixels.C
    ├── calibrate_pixel_energy.C
    ├── check_pixel_uniformity.C
    ├── strip_map.h
    ├── channel_map.txt        # board/channel → X/Y strip map
    ├── run_analyze.sh
    └── run_uniformity.sh
```

Generated plots, histograms, and peak tables are written under local output
folders and are ignored by git (see `.gitignore`).

---

## X6 calibration

Builds 2-junction + 1-ohmic coincidences, fits calibration peaks, and produces
linear coefficients `ADC = A + B · E[MeV]`.

### Configure

Edit `X6/config/config.txt`:

- `input_pattern` — glob of CoMPASS RAW files  
- `mapping_file` — channel map (`config/mapping.txt`)  
- `adc_threshold` — ADC cut  
- `n_peaks` / `peak_energies_mev` — number of peaks and energies in MeV  
- `coincidence_window` — timestamp window (`-1` = auto from timing fit)

Edit `X6/config/mapping.txt` for your cabling:

```
# board channel detector J/O det_ana strip pin
```

- Junction strips: `strip` 0–7, `pin` 0/1  
- Ohmic: `strip` 0–3, `pin` -1  

### Run

```bash
cd X6    # package symlink inside your CoMPASS project
root -l -b -q run_calibration.C
# or a custom config:
root -l -b -q 'run_calibration.C("config/config.txt")'
```

### Outputs (`X6/output/`)

| File | Contents |
|------|----------|
| `dt_timing.pdf` | Timing diagnostics |
| `ohmic_peaks.pdf` / `junction_peaks.pdf` | Peak spectra + fits (one detector/page) |
| `junction_epin0_epin1.pdf` | Junction \(E_{\mathrm{pin0}}\) vs \(E_{\mathrm{pin1}}\) (one detector/page) |
| `junction_position_energy.pdf` | Position \(\eta=(E_0-E_1)/(E_0+E_1)\) vs energy (one detector/page) |
| `calibration_linear_fits.pdf` | ADC vs E linear fits (one detector/page) |
| `x6_calibration_summary.txt` | Timing + peak info |
| `x6_calibration_coeffs.txt` | Calibration coefficients A, B |
| `x6_calibration_histograms.root` | Histograms / fits |

---

## W1 thickness measurement

Pipeline for a 16×16 double-sided strip detector (pixel = X∧Y coincidence).

### 1. Analyze a run (pixel ADC peaks)

```bash
cd W1_ThickMeas    # package symlink inside your CoMPASS project
./run_analyze.sh         # latest thick_meas* under ../DAQ
./run_analyze.sh 4       # ../DAQ/thick_meas_4
```

The script links `../DAQ/.../RAW/DataR_*.root` → `current.root` and runs
`analyze_w1_pixels.C`.

Outputs under `AnaW1Pix/`:

- `*_w1_pixel_summary.pdf`
- `*_w1_pixel_peaks.txt`

### 2. Energy calibration (two peak sets)

Use peak tables from two known energies (e.g. two source runs):

```bash
root -l -b -q 'calibrate_pixel_energy.C(
  "AnaW1Pix/thick_meas_1_w1_pixel_peaks.txt",
  "AnaW1Pix/thick_meas_4_w1_pixel_peaks.txt",
  E1_keV, E2_keV,
  "CalPixE/pixel_calib.txt",
  "CalPixE/pixel_calib.pdf")'
```

Replace `E1_keV` / `E2_keV` with the true energies for those runs.

### 3. Pixel energy uniformity

```bash
./run_uniformity.sh           # uses current.root
./run_uniformity.sh 1 4       # compare Mg_thick_meas_1 and _4
```

Requires `CalPixE/pixel_calib.txt`. Writes PDFs under `UnifPix/`.

### Channel map

Edit `W1_ThickMeas/channel_map.txt` if strip order / board assignment differs:

```
# board  channel  axis  strip
# axis = X or Y ; strip = 0 .. 15
```

---

## Notes

- Timestamp units follow the digitizer / CoMPASS export (often picoseconds for CAEN V274x).
- Run via the package symlink inside a CoMPASS project so `../DAQ` finds that project’s data.
- These are **interpreted ROOT macros** — no build step is required.
