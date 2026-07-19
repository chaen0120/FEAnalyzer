#!/usr/bin/env bash
set -euo pipefail

# Apply CalPixE/pixel_calib.txt to one or more runs and check energy uniformity.
# Usage:
#   ./run_uniformity.sh 1
#   ./run_uniformity.sh 1 4
#   ./run_uniformity.sh            # uses current.root only

RUN_PREFIX="Mg_thick_meas"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DAQ_DIR="$(cd "${SCRIPT_DIR}/../DAQ" && pwd)"
cd "${SCRIPT_DIR}"

mkdir -p UnifPix

CALIB_FILE="CalPixE/pixel_calib.txt"
if [[ ! -f "${CALIB_FILE}" ]]; then
  echo "Error: missing ${CALIB_FILE}" >&2
  exit 1
fi

input_files=()
label_parts=()

if [[ $# -eq 0 ]]; then
  if [[ ! -e current.root ]]; then
    echo "Error: current.root not found; pass run numbers or create the symlink" >&2
    exit 1
  fi
  input_files+=("current.root")
  label_parts+=("current")
else
  for run_number in "$@"; do
    if [[ ! "${run_number}" =~ ^[0-9]+$ ]]; then
      echo "Usage: $0 [run_number ...]" >&2
      exit 1
    fi
    target_dir="${DAQ_DIR}/${RUN_PREFIX}_${run_number}"
    mapfile -t root_files < <(compgen -G "${target_dir}/RAW/DataR_*.root" || true)
    if [[ ${#root_files[@]} -eq 0 ]]; then
      echo "Error: no DataR_*.root under ${target_dir}/RAW" >&2
      exit 1
    fi
    input_files+=("${root_files[0]}")
    label_parts+=("${run_number}")
  done
fi

ifs_join() {
  local IFS=','
  echo "$*"
}

label="$(IFS=_; echo "${label_parts[*]}")"
summary_pdf="UnifPix/pixel_uniformity_${label}.pdf"
input_csv="$(ifs_join "${input_files[@]}")"

echo "Calib:   ${CALIB_FILE}"
echo "Inputs:  ${input_csv}"
echo "PDF:     ${summary_pdf}"

root -l -b -q "check_pixel_uniformity.C(\"${CALIB_FILE}\", \"${input_csv}\", 0.0, 50.0, 0.0, 7000.0, 1.0, 500, \"channel_map.txt\", \"${summary_pdf}\")"
#int check_pixel_uniformity(calibFile, inputFiles, coincidenceGatePs, adcThreshold, energyMinKeV, energyMaxKeV, energyBinWidthKeV, minPixelEntries, channelMapFile, summaryPdf)
