#!/usr/bin/env bash
set -euo pipefail

# Target run prefix under ../DAQ (directories like thick_meas, thick_meas_1, ...)
RUN_PREFIX="thick_meas"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DAQ_DIR="$(cd "${SCRIPT_DIR}/../DAQ" && pwd)"
cd "${SCRIPT_DIR}"

usage() {
  echo "Usage: $0 [run_number]" >&2
  echo "  With no argument: use the latest ${RUN_PREFIX}* under DAQ." >&2
  echo "  With run_number N: use ${RUN_PREFIX}_N (e.g. $0 2 -> ${RUN_PREFIX}_2)." >&2
}

run_number=""
if [[ $# -ge 1 ]]; then
  if [[ $# -gt 1 ]] || [[ ! "$1" =~ ^[0-9]+$ ]]; then
    usage
    exit 1
  fi
  run_number="$1"
  target_dir="${DAQ_DIR}/${RUN_PREFIX}_${run_number}"
  if [[ ! -d "${target_dir}" ]]; then
    echo "Error: target directory not found: ${target_dir}" >&2
    exit 1
  fi
else
  # Latest matching DAQ directory by modification time
  mapfile -t candidates < <(
    find "${DAQ_DIR}" -mindepth 1 -maxdepth 1 -type d \( \
      -name "${RUN_PREFIX}" -o -name "${RUN_PREFIX}_*" \
    \) -printf '%T@ %p\n' | sort -nr | awk '{print $2}'
  )

  if [[ ${#candidates[@]} -eq 0 ]]; then
    echo "Error: no directories matching '${RUN_PREFIX}' under ${DAQ_DIR}" >&2
    exit 1
  fi

  target_dir="${candidates[0]}"
  target_base="$(basename "${target_dir}")"
  if [[ "${target_base}" =~ ^${RUN_PREFIX}_([0-9]+)$ ]]; then
    run_number="${BASH_REMATCH[1]}"
  fi
fi

mapfile -t root_files < <(compgen -G "${target_dir}/RAW/DataR_*.root" || true)
if [[ ${#root_files[@]} -eq 0 ]]; then
  echo "Error: no DataR_*.root under ${target_dir}/RAW" >&2
  exit 1
fi
if [[ ${#root_files[@]} -gt 1 ]]; then
  echo "Warning: multiple DataR_*.root found; using ${root_files[0]}" >&2
fi
root_file="${root_files[0]}"

link_name="current.root"
if [[ -L "${link_name}" ]]; then
  unlink "${link_name}"
elif [[ -e "${link_name}" ]]; then
  echo "Error: ${link_name} exists and is not a symlink" >&2
  exit 1
fi
ln -s "${root_file}" "${link_name}"

if [[ -n "${run_number}" ]]; then
  summary_pdf="AnaW1Pix/${RUN_PREFIX}_${run_number}_w1_pixel_summary.pdf"
else
  summary_pdf="AnaW1Pix/${RUN_PREFIX}_w1_pixel_summary.pdf"
fi

echo "Target directory: ${target_dir}"
echo "ROOT file:        ${root_file}"
echo "Symlink:          ${link_name} -> $(readlink "${link_name}")"
echo "Summary PDF:      ${summary_pdf}"
echo "Peaks TXT:        ${summary_pdf%_summary.pdf}_peaks.txt"
echo "Running analysis..."

root -l -b -q "analyze_w1_pixels.C(\"current.root\", 0.0, 50.0, 0.0, 1500.0, 1.0, \"channel_map.txt\", \"${summary_pdf}\")"
# int analyze_w1_pixels(inputFile, coincidenceGatePs, adcThreshold, energyMin, energyMax, energyBinWidth, channelMapFile, summaryPdf)
# Writes companion peaks file: <stem>_peaks.txt next to the summary PDF.
