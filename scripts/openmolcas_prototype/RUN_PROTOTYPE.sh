#!/usr/bin/env bash
set -euo pipefail

: "${MOLCAS_ROOT:?Set MOLCAS_ROOT to the OpenMolcas build or installation directory}"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
package_dir="$(cd "${script_dir}/.." && pwd)"
run_root="${1:-${package_dir}/reproduction_runs}"

mkdir -p "${run_root}/h2/work" "${run_root}/dy3/work"

env MOLCAS="${MOLCAS_ROOT}" \
  MOLCAS_WORKDIR="${run_root}/h2/work" \
  MOLCAS_MEM=2000 MOLCAS_NPROCS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  "${MOLCAS_ROOT}/pymolcas" "${package_dir}/inputs/h2_casscf_caspt2.input" \
  > "${run_root}/h2/h2_console.log" 2>&1

env MOLCAS="${MOLCAS_ROOT}" \
  MOLCAS_WORKDIR="${run_root}/dy3/work" \
  MOLCAS_MEM=6000 MOLCAS_DISK=20000 MOLCAS_NPROCS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  "${MOLCAS_ROOT}/pymolcas" "${package_dir}/inputs/dy3_cas97_soc_single_aniso_core.input" \
  > "${run_root}/dy3/dy3_console.log" 2>&1

grep -q "Happy landing" "${run_root}/h2/h2_console.log"
grep -q "Happy landing" "${run_root}/dy3/dy3_console.log"
printf 'Both prototype runs completed with Happy landing.\n'
