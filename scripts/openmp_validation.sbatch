#!/usr/bin/env bash
#SBATCH --job-name=pb-openmp-val
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=validation/slurm-%j.out
#SBATCH --error=validation/slurm-%j.err

set -uo pipefail

# Compare the original upstream code in serial against the optimized code in
# both serial and OpenMP modes. This lets us confirm that threaded runs still
# match the validated serial baseline while reporting speedups.
#
# By default this script runs all bundled case studies and preserves the
# case-study defaults exactly as shipped in case_studies/<case>/defaults.dat.
# That means surface-area calculations remain enabled wherever the example
# inputs enable them.

if [ -n "${REPO_DIR:-}" ]; then
  REPO_DIR=$(cd "$REPO_DIR" && pwd)
elif [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
  REPO_DIR=$(cd "$SLURM_SUBMIT_DIR" && pwd)
else
  SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
  REPO_DIR=$(cd "$SCRIPT_DIR/.." && pwd)
fi

OLD_REF=${OLD_REF:-upstream/main}
NEW_REF=${NEW_REF:-WORKTREE}
OUT_ROOT=${OUT_ROOT:-"$REPO_DIR/validation"}
FORTRAN_COMPILER=${FORTRAN_COMPILER:-gfortran}
JOB_ID=${SLURM_JOB_ID:-manual}
RUN_ROOT=${RUN_ROOT:-"${SLURM_TMPDIR:-/tmp}/poreblazer_openmp_${JOB_ID}_$$"}
KEEP_FULL_OUTPUTS=${KEEP_FULL_OUTPUTS:-0}
OLD_THREADS=${OLD_THREADS:-1}
NEW_SERIAL_THREADS=${NEW_SERIAL_THREADS:-1}
NEW_OPENMP_THREADS=${NEW_OPENMP_THREADS:-${SLURM_CPUS_PER_TASK:-16}}
WORKTREE_PATCH="$RUN_ROOT/new_worktree.patch"

CASES_STRING=${CASES:-"BHP CD121 CLO HKUST1 IRMOF1 LOV MIL101 MIL47V MOF180 ROG RON SLIT STO WEI ZIF8"}
read -r -a CASES <<< "$CASES_STRING"

mkdir -p "$OUT_ROOT"
OUT_DIR="$OUT_ROOT/${JOB_ID}_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUT_DIR"

OLD_SHA=$(git -C "$REPO_DIR" rev-parse --verify "$OLD_REF^{commit}")
if [ "$NEW_REF" = "WORKTREE" ]; then
  NEW_SHA=$(git -C "$REPO_DIR" rev-parse --verify "HEAD^{commit}")
else
  NEW_SHA=$(git -C "$REPO_DIR" rev-parse --verify "$NEW_REF^{commit}")
fi

OLD_TREE="$RUN_ROOT/source_old"
NEW_TREE="$RUN_ROOT/source_new"
mkdir -p "$RUN_ROOT"

cleanup() {
  git -C "$REPO_DIR" worktree remove --force "$OLD_TREE" >/dev/null 2>&1 || true
  git -C "$REPO_DIR" worktree remove --force "$NEW_TREE" >/dev/null 2>&1 || true
}
trap cleanup EXIT

echo "PoreBlazer OpenMP validation"
echo "Repository: $REPO_DIR"
echo "Cases: ${CASES[*]}"
echo "Old ref: $OLD_REF ($OLD_SHA), threads=$OLD_THREADS"
if [ "$NEW_REF" = "WORKTREE" ]; then
  echo "New ref: WORKTREE (base $NEW_SHA + local tracked diff)"
else
  echo "New ref: $NEW_REF ($NEW_SHA)"
fi
echo "New serial threads: $NEW_SERIAL_THREADS"
echo "New OpenMP threads: $NEW_OPENMP_THREADS"
echo "Fortran compiler: $FORTRAN_COMPILER"
echo "Run root: $RUN_ROOT"
echo "Output: $OUT_DIR"

git -C "$REPO_DIR" worktree add --detach "$OLD_TREE" "$OLD_SHA" || exit 1
git -C "$REPO_DIR" worktree add --detach "$NEW_TREE" "$NEW_SHA" || exit 1

if [ "$NEW_REF" = "WORKTREE" ]; then
  git -C "$REPO_DIR" diff --binary HEAD -- README.md .gitignore \
    'src/*.f90' 'src/Makefile*' 'src/README_OPENMP.txt' > "$WORKTREE_PATCH"
  if [ -s "$WORKTREE_PATCH" ]; then
    git -C "$NEW_TREE" apply "$WORKTREE_PATCH" || exit 1
  fi
fi

build_tree() {
  local label=$1
  local tree=$2
  local log_dir="$OUT_DIR/build_$label"

  mkdir -p "$log_dir"
  echo "Building $label in $tree"
  (
    cd "$tree/src" || exit 1
    make clean >"$log_dir/make-clean.log" 2>&1 || true
    make FORTRAN_COMPILER="$FORTRAN_COMPILER" >"$log_dir/make.log" 2>&1
  )
}

run_case() {
  local label=$1
  local tree=$2
  local case_name=$3
  local threads=$4
  local case_src="$REPO_DIR/case_studies/$case_name"
  local run_dir="$RUN_ROOT/runs/$label/$case_name"
  local out_dir="$OUT_DIR/$case_name/$label"
  local status=0

  rm -rf "$run_dir"
  mkdir -p "$run_dir" "$out_dir"
  cp -a "$case_src"/. "$run_dir"/

  echo "Running $label/$case_name with OMP_NUM_THREADS=$threads"
  (
    cd "$run_dir" || exit 1
    export OMP_NUM_THREADS="$threads"
    if [ "$threads" -gt 1 ]; then
      export OMP_PROC_BIND=close
      export OMP_PLACES=cores
    fi
    /usr/bin/time -p -o "$out_dir/time.txt" "$tree/src/poreblazer.exe" < input.dat > "$out_dir/stdout.txt" 2> "$out_dir/stderr.txt"
  )
  status=$?

  printf "%s\n" "$status" > "$out_dir/status.txt"
  if [ -f "$run_dir/summary.dat" ]; then
    cp "$run_dir/summary.dat" "$out_dir/summary.dat"
  fi

  if [ "$KEEP_FULL_OUTPUTS" = "1" ]; then
    mkdir -p "$out_dir/full"
    cp -a "$run_dir"/. "$out_dir/full"/
  fi

  return "$status"
}

safe_real() {
  awk '/^real / {print $2}' "$1" 2>/dev/null || true
}

safe_speedup() {
  python - "$1" "$2" <<'PY'
import sys
base = sys.argv[1]
trial = sys.argv[2]
try:
    b = float(base)
    t = float(trial)
    if t == 0:
        print("NA")
    else:
        print(f"{b/t:.2f}")
except Exception:
    print("NA")
PY
}

summary_without_metadata() {
  local summary_file=$1
  awk '!/^D_axes[[:space:]]/' "$summary_file"
}

build_tree old "$OLD_TREE" || {
  echo "Old-code build failed. See $OUT_DIR/build_old/make.log"
  exit 1
}

build_tree new "$NEW_TREE" || {
  echo "New-code build failed. See $OUT_DIR/build_new/make.log"
  exit 1
}

SUMMARY_CSV="$OUT_DIR/openmp_validation.csv"
printf "case,old_status,new_serial_status,new_openmp_status,old_real_s,new_serial_real_s,new_openmp_real_s,serial_diff_lines,openmp_diff_lines,openmp_vs_serial_diff_lines,openmp_speedup_vs_old,openmp_speedup_vs_new_serial\n" > "$SUMMARY_CSV"

for case_name in "${CASES[@]}"; do
  old_status=0
  new_serial_status=0
  new_openmp_status=0

  run_case old "$OLD_TREE" "$case_name" "$OLD_THREADS" || old_status=$?
  run_case new_serial "$NEW_TREE" "$case_name" "$NEW_SERIAL_THREADS" || new_serial_status=$?
  run_case new_openmp "$NEW_TREE" "$case_name" "$NEW_OPENMP_THREADS" || new_openmp_status=$?

  old_summary="$OUT_DIR/$case_name/old/summary.dat"
  new_serial_summary="$OUT_DIR/$case_name/new_serial/summary.dat"
  new_openmp_summary="$OUT_DIR/$case_name/new_openmp/summary.dat"

  serial_diff_file="$OUT_DIR/$case_name/serial.diff"
  openmp_diff_file="$OUT_DIR/$case_name/openmp.diff"
  openmp_vs_serial_diff_file="$OUT_DIR/$case_name/openmp_vs_serial.diff"

  if [ -f "$old_summary" ] && [ -f "$new_serial_summary" ]; then
    diff -u \
      <(summary_without_metadata "$old_summary") \
      <(summary_without_metadata "$new_serial_summary") \
      > "$serial_diff_file" || true
    serial_diff_lines=$(wc -l < "$serial_diff_file")
  else
    printf "summary.dat missing for old or new_serial run\n" > "$serial_diff_file"
    serial_diff_lines=-1
  fi

  if [ -f "$old_summary" ] && [ -f "$new_openmp_summary" ]; then
    diff -u \
      <(summary_without_metadata "$old_summary") \
      <(summary_without_metadata "$new_openmp_summary") \
      > "$openmp_diff_file" || true
    openmp_diff_lines=$(wc -l < "$openmp_diff_file")
  else
    printf "summary.dat missing for old or new_openmp run\n" > "$openmp_diff_file"
    openmp_diff_lines=-1
  fi

  if [ -f "$new_serial_summary" ] && [ -f "$new_openmp_summary" ]; then
    diff -u "$new_serial_summary" "$new_openmp_summary" > "$openmp_vs_serial_diff_file" || true
    openmp_vs_serial_diff_lines=$(wc -l < "$openmp_vs_serial_diff_file")
  else
    printf "summary.dat missing for new_serial or new_openmp run\n" > "$openmp_vs_serial_diff_file"
    openmp_vs_serial_diff_lines=-1
  fi

  old_real=$(safe_real "$OUT_DIR/$case_name/old/time.txt")
  new_serial_real=$(safe_real "$OUT_DIR/$case_name/new_serial/time.txt")
  new_openmp_real=$(safe_real "$OUT_DIR/$case_name/new_openmp/time.txt")

  openmp_speedup_vs_old=$(safe_speedup "${old_real:-NA}" "${new_openmp_real:-NA}")
  openmp_speedup_vs_new_serial=$(safe_speedup "${new_serial_real:-NA}" "${new_openmp_real:-NA}")

  printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" \
    "$case_name" "$old_status" "$new_serial_status" "$new_openmp_status" \
    "${old_real:-NA}" "${new_serial_real:-NA}" "${new_openmp_real:-NA}" \
    "$serial_diff_lines" "$openmp_diff_lines" "$openmp_vs_serial_diff_lines" \
    "$openmp_speedup_vs_old" "$openmp_speedup_vs_new_serial" \
    >> "$SUMMARY_CSV"
done

echo "OpenMP validation complete"
echo "Summary: $SUMMARY_CSV"
