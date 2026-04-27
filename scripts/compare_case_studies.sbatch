#!/usr/bin/env bash
#SBATCH --job-name=pb-compare-cases
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=validation/slurm-%j.out
#SBATCH --error=validation/slurm-%j.err

set -uo pipefail

# Compare all bundled case studies between the old upstream code and the
# current optimized code without modifying tracked case-study outputs.

if [ -n "${REPO_DIR:-}" ]; then
  REPO_DIR=$(cd "$REPO_DIR" && pwd)
elif [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
  REPO_DIR=$(cd "$SLURM_SUBMIT_DIR" && pwd)
else
  SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
  REPO_DIR=$(cd "$SCRIPT_DIR/.." && pwd)
fi

OLD_REF=${OLD_REF:-upstream/main}
NEW_REF=${NEW_REF:-HEAD}
OUT_ROOT=${OUT_ROOT:-"$REPO_DIR/validation"}
THREADS=${SLURM_CPUS_PER_TASK:-1}
FORTRAN_COMPILER=${FORTRAN_COMPILER:-gfortran}
JOB_ID=${SLURM_JOB_ID:-manual}
RUN_ROOT=${RUN_ROOT:-"${SLURM_TMPDIR:-/tmp}/poreblazer_compare_${JOB_ID}_$$"}
KEEP_FULL_OUTPUTS=${KEEP_FULL_OUTPUTS:-0}
WORKTREE_PATCH="$RUN_ROOT/new_worktree.patch"

CASES=(
  BHP
  CD121
  CLO
  HKUST1
  IRMOF1
  LOV
  MIL101
  MIL47V
  MOF180
  ROG
  RON
  SLIT
  STO
  WEI
  ZIF8
)

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

echo "PoreBlazer case-study comparison"
echo "Repository: $REPO_DIR"
echo "Old ref: $OLD_REF ($OLD_SHA)"
if [ "$NEW_REF" = "WORKTREE" ]; then
  echo "New ref: WORKTREE (base $NEW_SHA + local tracked diff)"
else
  echo "New ref: $NEW_REF ($NEW_SHA)"
fi
echo "Threads: $THREADS"
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
  local case_src="$REPO_DIR/case_studies/$case_name"
  local run_dir="$RUN_ROOT/runs/$label/$case_name"
  local out_dir="$OUT_DIR/$case_name/$label"
  local status=0

  mkdir -p "$run_dir" "$out_dir"
  cp -a "$case_src"/. "$run_dir"/

  echo "Running $label/$case_name"
  (
    cd "$run_dir" || exit 1
    export OMP_NUM_THREADS="$THREADS"
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

SUMMARY_CSV="$OUT_DIR/summary.csv"
printf "case,old_status,new_status,old_real_s,new_real_s,summary_diff_lines\n" > "$SUMMARY_CSV"

for case_name in "${CASES[@]}"; do
  old_status=0
  new_status=0

  run_case old "$OLD_TREE" "$case_name" || old_status=$?
  run_case new "$NEW_TREE" "$case_name" || new_status=$?

  old_summary="$OUT_DIR/$case_name/old/summary.dat"
  new_summary="$OUT_DIR/$case_name/new/summary.dat"
  diff_file="$OUT_DIR/$case_name/summary.diff"

  if [ -f "$old_summary" ] && [ -f "$new_summary" ]; then
    diff -u \
      <(summary_without_metadata "$old_summary") \
      <(summary_without_metadata "$new_summary") \
      > "$diff_file" || true
    diff_lines=$(wc -l < "$diff_file")
  else
    printf "summary.dat missing for old or new run\n" > "$diff_file"
    diff_lines=-1
  fi

  old_real=$(awk '/^real / {print $2}' "$OUT_DIR/$case_name/old/time.txt" 2>/dev/null || true)
  new_real=$(awk '/^real / {print $2}' "$OUT_DIR/$case_name/new/time.txt" 2>/dev/null || true)

  printf "%s,%s,%s,%s,%s,%s\n" \
    "$case_name" "$old_status" "$new_status" "${old_real:-NA}" "${new_real:-NA}" "$diff_lines" \
    >> "$SUMMARY_CSV"
done

echo "Comparison complete"
echo "Summary: $SUMMARY_CSV"
