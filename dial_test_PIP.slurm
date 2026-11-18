#!/bin/bash
#SBATCH --job-name=OT_PIP_test
#SBATCH -D .                               # run from the submission directory
#SBATCH -p devel                           # DIaL3 devel partition (max 4 nodes, 2 h)
#SBATCH -A dp446                           # DiRAC project account
#SBATCH --time=00:45:00
#SBATCH --nodes=2                          # 2 nodes x 128 ranks = 256 ranks
#SBATCH --ntasks-per-node=128
#SBATCH --exclusive                        # devel nodes are shared by default; exclusive for clean timings

# Test run for the 32768^2 PIP Orszag-Tang job (dial_OT_PIP_32768.slurm).
#
# Runs 8192^2 on 256 ranks (16 x 16). Each rank gets a 512 x 512 block, the same as
# 32768^2 on 4096 ranks, so the time per step here is what the full run will see.
# Checks: the build, the MPI launch across nodes, HDF5 snapshot output (t=0, 0.01, 0.02).
# Prints: measured time per step, dt, and a projection for the 32768^2 run.
#
# Expects in the submission directory: dial_env.sh, a.out (make -f Makefile.dial),
# setting_dial_OT_PIP_32768.txt (physics settings are taken from it unchanged).

source ./dial_env.sh

set -euo pipefail

NTEST=8192
RUNDIR=test_PIP_${NTEST}_${SLURM_JOB_ID}
mkdir -p "$RUNDIR"
cp a.out "$RUNDIR/"

# Same physics as the production run; only grid, layout, run length and output change.
sed -e "1s#.*#OT_PIP_${NTEST}/:outdir; output directory#" \
    -e "s/^32768:ix;/${NTEST}:ix;/" \
    -e "s/^32768:jx;/${NTEST}:jx;/" \
    -e "s/^1.0:tend;/0.02:tend;/" \
    -e "s/^0.1 :dtout;/0.01 :dtout;/" \
    -e "s/^1:esav;/0:esav;/" \
    -e "s/^64     :mpi_x;/16     :mpi_x;/" \
    -e "s/^64     :mpi_y;/16     :mpi_y;/" \
    setting_dial_OT_PIP_32768.txt > "$RUNDIR/setting.txt"
cd "$RUNDIR"

# Check every substitution took effect (fails loudly if the template changes).
for pat in "^OT_PIP_${NTEST}/:outdir;" "^PIP:flag_eqs;" "^${NTEST}:ix;" "^${NTEST}:jx;" \
           "^0.02:tend;" "^0.01 :dtout;" "^0:esav;" "^16     :mpi_x;" "^16     :mpi_y;"; do
    grep -q "$pat" setting.txt || { echo "ERROR: setting.txt is missing '$pat'" >&2; exit 1; }
done

OUTDIR=OT_PIP_${NTEST}
mkdir -p "$OUTDIR"

echo "Job $SLURM_JOB_ID on $SLURM_JOB_NUM_NODES nodes ($SLURM_JOB_NODELIST), $SLURM_NTASKS ranks, started $(date)"
echo "h5pfc: $(which h5pfc)   mpirun: $(which mpirun)"
START=$(date +%s)
set +e
mpirun -np "$SLURM_NTASKS" --bind-to core -x MALLOC_MMAP_MAX_ -x MALLOC_TRIM_THRESHOLD_ ./a.out > "$OUTDIR/run.log" 2>&1
STATUS=$?
set -e
echo "mpirun exit status $STATUS after $(( $(date +%s) - START )) s"

echo
echo "== Snapshots written"
ls -lh "$OUTDIR"/*.h5 2>/dev/null || echo "NONE - check HDF5 output in $RUNDIR/$OUTDIR/run.log"

echo
echo "== Problems in the log"
grep -iE "nan|error|abort|segmentation|killed" "$OUTDIR/run.log" | grep -v "cp: cannot create" | head -10 || echo "none found"

echo
echo "== Progress lines (time, dt, nout, nt, elapsed s)"
grep "Time,dt" "$OUTDIR/run.log" || true

# Timing from the last two snapshot lines: steps and wall time between them.
grep "Time,dt" "$OUTDIR/run.log" | sed 's/.*://' | tail -2 | awk -v n=$NTEST '
    NR==1 {t1=$1; nt1=$4; e1=$5}
    NR==2 {t2=$1; dt=$2; nt2=$4; e2=$5}
    END {
        if (NR < 2 || nt2 <= nt1) { print "\nNot enough progress lines to measure timing."; exit }
        sps   = (e2 - e1) / (nt2 - nt1)          # wall seconds per step (includes one snapshot write)
        dtbig = dt * n / 32768                    # dt scales with grid spacing
        steps = 1.0 / dtbig                       # to tend = 1, if dt stayed at its current value
        wall  = steps * sps / 3600
        printf "\n== Measured on %d^2\n", n
        printf "  seconds per step : %.3f   (%d steps between t=%.3f and t=%.3f)\n", sps, nt2-nt1, t1, t2
        printf "  dt               : %.3e\n", dt
        printf "\n== Projection for 32768^2 on 4096 ranks (32 nodes), same 512x512 block per rank\n"
        printf "  steps to t=1     : ~%.0f   (more if dt shrinks as shocks form)\n", steps
        printf "  wall time        : ~%.0f h  (job asks for 72 h)\n", wall
        printf "  core-hours       : ~%.0f   (allocation ~400000)\n", wall * 4096
    }'
