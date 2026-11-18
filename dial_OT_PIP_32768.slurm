#!/bin/bash
#SBATCH --job-name=OT_PIP_32768
#SBATCH -D .                               # run from the submission directory
#SBATCH -p slurm                           # DIaL3 standard partition (max 4 days, 82 nodes)
#SBATCH -A dp446                           # DiRAC project account
#SBATCH --time=72:00:00                    # estimate 23-49 h; emergency save at 70 h (esav in settings)
#SBATCH --nodes=32                         # 32 nodes x 128 ranks = 4096 ranks
#SBATCH --ntasks-per-node=128              # must give nodes*ntasks = mpi_x*mpi_y (64*64) in setting.txt
#SBATCH --exclusive
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=b.snow@exeter.ac.uk    # email on completion or failure

# Orszag-Tang partially ionised (PIP), 32768^2, SLW scheme.
# n_fraction=0.9, col=1, explicit collisions; no IR, radiation, resistivity, viscosity or conduction.
#
# Expects in the submission directory:
#   dial_env.sh                    environment (modules + parallel HDF5 from build_hdf5_dial.sh)
#   a.out                          built with: make -f Makefile.dial
#   setting_dial_OT_PIP_32768.txt  run settings (outdir, grid, mpi_x=64, mpi_y=64)
#
# The run happens in its own directory, RUNDIR, so the repo's setting.txt is never
# overwritten. If the job hits the wall time, the emergency save (esav, 70 h) will
# have written a snapshot. To restart from it, edit $RUNDIR/setting.txt:
#   flag_restart = 0                (continue from the last snapshot listed in config.txt)
#   indir        = OT_PIP_32768/    (restart data is read from indir, not outdir)
#   tend         = time remaining   (on restart the code sets tend = tend + current time)
#   esav         = 1                (re-arm the emergency save)
# then resubmit with:
#   sbatch --export=ALL,RESTART=1 dial_OT_PIP_32768.slurm

# Same compiler, MPI and HDF5 as the build (see dial_env.sh).
source ./dial_env.sh

set -euo pipefail

RUNDIR=run_OT_PIP_32768
SETTINGS=setting_dial_OT_PIP_32768.txt

if [ "${RESTART:-0}" != "1" ]; then
    # Fresh start: refuse to overwrite an earlier run's data.
    if [ -e "$RUNDIR" ]; then
        echo "ERROR: $RUNDIR already exists. Remove it, or resubmit with RESTART=1." >&2
        exit 1
    fi
    mkdir -p "$RUNDIR"
    cp "$SETTINGS" "$RUNDIR/setting.txt"
fi
cp a.out "$RUNDIR/"
cd "$RUNDIR"

# The code does not create its output directory itself.
OUTDIR=$(head -n1 setting.txt | cut -d: -f1)
mkdir -p "$OUTDIR"

# Check that the MPI layout in setting.txt matches the allocation. A mismatch
# silently falls back to a 1D strip decomposition.
MPIX=$(grep ':mpi_x;' setting.txt | cut -d: -f1 | tr -d ' ')
MPIY=$(grep ':mpi_y;' setting.txt | cut -d: -f1 | tr -d ' ')
MPIZ=$(grep ':mpi_z;' setting.txt | cut -d: -f1 | tr -d ' ')
if [ $((MPIX * MPIY * MPIZ)) -ne "$SLURM_NTASKS" ]; then
    echo "ERROR: mpi_x*mpi_y*mpi_z=$((MPIX*MPIY*MPIZ)) but SLURM_NTASKS=$SLURM_NTASKS" >&2
    exit 1
fi

echo "Job $SLURM_JOB_ID on $SLURM_JOB_NUM_NODES nodes, $SLURM_NTASKS ranks, started $(date)"
mpirun -np "$SLURM_NTASKS" --bind-to core -x MALLOC_MMAP_MAX_ -x MALLOC_TRIM_THRESHOLD_ ./a.out > "$OUTDIR/run.log" 2>&1
echo "Finished $(date)"
