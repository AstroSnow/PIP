#!/bin/bash
#SBATCH --job-name=OT_32768
#SBATCH --export=all                       # export all environment variables to the batch job
#SBATCH -D .                               # run from the submission directory
#SBATCH -p pq                              # parallel queue
#SBATCH --time=160:00:00                   # max wall time (estimate ~95 h on 496 ranks)
#SBATCH -A Research_Project-T108096        # research project to submit under
#SBATCH --nodes=31                         # 31 nodes x 16 ranks = 496 ranks (QOS normal caps users at 500 CPUs)
#SBATCH --ntasks-per-node=16               # must give nodes*ntasks = mpi_x*mpi_y (16*31) in setting.txt
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=b.snow@exeter.ac.uk    # email on completion or failure

# Orszag-Tang ideal MHD, 32768^2, SLW scheme, no physical resistivity/viscosity.
#
# Expects in the submission directory:
#   a.out                      built with: make -f Makefile.isca
#   setting_isca_OT_32768.txt  run settings (outdir, grid, mpi_x=16, mpi_y=31)
#
# The run happens in its own directory, RUNDIR, so the repo's setting.txt is never
# overwritten. To restart after a timeout, edit $RUNDIR/setting.txt:
#   flag_restart = 0                (continue from the last snapshot listed in config.txt)
#   indir        = OT_MHD_32768/    (restart data is read from indir, not outdir)
#   tend         = time remaining   (on restart the code sets tend = tend + current time)
# then resubmit with:
#   sbatch --export=ALL,RESTART=1 isca_OT_32768.slurm

module purge
module load foss/2019a
module load h5py/2.9.0-foss-2019a          # provides parallel HDF5 (h5pfc) built with foss/2019a

# Stop glibc returning large freed arrays (artvis allocates ~1 GB of automatic arrays per
# rank per call) to the OS every call. In the local 16384^2 run that churn cost 57% of
# CPU time as page faults; with these settings freed memory is reused (~1.8x faster per step).
export MALLOC_MMAP_MAX_=0
export MALLOC_TRIM_THRESHOLD_=68719476736
ulimit -s unlimited

set -euo pipefail

RUNDIR=run_OT_32768
SETTINGS=setting_isca_OT_32768.txt

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

# One line per node with an explicit slot count. Without "slots=", Open MPI assumes it can
# use every core it detects on a node, packs more than 16 ranks onto it, and then fails to
# bind ranks to cores outside this job's 16-core allocation (hwloc_set_cpubind error).
scontrol show hostnames "$SLURM_JOB_NODELIST" | sed "s/\$/ slots=$SLURM_NTASKS_PER_NODE/" > nodes.txt

echo "Job $SLURM_JOB_ID on $SLURM_JOB_NUM_NODES nodes, $SLURM_NTASKS ranks, started $(date)"
mpirun -x LD_LIBRARY_PATH -x PATH -x MALLOC_MMAP_MAX_ -x MALLOC_TRIM_THRESHOLD_ -machinefile nodes.txt --report-bindings \
       -np "$SLURM_NTASKS" --bind-to core ./a.out > "$OUTDIR/run.log" 2>&1
echo "Finished $(date)"
