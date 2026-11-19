# Environment for building and running PIP on DIaL3. Use with:  source dial_env.sh
# Used by build_hdf5_dial.sh, by "make -f Makefile.dial", and by dial_OT_PIP_32768.slurm,
# so the compiler, MPI and HDF5 are always the same.

# Keep conda's compilers/MPI/HDF5 off the PATH.
if type conda >/dev/null 2>&1; then
    conda deactivate 2>/dev/null || true
fi

module purge
# Intel compilers for the openmpi4/intel wrappers. If "mpif90 --version" fails after
# sourcing this file, try intel-parallel-studio/cluster.2019.5/aveeft here instead.
module load intel-oneapi/2022.1.2
module load openmpi4/intel/4.0.5

# openmpi4/intel/4.0.5 was built against PMIx 3 (/cm/shared/apps/cm-pmix3/3.1.4), but that
# directory no longer holds libpmix.so.2, so the loader falls back to the system PMIx 2.2.5
# in /lib64 (same soname) and mpirun aborts with "unexpected version of PMIx". Put the
# first PMIx 3 directory that really contains the library on LD_LIBRARY_PATH.
PMIX3_LIB=
for d in /cm/shared/apps/cm-pmix3/3.1.4/lib /cm/shared/apps/cm-pmix3/3.1.7/lib /cm/shared/apps/cm-pmix3/current/lib; do
    if [ -e "$d/libpmix.so.2" ]; then PMIX3_LIB=$d; break; fi
done
if [ -n "$PMIX3_LIB" ]; then
    export LD_LIBRARY_PATH=$PMIX3_LIB${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
else
    echo "dial_env.sh: WARNING: no PMIx 3 libpmix.so.2 found under /cm/shared/apps/cm-pmix3" >&2
fi

# Parallel HDF5 built by build_hdf5_dial.sh (DIaL3 has no parallel HDF5 module).
export HDF5_DIR=$HOME/sw/hdf5-1.14.6-ompi4-intel
export PATH=$HDF5_DIR/bin:$PATH
export LD_LIBRARY_PATH=$HDF5_DIR/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}

# Large automatic arrays (artvis etc. use ~1 GB per rank): Intel Fortran puts them on
# the stack, so lift the stack limit or the code segfaults.
ulimit -s unlimited
# Stop glibc returning large freed arrays to the OS every call (cost 57% of CPU time in
# the local 16384^2 run as page faults); freed memory is reused instead.
export MALLOC_MMAP_MAX_=0
export MALLOC_TRIM_THRESHOLD_=68719476736
