#!/bin/bash
# Build parallel HDF5 1.14.6 with Fortran on DIaL3, against openmpi4/intel/4.0.5.
# Run once on a login node, from the PIP_OT directory:   bash build_hdf5_dial.sh
# Takes about 20-30 minutes. Installs into $HDF5_DIR (set in dial_env.sh).

set -eo pipefail
source "$(dirname "$0")/dial_env.sh"
set -u

VER=1.14.6
SRC=$HOME/sw/src

echo "== Checking the compilers"
if ! mpicc --version >/dev/null 2>&1 || ! mpif90 --version >/dev/null 2>&1; then
    echo "ERROR: mpicc/mpif90 do not work. The Intel compiler module in dial_env.sh is probably" >&2
    echo "       wrong; edit its 'module load intel-...' line and try again." >&2
    exit 1
fi
mpif90 --version | head -1
which mpicc mpif90

if [ -x "$HDF5_DIR/bin/h5pfc" ]; then
    echo "HDF5 already installed at $HDF5_DIR; delete it to rebuild."
    exit 0
fi

mkdir -p "$SRC"
cd "$SRC"
TARBALL=hdf5-$VER.tar.gz
if [ ! -f "$TARBALL" ]; then
    echo "== Downloading HDF5 $VER"
    wget -O "$TARBALL" "https://github.com/HDFGroup/hdf5/releases/download/hdf5_$VER/$TARBALL" || {
        rm -f "$TARBALL"
        echo "ERROR: download failed. Download $TARBALL from https://www.hdfgroup.org/download-hdf5/" >&2
        echo "       and copy it to $SRC, then rerun this script." >&2
        exit 1
    }
fi

rm -rf "hdf5-$VER"
tar xzf "$TARBALL"
cd "hdf5-$VER"

echo "== Configuring"
CC=mpicc FC=mpif90 CFLAGS="-O2" FCFLAGS="-O2" \
    ./configure --prefix="$HDF5_DIR" --enable-parallel --enable-fortran \
                --disable-tests --disable-tools-tests 2>&1 | tail -30

echo "== Building"
make -j 16 > build.log 2>&1 || { tail -40 build.log; echo "ERROR: build failed, see $PWD/build.log" >&2; exit 1; }
make install > install.log 2>&1

echo "== Checking the install"
"$HDF5_DIR/bin/h5pfc" -showconfig | grep -E "HDF5 Version|Parallel HDF5|Fortran:|Fortran Compiler"
echo "Done. HDF5 is in $HDF5_DIR"
