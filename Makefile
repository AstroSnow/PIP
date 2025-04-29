TARGET = a.out
OBJECTS = Util_rot.o BMP_rot.o Matrix_rot.o \
	Globalvar.o Parameters.o MPI_rot.o Boundary_rot.o Scheme_rot.o \
	HC_rot.o \
	Res_rot.o Gra_rot.o PIP_rot.o Visc_rot.o\
	IOT_rot.o IO_rot.o   \
	HLL_rot.o Model_rot.o \
	Solver_rot.o\
	set_coordinate_NUG.o set_fold_grid.o \
	linear_wave.o  \
	shock_tube.o explosion.o currentsheet.o KH.o RT.o\
	Orszag_Tang.o FieldLoop.o ambipolar.o PNcoupling.o\
	HCtest.o HCtest_Tonly.o Flare.o Coalescence.o Sedov.o tearing.o \
	NUwave.o Reconnection.o BMP_density.o blob.o field_diffusion.o\
	CS_collapse.o asym_currentsheet.o CSC.o\
	Alfven_Damping.o cnd_tube.o shock_tube_ion.o\
	MRI.o disk_flare.o mass_load_prom.o KHnlev.o KHtube.o\
	Initial_rot.o  main.o relax_prom.o relax_prom2.o procedures.o hsstatic.o\
	resonator.o ionrectest.o  Complete_spectrum.o shock_tube_stab.o shock_tube_stab2.o \
	shock_tube_stab3.o kink_wave.o RMI.o kink_instability.o jetcrosssection.o HC_tests.o
MOD_FILES = util_rot.mod globalvar.mod parameters.mod \
	scheme_rot.mod hc_rot.mod \
	res_rot.mod gra_rot.mod pip_rot.mod visc_rot.mod\
	mpi_rot.mod io_rot.mod \
	hll_rot.mod model_rot.mod  \
	boundary_rot.mod solver_rot.mod matrix_rot.mod initial_rot.mod  procedures.mod

#FC = gfortran
FC = h5pfc
#FC = mpif90-mpich-mp -O2
#FC = mpif90-mpich-mp


FFLAGS = -O2
#FFLAGS = 
LDFLAGS =
LIB_DIR=.
#DEBUG = -g -pg
DEBUG= -ffpe-trap=invalid,zero,overflow -fbacktrace -fbounds-check -g
#DEBUG= -ffpe-trap=invalid,overflow -fbacktrace -fbounds-check -g
#DEBUG= -fbacktrace -fbounds-check -g
#DEBUG= -fbacktrace -g
DEBUG=  #-fallow-argument-mismatch -ffpe-trap=invalid,zero,overflow -fbacktrace -fbounds-check -g

.SUFFIXES : .o .f90
.f90.o:
	${FC} ${FFLAGS} ${DEBUG} -c $<


${TARGET} : ${OBJECTS}
	${FC} ${DEBUG} -L$(LIB_DIR) ${LDFLAGS} -o $@ ${OBJECTS}
#	rm ${OBJECTS} ${MOD_FILES}
clean:
	rm ${TARGET} ${OBJECTS} ${MOD_FILES}

datatidy:
	rm Data/*

print_vars:
	@echo "FC='${FC}'"
