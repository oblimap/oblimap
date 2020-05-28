#    Makefile for compiling the programs which are part of the of
#    the XXX package: an XXX model
#    Thomas Reerink

#    To execute the make file (in the src/ directory) e.g.:
#     make fortran_example_model_using_namelist_program  Notice: without the  .f90  extension
#     make                                               Short cut for the first program below
#     make all                                           Compiling all programs
#     make clean                                         Removing all object, module and executables files

#    For "make" syntax search on "make tutorial", or see:
#     http://www.gnu.org/software/make/manual/make.html

include Makefile.mpif90
#include Makefile.gfortran
#include Makefile.gfortran-mp-6-mac
#include Makefile.gfortran-lisa
#include Makefile.ifort-lisa
#include Makefile.ifort-staff.science

# The OBJ_PATH and the MOD_PATH are defined:
OBJ_PATH = object-files
MOD_PATH = module-files

# Pattern rules for creating object files from *.f90 (fortran 90) source files
$(OBJ_PATH)/%.o : %.f90
	$(F90) $(F90FLAGS) -c -o $@ $<

# This rule is used to build an executable from the object files:
% :
	$(F90) $(F90FLAGS)    -o $@ $^ $(LDFLAGS)



# ================================================================
#       Building the fortran_example_model_using_namelist_program
# ================================================================

# Adding this program to the list of all programs:
all_programs += fortran_mpi_win_allocate_shared_example_program

# List of object files for fortran_example_model_using_namelist_program:
fortran_mpi_win_allocate_shared_example_program_objects = \
        fortran_mpi_win_allocate_shared_example_program.o

# The executable (the target) depends on the object files:
fortran_mpi_win_allocate_shared_example_program : $(fortran_mpi_win_allocate_shared_example_program_objects:%=$(OBJ_PATH)/%)



# =============================================================
#       Cleaning up everything
# =============================================================

.PHONY : clean

clean :
	rm -fv $(all_programs) $(OBJ_PATH)/*.o $(MOD_PATH)/*.mod *.mod *.lst



# =============================================================
#       Building everything
# =============================================================

.PHONY : all

all : $(all_programs)
