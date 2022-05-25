#
# Fortran macros for Intel Fortran compiler
#
.SUFFIXES: .SUFFIXES .f .f90 .f95 .F .F90 .F95 .fpp
#
MOD_STD=$(FLIB_ROOT)/modules/
#
FC=ifort
# For some reason you have to swich off all optimisation to get ifort
# to link properly using Fedora Core 4!!!!
# I copied and pasted to below flags from a suggestion by a bloke
# called Andrew from some Intel mailing list
FFLAGS= -pthread -pad -save -O0 -extend-source

INC_PREFIX=-I
MOD_PREFIX=-I
LIB_PREFIX=-L
#
MOD_EXT=mod
MOD_SEARCH=$(MOD_PREFIX)$(MOD_STD)
#
# Experimental : the following deactivates an implicit rule
# which breaks havoc with the operation of this makefile
# It works at least with GNU make
%.o : %.mod
#
#.f90.o:
#	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS)   $<
%.o: %.f90
	$(FC) $(MOD_SEARCH) $(FFLAGS) -o $@ -c $<












