#
# Fortran macros for Intel Fortran compiler
#
.SUFFIXES: .SUFFIXES .f .f90 .f95 .F .F90 .F95 .fpp
#
LIB_STD=$(FLIB_ROOT)/lib/
MOD_STD=$(FLIB_ROOT)/modules/
INC_STD=$(FLIB_ROOT)/include/
#BIN_STD=$(FLIB_ROOT)/bin
#
FC=ifort
# For some reason you have to swich off all optimisation to get ifort
# to link properly using Fedora Core 4!!!!
# I copied and pasted to below flags from a suggestion by a bloke 
# called Andrew from some Intel mailing list
CFLAGS= -pthread -pad -Vaxlib -mp -save -O0 -w95 -extend_source
FFLAGS=     $(CFLAGS)
#FFLAGS_DEBUG= -g  
#FFLAGS_CHECK= -g 
LDFLAGS=      $(CFLAGS)
#
INC_PREFIX=-I
MOD_PREFIX=-I
LIB_PREFIX=-L
#
MOD_EXT=mod
MOD_SEARCH_STD= $(MOD_PREFIX)$(MOD_STD) 
MOD_SEARCH= $(MOD_SEARCH_STD) $(MOD_SEARCH_OTHER)
#INC_SEARCH= $(INC_PREFIX)$(INC_STD)
#
#
AR=ar
RANLIB=echo
#
#CPP=/bin/cpp -P
#COCO=$(BIN_STD)coco -I $(INC_STD)
#DEFS=
#
# Experimental : the following deactivates an implicit rule
# which breaks havoc with the operation of this makefile
# It works at least with GNU make
%.o : %.mod
#
#.F.o:
#	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS)  $(DEFS) $<
.f.o:
	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS)   $<
#.F90.o:
#	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS) $(DEFS) $<
.f90.o:
	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS)   $<
.F95.o:
	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS)   $<
#.f95.o:
#	cp $< $*.f90
#	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS)   $*.f90
#	rm -f $*.f90
#
#.fpp.f90:
#	$(COCO) $*
#









