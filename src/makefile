#
# Makefile for moldynmoncar
#

FCOM=ifort
LINK=ifort
FOPT=


#### link to Fortran XML library (xmlf90) ####

FLIB_ROOT=$(HOME)/xmlf90-1.2/macros
LIBS=-L$(FLIB_ROOT)/lib/ -lflib
MOD_SEARCH=-I$(FLIB_ROOT)/modules/

 
OBJS= structure.o fnc_constraint.o configuration.o  \
      common_block.o structure_reader.o fnc_constraint_reader.o \
      moldyn_control.o handler.o mdmc.o 

TARGET=mdmc


all: $(TARGET)


.f.o:
	$(FCOM) -c $(FOPT) $*.f
%.o: %.f90
	$(FCOM) $(MOD_SEARCH) -c $(FOPT) $*.f90
	

$(TARGET): $(OBJS)
	$(LINK) -o $(TARGET) $(OBJS) $(LIBS)


clean: 
	rm -f *.o
	rm -f *.mod
	rm -f $(TARGET)










