FCOM=ifort
LINK=ifort

FLIB_ROOT=xmlf90-1.2/macros
include $(FLIB_ROOT)/fortran.mk

#### link to Fortran XML library (xmlf90) ####

LIBS=-L$(FLIB_ROOT)/lib/ -lxmlf90_sax -lxmlf90_wxml
MOD_SEARCH=-I$(FLIB_ROOT)/modules/


OBJS= structure.o fnc_constraint.o configuration.o  \
      common_block.o structure_reader.o fnc_constraint_reader.o \
      md_control.o handler.o mdmc.o 

TARGET=mdmc


default: $(TARGET)


#.f.o:
#	$(FCOM) -c $(FOPT) $*.f
#%.o: %.f90
#	$(FCOM) $(MOD_SEARCH) -c $(FOPT) $*.f90


$(TARGET): $(OBJS)
	$(LINK) $(FFLAGS) -o $(TARGET) $(OBJS) $(LIBS)


clean: 
	rm -f *.o
	rm -f *.mod
	rm -f $(TARGET)










