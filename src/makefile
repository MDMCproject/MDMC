FLIB_ROOT=xmlf90-1.2/macros
include $(FLIB_ROOT)/fortran.mk

#### link to Fortran XML library (xmlf90) ####

LIBS=-L$(FLIB_ROOT)/lib/ -lxmlf90_sax -lxmlf90_wxml -lxmlf90_xpath
MOD_SEARCH=-I$(FLIB_ROOT)/modules/
#INC_SEARCH=-I./functions -I./algorithms -I./readers -I.
F90_FUNCTIONS=$(wildcard functions/*.f90)
OBJ_FUNCTIONS=$(patsubst %.f90,%.o,$(F90_FUNCTIONS))
F90_READERS=$(wildcard readers/*.f90)
OBJ_READERS=$(patsubst %.f90,%.o,$(F90_READERS))
F90_ALGORITHMS=$(wildcard algorithms/*.f90)
OBJ_ALGORITHMS=$(patsubst %.f90,%.o,$(F90_ALGORITHMS))


OBJS= various_constants.o structure_type_definition.o \
      structure_nn_methods.o tic_toc.o converters.o \
      structure.o fnc_constraint.o configuration.o \
      func_param.o \
      histogram.o rdf.o histogram_cutdown.o time_corr_hist_container.o s_q_time.o s_q_omega.o \
      $(OBJ_FUNCTIONS) function.o \
      common_block.o phasespace.o time_corr_algorithm.o fnc_constraint.o \
      configuration.o md_properties.o control_containers.o \
      $(OBJ_READERS) \
      func_params_wrapper.o $(OBJ_ALGORITHMS) \
      handler.o \
      fom_readers.o handler_structure_part.o mdmc.o

TARGET=mdmc


default: $(TARGET)

#%.o: %.f90
#	$(FC) $(MOD_SEARCH) $(FFLAGS) -o $@ -c $<


$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $(TARGET) $(OBJS) $(LIBS)


clean: 
	rm -f *.o
	rm -f *.mod
	rm -f $(TARGET)
	rm -f functions/*.o
	rm -f algorithms/*.o
	rm -f readers/*o
