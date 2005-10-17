#
# Macro file to compile and install modules in FLIB itself
#
include $(FLIB_ROOT)/fortran.mk
#
LIBRARY=libflib.a
#
CP=cp
install: $(OBJFILES)
	@echo "  ==> Updating $(LIBRARY) with $(OBJFILES)"
	$(AR) r $(LIB_STD)$(LIBRARY) $(OBJFILES)
	$(RANLIB) $(LIB_STD)$(LIBRARY)
	@echo "  ==> Installing modules: $(MODFILES)"
	@for i in $(MODFILES) ; do  \
	   if [ -f $$i.$(MOD_EXT) ] ; then \
           $(CP) $$i.$(MOD_EXT) $(MOD_STD) ; else  \
           $(CP) `echo $$i | tr [a-z] [A-Z]`.$(MOD_EXT) $(MOD_STD) ; fi ; done
#
# The convoluted logic above is there to support Intel Fortran's 
# convention of producing .mod files with uppercase names. 








