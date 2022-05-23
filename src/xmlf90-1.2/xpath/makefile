OBJFILES= m_path.o flib_xpath.o
MODFILES=$(OBJFILES:.o=)

LIBRARY=libxmlf90_xpath.a
FLIB_ROOT=../macros
LIBS=-L$(FLIB_ROOT)/lib/ -lxmlf90_sax
MK=$(FLIB_ROOT)/lib.mk
include $(MK)

clean:
	rm -f *.o  
	rm -f $(MOD_STD)*.$(MOD_EXT)
	rm -f $(LIB_STD)$(LIBRARY)
