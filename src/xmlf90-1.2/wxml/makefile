OBJFILES= m_wxml_elstack.o m_wxml_dictionary.o m_wxml_buffer.o m_wxml_core.o \
          m_wxml_text.o flib_wxml.o
MODFILES=$(OBJFILES:.o=)

FLIB_ROOT=../macros
MK=$(FLIB_ROOT)/lib.mk
include $(MK)

clean:
	rm -f *.o  
	rm -f $(MOD_STD)*.$(MOD_EXT)
	rm -f $(LIB_STD)$(LIBRARY)

