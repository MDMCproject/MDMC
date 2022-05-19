OBJFILES= m_debug.o m_buffer.o m_elstack.o m_dictionary.o m_io.o \
          m_reader.o m_converters.o m_charset.o  m_entities.o \
          m_xml_error.o  m_fsm.o  m_xml_parser.o flib_sax.o
MODFILES=$(OBJFILES:.o=)

FLIB_ROOT=../macros
MK=$(FLIB_ROOT)/lib.mk
include $(MK)

clean:
	rm -f *.o  
	rm -f $(MOD_STD)*.$(MOD_EXT)
	rm -f $(LIB_STD)$(LIBRARY)

