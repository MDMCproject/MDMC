program moldynmoncar

  use flib_sax
  use handler_class

  type(xml_t) :: fxml
  integer     :: iostat


  call open_xmlfile("input/test.xml",fxml,iostat)
  if (iostat /= 0) stop "Cannot open file."

  call xml_parse(fxml,           &
                 begin_element_handler=begin_element, &
                 end_element_handler=end_element,     &
                 pcdata_chunk_handler=pcdata_chunk, &
                 start_document=start_document, &
                 end_document=end_document)

end program moldynmoncar


