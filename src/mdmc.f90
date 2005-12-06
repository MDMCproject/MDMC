program mdmc

  use flib_sax
  use handler_class
  use structure_reader_class
  use flib_xpath

  type(xml_t) :: fxml
  integer     :: iostat

  character(len=99) :: filename, structure_filename
  type (dictionary_t) :: structure_attributes
  
  write (*,*) "Enter pathname/filename"
  read *, filename 
    
  
  call open_xmlfile(trim(filename),fxml,iostat)
  if (iostat /= 0) stop "Cannot open file."
  
  call get_node(fxml, path="//structure",attributes=structure_attributes,status=iostat)
  
  call get_value(structure_attributes, "filename", structure_filename, iostat)
  if (iostat == 0) then
    call make_structure(trim(structure_filename))
  end if
  
  call close_xmlfile(fxml)
  
  
  call startup_handler(trim(filename))
  

end program mdmc


