program mdmc

  use flib_sax
  use handler_class, only :  startup_handler
  use handler_structure_part_class, only : handler_structure_part
  use structure_reader_class
  use flib_xpath
  use fom_readers_class
  use rdf_class

  type(xml_t) :: fxml
  integer     :: iostat

  character(len=99) :: filename, structure_filename
  type (dictionary_t) :: structure_attributes
  
  logical :: build_structure_from_model = .true.
  
  
  write (*,*) "Enter pathname/filename"
  read *, filename 
    
  ! Because of bug in otherwise very useful xmlf90 library 
  ! need to read in input file in bits...... 
  
  call open_xmlfile(trim(filename),fxml,iostat)
  if (iostat /= 0) stop "Cannot open file."
  
  call get_node(fxml, path="//structure",attributes=structure_attributes,status=iostat)
  
  call get_value(structure_attributes, "filename", structure_filename, iostat)
  if (iostat == 0) then
    build_structure_from_model = .false.
    call make_structure(trim(structure_filename))
  end if
   
  call close_xmlfile(fxml)
  
  
  if (build_structure_from_model) then
    call handler_structure_part(trim(filename))
  end if
  
  
  call open_xmlfile(trim(filename),fxml,iostat)
  
  call get_node(fxml, path="//fom/rdf-fom/rdf-data",attributes=structure_attributes,status=iostat)
  
  call get_value(structure_attributes, "filename", structure_filename, iostat)
  if (iostat == 0) then
    call make_rdf_fom(trim(structure_filename))
  end if
  
  call close_xmlfile(fxml)  
  
  !call save_rdf(target_rdf_fom%rdf_data) 
  !stop
  
  call startup_handler(trim(filename))
  

end program mdmc


