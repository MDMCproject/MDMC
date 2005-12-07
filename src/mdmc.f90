program mdmc

  use flib_sax
  use handler_class
  use structure_reader_class
  use flib_xpath
  use fom_readers_class
  use rdf_class

  type(xml_t) :: fxml
  integer     :: iostat

  character(len=99) :: filename, structure_filename
  type (dictionary_t) :: structure_attributes
  
  write (*,*) "Enter pathname/filename"
  read *, filename 
    
  
  ! Because of bug in otherwise very useful xmlf90 library 
  ! need to read in input file in bits...... 
  
  call open_xmlfile(trim(filename),fxml,iostat)
  if (iostat /= 0) stop "Cannot open file."
  
  call get_node(fxml, path="//structure",attributes=structure_attributes,status=iostat)
  
  call get_value(structure_attributes, "filename", structure_filename, iostat)
  if (iostat == 0) then
    call make_structure(trim(structure_filename))
  !else
  !  call get_value(structure_attributes, "what-init-structure-to-build", &
  !                 structure_filename, iostat)
                   
                   
  end if
   
  
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


