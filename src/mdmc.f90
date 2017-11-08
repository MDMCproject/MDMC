program mdmc

  use flib_sax
  use handler_class, only :  startup_handler
  use handler_structure_part_class, only : handler_structure_part
  use structure_reader_class
  use flib_xpath
  use fom_readers_class
  use g_d_readers_class
  use s_qt_reader_class
  use s_qo_reader_class
  use rdf_class
  use alloc_q_omega_arrays_class

  type(xml_t) :: fxml
  integer     :: iostat

  character(len=99) :: filename, structure_filename
  type (dictionary_t) :: structure_attributes
  
  logical :: build_structure_from_model = .true.  
  
  ! set build in random generator function seeds
  ! note only need to do this once 
  !
  ! Note I can't seem to force random_number to always
  ! produce the same sequence of random numbers...
  
  integer :: n_seeds  
  logical :: output_dir_exist
  INTEGER, ALLOCATABLE :: new (:)
  call random_seed(size=n_seeds)
  ALLOCATE (new(I))
  new = 5
  CALL RANDOM_SEED (PUT=new(1:I))
  
  ! Ask user for a MDMC job file
  
  write (*,*) "Enter the filename of a MDMC job file: "
  read *, filename 
  
  ! create output directory - at this point in time hardcoded to 'output'
  
  inquire( directory="output", exist=output_dir_exist )
  if ( output_dir_exist == .false. ) then
    write(*,*), "create output directory"
    call system('mkdir output')
  end if    
  
  ! Because of bug in otherwise very useful xmlf90 library 
  ! need to read in input file in bits...... 
  
  ! Assume <structure> element is before the <fom> element
  ! and make structure differently depending on whether the
  ! <structure> element has an attribute "filename" or not
  
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
  
  ! Same as above but for the <fom> element. 
  
  call open_xmlfile(trim(filename),fxml,iostat)
  call get_node(fxml, path="//fom/rdf-fom/rdf-data",attributes=structure_attributes,status=iostat)
  call get_value(structure_attributes, "filename", structure_filename, iostat)
  call close_xmlfile(fxml)  
  
  if (iostat == 0) then
    call make_rdf_fom(trim(structure_filename))
  end if
  
  ! Similar to above but for the <control-object/q-values> and <control-object/omega-values> elements.
  ! As of this writing there is not an option for specifying q and omega values from a data file
  
  call alloc_q_omega_arrays(trim(filename))

  ! open main input file to check for g-d-rt-fom element
  ! otherwise this code is here for the same reason as above
  
  call open_xmlfile(trim(filename),fxml,iostat)
  call get_node(fxml, path="//fom/g-d-rt-fom/data-file",attributes=structure_attributes,status=iostat)
  call get_value(structure_attributes, "filename", structure_filename, iostat)
  call close_xmlfile(fxml)  
  
  if (iostat == 0) then
    call make_g_d_rt_fom_container(trim(structure_filename))
  end if  
    
  ! open main input file to check for s-qt-fom element
  ! otherwise this code is here for the same reason as above
  
  call open_xmlfile(trim(filename),fxml,iostat)
  call get_node(fxml, path="//fom/s-qt-fom/data-file",attributes=structure_attributes,status=iostat)
  call get_value(structure_attributes, "filename", structure_filename, iostat)
  call close_xmlfile(fxml)  
  
  if (iostat == 0) then
    call make_s_qt_fom_container(trim(structure_filename))
  end if    
  
  ! open main input file to check for s-qo-fom element
  ! otherwise this code is here for the same reason as above
  
  call open_xmlfile(trim(filename),fxml,iostat)
  call get_node(fxml, path="//fom/s-qo-fom/data-file",attributes=structure_attributes,status=iostat)
  call get_value(structure_attributes, "filename", structure_filename, iostat)
  call close_xmlfile(fxml)  
  
  if (iostat == 0) then
    call make_s_qo_fom_container(trim(structure_filename))
  end if    

  ! Now finally read in the rest of the job file and run the 
  ! simulation specified in the job file 
  
  call startup_handler(trim(filename))
  
end program mdmc


