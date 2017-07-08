module alloc_q_omega_arrays_class
use common_block_class
use various_constants_class
use converters_class
use control_containers_class
use time_corr_algorithm_class

  implicit none
  
  private :: begin_element, end_element, pcdata_chunk
  private :: start_document, end_document                               


contains


  ! If control-object/q-values and control-object/omega-values exist in job file then populate 
  ! setup_mdmc_control_params%omega_values and setup_mdmc_control_params%q_values. 
  ! Done it this way rather than directly in handler because of convenience.

  subroutine alloc_q_omega_arrays(filename)
    use flib_sax  
    use flib_xpath    
    character(len=*), intent(in) :: filename

    type(xml_t) :: fxml
    integer     :: iostat
    
    character(len=40) :: read_db, read_int
    type (dictionary_t) :: structure_attributes
    
    real(db) :: r_first, r_last, r_max
    real(db) :: t_last
    real(db) :: start_value, step_size
    
    integer :: q_index, omega_index
    integer :: i
    

    
    integer :: num_q_values = 0  ! Count the number of q values before allocating q array
    integer :: num_omega_values = 0  ! Count the number of omega values before allocating omega array
    integer :: n_step
    
    
    call open_xmlfile(trim(filename),fxml,iostat)
    if (iostat /= 0) stop "Cannot open file which is suppose to contain q and omega values."
    
    
    ! count number of q-values
         
    do
      call get_node(fxml, path="//job/control-object/q-values/q", &
        attributes=structure_attributes, status=iostat)
      if (iostat < 0) exit
      
      call get_value(structure_attributes,"n-step",read_int,status=iostat)
      num_q_values = num_q_values + string_to_int(read_int) + 1
    end do
      
    call close_xmlfile(fxml)  
    call open_xmlfile(trim(filename),fxml,iostat)
    
    
    ! check for omega-values
                            
    do
      call get_node(fxml, path="//job/control-object/omega-values/omega", &
        attributes=structure_attributes, status=iostat)
      if (iostat < 0) exit
      
      call get_value(structure_attributes,"n-step",read_int,status=iostat)
      num_omega_values = num_omega_values + string_to_int(read_int) + 1
    end do
      
    call close_xmlfile(fxml)  
    call open_xmlfile(trim(filename),fxml,iostat)
    
    
    ! allocate q and omega arrays
    
    allocate(setup_mdmc_control_params%q_values(num_q_values))
    allocate(setup_mdmc_control_params%omega_values(num_omega_values))
    
    ! read in q-values
         
    q_index = 1
    do
      call get_node(fxml, path="//job/control-object/q-values/q", &
        attributes=structure_attributes, status=iostat)
      if (iostat < 0) exit
      
      call get_value(structure_attributes,"n-step",read_int,status=iostat)
      n_step = string_to_int(read_int)
      call get_value(structure_attributes,"start",read_db,status=iostat)
      start_value = string_to_db(read_db)      
      call get_value(structure_attributes,"step",read_db,status=iostat)
      step_size = string_to_db(read_db)        
      
      do i = 0, n_step
        setup_mdmc_control_params%q_values(q_index) = start_value + i*step_size
        q_index = q_index + 1
      end do
    end do
      
    call close_xmlfile(fxml)  
    call open_xmlfile(trim(filename),fxml,iostat)
    
    
    ! read in omega-values
         
    omega_index = 1
    do
      call get_node(fxml, path="//job/control-object/omega-values/omega", &
        attributes=structure_attributes, status=iostat)
      if (iostat < 0) exit
      
      call get_value(structure_attributes,"n-step",read_int,status=iostat)
      n_step = string_to_int(read_int)
      call get_value(structure_attributes,"start",read_db,status=iostat)
      start_value = string_to_db(read_db)      
      call get_value(structure_attributes,"step",read_db,status=iostat)
      step_size = string_to_db(read_db)        
      
      do i = 0, n_step
        setup_mdmc_control_params%omega_values(omega_index) = start_value + i*step_size
        omega_index = omega_index + 1
      end do
    end do
      
    call close_xmlfile(fxml)  


  end subroutine alloc_q_omega_arrays

  
  
  !START_DOCUMENT
  subroutine start_document()
    use flib_sax
  end subroutine start_document


  ! BEGIN_ELEMENT
  subroutine begin_element(name,attributes)
    use flib_sax  
    character(len=*), intent(in)   :: name
    type(dictionary_t), intent(in) :: attributes
    
    integer :: status
    real(db) :: number_db, number_db2
    character(len=40) :: read_db, read_int    
    

    select case(name)                       
      !case("G-d")
	
    end select


  end subroutine begin_element


  ! PCDATA_CHUNK
  subroutine pcdata_chunk(chunk)
    use flib_sax
    character(len=*), intent(in) :: chunk

  end subroutine pcdata_chunk


  ! END_ELEMENT
  subroutine end_element(name)
    use flib_sax  
    character(len=*), intent(in)   :: name

  end subroutine end_element


  !END_DOCUMENT
  subroutine end_document()
    use flib_sax  
  end subroutine end_document

end module alloc_q_omega_arrays_class
