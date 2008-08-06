! Please read Umbrello for how this code perhaps should be 
! changed in the future

module g_d_readers_class
use common_block_class
use various_constants_class
use converters_class
use control_containers_class
use time_corr_algorithm_class

  implicit none
  
  private :: begin_element, end_element, pcdata_chunk
  private :: start_document, end_document                               


contains


  ! read G_d(r,t) data from file and store in target_g_d_rt_fom

  subroutine make_g_d_rt_fom_container(filename)
    use flib_sax  
    use flib_xpath    
    character(len=*), intent(in) :: filename

    type(xml_t) :: fxml
    integer     :: iostat
    
    character(len=40) :: read_db, read_int
    type (dictionary_t) :: structure_attributes
    
    real(db) :: r_first, r_last, r_max
    real(db) :: t_last
    real(db) :: r_bin, t_bin
    
    integer :: n_r_bin, n_time_bin
    integer :: i, j
    
    integer :: number_of_g_d  ! counts # of G-d elements
    
    ! before reading the rdf data file from a-to-z determine
    ! separately how many g-d elements it contains
    
    call open_xmlfile(trim(filename),fxml,iostat)
    if (iostat /= 0) stop "Cannot open g_d data file."
    
    
    ! get r_bin, density and number of atoms
    
    call get_node(fxml, path="//G_d-space-time-pair-correlation-function", attributes=structure_attributes, status=iostat)       
    call get_value(structure_attributes,"bin-length",read_db,status=iostat)
    if (iostat /= 0) stop "No bin-length attribute in g_d data file."
    r_bin = string_to_db(read_db)   
    
    !print *, "r_bin = ", r_bin
    
    call get_value(structure_attributes,"density",read_db,status=iostat)
    if (iostat /= 0) stop "No density attribute in g_d data file."
    target_g_d_rt_fom%density = string_to_db(read_db)  
    
    call get_value(structure_attributes,"n-atom",read_int,status=iostat)
    if (iostat /= 0) stop "No n-atom attribute in g_d data file."
    target_g_d_rt_fom%n_atom = string_to_int(read_int)
        
    !print *, target_g_d_rt_fom%n_atom   
        
    !stop    
        
     
    
    ! count first g_d element
    
    number_of_g_d = 1
    
    call get_node(fxml, path="//G_d-space-time-pair-correlation-function/G-d", attributes=structure_attributes, status=iostat)    
    
    call get_value(structure_attributes,"r",read_db,status=iostat)
    r_first = string_to_db(read_db)
    
    
    ! count remainder g_d element
    
    do 
      call get_node(fxml, path="//G_d-space-time-pair-correlation-function/G-d", &
        attributes=structure_attributes, status=iostat)
      if (iostat < 0) exit
      
      call get_value(structure_attributes,"r",read_db,status=iostat)
      r_last = string_to_db(read_db)
      
      call get_value(structure_attributes,"t",read_db,status=iostat)
      t_last = string_to_db(read_db)            
      
      number_of_g_d = number_of_g_d + 1
    end do


    r_max = r_last+r_first
    n_r_bin = nint( (r_last+r_first) / r_bin )
    n_time_bin = number_of_g_d / n_r_bin
    t_bin = t_last / (n_time_bin-1)
    
    
    target_g_d_rt_fom%n_r_bin = n_r_bin
    target_g_d_rt_fom%r_bin = r_bin
    target_g_d_rt_fom%n_t_bin = n_time_bin
    target_g_d_rt_fom%t_bin = t_bin


    allocate(target_g_d_rt_fom%obs(floor(r_max/r_bin), n_time_bin))
  
  
    call close_xmlfile(fxml)  
    

    call open_xmlfile(trim(filename),fxml,iostat)
    
    
    do j = 1, n_time_bin
      do i = 1, n_r_bin
      
        call get_node(fxml, path="//G_d-space-time-pair-correlation-function/G-d", &
          attributes=structure_attributes, status=iostat) 
        call get_value(structure_attributes,"G",read_db,status=iostat)
        target_g_d_rt_fom%obs(i, j) = string_to_db(read_db)         
        
        !print *, j, i 
      end do
    end do

!stop

!    ! now begin reading rdf data file from a-to-z
!
!    call open_xmlfile(filename,fxml,iostat)
!    if (iostat /= 0) stop "Cannot open file."
!
!    call xml_parse(fxml,           &
!                 begin_element_handler=begin_element, &
!                 end_element_handler=end_element,     &
!                 pcdata_chunk_handler=pcdata_chunk, &
!                 start_document=start_document, &
!                 end_document=end_document)

    
    call close_xmlfile(fxml)                             

  end subroutine make_g_d_rt_fom_container

  
  
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

end module g_d_readers_class
