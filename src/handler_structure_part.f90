module handler_structure_part_class
use flib_sax
use structure_reader_class

  implicit none
  
  
  private :: begin_element, end_element, pcdata_chunk
  private :: start_document, end_document
  
  public :: handler_structure_part
  

  logical, private  :: in_structure = .false. 
  
  ! both to be passed to make_simple_cubic_structure or make_fcc_structure
  real(db), private :: density   
  integer, dimension(ndim), private :: unitcells_xyz
                                             
  
  character(len=99), private :: what_init_structure_to_build
contains
  
  subroutine handler_structure_part(filename)
    character(len=*), intent(in) :: filename
  
    type(xml_t) :: fxml
    integer     :: iostat
  
    call open_xmlfile(trim(filename),fxml,iostat)
    if (iostat /= 0) stop "Cannot open file."

    call xml_parse(fxml,           &
                 begin_element_handler=begin_element, &
                 end_element_handler=end_element,     &
                 pcdata_chunk_handler=pcdata_chunk, &
                 start_document=start_document, &
                 end_document=end_document)
  
    call close_xmlfile(fxml)
  
  end subroutine handler_structure_part


  ! START_DOCUMENT
  subroutine start_document()
  end subroutine start_document


  ! BEGIN_ELEMENT
  subroutine begin_element(name,attributes)
    character(len=*), intent(in)   :: name
    type(dictionary_t), intent(in) :: attributes
    
    integer :: status, n
    integer :: number_int
    character(len=40) :: read_db, read_int
    character(len=40) :: units
    character(len=120) :: filename

   ! write(*,*) name
    select case(trim(name))
      case("structure")
        call get_value(attributes,"filename",filename,status)
        if (status == 0) then
          print *, "ERROR: not allowed to read in structure"
          print *, "from handler_structure_part"
          stop
        else
          in_structure = .true.
          
          
          ! which initial structure to create
          
          call get_value(attributes,"what-init-structure-to-build", &
                         what_init_structure_to_build, status)
          if (status /= 0) then
            write (*,*) "ERROR Initial structure type not specified in input file"
            stop
          end if
            
        end if
	
    end select


    ! if in structure
    
    if (in_structure) then
      select case(name)
        case("density")          
          ! read density value into read_dp(1)
          call get_value(attributes,"val",read_db,status)
          
          ! read units
          call get_value(attributes,"units",units,status)

          if (units == "atom/AA3") then
            density = string_to_db(read_db)
          else
            write(*,*) "ERROR, unit of density currently must be in atom/AA3"
            stop            
          end if
          
        case("number-of-unit-cells")
          ! check that the number of attributes of this element equals ndim
          number_int = number_of_entries(attributes)
          if (number_int /= ndim) then
            write(*,'(a,i2)') "ndim ", ndim
            write(*,'(a,i2)') "number of ar attributes ", number_int
            write(*,*) "ERROR, ndim not equal number of ar attr"
            stop
          end if
          
          ! read number of atoms into read_int(1)
          call get_value(attributes,"nx",read_int,status)         
          unitcells_xyz(1) = string_to_int(read_int)
          if (ndim > 1) then
            call get_value(attributes,"ny",read_int,status)         
            unitcells_xyz(2) = string_to_int(read_int)
          end if
          
          if (ndim > 2) then
            call get_value(attributes,"nz",read_int,status)         
            unitcells_xyz(3) = string_to_int(read_int)
          end if
          
      end select
    end if
    
  end subroutine begin_element


  ! PCDATA_CHUNK
  subroutine pcdata_chunk(chunk)
    character(len=*), intent(in) :: chunk

  end subroutine pcdata_chunk


  ! END_ELEMENT
  subroutine end_element(name)
    character(len=*), intent(in)   :: name
  
    if (name == "structure") then
      select case(trim(what_init_structure_to_build))
        case("simple-cubic")
          call make_simple_cubic_structure(density, unitcells_xyz)
          
        case("fcc")
          call make_fcc_structure(density, unitcells_xyz)
          
        case default
          write (*,*) "ERROR initial what-init-structure-to-build attribute"
          write (*,*) "in input file not recognized"
          write (*,'(2a)') "what-init-structure-to-build = ", what_init_structure_to_build
          stop         
      end select
    end if

  end subroutine end_element


  !END_DOCUMENT
  subroutine end_document()
  end subroutine end_document

end module handler_structure_part_class
