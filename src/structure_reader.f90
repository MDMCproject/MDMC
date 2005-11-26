module structure_reader_class
use flib_sax
use common_block_class, only : common_config
use various_constants_class


  implicit none
  
  private :: begin_element, end_element, pcdata_chunk
  private :: start_document, end_document
  
  integer, private :: count_number_atoms

contains

  subroutine make_structure(filename)
    character(len=*), intent(in) :: filename

    type(xml_t) :: fxml
    integer     :: iostat


    call open_xmlfile(filename,fxml,iostat)
    if (iostat /= 0) stop "Cannot open file."

    call xml_parse(fxml,           &
                 begin_element_handler=begin_element, &
                 end_element_handler=end_element,     &
                 pcdata_chunk_handler=pcdata_chunk, &
                 start_document=start_document, &
                 end_document=end_document)

  end subroutine make_structure

  
  subroutine make_simple_cubic_structure(density, n_atoms)
    real(db), intent(in) :: density
    integer, dimension(ndim) :: n_atoms
    
    real(db), dimension(ndim) :: edges
    real(db), dimension(ndim) :: gap, put_atom_at ! gap between atoms
    
    integer :: n = 1, nx, ny, nz, n_tot
    
    
    ! calculate box edge and gap

    edges = n_atoms / density**(1.0_db/3.0_db)
    common_config%str%box_edges = edges
    gap = edges / n_atoms
    
    
    n_tot = product(n_atoms)
    allocate(common_config%str%atoms(n_tot))   ! allocate mass, name...
    allocate(common_config%str%r(n_tot,ndim))  ! allocate coordinates
    
    do nz = 1, n_atoms(3)
      do ny = 1, n_atoms(2)
        do nx = 1, n_atoms(1)
          put_atom_at(1) = nx - 0.5_db
          put_atom_at(2) = ny - 0.5_db
          put_atom_at(3) = nz - 0.5_db
          put_atom_at = put_atom_at * gap - 0.5_db * edges

          common_config%str%r(n,:) = put_atom_at
          n = n + 1
        end do
      end do
    end do
    
  end subroutine make_simple_cubic_structure
  
  
  !START_DOCUMENT
  subroutine start_document()
  end subroutine start_document


  ! BEGIN_ELEMENT
  subroutine begin_element(name,attributes)
    character(len=*), intent(in)   :: name
    type(dictionary_t), intent(in) :: attributes
    
    integer :: status, ndata=0
    real(db) :: number_db(1)
    character(len=40) :: read_db
    character(len=40) :: control_object_name, read_number_atoms
    integer :: number_atoms(1);
    character(len=120) :: filename

    select case(name)
      case("molecule")
        call get_value(attributes,"number", read_number_atoms,status)
        call build_data_array(read_number_atoms, number_atoms, ndata)
        write(*,*) number_atoms
        allocate(common_config%str%atoms(number_atoms(1)))  ! allocate mass, name...
        allocate(common_config%str%r(number_atoms(1), ndim))  ! allocate coordinates
        count_number_atoms = 1

      case("atom")
        call get_value(attributes,"x3", read_db,status)
        ndata = 0
        call build_data_array(read_db, number_db, ndata)
        common_config%str%r(count_number_atoms,1) = number_db(1)

        call get_value(attributes,"y3", read_db,status)
        ndata = 0
        call build_data_array(read_db, number_db, ndata)
        common_config%str%r(count_number_atoms,2) = number_db(1)

        call get_value(attributes,"z3", read_db,status)
        ndata = 0
        call build_data_array(read_db, number_db, ndata)
        common_config%str%r(count_number_atoms,3) = number_db(1)

        write (*,*) count_number_atoms, &
          common_config%str%r(count_number_atoms,1), &
          common_config%str%r(count_number_atoms,2), &
          common_config%str%r(count_number_atoms,3)
        count_number_atoms = count_number_atoms + 1        

	
    end select


  end subroutine begin_element


  ! PCDATA_CHUNK
  subroutine pcdata_chunk(chunk)
    character(len=*), intent(in) :: chunk

  end subroutine pcdata_chunk


  ! END_ELEMENT
  subroutine end_element(name)
    character(len=*), intent(in)   :: name



  end subroutine end_element


  !END_DOCUMENT
  subroutine end_document()
  end subroutine end_document

end module structure_reader_class
