module structure_reader_class
use common_block_class, only : common_config
use various_constants_class
use converters_class


  implicit none
  
  private :: begin_element, end_element, pcdata_chunk
  private :: start_document, end_document
  
  integer, private :: count_number_atoms = 1

contains

  subroutine make_structure(filename)
    use flib_sax  
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

    
    call close_xmlfile(fxml)                             

  end subroutine make_structure


  ! returns a structure which is locally simple cubic. However, this does not
  ! imply that the resulting box containing the atoms is cubic - this
  ! is only the case when unitcells_xyz(1:3) contains three identical integers
  ! for the number of unit cells in the x,y and z directions. 
  ! Further the resulting structure is returned with its center at the origin.
  
  subroutine make_simple_cubic_structure(density, unitcells_xyz)
    real(db), intent(in) :: density
    integer, dimension(ndim) :: unitcells_xyz
    
    real(db), dimension(ndim) :: edges
    real(db), dimension(ndim) :: gap, put_atom_at ! gap between atoms
    
    integer :: n = 1, nx, ny, nz, n_tot
    
    
    ! calculate box edges and gap between atoms

    edges = unitcells_xyz / density**(1.0_db/3.0_db)
    common_config%str%box_edges = edges
    
    
    ! make the gap between atoms the same in x, y, and z. This means the
    ! returning structure is locally simple cubic
         
    gap = edges / unitcells_xyz       
    
    
    n_tot = product(unitcells_xyz)
    allocate(common_config%str%atoms(n_tot))   ! allocate mass, name...
    allocate(common_config%str%r(n_tot,ndim))  ! allocate coordinates
    
    do nz = 1, unitcells_xyz(3)
      do ny = 1, unitcells_xyz(2)
        do nx = 1, unitcells_xyz(1)
          put_atom_at(1) = nx - 0.5_db
          put_atom_at(2) = ny - 0.5_db
          put_atom_at(3) = nz - 0.5_db
          put_atom_at = put_atom_at * gap - 0.5_db * edges

          common_config%str%r(n,:) = put_atom_at
          n = n + 1
        end do
      end do
    end do
    
    
    ! print to screen
    
    write(*, *) " "
    print *, 'Simple cubic structure created: '
    write(*,'(a,3f12.6)') "   Box size = ", edges
    print *, '  total number of atoms = ', product(unitcells_xyz)  
    write(*, *) " "
        
  end subroutine make_simple_cubic_structure
  
  
  ! returns a structure which is locally fcc. However, this does not
  ! imply that the resulting box containing the atoms is cubic - this
  ! is only the case when unitcells_xyz(1:3) contains three identical integers
  ! for the number of unit cells in the x,y and z directions. 
  ! Notice the total number of atoms in the returning structure is
  ! 4 * product(unitcells_xyz). 
  ! Further the resulting structure is returned with its center at the origin.
    
  subroutine make_fcc_structure(density, unitcells_xyz)
    real(db), intent(in) :: density
    integer, dimension(ndim) :: unitcells_xyz
    
    real(db), dimension(ndim) :: edges
    real(db), dimension(ndim) :: gap, put_atom_at ! gap between atoms
    
    integer :: n = 1, nx, ny, nz, n_tot, j
    
    
    ! calculate box edges and gap. Notice to accummodate product(unitcells_xyz) 'unit cells'
    ! the box is blown up by 4^(1/3) in each direction. The resulting number of atoms
    ! is therefore 4 * product(unitcells_xyz)

    edges = unitcells_xyz / (density/4.0_db)**(1.0_db/3.0_db)
    common_config%str%box_edges = edges
    gap = edges / unitcells_xyz
    
    
    n_tot = 4*product(unitcells_xyz)
    allocate(common_config%str%atoms(n_tot))   ! allocate mass, name...
    allocate(common_config%str%r(n_tot,ndim))  ! allocate coordinates
    
    do nz = 0, unitcells_xyz(3)-1
      do ny = 0, unitcells_xyz(2)-1
        do nx = 0, unitcells_xyz(1)-1
          put_atom_at(1) = nx + 0.25_db
          put_atom_at(2) = ny + 0.25_db
          put_atom_at(3) = nz + 0.25_db
          put_atom_at = put_atom_at * gap - 0.5_db * edges
          
          do j = 0, 3
            common_config%str%r(n,:) = put_atom_at
            if (j /= 3) then
              if (j /= 0) then 
                common_config%str%r(n,1) = common_config%str%r(n,1) + 0.5_db * gap(1)
              end if
              if (j /= 1) then 
                common_config%str%r(n,2) = common_config%str%r(n,2) + 0.5_db * gap(2)
              end if
              if (j /= 2) then 
                common_config%str%r(n,3) = common_config%str%r(n,3) + 0.5_db * gap(3)
              end if                            
            end if
            n = n + 1
          end do
        end do
      end do
    end do
    
    ! print to screen
    
    write(*, *) " "
    print *, 'FCC structure created: '
    write(*,'(a,3f12.6)') "   Box size = ", edges
    print *, '  total number of atoms = ', 4*product(unitcells_xyz)  
    write(*, *) " "
    
  end subroutine make_fcc_structure  
  
  
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
    character(len=40) :: read_db
    character(len=40) :: control_object_name, read_number_atoms
    integer :: number_atoms
    
    character(len=120) :: filename

    select case(name)
      case("molecule")
        call get_value(attributes,"title", common_config%str%title, status)
        
      case("box-edges")
        call get_value(attributes,"x", read_db,status)
        common_config%str%box_edges(1) = string_to_db(read_db)

        call get_value(attributes,"y", read_db,status)
        common_config%str%box_edges(2) = string_to_db(read_db)

        call get_value(attributes,"z", read_db,status)
        common_config%str%box_edges(3) = string_to_db(read_db)
        
    
      case("atomArray")
        call get_value(attributes,"number", read_number_atoms,status)
        
        if (status /= 0) then
          write (*,*) "ERROR: number attribute missing for atomArray element"
          stop
        end if       
        
        number_atoms = string_to_int(read_number_atoms)
        allocate(common_config%str%atoms(number_atoms))  ! allocate mass, name...
        allocate(common_config%str%r(number_atoms, ndim))  ! allocate coordinates
        count_number_atoms = 1

      case("atom")
        call get_value(attributes,"x3", read_db,status)
        common_config%str%r(count_number_atoms,1) = string_to_db(read_db)

        call get_value(attributes,"y3", read_db,status)
        common_config%str%r(count_number_atoms,2) = string_to_db(read_db)

        call get_value(attributes,"z3", read_db,status)
        common_config%str%r(count_number_atoms,3) = string_to_db(read_db)

        count_number_atoms = count_number_atoms + 1        

	
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

    select case(name)
      case("atomArray")
        if (count_number_atoms /= size(common_config%str%atoms)+1) then
          write (*,*) "ERROR: number of atoms read in not equal to element attribute"
          write (*,*) count_number_atoms
          write (*,*) size(common_config%str%atoms)+1
          stop
        end if

    end select

  end subroutine end_element


  !END_DOCUMENT
  subroutine end_document()
    use flib_sax  
  end subroutine end_document

end module structure_reader_class
