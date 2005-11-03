module handler_class
use flib_sax
use common_block_class 
use md_control_class
use structure_reader_class
use gpe_class

  implicit none
  
  public :: begin_element, end_element, pcdata_chunk
  public :: start_document, end_document
  

  logical, private  :: in_constraints = .false., in_structure = .false. 
	logical, private  :: in_gpe = .false.
  
  real(db), private :: density   
  integer, dimension(3), private :: n_atoms  ! both to be passed to 
                                             ! make_simple_cubic_structure

contains

  ! START_DOCUMENT
  subroutine start_document()
  end subroutine start_document


  ! BEGIN_ELEMENT
  subroutine begin_element(name,attributes)
    character(len=*), intent(in)   :: name
    type(dictionary_t), intent(in) :: attributes
    
    integer :: status, ndata=0
    real(db) :: number_db(1)
    integer :: number_int(1)
    character(len=40) :: read_db, read_int
    character(len=40) :: control_object_name, units
    character(len=120) :: filename

    !write(*,*) name
    select case(name)
      case("structure")
        call get_value(attributes,"filename",filename,status)
        if (status == 0) then
          call make_structure(filename)
        else
          in_structure = .true.
        end if

      case("constraints")
        in_constraints = .true.
				
		  case("gpe")
			  in_gpe = .true.

      case("control-object")
	      call get_value(attributes,"name",control_object_name,status)
	
        select case(control_object_name)
          case("moldyn_control")
            call run_md_control(common_configuration)	  

        end select
	
    end select

		
    ! if in gpe

    if (in_gpe) then
      select case(name)
        case("lj-potential")
          write(*,*) "in lj-potential"        
          call get_value(attributes,"num-variable",read_int,status)
					ndata = 0
          call build_data_array(read_int, number_int, ndata)
					common_pe_list(1)%name = "lj-potential"
					allocate(common_pe_list(1)%vars(number_int(1)))
					write (*,*) "number_int ", size(common_pe_list(1)%vars)
          
        case("sigma")
          call get_value(attributes,"val",read_db,status)
					ndata = 0
          call build_data_array(read_db, number_db, ndata)
          common_pe_list(1)%vars(1)%name = "sigma"
          common_pe_list(1)%vars(1)%val = number_db(1)
          write(*,*) common_pe_list(1)%vars(1)%name, common_pe_list(1)%vars(1)%val

        case("epsilon")
          call get_value(attributes,"val",read_db,status)
					ndata = 0
          call build_data_array(read_db, number_db, ndata)
          common_pe_list(1)%vars(2)%name = "epsilon"
          common_pe_list(1)%vars(2)%val = number_db(1)
          write(*,*) common_pe_list(1)%vars(2)%name, common_pe_list(1)%vars(2)%val
          
      end select
    end if
		
		
    ! if in constraints

    if (in_constraints) then
      select case(name)
        case("fnc_constraint")
          call get_value(attributes,"filename",filename,status)
          call make_fnc_constraint(filename)
          call add_constraint(common_configuration, my_fnc_constraint)
      end select
    end if

    
    ! if in structure
    
    if (in_structure) then
      select case(name)
        case("density")          
          ! read density value into read_dp(1)
          call get_value(attributes,"value",read_db,status)
          ndata = 0
          call build_data_array(read_db, number_db, ndata)
          
          ! read units
          call get_value(attributes,"units",units,status)
          
          if (units == "atom/AA3") then
            density = number_db(1)
          else
            ! convert other units to atom/AA3 somehow
          end if
          
        case("ar")
          ! read number of atoms into read_int(1)
          call get_value(attributes,"nx",read_int,status)
          ndata = 0
          call build_data_array(read_int, number_int, ndata)          
          n_atoms(1) = number_int(1)

          call get_value(attributes,"ny",read_int,status)
          ndata = 0
          call build_data_array(read_int, number_int, ndata)          
          n_atoms(2) = number_int(1)
          
          call get_value(attributes,"nz",read_int,status)
          ndata = 0
          call build_data_array(read_int, number_int, ndata)          
          n_atoms(3) = number_int(1)
          
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

    select case(name)
      case("constraints")
        in_constraints = .false.
        
      case("gpe")
        in_gpe = .false.        
        
      case("structure")
        call make_simple_cubic_structure(density, n_atoms)
        in_structure = .false.

    end select

  end subroutine end_element


  !END_DOCUMENT
  subroutine end_document()
  end subroutine end_document

end module handler_class
