module handler_class
use flib_sax
use common_block_class
use md_control_class
use mdmc_control_class
use structure_reader_class
use control_containers_class

  implicit none
  
  
  private :: begin_element, end_element, pcdata_chunk
  private :: start_document, end_document
  
  public :: startup_handler
  

  logical, private  :: in_constraints = .false.
	logical, private  :: in_gpe = .false., in_fom = .false.
	logical, private  :: in_md_control = .false., in_mdmc_control = .false.
  
  !real(db), private :: density   
  !integer, dimension(ndim), private :: n_atoms  ! both to be passed to 
                                             ! make_simple_cubic_structure
                                             
  integer, dimension(ndim), private :: n_param
  
  !character(len=99), private :: what_init_structure_to_build
contains
  
  subroutine startup_handler(filename)
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
  
  end subroutine startup_handler


  ! START_DOCUMENT
  subroutine start_document()
  end subroutine start_document


  ! BEGIN_ELEMENT
  subroutine begin_element(name,attributes)
    character(len=*), intent(in)   :: name
    type(dictionary_t), intent(in) :: attributes
    
    integer :: status, ndata=0, n
    real(db) :: number_db(1)
    integer :: number_int(1)
    character(len=40) :: read_db, read_int
    character(len=40) :: control_object_name, units
    character(len=120) :: filename


    select case(trim(name))
      case("constraints")
        in_constraints = .true.
				
		  case("gpe")
			  in_gpe = .true.
			  
		  case("fom")
			  in_fom = .true.			  
        
      case("control-object")
	      call get_value(attributes,"name",control_object_name,status)
	
	
        select case(control_object_name)
          case("md_control")
            in_md_control = .true.
            
          case("mdmc_control")
            in_mdmc_control = .true.
        end select
	
    end select

    
    ! if in md_control
    if (in_md_control) then
      select case(name)
      
        !case("use-near-neighbour-method")
        !  setup_md_control_params%use_near_neighbour_method = .true.
      
        case("delta-r")       
          call get_value(attributes,"val",read_db,status)
					ndata = 0
          call build_data_array(read_db, number_db, ndata)
          setup_md_control_params%delta_r = number_db(1)
          
        case("r-cut")       
          call get_value(attributes,"val",read_db,status)
					ndata = 0
          call build_data_array(read_db, number_db, ndata)
          setup_md_control_params%r_cut = number_db(1)
          
        case("step-limit")       
          call get_value(attributes,"number",read_int,status)
					ndata = 0
          call build_data_array(read_int, number_int, ndata)
          setup_md_control_params%step_limit = number_int(1)
          
        case("average-over-this-many-step")       
          call get_value(attributes,"number",read_int,status)
					ndata = 0
          call build_data_array(read_int, number_int, ndata)
          setup_md_control_params%average_over_this_many_step = number_int(1)     
          
        !case("perform-initial-temperature-calibration")
        !  setup_md_control_params%perform_initial_temperature_calibration = .true.               
          
        case("total-step-temp-cali")       
          call get_value(attributes,"number",read_int,status)
					ndata = 0
          call build_data_array(read_int, number_int, ndata)
          setup_md_control_params%total_step_temp_cali = number_int(1)       
          
        case("adjust-temp-at-interval")       
          call get_value(attributes,"number",read_int,status)
					ndata = 0
          call build_data_array(read_int, number_int, ndata)
          setup_md_control_params%adjust_temp_at_interval = number_int(1)          
            
        case("temperature")       
          call get_value(attributes,"val",read_db,status)
					ndata = 0
          call build_data_array(read_db, number_db, ndata)
          
          ! to convert to dimensionless units multiply by 8.314/1000 K^-1
          setup_md_control_params%temperature = 0.008314 * number_db(1)            
            
        case("time-step")       
          call get_value(attributes,"val",read_db,status)
					ndata = 0
          call build_data_array(read_db, number_db, ndata)
          setup_md_control_params%time_step = number_db(1) 
               
        
        case("calculate-rdf")
          setup_md_control_params%calculate_rdf = .true.   
          
        case("r-max")       
          call get_value(attributes,"val",read_db,status)
					ndata = 0
          call build_data_array(read_db, number_db, ndata)
          setup_md_control_params%r_max = number_db(1)    
          
        case("number-bins")       
          call get_value(attributes,"number",read_int,status)
					ndata = 0
          call build_data_array(read_int, number_int, ndata)
          setup_md_control_params%number_bins = number_int(1)                 
                          
        case("cal-rdf-at-interval")       
          call get_value(attributes,"number",read_int,status)
					ndata = 0
          call build_data_array(read_int, number_int, ndata)
          setup_md_control_params%cal_rdf_at_interval = number_int(1)
          
        case("average-over-this-many-rdf")       
          call get_value(attributes,"number",read_int,status)
					ndata = 0
          call build_data_array(read_int, number_int, ndata)
          setup_md_control_params%average_over_this_many_rdf = number_int(1)
                               
      end select
    end if


    ! if in mdmc_control
    
    if (in_mdmc_control) then
      select case(name)
      
        case("delta-r")       
          call get_value(attributes,"val",read_db,status)
					ndata = 0
          call build_data_array(read_db, number_db, ndata)
          setup_mdmc_control_params%delta_r = number_db(1)
          
        case("r-cut")       
          call get_value(attributes,"val",read_db,status)
					ndata = 0
          call build_data_array(read_db, number_db, ndata)
          setup_mdmc_control_params%r_cut = number_db(1)
          
        case("md-steps-per-trajectory")       
          call get_value(attributes,"number",read_int,status)
					ndata = 0
          call build_data_array(read_int, number_int, ndata)
          setup_mdmc_control_params%md_steps_per_trajectory = number_int(1)
          
        case("mc-steps")       
          call get_value(attributes,"number",read_int,status)
					ndata = 0
          call build_data_array(read_int, number_int, ndata)
          setup_mdmc_control_params%mc_steps = number_int(1)           
            
        case("temperature")       
          call get_value(attributes,"val",read_db,status)
					ndata = 0
          call build_data_array(read_db, number_db, ndata)
          
          ! to convert to dimensionless units multiply by 8.314/1000 K^-1
          setup_mdmc_control_params%temperature = 0.008314 * number_db(1)            
            
        case("time-step")       
          call get_value(attributes,"val",read_db,status)
					ndata = 0
          call build_data_array(read_db, number_db, ndata)
          setup_mdmc_control_params%time_step = number_db(1) 
                               
      end select
    end if
        
		
    ! if in gpe

    if (in_gpe) then
      select case(name)
        case("lj-potential")       
          call add_potential(common_pe_list, target_lj_pe)
          
        case("sigma")
          call get_value(attributes,"val",read_db,status)
					ndata = 0
          call build_data_array(read_db, number_db, ndata)
          
          call add_func_param(target_lj_pe%params, "sigma", number_db(1))
          
        case("epsilon")
          call get_value(attributes,"val",read_db,status)
					ndata = 0
          call build_data_array(read_db, number_db, ndata)
          
          call add_func_param(target_lj_pe%params, "epsilon", number_db(1))
          
        case("r-cut")
          call get_value(attributes,"val",read_db,status)
					ndata = 0
          call build_data_array(read_db, number_db, ndata)
          
          call add_func_param(target_lj_pe%params, "r-cut", number_db(1))      
          
      end select
    end if


    ! if in fom

    if (in_fom) then
      select case(name)
        case("rdf-fom")       
          call add_potential(common_fom_list, target_rdf_fom)
      
          
        case("scale-factor")
          call get_value(attributes,"val",read_db,status)
					ndata = 0
          call build_data_array(read_db, number_db, ndata)
          
          target_rdf_fom%scale_factor = number_db(1)
          
        case("sigma")
          call get_value(attributes,"val",read_db,status)
					ndata = 0
          call build_data_array(read_db, number_db, ndata)
          
          target_rdf_fom%weight = 1 / (number_db(1)*number_db(1))  
          
      end select
    end if
		
		
    ! if in constraints

    if (in_constraints) then
      select case(name)
        case("fnc_constraint")
          call get_value(attributes,"filename",filename,status)
          call make_fnc_constraint(filename)
          call add_constraint(common_config, target_fnc_constr)
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
        
      case("control-object")
        if (in_md_control == .true.) then
          call run_md_control(common_config, setup_md_control_params)
        end if
        if (in_mdmc_control == .true.) then
          call run_mdmc_control(common_config, setup_mdmc_control_params)
        end if
        
    end select

  end subroutine end_element


  !END_DOCUMENT
  subroutine end_document()
  end subroutine end_document

end module handler_class
