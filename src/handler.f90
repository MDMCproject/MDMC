module handler_class
use flib_sax
use common_block_class
use md_control_class
use mdmc_control_class
use md_gridsearch_control_class
use mdmc_control_time_corr_class
use md_control_time_corr_class
use structure_reader_class
use control_containers_class
use converters_class

  implicit none
  
  
  private :: begin_element, end_element, pcdata_chunk
  private :: start_document, end_document
  
  public :: startup_handler
  

  logical, private  :: in_constraints = .false., in_use_near_neighbour_method = .false.
  logical, private  :: in_gpe = .false., in_fom = .false.
  logical, private  :: in_md_control = .false., in_mdmc_control = .false.
  
  ! Notice use the same setup param container for md_gridsearch_control as is used for
  ! mdmc_control - see code below
  logical, private :: in_md_gridsearch_control = .false.
  logical, private  :: in_mdmc_control_time_corr = .false. 
  logical, private  :: in_md_control_time_corr = .false.  
  
  ! nn_list_r_cut and nn_list_delta_r are the two attributes that are needed to setup
  ! a nearest neighbour list. 
  ! Notice that if nn_list_r_cut stays equal to 0.0 then nn_list%ignore_list stays 
  ! equal to .true. when attempting to make a nearest neighbour list with 
  ! make_near_neighb_list()
  ! nn_list_delta_r defines the buffer region so that the near-neighbour pairs includes
  ! all pairs with distances less than r_cut+delta_r
  real(db), private :: nn_list_r_cut = 0.0                                      
  real(db), private :: nn_list_delta_r       
                                             
  integer, dimension(ndim), private :: n_param
  
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
    
    integer :: status, n
    real(db) :: number_db
    integer :: size_of_rdf_cal_val_array  ! used in rdf-fom element
    character(len=40) :: read_db, read_int
    character(len=40) :: control_object_name, units
    character(len=120) :: filename


    select case(trim(name))
      case("constraints")
        in_constraints = .true.

      case("use-near-neighbour-method")
        in_use_near_neighbour_method = .true.
      
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
            
          case("md_gridsearch_control")
            in_md_gridsearch_control = .true.
            
          case("mdmc_control_time_corr")
            in_mdmc_control_time_corr = .true.  
            
          case("md_control_time_corr")
            in_md_control_time_corr = .true.            
        end select
	
    end select


    ! if in use-near-neighbour-method
    if (in_use_near_neighbour_method) then
      select case(name)
        case("delta-r")       
          call get_value(attributes,"val",read_db,status)
          nn_list_delta_r = string_to_db(read_db)
          
        case("r-cut")       
          call get_value(attributes,"val",read_db,status)
          nn_list_r_cut = string_to_db(read_db)      
      end select
    end if    
    
         
    ! if in md_control
    if (in_md_control) then
      select case(name)
        case("step-limit")       
          call get_value(attributes,"number",read_int,status)
          setup_md_control_params%step_limit = string_to_int(read_int)
          
        case("average-over-this-many-step")       
          call get_value(attributes,"number",read_int,status)
          setup_md_control_params%average_over_this_many_step = string_to_int(read_int)
          
        !case("perform-initial-temperature-calibration")
        !  setup_md_control_params%perform_initial_temperature_calibration = .true.               
          
        case("total-step-temp-cali")       
          call get_value(attributes,"number",read_int,status)
          setup_md_control_params%total_step_temp_cali = string_to_int(read_int)
          
        case("adjust-temp-at-interval")       
          call get_value(attributes,"number",read_int,status)
          setup_md_control_params%adjust_temp_at_interval = string_to_int(read_int)
            
        case("temperature")       
          call get_value(attributes,"val",read_db,status)
          
          ! to convert to dimensionless units
          setup_md_control_params%temperature = string_to_db(read_db) / T_unit
            
        case("time-step")       
          call get_value(attributes,"val",read_db,status)
          setup_md_control_params%time_step = string_to_db(read_db)
               
        
        case("calculate-rdf")
          setup_md_control_params%calculate_rdf = .true.   
        
        case("r-max")       
          call get_value(attributes,"val",read_db,status)
          
          ! it is assumed for now that r-max should be smaller than L/2 here
          
          number_db = string_to_db(read_db)
          if (number_db > minval(common_config%str%box_edges)/2.0) then
            number_db = minval(common_config%str%box_edges)/2.0
          end if
          
          setup_md_control_params%r_max = number_db    
          
        case("bin-length")       
          call get_value(attributes,"val",read_db,status)
          setup_md_control_params%bin_length = string_to_db(read_db)                                           
                          
        case("cal-rdf-at-interval")       
          call get_value(attributes,"number",read_int,status)
          setup_md_control_params%cal_rdf_at_interval = string_to_int(read_int)
          
        case("average-over-this-many-rdf")       
          call get_value(attributes,"number",read_int,status)
          setup_md_control_params%average_over_this_many_rdf = string_to_int(read_int)                  
          
      end select
    end if


    ! if in mdmc_control
    
    if (in_mdmc_control .or. in_md_gridsearch_control .or. in_mdmc_control_time_corr &
        .or. in_md_control_time_corr) then
      select case(name)
 
        case("total-steps-initial-equilibration")       
          call get_value(attributes,"number",read_int,status)
          setup_mdmc_control_params%total_steps_initial_equilibration = string_to_int(read_int)  
          
          call get_value(attributes,"average-over",read_int,status)
          setup_mdmc_control_params%average_over_initial_equilibration = string_to_int(read_int)                
          
        case("total-step-temp-cali")       
          call get_value(attributes,"number",read_int,status)
          setup_mdmc_control_params%total_step_temp_cali = string_to_int(read_int)
          
        case("adjust-temp-at-interval")       
          call get_value(attributes,"number",read_int,status)
          setup_mdmc_control_params%adjust_temp_at_interval = string_to_int(read_int)     
          
         case("total-step-temp-cali-repeated")       
          call get_value(attributes,"number",read_int,status)
          setup_mdmc_control_params%total_step_temp_cali_repeated = string_to_int(read_int)
          
        case("adjust-temp-at-interval-repeated")       
          call get_value(attributes,"number",read_int,status)
          setup_mdmc_control_params%adjust_temp_at_interval_repeated = string_to_int(read_int)               
                            
        case("md-steps-repeated-equilibration")       
          call get_value(attributes,"number",read_int,status)
          setup_mdmc_control_params%md_steps_repeated_equilibration = string_to_int(read_int)
          
          call get_value(attributes,"average-over",read_int,status)
          setup_mdmc_control_params%average_over_repeated_equilibration = string_to_int(read_int)                        
                          
        case("cal-rdf-at-interval")       
          call get_value(attributes,"number",read_int,status)
          setup_mdmc_control_params%cal_rdf_at_interval = string_to_int(read_int)
          
        case("average-over-this-many-rdf")       
          call get_value(attributes,"number",read_int,status)
          setup_mdmc_control_params%average_over_this_many_rdf = string_to_int(read_int)
                                  
          
        case("mc-steps")       
          call get_value(attributes,"number",read_int,status)
          setup_mdmc_control_params%mc_steps = string_to_int(read_int)            
          
            
        case("temperature")       
          call get_value(attributes,"val",read_db,status)
          
          ! to convert to dimensionless units
          setup_mdmc_control_params%temperature = string_to_db(read_db) / T_unit           
            
        case("time-step")       
          call get_value(attributes,"val",read_db,status)
          setup_mdmc_control_params%time_step = string_to_db(read_db)

        case("temperature-mc")       
          call get_value(attributes,"val",read_db,status)
          setup_mdmc_control_params%temperature_mc = string_to_db(read_db)
          
          
        ! r-max and bin-length are here only used when wanting to save g(r) to file
          
        case("r-max")       
          call get_value(attributes,"val",read_db,status)
          
          ! it is assumed for now that r-max should be smaller than L/2 here
          
          number_db = string_to_db(read_db)
          if (number_db > minval(common_config%str%box_edges)/2.0) then
            number_db = minval(common_config%str%box_edges)/2.0
          end if
          
          setup_mdmc_control_params%r_max = number_db    
          
        case("bin-length")       
          call get_value(attributes,"val",read_db,status)
          setup_mdmc_control_params%bin_length = string_to_db(read_db)        
          
          
        ! time correlation stuff  
          
        case("n-time-buffers")       
          call get_value(attributes,"number",read_int,status)
          setup_mdmc_control_params%n_time_buffers = string_to_int(read_int)    
           
        case("n-time-evals")       
          call get_value(attributes,"number",read_int,status)
          setup_mdmc_control_params%n_time_evals = string_to_int(read_int)
          
        case("n-buffer-average-over")       
          call get_value(attributes,"number",read_int,status)
          setup_mdmc_control_params%n_buffer_average_over = string_to_int(read_int)
          
        case("n-delta-t")       
          call get_value(attributes,"number",read_int,status)
          setup_mdmc_control_params%n_delta_t = string_to_int(read_int)                              
               
          
      end select
    end if
        
		
    ! if in gpe

    if (in_gpe) then
      select case(name)
        case("lj-potential")       
          call add_function(common_pe_list, target_lj_pe)
          
        case("sigma")
          call get_value(attributes,"val",read_db,status)
          call add_func_param(target_lj_pe%params, "sigma", string_to_db(read_db))
          
          call get_value(attributes,"fixed",read_db,status)
          
          if (status == 0) then 
            if (trim(read_db) == "no") then
              call set_func_param_fixed(target_lj_pe%params, "sigma", .false.)
              
              call get_value(attributes,"min",read_db,status)
              call set_func_param_min(target_lj_pe%params, "sigma", &
                                      string_to_db(read_db))
                                      
              call get_value(attributes,"max",read_db,status)
              call set_func_param_max(target_lj_pe%params, "sigma", &
                                      string_to_db(read_db)) 
                                      
              call get_value(attributes,"max-move",read_db,status)
              call set_func_param_max_move(target_lj_pe%params, "sigma", &
                                      string_to_db(read_db))   
                                      
                                                                       
            end if
          end if
          
        case("epsilon")
          call get_value(attributes,"val",read_db,status)
          call add_func_param(target_lj_pe%params, "epsilon", string_to_db(read_db))
           
          call get_value(attributes,"fixed",read_db,status)
          
          if (status == 0) then 
            if (trim(read_db) == "no") then
              call set_func_param_fixed(target_lj_pe%params, "epsilon", .false.)
              
              call get_value(attributes,"min",read_db,status)
              call set_func_param_min(target_lj_pe%params, "epsilon", &
                                      string_to_db(read_db))
                                      
              call get_value(attributes,"max",read_db,status)
              call set_func_param_max(target_lj_pe%params, "epsilon", &
                                      string_to_db(read_db)) 
                                      
              call get_value(attributes,"max-move",read_db,status)
              call set_func_param_max_move(target_lj_pe%params, "epsilon", &
                                      string_to_db(read_db))                                                            
            end if
          end if          
          
          
        case("r-cut")
          call get_value(attributes,"val",read_db,status)

          call add_func_param(target_lj_pe%params, "r-cut", string_to_db(read_db))      
          
      end select
    end if


    ! if in fom

    if (in_fom) then
      select case(name)
        case("rdf-fom")       
          call add_function(common_fom_list, target_rdf_fom)
          
        case("r-max")
          call get_value(attributes,"val",read_db,status)
          number_db = string_to_db(read_db) 
          
          ! It is assumed here that rdf-data is positioned
          ! above r-max in rdf-fom
          
          size_of_rdf_cal_val_array = floor(number_db/target_rdf_fom%rdf_data%bin_length)
          
          ! this number cannot be larger the number of data points
          
          if (size_of_rdf_cal_val_array > size(target_rdf_fom%rdf_data%val)) then
            size_of_rdf_cal_val_array = size(target_rdf_fom%rdf_data%val)
          end if
          
          target_rdf_fom%rdf_cal = make_rdf(product(common_config%str%box_edges), &
              size(common_config%str%atoms), size_of_rdf_cal_val_array, &
              target_rdf_fom%rdf_data%bin_length)  
          target_rdf_fom%hist = make_histogram(number_db, &
              target_rdf_fom%rdf_data%bin_length)          
          
        case("scale-factor")
          call get_value(attributes,"val",read_db,status)
            target_rdf_fom%scale_factor = string_to_db(read_db)
          
        case("sigma")
          call get_value(attributes,"val",read_db,status)
          number_db = string_to_db(read_db)
          target_rdf_fom%weight = 1 / (number_db*number_db)  
          
      end select
    end if
		
		
    ! if in constraints

    if (in_constraints) then
      select case(name)
        case("fnc_constraint")
          call get_value(attributes,"filename",filename,status)
          !call make_fnc_constraint(filename)
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
        if (in_md_gridsearch_control == .true.) then
          call run_md_gridsearch_control(common_config, setup_mdmc_control_params)           
        end if 
        if (in_md_control_time_corr == .true.) then
          ! check if time correlation parameters are ok
          
          if (setup_mdmc_control_params%n_time_buffers > &
              setup_mdmc_control_params%n_time_evals) then
            write(*,*) " "
            write(*,*) "ERROR in handler.f90"
            write(*,*) "n_time_buffers > n_time_evals"
            stop
          end if
          
          
          ! force n_time_evals/n_time_buffers to have no remainder since this make
          ! the buffer algorithm better work better
          
          if (mod(setup_mdmc_control_params%n_time_evals, &
                  setup_mdmc_control_params%n_time_buffers)) then
            setup_mdmc_control_params%n_time_evals = setup_mdmc_control_params%n_time_buffers &
              * (setup_mdmc_control_params%n_time_evals/setup_mdmc_control_params%n_time_buffers)
            print *, "n-time-evals adjusted reduced to ", setup_mdmc_control_params%n_time_evals
          end if
          
          
          call run_md_control_time_corr(common_config, setup_mdmc_control_params)           
        end if 
        if (in_mdmc_control_time_corr == .true.) then
          ! check if time correlation parameters are ok
          
          if (setup_mdmc_control_params%n_time_buffers > &
              setup_mdmc_control_params%n_time_evals) then
            write(*,*) " "
            write(*,*) "ERROR in handler.f90"
            write(*,*) "n_time_buffers > n_time_evals"
            stop
          end if
          
          
          ! force n_time_evals/n_time_buffers to have no remainder since this make
          ! the buffer algorithm better work better
          
          if (mod(setup_mdmc_control_params%n_time_evals, &
                  setup_mdmc_control_params%n_time_buffers)) then
            setup_mdmc_control_params%n_time_evals = setup_mdmc_control_params%n_time_buffers &
              * (setup_mdmc_control_params%n_time_evals/setup_mdmc_control_params%n_time_buffers)
            print *, "n-time-evals adjusted reduced to ", setup_mdmc_control_params%n_time_evals
          end if
          
          
          call run_mdmc_control_time_corr(common_config, setup_mdmc_control_params)           
        end if        
        
      case("fom")
        in_fom = .false.                
        
      case("use-near-neighbour-method")
      
        ! it is assumed that r-cut+delta-r should be smaller than L/2, hence
        ! check for this below and if necessary change r-cut accordingly
        
        if (nn_list_r_cut /= 0.0) then
          if (nn_list_r_cut+nn_list_delta_r > minval(common_config%str%box_edges)/2.0) then
            nn_list_r_cut = minval(common_config%str%box_edges)/2.0 - nn_list_delta_r
          end if
          
          ! Print to screen
    
          write(*, *) " "
          print *, 'Nearest neighbour method will be used with: '
          write(*,'(a,f12.6)') "   r-cut = ", nn_list_r_cut 
          write(*,'(a,f12.6)') "   delta-r = ", nn_list_delta_r
          write(*, *) " "
        end if
   
        ! make nn-list. This is regardless of whether it will be used or not because at
        ! the moment it is in make_near_neighb_list() that str%nn_list%ignore_list is 
        ! set equal to .false. if it is decided that the nearest neighbour method should 
        ! be used
        
        call make_near_neighb_list(common_config%str, nn_list_r_cut, nn_list_delta_r)
   
    end select

  end subroutine end_element


  !END_DOCUMENT
  subroutine end_document()
  end subroutine end_document

end module handler_class
