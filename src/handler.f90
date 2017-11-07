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
use time_corr_algorithm_class

  implicit none
  
  
  private :: begin_element, end_element, pcdata_chunk
  private :: start_document, end_document
  
  private :: check_when_temperature_cali
  
  public :: startup_handler
  
  ! The logical variables named 'in_ + something' are used to keep track of when the 
  ! XML reader is inside or not inside an XML element
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
  
  integer, private :: n_md_step_between_buffers  ! number of md steps between buffers
  integer, private :: n_buffers  ! number of buffers
  
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
    
    integer :: status, n, index
    real(db) :: number_db
    real(db) :: l_start, l_end, l_step  ! used for reading in q and omega values
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
          ! check unit
          call get_value(attributes,"units",units,status)
          if (units /= "AA") then
            write(*,*) "ERROR, unit of delta-r must be in AA"
            stop
          end if 

          call get_value(attributes,"val",read_db,status)
          nn_list_delta_r = string_to_db(read_db)
          
        case("r-cut")
          ! check unit
          call get_value(attributes,"units",units,status)
          if (units /= "AA") then
            write(*,*) "ERROR, unit of r-cut must be in AA"
            stop            
          end if
          
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
          
        case("total-step-temp-cali")       
          call get_value(attributes,"number",read_int,status)
          setup_md_control_params%total_step_temp_cali = string_to_int(read_int)
          
        case("adjust-temp-at-interval")       
          call get_value(attributes,"number",read_int,status)
          setup_md_control_params%adjust_temp_at_interval = string_to_int(read_int)
            
        case("temperature")       
          ! check unit
          call get_value(attributes,"units",units,status)
          if (units /= "Kelvin") then
            write(*,*) "ERROR, unit of temperature must be in Kelvin"
            stop            
          end if            
            
          ! read temperature and convert it to internal temperature unit
          call get_value(attributes,"val",read_db,status)
          setup_md_control_params%temperature = string_to_db(read_db) / T_unit
            
        case("time-step")       
          ! check unit
          call get_value(attributes,"units",units,status)
          if (units /= "10e-13s") then
            write(*,*) "ERROR, unit of time-step must be in 10e-13s"
            stop            
          end if 
          
          call get_value(attributes,"val",read_db,status)
          setup_md_control_params%md_delta_t = string_to_db(read_db)
        
        case("calculate-rdf")
          setup_md_control_params%calculate_rdf = .true.   
        
        case("r-max")
          ! check unit
          call get_value(attributes,"units",units,status)
          if (units /= "AA") then
            write(*,*) "ERROR, unit of r-max must be in AA"
            stop            
          end if
          
          call get_value(attributes,"val",read_db,status)
          
          ! it is assumed for now that r-max should be smaller than half
          ! the shortest box edge
          
          number_db = string_to_db(read_db)
          if (number_db > minval(common_config%str%box_edges)/2.05) then
            number_db = minval(common_config%str%box_edges)/2.05
            print *, "Job file input r-max adjusted to be less than half of the shortest box edge."            
            write(*,'(a,f10.4)') " r-max adjusted to: ", number_db
            write(*, *) " "
          end if
          
          setup_md_control_params%r_max = number_db    
          
        case("bin-length")
          ! check unit
          call get_value(attributes,"units",units,status)
          if (units /= "AA") then
            write(*,*) "ERROR, unit of bin-length must be in AA"
            stop            
          end if
          
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


    ! If in mdmc_control, md_gridsearch_control, mdmc_control_time_corr or md_control_time_corr.
    ! All these algorithms store user input in the same container
    
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
          
          ! read and to convert to internal temperature units
          setup_mdmc_control_params%temperature = string_to_db(read_db) / T_unit           
            
        case("time-step")       
          call get_value(attributes,"val",read_db,status)
          setup_mdmc_control_params%md_delta_t = string_to_db(read_db)          
          
        case("temperature-mc")       
          call get_value(attributes,"val",read_db,status)
          setup_mdmc_control_params%temperature_mc = string_to_db(read_db)
          
          
        ! r-max and bin-length are here only used when wanting to save g(r) to file
          
        case("r-max")
          call get_value(attributes,"val",read_db,status)
          
          ! it is assumed for now that r-max should be smaller than half
          ! the shortest box edge
          
          number_db = string_to_db(read_db)
          if (number_db > minval(common_config%str%box_edges)/2.05) then
            number_db = minval(common_config%str%box_edges)/2.05
            print *, "Job file input r-max adjusted to be less than half of the shortest box edge."            
            write(*,'(a,f10.4)') " r-max adjusted to: ", number_db
            write(*, *) " "
          end if
          
          setup_mdmc_control_params%r_max = number_db    
          
        case("bin-length")       
          call get_value(attributes,"val",read_db,status)
          setup_mdmc_control_params%bin_length = string_to_db(read_db)        
          
          
        ! time correlation stuff, i.e. that contained within the <time-correlation> element
          
        case("n-md-step-between-buffers")       
          ! This parameter is used to add a time delay between
          ! buffers, where one buffer calculates one g(r,t). See also the
          ! comments elsewhere in this file where n_md_step_between_buffers is used
          call get_value(attributes,"number",read_int,status)
          n_md_step_between_buffers = string_to_int(read_int);                   

        case("n-time-bin")       
          call get_value(attributes,"number",read_int,status)
          setup_mdmc_control_params%n_time_bin = string_to_int(read_int)

        case("n-g-r-t-to-average-over")       
          call get_value(attributes,"number",read_int,status)
          
          ! set this time_corr_container private attribute directly
          call set_n_g_r_t_to_average_over(string_to_int(read_int))
          
        case("md-per-time-bin")       
          call get_value(attributes,"number",read_int,status) 
          setup_mdmc_control_params%md_per_time_bin = string_to_int(read_int)                              
               
          
      end select
    end if
        
		
    ! if in gpe

    if (in_gpe) then
      select case(name)
        case("lj-potential")       
          call add_function(common_pe_list, target_lj_pe)
          
        case("sigma")
          ! check unit
          call get_value(attributes,"units",units,status)
          if (units /= "AA") then
            write(*,*) "ERROR, unit of sigma must be in AA"
            stop            
          end if
          
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
          ! check unit
          call get_value(attributes,"units",units,status)
          if (units /= "KJ/mole") then
            write(*,*) "ERROR, unit of epsilon must be in KJ/mole"
            stop            
          end if
          
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
          ! check unit
          call get_value(attributes,"units",units,status)
          if (units /= "AA") then
            write(*,*) "ERROR, unit of r-cut must be in AA"
            stop            
          end if
          
          call get_value(attributes,"val",read_db,status)
          call add_func_param(target_lj_pe%params, "r-cut", string_to_db(read_db))      
          
      end select
    end if


    ! if in fom

    if (in_fom) then
      select case(name)
        case("rdf-fom")       
          call add_function(common_fom_list, target_rdf_fom)
          
        case("g-d-rt-fom")       
          call add_function(common_fom_list, target_g_d_rt_fom)          

        case("s-qt-fom")       
          call add_function(common_fom_list, target_s_qt_fom)   

        case("s-qo-fom")       
          call add_function(common_fom_list, target_s_qo_fom)   
                              
        case("r-max")
          ! Within the job file and within the <rdf-fom> element 
          ! it is required that the <rdf-data> element is positioned above 
          ! the <r-max> element. That is that the rdf-data has been read
          ! and processed at this point.
          ! TODO: fix this so this is no longer required (see related xmlf90 comments in mdmc.f90)
          
          if (target_rdf_fom%rdf_data%bin_length == 0.0) then
            print *, " "
            print *, "Within the job file and within the <rdf-fom> element"
            print *, "it is required that the <rdf-data> element is positioned"
            print *, "above the <r-max> element"
            stop
          end if
            
          call get_value(attributes,"val",read_db,status)
          number_db = string_to_db(read_db)           

          ! Set the size of the calculated rdf to be <r-max> divided the data bin lenght
          
          size_of_rdf_cal_val_array = floor(number_db/target_rdf_fom%rdf_data%bin_length) 
          
          ! this number is here restricted to be smaller than or equal to the number of 
          ! data points. This is done to avoid wasting time calculating g(r) for points
          ! for which no data are available
          
          if (size_of_rdf_cal_val_array > size(target_rdf_fom%rdf_data%val)) then
            size_of_rdf_cal_val_array = size(target_rdf_fom%rdf_data%val)
          end if
          
          ! Allocate space for rdf_cal
          
          target_rdf_fom%rdf_cal = make_rdf(product(common_config%str%box_edges), &
              size(common_config%str%atoms), size_of_rdf_cal_val_array, &
              target_rdf_fom%rdf_data%bin_length)       
 
          ! set stuff which enable user to create histogram suitable for calculating rdf
          
          setup_mdmc_control_params%n_r_bin_cal = size_of_rdf_cal_val_array
          setup_mdmc_control_params%bin_length_cal = target_rdf_fom%rdf_data%bin_length      
          
        case("scale-factor")
          call get_value(attributes,"val",read_db,status)
            target_rdf_fom%scale_factor = string_to_db(read_db)
            target_s_qo_fom%scale_factor = string_to_db(read_db)
            target_s_qt_fom%scale_factor = string_to_db(read_db)
          
        case("sigma")
          call get_value(attributes,"val",read_db,status)
          number_db = string_to_db(read_db)
          target_rdf_fom%weight = 1 / (number_db*number_db)
          
        case("ignore-errors")
          target_s_qo_fom%ignore_errors = .true.         
          
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

    integer :: n_time_bin_between_buffers

    select case(name)
      case("constraints")
        in_constraints = .false.
        
      case("gpe")
        in_gpe = .false.        
        
      case("control-object")
        
        ! Do checks of MD steps values and temperature calibration values
        
        if ( setup_md_control_params%step_limit > 0 ) then
          call check_when_temperature_cali(setup_md_control_params%step_limit, & 
                                           setup_md_control_params%total_step_temp_cali, &
                                           setup_md_control_params%average_over_this_many_step, &
                                           setup_md_control_params%adjust_temp_at_interval)
        end if
        
        if ( setup_mdmc_control_params%total_steps_initial_equilibration > 0 ) then
          call check_when_temperature_cali(setup_mdmc_control_params%total_steps_initial_equilibration, & 
                                           setup_mdmc_control_params%total_step_temp_cali, &
                                           setup_mdmc_control_params%average_over_initial_equilibration, &
                                           setup_mdmc_control_params%adjust_temp_at_interval)            
        end if
        
        if ( setup_mdmc_control_params%md_steps_repeated_equilibration > 0 ) then
          call check_when_temperature_cali(setup_mdmc_control_params%md_steps_repeated_equilibration, & 
                                           setup_mdmc_control_params%total_step_temp_cali_repeated, &
                                           setup_mdmc_control_params%average_over_repeated_equilibration, &
                                           setup_mdmc_control_params%adjust_temp_at_interval_repeated)              
        end if
        
        
        ! Start running an algorithm 
          
        if (in_md_control == .true.) then
          call run_md_control(common_config, setup_md_control_params)
        end if
        
        if (in_mdmc_control == .true.) then
          call run_mdmc_control(common_config, setup_mdmc_control_params)
        end if
        
        if (in_md_gridsearch_control == .true.) then
          call run_md_gridsearch_control(common_config, setup_mdmc_control_params)           
        end if
        
        
        ! Before calling algorithms that compare against dynamical structure factor information
        ! do the following addition checks and adjustments
        
        if (in_md_control_time_corr == .true. .or. in_mdmc_control_time_corr == .true.) then
          
          ! Within a 'buffer' a full g(r,t) is calculated.
          ! A time bin for g(r,t) has the length: MD time step * md_per_time_bin.
          ! Determine the number of time bins to use between when buffers individual 
          ! buffers.
          
          if (n_md_step_between_buffers > setup_mdmc_control_params%md_per_time_bin * &
                                          setup_mdmc_control_params%n_time_bin) then
            n_md_step_between_buffers = setup_mdmc_control_params%md_per_time_bin * setup_mdmc_control_params%n_time_bin
            print *, "Job file input n-md-step-between-buffers too large."
            print *, "Must to smaller than or equal to md-per-time-bin * n-time-bin."
            print *, "n-md-step-between-buffers adjusted to: ", n_md_step_between_buffers
            print *, " "
          end if

          n_time_bin_between_buffers = nint( dble(n_md_step_between_buffers) &
            / dble(setup_mdmc_control_params%md_per_time_bin) )
                    
          print *, "Number of time bins between buffers for calculating g(r,t) = ", n_time_bin_between_buffers
          print *, " "

          ! For computing g(r,t), one per buffer, it is computationally most efficient if the number
          ! of time bins used for calculating a g(r,t) is a multiple of the time bins between when 
          ! individual buffers starts 
          ! Therefore, adjust n_time_bin:
          
          if ( mod(setup_mdmc_control_params%n_time_bin, n_time_bin_between_buffers) /= 0) then
            setup_mdmc_control_params%n_time_bin = n_time_bin_between_buffers &
              * ceiling(dble(setup_mdmc_control_params%n_time_bin)/dble(n_time_bin_between_buffers))           
            print *, "Job file n-time-bin adjusted to be a multiple of"
            print *, "md-per-time-bin * n-md-step-between-buffers."
            print *, "Job file n-time-bin adjusted to: ", setup_mdmc_control_params%n_time_bin
            print *, " "
          end if      
          
          ! A sanity check of the above code
          
          if (n_time_bin_between_buffers > &
              setup_mdmc_control_params%n_time_bin) then
            write(*,*) " "
            write(*,*) "ERROR in handler.f90"
            write(*,*) "n_time_bin_between_buffers > setup_mdmc_control_params%n_time_bin"
            stop
          end if
    
          ! Further it is a waste of memory (and computation) if the number of buffers
          ! larger than the number of g(r,t) that needs calculating, which is n-g-r-t-to-average-over
          ! TODO: adjust this automatically for the user instead
          
          n_buffers = setup_mdmc_control_params%n_time_bin/n_time_bin_between_buffers
          if (n_buffers > get_n_g_r_t_to_average_over()) then
            write(*,*) " "
            write(*,*) "In Job file please increase n-md-step-between-buffers or n-g-r-t-to-average-over"
            write(*,*) "otherwise when g(r,t) more buffers are created than needed."
            write(*,*) "That is more buffers than n-g-r-t-to-average-over which is inefficient."
            stop            
          end if
              
          ! Set the private attribute n_buffer in time_corr_container. n_buffer is
          ! the number of time bins between when g(r,t)s are calculated
          
          call set_n_buffers( n_buffers )
          
          
          if (in_mdmc_control_time_corr == .true.) then
            call run_mdmc_control_time_corr(common_config, setup_mdmc_control_params)
          else
            call run_md_control_time_corr(common_config, setup_mdmc_control_params)    
          end if             
        end if 
        
      case("fom")
        in_fom = .false.                
        
      case("use-near-neighbour-method")
      
        ! it is assumed that r-cut+delta-r should be smaller than L/2, hence
        ! check for this below and if necessary change r-cut accordingly
        
        if (nn_list_r_cut /= 0.0) then
          if (nn_list_r_cut+nn_list_delta_r > minval(common_config%str%box_edges)/2.05) then
            nn_list_r_cut = minval(common_config%str%box_edges)/2.05 - nn_list_delta_r
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
  
  ! Additional user input check when temperature calibration is first part of 
  ! a MD simulation 
  
  subroutine check_when_temperature_cali(total_steps, temperature_steps, average_over, cali_when)
    integer, intent(in) :: total_steps
    integer, intent(in) :: temperature_steps
    integer, intent(in) :: average_over
    integer, intent(in) :: cali_when    
    
    ! To ensure that no MD steps are vasted then insist
    
    if ( mod(total_steps, average_over) /= 0) then        
      write(*,*) "ERROR in input in job file"
      write(*,*) "To ensure that no MD steps are vasted then: "      
      write(*,*) "The total number of steps in MD simulation must be a multiple of "
      write(*,*) "the average-over XML attribute."
      stop
    end if      
    
    
    if ( temperature_steps > 0 ) then
      return    
    end if
    
    ! To ensure that no MD steps are vasted then insist
    
    if ( mod(temperature_steps, cali_when) /= 0) then        
      write(*,*) "ERROR in input in job file"
      write(*,*) "To ensure that no MD steps are vasted then: "      
      write(*,*) "The total number of MD steps temperature calibration must be a multiple of "
      write(*,*) "the number attribute of the adjust-temp-at-interval-repeated XML element."
      stop
    end if            
    
    ! Ensure that temperature calibration can fit simulation
    
    if ( temperature_steps < total_steps ) then        
      write(*,*) "ERROR in input in job file"
      write(*,*) "Number of steps used for temperature calibration must be smaller than "      
      write(*,*) "the number of steps assigned for the MD simulation."
      stop
    end if      
    
    ! Required to ensure enough MD steps to check if total energy has drifted since last 
    ! temperature calibration of velocities
    
    if ( mod(total_steps-temperature_steps, average_over) /= 0) then        
      write(*,*) "ERROR in input in job file"
      write(*,*) "There must be a gap between end of temperature calibration and end of "      
      write(*,*) "simulation is a multiple of the average-over XML attribute."
      stop
    end if    
  

  end subroutine check_when_temperature_cali 

end module handler_class
