module mdmc_control_time_corr_class
use configuration_class
use function_class
use common_block_class, only : common_config, common_pe_list, common_fom_list, target_rdf_fom
use phasespace_class
use md_properties_class
use control_containers_class
use tic_toc_class
use rdf_class
use structure_reader_class
use func_params_wrapper_class
use time_corr_algorithm_class
use time_corr_hist_container_class
use s_q_time_class
use s_q_omega_class

  implicit none

  public :: run_mdmc_control_time_corr
  
  private :: acceptable_temperature
  private :: acceptable_energy

contains

  subroutine run_mdmc_control_time_corr(a_config, c)
    use flib_wxml  
    type (configuration), intent(inout) :: a_config
    type (mdmc_control_container) :: c
    
    real (db) :: time_now = 0.0    
    integer :: i, i_md, j
    real(db) :: density  ! used for cal input argument to save_rdf
    
    real(db) :: sum_kin_energy
    real(db) :: temp_adjust_factor
    
    real(db) :: average_energy_end_of_temp_calibration
  
    real(db) :: fom_val, fom_old

    real(db) :: fom_best_including_rejected  ! used for keeping track of best FOM including rejected
    real(db) :: fom_best_accepted  ! used for keeping track of best accepted FOM (as a sanity check)
    type (phasespace) :: my_ps_best ! to store status of best phasespace

    type (phasespace) :: my_ps, my_ps_old
    type (md_properties) :: my_props
    real(db) :: pressure_comp = 0.0, pot_energy = 0.0

    type(time_corr_hist_container) :: my_time_corr_container
    
    type (s_q_time) :: my_s_q_time
    type (s_q_omega) :: my_s_q_omega    
	
    integer :: print_to_file = 555
    integer :: print_to_screen = 0
    
    real (db) :: ran_num  ! for random generator
    real (db) :: delta_fom ! for metropolis
    logical :: accept_parameters ! for metropolis
    
    type (xmlf_t) :: xf
    call xml_OpenFile("output/mdmc_results.xml", xf, indent=.true.)
    call xml_AddXMLDeclaration(xf, "UTF-8")
    call xml_NewElement(xf, "mdmc-control-results")
    call xml_AddAttribute(xf, "title", "rho = " // str(density, format="(f10.5)") &
                                         // "atoms/AA-3")
    call xml_NewElement(xf, "this-file-was-created")
    call xml_AddAttribute(xf, "when", get_current_date_and_time())
    call xml_EndElement(xf, "this-file-was-created")
	  
    if (print_to_file /= 0) then
      open(print_to_file, file="output/job_log.txt")
    end if
		
    write(print_to_file,*) "In run_mdmc_control_time_corr"

    call tic

    ! prepare phasespaces
    
    my_ps = make_phasespace(a_config%str, c%temperature)
    my_ps_old = copy_phasespace(my_ps)
    my_ps_best = copy_phasespace(my_ps)

    ! prepare various other containers

    my_time_corr_container = make_time_corr_hist_container(c%r_max, c%bin_length, c%n_time_bin, &
                             c%md_per_time_bin * c%md_delta_t)                       
    my_s_q_time = make_s_q_time(c%q_values, c%md_per_time_bin * c%md_delta_t, c%n_time_bin)
    my_s_q_omega = make_s_q_omega(c%q_values, c%omega_values)
    
    ! Setting various variables that will be accumulated to zero
    
    call clear_time_corr_hist_container(my_time_corr_container)
                                         
! -------------- initial equilibration ---------------- !

    sum_kin_energy = 0.0
    
    do i = 1, c%total_steps_initial_equilibration
      time_now = c%md_delta_t * i   ! perhaps print this one out 
      
      ! Do one trajectory of length = 1 where pressure_comp and pot_energy is also
      ! calculated
      
      call trajectory_in_phasespace(my_ps, common_pe_list, 1, c%md_delta_t, & 
                                    pressure_comp, pot_energy)
      
      ! Calculate proproperties from the phase-space configuration like the kinetic energy
      ! including moving averages of these
      
      call md_cal_properties(my_ps, my_props, common_pe_list, pressure_comp, pot_energy)
      
      ! Optionally adjust the temperature in the initial stage of this first simulation

      if (i < c%total_step_temp_cali) then
        ! Needs the average kin energy for adjusting the temperature.
        ! This could done more clean (code) by introducing a second instance of my_props
        ! i.e. not keep track of the 
          
        sum_kin_energy = sum_kin_energy + my_props%kin_energy%val
        
        if (mod(i,c%adjust_temp_at_interval) == 0) then
          ! see src/algorithms/md_gridsearch_control.doc for explanation of the expression below

          temp_adjust_factor = sqrt(c%adjust_temp_at_interval * 1.5 * c%temperature / &
              sum_kin_energy * (size(my_ps%str%atoms)-1.0) / size(my_ps%str%atoms))
          my_ps%p = my_ps%p * temp_adjust_factor
          sum_kin_energy = 0.0
        end if
      end if
      
      ! The idea is to keep a record of what the average total energy is after the
      ! temperature calibration and then compare this value at end of the entire
      ! initial calibration to check if the total energy has not drifted due to 
      ! numerical errors. Not sure if this is the best point to record this total average
      
      if (i == c%total_step_temp_cali) then
        average_energy_end_of_temp_calibration = my_props.tot_energy.ave 
      end if      

  
      ! print out stuff at interval = average_over_initial_equilibration
        
      if (mod(i,c%average_over_initial_equilibration) == 0) then 
        call md_print_properties(print_to_file, my_props)
        
        if (sum(sum(my_ps%p,1)) > 0.0001) then
          write(*,*) "ERROR:"
          write(*,'(a,3f12.6)') "total momentum ", sum(my_ps%p,1)
          write(*,*) "Serious problem - total momentum different from zero"
          stop
        end if
        
        ! Reset calculation of rolling averages
        
        call md_reset_properties(my_props)
        
        write(*, '(a,i8,a,f12.4,a)') "MD steps = ", i, " MD run-time = ", time_now, "*10e-13"
      end if
      
    end do
    
    write(print_to_file, *) " "
    write(print_to_file, '(a,f12.4,a)') "Taken ", toc(), " seconds to execute initial equilibration."
    write(print_to_file, *) " "


    ! Determine if the MD simulation reached an acceptable equilibrium
    !
    ! Note the (2.0/ndim) factor is to convert from dimensionless kin_energy per atom to 
    ! dimensionless temperature
    !
    ! First up check the temperature and the end of the simulation is OK
    
    if ( acceptable_temperature((2.0/ndim)*my_props.kin_energy.ave, &
         c%temperature, 0.2d+0) == .false.) then
         write(print_to_screen, *) "Initial equilibration did not reach equilibrium"
         write(print_to_screen, *) "Temperature outside acceptable value"
         stop
    end if   


    ! Secondly check that the total energy has not drifted since the temperature
    ! calibration

    if ( acceptable_energy(average_energy_end_of_temp_calibration, &
         my_props.tot_energy.ave, 0.1d+0) == .false.) then
         write(print_to_screen, *) "Initial equilibration did not reach equilibrium - STOP"
         write(print_to_screen, *) "Total energy drift is too high since temperature calibration ended"
         stop
    end if     
    
 ! ---------------  finished initial equilibration -------------- !      

 
 ! -------- calculate first FOM ------------------------- !

  density = size(my_ps%str%atoms) / product(my_ps%str%box_edges) ! for printing 

  ! cal 1st time-dependent correlation container - i.e. move phase-space forward
  ! as required for this and average over 'x' number of these containers as specified in the MDMC job file
  
  call cal_time_corr_container(my_time_corr_container, my_ps, common_pe_list, c%md_per_time_bin, c%md_delta_t)
  
  ! from container calculate FOM etc  
  
  call print_g_d(my_time_corr_container, product(my_ps%str%box_edges), size(my_ps%str%atoms), &
                 c%temperature, "output/first_g_d.xml") 
  call cal_s_q_time(my_time_corr_container, my_ps%str, my_s_q_time)
  call print_s_q_time(my_s_q_time, density, c%temperature, "output/first_s_q_time.xml")
  call cal_s_q_omega(my_s_q_time, my_ps%str, my_s_q_omega)
  call print_s_q_omega(my_s_q_omega, density, c%temperature, "output/first_s_q_omega.xml")
  fom_val = func_val(my_s_q_omega, common_fom_list)
  call clear_time_corr_hist_container(my_time_corr_container)
  
  write(print_to_file,'(a,f14.5)') "1st FOM = ", fom_val
  write(print_to_file, '(a,f14.5)') "Finished cal 1st FOM. Time: ", toc()
  write(print_to_file, *) " "
  write(print_to_screen,'(a,f14.5)') "1st FOM = ", fom_val
  write(print_to_screen, '(a,f14.5)') "Finished cal 1st FOM. Time: ", toc()

  ! print FOM to xml file
  
  call xml_NewElement(xf, "accept")
  call add_xml_attribute_func_params(xf, common_pe_list)
  call xml_AddAttribute(xf, "val", str(fom_val, format="(f14.5)"))
  call xml_EndElement(xf, "accept")    

  ! store best solution so far
  
  fom_best_accepted = fom_val
  fom_best_including_rejected = fom_val
  call backup_best_func_params(common_pe_list)
  call shallow_copy_phasespace(my_ps, my_ps_best)      
 
 ! -------- finished calculating first FOM -------------------------- !



 ! ---------------------- core mdmc part ------------------------ !              
               
  
    ! save state
    
    fom_old = fom_val
    call backup_func_params(common_pe_list)
    call shallow_copy_phasespace(my_ps, my_ps_old)       
                    
    do i = 1, c%mc_steps

      write (print_to_screen, *) "Begin MC step number ", i
      write (print_to_file, *) "Begin MC step number ", i
      write (print_to_file, *) " "

      call move_random_func_params(common_pe_list)
      
      ! do MD simulation
      
      sum_kin_energy = 0.0
      
      do i_md = 1, c%md_steps_repeated_equilibration

        call trajectory_in_phasespace(my_ps, common_pe_list, 1, c%md_delta_t, & 
                                  pressure_comp, pot_energy)

        call md_cal_properties(my_ps, my_props, common_pe_list, pressure_comp, pot_energy)       
      
        ! As for initial MD simulation here also optionally adjust the temperature
        
        if (i_md < c%total_step_temp_cali_repeated) then
          sum_kin_energy = sum_kin_energy + my_props%kin_energy%val
          if (mod(i_md,c%adjust_temp_at_interval_repeated) == 0) then
            temp_adjust_factor = sqrt(c%adjust_temp_at_interval_repeated * 1.5 * c%temperature / &
                sum_kin_energy * (size(my_ps%str%atoms)-1.0) / size(my_ps%str%atoms))
            my_ps%p = my_ps%p * temp_adjust_factor
            sum_kin_energy = 0.0
          end if
        end if    
        
        ! store total energy after temperature calibration to later check if
        ! this energy has drifted too much due to numerical errors
        
        if (i_md == c%total_step_temp_cali_repeated) then
          average_energy_end_of_temp_calibration = my_props.tot_energy.ave 
        end if      
        
        if (mod(i_md,c%average_over_repeated_equilibration) == 0) then 
          call md_print_properties(print_to_file, my_props)

          call md_reset_properties(my_props)
        end if
              
      end do
         
      write(print_to_file, '(a,f14.5)') "Finished MD simulation within MC loop. Time: ", toc() 
      write(print_to_file, *) " "
      write(print_to_screen, '(a,f14.5)') "Finished MD simulation within MC loop. Time: ", toc()      
      

      ! Determine if the MD simulation reached an acceptable equilibrium.
      ! If not then reset phase-space and start a new MC cycle
      !
      ! Note the (2.0/ndim) factor is to convert from dimensionless kin_energy per atom to 
      ! dimensionless temperature
      !
      ! Check temperature OK and if energy has been drifting too much
      
      if ( acceptable_temperature((2.0/ndim)*my_props.kin_energy.ave, &
          c%temperature, 0.2d+0) == .false. .or. acceptable_energy(average_energy_end_of_temp_calibration, &
          my_props.tot_energy.ave, 0.1d+0) == .false.) then
        write(print_to_screen, *) "Did not reach equilibrium in MD run within MC loop"
        write(print_to_file, *) "Did not reach equilibrium in MD run within MC loop"

        ! print out to mdmc_results.xml that that MD simulation failed to reach
        ! equilibrium for this PE parameter values choice
        
        call xml_NewElement(xf, "failed-md-run")
        call add_xml_attribute_func_params(xf, common_pe_list)
        call xml_EndElement(xf, "failed-md-run")   

        call restore_func_params(common_pe_list)
        call shallow_copy_phasespace(my_ps_old, my_ps)         

        cycle
      end if   

      ! cal time-dependent container, FOM etc
    
      call cal_time_corr_container(my_time_corr_container, my_ps, common_pe_list, c%md_per_time_bin, c%md_delta_t)   
      call cal_s_q_time(my_time_corr_container, my_ps%str, my_s_q_time)
      call cal_s_q_omega(my_s_q_time, my_ps%str, my_s_q_omega)
      fom_val = func_val(my_s_q_omega, common_fom_list)
      call clear_time_corr_hist_container(my_time_corr_container)

      write(print_to_file,'(a,f14.5)') "FOM = ", fom_val
      write(print_to_file, '(a,f14.5)') "Finished cal FOM. Time: ", toc()
      write(print_to_file, *) " "
      write(print_to_screen,'(a,f14.5)') "FOM = ", fom_val
      write(print_to_screen, '(a,f14.5)') "Finished cal FOM. Time: ", toc()
      
      
      ! Check if new best FOM
      ! Here whether or not this parameter move will be rejected or accepted 
      ! If the new FOM value is lower it should always be accepted, but as a
      ! sanity check if this is indeed the case by tracking two fom_best...
      
      if (fom_val < fom_best_including_rejected) then
        fom_best_including_rejected = fom_val
        call backup_best_func_params(common_pe_list)
        call shallow_copy_phasespace(my_ps, my_ps_best)      
      end if
     
      
      ! Metropolis check

      delta_fom = fom_val - fom_old

      accept_parameters = .false.

      if(delta_fom <= 0.0) then
          ! new parameters accepted
        accept_parameters = .true.    
      else
        ! get random number between 0 and 1
        call random_number(ran_num)
          
        if(exp(- delta_fom / c%temperature_mc) > ran_num) then
          accept_parameters = .true.
        end if
      end if
 
        
      if (accept_parameters) then
      
        call xml_NewElement(xf, "accept")
        call xml_AddAttribute(xf, "N", str(i, format="(i)"))
        call add_xml_attribute_func_params(xf, common_pe_list)
        call xml_AddAttribute(xf, "val", str(fom_val, format="(f14.5)"))
        call xml_EndElement(xf, "accept")       
     
      
        ! save state
    
        fom_old = fom_val
        call backup_func_params(common_pe_list)
        call shallow_copy_phasespace(my_ps, my_ps_old)
        
        ! For a sanity check later, here keep record of the best accepted FOM found 
      
        if (fom_val < fom_best_accepted) then
          fom_best_accepted = fom_val     
        end if        
        
      else
      
        call xml_NewElement(xf, "rejected")
        call xml_AddAttribute(xf, "N", str(i, format="(i)"))
        call add_xml_attribute_func_params(xf, common_pe_list)
        call xml_AddAttribute(xf, "val", str(fom_val, format="(f14.5)"))
        call xml_EndElement(xf, "rejected")
        
        
        ! restore state
    
        call restore_func_params(common_pe_list)
        call shallow_copy_phasespace(my_ps_old, my_ps)          
              
      end if  
      
    end do
    
    ! Another sanity check
    
    if ( fom_best_including_rejected /= fom_best_accepted ) then
         write(print_to_screen, *) "Fundamental error - STOP"
         write(print_to_screen, *) "The best FOM found among accepted and rejected parameter moves"
         write(print_to_screen, *) "different from the best FOM found among accepted moves only."
         write(print_to_screen, *) "Sound not be possible since all FOM found smaller than previous"
         write(print_to_screen, *) "should always be accepted."
         stop
    end if    
     
 ! ---------------- finished core mdmc part ------------------------ !    


    ! print out the best solution found. Note my_ps_best has a copy of the phase-space
    ! and the end point of that FOM calculation, not the start of it. Therefore the 
    ! re-calculated FOM although hopefully very simular will in general be different
    
    call restore_best_func_params(common_pe_list)
    call shallow_copy_phasespace(my_ps_best, my_ps)  
        
    call cal_time_corr_container(my_time_corr_container, my_ps, common_pe_list, c%md_per_time_bin, c%md_delta_t)   
    call cal_s_q_time(my_time_corr_container, my_ps%str, my_s_q_time)
    call print_s_q_time(my_s_q_time, density, c%temperature, "output/best_s_q_time.xml")
    call cal_s_q_omega(my_s_q_time, my_ps%str, my_s_q_omega)
    call print_s_q_omega(my_s_q_omega, density, c%temperature, "output/best_s_q_omega.xml")
    fom_val = func_val(my_s_q_omega, common_fom_list)
    call print_g_d(my_time_corr_container, product(my_ps%str%box_edges), size(my_ps%str%atoms), &
                   c%temperature, "output/best_g_d.xml")   
    call clear_time_corr_hist_container(my_time_corr_container)
      
    write(print_to_file,'(a,f14.5)') "Best FOM found = ", fom_best_including_rejected
    write(print_to_file,'(a,f14.5)') "Best FOM found (re-calculated from different phase-space point) = ", fom_val   
    write(print_to_file,'(a)') "With: "
    call print_all_func_params(print_to_file, common_pe_list)
    write(print_to_file, '(a,f14.5)') "Finished printing best fom. Time: ", toc()
    write(print_to_file, *) " "
    
    print *, ' '
    write(print_to_screen,'(a,f11.2,a)') 'Job took ', toc(), ' seconds to execute.'

    if (print_to_file /= 0) then
	    close(print_to_file)
    end if    
	  
    call xml_EndElement(xf, "mdmc-control-results")
    call xml_Close(xf)
    
  end subroutine run_mdmc_control_time_corr
  
  
  ! check to see if temperature is within percentage (where 1 is 100%) of t_target
  !
  function acceptable_temperature(t, t_target, percentage) result (yes_or_no)
    real (db), intent(in) :: t, t_target, percentage
    logical :: yes_or_no 
    
    yes_or_no = .true.
    
    if (abs((t - t_target)/t_target) > percentage) then
      yes_or_no = .false.
    end if
  end function acceptable_temperature
  
  
  ! check to see if energy is within percentage (where 1 is 100%) of e_target
  !
  function acceptable_energy(e, e_target, percentage) result (yes_or_no)
    real (db), intent(in) :: e, e_target, percentage
    logical :: yes_or_no 
    
    yes_or_no = .true.
    
    if (abs((e - e_target)/e_target) > percentage) then
      yes_or_no = .false.
    end if
  end function acceptable_energy 

end module mdmc_control_time_corr_class  