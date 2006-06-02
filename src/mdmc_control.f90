module mdmc_control_class
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

  implicit none

  public :: run_mdmc_control
  
  private :: acceptable_temperature
  private :: acceptable_energy  

contains

  subroutine run_mdmc_control(a_config, c)
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

    real(db) :: fom_best  ! used for keeping track of best fom
    type (phasespace) :: my_ps_best ! to store status of best phasespace

    type (phasespace) :: my_ps, my_ps_old
    type (md_properties) :: my_props
    real(db) :: pressure_comp = 0.0, pot_energy = 0.0
    type (rdf) :: my_rdf
    type (histogram) :: my_histogram
	
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
      open(print_to_file, file="output/job_summary.txt")
    end if
	  		
    write(print_to_file,*) "In run_md_control"

    call tic

    ! initiate phasespace
    
    my_ps = make_phasespace(a_config%str, c%temperature)
    
    my_ps_old = copy_phasespace(my_ps)
    
    my_ps_best = copy_phasespace(my_ps)
                        
                        
    ! to print out g(r) to file (otherwise neither my_histogram nor my_rdf needed)
    
    my_histogram = make_histogram(c%r_max, c%bin_length)
                        
    my_rdf = make_rdf(product(a_config%str%box_edges), size(a_config%str%atoms), &
                       floor(c%r_max/c%bin_length), c%bin_length)                    
                                               
                                               
                        
               
! -------------- initial equilibration ---------------- !

    sum_kin_energy = 0.0
    
    do i = 1, c%total_steps_initial_equilibration
      time_now = c%time_step * i   ! perhaps print this one out 
      
      ! do one trajectory of length = 1 where pressure_comp and pot_energy is also
      ! calculated
      
      !call trajectory_in_phasespace(my_ps, common_pe_list, 1, c%time_step)
      call trajectory_in_phasespace(my_ps, common_pe_list, 1, c%time_step, & 
                                    pressure_comp, pot_energy)
      
      !call md_cal_properties(my_ps, my_props, common_pe_list)
      call md_cal_properties(my_ps, my_props, common_pe_list, pressure_comp, pot_energy)
      
      ! case you want to adject the temperature in the initial stages of the MD simulation
      ! (notice c%total_step_temp_cali = 0 if <perform-initial-temperature-calibration> 
      ! element not specified in input file)

      if (i < c%total_step_temp_cali) then
        sum_kin_energy = sum_kin_energy + my_props%kin_energy%val
        if (mod(i,c%adjust_temp_at_interval) == 0) then
          temp_adjust_factor = sqrt(c%adjust_temp_at_interval * 1.5 * c%temperature / &
              sum_kin_energy * (size(my_ps%str%atoms)-1.0) / size(my_ps%str%atoms))
          my_ps%p = my_ps%p * temp_adjust_factor
          sum_kin_energy = 0.0
        end if
      end if
     
      

      ! accumulate the calculated MD property values
        
      call md_accum_properties(my_props)
        
        
      ! print out stuff and interval = average_over_this_many_steps
        
      if (mod(i,c%average_over_initial_equilibration) == 0) then 
        call md_print_properties(print_to_file, my_props)
        
        if (sum(sum(my_ps%p,1)) > 0.0001) then
          write(*,*) "ERROR:"
          write(*,'(a,3f12.6)') "total momentum ", sum(my_ps%p,1)
          write(*,*) "Serious problem - total momentum different from zero"
          stop
        end if
        
        call md_reset_properties(my_props)
        write(*, '(a,i8,a,f12.4,a)') "MD steps = ", i, " MD run-time = ", time_now, "*10e-13"
      end if
      
      
      ! Store the average energy at the point when finished the temperature calibration.
      ! Note this must be done after md_print_properties has been called - since only this
      ! subfunction alters the my_props.ave value.
      
      if (i == c%total_step_temp_cali) then
        average_energy_end_of_temp_calibration = my_props.tot_energy.ave 
      end if

    end do
    
    write(print_to_file, *) " "
    write(print_to_file, '(a,f12.4,a)') "Taken ", toc(), " seconds to execute initial equilibration."
    write(print_to_file, *) " "


    ! Determine if equilibrium was reached
    !
    ! Note the (2.0/ndim) factor is to convert from dimensionless kin_energy per atom to 
    ! dimensionless temperature
    
    if ( acceptable_temperature((2.0/ndim)*my_props.kin_energy.ave, &
         c%temperature) == .false.) then
         write(print_to_screen, *) "Initial equilibration did not reach equilibrium"
         write(print_to_screen, *) "Temperature outside acceptable value - STOP"
         stop
    end if   
    
    if ( acceptable_energy(average_energy_end_of_temp_calibration, &
         my_props.tot_energy.ave) == .false.) then
         write(print_to_screen, *) "Initial equilibration did not reach equilibrium - STOP"
         write(print_to_screen, *) "Energy outside acceptable value - STOP"
         stop
    end if     
    
 ! ---------------  finished initial equilibration -------------- !      
 
 
 
 ! ----------- calculate first FOM ------------- !        
    
    do j = 1, c%average_over_this_many_rdf

      call trajectory_in_phasespace(my_ps, common_pe_list, c%cal_rdf_at_interval, c%time_step)
      
      call func_accum_histogram(my_ps%str, common_fom_list)
      
      ! to print out rdf
      
      call accum_histogram(my_histogram, my_ps%str)
      
    end do 
    
    !! print also out what the param values are
    fom_val = func_val(my_ps%str, common_fom_list)
    call func_clear_histogram(my_ps%str, common_fom_list)
    write(print_to_file,'(a,f12.4)') "1st FOM = ", fom_val
    write(print_to_file, '(a,f12.4)') "Finished cal 1st FOM. Time: ", toc()
    write(print_to_file, *) " "
    write(print_to_screen, '(a,f12.4)') "Finished cal 1st FOM. Time: ", toc()
    
    
    call xml_NewElement(xf, "accept")
    call add_xml_attribute_func_params(xf, common_pe_list)
    call xml_AddAttribute(xf, "val", str(fom_val, format="(f10.5)"))
    call xml_EndElement(xf, "accept")    


    ! to print out rdf 
      
    call cal_rdf(my_rdf, my_histogram)
            
    density = size(my_ps%str%atoms) / product(my_ps%str%box_edges)
    call save_rdf(my_rdf, c%temperature, density)
    call clear_histogram(my_histogram)


    ! store best solution so far
    
    fom_best = fom_val
    call backup_best_func_params(common_pe_list)
    call shallow_copy_phasespace(my_ps, my_ps_best)



 ! ----------- finished calculating first FOM ------------- !




 ! ---------------------- mdmc part ------------------------ !              
               
  
    ! save state
    
    fom_old = fom_val
    call backup_func_params(common_pe_list)
    call shallow_copy_phasespace(my_ps, my_ps_old)       
      
              
    do i = 1, c%mc_steps

      write (print_to_screen, *) "Begin MC step number ", i
      
      write (print_to_file, *) "Begin MC step number ", i
      write (print_to_file, *) " "

      call move_random_func_params(common_pe_list)
      

      ! do MD equilibration
      
      sum_kin_energy = 0.0
      
      do i_md = 1, c%md_steps_repeated_equilibration

        call trajectory_in_phasespace(my_ps, common_pe_list, 1, c%time_step, & 
                                  pressure_comp, pot_energy)

        call md_cal_properties(my_ps, my_props, common_pe_list, pressure_comp, pot_energy)       
      
      
        if (i_md < c%total_step_temp_cali_repeated) then
          sum_kin_energy = sum_kin_energy + my_props%kin_energy%val
          if (mod(i_md,c%adjust_temp_at_interval_repeated) == 0) then
            temp_adjust_factor = sqrt(c%adjust_temp_at_interval_repeated * 1.5 * c%temperature / &
                sum_kin_energy * (size(my_ps%str%atoms)-1.0) / size(my_ps%str%atoms))
            my_ps%p = my_ps%p * temp_adjust_factor
            sum_kin_energy = 0.0
          end if
        end if        
      
      
        call md_accum_properties(my_props)
        
        if (mod(i_md,c%average_over_repeated_equilibration) == 0) then 
          call md_print_properties(print_to_file, my_props)

          call md_reset_properties(my_props)
        end if
        
        ! Store the average energy at the point when finished the temperature calibration.
        ! Note this must be done after md_print_properties has been called - since only this
        ! subfunction alters the my_props.ave value.
      
        if (i_md == c%total_step_temp_cali_repeated) then
          average_energy_end_of_temp_calibration = my_props.tot_energy.ave 
        end if
              
      end do
         
      write(print_to_file, '(a,f12.4)') "Finished repeated MD trajectory. Time: ", toc() 
      write(print_to_file, *) " "
      write(print_to_screen, '(a,f12.4)') "Finished repeated MD trajectory. Time: ", toc()      
      

      ! Determine if equilibrium was reached
      !
      ! Note the (2.0/ndim) factor is to convert from dimensionless kin_energy per atom to 
      ! dimensionless temperature
      
      if ( acceptable_temperature((2.0/ndim)*my_props.kin_energy.ave, &
          c%temperature) == .false. .or. acceptable_energy(average_energy_end_of_temp_calibration, &
          my_props.tot_energy.ave) == .false.) then
        write(print_to_screen, *) "Did not reach equilibrium in MD run"
        write(print_to_file, *) "Did not reach equilibrium in MD run"

        call xml_NewElement(xf, "failed-md-run")
        call add_xml_attribute_func_params(xf, common_pe_list)
        call xml_EndElement(xf, "failed-md-run")   

        call restore_func_params(common_pe_list)
        call shallow_copy_phasespace(my_ps_old, my_ps)         

        cycle
      end if   


      ! cal averaged rdf and FOM
      
      do j = 1, c%average_over_this_many_rdf

        call trajectory_in_phasespace(my_ps, common_pe_list, c%cal_rdf_at_interval, c%time_step)
       
        call func_accum_histogram(my_ps%str, common_fom_list)
        
        ! to print out rdf
        
        !if (i == c%mc_steps) call accum_histogram(my_histogram, my_ps%str)
        
      end do 
      
      !! print also out what the param values are
      fom_val = func_val(my_ps%str, common_fom_list)
      call func_clear_histogram(my_ps%str, common_fom_list)
      write(print_to_file,'(a,f12.4)') "FOM = ", fom_val
      write(print_to_file, '(a,f12.4)') "Finished cal FOM. Time: ", toc()
      write(print_to_file, *) " "
      write(print_to_screen, '(a,f12.4)') "Finished cal FOM. Time: ", toc()
      
      
      ! check if new best fom
      
      if (fom_val < fom_best) then
        fom_best = fom_val
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
        call random_number(ran_num)
          
        if(exp(- delta_fom / c%temperature_mc) > ran_num) then
          accept_parameters = .true.
        end if
      end if
 
        
      if (accept_parameters) then
      
        call xml_NewElement(xf, "accept")
        call add_xml_attribute_func_params(xf, common_pe_list)
        call xml_AddAttribute(xf, "val", str(fom_val, format="(f10.5)"))
        call xml_EndElement(xf, "accept")       
     
      
        ! save state
    
        fom_old = fom_val
        call backup_func_params(common_pe_list)
        call shallow_copy_phasespace(my_ps, my_ps_old) 
        
      else
      
        call xml_NewElement(xf, "rejected")
        call add_xml_attribute_func_params(xf, common_pe_list)
        call xml_AddAttribute(xf, "val", str(fom_val, format="(f10.5)"))
        call xml_EndElement(xf, "rejected")
        
        
        ! restore state
    
        call restore_func_params(common_pe_list)
        call shallow_copy_phasespace(my_ps_old, my_ps)          
              
      end if  
      
    end do
     
 ! ---------------- finished mdmc part ------------------------ !      
     
    
     
    ! print out rdf for best fom
    
    call restore_best_func_params(common_pe_list)
    call shallow_copy_phasespace(my_ps_best, my_ps)
        
    do i = 1, 5    
        
    do j = 1, c%average_over_this_many_rdf

      call trajectory_in_phasespace(my_ps, common_pe_list, c%cal_rdf_at_interval, c%time_step)
       
      call func_accum_histogram(my_ps%str, common_fom_list)
        
      ! to print out rdf
        
      call accum_histogram(my_histogram, my_ps%str)
        
    end do 
      
    fom_val = func_val(my_ps%str, common_fom_list)
    call func_clear_histogram(my_ps%str, common_fom_list)
    write(print_to_file,'(a,f12.4)') "BEST FOM = ", fom_val
    write(print_to_file,'(a)') "WITH: "
    call print_all_func_params(print_to_file, common_pe_list)
    write(print_to_file, '(a,f12.4)') "Finished cal FOM. Time: ", toc()
    write(print_to_file, *) " "
      
    call cal_rdf(my_rdf, my_histogram)
            
    density = size(my_ps%str%atoms) / product(my_ps%str%box_edges)
    call save_rdf(my_rdf, c%temperature, density)
    call clear_histogram(my_histogram)
    
    end do   
                      
    !call save_structure(my_ps%str, "output/argon_structure.xml")
    
    print *, ' '
    print *, 'Job took ', toc(), ' seconds to execute.'
    
    if (print_to_file /= 0) then
	    close(print_to_file)
	  end if    
	  
	  
	  call xml_EndElement(xf, "mdmc-control-results")
    call xml_Close(xf)
    
  end subroutine run_mdmc_control
  
  
  ! check to see if temperature is within certain limits of t_target
  function acceptable_temperature(t, t_target) result (yes_or_no)
    real (db), intent(in) :: t, t_target
    logical :: yes_or_no 
    
    yes_or_no = .true.
    
    if (abs((t - t_target)/t_target) > 0.2) then
      yes_or_no = .false.
    end if
  end function acceptable_temperature
  
  ! check to see if energy is within certain limits of t_target
  function acceptable_energy(e, e_target) result (yes_or_no)
    real (db), intent(in) :: e, e_target
    logical :: yes_or_no 
    
    yes_or_no = .true.
    
    if (abs((e - e_target)/e_target) > 0.1) then
      yes_or_no = .false.
    end if
  end function acceptable_energy

end module mdmc_control_class
