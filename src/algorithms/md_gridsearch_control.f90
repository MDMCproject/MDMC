module md_gridsearch_control_class
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

contains

  ! This is a test run-control subroutine - which assumes a known
  ! MD potential is used.
  
  ! This subroutine runs a series of MD runs with various Epot
  ! param values as defined on a grid.
  !
  ! It start off with an initial equilibrium run. Then a run
  ! to calculate g(r) (which may be an average of a number of
  ! of g(r) on a path) and at the end of this run the FOM. Then
  ! Epot param(s) are changed and a MD run to equilibrium is
  ! executed followed by the same run for calculating g(r) again
  ! and so on. 

  subroutine run_md_gridsearch_control(a_config, c)
    type (configuration), intent(inout) :: a_config
    type (mdmc_control_container) :: c
    
    real (db) :: time_now = 0.0    
    integer :: i, i_md, j
    real(db) :: density  ! used for cal input argument to save_rdf
    
    real(db) :: sum_kin_energy = 0.0
    real(db) :: temp_adjust_factor
  
    ! md storage containers

    type (phasespace) :: my_ps
    type (md_properties) :: my_props
    real(db) :: pressure_comp = 0.0, pot_energy = 0.0
    
    
    ! For printing out calculated rdf. Notice that the binning of 
    ! the rdf_printout_histogram is controlled by sub XML elements 
    ! of <control-object> and independent of the binning of any data
    
    type (rdf) :: rdf_printout
    type (histogram) :: rdf_printout_histogram 
    type (histogram) :: rdf_cal_histogram
    
    real (db), dimension(6) :: sigma_search  ! to hold what sigmas
                                             ! loop over
	
    integer :: print_to_file = 222  ! if set to zero then always
                                    ! print to screen only
    integer :: print_to_screen = 0
	  
    if (print_to_file /= 0) then
      open(print_to_file, file="output/job_log.txt")
    end if
	  		
    ! specify sigmas to loop over		
    sigma_search = (/ (2.2+0.4*float(i), i = 0, 5) /)				
	  		
	  		
    write(print_to_file,*) "In run_md_gridsearch_control"

    call tic

    ! initiate phasespace
    
    my_ps = make_phasespace(a_config%str, c%temperature)
                        
    
    ! to print out g(r) to file 
    
    rdf_printout_histogram = make_histogram(c%r_max, c%bin_length)
    
    rdf_cal_histogram = make_histogram(c%n_r_bin_cal, c%bin_length_cal)
                        
    rdf_printout = make_rdf(product(a_config%str%box_edges), size(a_config%str%atoms), &
                       floor(c%r_max/c%bin_length), c%bin_length)                    
                                               
                        
               
! -------------- settling the system ---------------- !

    if ( associated (common_pe_list%pt_lj_pe) ) then
      call set_func_param_val(common_pe_list%pt_lj_pe%params, "sigma", sigma_search(1))
      
      write(print_to_file, '(a,f12.4,a,f12.4)') "Starting Epot values sigma = ", &
        get_func_param_val(common_pe_list%pt_lj_pe%params, "sigma"), " and epsilon = ", &
        get_func_param_val(common_pe_list%pt_lj_pe%params, "epsilon")
    else
      print *, "ASSUME LJ POTENTIAL IN MD_GRIDSEARCH_CONTROL.F90"
      stop
    end if


    do i = 1, c%total_steps_initial_equilibration
      time_now = c%md_delta_t * i   ! perhaps print this one out 
      
      ! do one trajectory of length = 1 where pressure_comp and pot_energy is also
      ! calculated
      
      call trajectory_in_phasespace(my_ps, common_pe_list, 1, c%md_delta_t, & 
                                    pressure_comp, pot_energy)
      
      call md_cal_properties(my_ps, my_props, common_pe_list, pressure_comp, pot_energy)
      
      
      ! Optionally adjust the temperature

      if (i < c%total_step_temp_cali) then
        sum_kin_energy = sum_kin_energy + my_props%kin_energy%val
        if (mod(i,c%adjust_temp_at_interval) == 0) then
          ! see src/algorithms/md_gridsearch_control.doc for explanation of the expression below
          !
          temp_adjust_factor = sqrt( real(c%adjust_temp_at_interval,db) * 1.5 * c%temperature / &
              sum_kin_energy * (size(my_ps%str%atoms)-1.0) / real(size(my_ps%str%atoms),db))
          my_ps%p = my_ps%p * temp_adjust_factor
          sum_kin_energy = 0.0
        end if
      end if
   
        
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
      
      
    end do
    
    write(print_to_file, *) " "
    write(print_to_file, '(a,f12.4,a)') "Taken ", toc(), " seconds to execute initial equilibration."
    write(print_to_file, *) " "
    
 ! ---------------  finished settling the system -------------- !               
               
               
 ! --- calculate equilibrium g(r)'s for different Epot param values --- !              
               
    do i = 1, size(sigma_search)
    
      if ( associated (common_pe_list%pt_lj_pe) ) then
        call set_func_param_val(common_pe_list%pt_lj_pe%params, "sigma", sigma_search(i))
      
        write(print_to_file, '(a,f12.4,a,f12.4)') "Now Epot values sigma = ", &
          get_func_param_val(common_pe_list%pt_lj_pe%params, "sigma"), " and epsilon = ", &
          get_func_param_val(common_pe_list%pt_lj_pe%params, "epsilon")
        write(print_to_file, *) " "  
      else
        print *, "ASSUME LJ POTENTIAL IN MD_GRIDSEARCH_CONTROL.F90"
        stop
      end if
      
    
      if (i > 1) then
        ! do repeated MD equilibration
        write(print_to_file, '(a,f12.4)') "Start a repeated MD run with sigma = ", sigma_search(i)
        write(print_to_file, *) " "
        
        do i_md = 1, c%md_steps_repeated_equilibration

          call trajectory_in_phasespace(my_ps, common_pe_list, 1, c%md_delta_t, & 
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
          
          
          if (mod(i_md,c%average_over_repeated_equilibration) == 0) then 
            call md_print_properties(print_to_file, my_props)

            call md_reset_properties(my_props)
          end if
          
        end do
         
        write(print_to_file, '(a,f12.4)') "Finished repeated MD trajectory. Time: ", toc() 
        write(print_to_file, *) " "
        write(print_to_screen, '(a,f12.4)') "Finished repeated MD trajectory. Time: ", toc()
      end if       
      

      ! cal averaged rdf and FOM
      
      write(print_to_file, '(a,f12.4)') "Start cal FOM with sigma = ", sigma_search(i)
      write(print_to_file, *) " "
      
      do j = 1, c%average_over_this_many_rdf
        

        call trajectory_in_phasespace(my_ps, common_pe_list, c%cal_rdf_at_interval, c%md_delta_t)

       
        call accum_histogram(rdf_cal_histogram, my_ps%str)
        
        ! to print out rdf
        
        call accum_histogram(rdf_printout_histogram, my_ps%str)
        
      end do 
      
      write(print_to_file,'(a,f12.4)') "FOM = ", func_val(rdf_cal_histogram, common_fom_list)
      write(print_to_file, '(a,f12.4)') "Finished cal FOM. Time: ", toc()
      write(print_to_file, *) " "
      write(print_to_screen, '(a,f12.4)') "Finished cal FOM. Time: ", toc()
      
      
      ! print out rdf
      
      call cal_rdf(rdf_printout, rdf_printout_histogram)
            
      density = size(my_ps%str%atoms) / product(my_ps%str%box_edges)
      call save_rdf(rdf_printout, c%temperature, density)
      call clear_histogram(rdf_printout_histogram)
      
      call clear_histogram(rdf_cal_histogram)
          
    end do ! end of sigma search
                             
 
   
    !call save_structure(my_ps%str, "output/argon_structure.xml")
    
    print *, ' '
    print *, 'Job took ', toc(), ' seconds to execute.'
    
    if (print_to_file /= 0) then
	    close(print_to_file)
	  end if    
    
  end subroutine run_md_gridsearch_control

end module md_gridsearch_control_class
