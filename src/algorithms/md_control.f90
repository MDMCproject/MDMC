module md_control_class
use configuration_class
use function_class
use common_block_class, only : common_config, common_pe_list
use phasespace_class
use md_properties_class
use control_containers_class
use tic_toc_class
use rdf_class
use structure_reader_class

  implicit none

contains

  subroutine run_md_control(a_config, c)
    type (configuration), intent(inout) :: a_config
    type (md_control_container) :: c
    
    real (db) :: time_now = 0.0    
    integer :: i
    real(db) :: density  ! used for cal input argument to save_rdf
    
    real(db) :: sum_kin_energy = 0.0
    real(db) :: temp_adjust_factor
    
    ! to calculate pressure correction
    
    real(db) :: pressure_corr
    real(db) :: sig_r_cut3
    real(db) :: sigma, epsilon
  
  
    ! md storage containers

    type (phasespace) :: my_ps
    type (md_properties) :: my_props
    real(db) :: pressure_comp = 0.0, pot_energy = 0.0
    type (rdf) :: my_rdf
    type (histogram) :: my_histogram
    
    integer :: print_to_file = 111
    integer :: print_to_screen = 0
	  
	  if (print_to_file /= 0) then
	    open(print_to_file, file="output/job_log.txt")
	  end if
	  
		
    write(print_to_file, *) "In run_md_control"


    call tic

    ! initiate
    
    my_ps = make_phasespace(a_config%str, c%temperature)                      
    
    my_histogram = make_histogram(c%r_max, c%bin_length)
                        
    my_rdf = make_rdf(product(a_config%str%box_edges), size(a_config%str%atoms), &
                      floor(c%r_max/c%bin_length), c%bin_length)
 
                          
    
    do i = 1, c%step_limit
      time_now = c%md_delta_t * i   ! perhaps print this one out 
      
      ! do one trajectory of length = 1 where pressure_comp and pot_energy is also
      ! calculated
      
      call trajectory_in_phasespace(my_ps, common_pe_list, 1, c%md_delta_t, & 
                                    pressure_comp, pot_energy)
      
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
     
        
      ! print out stuff and interval = average_over_this_many_steps
        
      if (mod(i,c%average_over_this_many_step) == 0) then 
        call md_print_properties(print_to_file, my_props)
        
        if (sum(sum(my_ps%p,1)) > 0.0001) then
          write(*,*) "ERROR:"
          write(*,'(a,3f12.6)') "total momentum ", sum(my_ps%p,1)
          write(*,*) "Serious problem - total momentum different from zero"
          stop
        end if
        
        
        ! Print out pressure correction term 
        ! Notice, the calculation is limited for simplicity to the case where
        ! we use the neighest nearbour method since some r_cut value is needed
        ! to cal sig_r_cut3 below. For the case where nn-method not use we could
        ! perhaps define a r_cut value to L/2 or something like that
        
        if ( associated (common_pe_list%pt_lj_pe) ) then
          if (a_config%str%nn_list%r_cut /= 0.0) then
            density = size(my_ps%str%atoms)/product(my_ps%str%box_edges)
            sigma = get_func_param_val(common_pe_list%pt_lj_pe%params, "sigma")  
		        epsilon = get_func_param_val(common_pe_list%pt_lj_pe%params, "epsilon")
            sig_r_cut3 = (sigma / a_config%str%nn_list%r_cut)**3
            pressure_corr = 32.0 * pi_value* density**2 * sigma**3 * epsilon * &
                            (sig_r_cut3**3 - 1.5*sig_r_cut3) / 9.0
                            
            ! convert to atm
            
            pressure_corr = pressure_corr * P_unit
            
            write(print_to_file,'(a,f12.6)') "P_corr(atm) = ", pressure_corr
          end if
        end if
        
        
        call md_reset_properties(my_props)
        write(print_to_file, '(a,i8,a,f12.4,a)') "MD steps = ", i, " MD run-time = ", time_now, "*10e-13"
        write(print_to_screen, '(a,i8,a,f12.4,a)') "MD steps = ", i, " MD run-time = ", time_now, "*10e-13"
      end if
      
      
      ! if rdf is set to be calculated
      
      if (c%calculate_rdf) then
        if (i > c%total_step_temp_cali) then
          if ( mod(i-c%total_step_temp_cali, c%cal_rdf_at_interval) == 0 ) then
            call accum_histogram(my_histogram, my_ps%str)
          end if
          
          
          ! when average_over_this_many_rdf rdf's have been calculated then print out the average
          ! of these and reset rdf sum container
          
          if ( mod(i-c%total_step_temp_cali, c%cal_rdf_at_interval*c%average_over_this_many_rdf) == 0 ) then
            
            call cal_rdf(my_rdf, my_histogram)
            
            density = size(my_ps%str%atoms) / product(my_ps%str%box_edges)
            call save_rdf(my_rdf, c%temperature, density)
            
            call clear_histogram(my_histogram)
          end if
        
        end if
        
        
      end if
      
    end do
   
    call save_structure(my_ps%str, "output/argon_structure.xml")
    
    write(print_to_file, *) " "
    write(print_to_file, '(a,f12.4,a)') "Job took ", toc(), " seconds to execute."
    
    
    if (print_to_file /= 0) then
	    close(print_to_file)
	  end if

  end subroutine run_md_control

end module md_control_class
