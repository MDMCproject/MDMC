module md_control_class
use configuration_class
use gpe_class
use common_block_class, only : common_config, common_pe_list
!use common_potential_block_class
use phasespace_class
use md_properties_class
use control_containers_class
use tic_toc_class

  implicit none

contains

  subroutine run_md_control(a_config, c)
    type (configuration), intent(inout) :: a_config
    type (md_control_container) :: c
    
    real (db) :: time_now = 0.0    
    integer :: i
    
    real(db) :: sum_kin_energy = 0.0
    real(db) :: temp_adjust_factor
  
    ! md storage containers

    type (phasespace) :: my_ps
    type (md_properties) :: my_props
    real(db) :: pressure_comp = 0.0, pot_energy = 0.0
		
		
    write(*,*) "In run_md_control"

    call tic

    ! initiate phasespace
    
    my_ps = make_phasespace(a_config%str, c%r_cut, c%delta_r, &
                        c%temperature)
    
    
    do i = 1, c%step_limit
      time_now = c%time_step * i   ! perhaps print this one out 
      
      ! do one trajectory of length = 1 where pressure_comp and pot_energy is also
      ! calculated
      
      !call trajectory_in_phasespace(my_ps, common_pe_list, 1, c%time_step)
      call trajectory_in_phasespace(my_ps, common_pe_list, 1, c%time_step, & 
                                    pressure_comp, pot_energy)
      
      !call md_cal_properties(my_ps, my_props, common_pe_list)
      call md_cal_properties(my_ps, my_props, common_pe_list, pressure_comp, pot_energy)
      
      ! case you want to adject the temperature in the initial stages of the MD simulation
        
      if (c%perform_initial_temperature_calibration) then
        if (i < c%total_step) then
          sum_kin_energy = sum_kin_energy + my_props%kin_energy%val
          if (mod(i,c%adjust_temp_at_interval) == 0) then
            temp_adjust_factor = sqrt(c%adjust_temp_at_interval * 1.5 * c%temperature / &
               sum_kin_energy * (size(my_ps%str%atoms)-1.0) / size(my_ps%str%atoms))
            !write(*,*) sum_kin_energy / c%adjust_temp_at_interval   
            my_ps%p = my_ps%p * temp_adjust_factor
            !write(*,*) sum(my_ps%p * my_ps%p) * 0.5 / size(my_ps%str%atoms)
            sum_kin_energy = 0.0
          end if
        end if
      end if      
      

      ! accumulate the calculated MD property values
        
      call md_accum_properties(my_props)
        
        
      ! print out stuff and interval = average_over_this_many_steps
        
      if (mod(i,c%average_over_this_many_step) == 0) then
        call md_print_properties(my_props)
        write(*,'(a,3f12.6)') "total momentum ", sum(my_ps%p,1)
        call md_reset_properties(my_props)
      end if
      
    end do
    
    print *, ' '
    print *, 'Job took ', toc(), ' seconds to execute.'
  end subroutine run_md_control

end module md_control_class
