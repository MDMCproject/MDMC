module md_control_class
use configuration_class
use gpe_class
use common_block_class, only : common_config
use common_potential_block_class
use phasespace_class
use md_properties_class

  implicit none

contains

  subroutine run_md_control(a_config, r_cut, delta_r)
    type (configuration), intent(inout) :: a_config
    real (db), intent(in) :: r_cut, delta_r

		real (db) :: delta_t, temperature
    real (db) :: time_now = 0.0
    integer :: n_md = 100
    integer :: i
    integer :: steps_to_equilibrium = 0
    integer :: adjust_temp_after_this_many_steps = 100
    integer :: average_over_this_many_steps = 20
    
    type (phasespace) :: my_ps
    type (md_properties) :: my_props
    
    real(db) :: pressure_comp, pot_energy
		
		
    write(*,*) "In run_md_control"

    call md_reset_properties(my_props)
    
    delta_t = 0.005
    temperature = 1.0
    
    my_ps = make_phasespace(a_config%str, r_cut, delta_r, &
                        temperature)
    
    do i = 1, n_md
      time_now = delta_t * i
      
      call trajectory_in_phasespace_extra(my_ps, 1, delta_t, pressure_comp, pot_energy)

      ! calculate e_kin, e_tot, pressure
      call md_cal_properties_extra(my_ps, my_props, common_gpe, pressure_comp, pot_energy)
      if (i < steps_to_equilibrium) then
        
      else
        call md_accum_properties(my_props)
        if (mod(i,average_over_this_many_steps) == 0) then
          call md_print_properties(my_props)
          call md_reset_properties(my_props)
        end if
      end if
      
    end do
  end subroutine run_md_control

end module md_control_class
