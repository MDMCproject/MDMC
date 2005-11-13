module md_control_class
use configuration_class
use gpe_class
use common_block_class, only : common_config
use common_potential_block_class
use phasespace_class

  implicit none

contains

  subroutine run_md_control(a_config, r_cut, delta_r)
    type (configuration), intent(inout) :: a_config
    real (db), intent(in) :: r_cut, delta_r

		real (db) :: pot_energy, delta_t, temperature
    real (db) :: time_now = 0.0
    integer :: n_md = 10
    integer :: i
    integer :: steps_to_equilibrium = 0
    integer :: adjust_temp_after_this_many_steps = 100
    integer :: average_over_this_many_steps = 2000
    
    type (phasespace) :: my_ps
		
    write(*,*) "In run_md_control"

    delta_t = 0.005
    temperature = 1.0
    
    my_ps = make_phasespace(a_config%str, r_cut, delta_r, &
                        temperature)
    
    do i = 1, n_md
      time_now = delta_t * i
      
      call trajectory_in_phasespace(my_ps, 1, delta_t)

      ! calculate e_kin, e_tot, pressure
      
      !if (i < steps_to_equilibrium) then
        
      !else
        ! accum properties
      !  if (mod(i,average_over_this_many_steps) == 0) then
          ! print mean s.d. and reset properties
      !  end if
      !end if
      
    end do
    ! md_cal_properties(phasespace..., props)
    ! md_print_properties(props)
    
  end subroutine run_md_control

end module md_control_class
