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
    integer :: n_md
    type (phasespace) :: phasespace_backup, ps_temp
		
    write(*,*) "In run_md_control"

    n_md = 1
    delta_t = 0.005
    temperature = 1.0
    
    phasespace_backup = make_phasespace(a_config%str, r_cut, delta_r, &
                        temperature)
    
    call trajectory_in_phasespace(phasespace_backup, n_md, delta_t)

    ! md_cal_properties(phasespace..., props)
    ! md_print_properties(props)
    
  end subroutine run_md_control

end module md_control_class
