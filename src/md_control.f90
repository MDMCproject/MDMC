module md_control_class
use configuration_class
use gpe_class
use common_block_class, only : common_config
use phasespace_class

  implicit none

contains

  subroutine run_md_control(a_config, r_cut, delta_r)
    type (configuration), intent(inout) :: a_config
    real (db), intent(in) :: r_cut, delta_r

		real (db) :: pot_energy
    type (phasespace) :: phasespace_backup
		
    write(*,*) "In run_md_control"

    
   phasespace_backup = make_phasespace(a_config%str, r_cut, delta_r)
		
		pot_energy = gpe_val(common_config%str, &
		                     common_pe_list(1)%vars)
                         
    write (*, '(a,f10.3)') "Epot = ", pot_energy
  end subroutine

end module md_control_class
