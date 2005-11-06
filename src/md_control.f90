module md_control_class
use configuration_class
use gpe_class
use common_block_class, only : common_config

  implicit none

contains

  subroutine run_md_control(a_config)
    type (configuration), intent(inout) :: a_config

		real (db) :: pot_energy
		
    write(*,*) "In run_moldyn_control"

		
		pot_energy = gpe_val(common_config%str, &
		                     common_pe_list(1)%vars)
                         
    write (*, '(a,f10.3)') "Epot = ", pot_energy
  end subroutine

end module md_control_class
