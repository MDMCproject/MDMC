module md_control_class
use configuration_class
use gpe_class
use common_block_class, only : common_configuration

  implicit none

contains

  subroutine run_md_control(a_config)
    type (configuration), intent(inout) :: a_config

		real (db) :: pot_energy
		
    write(*,*) "In run_moldyn_control"

		
		pot_energy = gpe_val(common_configuration%cf_structure, &
		                     common_pe_list(1)%vars)
                         
    write(*,*) "pot_energy ", pot_energy
  end subroutine

end module md_control_class
