module mc_control_class
use configuration_class
use function_class
use common_block_class, only : common_config, common_pe_list
!use common_potential_block_class
use phasespace_class
use md_properties_class
use control_containers_class
use tic_toc_class
use rdf_class
use structure_reader_class

  implicit none

contains

  subroutine run_mc_control(a_config, c)
    type (configuration), intent(inout) :: a_config
    type (md_control_container) :: c
    
   ! PSUDOCODE ONLY FOR NOW
   
   ! slow version
   ! move_one_atom(config_new)
   ! is_valid_structure(config_new)
   ! func_val(str, list)
   
   ! fast version
   ! is_valid_atom_move(struct_new,which_atom,one_atom_dists)
   ! nn_update_histogram_one_atom(....)
   ! func_val_nn(str, list, nn)

  end subroutine run_mc_control

end module mc_control_class
