module md_control_class
use configuration_class

  implicit none

contains

  subroutine run_md_control(a_config)
    type (configuration), intent(inout) :: a_config

    write(*,*) "In run_moldyn_control"
  end subroutine

end module md_control_class
