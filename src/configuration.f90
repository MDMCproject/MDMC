module configuration_class
!use structure_class
use fnc_constraint_class

  implicit none


  ! Holds pointers to all possible constraint types

  type constraints
    private
    type (fnc_constraint), pointer :: pt_fnc_constraint => null()
  !  type (cutoff_constraint), pointer :: pt_cutoff_constraint 
  end type constraints


  ! This type holds a structure and constraints that are associated
  ! with this structure. Note that, the structure are not required
  ! to satisfy these constraints, but you may ask the structure IF
  ! it does satisfy them  

  type configuration
    type (structure) :: cf_structure
    type (constraints) :: cf_constraints
  end type configuration


  interface add_constraint
    module procedure add_fnc_constraint
  end interface

  !type (configuration)  :: my_configuration

contains


  subroutine add_fnc_constraint(a_config, fnc_target)
    type (configuration), intent(inout) :: a_config
    type (fnc_constraint), target, intent(in) :: fnc_target
    
    a_config%cf_constraints%pt_fnc_constraint => fnc_target

  end subroutine add_fnc_constraint


  function is_valid_structure(a_config) result(flag)
    type (configuration), intent(in) :: a_config
    logical :: flag

    if ( associated (a_config%cf_constraints%pt_fnc_constraint) ) then
      flag = validate_fnc_constraint (a_config%cf_structure)
    !  if (.not. flag) exit this functions
    !else if
    !  etc....
    end if
  end function is_valid_structure


end module configuration_class
