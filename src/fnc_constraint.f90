module fnc_constraint_class
use structure_class

implicit none


  type fnc_group
    integer :: atomindex
    integer, dimension(:), pointer :: neighbour_index
    integer, dimension(:), pointer :: constraint_type_id
  end type fnc_group

  type fnc_constraint
    type (fnc_group), dimension(:), pointer :: fnc_groups
    real(db), dimension(:), pointer :: rmin
    real(db), dimension(:), pointer :: rmax
  end type fnc_constraint

contains

  function validate_fnc_constraint(a_structure) result (flag)
    type (structure), intent(in)  :: a_structure
    logical  :: flag

  end function validate_fnc_constraint

  function make_fnc_constraint(filename) result (new_fnc_constraint)
    character(*)  :: filename
    type (fnc_constraint) :: new_fnc_constraint
 
  end function make_fnc_constraint

end module fnc_constraint_class
