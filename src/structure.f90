module structure_class

implicit none


  type atom
    character(len=2)  :: element_type
    double precision, dimension(3) :: position
  end type atom

  type structure
    type (atom), dimension(:), allocatable :: atoms
    double precision :: half_cell_edge
    character(len=120) :: title
  end type structure


contains


  subroutine swap_atoms(a_structure, atom_number1, atom_number2)
    type (structure), intent(inout) :: a_structure
    integer, intent(in) :: atom_number1, atom_number2

  end subroutine swap_atoms

end module structure_class
