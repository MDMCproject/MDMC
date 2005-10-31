module structure_class
use various_constants_class

implicit none


  type atom
    character(len=2)  :: element_type
    real(db), dimension(3) :: r   ! cartesian coor in units AA
    real(db) :: mass=39.95        ! in units of amu
    real(db) :: inv_mass=0.025    ! to save comp time
  end type atom

  type structure
    type (atom), dimension(:), allocatable :: atoms
    
    ! box_edges in units of AA. Atoms are forces into a box with dimensions
    ! defined in box_edges and with the box centered around the origin. That 
    ! means e.g. along the x-direction the atoms are placed
    ! within -box_edge(1)/2 and box_edge(1)/2
    real(db), dimension(3) :: box_edges  
                           
    real(db) :: density    ! in units of atom/AA3
    character(len=120) :: title
  end type structure

contains


  subroutine swap_atoms(a_structure, atom_number1, atom_number2)
    type (structure), intent(inout) :: a_structure
    integer, intent(in) :: atom_number1, atom_number2

  end subroutine swap_atoms

end module structure_class
