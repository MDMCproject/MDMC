module structure_class
use various_constants_class

implicit none


  type atom
    character(len=2)  :: element_type
    real(db), dimension(ndim) :: r   ! cartesian coor in units AA
    real(db) :: mass=1.0 !39.95        ! in units of amu
    real(db) :: inv_mass=1.0 !0.025    ! to save comp time
  end type atom

  type structure
    type (atom), dimension(:), allocatable :: atoms
    
    ! box_edges in units of AA. Atoms are forces into a box with dimensions
    ! defined in box_edges and with the box centered around the origin. That 
    ! means e.g. along the x-direction the atoms are placed
    ! within -box_edge(1)/2 and box_edge(1)/2
    real(db), dimension(ndim) :: box_edges  
                           
    real(db) :: density    ! in units of atom/AA3
    character(len=120) :: title
  end type structure

contains

  function copy_structure(str_in) result(str_out)
    type (structure), intent(in) :: str_in
    type (structure) :: str_out
    
    integer :: i
    real(db) :: sum_out, sum_in
    
    str_out%box_edges = str_in%box_edges
    str_out%density = str_in%density
    str_out%title = str_in%title
    
    allocate(str_out%atoms(size(str_in%atoms)))
    
    sum_out = 0.0
    sum_in = 0.0
    do i = 1, size(str_in%atoms)
      str_out%atoms(i)%r = str_in%atoms(i)%r
      str_out%atoms(i)%mass = str_in%atoms(i)%mass
      str_out%atoms(i)%inv_mass = str_in%atoms(i)%inv_mass
      str_out%atoms(i)%element_type = str_in%atoms(i)%element_type
      
      sum_out = sum_out + sum(str_out%atoms(i)%r)
      sum_in = sum_in + sum(str_in%atoms(i)%r)
    end do
    
    if (sum_out /= sum_in) then
      write(*,*) "ERROR: problem with copying structure"
      stop
    end if
  
  end function copy_structure

  subroutine swap_atoms(a_structure, atom_number1, atom_number2)
    type (structure), intent(inout) :: a_structure
    integer, intent(in) :: atom_number1, atom_number2

  end subroutine swap_atoms

end module structure_class
