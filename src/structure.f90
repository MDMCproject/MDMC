module structure_class
use various_constants_class

implicit none

  public :: apply_boundary_condition
  public :: copy_structure
  public :: swap_atoms

  type atom
    character(len=2)  :: element_type
    real(db) :: mass=1.0 !39.95        ! in units of amu
    real(db) :: inv_mass=1.0 !0.025    ! to save comp time
  end type atom

  type structure
    type (atom), dimension(:), allocatable :: atoms
    
    ! Fortran90 have many utilities for handling arrays hence
    ! the atom coordinates are here stored in ONE array rather
    ! than as atoms%r(1:ndim) for each atom as was originally done
    !
    ! Index 1 : atom number
    ! Index 2 : atom coordinates
    
    real(db), dimension(:,:), allocatable :: r  ! cartesian coor in units AA
    
    ! box_edges in units of AA. Atoms are forces into a box with dimensions
    ! defined in box_edges and with the box centered around the origin. That 
    ! means e.g. along the x-direction the atoms are placed
    ! within -box_edge(1)/2 and box_edge(1)/2
    real(db), dimension(ndim) :: box_edges  
                           
    real(db) :: density    ! in units of atom/AA3
    character(len=120) :: title
  end type structure

contains

  ! apply periodic boundary conditions, however assume here that atoms have
  ! not moved out more that box_edge length away from the region
  ! -box_edge/2 : box_edge/2
  subroutine apply_boundary_condition(str)
    type (structure), intent(inout) :: str
    
    integer :: i_a, i
    
    do i_a = 1, size(str%atoms)
      
      do i = 1, ndim
        if (str%r(i_a,i) >= 0.5 * str%box_edges(i)) then
          str%r(i_a,i) = str%r(i_a,i) - str%box_edges(i)
        end if
        if (str%r(i_a,i) < -0.5 * str%box_edges(i)) then
          str%r(i_a,i) = str%r(i_a,i) + str%box_edges(i)
        end if       
      end do
      
      
      ! check that atoms have not moved further than one box length!!
      ! This might be caused by 1) crap initial structure where two atoms starts
      ! off being far too close to each other 2) when using the nearest neighbour
      ! method be careful not to wait too long before updating which can cause
      ! the effect described in 1)
      
      do i = 1, ndim
        if (str%r(i_a,i) >= 0.5 * str%box_edges(i)) then
          write (*,*) "atoms move further than the length of the box within one step!!"
          write (*,*) "ERROR written from apply_boundary_condition"
          stop
        end if
        if (str%r(i_a,i) < -0.5 * str%box_edges(i)) then
          write (*,*) "atoms move further than the length of the box within one step!!"
          write (*,*) "ERROR written from apply_boundary_condition"
          stop
        end if       
      end do

    end do
  end subroutine apply_boundary_condition
  
  
  function copy_structure(str_in) result(str_out)
    type (structure), intent(in) :: str_in
    type (structure) :: str_out
    
    integer :: i
    real(db) :: sum_out, sum_in
    
    str_out%box_edges = str_in%box_edges
    str_out%density = str_in%density
    str_out%title = str_in%title
    
    allocate(str_out%atoms(size(str_in%atoms)))  ! allocate mass, name...
    allocate(str_out%r(size(str_in%atoms),ndim)) ! allocate coordinates

    do i = 1, size(str_in%atoms)
      str_out%atoms(i)%mass = str_in%atoms(i)%mass
      str_out%atoms(i)%inv_mass = str_in%atoms(i)%inv_mass
      str_out%atoms(i)%element_type = str_in%atoms(i)%element_type
    end do
    
    str_out%r = str_in%r
  
  end function copy_structure

  subroutine swap_atoms(a_structure, atom_number1, atom_number2)
    type (structure), intent(inout) :: a_structure
    integer, intent(in) :: atom_number1, atom_number2

  end subroutine swap_atoms

end module structure_class
