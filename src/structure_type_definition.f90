module structure_type_definition_class
use various_constants_class

implicit none

  ! moved here because used in structure_nn_methods 
  public :: apply_boundary_condition_one_atom


  ! the below is the fastest and most compressed collection of single link
  ! list I have ever seen. Its drawback is that it is a pain to handle but
  ! luckily it only needs to be used in some private subroutines to do with
  ! calculating attributes in near_neighb_list
  
  type cell_list
    logical :: ignore_cell_method = .true.
    integer, dimension(ndim) :: num_cells
    integer, dimension(:), allocatable :: list
  end type cell_list


  ! This type is only include to speed up the calculation. Most importantly by
  ! replacing the double summation over atoms by a sum over nearest neighbours
  ! and secondly by storing distances that are used as input to values to the
  ! various potential and FOM functions. Therefore make sure that this stored
  ! values are always up-to-date; this is the after a near-neighb-list
  ! has just been build and after each call to cal_nn_distances()

  type near_neighb_list
    ! for the case where it is either decide
    ! that this list is no-point or requires more
    ! memory than what can be allocated or otherwise
    ! what don't what to use near neighbour method
    logical :: ignore_list = .true.  
  
    ! the list - for now just pairs of j1, j2
    ! Notice the reason for the n_pairs attribute is that there may be more
    ! space allocated then there are nearest neigbour pairs for a given structure
    ! For this reason also when looping over pairs you most not use pairs_allocated
    ! or size(pairs) but always n_pairs
    integer, dimension(:), allocatable :: pairs
    integer :: pairs_allocated
    integer :: n_pairs
  
    ! stuff decide if list needs updating and with rebuilding
    logical :: needs_updating = .true.
    real(db) :: delta_r, r_cut
    real(db) :: max_an_atom_has_moved = 0.0_db
  
    ! storage of numbers for potential later use
    character(len=2) :: what_is_stored
    real(db), dimension(:), allocatable :: dists
    
    type (cell_list) :: cells
  end type near_neighb_list
  
  
  type atom
    character(len=2)  :: element_type = "Ar"
    real(db) :: mass = 39.95        ! in units of amu
    real(db) :: inv_mass = 0.025031289111  ! to save comp time
  end type atom


  ! Associate a near_neighb_list attribute the the type structure.
  ! If the near_neighb_list attribute ignore_list = .false. then
  ! the near_neighb_list is assumed to match the resulting 
  ! structure. This attribute is put in order to try and 'hide' it.
  ! Meaning that the nearest neighbour method may be used behind 
  ! the scene but would not be explicitely accessed from any of the
  ! run-control-f90 code-scripts.

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
      
    character(len=120) :: title = " "
    
    type (near_neighb_list) :: nn_list
  end type structure
  
contains

  subroutine apply_boundary_condition_one_atom(vec, box_edges)
    real (db), dimension(ndim), intent(inout) :: vec
    real (db), dimension(ndim), intent(in) :: box_edges
    
    integer :: i
    
    do i = 1, ndim
      if (vec(i) >= 0.5 * box_edges(i)) then
        vec(i) = vec(i) - box_edges(i)
      end if
      if (vec(i) < -0.5 * box_edges(i)) then
        vec(i) = vec(i) + box_edges(i)
      end if       
    end do
    
  end subroutine apply_boundary_condition_one_atom  
  
end module structure_type_definition_class