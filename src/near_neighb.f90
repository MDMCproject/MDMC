module near_neighb_class
use structure_class

  implicit none

  ! used for the cell method (use rapaport's notation)
  integer, parameter, private :: n_offset = 14

  
  private :: get_num_near_neighb 
  private :: get_num_near_neighb_with_cell  ! using cell method
  private :: get_num_near_neighb_without_cell  ! not using cell method
  
  ! the below is the fastest and most compressed collection of single link
  ! list I have ever seen. Its drawback is that it is a pain to handle but
  ! luckily it only needs to be used in this module
  type cell_list
    private
    logical :: ignore_cell_method = .true.
    integer, dimension(ndim) :: num_cells
    integer, dimension(:), allocatable :: list
  end type cell_list

  ! This type is only include to speed up the calculation. Most importantly by
  ! replacing the double summation over atoms by a sum over nearest neighboors
  ! and secondly by storing distances that may be used at a later stage.
  ! Notice that this distances are here stored as REAL rather then REAL (DB)
  ! to save memory, since this type will be the single most memory hungry 
  ! component of this code and knowing these distances to 6 digits accuracy
  ! is e.g. definitely good enough for calculating a histogram, pdf etc from

  type near_neighb_list
    ! the list - for now just pairs of j1, j2
    integer, dimension(:), allocatable :: pairs
    integer :: pairs_allocated
    integer :: n_pairs
  
    ! stuff decide if list needs updating and with rebuilding
    logical :: needs_updating = .true.
    real(db) :: delta_r, r_cut
    real(db) :: max_an_atom_has_moved
  
    ! for the case where it is either decide
    ! that this list is no-point or requires more
    ! memory than what can be allocated
    logical :: ignore_list = .true.
  
    ! storage of numbers for potential later use
    character(len=2) :: what_is_stored
    real, dimension(:), allocatable :: dists
    
    type (cell_list) :: cells
  end type near_neighb_list



contains

  function get_num_near_neighb(cells, str, r_cut_neighbour) &
         result (num)
    type (cell_list), intent(inout) :: cells
    type (structure), intent(in) :: str
    real(db), intent(in) :: r_cut_neighbour
    integer :: num

    if (cells%ignore_cell_method) then
      num = get_num_near_neighb_without_cell(str, r_cut_neighbour)
    else
      num = get_num_near_neighb_with_cell(cells, str, r_cut_neighbour)
    end if
  
  end function get_num_near_neighb
  
  
  function get_num_near_neighb_without_cell(str, r_cut_neighbour) &
         result (num)
    type (structure), intent(in) :: str
    real(db), intent(in) :: r_cut_neighbour
    integer :: num

    real(db), dimension(ndim) :: diff_vec  ! vector between two atoms    
    integer :: i, i1, i2, n_tot
    real(db) :: rr_cut_neighbour    
    
    n_tot = size(str%atoms)
   
    rr_cut_neighbour = r_cut_neighbour*r_cut_neighbour    
    
    num = 0
		do i1 = 1, n_tot
		  do i2 = i1+1, n_tot
        diff_vec = str%atoms(i1)%r - str%atoms(i2)%r
        
        do i = 1, ndim
          if (diff_vec(i) >= 0.5 * str%box_edges(i)) then
            diff_vec(i) = diff_vec(i) - str%box_edges(i)
          end if
          if (diff_vec(i) < -0.5 * str%box_edges(i)) then
            diff_vec(i) = diff_vec(i) + str%box_edges(i)
          end if       
        end do
        
        if (sum(diff_vec*diff_vec) < rr_cut_neighbour) then
          num = num + 1
        end if
			
			end do
		end do
  
  end function get_num_near_neighb_without_cell  

  
  function get_num_near_neighb_with_cell(cells, str, r_cut_neighbour) &
           result (num)
    type (cell_list), intent(inout) :: cells
    type (structure), intent(in) :: str
    real(db), intent(in) :: r_cut_neighbour
    integer :: num
    
  
    real(db) :: rr_cut_neighbour
    real(db), dimension(ndim) :: inverse_cell_width, dummy
    real(db), dimension(ndim) :: shift  ! atom shift caused by cell method
    real(db), dimension(ndim) :: diff_vec  ! vector between two atoms
    integer, dimension(ndim) :: cp_3d, cp1_3d, cp2_3d  ! 3D cell positions
    integer :: cp_1d, cp1_1d, cp2_1d  ! 1D cell positions
    integer :: n_tot, i, j, cp1x, cp1y, cp1z, j1, j2
    
    integer, dimension(ndim,14) :: offset_vals
  
    offset_vals = reshape( (/0,0,0, 1,0,0, 1,1,0, 0,1,0, -1,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1, -1,1,1, -1,0,1, -1,-1,1, 0,-1,1, 1,-1,1/), & 
                  shape(offset_vals))
    
    
    
    n_tot = size(str%atoms)
   
    rr_cut_neighbour = r_cut_neighbour*r_cut_neighbour
    
    inverse_cell_width = cells%num_cells / str%box_edges
    
    write(*, '(a,3f10.3)') "inverse_cell_witdh ", inverse_cell_width
    
    do i = 1, product(cells%num_cells)
      cells%list(n_tot+i) = -1
    end do
    
    do i = 1, n_tot
      dummy = str%atoms(i)%r + 0.5 * str%box_edges
      
      ! find 3D cell position
      cp_3d = floor(dummy * inverse_cell_width)
      
      ! here assume that we are in 3D
      ! locating the cell position of an atom in a 1D array by having the fastest
      ! running dimension to be (1), second fastest (2) and slowest (3). 
      ! that means (cp_3d(1),cp_3d(2), cp_3d(3)) has the following index in a 
      ! 1D (which first entry is index=0):
      !
      ! index = cp_3d(3) * num_cells(1)*num_cells(2) 
      !         + cp_3d(2) * num_cells(1)
      !         + cp_3d(1)
      ! which is what is calculated below
      ! cp_1d stands for cell position 1D
      cp_1d = (cp_3d(3) * cells%num_cells(2) + cp_3d(2)) * cells%num_cells(1) &
              + cp_3d(1)
      
      ! in the speciel cells the 'first' element in each cell is stored in 
      ! the last n_tot+1, n_tot+2, ... n_tot+product(num_cells) elements
      cp_1d = cp_1d + n_tot + 1
      
      ! put the previous atom listed in cell-list cells%list(cp_1d) into
      ! cells%list(i) (which was also the 1st atom listed in that list)
      cells%list(i) = cells%list(cp_1d)
      
      ! now let this cell-list point to the current atom number (here i). This 
      ! then takes on the significance of the first atom in this cell-list; 
      ! the content of that element then points to the next atom and so on 
      ! until -1 is encountered
      cells%list(cp_1d) = i
      
      !write(*, '(a,i10)') "i = ", i
      !write(*, '(a,3f10.3)') "dummy ", dummy
      !write(*, '(a,3i10)') "cp_3d ", cp_3d
      !write(*, '(a,i10)') "cp_1d ",cp_1d

    end do
    
    !do i = 1, n_tot+product(cells%num_cells)
    !  write(*,'(i10, i10)') i, cells%list(i)
    !end do
    
  
    num = 0
    do cp1z = 0, cells%num_cells(3)-1
      do cp1y = 0, cells%num_cells(2)-1
        do cp1x = 0, cells%num_cells(1)-1
          cp1_3d(1)=cp1x; cp1_3d(2)=cp1y; cp1_3d(3)=cp1z
    !      
          cp1_1d = (cp1_3d(3) * cells%num_cells(2) + cp1_3d(2)) * & 
                   cells%num_cells(1) + cp1_3d(1) + n_tot + 1
          
          do i = 1, n_offset
            cp2_3d = cp1_3d + offset_vals(:,i)
            shift = 0.0
            
            ! wrap around
            do j = 1, ndim
              if (cp2_3d(j) >= cells%num_cells(j)) then
                cp2_3d(j) = 0
                shift(j) = str%box_edges(j)
              else if (cp2_3d(j) < 0) then
                cp2_3d(j) = cells%num_cells(j) - 1
                shift(j) = - str%box_edges(j)
              end if
            end do
            
            cp2_1d = (cp2_3d(3) * cells%num_cells(2) + cp2_3d(2)) * & 
                     cells%num_cells(1) + cp2_3d(1) + n_tot + 1
            
            j1 = cells%list(cp1_1d)
            
          
            do while (j1 >= 0)
              j2 = cells%list(cp2_1d)
              do while (j2 >= 0)
              
                ! remember to not counting twice distances within same cell
                ! i.e. when cp1_1d==cp2_1d
                if (cp1_1d /= cp2_1d .or. j2 < j1) then
                  diff_vec = str%atoms(j1)%r - str%atoms(j2)%r
                  diff_vec = diff_vec - shift
                  
                  if (sum(diff_vec*diff_vec) < rr_cut_neighbour) then
                    num = num + 1
                  
                  end if
                end if
              
                j2 = cells%list(j2)
              end do
              j1 = cells%list(j1)
            end do
          
            
          end do
        end do
      end do
    end do

  
  end function get_num_near_neighb_with_cell

  function make_near_neighb_list(str, r_cut, delta_r) result (nn_list)
    type (structure), intent(in) :: str
    real (db), intent(in) :: r_cut, delta_r
    type (near_neighb_list) :: nn_list
    
    integer :: i
    integer :: ignore_cell_method = .false.  ! used for determining whether 
                                             ! to ignore cell method
    integer :: n_tot
    
    n_tot = size(str%atoms)
    
    ! stuff to do with when to update neighbor list
    
    nn_list%needs_updating = .true.
    nn_list%delta_r = delta_r
    nn_list%r_cut = r_cut   
    nn_list%max_an_atom_has_moved = 0.0
    
                                              
    ! do stuff associated with the cell method
                                             
    nn_list%cells%num_cells = floor(str%box_edges / (r_cut+delta_r))
    
    do i = 1, ndim
      if (nn_list%cells%num_cells(i) < 4) then
        ignore_cell_method = .true.
      end if
    end do
    
    nn_list%cells%ignore_cell_method = ignore_cell_method
    
    allocate(nn_list%cells%list(n_tot+product(nn_list%cells%num_cells)))
    
    write(*,'(a,3f10.3)') "box_edges ", str%box_edges
    write(*,'(a,3i4)') "num_cells ", nn_list%cells%num_cells
    
    
    ! the list        
    nn_list%n_pairs = get_num_near_neighb(nn_list%cells, &
                         str, r_cut+delta_r)

    nn_list%pairs_allocated = nn_list%n_pairs+n_tot
    allocate(nn_list%pairs( 2*nn_list%pairs_allocated ))
    allocate(nn_list%dists( nn_list%pairs_allocated ))
    
    write(*,'(a,i7)') "num_near_neighbours ", nn_list%n_pairs
    

  end function make_near_neighb_list
  
  
   subroutine build_near_neighb_with_cell(cells, str, nn_list, r_cut_neighbour)
     type (cell_list), intent(inout) :: cells
     type (structure), intent(in) :: str
     type (near_neighb_list), intent(out) :: nn_list
     real(db), intent(in) :: r_cut_neighbour
    

    integer :: num  ! number of nearest neighbour pairs
      
    real(db) :: rr_cut_neighbour
    real(db), dimension(ndim) :: inverse_cell_width, dummy
    real(db), dimension(ndim) :: shift  ! atom shift caused by cell method
    real(db), dimension(ndim) :: diff_vec  ! vector between two atoms
    integer, dimension(ndim) :: cp_3d, cp1_3d, cp2_3d  ! 3D cell positions
    integer :: cp_1d, cp1_1d, cp2_1d  ! 1D cell positions
    integer :: n_tot, i, j, cp1x, cp1y, cp1z, j1, j2
    
    integer, dimension(ndim,14) :: offset_vals
  
    offset_vals = reshape( (/0,0,0, 1,0,0, 1,1,0, 0,1,0, -1,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1, -1,1,1, -1,0,1, -1,-1,1, 0,-1,1, 1,-1,1/), & 
                  shape(offset_vals))
    
    
    
    n_tot = size(str%atoms)
   
    rr_cut_neighbour = r_cut_neighbour*r_cut_neighbour
    
    inverse_cell_width = cells%num_cells / str%box_edges
    
    write(*, '(a,3f10.3)') "inverse_cell_witdh ", inverse_cell_width
    
    do i = 1, product(cells%num_cells)
      cells%list(n_tot+i) = -1
    end do
    
    do i = 1, n_tot
      dummy = str%atoms(i)%r + 0.5 * str%box_edges
      
      ! find 3D cell position
      cp_3d = floor(dummy * inverse_cell_width)
      
      ! here assume that we are in 3D
      ! locating the cell position of an atom in a 1D array by having the fastest
      ! running dimension to be (1), second fastest (2) and slowest (3). 
      ! that means (cp_3d(1),cp_3d(2), cp_3d(3)) has the following index in a 
      ! 1D (which first entry is index=0):
      !
      ! index = cp_3d(3) * num_cells(1)*num_cells(2) 
      !         + cp_3d(2) * num_cells(1)
      !         + cp_3d(1)
      ! which is what is calculated below
      ! cp_1d stands for cell position 1D
      cp_1d = (cp_3d(3) * cells%num_cells(2) + cp_3d(2)) * cells%num_cells(1) &
              + cp_3d(1)
      
      ! in the speciel cells the 'first' element in each cell is stored in 
      ! the last n_tot+1, n_tot+2, ... n_tot+product(num_cells) elements
      cp_1d = cp_1d + n_tot + 1
      
      ! put the previous atom listed in cell-list cells%list(cp_1d) into
      ! cells%list(i) (which was also the 1st atom listed in that list)
      cells%list(i) = cells%list(cp_1d)
      
      ! now let this cell-list point to the current atom number (here i). This 
      ! then takes on the significance of the first atom in this cell-list; 
      ! the content of that element then points to the next atom and so on 
      ! until -1 is encountered
      cells%list(cp_1d) = i
      
      !write(*, '(a,i10)') "i = ", i
      !write(*, '(a,3f10.3)') "dummy ", dummy
      !write(*, '(a,3i10)') "cp_3d ", cp_3d
      !write(*, '(a,i10)') "cp_1d ",cp_1d

    end do
    
    !do i = 1, n_tot+product(cells%num_cells)
    !  write(*,'(i10, i10)') i, cells%list(i)
    !end do
    
  
    num = 0
    do cp1z = 0, cells%num_cells(3)-1
      do cp1y = 0, cells%num_cells(2)-1
        do cp1x = 0, cells%num_cells(1)-1
          cp1_3d(1)=cp1x; cp1_3d(2)=cp1y; cp1_3d(3)=cp1z
    !      
          cp1_1d = (cp1_3d(3) * cells%num_cells(2) + cp1_3d(2)) * & 
                   cells%num_cells(1) + cp1_3d(1) + n_tot + 1
          
          do i = 1, n_offset
            cp2_3d = cp1_3d + offset_vals(:,i)
            shift = 0.0
            
            ! wrap around
            do j = 1, ndim
              if (cp2_3d(j) >= cells%num_cells(j)) then
                cp2_3d(j) = 0
                shift(j) = str%box_edges(j)
              else if (cp2_3d(j) < 0) then
                cp2_3d(j) = cells%num_cells(j) - 1
                shift(j) = - str%box_edges(j)
              end if
            end do
            
            cp2_1d = (cp2_3d(3) * cells%num_cells(2) + cp2_3d(2)) * & 
                     cells%num_cells(1) + cp2_3d(1) + n_tot + 1
            
            j1 = cells%list(cp1_1d)
            
          
            do while (j1 >= 0)
              j2 = cells%list(cp2_1d)
              do while (j2 >= 0)
              
                ! remember to not counting twice distances within same cell
                ! i.e. when cp1_1d==cp2_1d
                if (cp1_1d /= cp2_1d .or. j2 < j1) then
                  diff_vec = str%atoms(j1)%r - str%atoms(j2)%r
                  diff_vec = diff_vec - shift
                  
                  if (sum(diff_vec*diff_vec) < rr_cut_neighbour) then
                    
                    ! this number must never be bigger than
                    ! nn_list%pairs_allocated
                    num = num + 1
                    
                    if (num > nn_list%pairs_allocated) then
                      write (*, '(a)') "ERROR in build_near_neighb_with_cell"
                      stop
                    end if
                    
                    nn_list%pairs(2*num-1) = j1
                    nn_list%pairs(2*num) = j2
                  
                  end if
                end if
              
                j2 = cells%list(j2)
              end do
              j1 = cells%list(j1)
            end do
          
            
          end do
        end do
      end do
    end do
    
 end subroutine build_near_neighb_with_cell


end module near_neighb_class
