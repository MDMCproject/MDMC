module structure_nn_methods_class
use various_constants_class
use structure_type_definition_class

implicit none

  public :: copy_near_neighb_list
  public :: make_near_neighb_list
  public :: build_near_neighb
  public :: update_nn_list_flags


  ! used for the cell method (use rapaport's notation)
  integer, parameter, private :: n_offset = 14
  
  private :: get_num_near_neighb 
  private :: get_num_near_neighb_with_cell  ! using cell method
  private :: get_num_near_neighb_without_cell  ! not using cell method
  private :: build_near_neighb_with_cell    
  private :: build_near_neighb_without_cell  
  private :: copy_cell_list
  private :: cal_nn_distances     ! populates the dists attribute of the near_neighb_list only  

contains

  function copy_cell_list(cell_in) result(cell_out)
    type (cell_list), intent(in) :: cell_in
    type (cell_list) :: cell_out
    
    integer :: i
    
    cell_out%ignore_cell_method = cell_in%ignore_cell_method
    cell_out%num_cells = cell_in%num_cells
    
    allocate(cell_out%list(size(cell_in%list))) 
    cell_out%list = cell_in%list
    
  end function copy_cell_list


  function copy_near_neighb_list(nn_list_in) result(nn_list_out)
    type (near_neighb_list), intent(in) :: nn_list_in
    type (near_neighb_list) :: nn_list_out
  
    integer :: i
    
    nn_list_out%ignore_list = nn_list_in%ignore_list 
    
    if ( nn_list_in%ignore_list == .true. ) then
      return
    end if 
  
    nn_list_out%pairs_allocated = nn_list_in%pairs_allocated
    nn_list_out%n_pairs = nn_list_in%n_pairs
    nn_list_out%needs_updating = nn_list_in%needs_updating
    nn_list_out%delta_r = nn_list_in%delta_r
    nn_list_out%r_cut = nn_list_in%r_cut
    nn_list_out%max_an_atom_has_moved = nn_list_in%max_an_atom_has_moved
    nn_list_out%what_is_stored = nn_list_in%what_is_stored
  
  
    allocate(nn_list_out%pairs(size(nn_list_in%pairs)))
    nn_list_out%pairs = nn_list_in%pairs

    
    allocate(nn_list_out%dists(size(nn_list_in%dists)))
    nn_list_out%dists = nn_list_in%dists
    
    nn_list_out%cells = copy_cell_list(nn_list_in%cells)
    
  end function copy_near_neighb_list

  
  subroutine update_nn_list_flags(max_an_atom_has_moved, nn_list)
    real(db), intent(in) :: max_an_atom_has_moved
    type (near_neighb_list), intent(inout) :: nn_list
    
    nn_list%max_an_atom_has_moved = nn_list%max_an_atom_has_moved &
                             + max_an_atom_has_moved  
                             
    if (nn_list%max_an_atom_has_moved > 0.5 * nn_list%delta_r) then
      nn_list%needs_updating = .true.
      nn_list%max_an_atom_has_moved = 0.0
    end if 
    
  end subroutine update_nn_list_flags
  

  subroutine make_near_neighb_list(str, r_cut, delta_r)
    type (structure), intent(inout) :: str
    real (db), intent(in) :: r_cut, delta_r
    
    integer :: i
    integer :: ignore_cell_method    ! used for determining whether 
                                     ! to ignore cell method
    integer :: n_tot
    
    ! first of all determine whether or not to completely ignore 
    ! using the nearest neighbour method. Perhaps the only reason 
    ! why you want to do this is if the memory storage required out-weights
    ! the other benifits - for now simply just assume that everything
    ! is fine and only if there is no use-near-neighbour-method element 
    ! specified in the input file (in which case this causes r_cut=0.0)
    ! will the nn_list be ignored completely.
    
    if (r_cut /= 0.0) then
      str%nn_list%ignore_list = .false.
    else
      return  
    end if
    
    n_tot = size(str%atoms)
    
    
    ! stuff to do with when to update neighbor list
    
    str%nn_list%needs_updating = .true.
    str%nn_list%delta_r = delta_r
    str%nn_list%r_cut = r_cut   
    str%nn_list%max_an_atom_has_moved = 0.0
    
                                              
    ! Determine whether or not to use cell method
                                             
    str%nn_list%cells%num_cells = floor(str%box_edges / (r_cut+delta_r))
    
    ignore_cell_method = .false.
    
    do i = 1, ndim
      if (str%nn_list%cells%num_cells(i) < 4) then
        ignore_cell_method = .true.
        exit
      end if
    end do
    
    str%nn_list%cells%ignore_cell_method = ignore_cell_method
    
    allocate(str%nn_list%cells%list(n_tot+product(str%nn_list%cells%num_cells)))
    
    
    ! get number of nn atom pairs required for the structure: str         
    
    str%nn_list%n_pairs = get_num_near_neighb(str%nn_list%cells, &
                         str, r_cut+delta_r)

    ! nn_list cannot possibly be bigger than N(N-1)/2, where N is the total
    ! number of atoms in the box. 
    ! Proper exception handling should have been used below but for now simply 
    ! use brute force
    
    if (str%nn_list%n_pairs > n_tot*(n_tot-1) / 2 ) then
      write(*,*) " "
      write(*,*) "ERROR in make_near_neighb_list"
      write(*,'(a,i6)') "n_pairs = ", str%nn_list%n_pairs
      write(*,'(a,i6)') "and number of possible pairs are ", n_tot*(n_tot-1) / 2
      stop
    end if
    

    str%nn_list%pairs_allocated = str%nn_list%n_pairs+10*n_tot
    
    if (str%nn_list%pairs_allocated > n_tot*(n_tot-1) / 2 ) then
      str%nn_list%pairs_allocated = n_tot*(n_tot-1) / 2
    end if
    
    allocate(str%nn_list%pairs( 2*str%nn_list%pairs_allocated ))
    allocate(str%nn_list%dists( str%nn_list%pairs_allocated ))
    
    ! general warning since no point in using nearest neighbour method
    ! if nn_list close to max number of allowed and this is likely to 
    ! be the case when
    
    if ( r_cut+delta_r > maxval(str%box_edges) ) then
      write(*,*) " "
      write(*,*) "WARNING: using nearest neighbour method may be inefficient"
      write(*,'(a,i6)') "pairs allocated = ", str%nn_list%n_pairs
      write(*,'(a,i6)') "and number of possible pairs are ", n_tot*(n_tot-1) / 2
      write(*,*) " "
    end if
    
    
    ! a nn-list is here a member of a structure type. It is implicitely assumed
    ! this nn-list is always in sync with the structure, therefore:
    
    call build_near_neighb(str)

  end subroutine make_near_neighb_list


  subroutine build_near_neighb(str)
    type (structure), intent(inout) :: str
    
    ! For now you are not allowed to call this function
    ! when nn_list%ignore_list == .true.
    ! This restriction will probably be removed later
    
    if (str%nn_list%ignore_list == .true.) then
      write(*,'(a)') "ERROR in build_near_neighb"
      write(*,'(a)') "Not allowed to rebuild nn_list when"
      write(*,'(a)') "ignore_list == .true."
      stop
    end if
    
    ! if needs_updating = .false. then the nn-list is assumed to still
    ! be ok - hence no need building a new one. However, the atoms may
    ! have moved and the nn_list%dists still needs to be updated to these
    ! new atomic positions
    !
    ! Perhaps at some point will somehow make build_near_neighb or another
    ! function a 'permament' member of structure such that whenever the 
    ! atomic position coordinates are changed the nn_list is updated
    ! accordingly
    
    if (str%nn_list%needs_updating == .false.) then
      call cal_nn_distances(str)
      return
    end if
    
    
    if (str%nn_list%cells%ignore_cell_method == .true.) then
      call build_near_neighb_without_cell(str)
    else
      call build_near_neighb_with_cell(str)
    end if
    
    
    ! at this point a fresh nn-list has just been build, hence
    
    str%nn_list%needs_updating = .false.
    
  end subroutine build_near_neighb


!!!!!!!!!!!!!!!!!!!!!!! private functions/subroutines !!!!!!!!!!!!!!!!!!!!!!! 
  
  ! build new nn-list and populate dists array with, for now, r2 values
  
  subroutine build_near_neighb_without_cell(str)
    type (structure), intent(inout) :: str
    
    integer :: num  ! number of nearest neighbour pairs     
    
    real(db), dimension(ndim) :: diff_vec  ! vector between two atoms    
    integer :: i, i1, i2, n_tot
    real(db) :: r_cut_neighbour, rr_cut_neighbour    
    real(db) :: rr  ! used when storing |r|^2    
    
    n_tot = size(str%atoms)
   
    r_cut_neighbour = str%nn_list%r_cut + str%nn_list%delta_r
    rr_cut_neighbour = r_cut_neighbour*r_cut_neighbour    
    
    num = 0
    do i1 = 1, n_tot
      do i2 = i1+1, n_tot
        diff_vec = str%r(i1,:) - str%r(i2,:)
      
        call apply_boundary_condition_one_atom(diff_vec, str%box_edges)
      
        rr = sum(diff_vec*diff_vec)
        if (rr < rr_cut_neighbour) then
          num = num + 1
          
          if (num > str%nn_list%pairs_allocated) then
            write (*, '(a)') "ERROR in build_near_neighb_with_cell"
            stop
          end if
                    
          str%nn_list%pairs(2*num-1) = i1
          str%nn_list%pairs(2*num) = i2          
          str%nn_list%dists(num) = rr
        end if
  		
      end do
    end do
    
    ! update number of pairs
    str%nn_list%n_pairs = num
    
  end subroutine build_near_neighb_without_cell
  
  
  ! same as build_near_neighb_without_cell but using cell method
  
  subroutine build_near_neighb_with_cell(str)
    type (structure), intent(inout) :: str
    
    integer :: num  ! number of nearest neighbour pairs     
    real(db) :: r_cut_neighbour, rr_cut_neighbour
    real(db), dimension(ndim) :: inverse_cell_width, dummy
    real(db), dimension(ndim) :: shift  ! atom shift caused by cell method
    real(db), dimension(ndim) :: diff_vec  ! vector between two atoms
    integer, dimension(ndim) :: cp_3d, cp1_3d, cp2_3d  ! 3D cell positions
    integer :: cp_1d, cp1_1d, cp2_1d  ! 1D cell positions
    integer :: n_tot, i, j, cp1x, cp1y, cp1z, j1, j2
    real(db) :: rr  ! used when storing |r|^2     
    
    integer, dimension(ndim,14) :: offset_vals
  
    offset_vals = reshape( (/0,0,0, 1,0,0, 1,1,0, 0,1,0, -1,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1, -1,1,1, -1,0,1, -1,-1,1, 0,-1,1, 1,-1,1/), & 
                  shape(offset_vals))
    
    n_tot = size(str%atoms)
    r_cut_neighbour = str%nn_list%r_cut + str%nn_list%delta_r
    rr_cut_neighbour = r_cut_neighbour*r_cut_neighbour
    
    inverse_cell_width = str%nn_list%cells%num_cells / str%box_edges
        
    do i = 1, product(str%nn_list%cells%num_cells)
      str%nn_list%cells%list(n_tot+i) = -1
    end do
    
    do i = 1, n_tot
      dummy = str%r(i,:) + 0.5 * str%box_edges
      
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
      cp_1d = (cp_3d(3) * str%nn_list%cells%num_cells(2) + cp_3d(2)) * & 
              str%nn_list%cells%num_cells(1) + cp_3d(1)
      
      
      ! Here I could check whether I have remembered to do wrap around 
      
      !if (cp_1d < 0) then
      !  write(*,*) "ERROR - serious error in simulation in build_near_neighb_with_cell"
      !  write(*,*) "Likely causes are that you have forgot to do wrap around"
      !  write(*,'(a,i6)') "Value of cp_1d = ", cp_1d
      !  stop
      !end if 
      
      
      ! in the speciel cells the 'first' element in each cell is stored in 
      ! the last n_tot+1, n_tot+2, ... n_tot+product(num_cells) elements
      cp_1d = cp_1d + n_tot + 1
      
      ! put the previous atom listed in cell-list cells%list(cp_1d) into
      ! cells%list(i) (which was also the 1st atom listed in that list)
      str%nn_list%cells%list(i) = str%nn_list%cells%list(cp_1d)
      
      ! now let this cell-list point to the current atom number (here i). This 
      ! then takes on the significance of the first atom in this cell-list; 
      ! the content of that element then points to the next atom and so on 
      ! until -1 is encountered
      str%nn_list%cells%list(cp_1d) = i
      

    end do
    
    !do i = 1, n_tot+product(str%nn_list%cells%num_cells)
    !  write(*,'(a,i10, a,i10)') "n=",i, " cells(n)=", str%nn_list%cells%list(i)
    !end do
    !stop
    
  
    num = 0
    do cp1z = 0, str%nn_list%cells%num_cells(3)-1
      do cp1y = 0, str%nn_list%cells%num_cells(2)-1
        do cp1x = 0, str%nn_list%cells%num_cells(1)-1
          cp1_3d(1)=cp1x; cp1_3d(2)=cp1y; cp1_3d(3)=cp1z
    !      
          cp1_1d = (cp1_3d(3) * str%nn_list%cells%num_cells(2) + cp1_3d(2)) * & 
                   str%nn_list%cells%num_cells(1) + cp1_3d(1) + n_tot + 1
          
          do i = 1, n_offset
            cp2_3d = cp1_3d + offset_vals(:,i)
            shift = 0.0
            
            ! wrap around
            do j = 1, ndim
              if (cp2_3d(j) >= str%nn_list%cells%num_cells(j)) then
                cp2_3d(j) = 0
                shift(j) = str%box_edges(j)
              else if (cp2_3d(j) < 0) then
                cp2_3d(j) = str%nn_list%cells%num_cells(j) - 1
                shift(j) = - str%box_edges(j)
              end if
            end do
            
            cp2_1d = (cp2_3d(3) * str%nn_list%cells%num_cells(2) + cp2_3d(2)) * & 
                     str%nn_list%cells%num_cells(1) + cp2_3d(1) + n_tot + 1
            
            j1 = str%nn_list%cells%list(cp1_1d)
            
            
            do while (j1 >= 0)
              j2 = str%nn_list%cells%list(cp2_1d)
              do while (j2 >= 0)
                !write(*,'(a,i6,a,i6)') "j1=", j1-1, " j2=", j2-1
                ! remember to not counting twice distances within same cell
                ! i.e. when cp1_1d==cp2_1d
                if (cp1_1d /= cp2_1d .or. j2 < j1) then
                  diff_vec = str%r(j1,:) - str%r(j2,:)
                  diff_vec = diff_vec - shift
                  !write(*,'(a,i6,a,i6,a,i6)') "j1=", j1-1, " j2=", j2-1, " num=", num
                  !write(*,'(3f12.6)') str%r(j1,:)
                  !write(*,'(3f12.6)') str%r(j2,:)
                  !write(*,'(3f12.6)') diff_vec
                  !write(*,'(3f12.6)') shift
                  rr = sum(diff_vec*diff_vec)
                  if (rr < rr_cut_neighbour) then
                    
                    ! this number must never be bigger than
                    ! str%nn_list%pairs_allocated
                    num = num + 1
                    
                    if (num > str%nn_list%pairs_allocated) then
                      write (*, '(a)') "ERROR in build_near_neighb_with_cell"
                      stop
                    end if
                    
                    str%nn_list%pairs(2*num-1) = j1
                    str%nn_list%pairs(2*num) = j2
                    str%nn_list%dists(num) = rr
                    !write(*,'(a,i6,a,i6,a,i6)') "j1=", j1-1, " j2=", j2-1, " num=", num
                    !stop
                  end if
                end if
              
                j2 = str%nn_list%cells%list(j2)
              end do
              j1 = str%nn_list%cells%list(j1)
            end do
          
            
          end do
        end do
      end do
    end do
    
    
    ! update number of pairs
    str%nn_list%n_pairs = num
    
  end subroutine build_near_neighb_with_cell


  ! populate the dists array only, for now, with r2 values. This function
  ! does not 'build' a new nn-list but assumes the existing nn-list is still
  ! ok
  
  subroutine cal_nn_distances(str)
    type (structure), intent(inout) :: str
    
    integer :: i1, i2, j, i
    real(db), dimension(3) :: diff_vec
    
    do j = 1, str%nn_list%n_pairs
  	  i1 = str%nn_list%pairs(2*j-1)
  	  i2 = str%nn_list%pairs(2*j)
      
      diff_vec = str%r(i1,:) - str%r(i2,:)
        
      call apply_boundary_condition_one_atom(diff_vec, str%box_edges)
        
      str%nn_list%dists(j) = sum(diff_vec*diff_vec)
    end do
    
      str%nn_list%what_is_stored = "r2"
      
  end subroutine cal_nn_distances
  
  
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
      diff_vec = str%r(i1,:) - str%r(i2,:)
      
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
    
    !write(*, '(a,3f10.3)') "inverse_cell_witdh ", inverse_cell_width
    
    do i = 1, product(cells%num_cells)
      cells%list(n_tot+i) = -1
    end do
    
    do i = 1, n_tot
      dummy = str%r(i,:) + 0.5 * str%box_edges
      
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
                  diff_vec = str%r(j1,:) - str%r(j2,:)
                  diff_vec = diff_vec - shift
                  
                  if (sum(diff_vec*diff_vec) < rr_cut_neighbour) then
                    num = num + 1
                  
                  
                   ! write(*, '(3i6)') j1, j2, num
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


end module structure_nn_methods_class