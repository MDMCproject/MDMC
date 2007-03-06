module histogram_class
use various_constants_class
use structure_class

implicit none


  public :: make_histogram
  public :: copy_histogram
  public :: cal_histogram
  public :: accum_histogram, clear_histogram


  ! it is assumed that histogram stores distances in bins of size
  ! bin_length; with the first bin representing distances between
  ! zero and bin_length and the last bin representing distances
  ! between (n-bins-1)*bin_length and n_bins*bin_length
  
  type histogram
    real(db) :: bin_length   
    
    ! used when calculating a single histogram
    integer, dimension(:), allocatable :: val
    
    ! used when using accumulating histograms
    integer, dimension(:), allocatable :: sum 
    integer :: n_accum  
    
    ! Notice in the current design %val and %sum are not
    ! in sync. 
    
    ! ....
    
  end type histogram

contains


  ! Same as cal_histogram but in additional this function also adds the
  ! calculated histogram to the hist%sum attribute and increase hist%n_accum

  subroutine accum_histogram(hist, str)
    type (histogram), intent(inout) :: hist
    type (structure), intent(in) :: str    
    
    ! fill up the hist%val array
    
    call cal_histogram(hist, str)
    
    
    ! add the hist%val array to the hist%sum array
    
    hist%n_accum = hist%n_accum + 1
    hist%sum = hist%sum + hist%val
  
  end subroutine accum_histogram
  
  
  ! This one is used together with accum_histogram. 
  ! It resets the hist%n_accum and hist%sum attributes
  
  subroutine clear_histogram(hist)
    type (histogram), intent(inout) :: hist  
  
    hist%val = 0
    hist%n_accum = 0
    hist%sum = 0
  
  end subroutine clear_histogram
  
  
  function make_histogram(r_max, bin_length) result(hist)
    real(db), intent(in) :: r_max, bin_length
    type (histogram) :: hist
    
    
    hist%bin_length = bin_length 
    
    allocate(hist%val(floor(r_max/bin_length)))
    allocate(hist%sum(floor(r_max/bin_length)))
    
    ! initialise
    
    hist%val = 0
    hist%sum = 0
    hist%n_accum = 0
  
  end function make_histogram
  
  
  function copy_histogram(hist_in) result(hist_out)
    type (histogram), intent(in) :: hist_in
    type (histogram) :: hist_out
    
    integer :: n_bin 
    
    n_bin = size(hist_in%val)
    
    if (n_bin < 1) then
      print *, "ERROR in copy histogram"
      print *, "Try to copy empty histogram"
      stop
    end if
    
    hist_out%bin_length = hist_in%bin_length
    
    allocate(hist_out%val(n_bin))
    allocate(hist_out%sum(n_bin))
    
    hist_out%val = hist_in%val
    hist_out%sum = hist_in%sum
    hist_out%n_accum = hist_in%n_accum
  
  end function copy_histogram  
  
  
  ! Populate the histogram%val array (but doesn't alter %sum and %n_accum)
  
  subroutine cal_histogram(hist, str)
    type (histogram), intent(inout) :: hist
    type (structure), intent(in) :: str
    
    real(db) :: rr_max, r_max
    integer :: n_tot        ! number of atoms
    integer :: i, i1, i2
    integer :: which_bin    ! where which_bin=1 is the first bin: [0:bin_length]
    real (db) :: rr
    real (db), dimension(ndim) :: diff_vec   
		
		
	r_max = hist%bin_length * size(hist%val)
    rr_max = r_max * r_max
    
    n_tot = size(str%atoms)   ! number of atoms
      
    hist%val = 0

    ! Do different summations depending on whether the nearest neighbour list is
    ! in use    
    
    if (str%nn_list%ignore_list == .true. .or. str%nn_list%r_cut < r_max) then
    
      do i1 = 1, n_tot
        do i2 = i1+1, n_tot
          diff_vec = str%r(i1,:) - str%r(i2,:)
          
          call apply_boundary_condition_to_vector_expensive(diff_vec, str%box_edges)
          
          
          rr = sum(diff_vec*diff_vec)
          
          if (rr < rr_max) then
            which_bin = ceiling(sqrt(rr)/hist%bin_length) 
            hist%val(which_bin) = hist%val(which_bin) + 2 ! two to count both i1<i2 and i1>i2
          end if
  			
        end do
      end do    
		
    else
      ! for now also assumes that distances are stored as r2
      ! and I have put this if statement in so that I don't forget to
      ! change the code in this function when I no longer assume this
    
      if (str%nn_list%what_is_stored /= "r2") then
        print *, "ERROR: what_is_stored is different from r2"
        print *, "Time to update code in nn_update_histogram subroutine"
        stop
      end if
     
      do i = 1, str%nn_list%n_pairs
        rr = str%nn_list%dists(i) 
        if (rr < rr_max) then
          which_bin = ceiling(sqrt(rr)/hist%bin_length) 
          hist%val(which_bin) = hist%val(which_bin) + 2 ! two to count both i1<i2 and i1>i2
        end if
      end do

		end if

  end subroutine cal_histogram


end module histogram_class