module histogram_class
use various_constants_class
use structure_class
!use near_neighb_class

implicit none


  ! it is assumed that histogram stores distances in bins of size
  ! bin_length; with the first bin representing distances between
  ! zero and bin_length and the last bin representing distances
  ! between (n-bins-1)*bin_length and n_bins*bin_length
  
  type histogram
    real(db) :: bin_length    ! length of a bin
    integer, dimension(:), allocatable :: h    ! the histogram 
    
    real(db) :: r_max  ! included for convenience
  end type histogram

contains

  function make_histogram(r_max, n_bin) result(hist)
    real(db), intent(in) :: r_max
    integer, intent(in) :: n_bin
    type (histogram) :: hist
    
    hist%bin_length = r_max / n_bin 
    hist%r_max = r_max
    
    allocate(hist%h(n_bin))
  
  end function make_histogram
  
  function copy_histogram(hist_in) result(hist_out)
    type (histogram), intent(in) :: hist_in
    type (histogram) :: hist_out
    
    integer :: n_bin
    
    n_bin = size(hist_in%h)
    
    if (n_bin < 1) then
      print *, "ERROR in copy histogram"
      print *, "Try to copy empty histogram"
      stop
    end if
    
    hist_out%r_max = hist_in%r_max
    hist_out%bin_length = hist_in%bin_length
    
    allocate(hist_out%h(n_bin))
    hist_out%h = hist_in%h
  
  end function copy_histogram  
  
  
!  subroutine cal_histogram_nn(hist, nn)
!    type (histogram), intent(inout) :: hist
!    type (near_neighb_list), intent(in) :: nn
    
!    if (nn%what_is_stored == "r2") then
       
    !
!    else if (nn%what_is_stored == "r") then
    !
    
!    else
    !
    
!    end if
    
!  end subroutine cal_histogram_nn  
  
  
  subroutine cal_histogram(hist, str)
    type (histogram), intent(inout) :: hist
    type (structure), intent(in) :: str
    
    real(db) :: rr_max
    integer :: n_tot        ! number of atoms
    integer :: i, i1, i2
    integer :: which_bin    ! where which_bin=1 is the first bin: [0:bin_length]
    real (db) :: rr
		real (db), dimension(ndim) :: diff_vec   
		
		
    rr_max = hist%r_max * hist%r_max
    
		n_tot = size(str%atoms)   ! number of atoms
    
    
    hist%h = 0
    
		do i1 = 1, n_tot
		  do i2 = i1+1, n_tot
        diff_vec = str%r(i1,:) - str%r(i2,:)
        
        do i = 1, ndim
          if (diff_vec(i) >= 0.5_db * str%box_edges(i)) then
            diff_vec(i) = diff_vec(i) - str%box_edges(i)
          end if
          if (diff_vec(i) < -0.5_db * str%box_edges(i)) then
            diff_vec(i) = diff_vec(i) + str%box_edges(i)
          end if       
        end do
        
        rr = sum(diff_vec*diff_vec)
        
        if (rr < rr_max) then
          which_bin = ceiling(sqrt(rr)/hist%bin_length) 
          hist%h(which_bin) = hist%h(which_bin) + 1
        end if
			
			end do
		end do    

  end subroutine cal_histogram


end module histogram_class