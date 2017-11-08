module histogram_cutdown_class
use various_constants_class
use structure_class

implicit none

  public :: make_histogram_cutdown
  public :: cal_g_s_histogram, cal_g_d_histogram 
  public :: cal_g_s_histogram_no_wrap_around 

  type histogram_cutdown
    real(db) :: bin_length   
    integer, dimension(:), allocatable :: val
  end type histogram_cutdown  

contains

  function make_histogram_cutdown(r_max, bin_length) result(hist)
    real(db), intent(in) :: r_max, bin_length
    type (histogram_cutdown) :: hist
    
    hist%bin_length = bin_length 
    allocate(hist%val(floor(r_max/bin_length)))
  end function make_histogram_cutdown

  
  ! return sum of squared differences
  ! Applies periodic boundary wrap-around for act_r
  !
  function cal_g_s_histogram(g_s_hist, org_r, act_r, box_edges) result(rr_sum)
    type (histogram_cutdown), intent(inout) :: g_s_hist
    real(db), dimension(:,:), intent(in) :: org_r, act_r
    real (db), dimension(ndim), intent(in) :: box_edges    
    real(db) :: rr_sum
  
    real(db) :: rr_max, r_max
    integer :: n_tot        ! number of atoms
    integer :: i, i1
    integer :: which_bin    ! where which_bin=1 is the first bin: [0:bin_length]
    real (db) :: rr
    real (db), dimension(ndim) :: diff_vec   
		
    r_max = g_s_hist%bin_length * size(g_s_hist%val)
    rr_max = r_max * r_max
    
    n_tot = size(org_r(:,1))   ! number of atoms
      
    g_s_hist%val = 0  
    rr_sum = 0.0
  
    do i1 = 1, n_tot
      diff_vec = act_r(i1,:) - org_r(i1,:)
 
      ! The command below does not actually seem to make that much difference
      ! to the execution speed 
          
      call apply_boundary_condition_to_vector_expensive(diff_vec, box_edges)
                  
      rr = sum(diff_vec*diff_vec)
      
      if (rr < rr_max) then
        which_bin = ceiling(sqrt(rr)/g_s_hist%bin_length) 
        g_s_hist%val(which_bin) = g_s_hist%val(which_bin) + 1
      end if
      
      rr_sum = rr_sum + rr  
    end do      
  
  end function cal_g_s_histogram


  ! same as cal_g_s_histogram but no wrap-around
  ! used for calculating g_s that do not suffer from periodic boundary effect...
  !
  function cal_g_s_histogram_no_wrap_around(g_s_hist, org_r, act_r, box_edges) result(rr_sum)
    type (histogram_cutdown), intent(inout) :: g_s_hist
    real(db), dimension(:,:), intent(in) :: org_r, act_r
    real (db), dimension(ndim), intent(in) :: box_edges    
    real(db) :: rr_sum
  
    real(db) :: rr_max, r_max
    integer :: n_tot        ! number of atoms
    integer :: i, i1
    integer :: which_bin    ! where which_bin=1 is the first bin: [0:bin_length]
    real (db) :: rr
    real (db), dimension(ndim) :: diff_vec   
		
    r_max = g_s_hist%bin_length * size(g_s_hist%val)
    rr_max = r_max * r_max
    
    n_tot = size(org_r(:,1))   ! number of atoms
      
    g_s_hist%val = 0  
    rr_sum = 0.0
  
    do i1 = 1, n_tot
      diff_vec = act_r(i1,:) - org_r(i1,:)
 
      ! The command below does not actually seem to make that much difference
      ! to the execution speed 
                    
      rr = sum(diff_vec*diff_vec)
      
      if (rr < rr_max) then
        which_bin = ceiling(sqrt(rr)/g_s_hist%bin_length) 
        g_s_hist%val(which_bin) = g_s_hist%val(which_bin) + 1
      end if
      
      rr_sum = rr_sum + rr
    end do      
  
  end function cal_g_s_histogram_no_wrap_around
    

  subroutine cal_g_d_histogram(g_d_hist, org_r, act_r, box_edges)
    type (histogram_cutdown), intent(inout) :: g_d_hist
    real(db), dimension(:,:), intent(in) :: org_r, act_r
    real (db), dimension(ndim), intent(in) :: box_edges
    real(db) :: rr_sum
  
    real(db) :: rr_max, r_max
    integer :: n_tot        ! number of atoms
    integer :: i, i1, i2
    integer :: which_bin    ! where which_bin=1 is the first bin: [0:bin_length]
    real (db) :: rr
    real (db), dimension(ndim) :: diff_vec   
		
    r_max = g_d_hist%bin_length * size(g_d_hist%val)
    rr_max = r_max * r_max
    
    n_tot = size(org_r(:,1))   ! number of atoms
      
    g_d_hist%val = 0  
    rr_sum = 0.0
  
    do i1 = 1, n_tot
      do i2 = 1, n_tot
        if (i1 /= i2) then
          diff_vec = act_r(i1,:) - org_r(i2,:)
      
          ! The command below does not actually seem to make that much difference
          ! to the execution speed 
          
          call apply_boundary_condition_to_vector_expensive(diff_vec, box_edges)
          
          rr = sum(diff_vec*diff_vec)
      
          if (rr < rr_max) then
            which_bin = ceiling(sqrt(rr)/g_d_hist%bin_length) 
            g_d_hist%val(which_bin) = g_d_hist%val(which_bin) + 1
          end if
        end if
  		end do  
    end do      
  
  end subroutine cal_g_d_histogram    

end module histogram_cutdown_class