module time_corr_algorithm_class
use various_constants_class         
use structure_class
use tic_toc_class
use histogram_cutdown_class
use time_corr_hist_container_class
use phasespace_class

implicit none

  public :: cal_time_corr_container
  public :: set_n_g_r_t_to_average_over, set_n_buffers  ! used in handler only

  private :: update_act_r 
  private :: allocate_buffer, reset_buffer
  private :: do_time_correlation


  type buffer
    ! time_count may be negative, but buffers are only calculated
    ! for time_count values: 0,1,2,...,N-1, where N equals the XML 
    ! element <n-time-bin> 
    
    integer :: time_count
 
 
    ! org_r holds the atomic coordinates at time_count=0 and
    ! act_r at any later time. 
    ! Both have dimension n_atoms x ndim
    
    real(db), dimension(:,:), allocatable :: org_r, act_r
    
    
    ! sum of square differences between org_r and act_r for
    ! all the atoms. This array has dimension <n-time-bin>
    
    real(db), dimension(:), allocatable :: sum_square_diffs
   
    
    ! histograms for g_s and g_d for various values of time_count.
    ! These arrays have dimension <n-time-bin>, i.e. that many histograms.
    ! So effectively g_s and g_d are functions of r and t.
    
    type (histogram_cutdown), dimension(:), allocatable :: g_s_hists
    type (histogram_cutdown), dimension(:), allocatable :: g_d_hists
  end type buffer
  
  private buffer


  ! bufs holds all the buffers

  type (buffer), dimension(:), allocatable, private :: bufs  
    
    
  integer, private :: num_buffs_cal_thus_far

  integer, private :: n_g_r_t_to_average_over
  integer, private :: n_buffers

contains

  ! calculate the average of n_g_r_t_to_average_over g_s_hist(r,t)'s and g_d_hist(r,t)'s

  subroutine cal_time_corr_container(container, ps, f_list, md_per_time_bin, time_step)
    type (time_corr_hist_container), intent(inout) :: container
    type (phasespace), intent(inout) :: ps
    type (func_list), intent(in) :: f_list
    integer, intent(in) :: md_per_time_bin
    real(db), intent(in) :: time_step
    
    integer n_time_bin

    call check_if_time_corr_hist_container_is_allocated(container)

    n_time_bin = get_time_corr_hist_n_time_bin(container)

    if (allocated(bufs) == .false.) then
      call allocate_buffer(n_time_bin, &
                           size(ps%str%atoms), &
                           get_time_corr_hist_r_max(container), &
                           get_time_corr_hist_r_bin(container))
    end if
    
    call reset_buffer(n_time_bin)
    
    call check_if_time_corr_hist_container_is_allocated(container)


!    character(len=80) :: filename
    
    ! here keep on calling do_time_correlation until enough buffers have 
    ! been calculated
    
    do while (do_time_correlation(container, ps%str) == .false.)
    
! TO PRINT OUT FOR MOVIE - UNCOMMENT STUFF BELOW    
!      if (filename_number < 500) then
!      
!        filename = ""
!        if (filename_number < 10) then
!          write(filename, '(i1)') filename_number
!        else if (filename_number < 100) then
!          write(filename, '(i2)') filename_number
!        else if (filename_number < 1000) then
!          write(filename, '(i3)') filename_number
!        else
!          write(*,*) "ERROR: in cal_full_time_correlation"
!          stop
!        end if
!        
!        filename = md_prefix // trim(filename) // ".xml"
!        filename_number = filename_number + 1
!    
!        call save_structure(ps%str, filename)
!      
!      end if
    
    
      call trajectory_in_phasespace(ps, f_list, md_per_time_bin, time_step)
    end do
     
  end subroutine cal_time_corr_container



  subroutine set_n_buffers(n)
    integer, intent(in) :: n
    
    n_buffers = n
  end subroutine set_n_buffers



  subroutine set_n_g_r_t_to_average_over(n)
    integer, intent(in) :: n
    
    n_g_r_t_to_average_over = n
  end subroutine set_n_g_r_t_to_average_over

  
  
  ! update act_r to str%r but without wrap-around

  subroutine update_act_r(act_r, str)
    real(db), dimension(:,:), intent(inout) :: act_r
    type (structure), intent(in) :: str

    integer :: i

    do i = 1, size(str%r(:,1))
    
      ! use formula (5.2.5) in Rapaport book
      
      act_r(i,:) = str%r(i,:) + &
        nint( (act_r(i,:) - str%r(i,:))/str%box_edges ) * str%box_edges
    end do
    
  end subroutine update_act_r


  ! This function updates the bufs container after the MD has moved the atoms
  ! along by md_per_time_bin*time_step.
  ! It returns .true. when enough data have been simulated to calculate the average
  ! of n_g_r_t_to_average_over g_s_hist(r,t)'s and g_d_hist(r,t)'s. 

  function do_time_correlation(container, str) result(job_done)
    type (time_corr_hist_container), intent(inout) :: container
    type (structure), intent(in) :: str
    
    logical :: job_done
  
  
    integer :: i_buf, i_time
    integer :: tc    ! shorthand for time_count
    integer :: n_time_bin
    real(db) :: fac
   
    
    ! bufs holds all the buffers
    
    n_time_bin = size(bufs(1)%sum_square_diffs)
  
  
    ! sum over all the buffers
    
    do i_buf = 1, size(bufs)
      tc = bufs(i_buf)%time_count
      
      if (tc == 0) then
        bufs(i_buf)%org_r = str%r
        bufs(i_buf)%act_r = str%r


        ! calculate g_d(r,t=0) = g(r)
        
        call cal_g_d_histogram(bufs(i_buf)%g_d_hists(tc+1), & 
          bufs(i_buf)%org_r, bufs(i_buf)%act_r, str%box_edges)
                   
      elseif (tc > 0) then   
      
        ! update actual act_r positions relative to org_r 
           
        call update_act_r(bufs(i_buf)%act_r, str)
        
        
        ! calculate g_s(r,t) and the sum of square differences
        
        bufs(i_buf)%sum_square_diffs(tc+1) = cal_g_s_histogram_no_wrap_around( & 
          bufs(i_buf)%g_s_hists(tc+1), bufs(i_buf)%org_r, bufs(i_buf)%act_r, str%box_edges)
        
        
        ! calculate g_d(r,t)
          
        call cal_g_d_histogram(bufs(i_buf)%g_d_hists(tc+1), bufs(i_buf)%org_r, bufs(i_buf)%act_r, str%box_edges)  
      end if
      
      
      ! check if this buffer reached time series end point
      
      if (tc + 1 == n_time_bin) then
      
        ! reset time_count
        
        bufs(i_buf)%time_count = 0
      
        
        ! notice at this point einstein_diffuse_exp holds the sum of 'sum
        ! of square differences' 
        
        container%einstein_diffuse_exp = container%einstein_diffuse_exp + bufs(i_buf)%sum_square_diffs
        
        
        ! sum up G_s(r,t) and G_d(r,t) histograms for each n_time_bin
        
        do i_time = 2, n_time_bin
          container%g_s_hists_sum(i_time)%val = container%g_s_hists_sum(i_time)%val + bufs(i_buf)%g_s_hists(i_time)%val
        end do

        do i_time = 1, n_time_bin
          container%g_d_hists_sum(i_time)%val = container%g_d_hists_sum(i_time)%val + bufs(i_buf)%g_d_hists(i_time)%val
        end do
        
        
        ! number of buffers calculated thus far
        
        num_buffs_cal_thus_far = num_buffs_cal_thus_far + 1
        
        print *, "num_buffs_cal_thus_far ", num_buffs_cal_thus_far
        
        container%n_accum = num_buffs_cal_thus_far
        
        
        ! check if enough buffers have been completed to calculate average
        
        if (num_buffs_cal_thus_far == n_g_r_t_to_average_over) then
        
          ! calculate einstein diff. expression according to (5.2.4) Rapaport book
        
          fac = 1 / (ndim*2*container%time_bin*n_g_r_t_to_average_over*get_n_atom(str))
          do i_time = 2, n_time_bin             
            container%einstein_diffuse_exp(i_time) = container%einstein_diffuse_exp(i_time) * fac / (i_time-1)
          end do           
          
          
          ! finish what needs to be calculated
          
          job_done = .true.
          return
        end if
      
      else
        bufs(i_buf)%time_count = bufs(i_buf)%time_count + 1
      end if
    
    end do ! end i_buf
    
    
    ! have not completed enough buffers yet
    
    job_done = .false.
    
  end function do_time_correlation

  
  ! Allocate bufs, which has dimensions n_buffers * number of time bins * number of r bins.
  ! It also has dimensions n_buffers * n_atoms * ndim.

  subroutine allocate_buffer(n_time_bin, n_atoms, r_max, bin_length)
    integer, intent(in) :: n_time_bin, n_atoms
    real(db), intent(in) :: r_max, bin_length
 
    integer :: i, j, n_bin
    real(db) :: r_i, temp


    allocate(bufs(n_buffers))

    do i = 0, n_buffers-1
      
      allocate(bufs(i+1)%org_r(n_atoms, ndim)) 
      allocate(bufs(i+1)%act_r(n_atoms, ndim))
      allocate(bufs(i+1)%sum_square_diffs(n_time_bin))
      
      allocate(bufs(i+1)%g_s_hists(n_time_bin))
      allocate(bufs(i+1)%g_d_hists(n_time_bin))
      
      do j = 1, n_time_bin
        bufs(i+1)%g_s_hists(j) = make_histogram_cutdown(r_max, bin_length)
        bufs(i+1)%g_d_hists(j) = make_histogram_cutdown(r_max, bin_length)
      end do
      
    end do
    
  end subroutine allocate_buffer
   

  
  subroutine reset_buffer(n_time_bin)
    integer, intent(in) :: n_time_bin
  
    integer :: i
  
    do i = 0, n_buffers-1
      bufs(i+1)%time_count = - i * n_time_bin / n_buffers
    end do  
  
    num_buffs_cal_thus_far = 0

  end subroutine reset_buffer
    
end module time_corr_algorithm_class
