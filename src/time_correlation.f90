module time_correlation_class
use various_constants_class         
use structure_class
use tic_toc_class
use histogram_cutdown_class

implicit none

  public :: init_time_correlation, do_time_correlation
  public :: set_n_buffer_average_over, set_n_time_buffers
  public :: print_g_d, print_g_s, print_einstein_diffuse_exp
  public :: clear_time_correlation

  private :: update_act_r 




  type buffer
    ! time_count may be negative, but buffers are only calculated
    ! for time_count values: 0,1,2,...,N-1, where N equals the XML 
    ! element <n-time-evals> 
    
    integer :: time_count
 
 
    ! org_r holds the atomic coordinates at time_count=0 and
    ! act_r at any later time. 
    ! Both have dimension n_atoms x ndim
    
    real(db), dimension(:,:), allocatable :: org_r, act_r
    
    
    ! sum of square differences between org_r and act_r for
    ! all the atoms. This array has dimension <n-time-evals>
    
    real(db), dimension(:), allocatable :: sum_square_diffs
   
    
    ! histograms for g_s and g_d for various values of time_count.
    ! These arrayes have dimension <n-time-evals>
    
    type (histogram_cutdown), dimension(:), allocatable :: g_s_hists
    type (histogram_cutdown), dimension(:), allocatable :: g_d_hists
  end type buffer
  
  private buffer


  ! bufs holds all the buffers

  type (buffer), dimension(:), allocatable, private :: bufs  
    
  real(db), dimension(:), allocatable, private :: einstein_diffuse_exp
  
  
  ! g_s_hists_sum used to hold the sum of histograms and g_s_pair_func to hold
  ! the G_s(r,t) space-time pair correlation function. g_s_pair_func has dimensions
  ! n_bin x n_time_evals. g_prefac is to convert from histogram values to g values.
  
  type (histogram_cutdown), dimension(:), allocatable, private :: g_s_hists_sum
  type (histogram_cutdown), dimension(:), allocatable, private :: g_d_hists_sum
  !real(db), dimension(:,:,:), allocatable, private :: g_s_pair_func
  real(db), dimension(:), allocatable, private :: g_prefac 
  
  
  ! Temporary container for reading in G_d dataset from file
  
  real(db), dimension(:,:), allocatable :: g_d_data
    
    
  integer, private :: num_buffs_cal_thus_far




  integer, private :: n_buffer_average_over
  integer, private :: n_time_buffers

contains


  function g_d_fom_val() result(val)
    real(db) :: val
    
    integer :: i, i_bin
    
    integer n_eval_times, n_bin
    
    n_eval_times = size(g_d_data, 2)
    n_bin = size(g_d_data, 1)
  
    val = 0.0
    
    g_prefac = g_prefac / n_buffer_average_over

    do i = 1, n_eval_times    
      do i_bin = 20, n_bin
         val = val + (g_d_data(i_bin, i) - g_d_hists_sum(i)%val(i_bin)*g_prefac(i_bin))**2
      end do 
    end do
    
    g_prefac = g_prefac * n_buffer_average_over  
  
  end function g_d_fom_val


  subroutine set_n_time_buffers(n)
    integer, intent(in) :: n
    
    n_time_buffers = n
  end subroutine set_n_time_buffers

  subroutine set_n_buffer_average_over(n)
    integer, intent(in) :: n
    
    n_buffer_average_over = n
  end subroutine set_n_buffer_average_over


 
  
  
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



  ! time_step is here equal to real MD delta_t times n_delta_t

  function do_time_correlation(str, time_step) result(job_done)
    type (structure), intent(in) :: str
    real(db), intent(in) :: time_step   
    
    logical :: job_done
  
  
    integer :: i_buf, i_time
    integer :: tc    ! shorthand for time_count
    integer :: n_time_evals
    real(db) :: fac
   
    
    ! bufs holds all the buffers
    
    n_time_evals = size(bufs(1)%sum_square_diffs)
  
  
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
        
        bufs(i_buf)%sum_square_diffs(tc+1) = cal_g_s_histogram( & 
          bufs(i_buf)%g_s_hists(tc+1), bufs(i_buf)%org_r, bufs(i_buf)%act_r, str%box_edges)
        
        
        ! calculate g_d(r,t)
          
        call cal_g_d_histogram(bufs(i_buf)%g_d_hists(tc+1), bufs(i_buf)%org_r, bufs(i_buf)%act_r, str%box_edges)  
      end if
      
      
      ! check if this buffer reached time series end point
      
      if (tc + 1 == n_time_evals) then
      
        ! reset time_count
        
        bufs(i_buf)%time_count = 0
      
        
        ! notice at this point einstein_diffuse_exp holds the sum of 'sum
        ! of square differences' 
        
        einstein_diffuse_exp = einstein_diffuse_exp + bufs(i_buf)%sum_square_diffs
        
        
        ! sum up G_s(r,t) and G_d(r,t) histograms for each n_time_evals
        
        do i_time = 2, n_time_evals
          g_s_hists_sum(i_time)%val = g_s_hists_sum(i_time)%val + bufs(i_buf)%g_s_hists(i_time)%val
        end do

        do i_time = 1, n_time_evals
          g_d_hists_sum(i_time)%val = g_d_hists_sum(i_time)%val + bufs(i_buf)%g_d_hists(i_time)%val
        end do
        
        
        ! number of buffers calculated thus far
        
        num_buffs_cal_thus_far = num_buffs_cal_thus_far + 1
        
        print *, "num_buffs_cal_thus_far ", num_buffs_cal_thus_far
        
        
        ! check if enough buffers have been completed to calculate average
        
        if (num_buffs_cal_thus_far == n_buffer_average_over) then
        
          ! calculate einstein diff. expression according to (5.2.4) Rapaport book
        
          fac = 1 / (ndim*2*time_step*n_buffer_average_over)
          do i_time = 2, n_time_evals             
            einstein_diffuse_exp(i_time) = einstein_diffuse_exp(i_time) * fac / (i_time-1)
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
  

  subroutine init_time_correlation(n_time_evals, n_atoms, r_max, bin_length)
    integer, intent(in) :: n_time_evals, n_atoms
    real(db), intent(in) :: r_max, bin_length
 
    integer :: i, j, n_bin
    real(db) :: r_i, temp


    allocate(bufs(n_time_buffers))

    do i = 0, n_time_buffers-1
      
      allocate(bufs(i+1)%org_r(n_atoms, ndim)) 
      allocate(bufs(i+1)%act_r(n_atoms, ndim))
      allocate(bufs(i+1)%sum_square_diffs(n_time_evals))
      
      allocate(bufs(i+1)%g_s_hists(n_time_evals))
      allocate(bufs(i+1)%g_d_hists(n_time_evals))
      
      do j = 1, n_time_evals
        bufs(i+1)%g_s_hists(j) = make_histogram_cutdown(r_max, bin_length)
        bufs(i+1)%g_d_hists(j) = make_histogram_cutdown(r_max, bin_length)
      end do
      
    end do

    
    call clear_time_correlation(n_time_evals)
    
  end subroutine init_time_correlation
   
  
  subroutine clear_time_correlation(n_time_evals)
    integer, intent(in) :: n_time_evals
  
    integer :: i
  
    do i = 0, n_time_buffers-1
      bufs(i+1)%time_count = - i * n_time_evals / n_time_buffers
    end do  
  
    num_buffs_cal_thus_far = 0

  end subroutine clear_time_correlation
  


    
end module time_correlation_class
