module time_correlation_class
use various_constants_class         
use structure_class
use tic_toc_class

implicit none

  public :: init_time_correlation, do_time_correlation

  private :: cal_g_s_histogram
  private :: clear_time_correlation
  private :: update_act_r
  private :: print_g_s, print_einstein_diffuse_exp

  private :: make_histogram_cutdown


  type histogram_cutdown
    real(db) :: bin_length   
    integer, dimension(:), allocatable :: val
  end type histogram_cutdown


  type buffer
    integer :: time_count
    real(db), dimension(:,:), allocatable :: org_r, act_r
    real(db), dimension(:), allocatable :: sum_square_diffs
    
    type (histogram_cutdown), dimension(:), allocatable :: g_s_hists
    type (histogram_cutdown), dimension(:), allocatable :: g_d_hists
  end type buffer


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
    
  integer, private :: num_buffs_cal_thus_far


  character(len=27), parameter, private :: filename_prefix1 = "output/einstein_diffuse_exp"  
  integer, private :: filname_number1 = 1  ! so first saved rdf file will be called einstein_diffuse_exp1.xml
  
  character(len=10), parameter, private :: filename_prefix2 = "output/g_s"  
  integer, private :: filname_number2 = 1  
  
  character(len=10), parameter, private :: filename_prefix3 = "output/g_d"  
  integer, private :: filname_number3 = 1    


contains

  ! The temperature is required to passed in dimensional units here.

  subroutine print_einstein_diffuse_exp(temperature, density, time_step)
    use flib_wxml
    real(db), optional, intent(in) :: temperature, density, time_step
    
    type (xmlf_t) :: xf
    integer :: i, n_eval_times
    character(len=50) :: filename
    
    if (filname_number1 < 10) then
      write(filename, '(i1)') filname_number1
    else if (filname_number1 < 100) then
      write(filename, '(i2)') filname_number1
    else if (filname_number1 < 1000) then
      write(filename, '(i3)') filname_number1
    else
      write(*,*) "ERROR: in save_rdf"
      write(*,*) "It is assumed that you did not intend to write"
      write(*,*) "to disk 1000 rdf xml files!!!!"
      stop
    end if
    
    filname_number1 = filname_number1 + 1
    
    filename = filename_prefix1 // trim(filename) // ".xml"
    
    write(*,'(3a)') "Write ", trim(filename), " to disk"
    
    n_eval_times = size(einstein_diffuse_exp)    
    
    call xml_OpenFile(filename, xf, indent=.true.)
    
    call xml_AddXMLDeclaration(xf, "UTF-8")
    call xml_NewElement(xf, "einstein-diffuse-exp")
    
    ! notice convert units of temperature from dimensionless to K  
    if (present(temperature) .and. present(density)) then
      call xml_AddAttribute(xf, "title", "T = " // trim(str(temperature * T_unit, format="(f10.5)")) // &
                                         " K: rho = " // trim(str(density, format="(f10.5)")) &
                                         // " atoms/AA-3")
    else if (present(temperature)) then 
      call xml_AddAttribute(xf, "title", "T = " // str(temperature, format="(f10.5)") // &
                                         "K")
    else  if (present(density)) then 
      call xml_AddAttribute(xf, "title", "rho = " // str(density, format="(f10.5)") &
                                         // "atoms/AA-3")
    end if
    
    call xml_AddAttribute(xf, "time-unit", "10^-13 s")
    call xml_AddAttribute(xf, "diffuse-units", "10^13 AA^2 s^-1")
    
    
    call xml_NewElement(xf, "this-file-was-created")
    call xml_AddAttribute(xf, "when", get_current_date_and_time())
    call xml_EndElement(xf, "this-file-was-created")
    
    do i = 1, n_eval_times
      call xml_NewElement(xf, "D-of-t")
      call xml_AddAttribute(xf, "t", str((i-1)*time_step, format="(f10.5)"))
      call xml_AddAttribute(xf, "D", str(einstein_diffuse_exp(i), format="(f10.5)"))
      call xml_EndElement(xf, "D-of-t")
    end do 
    
    call xml_EndElement(xf, "einstein-diffuse-exp")
    
    call xml_Close(xf)    
  
  end subroutine print_einstein_diffuse_exp


  ! The temperature is required to passed in dimensional units here.

  subroutine print_g_s(temperature, density, time_step, n_buffer_average_over)
    use flib_wxml
    real(db), optional, intent(in) :: temperature, density, time_step
    integer, intent(in) :: n_buffer_average_over
    
    type (xmlf_t) :: xf
    integer :: i, n_eval_times, n_bin, i_bin
    real(db) :: bin_length
    character(len=50) :: filename
    
    if (filname_number2 < 10) then
      write(filename, '(i1)') filname_number2
    else if (filname_number2 < 100) then
      write(filename, '(i2)') filname_number2
    else if (filname_number2 < 1000) then
      write(filename, '(i3)') filname_number2
    else
      write(*,*) "ERROR: in save_rdf"
      write(*,*) "It is assumed that you did not intend to write"
      write(*,*) "to disk 1000 rdf xml files!!!!"
      stop
    end if
    
    filname_number2 = filname_number2 + 1
    
    filename = filename_prefix2 // trim(filename) // ".xml"
    
    write(*,'(3a)') "Write ", trim(filename), " to disk"
    
    
    n_eval_times = size(einstein_diffuse_exp)    
    n_bin = size(g_s_hists_sum(1)%val)
    bin_length = g_s_hists_sum(1)%bin_length
    
    call xml_OpenFile(filename, xf, indent=.true.)
    
    call xml_AddXMLDeclaration(xf, "UTF-8")
    call xml_NewElement(xf, "G_s-space-time-pair-correlation-function")
    
    ! notice convert units of temperature from dimensionless to K  
    if (present(temperature) .and. present(density)) then
      call xml_AddAttribute(xf, "title", "T = " // trim(str(temperature * T_unit, format="(f10.5)")) // &
                                         " K: rho = " // trim(str(density, format="(f10.5)")) &
                                         // " atoms/AA-3")
    else if (present(temperature)) then 
      call xml_AddAttribute(xf, "title", "T = " // str(temperature, format="(f10.5)") // &
                                         "K")
    else  if (present(density)) then 
      call xml_AddAttribute(xf, "title", "rho = " // str(density, format="(f10.5)") &
                                         // "atoms/AA-3")
    end if
    
    call xml_AddAttribute(xf, "bin-length", str(bin_length, format="(f10.5)"))
    
    call xml_AddAttribute(xf, "time-unit", "10^-13 s")
    call xml_AddAttribute(xf, "r-units", "AA")
    
    
    call xml_NewElement(xf, "this-file-was-created")
    call xml_AddAttribute(xf, "when", get_current_date_and_time())
    call xml_EndElement(xf, "this-file-was-created")
    
    
    g_prefac = g_prefac / n_buffer_average_over

    do i = 1, n_eval_times    
      do i_bin = 1, n_bin
        call xml_NewElement(xf, "G-s")
        call xml_AddAttribute(xf, "r", str((i_bin-0.5)*bin_length, format="(f10.5)"))
        call xml_AddAttribute(xf, "t", str((i-1)*time_step, format="(f10.5)"))
        call xml_AddAttribute(xf, "G", str(g_s_hists_sum(i)%val(i_bin)*g_prefac(i_bin), format="(f15.5)"))
        call xml_EndElement(xf, "G-s")
      end do 
    end do
    
    g_prefac = g_prefac * n_buffer_average_over
    
    call xml_EndElement(xf, "G_s-space-time-pair-correlation-function")
    
    call xml_Close(xf)    
  
  end subroutine print_g_s


  ! The temperature is required to passed in dimensional units here.

  subroutine print_g_d(temperature, density, time_step, n_buffer_average_over)
    use flib_wxml
    real(db), optional, intent(in) :: temperature, density, time_step
    integer, intent(in) :: n_buffer_average_over
    
    type (xmlf_t) :: xf
    integer :: i, n_eval_times, n_bin, i_bin
    real(db) :: bin_length
    character(len=50) :: filename
    
    if (filname_number3 < 10) then
      write(filename, '(i1)') filname_number3
    else if (filname_number3 < 100) then
      write(filename, '(i2)') filname_number3
    else if (filname_number3 < 1000) then
      write(filename, '(i3)') filname_number3
    else
      write(*,*) "ERROR: in save_rdf"
      write(*,*) "It is assumed that you did not intend to write"
      write(*,*) "to disk 1000 rdf xml files!!!!"
      stop
    end if
    
    filname_number3 = filname_number3 + 1
    
    filename = filename_prefix3 // trim(filename) // ".xml"
    
    write(*,'(3a)') "Write ", trim(filename), " to disk"
    
    
    n_eval_times = size(einstein_diffuse_exp)    
    n_bin = size(g_s_hists_sum(1)%val)
    bin_length = g_s_hists_sum(1)%bin_length
    
    call xml_OpenFile(filename, xf, indent=.true.)
    
    call xml_AddXMLDeclaration(xf, "UTF-8")
    call xml_NewElement(xf, "G_d-space-time-pair-correlation-function")
    
    ! notice convert units of temperature from dimensionless to K  
    if (present(temperature) .and. present(density)) then
      call xml_AddAttribute(xf, "title", "T = " // trim(str(temperature * T_unit, format="(f10.5)")) // &
                                         " K: rho = " // trim(str(density, format="(f10.5)")) &
                                         // " atoms/AA-3")
    else if (present(temperature)) then 
      call xml_AddAttribute(xf, "title", "T = " // str(temperature, format="(f10.5)") // &
                                         "K")
    else  if (present(density)) then 
      call xml_AddAttribute(xf, "title", "rho = " // str(density, format="(f10.5)") &
                                         // "atoms/AA-3")
    end if
    
    call xml_AddAttribute(xf, "bin-length", str(bin_length, format="(f10.5)"))
    
    call xml_AddAttribute(xf, "time-unit", "10^-13 s")
    call xml_AddAttribute(xf, "r-units", "AA")
    
    
    call xml_NewElement(xf, "this-file-was-created")
    call xml_AddAttribute(xf, "when", get_current_date_and_time())
    call xml_EndElement(xf, "this-file-was-created")
    
    
    g_prefac = g_prefac / n_buffer_average_over

    do i = 1, n_eval_times    
      do i_bin = 1, n_bin
        call xml_NewElement(xf, "G-d")
        call xml_AddAttribute(xf, "r", str((i_bin-0.5)*bin_length, format="(f10.5)"))
        call xml_AddAttribute(xf, "t", str((i-1)*time_step, format="(f10.5)"))
        call xml_AddAttribute(xf, "G", str(g_d_hists_sum(i)%val(i_bin)*g_prefac(i_bin), format="(f15.5)"))
        call xml_EndElement(xf, "G-d")
      end do 
    end do
    
    g_prefac = g_prefac * n_buffer_average_over
    
    call xml_EndElement(xf, "G_d-space-time-pair-correlation-function")
    
    call xml_Close(xf)    
  
  end subroutine print_g_d  
  
  
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

  subroutine do_time_correlation(str, n_buffer_average_over, time_step, temperature, density)
    type (structure), intent(in) :: str
    integer, intent(in) :: n_buffer_average_over
    real(db), intent(in) :: time_step
    real(db), optional, intent(in) :: temperature, density    
  
    integer :: i_buf, i_time
    integer :: tc    ! shorthand for time_count
    integer :: n_time_evals
    real(db) :: fac
    
    n_time_evals = size(bufs(1)%sum_square_diffs)
  
    do i_buf = 1, size(bufs)
      tc = bufs(i_buf)%time_count
      
      if (tc == 0) then
        bufs(i_buf)%org_r = str%r
        bufs(i_buf)%act_r = str%r
        
        call cal_g_d_histogram(bufs(i_buf)%g_d_hists(tc+1), & 
          bufs(i_buf)%org_r, bufs(i_buf)%act_r )         
      elseif (tc > 0) then      
        call update_act_r(bufs(i_buf)%act_r, str)
        
        bufs(i_buf)%sum_square_diffs(tc+1) = cal_g_s_histogram( & 
          bufs(i_buf)%g_s_hists(tc+1), bufs(i_buf)%org_r, bufs(i_buf)%act_r )
          
        call cal_g_d_histogram(bufs(i_buf)%g_d_hists(tc+1), bufs(i_buf)%org_r, bufs(i_buf)%act_r )  
      end if
      
      ! check if this buffer reached time series end point
      
      if (tc + 1 == n_time_evals) then
        bufs(i_buf)%time_count = 0
        
        ! notice at this point einstein_diffuse_exp holds just the sum of sums
        ! of square differences 
        
        einstein_diffuse_exp = einstein_diffuse_exp + bufs(i_buf)%sum_square_diffs
        
        
        ! sum up histograms for the calculation of G_s(r,t)
        
        do i_time = 2, n_time_evals
          g_s_hists_sum(i_time)%val = g_s_hists_sum(i_time)%val + bufs(i_buf)%g_s_hists(i_time)%val
        end do

        do i_time = 1, n_time_evals
          g_d_hists_sum(i_time)%val = g_d_hists_sum(i_time)%val + bufs(i_buf)%g_d_hists(i_time)%val
        end do
        
        
        num_buffs_cal_thus_far = num_buffs_cal_thus_far + 1
        
        print *, "num_buffs_cal_thus_far ", num_buffs_cal_thus_far
        
        
        ! check if enough buffers have been completed to calculate averages
        
        if (num_buffs_cal_thus_far == n_buffer_average_over) then
          fac = 1 / (ndim*2*time_step*n_buffer_average_over)
          do i_time = 2, n_time_evals
          
            ! expression (5.2.4) Rapaport book 
             
            einstein_diffuse_exp(i_time) = einstein_diffuse_exp(i_time) * fac / (i_time-1)
          end do 
          call print_einstein_diffuse_exp(temperature, density, time_step)
          
          call print_g_s(temperature, density, time_step, n_buffer_average_over)
          call print_g_d(temperature, density, time_step, n_buffer_average_over)
          
          call clear_time_correlation()
        end if
      
      else
        bufs(i_buf)%time_count = bufs(i_buf)%time_count + 1
      end if
    
    end do ! end i_buf
    
    
  
  end subroutine do_time_correlation


  subroutine init_time_correlation(n_time_buffers, n_time_evals, n_atoms, r_max, bin_length)
    integer, intent(in) :: n_time_buffers, n_time_evals, n_atoms
    real(db), intent(in) :: r_max, bin_length
 
    integer :: i, j, n_bin
    real(db) :: r_i, temp


    allocate(bufs(n_time_buffers))
    allocate(einstein_diffuse_exp(n_time_evals))

    do i = 0, n_time_buffers-1
      bufs(i+1)%time_count = - i * n_time_evals / n_time_buffers
      
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
    
    allocate(g_s_hists_sum(n_time_evals))
    allocate(g_d_hists_sum(n_time_evals))
    
    do i = 1, n_time_evals
      g_s_hists_sum(i) = make_histogram_cutdown(r_max, bin_length)
      g_d_hists_sum(i) = make_histogram_cutdown(r_max, bin_length)
    end do
    
    
    n_bin = size(g_s_hists_sum(1)%val)
    
    !allocate(g_s_pair_func(n_bin, n_time_evals))
    
    
    allocate(g_prefac(n_bin))
    
    temp = 1.0 / ( 4.0 * pi_value * n_atoms * bin_length**3 )
    
    do i = 1, n_bin
      r_i = i - 0.5
      g_prefac(i) = temp / (r_i**2)
    end do    
    
    
    call clear_time_correlation()
    
  end subroutine init_time_correlation
   
  
  subroutine clear_time_correlation()
    integer :: i
  
    num_buffs_cal_thus_far = 0
    einstein_diffuse_exp = 0.0
    
    do i = 1, size(g_s_hists_sum)
      g_s_hists_sum(i)%val = 0
      g_d_hists_sum(i)%val = 0
    end do 
    
  end subroutine clear_time_correlation
  
  
  function make_histogram_cutdown(r_max, bin_length) result(hist)
    real(db), intent(in) :: r_max, bin_length
    type (histogram_cutdown) :: hist
    
    hist%bin_length = bin_length 
    allocate(hist%val(floor(r_max/bin_length)))
  end function make_histogram_cutdown
  

  function cal_g_s_histogram(g_s_hist, org_r, act_r) result(rr_sum)
    type (histogram_cutdown), intent(inout) :: g_s_hist
    real(db), dimension(:,:), intent(in) :: org_r, act_r
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
      
      rr = sum(diff_vec*diff_vec)
      
      if (rr < rr_max) then
        which_bin = ceiling(sqrt(rr)/g_s_hist%bin_length) 
        g_s_hist%val(which_bin) = g_s_hist%val(which_bin) + 1
      end if
      
      rr_sum = rr_sum + rr
  		  
    end do      
  
  end function cal_g_s_histogram
    

  subroutine cal_g_d_histogram(g_d_hist, org_r, act_r)
    type (histogram_cutdown), intent(inout) :: g_d_hist
    real(db), dimension(:,:), intent(in) :: org_r, act_r
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
      
          rr = sum(diff_vec*diff_vec)
      
          if (rr < rr_max) then
            which_bin = ceiling(sqrt(rr)/g_d_hist%bin_length) 
            g_d_hist%val(which_bin) = g_d_hist%val(which_bin) + 1
          end if
        end if
  		end do  
    end do      
  
  end subroutine cal_g_d_histogram    
    
end module time_correlation_class
