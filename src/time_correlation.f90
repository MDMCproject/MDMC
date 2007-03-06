module time_correlation_class
use various_constants_class         
use structure_class
use tic_toc_class

implicit none

  public :: init_time_correlation, do_time_correlation
  public :: set_n_buffer_average_over, set_n_time_buffers
  public :: print_g_d, print_g_s, print_einstein_diffuse_exp
  public :: clear_time_correlation

  private :: cal_g_s_histogram, cal_g_d_histogram
  private :: update_act_r 

  private :: make_histogram_cutdown


  type histogram_cutdown
    real(db) :: bin_length   
    integer, dimension(:), allocatable :: val
  end type histogram_cutdown
  
  private histogram_cutdown


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


  character(len=27), parameter, private :: filename_prefix1 = "output/einstein_diffuse_exp"  
  integer, private :: filname_number1 = 1  ! so first saved rdf file will be called einstein_diffuse_exp1.xml
  
  character(len=10), parameter, private :: filename_prefix2 = "output/g_s"  
  integer, private :: filname_number2 = 1  
  
  character(len=10), parameter, private :: filename_prefix3 = "output/g_d"  
  integer, private :: filname_number3 = 1    

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

  subroutine print_g_s(temperature, density, time_step)
    use flib_wxml
    real(db), optional, intent(in) :: temperature, density, time_step
    
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
    
    
    g_prefac = g_prefac / (n_buffer_average_over*density)

    do i = 1, n_eval_times    
      do i_bin = 1, n_bin
        call xml_NewElement(xf, "G-s")
        call xml_AddAttribute(xf, "r", str((i_bin-0.5)*bin_length, format="(f10.5)"))
        call xml_AddAttribute(xf, "t", str((i-1)*time_step, format="(f10.5)"))
        call xml_AddAttribute(xf, "G", str(g_s_hists_sum(i)%val(i_bin)*g_prefac(i_bin), format="(f15.5)"))
        call xml_EndElement(xf, "G-s")
      end do 
    end do
    
    g_prefac = g_prefac * (n_buffer_average_over*density)
    
    call xml_EndElement(xf, "G_s-space-time-pair-correlation-function")
    
    call xml_Close(xf)    
  
  end subroutine print_g_s


  ! The temperature is required to passed in dimensional units here.

  subroutine print_g_d(temperature, volume, n_atom, time_step)
    use flib_wxml
    real(db), optional, intent(in) :: temperature, volume, time_step
    integer, optional, intent(in) :: n_atom
    
    type (xmlf_t) :: xf
    integer :: i, n_eval_times, n_bin, i_bin
    real(db) :: bin_length
    character(len=50) :: filename
    real(db) :: density
    
    density = n_atom / volume
    
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
    if (present(temperature) .and. present(volume) .and. present(n_atom)) then
      call xml_AddAttribute(xf, "title", "T = " // trim(str(temperature * T_unit, format="(f10.5)")) // &
                                         " K: rho = " // trim(str(volume, format="(f10.5)")) &
                                         // " atoms/AA-3")
    else if (present(temperature)) then 
      call xml_AddAttribute(xf, "title", "T = " // str(temperature, format="(f10.5)") // &
                                         "K")
    else  if (present(volume) .and. present(n_atom)) then 
      call xml_AddAttribute(xf, "title", "rho = " // str(volume, format="(f10.5)") &
                                         // "atoms/AA-3")
    end if
    
    call xml_AddAttribute(xf, "bin-length", str(bin_length, format="(f10.5)"))
    
    call xml_AddAttribute(xf, "time-unit", "10^-13 s")
    call xml_AddAttribute(xf, "r-units", "AA")
    
    
    call xml_NewElement(xf, "this-file-was-created")
    call xml_AddAttribute(xf, "when", get_current_date_and_time())
    call xml_EndElement(xf, "this-file-was-created")
    
    
    g_prefac = g_prefac / (n_buffer_average_over*density*(n_atom-1))

    do i = 1, n_eval_times    
      do i_bin = 1, n_bin
        call xml_NewElement(xf, "G-d")
        call xml_AddAttribute(xf, "r", str((i_bin-0.5)*bin_length, format="(f10.5)"))
        call xml_AddAttribute(xf, "t", str((i-1)*time_step, format="(f10.5)"))
        call xml_AddAttribute(xf, "G", str(g_d_hists_sum(i)%val(i_bin)*g_prefac(i_bin), format="(f15.5)"))
        call xml_EndElement(xf, "G-d")
      end do 
    end do
    
    g_prefac = g_prefac * (n_buffer_average_over*density*(n_atom-1))
    
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
    allocate(einstein_diffuse_exp(n_time_evals))

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
    
    allocate(g_s_hists_sum(n_time_evals))
    allocate(g_d_hists_sum(n_time_evals))
    
    do i = 1, n_time_evals
      g_s_hists_sum(i) = make_histogram_cutdown(r_max, bin_length)
      g_d_hists_sum(i) = make_histogram_cutdown(r_max, bin_length)
    end do
    
    
    n_bin = size(g_s_hists_sum(1)%val)
    
    
    allocate(g_prefac(n_bin))
    
    temp = 1.0 / ( 4.0 * pi_value * bin_length**3 )
    
    do i = 1, n_bin
      r_i = i - 0.5
      g_prefac(i) = temp / (r_i**2)
      
      !write (*,*) i, g_prefac(i)
    end do    
    !stop
    
    call clear_time_correlation(n_time_evals)
    
  end subroutine init_time_correlation
   
  
  subroutine clear_time_correlation(n_time_evals)
    integer, intent(in) :: n_time_evals
  
    integer :: i
  
  
    do i = 0, n_time_buffers-1
      bufs(i+1)%time_count = - i * n_time_evals / n_time_buffers
    end do  
  
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
    
end module time_correlation_class
