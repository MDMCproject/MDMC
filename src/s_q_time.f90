module s_q_time_class
use various_constants_class         
use structure_type_definition_class 
use structure_nn_methods_class
use tic_toc_class
use time_corr_hist_container_class

implicit none

  public :: cal_s_q_time
  public :: make_s_q_time
  public :: check_if_s_q_time_allocated
  public :: print_s_q_time
  
  private :: make_and_cal_integrat_over_r_precal
  private :: is_integrat_over_r_precal_allocated
  
  type s_q_time
    real(db), dimension(:), allocatable :: q
    real(db) :: delta_time  ! this number together with the lenght of
                            ! the time dimension determines the time binning
                            
    real(db), dimension(:,:), allocatable :: self   ! dimensions (n_q, n_time)
    real(db), dimension(:,:), allocatable :: diff   ! dimensions (n_q, n_time)
  end type s_q_time

  type integrat_over_r_precal
    real(db), dimension(:), allocatable :: q  ! used perhaps in the future to figure out
                                             ! if this container needs recalculating
                                             
    real(db), dimension(:,:), allocatable :: val  ! dimension (n_r, n_q)
  end type integrat_over_r_precal
  
  private integrat_over_r_precal
  
  
  type (integrat_over_r_precal), private :: local_integrat_over_r_precal  
  
  
  
  character(len=15), parameter, private :: filename_prefix = "output/s_q_time"    
  integer, private :: filname_number = 1    
  
contains

  function get_s_q_time_n_time(container) result(n_time)
    type(s_q_time), intent(in) :: container
    integer :: n_time 
    
    call check_if_s_q_time_allocated(container)     

    n_time = size(container%self, 2)
  
  end function get_s_q_time_n_time 



  subroutine cal_s_q_time(g_r_t, str, s_q_t)
    type(time_corr_hist_container), intent(in) :: g_r_t
    type(structure), intent(in) :: str
    type(s_q_time), intent(inout) :: s_q_t
  
    integer :: n_q, n_t, n_r, i_q, i_t, i_r, n_atom
    real(db) :: density, pre_s_d, pre_s_s
    real(db), dimension(:), allocatable :: prefac_s_d, prefac_s_s
  
    ! some checks
       
    call check_if_time_corr_hist_container_is_allocated(g_r_t)
    call check_if_s_q_time_allocated(s_q_t)
    
    
    ! Has local_integrat_over_r_precal been allocated and populated
    ! Notice should also compare the Q array in that container with 
    ! Q array in the s_q_t container  
    
    n_r = get_time_corr_hist_n_r_bin(g_r_t)
    
    if (is_integrat_over_r_precal_allocated(local_integrat_over_r_precal) == .false.) then
      local_integrat_over_r_precal = make_and_cal_integrat_over_r_precal( &
        n_r, get_time_corr_hist_r_bin(g_r_t) &
        , s_q_t%q)
    end if
  

    ! calculate 
    
    n_atom = get_n_atom(str)
    density = get_density(str)
    
    n_q = size(s_q_t%q)
    n_t = get_time_corr_hist_n_time_bin(g_r_t)
    
    
    allocate(prefac_s_d(n_r))
    allocate(prefac_s_s(n_r))
    
    prefac_s_d =  g_r_t%volume_prefac / (density*(n_atom-1)*g_r_t%n_accum)
    prefac_s_s =  g_r_t%volume_prefac / (density*g_r_t%n_accum)
    
   
    ! notice the was the calculation below could be done faster
    
    s_q_t%diff = 0
    
    pre_s_d = density * (n_atom-1)/n_atom
    pre_s_s = density / n_atom
    
    do i_q = 1, n_q
      do i_t = 1, n_t
      
        do i_r = 1, n_r
          s_q_t%diff(i_q, i_t) = s_q_t%diff(i_q, i_t) + pre_s_d * &
            (prefac_s_d(i_r)*g_r_t%g_d_hists_sum(i_t)%val(i_r)-1) * &
                            local_integrat_over_r_precal%val(i_r,i_q)
          s_q_t%self(i_q, i_t) = s_q_t%self(i_q, i_t) + pre_s_s * &
             (prefac_s_s(i_r)*g_r_t%g_s_hists_sum(i_t)%val(i_r)) * &
                            local_integrat_over_r_precal%val(i_r,i_q)
        end do                          
      end do
    end do


    ! force the correct exact values for 'self' when t = 0
    
    do i_q = 1, n_q
       s_q_t%self(i_q, 1) = 1.0;
    end do


    deallocate(prefac_s_d)
    deallocate(prefac_s_s)

  end subroutine cal_s_q_time



  ! The temperature is assumed to be in dimensionless units.
  ! density should perhaps also be made an optional argument 

  subroutine print_s_q_time(container, density, temperature)
    use flib_wxml
    type(s_q_time), intent(in) :: container  
    real(db), intent(in) :: density    
    real(db), optional, intent(in) :: temperature
    
    type (xmlf_t) :: xf
    integer :: i_q, i_t, n_time
    
    character(len=50) :: filename
    
    if (filname_number < 10) then
      write(filename, '(i1)') filname_number
    else if (filname_number < 100) then
      write(filename, '(i2)') filname_number
    else if (filname_number < 1000) then
      write(filename, '(i3)') filname_number
    else
      write(*,*) "ERROR: in save_rdf"
      write(*,*) "It is assumed that you did not intend to write"
      write(*,*) "to disk 1000 rdf xml files!!!!"
      stop
    end if
    
    filname_number = filname_number + 1
    
    filename = filename_prefix // trim(filename) // ".xml" ! here filename is just a number on rhs
    
    write(*,'(3a)') "Write ", trim(filename), " to disk"
    
    
    call xml_OpenFile(filename, xf, indent=.true.)
    
    call xml_AddXMLDeclaration(xf, "UTF-8")
    call xml_NewElement(xf, "s-q-time")
    
    ! notice convert units of temperature from dimensionless to K  
    if (present(temperature)) then
      call xml_AddAttribute(xf, "title", "T = " // trim(str(temperature * T_unit, format="(f10.5)")) // &
                                         " K: rho = " // trim(str(density, format="(f10.5)")) &
                                         // " atoms/AA-3")
    else 
      call xml_AddAttribute(xf, "title", "rho = " // str(density, format="(f10.5)") &
                                         // "atoms/AA-3")                                                                        
    end if
    
    !call xml_AddAttribute(xf, "bin-length", str(bin_length, format="(f10.5)"))
    
    call xml_AddAttribute(xf, "q-units", "AA^-1")
    call xml_AddAttribute(xf, "time-unit", "10^-13 s")        
    
    call xml_NewElement(xf, "this-file-was-created")
    call xml_AddAttribute(xf, "when", get_current_date_and_time())
    call xml_EndElement(xf, "this-file-was-created")  
    
    
    n_time = size(container%self, 2)

    do i_t = 1, n_time   
      do i_q = 1, size(container%q)
        call xml_NewElement(xf, "SQt")
        call xml_AddAttribute(xf, "q", str(container%q(i_q), format="(f10.5)"))
        call xml_AddAttribute(xf, "t", str((i_t-1)*container%delta_time, format="(f10.5)"))
        call xml_AddAttribute(xf, "S-self", str(container%self(i_q, i_t), format="(f15.5)"))
        call xml_AddAttribute(xf, "S-diff", str(container%diff(i_q, i_t), format="(f15.5)"))        
        call xml_EndElement(xf, "SQt")
      end do 
    end do
    
    call xml_EndElement(xf, "s-q-time")
    
    call xml_Close(xf)    
  
  end subroutine print_s_q_time


  subroutine check_if_s_q_time_allocated(container)
    type(s_q_time), intent(in) :: container
  
    if (allocated(container%q) == .false.) then
      write(*,*) " "
      write(*,*) "ERROR in s_q_time.f90"
      write(*,*) "Forgot to allocate s_q_time"
      stop
    end if  
  end subroutine check_if_s_q_time_allocated


  function is_integrat_over_r_precal_allocated(container) result(boolean)
    type(integrat_over_r_precal), intent(in) :: container
    logical :: boolean
  
    if (allocated(container%val) == .false.) then
      boolean = .false.
    else
      boolean = .true.
    end if  
    
  end function is_integrat_over_r_precal_allocated



  function make_and_cal_integrat_over_r_precal(n_r, bin_length, q) result(container)
    real(db), dimension(:), intent(in) :: q
    real(db), intent(in) :: bin_length
    integer, intent(in) :: n_r
    type (integrat_over_r_precal) :: container
    
    integer :: i_q, i_r, n_q
    real(db) :: lower_lim, upper_lim, prefac, r
    
    n_q = size(q)
    
    allocate(container%q(n_q))
    
    container%q = q
    
    allocate(container%val(n_r, n_q)) 
    
    do i_q = 1, n_q
      prefac = 4 * pi_value / q(i_q)**2
      
      lower_lim = 0  ! since sin(Q*r) = 0 for r=0
      
      do i_r = 1, n_r
        r = i_r * bin_length
        
        upper_lim = sin(q(i_q)*r) / q(i_q) - r * cos(q(i_q)*r)
        
        container%val(i_r,i_q) = (upper_lim - lower_lim) * prefac
        
        lower_lim = upper_lim  
      end do
    end do
    
  end function make_and_cal_integrat_over_r_precal



  ! delta_time and n_time are needed to pass on the time-binning info to this container

  function make_s_q_time(q, delta_time, n_time) result(container)
    real(db), dimension(:), intent(in) :: q
    real(db), intent(in) :: delta_time
    integer, intent(in) :: n_time
    type (s_q_time) :: container
    
    integer :: n_q
    
    n_q = size(q)
    
    allocate(container%q(n_q))
    
    container%q = q
    container%delta_time = delta_time
    
    allocate(container%self(n_q,n_time)) 
    allocate(container%diff(n_q,n_time))
    
    container%self = 0
    container%diff = 0
    
  end function make_s_q_time




end module s_q_time_class