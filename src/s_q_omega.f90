module s_q_omega_class
use various_constants_class         
use structure_type_definition_class 
use structure_nn_methods_class
use tic_toc_class
use time_corr_hist_container_class
use s_q_time_class

implicit none

  public :: cal_s_q_omega
  public :: make_s_q_omega
  public :: check_if_s_q_omega_allocated
  public :: print_s_q_omega
  
  private :: make_and_cal_precal_cos_values
  private :: is_precal_cos_values_allocated
  
  type s_q_omega
    real(db), dimension(:), allocatable :: q
    real(db), dimension(:), allocatable :: omega
                            
    real(db), dimension(:,:), allocatable :: self   ! dimensions (n_q, n_omega)
    real(db), dimension(:,:), allocatable :: diff   ! dimensions (n_q, n_omega)
  end type s_q_omega

  type precal_cos_values
    real(db), dimension(:), allocatable :: omega  ! used perhaps in the future to figure out
                                                  ! if this container needs recalculating
                                             
    real(db), dimension(:,:), allocatable :: val  ! dimension (n_t, n_omega)
  end type precal_cos_values
  
  private precal_cos_values
  
  
  type (precal_cos_values), private :: local_precal_cos_values  
  
  
  character(len=16), parameter, private :: filename_prefix = "output/s_q_omega"    
  integer, private :: filname_number = 1    
  
contains

  ! Calculate dynamical structure factors from intermediate scattering function as
  !
  ! S(q,omega) = (1/pi) * int_0^infinity ( cos(omega*t) * S(q,t)) dt
  !
  ! where this integral is approximated according to method make_and_cal_precal_cos_values()
  !
  subroutine cal_s_q_omega(s_qt, str, s_qo)
    type(s_q_time), intent(in) :: s_qt
    type(structure), intent(in) :: str
    type(s_q_omega), intent(inout) :: s_qo
  
    integer :: n_q, n_omega, n_t, i_q, i_omega, n_atom
    real(db) :: density, pre_s_d, pre_s_s
    real(db), dimension(:), allocatable :: prefac_s_d, prefac_s_s
  
    ! some checks
       
    call check_if_s_q_time_allocated(s_qt)
    call check_if_s_q_omega_allocated(s_qo)
    
    
    ! Has local_precal_cos_values been allocated and populated
    ! Notice should also compare the omega array in that container with 
    ! omega array in the s_q_omega container  
    
    n_t = get_s_q_time_n_t_bin(s_qt)
    
    if (is_precal_cos_values_allocated(local_precal_cos_values) == .false.) then
      local_precal_cos_values = make_and_cal_precal_cos_values( &
        n_t, s_qt%t_bin, s_qo%omega)
    end if
  

    ! calculate S(Q,omega)
  
    s_qo%self = matmul(s_qt%self, local_precal_cos_values%val)
    s_qo%diff = matmul(s_qt%diff, local_precal_cos_values%val) 

  end subroutine cal_s_q_omega


  ! The temperature is assumed to be in dimensionless units.
  ! Print S(q, omega) to either filename given by name_of_file or
  ! if this argument is not specified the name given by module 
  ! attribute filename_prefix + number
  !
  subroutine print_s_q_omega(container, density, temperature, name_of_file)
    use flib_wxml
    type(s_q_omega), intent(in) :: container  
    real(db), intent(in) :: density    
    real(db), intent(in) :: temperature
    character(len=*), optional, intent(in) :: name_of_file
    
    type (xmlf_t) :: xf
    integer :: i_q, i_o, n_omega
    
    character(len=50) :: filename
    
    if ( present(name_of_file) ) then
      filename = name_of_file
    else
      if (filname_number < 10) then
        write(filename, '(i1)') filname_number
      else if (filname_number < 100) then
        write(filename, '(i2)') filname_number
      else if (filname_number < 1000) then
        write(filename, '(i3)') filname_number
      else
        write(*,*) "ERROR: in print_s_q_omega"
        write(*,*) "It is assumed that you did not intend to write"
        write(*,*) "to disk 1000 s(q,omega) xml files!"
        stop
      end if
    
      filname_number = filname_number + 1
      filename = filename_prefix // trim(filename) // ".xml" ! here filename is just a number on rhs  
    end if
    
    write(*,'(3a)') "Write ", trim(filename), " to disk"
    
    
    call xml_OpenFile(filename, xf, indent=.true.)
    
    call xml_AddXMLDeclaration(xf, "UTF-8")
    call xml_NewElement(xf, "s-q-omega")
    
    ! notice convert units of temperature from dimensionless to K  
    call xml_AddAttribute(xf, "title", "T = " // trim(str(temperature * T_unit, format="(f10.5)")) // &
                                         " K: rho = " // trim(str(density, format="(f10.5)")) &
                                         // " atoms/AA-3")

    call xml_AddAttribute(xf, "n-omega-points", str(size(container%omega), format="(i)") )
    call xml_AddAttribute(xf, "n-q-points", str(size(container%q), format="(i)") )
    call xml_AddAttribute(xf, "q-units", "AA^-1")
    call xml_AddAttribute(xf, "omega-unit", "1/[10^-13 s]")        
    
    call xml_NewElement(xf, "this-file-was-created")
    call xml_AddAttribute(xf, "when", get_current_date_and_time())
    call xml_EndElement(xf, "this-file-was-created")  
    
    
    n_omega = size(container%self, 2)

    do i_o = 1, n_omega   
      do i_q = 1, size(container%q)
        call xml_NewElement(xf, "SQomega")
        call xml_AddAttribute(xf, "q", str(container%q(i_q), format="(f10.5)"))
        call xml_AddAttribute(xf, "omega", str(container%omega(i_o), format="(f10.5)"))
        call xml_AddAttribute(xf, "S-self", str(container%self(i_q, i_o), format="(f15.5)"))
        call xml_AddAttribute(xf, "S-diff", str(container%diff(i_q, i_o), format="(f15.5)"))        
        call xml_EndElement(xf, "SQomega")
      end do 
    end do
    
    call xml_EndElement(xf, "s-q-omega")
    
    call xml_Close(xf)    
  
  end subroutine print_s_q_omega


  subroutine check_if_s_q_omega_allocated(container)
    type(s_q_omega), intent(in) :: container
  
    if (allocated(container%q) == .false.) then
      write(*,*) " "
      write(*,*) "ERROR in s_q_omega.f90"
      write(*,*) "Forgot to allocate s_q_omega"
      stop
    end if  
  end subroutine check_if_s_q_omega_allocated


  function is_precal_cos_values_allocated(container) result(boolean)
    type(precal_cos_values), intent(in) :: container
    logical :: boolean
  
    if (allocated(container%val) == .false.) then
      boolean = .false.
    else
      boolean = .true.
    end if  
    
  end function is_precal_cos_values_allocated


  ! Calculate factors to convert S(q,t) to S(q,omega) that is:
  ! S(q, omega) = precal(t_i)*S(q,t_i))
  ! where this is an approximation to the fourier integral over time to get to S(q,omega)
  ! This approximation is here done by using method 2 page 31 in my handwritten notes.
  ! 
  function make_and_cal_precal_cos_values(n_t, time_length, omega) result(container)
    real(db), dimension(:), intent(in) :: omega
    real(db), intent(in) :: time_length
    integer, intent(in) :: n_t
    type (precal_cos_values) :: container
    
    integer :: i_omega, i_t, n_omega
    real(db) :: t, prefac, lower_lim, upper_lim
    
    n_omega = size(omega)
    
    allocate(container%omega(n_omega))
    
    container%omega = omega
    
    allocate(container%val(n_t, n_omega)) 
    
    
    do i_omega = 1, n_omega
      ! when calculating integral this way (method 2 in my notes page 31)
      ! I need to calculate omega=0 separately
      
      if ( omega(i_omega) == 0.0 ) then
        lower_lim = 0  ! since first bin is assumed only to be a half delta-time
        do i_t = 0, n_t-1
          upper_lim = (i_t + 0.5) * time_length  ! all subsequent bins are a full delta-time
          
          container%val(i_t+1,i_omega) = (upper_lim - lower_lim) / pi_value
          lower_lim = upper_lim 
        end do
      else
        prefac = 1 / ( pi_value * omega(i_omega) )
      
        lower_lim = 0  ! since sin(omega*t) = 0 for t=0
      
        do i_t = 0, n_t-1
          t = (i_t + 0.5) * time_length
        
          upper_lim = sin(omega(i_omega)*t)
        
          container%val(i_t+1,i_omega) = (upper_lim - lower_lim) * prefac
          lower_lim = upper_lim  
        end do
      end if
    end do
    
  end function make_and_cal_precal_cos_values



  ! t_bin and n_time are needed to pass on the time-binning info to this container

  function make_s_q_omega(q, omega) result(container)
    real(db), dimension(:), intent(in) :: q
    real(db), dimension(:), intent(in) :: omega
    type (s_q_omega) :: container
    
    integer :: n_q
    integer :: n_omega
    
    n_q = size(q)
    n_omega = size(omega)
    
    allocate(container%q(n_q))
    allocate(container%omega(n_omega))
    
    container%q = q
    container%omega = omega
    
    allocate(container%self(n_q,n_omega)) 
    allocate(container%diff(n_q,n_omega))
    
    container%self = 0
    container%diff = 0
    
  end function make_s_q_omega

end module s_q_omega_class