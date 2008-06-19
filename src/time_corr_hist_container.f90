module time_corr_hist_container_class
use various_constants_class
use structure_class
use histogram_cutdown_class
use tic_toc_class

implicit none


  public :: clear_time_corr_hist_container
  public :: make_time_corr_hist_container 

  public :: print_g_d, print_g_s, print_einstein_diffuse_exp
  
  public :: check_if_time_corr_hist_container_is_allocated
  
  public :: get_time_corr_hist_r_max, get_time_corr_hist_bin_length, &
            get_time_corr_hist_n_time_evals, get_time_corr_hist_n_bin
            

  ! for debugging
  
  public :: print_h_s_hist, print_h_d_hist


  type time_corr_hist_container
    integer :: n_accum

    real(db) :: time_step  ! otherwise time-binning not fully defined.
                           ! time between different cals of g_s and g_d etc.
                           !(r-binning stuff is in histogram_cutdown) 

    ! dimensions are n_bin x n_time
    type (histogram_cutdown), dimension(:), allocatable :: g_s_hists_sum
    type (histogram_cutdown), dimension(:), allocatable :: g_d_hists_sum    
    
    ! dimension is n_time   
    real(db), dimension(:), allocatable :: einstein_diffuse_exp
  
    ! dimension is n_bin. Store 1/(4*pi*r^2*dr).
    real(db), dimension(:), allocatable :: volume_prefac 
    
  end type time_corr_hist_container
  
  
  character(len=27), parameter, private :: filename_prefix1 = "output/einstein_diffuse_exp"  
  integer, private :: filname_number1 = 1  ! so first saved rdf file will be called einstein_diffuse_exp1.xml
  
  character(len=10), parameter, private :: filename_prefix2 = "output/g_s"  
  integer, private :: filname_number2 = 1  
  
  character(len=10), parameter, private :: filename_prefix3 = "output/g_d"  
  integer, private :: filname_number3 = 1      

contains

  subroutine check_if_time_corr_hist_container_is_allocated(container)
    type(time_corr_hist_container), intent(in) :: container
  
    if (allocated(container%einstein_diffuse_exp) == .false.) then
      write(*,*) " "
      write(*,*) "ERROR in time_corr_hist_container.f90"
      write(*,*) "Forgot to allocate time_corr_hist_container"
      stop
    end if  
  end subroutine check_if_time_corr_hist_container_is_allocated




  subroutine clear_time_corr_hist_container(container)
    type(time_corr_hist_container), intent(inout) :: container
    
    integer :: i
    
    call check_if_time_corr_hist_container_is_allocated(container)
     
    container%n_accum = 0
    container%einstein_diffuse_exp = 0.0
    
    do i = 1, size(container%g_s_hists_sum)
      container%g_s_hists_sum(i)%val = 0
      container%g_d_hists_sum(i)%val = 0
    end do 
    
  end subroutine clear_time_corr_hist_container
  
  
  ! notice delta_time here is the time between to time calculations of g_d, g_s etc.

  function make_time_corr_hist_container(r_max, bin_length, n_time_evals, delta_time) result(container)
    real(db), intent(in) :: r_max, bin_length, delta_time
    integer, intent(in) :: n_time_evals
    
    type (time_corr_hist_container) :: container
    integer :: n_bin, i
    real(db) :: temp, r_i
    
    container%n_accum = 0
    
    container%time_step = delta_time 
    
    allocate(container%einstein_diffuse_exp(n_time_evals))
    
    allocate(container%g_s_hists_sum(n_time_evals))
    allocate(container%g_d_hists_sum(n_time_evals))
    
    do i = 1, n_time_evals
      container%g_s_hists_sum(i) = make_histogram_cutdown(r_max, bin_length)
      container%g_d_hists_sum(i) = make_histogram_cutdown(r_max, bin_length)
    end do
    
    
    n_bin = size(container%g_s_hists_sum(1)%val)
    
    
    allocate(container%volume_prefac(n_bin))
    
    
    ! cal volume element as (4*pi/3)*delta_r^3*(i^3 - (i-1)^3), where 
    ! i = 1, 2, ... n_bin
    
    temp = 3.0 / ( 4.0 * pi_value * bin_length**3 )
    
    do i = 1, n_bin
      !r_i = i - 0.5
      container%volume_prefac(i) = temp / ( i**3-(i-1)**3 )
    end do    
    
  end function make_time_corr_hist_container


  function get_time_corr_hist_n_time_evals(container) result(n_time_evals)
    type(time_corr_hist_container), intent(in) :: container
    integer :: n_time_evals
      
     call check_if_time_corr_hist_container_is_allocated(container) 
     
     n_time_evals = size(container%einstein_diffuse_exp)
  end function get_time_corr_hist_n_time_evals
  
  
  function get_time_corr_hist_bin_length(container) result(bin_length)
    type(time_corr_hist_container), intent(in) :: container
    real(db) :: bin_length  
    
    call check_if_time_corr_hist_container_is_allocated(container)     
    
    bin_length = container%g_s_hists_sum(1)%bin_length
  
  end function get_time_corr_hist_bin_length


  function get_time_corr_hist_r_max(container) result(r_max)
    type(time_corr_hist_container), intent(in) :: container
    real(db) :: r_max  
    
    integer n_bin
    real(db) :: bin_length
    
    call check_if_time_corr_hist_container_is_allocated(container)     
    
    bin_length = container%g_s_hists_sum(1)%bin_length
    n_bin = size(container%g_s_hists_sum(1)%val)
    
    r_max = bin_length * n_bin
  
  end function get_time_corr_hist_r_max
  
  
  function get_time_corr_hist_n_bin(container) result(n_bin)
    type(time_corr_hist_container), intent(in) :: container
    integer :: n_bin  
    
    call check_if_time_corr_hist_container_is_allocated(container)     

    n_bin = size(container%g_s_hists_sum(1)%val)
  
  end function get_time_corr_hist_n_bin  


 ! The temperature is assumed to be in dimensionless units.

  subroutine print_einstein_diffuse_exp(container, density, temperature)
    use flib_wxml
    type(time_corr_hist_container), intent(in) :: container
    real(db), optional, intent(in) :: temperature, density
    
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
    
    n_eval_times = size(container%einstein_diffuse_exp)    
    
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
      call xml_AddAttribute(xf, "t", str((i-1)*container%time_step, format="(f10.5)"))
      call xml_AddAttribute(xf, "D", str(container%einstein_diffuse_exp(i), format="(f10.5)"))
      call xml_EndElement(xf, "D-of-t")
    end do 
    
    call xml_EndElement(xf, "einstein-diffuse-exp")
    
    call xml_Close(xf)    
  
  end subroutine print_einstein_diffuse_exp


  ! The temperature is assumed to be in dimensionless units.

  subroutine print_g_s(container, density, temperature)
    use flib_wxml
    type(time_corr_hist_container), intent(inout) :: container  ! inout because volume_prefac modified
    real(db), intent(in) :: density    
    real(db), optional, intent(in) :: temperature
    
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
      write(*,*) "ERROR: in print_g_s"
      write(*,*) "It is assumed that you did not intend to write"
      write(*,*) "to disk 1000 g_s xml files!!!!"
      stop
    end if
    
    filname_number2 = filname_number2 + 1
    
    filename = filename_prefix2 // trim(filename) // ".xml"  ! here filename is just a number on rhs
    
    write(*,'(3a)') "Write ", trim(filename), " to disk"
    
    
    n_eval_times = size(container%einstein_diffuse_exp)    
    n_bin = size(container%g_s_hists_sum(1)%val)
    bin_length = container%g_s_hists_sum(1)%bin_length
    
    call xml_OpenFile(filename, xf, indent=.true.)
    
    call xml_AddXMLDeclaration(xf, "UTF-8")
    call xml_NewElement(xf, "G_s-space-time-pair-correlation-function")
    
    ! notice convert units of temperature from dimensionless to K  
    if (present(temperature)) then
      call xml_AddAttribute(xf, "title", "T = " // trim(str(temperature * T_unit, format="(f10.5)")) // &
                                         " K: rho = " // trim(str(density, format="(f10.5)")) &
                                         // " atoms/AA-3")
    else 
      call xml_AddAttribute(xf, "title", "rho = " // str(density, format="(f10.5)") &
                                         // "atoms/AA-3")
    end if
    
    call xml_AddAttribute(xf, "bin-length", str(bin_length, format="(f10.5)"))
    
    call xml_AddAttribute(xf, "time-unit", "10^-13 s")
    call xml_AddAttribute(xf, "r-units", "AA")
    
    
    call xml_NewElement(xf, "this-file-was-created")
    call xml_AddAttribute(xf, "when", get_current_date_and_time())
    call xml_EndElement(xf, "this-file-was-created")
    
    
    container%volume_prefac = container%volume_prefac / (container%n_accum*density)

    do i = 1, n_eval_times    
      do i_bin = 1, n_bin
        call xml_NewElement(xf, "G-s")
        call xml_AddAttribute(xf, "r", str((i_bin-0.5)*bin_length, format="(f10.5)"))
        call xml_AddAttribute(xf, "t", str((i-1)*container%time_step, format="(f10.5)"))
        call xml_AddAttribute(xf, "G", str(container%g_s_hists_sum(i)%val(i_bin) * &
                              container%volume_prefac(i_bin), format="(f15.5)"))
        call xml_EndElement(xf, "G-s")
      end do 
    end do
    
    container%volume_prefac = container%volume_prefac * (container%n_accum*density)
    
    call xml_EndElement(xf, "G_s-space-time-pair-correlation-function")
    
    call xml_Close(xf)    
  
  end subroutine print_g_s


  ! The temperature is assumed to be in dimensionless units.

  subroutine print_g_d(container, volume, n_atom, temperature)
    use flib_wxml
    type(time_corr_hist_container), intent(inout) :: container  ! inout because volume_prefac modified  
    real(db), intent(in) :: volume 
    real(db), optional, intent(in) :: temperature
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
      write(*,*) "ERROR: in print_g_d"
      write(*,*) "It is assumed that you did not intend to write"
      write(*,*) "to disk 1000 g_d xml files!!!!"
      stop
    end if
    
    filname_number3 = filname_number3 + 1
    
    filename = filename_prefix3 // trim(filename) // ".xml"
    
    write(*,'(3a)') "Write ", trim(filename), " to disk"
    
    
    n_eval_times = size(container%einstein_diffuse_exp)    
    n_bin = size(container%g_s_hists_sum(1)%val)
    bin_length = container%g_s_hists_sum(1)%bin_length
    
    call xml_OpenFile(filename, xf, indent=.true.)
    
    call xml_AddXMLDeclaration(xf, "UTF-8")
    call xml_NewElement(xf, "G_d-space-time-pair-correlation-function")
    
    ! notice convert units of temperature from dimensionless to K  
    if (present(temperature)) then
      call xml_AddAttribute(xf, "title", "T = " // trim(str(temperature * T_unit, format="(f10.5)")) // &
                                         " K: rho = " // trim(str(density, format="(f10.5)")) &
                                         // " atoms/AA-3")
    else 
      call xml_AddAttribute(xf, "title", "rho = " // str(density, format="(f10.5)") &
                                         // "atoms/AA-3")
    end if
    
    call xml_AddAttribute(xf, "bin-length", str(bin_length, format="(f10.5)"))
    
    call xml_AddAttribute(xf, "time-unit", "10^-13 s")
    call xml_AddAttribute(xf, "r-units", "AA")
    
    
    call xml_NewElement(xf, "this-file-was-created")
    call xml_AddAttribute(xf, "when", get_current_date_and_time())
    call xml_EndElement(xf, "this-file-was-created")
    
    
    container%volume_prefac = container%volume_prefac / (container%n_accum*density*(n_atom-1))

    do i = 1, n_eval_times    
      do i_bin = 1, n_bin
        call xml_NewElement(xf, "G-d")
        call xml_AddAttribute(xf, "r", str((i_bin-0.5)*bin_length, format="(f10.5)"))
        call xml_AddAttribute(xf, "t", str((i-1)*container%time_step, format="(f10.5)"))
        call xml_AddAttribute(xf, "G", str(container%g_d_hists_sum(i)%val(i_bin) * &
                                           container%volume_prefac(i_bin), format="(f15.5)"))
        call xml_EndElement(xf, "G-d")
      end do 
    end do
    
    container%volume_prefac = container%volume_prefac * (container%n_accum*density*(n_atom-1))
    
    call xml_EndElement(xf, "G_d-space-time-pair-correlation-function")
    
    call xml_Close(xf)    
  
  end subroutine print_g_d  


  ! 1st version of this function was written for debugging purposes

  subroutine print_h_s_hist(container, density, temperature)
    use flib_wxml
    type(time_corr_hist_container), intent(inout) :: container  ! inout because volume_prefac modified
    real(db), intent(in) :: density    
    real(db), optional, intent(in) :: temperature
    
    type (xmlf_t) :: xf
    integer :: i, n_eval_times, n_bin, i_bin
    real(db) :: bin_length
    character(len=50) :: filename
    
!    if (filname_number2 < 10) then
!      write(filename, '(i1)') filname_number2
!    else if (filname_number2 < 100) then
!      write(filename, '(i2)') filname_number2
!    else if (filname_number2 < 1000) then
!      write(filename, '(i3)') filname_number2
!    else
!      write(*,*) "ERROR: in save_rdf"
!      write(*,*) "It is assumed that you did not intend to write"
!      write(*,*) "to disk 1000 rdf xml files!!!!"
!      stop
!    end if
!    
!    filname_number2 = filname_number2 + 1
!    
!    filename = filename_prefix2 // trim(filename) // ".xml"  ! here filename is just a number on rhs
    
    filename = "output/h_s_histogram.xml"
    
    write(*,'(3a)') "Write ", trim(filename), " to disk"
    
    
    n_eval_times = size(container%einstein_diffuse_exp)    
    n_bin = size(container%g_s_hists_sum(1)%val)
    bin_length = container%g_s_hists_sum(1)%bin_length
    
    call xml_OpenFile(filename, xf, indent=.true.)
    
    call xml_AddXMLDeclaration(xf, "UTF-8")
    call xml_NewElement(xf, "normalised-h-s-histogram")
    
    ! notice convert units of temperature from dimensionless to K  
    if (present(temperature)) then
      call xml_AddAttribute(xf, "title", "T = " // trim(str(temperature * T_unit, format="(f10.5)")) // &
                                         " K: rho = " // trim(str(density, format="(f10.5)")) &
                                         // " atoms/AA-3")
    else 
      call xml_AddAttribute(xf, "title", "rho = " // str(density, format="(f10.5)") &
                                         // "atoms/AA-3")
    end if
    
    call xml_AddAttribute(xf, "bin-length", str(bin_length, format="(f10.5)"))
    
    call xml_AddAttribute(xf, "time-unit", "10^-13 s")
    call xml_AddAttribute(xf, "r-units", "AA")
    
    
    call xml_NewElement(xf, "this-file-was-created")
    call xml_AddAttribute(xf, "when", get_current_date_and_time())
    call xml_EndElement(xf, "this-file-was-created")
    

    do i = 1, n_eval_times    
      do i_bin = 1, n_bin
        call xml_NewElement(xf, "h-s")
        call xml_AddAttribute(xf, "r", str((i_bin-0.5)*bin_length, format="(f10.5)"))
        call xml_AddAttribute(xf, "t", str((i-1)*container%time_step, format="(f10.5)"))
        call xml_AddAttribute(xf, "h", str(container%g_s_hists_sum(i)%val(i_bin) / &
                              container%n_accum, format="(i10)"))
        call xml_EndElement(xf, "h-s")
      end do 
    end do
    
    call xml_EndElement(xf, "normalised-h-s-histogram")
    
    call xml_Close(xf)    
  
  end subroutine print_h_s_hist


  ! 1st version of this function was written for debugging purposes

  subroutine print_h_d_hist(container, density, temperature)
    use flib_wxml
    type(time_corr_hist_container), intent(inout) :: container  ! inout because volume_prefac modified
    real(db), intent(in) :: density    
    real(db), optional, intent(in) :: temperature
    
    type (xmlf_t) :: xf
    integer :: i, n_eval_times, n_bin, i_bin
    real(db) :: bin_length
    character(len=50) :: filename
    
!    if (filname_number2 < 10) then
!      write(filename, '(i1)') filname_number2
!    else if (filname_number2 < 100) then
!      write(filename, '(i2)') filname_number2
!    else if (filname_number2 < 1000) then
!      write(filename, '(i3)') filname_number2
!    else
!      write(*,*) "ERROR: in save_rdf"
!      write(*,*) "It is assumed that you did not intend to write"
!      write(*,*) "to disk 1000 rdf xml files!!!!"
!      stop
!    end if
!    
!    filname_number2 = filname_number2 + 1
!    
!    filename = filename_prefix2 // trim(filename) // ".xml"  ! here filename is just a number on rhs
    
    filename = "output/h_d_histogram.xml"
    
    write(*,'(3a)') "Write ", trim(filename), " to disk"
    
    
    n_eval_times = size(container%einstein_diffuse_exp)    
    n_bin = size(container%g_d_hists_sum(1)%val)
    bin_length = container%g_d_hists_sum(1)%bin_length
    
    call xml_OpenFile(filename, xf, indent=.true.)
    
    call xml_AddXMLDeclaration(xf, "UTF-8")
    call xml_NewElement(xf, "normalised-h-d-histogram")
    
    ! notice convert units of temperature from dimensionless to K  
    if (present(temperature)) then
      call xml_AddAttribute(xf, "title", "T = " // trim(str(temperature * T_unit, format="(f10.5)")) // &
                                         " K: rho = " // trim(str(density, format="(f10.5)")) &
                                         // " atoms/AA-3")
    else 
      call xml_AddAttribute(xf, "title", "rho = " // str(density, format="(f10.5)") &
                                         // "atoms/AA-3")
    end if
    
    call xml_AddAttribute(xf, "bin-length", str(bin_length, format="(f10.5)"))
    
    call xml_AddAttribute(xf, "time-unit", "10^-13 s")
    call xml_AddAttribute(xf, "r-units", "AA")
    
    
    call xml_NewElement(xf, "this-file-was-created")
    call xml_AddAttribute(xf, "when", get_current_date_and_time())
    call xml_EndElement(xf, "this-file-was-created")
    

    do i = 1, n_eval_times    
      do i_bin = 1, n_bin
        call xml_NewElement(xf, "h-d")
        call xml_AddAttribute(xf, "r", str((i_bin-0.5)*bin_length, format="(f10.5)"))
        call xml_AddAttribute(xf, "t", str((i-1)*container%time_step, format="(f10.5)"))
        call xml_AddAttribute(xf, "h", str(container%g_d_hists_sum(i)%val(i_bin) / &
                              container%n_accum, format="(i10)"))
        call xml_EndElement(xf, "h-d")
      end do 
    end do
    
    call xml_EndElement(xf, "normalised-h-d-histogram")
    
    call xml_Close(xf)    
  
  end subroutine print_h_d_hist

end module time_corr_hist_container_class