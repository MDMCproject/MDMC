module time_corr_hist_container_class
use various_constants_class
use structure_class
use histogram_cutdown_class

implicit none


  public :: clear_time_corr_hist_container
  public :: make_time_corr_hist_container 


  type time_corr_hist_container
    integer :: n_accum

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

  subroutine clear_time_corr_hist_container(container)
    type(time_corr_hist_container), intent(inout) :: container
    
    integer :: i
    
    if (allocated(container%einstein_diffuse_exp) == .false.) then
      write(*,*) " "
      write(*,*) "ERROR in time_corr_hist_container.f90"
      write(*,*) "Forgot to allocate time_corr_hist_container"
      stop
    end if
     
    container%n_accum = 0
    container%einstein_diffuse_exp = 0.0
    
    do i = 1, size(container%g_s_hists_sum)
      container%g_s_hists_sum(i)%val = 0
      container%g_d_hists_sum(i)%val = 0
    end do 
    
  end subroutine clear_time_corr_hist_container
  

  function make_time_corr_hist_container(r_max, bin_length, n_time_evals, delta_time) result(container)
    real(db), intent(in) :: r_max, bin_length, delta_time
    integer, intent(in) :: n_time_evals
    
    type (time_corr_hist_container) :: container
    integer :: n_bin, i
    real(db) :: temp, r_i
    
    container%n_accum = 0
    
    allocate(container%einstein_diffuse_exp(n_time_evals))
    
    allocate(container%g_s_hists_sum(n_time_evals))
    allocate(container%g_d_hists_sum(n_time_evals))
    
    do i = 1, n_time_evals
      container%g_s_hists_sum(i) = make_histogram_cutdown(r_max, bin_length)
      container%g_d_hists_sum(i) = make_histogram_cutdown(r_max, bin_length)
    end do
    
    
    n_bin = size(container%g_s_hists_sum(1)%val)
    
    
    allocate(container%volume_prefac(n_bin))
    
    temp = 1.0 / ( 4.0 * pi_value * bin_length**3 )
    
    do i = 1, n_bin
      r_i = i - 0.5
      container%volume_prefac(i) = temp / (r_i**2)
    end do    
    
  end function make_time_corr_hist_container


 ! The temperature is required to passed in dimensional units here.

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

end module time_corr_hist_container_class