module rdf_class
use various_constants_class
use structure_class
use histogram_class
use tic_toc_class

implicit none

  public :: make_rdf, cal_rdf, save_rdf

  ! it is assumed that histogram stores distances in bins of size
  ! bin_length; with the first bin representing distances between
  ! zero and bin_length and the last bin representing distances
  ! between (n-bins-1)*bin_length and n_bins*bin_length
  
  type rdf
    ! bin_length is e.g. used in cal_rdf
    ! to make sure that the histogram passed to that function is 
    ! compatable with how g(r) is defined here
    
    real(db) :: bin_length = 0.0
    
    
    ! Notice the val values are always to be the g(r) values at the
    ! center of each bin. I.e. at 0.5*bin_length, 1.5*bin_length etc.... 
    
    real(db), dimension(:), allocatable :: val
    real(db), dimension(:), allocatable :: prefac   ! g(r)=histogram*prefac
  end type rdf
  
  character(len=10), parameter, private :: rdf_filename_prefix = "output/rdf"  
  integer, private :: filename_number = 1  ! so first saved rdf file will be called rdf1.xml

contains

  function make_rdf(volume, n_atoms, n_bins, bin_length) result(a_rdf)
    real(db), intent(in) :: volume, bin_length
    integer, intent(in) :: n_atoms, n_bins
    type (rdf) :: a_rdf
    
    real(db) :: temp, r_i
    integer :: i, n_bin
    
    a_rdf%bin_length = bin_length
 
    ! make val attribute
    
    allocate(a_rdf%val(n_bins))
    a_rdf%val = 0.0
 
    ! make prefac attribute
    
    allocate(a_rdf%prefac(n_bins))
    
    ! On pages 19a and 19aa of my handwritten notes I consider different 
    ! possible choices for defining g(r). Here uses
    ! 
    !   g(r_i) = V * hist(r_i) / (n_atoms^2*4*pi*(r_i)^2*dr
    
    ! Calculate prefac = V / (n_atoms^2*4*pi*(r_i)^2*dr
    ! Note 4*pi*(r_i)^2*dr = 4*pi*(i-1/2)^2*dr^3 
    ! with width dr and center position r_i=dr*(i-1/2).

    temp = volume / (4.0 * pi_value * (n_atoms**2) * bin_length**3 )
    
    do i = 1, n_bins
      r_i = i - 0.5_db
      a_rdf%prefac(i) = temp / (r_i**2)
    end do
  
  end function make_rdf
  

  ! Convert a histogram or accumulated histogram into an rdf by applying the 
  ! prefactors calculated in make_rdf. Here assumed output rdf and input hist 
  ! have the same pre-allocated array sizes
  !
  subroutine cal_rdf(a_rdf, hist)
    type (rdf), intent(inout) :: a_rdf
    type (histogram), intent(in) :: hist
    
    ! TODO: check that histogram has the same size as that of the rdf
    
    if (hist%n_accum > 0) then
      a_rdf%val = a_rdf%prefac * hist%sum / hist%n_accum
    else
      a_rdf%val = a_rdf%prefac * hist%val
    end if

  end subroutine cal_rdf


  ! Notice the a_rdf%val array are assumed to represent g(r) at the center
  ! of each bin. That is at: 0.5*bin_length, 1.5*bin_length etc......
  ! The temperature is required to passed in dimensional units here.
  !
  subroutine save_rdf(a_rdf, temperature, density)
    use flib_wxml
    type (rdf), intent(in) :: a_rdf
    real(db), optional, intent(in) :: temperature, density

    type (xmlf_t) :: xf
    integer :: i, n_bin
    character(len=50) :: filename
    
    if (filename_number < 10) then
      write(filename, '(i1)') filename_number
    else if (filename_number < 100) then
      write(filename, '(i2)') filename_number
    else if (filename_number < 1000) then
      write(filename, '(i3)') filename_number
    else
      write(*,*) "ERROR: in save_rdf"
      write(*,*) "It is assumed that you did not intend to write"
      write(*,*) "to disk 1000 rdf xml files!!!!"
      stop
    end if
    
    filename_number = filename_number + 1
    
    filename = rdf_filename_prefix // trim(filename) // ".xml"
    
    write(*,'(3a)') "Write ", trim(filename), " to disk"
    
    n_bin = size(a_rdf%val)    
    
    call xml_OpenFile(filename, xf, indent=.true.)
    
    call xml_AddXMLDeclaration(xf, "UTF-8")
    call xml_NewElement(xf, "rdf")
    
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
    
    call xml_AddAttribute(xf, "units", "AA")
    call xml_AddAttribute(xf, "bin-length", str(a_rdf%bin_length, format="(f10.5)"))
    call xml_NewElement(xf, "this-file-was-created")
    call xml_AddAttribute(xf, "when", get_current_date_and_time())
    call xml_EndElement(xf, "this-file-was-created")
    
    do i = 1, n_bin
      call xml_NewElement(xf, "g-of-r")
      call xml_AddAttribute(xf, "r", str((i-0.5)*a_rdf%bin_length, format="(f10.5)"))
      call xml_AddAttribute(xf, "g", str(a_rdf%val(i), format="(f10.5)"))
      call xml_EndElement(xf, "g-of-r")
    end do 
    
    call xml_EndElement(xf, "rdf")
    call xml_Close(xf)
    
  end subroutine save_rdf

end module rdf_class