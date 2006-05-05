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
    ! Notice r_max = size(h) * bin_length. Hence strictly only one
    ! of bin_length and r_max are needed, but both are included for
    ! convenience
    ! Also, notice the bin_length and r_max are e.g. used in cal_rdf
    ! to make sure that the histogram passed to that function is 
    ! compatable with how g(r) is defined in this type
    
    real(db) :: bin_length  
    real(db) :: r_max  
    
    real(db), dimension(:), allocatable :: val
    real(db), dimension(:), allocatable :: prefac   ! g(r)=histogram*prefac
  end type rdf
  
  
  character(len=10), parameter, private :: rdf_filename_prefix = "output/rdf"  
  integer, private :: filname_number = 1  ! so first saved rdf file will be called rdf1.xml

contains

  function make_rdf(volume, n_atoms, r_max, n_bin) result(a_rdf)
    real(db), intent(in) :: volume, r_max
    integer, intent(in) :: n_atoms, n_bin
    type (rdf) :: a_rdf
    
    real(db) :: temp, r_i
    integer :: i
    
    a_rdf%r_max = r_max
    a_rdf%bin_length = r_max / n_bin
 
 
    ! make val attribute
    
    allocate(a_rdf%val(n_bin))
    a_rdf%val = 0.0
 
 
    ! make prefac attribute
    
    allocate(a_rdf%prefac(n_bin))
    

    temp = volume / (2.0 * pi_value * (n_atoms**2) * a_rdf%bin_length**3 )
    
    do i = 1, n_bin
      r_i = i - 0.5_db
      a_rdf%prefac(i) = temp / (r_i**2)
    end do
  
  end function make_rdf
  
  
  
  subroutine cal_rdf(a_rdf, hist)
    type (rdf), intent(inout) :: a_rdf
    type (histogram), intent(in) :: hist
    
    ! check that histogram has the same format as that required to calculate
    ! the rdf
    
    ! PUT CODE IN HERE LATER 
    
    if (hist%n_accum > 0) then
      a_rdf%val = a_rdf%prefac * hist%sum / hist%n_accum
    else
      a_rdf%val = a_rdf%prefac * hist%val
    end if

  end subroutine cal_rdf


  subroutine save_rdf(a_rdf, temperature, density)
    use flib_wxml
    type (rdf), intent(in) :: a_rdf
    real(db), optional, intent(in) :: temperature, density

    type (xmlf_t) :: xf
    integer :: i, n_bin
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
    
    filename = rdf_filename_prefix // trim(filename) // ".xml"
    
    write(*,'(3a)') "Write ", trim(filename), " to disk"
    
    n_bin = size(a_rdf%val)    
    
    call xml_OpenFile(filename, xf, indent=.true.)
    
    call xml_AddXMLDeclaration(xf, "UTF-8")
    call xml_NewElement(xf, "rdf")
    
    if (present(temperature) .and. present(density)) then
      call xml_AddAttribute(xf, "title", "T = " // trim(str(temperature, format="(f10.5)")) // &
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
    call xml_AddAttribute(xf, "n-bin", str(n_bin))
    call xml_AddAttribute(xf, "r-max", str(a_rdf%r_max, format="(f10.5)"))
    
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