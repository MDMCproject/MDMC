module rdf_class
use various_constants_class
use structure_class
use near_neighb_class
use histogram_class
use tic_toc_class

implicit none


  ! it is assumed that histogram stores distances in bins of size
  ! bin_length; with the first bin representing distances between
  ! zero and bin_length and the last bin representing distances
  ! between (n-bins-1)*bin_length and n_bins*bin_length
  
  type rdf
    real(db), dimension(:), allocatable :: g_of_r    ! g(r)
    real(db), dimension(:), allocatable :: prefac   ! g(r)=hist%h*prefac
    type (histogram) :: hist 
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
    
    a_rdf%hist = make_histogram(r_max, n_bin)
 
    allocate(a_rdf%g_of_r(n_bin))
    allocate(a_rdf%prefac(n_bin))
    

    temp = volume / (2.0 * pi_value * (n_atoms**2) * a_rdf%hist%bin_length**3 )
    
    do i = 1, n_bin
      r_i = i - 0.5_db
      a_rdf%prefac(i) = temp / (r_i**2)
    end do
  
  end function make_rdf
  
  
  subroutine cal_rdf_nn(hist, nn)
    type (histogram), intent(inout) :: hist
    type (near_neighb_list), intent(in) :: nn
    
    if (nn%what_is_stored == "r2") then
       
    !
    else if (nn%what_is_stored == "r") then
    !
    
    else
    !
    
    end if
    
  end subroutine cal_rdf_nn  
  
  
  subroutine cal_rdf(a_rdf, str)
    type (rdf), intent(inout) :: a_rdf
    type (structure), intent(in) :: str
    
    call cal_histogram(a_rdf%hist, str)
    
    a_rdf%g_of_r = a_rdf%prefac * a_rdf%hist%h

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
    
    n_bin = size(a_rdf%g_of_r)    
    
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
    call xml_AddAttribute(xf, "r-max", str(a_rdf%hist%r_max, format="(f10.5)"))
    
    call xml_NewElement(xf, "this-file-was-created")
    call xml_AddAttribute(xf, "when", get_current_date_and_time())
    call xml_EndElement(xf, "this-file-was-created")
    
    do i = 1, n_bin
      !write(*,'(2f12.6)') (i-0.5)*a_rdf%hist%bin_length, a_rdf%g_of_r(i)
      call xml_NewElement(xf, "g-of-r")
      call xml_AddAttribute(xf, "r", str((i-0.5)*a_rdf%hist%bin_length, format="(f10.5)"))
      call xml_AddAttribute(xf, "g", str(a_rdf%g_of_r(i), format="(f10.5)"))
      call xml_EndElement(xf, "g-of-r")
    end do 
    
    call xml_EndElement(xf, "rdf")
    
    call xml_Close(xf)
    
  end subroutine save_rdf

end module rdf_class