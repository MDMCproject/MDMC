module phasespace_class
use structure_class
use near_neighb_class
use common_potential_block_class

  implicit none


  ! This type is only include to speed up the calculation. Most importantly by
  ! replacing the double summation over atoms by a sum over nearest neighboors
  ! and secondly by storing distances that may be used at a later stage.
  ! Notice that this distances are here stored as REAL rather then REAL (DB)
  ! to save memory, since this type will be the single most memory hungry 
  ! component of this code and knowing these distances to 6 digits accuracy
  ! is e.g. definitely good enough for calculating a histogram, pdf etc from

  type phasespace
    type (structure) :: str
  
    type (near_neighb_list) :: neighb_list
  end type phasespace


contains

  subroutine trajectory_in_phasespace(ps, n_md, delta_t)
    type (phasespace), intent(inout) :: ps
    integer, intent(in) :: n_md
    real(db), intent(in) :: delta_t   
    
    real(db) :: pot_energy
    
    pot_energy = gpe_val(ps%str, common_gpe)
                         
    write (*, '(a,f10.3)') "Epot = ", pot_energy
  
  end subroutine trajectory_in_phasespace

  function make_phasespace(str, r_cut, delta_r) result (ps)
    type (structure), intent(in) :: str
    real (db), intent(in) :: r_cut, delta_r
    type (phasespace) :: ps

    integer :: n_tot
    
    n_tot = size(str%atoms)
    
    ps%neighb_list = make_near_neighb_list(str, r_cut, delta_r)
    
    ps%str = copy_structure(str)
    
  end function make_phasespace


end module phasespace_class
