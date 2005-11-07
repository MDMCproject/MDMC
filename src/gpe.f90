module gpe_class
use various_constants_class
use structure_class
use variable_class

implicit none

  type potential
    character(len=99)  :: name = "fisse"
    type (variable), dimension(:), allocatable :: vars
  end type potential

	integer, parameter, private    :: MAX_ITEMS = 64
	type (potential), dimension(MAX_ITEMS) :: common_pe_list

contains

  function gpe_val(str, vars) result (val)
    type (structure), intent(in) :: str
		type (variable), dimension(:), intent(in) :: vars
    real (db) :: val
    
    integer :: i1, i2, n_tot, i
    real (db) :: r_cut, rr_cut 
    real (db) :: sigma, epsilon, epsilon_times4
    real (db) :: rr, rri, rri3
		real (db), dimension(ndim) :: diff_vec
    
    
		sigma = vars(1)%val
		epsilon = vars(2)%val
    epsilon_times4 = 4*epsilon
    
    r_cut = 2**0.166666666666667 * sigma
    rr_cut = r_cut * r_cut
    
		n_tot = size(str%atoms)
    
    val = 0.0
    
		do i1 = 1, n_tot
		  do i2 = i1+1, n_tot
        diff_vec = str%atoms(i1)%r - str%atoms(i2)%r
        
        do i = 1, ndim
          if (diff_vec(i) >= 0.5 * str%box_edges(i)) then
            diff_vec(i) = diff_vec(i) - str%box_edges(i)
          end if
          if (diff_vec(i) < -0.5 * str%box_edges(i)) then
            diff_vec(i) = diff_vec(i) + str%box_edges(i)
          end if       
        end do
        
        rr = sum(diff_vec*diff_vec)
        
        if (rr < rr_cut) then
          rri = sigma*sigma / rr
          rri3 = rri * rri * rri
          val = val + epsilon_times4 * rri3 * (rri3 - 1.0) + 1.0   
        end if
			
			end do
		end do
    
  end function gpe_val
  
  
  subroutine gpe_deriv(str,deriv)
    type (structure), intent(in) :: str
    real(db), dimension(:,:), intent(out) :: deriv
    
    integer :: i
    

  end subroutine gpe_deriv
    
    
end module gpe_class
