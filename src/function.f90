module function_class
use various_constants_class
use structure_class
use func_param_class
use lennard_jones_class
use rdf_fom_class
use g_d_rt_fom_class
use s_q_time_fom_class

implicit none

  public :: func_val
  public :: func_deriv

  private :: add_lj_pe_container
  private :: add_rdf_fom_container

  type func_list
    !private
    type (lj_pe_container), pointer :: pt_lj_pe => null()
    type (rdf_fom_container), pointer :: pt_rdf_fom => null()
    type (g_d_rt_fom_container), pointer :: pt_g_d_rt_fom => null()
    type (s_qt_fom_container), pointer :: pt_s_qt_fom => null()
  end type func_list

  interface add_function
    module procedure add_lj_pe_container
    module procedure add_rdf_fom_container
    module procedure add_g_d_rt_fom_container
    module procedure add_s_qt_fom_container
  end interface
  
  interface func_val
    module procedure func_val_histogram
    module procedure func_val_structure
    module procedure func_val_time_corr
    module procedure func_val_s_q_time
  end interface  
  
contains

  subroutine add_lj_pe_container(list, pe_target)
    type (func_list), intent(inout) :: list
    type (lj_pe_container), target, intent(in) :: pe_target
    
    list%pt_lj_pe => pe_target
  end subroutine add_lj_pe_container
    
  subroutine add_rdf_fom_container(list, fom_target)
    type (func_list), intent(inout) :: list
    type (rdf_fom_container), target, intent(in) :: fom_target
    
    list%pt_rdf_fom => fom_target
  end subroutine add_rdf_fom_container  

  subroutine add_g_d_rt_fom_container(list, fom_target)
    type (func_list), intent(inout) :: list
    type (g_d_rt_fom_container), target, intent(in) :: fom_target
    
    list%pt_g_d_rt_fom => fom_target
  end subroutine add_g_d_rt_fom_container  
  
  subroutine add_s_qt_fom_container(list, fom_target)
    type (func_list), intent(inout) :: list
    type (s_qt_fom_container), target, intent(in) :: fom_target
    
    list%pt_s_qt_fom => fom_target
  end subroutine add_s_qt_fom_container    
  
    
  function func_val_structure(str, list) result (val)
    type (structure), intent(in) :: str
    type (func_list), intent(inout) :: list
    real (db) :: val

    if ( associated(list%pt_lj_pe) ) then
      val = lj_val(str, list%pt_lj_pe)
    else if ( associated(list%pt_rdf_fom) ) then
      val = rdf_fom_val(str, list%pt_rdf_fom)
    else
      print *, "Error in func_val_structure"
      stop
    end if
  end function func_val_structure

  function func_val_histogram(hist, list) result (val)
    type (histogram), intent(in) :: hist
    type (func_list), intent(inout) :: list
    real (db) :: val

    if ( associated(list%pt_lj_pe) ) then
      print *, "Problem in func_val_histogram"
      stop
    else if ( associated(list%pt_rdf_fom) ) then
      val = rdf_fom_val(hist, list%pt_rdf_fom)
    else
      print *, "Error in func_val_histogram"
      stop
    end if
  end function func_val_histogram  
  
  function func_val_time_corr(time_corr, list) result (val)
    type (time_corr_hist_container), intent(in) :: time_corr
    type (func_list), intent(inout) :: list
    real (db) :: val

    if ( associated(list%pt_lj_pe) ) then
      print *, "Problem in func_val_time_corr"
      stop
    else if ( associated(list%pt_rdf_fom) ) then
      print *, "Problem in func_val_time_corr"
      stop
    else if ( associated(list%pt_g_d_rt_fom) ) then
      val = g_d_rt_fom_val(time_corr, list%pt_g_d_rt_fom)
    else
      print *, "Error in func_val_time_corr"
      stop
    end if
  end function func_val_time_corr
  
  function func_val_s_q_time(sqtime, list) result (val)
    type (s_q_time), intent(in) :: sqtime
    type (func_list), intent(inout) :: list
    real (db) :: val

    if ( associated(list%pt_lj_pe) ) then
      print *, "Problem in func_val_s_q_time"
      stop
    else if ( associated(list%pt_rdf_fom) ) then
      print *, "Problem in func_val_s_q_time"
      stop
    else if ( associated(list%pt_g_d_rt_fom) ) then
      print *, "Problem in func_val_s_q_time"
      stop
    else if ( associated(list%pt_s_qt_fom) ) then
      val = s_qt_fom_val(sqtime, list%pt_s_qt_fom)      
    else
      print *, "Error in func_val_s_q_time"
      stop
    end if
  end function func_val_s_q_time  
  
    
  subroutine func_deriv(str, deriv, list, pressure_comp, pot_energy)
    type (structure), intent(in) :: str
    real(db), dimension(:,:), intent(out) :: deriv
    type (func_list), intent(in) :: list
 		real (db), optional, intent(out) :: pressure_comp, pot_energy
 		
 		! test for optional arguments
    
    if (present(pressure_comp) .and. present(pot_energy)) then
      call lj_deriv(str, deriv, list%pt_lj_pe, pressure_comp, pot_energy)
    else if (present(pressure_comp)==.false. .and. present(pot_energy)==.false.) then
      call lj_deriv(str, deriv, list%pt_lj_pe)
    else
      write(*,*) "ERROR in func_deriv"
      write(*,*) "Either both pressure_comp and pot_energy must be present"
      write(*,*) "none of them."
      stop    
    end if 
 		
  end subroutine func_deriv

      
end module function_class
