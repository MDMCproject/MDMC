! This class used to 'simulate' dynamical polymorphism
! f90 shortcoming (compared to a fully OO language) 

module common_block_class
use configuration_class
use function_class

  ! Used for reading in initial configuration from XML input
  ! and may than be altered in subsequent algorithms
  
  type(configuration) :: common_config

  type (fnc_constraint), target :: target_fnc_constr

  type(func_list) :: common_pe_list
  type(func_list) :: common_fom_list
  type (lj_pe_container), target :: target_lj_pe
  type (rdf_fom_container), target :: target_rdf_fom
  type (g_d_rt_fom_container), target :: target_g_d_rt_fom
  type (s_qt_fom_container), target :: target_s_qt_fom
  type (s_qo_fom_container), target :: target_s_qo_fom


end module common_block_class
