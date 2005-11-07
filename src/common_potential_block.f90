module common_potential_block_class
use gpe_class

  type(pe_list) :: common_gpe

  ! Used to 'simulate' dynamical polymorphism, see module
  ! configuration_class. I may not need to use targets here
  ! but could instead allocate memery directly from the
  ! constraint pointer in common_config rather than using
  ! target, but not convinced that this would work so for
  ! now do it using targets (I would be easy to change later)
  type (lj_pe_container), target :: target_lj_pe

end module common_potential_block_class

