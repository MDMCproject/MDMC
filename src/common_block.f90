module common_block_class
use configuration_class

  type(configuration) :: common_configuration

  ! Used to 'simulate' dynamical polymorphism, see module
  ! configuration_class. I may not need to use targets here
  ! but could instead allocate memery directly from the
  ! constraint pointer in common_configuration rather than using
  ! target, but not convinced that this would work so for
  ! now do it using targets (I would be easy to change later)
  type (fnc_constraint), target :: my_fnc_constraint

end module common_block_class
