module common_block_class
use configuration_class

  type(configuration) :: common_configuration
  type (fnc_constraint), target :: my_fnc_constraint

end module common_block_class
