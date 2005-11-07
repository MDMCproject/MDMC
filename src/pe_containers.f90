! include in this file all pe containers

module lj_pe_container_class
use variable_class

implicit none


  type lj_pe_container
    type (variable), dimension(:), allocatable :: vars
  end type lj_pe_container


end module lj_pe_container_class
