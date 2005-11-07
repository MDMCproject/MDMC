module various_constants_class

implicit none

  ! integer, parameter :: db = selected_real_kind(15,307)
  integer, parameter :: db = kind(1.d0)
  
  ! Notice there are places in the code where it is assumed that atoms move
  ! in 3D, in particular in structure_reader.f90
  integer, parameter :: ndim = 3
  

end module various_constants_class
