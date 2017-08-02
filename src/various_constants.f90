module various_constants_class

implicit none

  ! integer, parameter :: db = selected_real_kind(15,307)
  integer, parameter :: db = kind(1.d0)
  
  ! Notice there are places in the code where it is assumed that atoms move
  ! in 3D, in particular in structure_reader.f90
  integer, parameter :: ndim = 3
  
  
  ! this code does it calculations assuming:
  !
  !   m_unit = 1 AMU = 10^-3 KG/MOL
  !   r_unit = 10^-10 METER
  !   t_unit = 10^-13 SEC
  !
  ! which implies:
  !
  !   e_unit = 1 KJ/MOL
  !   T_unit = 1000 / 8.314 K
  !   P_unit = 16387.2 atm
  ! 
 
  real(db), parameter :: T_unit = 120.2790473899446716382_db
  real(db), parameter :: P_unit = 16387.2_db
  
  
  ! To specify when a no-data datapoint (likelihood of having the
  ! number specified below as a real datapoint is close to zero).
  
  real(db), parameter :: no_datapoint_available = 120279047389.944_db
  
    
  ! values below are copied from Ed Akin's figure 2.1
  
  real(db), parameter:: deg_per_rad  = 57.295779513082320876798155_db
  real(db), parameter:: rad_per_deg  = 0.017453292519943295769237_db

  real(db), parameter:: e_value      =  2.71828182845904523560287_db
  real(db), parameter:: e_recip      =  0.3678794411714423215955238_db
  real(db), parameter:: e_squared    =  7.389056098930650227230427_db
  real(db), parameter:: log10_of_e   =  0.4342944819032518276511289_db

  real(db), parameter:: euler        =  0.5772156649015328606_db
  real(db), parameter:: euler_log    = -0.5495393129816448223_db
  real(db), parameter:: gamma        =  0.577215664901532860606512_db
  real(db), parameter:: gamma_log    = -0.549539312981644822337662_db
  real(db), parameter:: golden_ratio =  1.618033988749894848_db

  real(db), parameter:: ln_2         =  0.6931471805599453094172321_db
  real(db), parameter:: ln_10        =  2.3025850929940456840179915_db
  real(db), parameter:: log10_of_2   =  0.3010299956639811952137389_db

  real(db), parameter:: pi_value     =  3.141592653589793238462643_db
  real(db), parameter:: pi_ln        =  1.144729885849400174143427_db
  real(db), parameter:: pi_log10     =  0.4971498726941338543512683_db
  real(db), parameter:: pi_over_2    =  1.570796326794896619231322_db
  real(db), parameter:: pi_over_3    =  1.047197551196597746154214_db
  real(db), parameter:: pi_over_4    =  0.7853981633974483096156608_db
  real(db), parameter:: pi_recip     =  0.3183098861837906715377675_db
  real(db), parameter:: pi_squared   =  9.869604401089358618834491_db
  real(db), parameter:: pi_sq_root   =  1.772453850905516027298167_db

  real(db), parameter:: sq_root_of_2 =  1.4142135623730950488_db
  real(db), parameter:: sq_root_of_3 =  1.7320508075688772935_db  

end module various_constants_class
