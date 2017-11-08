module control_containers_class
use various_constants_class

implicit none

  ! Contains values to control a MD simulation and comparing with rdf data

  type md_control_container
    integer :: step_limit = 0  ! max number of MD steps
    integer :: average_over_this_many_step   ! when cal properties find average over this many MD steps
    
    real(db) :: temperature   ! stored internally in units of 120.279K (=1000/8.314K)
    real(db) :: md_delta_t     ! stored internally in units of 10e-13seconds
                              
    integer :: total_step_temp_cali = 0    ! the max number of MD steps over which the temperature is adjusted
    integer :: adjust_temp_at_interval ! at what intervals should the temperature be adjusted
                              
    logical :: calculate_rdf = .false.                          
    real(db) :: r_max                  ! cal rdf for interval [0:r_max]
    real(db) :: bin_length             ! at positions 0.5*bin_length, 1.5*bin_length, etc ....
    integer :: cal_rdf_at_interval     ! at what intervals should rdf be calculated
    integer :: average_over_this_many_rdf  ! how many rdf's do you want to average over
  end type md_control_container

  
  ! Contains values to control a MD simulation, MC approach, rdf data as well as S(q, omega) data

  type mdmc_control_container
    ! to initially settle the system into equilibrium
    
    integer :: total_steps_initial_equilibration = 0 ! max number of MD steps
    integer :: average_over_initial_equilibration ! when cal MD properties: ave over this many steps
    integer :: total_step_temp_cali = 0    ! the max number of MD steps where the temperature is adjusted
    integer :: adjust_temp_at_interval ! at what intervals should the temperature be adjusted
                                       ! for this initial calibration stage
     
    ! After the system has settled into equilibrium then perform a number of repeated MC steps, each which
    ! consists of changing PE params and running 'shorter' MD runs with enough MD steps to calculate a 
    ! reasonable histogram and FOM
    
    integer :: md_steps_repeated_equilibration = 0
    integer :: average_over_repeated_equilibration ! when cal MD properties: ave over this many steps    
    integer :: total_step_temp_cali_repeated = 0  
    integer :: adjust_temp_at_interval_repeated    
    logical :: calculate_rdf = .true.
    integer :: cal_rdf_at_interval     ! at what intervals should rdf be calculated
    integer :: average_over_this_many_rdf  ! how many rdf's do you want to average over
    
    ! r-max and bin-length used when wanting to save g(r) to file for algorithm that 
    ! work with structural data only.
    ! For algorithms with compare with dynamical data these determine the r-binning of
    ! histograms created to calculated g(r,t), S(q,t) etc.
    
    real(db) :: r_max                  ! cal rdf for interval [0:r_max]
    real(db) :: bin_length             ! at positions 0.5*bin_length, 1.5*bin_length, etc ....    
    
    ! Used for creating r-binning for histograms to compare with g(r) data
    ! Not used to determine r-binning for histograms to calculate g(r,t)
    ! TODO: see comments else, but these could be removed
    
    integer :: n_r_bin_cal = -1000 
    real(db) :: bin_length_cal = -1000.0   
        
    integer :: mc_steps   ! number of Monte carlo steps  
    
    real(db) :: temperature   ! stored internally in units of 120.279K (=1000/8.314K)
    real(db) :: md_delta_t     ! MD delta-t, stored internally in units of 10e-13seconds
    real(db) :: temperature_mc   ! used in metropolis accept/reject probability
    
    ! Parameters that determine the time binning of histograms used for calculating g(r,t) etc
    ! A g(r,t) time bin has lenght: md_per_time_bin*md_delta_t
    ! Time correlation is calculated over the time = n_time_bin*md_per_time_bin*md_delta_t
    
    integer :: n_time_bin    ! number of time bins
    integer :: md_per_time_bin  ! number of md steps that makes up a time bin
    
    real(db), dimension(:), allocatable :: q_values
    real(db), dimension(:), allocatable :: omega_values
  end type mdmc_control_container


  type (md_control_container) :: setup_md_control_params
  type (mdmc_control_container) :: setup_mdmc_control_params

end module control_containers_class