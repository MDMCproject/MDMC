module control_containers_class
use various_constants_class

implicit none


  type md_control_container
    integer :: step_limit   ! max number of MD steps
    integer :: average_over_this_many_step   ! when cal properties find ave,esd over this many MD steps
    
    real(db) :: temperature   ! stored internally in units of 120.279K (=1000/8.314K)
    real(db) :: time_step     ! stored internally in units of 10e-13seconds
                              
    !logical :: perform_initial_temperature_calibration = .false.
    integer :: total_step_temp_cali = 0    ! the total number of MD steps where the temperature is adjusted
    integer :: adjust_temp_at_interval ! at what intervals should the temperature be adjusted
                              ! in this initial calibration stage. Notice that adjust_temp_at_interval
                              ! most be equal or smaller than total_step_temp_cali
                              
    logical :: calculate_rdf = .false.                          
    real(db) :: r_max                  ! cal rdf for interval [0:r_max]
    real(db) :: bin_length             ! at positions 0.5*bin_length, 1.5*bin_length, etc ....
    integer :: cal_rdf_at_interval     ! at what intervals should rdf be calculated
    integer :: average_over_this_many_rdf  ! how many rdf's do you want to average over
  end type md_control_container


  type mdmc_control_container
    ! to initially settle system
    integer :: total_steps_initial_equilibration
    integer :: average_over_initial_equilibration  ! when cal MD properties: ave over this many steps
    integer :: total_step_temp_cali = 0    ! the total number of MD steps where the temperature is adjusted
    integer :: adjust_temp_at_interval ! at what intervals should the temperature be adjusted
                              ! in this initial calibration stage. Notice that adjust_temp_at_interval
                              ! most be equal or smaller than total_step_temp_cali
     
    ! After the system has settle into equilibrium then perform a number of repeated MC steps, which
    ! consist of changing Epot params, run a 'short' MD run to get the new system into a state close
    ! to equilibrium, and finally enough MD steps to calculate a reasonable histogram and FOM.                          
    integer :: md_steps_repeated_equilibration
    integer :: average_over_repeated_equilibration ! when cal MD properties: ave over this many steps    
    integer :: total_step_temp_cali_repeated = 0    
    integer :: adjust_temp_at_interval_repeated    
    logical :: calculate_rdf = .true.
    integer :: cal_rdf_at_interval     ! at what intervals should rdf be calculated
    integer :: average_over_this_many_rdf  ! how many rdf's do you want to average over
    
    ! r-max and bin-length are here only used when wanting to save g(r) to file
    real(db) :: r_max                  ! cal rdf for interval [0:r_max]
    real(db) :: bin_length             ! at positions 0.5*bin_length, 1.5*bin_length, etc ....    
    
    ! Used for creating a histogram or time_corr_hist that has a binning and 
    ! length which is determined by the binning of data plus optionally some
    ! additional fom elements (like <r-max>)
    integer :: n_r_bin_cal = -1000 
    real(db) :: bin_length_cal = -1000.0   
        
    integer :: mc_steps   ! number of Monte carlo steps  
    
    real(db) :: temperature   ! stored internally in units of 120.279K (=1000/8.314K)
    real(db) :: time_step     ! stored internally in units of 10e-13seconds
    real(db) :: temperature_mc   ! used in metropolis
    
    ! parameters for time correlation run
    integer :: n_time_bin    ! number of time bins
    integer :: md_per_time_bin  ! number of md steps in a time bin. Robert say this one probably always 1 
    
    !real(db) :: g_d_data_time_step = 0.0 ! VERY BAD WHAT I AM DOING HERE
    
    real(db), dimension(:), allocatable :: q_values
    real(db), dimension(:), allocatable :: omega_values
    
  end type mdmc_control_container


  type (md_control_container) :: setup_md_control_params
  type (mdmc_control_container) :: setup_mdmc_control_params

end module control_containers_class