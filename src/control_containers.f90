module control_containers_class
use various_constants_class

implicit none


  type md_control_container
    integer :: step_limit   ! max number of MD steps
    integer :: average_over_this_many_step   ! when cal properties find ave,esd over this many MD steps
    
    real(db) :: temperature   ! stored internally in units of 120.279K (=1000/8.314K)
    real(db) :: time_step     ! stored internally in units of 10e-13seconds
    
    !logical :: use_near_neighbour_method = .false.
    real(db) :: r_cut = 0.0   ! cal the potential only to this value. If r_cut == 0.0 it is assumed that
                              ! use_near_neighbour_method = .false.
    real(db) :: delta_r       ! buffer region so that the near-neighbour pairs includes
                              ! all pairs with distances less than r_cut+delta_r
                              
    !logical :: perform_initial_temperature_calibration = .false.
    integer :: total_step_temp_cali = 0    ! the total number of MD steps where the temperature is adjusted
    integer :: adjust_temp_at_interval ! at what intervals should the temperature be adjusted
                              ! in this initial calibration stage. Notice that adjust_temp_at_interval
                              ! most be equal or smaller than total_step_temp_cali
                              
    logical :: calculate_rdf = .false.                          
    real(db) :: r_max                  ! cal rdf for interval [0:r_max]
    integer :: number_bins             ! over this many bins
    integer :: cal_rdf_at_interval     ! at what intervals should rdf be calculated
    integer :: average_over_this_many_rdf  ! how many rdf's do you want to average over
  end type md_control_container


  type mdmc_control_container
    ! to initially settle system
    integer :: total_steps_to_settle_system
    integer :: average_over_this_many_step
    integer :: total_step_temp_cali = 0    ! the total number of MD steps where the temperature is adjusted
    integer :: adjust_temp_at_interval ! at what intervals should the temperature be adjusted
                              ! in this initial calibration stage. Notice that adjust_temp_at_interval
                              ! most be equal or smaller than total_step_temp_cali
                              
    ! a MD trajectory is split up into two parts. First part consist of a number of
    ! steps to bring the system not necessarily to equilibrium but something perhaps
    ! close to equilibirum. The second part is where a rdf is calculated
    integer :: md_steps_first_part
    logical :: calculate_rdf = .true.
    real(db) :: r_max                  ! cal rdf for interval [0:r_max]
    integer :: number_bins             ! over this many bins
    integer :: cal_rdf_at_interval     ! at what intervals should rdf be calculated
    integer :: average_over_this_many_rdf  ! how many rdf's do you want to average over
        
    integer :: mc_steps   ! number of Monte carlo steps  
      
    integer :: md_steps_per_trajectory   ! number of MD steps per trajectory

    
    real(db) :: temperature   ! stored internally in units of 120.279K (=1000/8.314K)
    real(db) :: time_step     ! stored internally in units of 10e-13seconds
    real(db) :: temperature_mc   ! used in metropolis
    
    !logical :: use_near_neighbour_method = .false.
    real(db) :: r_cut = 0.0   ! cal the potential only to this value. If r_cut == 0.0 it is assumed that
                              ! use_near_neighbour_method = .false.
    real(db) :: delta_r       ! buffer region so that the near-neighbour pairs includes
                              ! all pairs with distances less than r_cut+delta_r
  end type mdmc_control_container


  type (md_control_container) :: setup_md_control_params
  type (mdmc_control_container) :: setup_mdmc_control_params

end module control_containers_class