module md_control_time_corr_class
use configuration_class
use function_class
use common_block_class, only : common_config, common_pe_list, common_fom_list, target_rdf_fom
use phasespace_class
use md_properties_class
use control_containers_class
use tic_toc_class
use rdf_class
use structure_reader_class
use func_params_wrapper_class
use time_corr_algorithm_class
use time_corr_hist_container_class
use s_q_time_class
use s_q_omega_class

  implicit none

  public :: run_md_control_time_corr
  
  private :: acceptable_temperature
  private :: acceptable_energy  
  
contains

  subroutine run_md_control_time_corr(a_config, c)
    use flib_wxml  
    type (configuration), intent(inout) :: a_config
    type (mdmc_control_container) :: c
    
    real (db) :: time_now = 0.0    
    integer :: i, i_md, j
    real(db) :: density  ! used for cal input argument to save_rdf
    
    real(db) :: sum_kin_energy
    real(db) :: temp_adjust_factor
    
    real(db) :: average_energy_end_of_temp_calibration
  
    real(db) :: fom_val, fom_old

    real(db) :: fom_best  ! used for keeping track of best fom
    type (phasespace) :: my_ps_best ! to store status of best phasespace

    type (phasespace) :: my_ps, my_ps_old
    type (md_properties) :: my_props
    real(db) :: pressure_comp = 0.0, pot_energy = 0.0

    type(time_corr_hist_container) :: my_time_corr_container
    
    type (s_q_time) :: my_s_q_time
    type (s_q_omega) :: my_s_q_omega
	
    integer :: print_to_file = 555
    integer :: print_to_screen = 0
    
    real (db) :: ran_num  ! for random generator
    real (db) :: delta_fom ! for metropolis
    logical :: accept_parameters ! for metropolis
    
    type (xmlf_t) :: xf
    call xml_OpenFile("output/mdmc_results.xml", xf, indent=.true.)
    call xml_AddXMLDeclaration(xf, "UTF-8")
    call xml_NewElement(xf, "md-control-time-corr-results")
    call xml_AddAttribute(xf, "title", "rho = " // str(density, format="(f10.5)") &
                                         // "atoms/AA-3")
    call xml_NewElement(xf, "this-file-was-created")
    call xml_AddAttribute(xf, "when", get_current_date_and_time())
    call xml_EndElement(xf, "this-file-was-created")
	  
    if (print_to_file /= 0) then
      open(print_to_file, file="output/job_log.txt")
    end if

    write(print_to_file,*) "In run_md_control"

    call tic

    ! prepare phasespaces
    
    my_ps = make_phasespace(a_config%str, c%temperature)
    my_ps_old = copy_phasespace(my_ps)
    my_ps_best = copy_phasespace(my_ps)
                                                                   
    
    ! prepare space and time correlation histogram container
    
    my_time_corr_container = make_time_corr_hist_container(c%r_max, c%bin_length, c%n_time_bin, &
                             c%md_per_time_bin * c%md_delta_t)
    call clear_time_corr_hist_container(my_time_corr_container)                         
    
    ! to cal and print out s_q_time and s_q_omega
    
    my_s_q_time = make_s_q_time(c%q_values, c%md_per_time_bin * c%md_delta_t, c%n_time_bin)
    my_s_q_omega = make_s_q_omega(c%q_values, c%omega_values)                                         
                        
               
! -------------- initial equilibration ---------------- !

    sum_kin_energy = 0.0
    
    do i = 1, c%total_steps_initial_equilibration
      time_now = c%md_delta_t * i 
      
      ! do one trajectory of length = 1 where pressure_comp and pot_energy is also
      ! calculated
      
      call trajectory_in_phasespace(my_ps, common_pe_list, 1, c%md_delta_t, & 
                                    pressure_comp, pot_energy)
      
      call md_cal_properties(my_ps, my_props, common_pe_list, pressure_comp, pot_energy)
      
      ! Optionally adjust the temperature

      if (i < c%total_step_temp_cali) then
        sum_kin_energy = sum_kin_energy + my_props%kin_energy%val
        if (mod(i,c%adjust_temp_at_interval) == 0) then
          ! see src/algorithms/md_gridsearch_control.doc for explanation of the expression below
          !        
          temp_adjust_factor = sqrt(c%adjust_temp_at_interval * 1.5 * c%temperature / &
              sum_kin_energy * (size(my_ps%str%atoms)-1.0) / size(my_ps%str%atoms))
          my_ps%p = my_ps%p * temp_adjust_factor
          sum_kin_energy = 0.0
        end if
      end if
        
      ! The idea is to keep a record of what the average total energy is after the
      ! temperature calibration and then compare this value at end of the entire
      ! initial calibration to check if the total energy has not drifted due to 
      ! numerical errors. Not sure if this is the best point to record this total average
      
      if (i == c%total_step_temp_cali) then
        average_energy_end_of_temp_calibration = my_props%tot_energy%ave 
      end if
      
      ! print out stuff at interval = average_over_this_many_steps
        
      if (mod(i,c%average_over_initial_equilibration) == 0) then 
        call md_print_properties(print_to_file, my_props)
        
        if (sum(sum(my_ps%p,1)) > 0.0001) then
          write(*,*) "ERROR:"
          write(*,'(a,3f12.6)') "total momentum ", sum(my_ps%p,1)
          write(*,*) "Serious problem - total momentum different from zero"
          stop
        end if
        
        call md_reset_properties(my_props)
        write(*, '(a,i8,a,f15.5,a)') "MD steps = ", i, " MD run-time = ", time_now, "*10e-13"
      end if
    end do
    
    write(print_to_file, *) " "
    write(print_to_file, '(a,f12.4,a)') "Taken ", toc(), " seconds to execute initial equilibration."
    write(print_to_file, *) " "

    ! Determine if equilibrium was reached
    !
    ! Note the (2.0/ndim) factor is to convert from dimensionless kin_energy per atom to 
    ! dimensionless temperature
    
    if (.not. acceptable_temperature((2.0/ndim)*my_props%kin_energy%ave,c%temperature)) then
         write(print_to_screen, *) "Initial equilibration did not reach equilibrium"
         write(print_to_screen, *) "Temperature outside acceptable value - STOP"
         stop
    end if   
    
    if (.not. acceptable_energy(average_energy_end_of_temp_calibration, &
         my_props%tot_energy%ave)) then
         write(print_to_screen, *) "Initial equilibration did not reach equilibrium - STOP"
         write(print_to_screen, *) "Energy outside acceptable value - STOP"
         stop
    end if     
 
 ! -------- calculate and save g(r,t), S(q,t) and S(q,omega) ------------------------- !

  density = size(my_ps%str%atoms) / product(my_ps%str%box_edges) ! for printing

  write(print_to_screen, *) "Start calculating space-time histograms"  
  call cal_time_corr_container(my_time_corr_container, my_ps, common_pe_list, c%md_per_time_bin, c%md_delta_t)
  write(print_to_screen, *) "Write space-time correlations functions to disk"  
  call print_g_d(my_time_corr_container, product(my_ps%str%box_edges), size(my_ps%str%atoms), c%temperature)
  call print_h_d_hist(my_time_corr_container, density, c%temperature)
  call print_g_s(my_time_corr_container, density, c%temperature)
  call print_h_s_hist(my_time_corr_container, density, c%temperature)
  call print_einstein_diffuse_exp(my_time_corr_container, size(my_ps%str%atoms), density, c%temperature)
  
  write(print_to_screen, *) "Start calculating intermediate scattering function" 
  call cal_s_q_time(my_time_corr_container, my_ps%str, my_s_q_time)
  write(print_to_screen, *) "Write intermediate scattering function to disk" 
  call print_s_q_time(my_s_q_time, density, c%temperature)
  
  write(print_to_screen, *) "Start calculating dynamical structure factor function"   
  call cal_s_q_omega(my_s_q_time, my_ps%str, my_s_q_omega)
  write(print_to_screen, *) "Write dynamical structure factors to disk"   
  call print_s_q_omega(my_s_q_omega, density, c%temperature)  
 
 ! -------- calculate and save g(r,t), S(q,t) and S(q,omega) -------------------------- !

    print *, ' '
    print *, 'Job took ', toc(), ' seconds to execute.'
    
    if (print_to_file /= 0) then
	    close(print_to_file)
    end if    
	  
    call xml_EndElement(xf, "md-control-time-corr-results")
    call xml_Close(xf)
    
  end subroutine run_md_control_time_corr
  
  
  ! check to see if temperature is within certain limits of t_target
  !
  function acceptable_temperature(t, t_target) result (yes_or_no)
    real (db), intent(in) :: t, t_target
    logical :: yes_or_no 
    
    yes_or_no = .true.
    
    if (abs((t - t_target)/t_target) > 0.2) then
      yes_or_no = .false.
    end if
  end function acceptable_temperature
  
  
  ! check to see if energy is within certain limits of t_target
  !
  function acceptable_energy(e, e_target) result (yes_or_no)
    real (db), intent(in) :: e, e_target
    logical :: yes_or_no 
    
    yes_or_no = .true.
    
    if (abs((e - e_target)/e_target) > 0.1) then
      yes_or_no = .false.
    end if
  end function acceptable_energy

end module md_control_time_corr_class
