module mdmc_control_class
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

  implicit none

contains

  subroutine run_mdmc_control(a_config, c)
    type (configuration), intent(inout) :: a_config
    type (mdmc_control_container) :: c
    
    real (db) :: time_now = 0.0    
    integer :: i, i_md
    real(db) :: density  ! used for cal input argument to save_rdf
    
    real(db) :: sum_kin_energy = 0.0
    real(db) :: temp_adjust_factor
  
    ! md storage containers

    type (phasespace) :: my_ps
    type (md_properties) :: my_props
    real(db) :: pressure_comp = 0.0, pot_energy = 0.0
    type (rdf) :: my_rdf, my_rdf_sum
	
		
    write(*,*) "In run_md_control"

    call tic

    ! initiate phasespace
    
    my_ps = make_phasespace(a_config%str, c%temperature)
                        
                        
    ! Notice a design issue: fix this later
    ! the problem is that make_phasespace does not currently properly initiate 
    ! the histogram (not a problem to fix later)
    
!    my_ps%neighb_list%hist = copy_histogram(target_rdf_fom%rdf_data%hist)
                                               
                        
               
! -------------- settling the system ---------------- !


    do i = 1, c%total_steps_to_settle_system
      time_now = c%time_step * i   ! perhaps print this one out 
      
      ! do one trajectory of length = 1 where pressure_comp and pot_energy is also
      ! calculated
      
      !call trajectory_in_phasespace(my_ps, common_pe_list, 1, c%time_step)
      call trajectory_in_phasespace(my_ps, common_pe_list, 1, c%time_step, & 
                                    pressure_comp, pot_energy)
      
      !call md_cal_properties(my_ps, my_props, common_pe_list)
      call md_cal_properties(my_ps, my_props, common_pe_list, pressure_comp, pot_energy)
      
      ! case you want to adject the temperature in the initial stages of the MD simulation
      ! (notice c%total_step_temp_cali = 0 if <perform-initial-temperature-calibration> 
      ! element not specified in input file)

      if (i < c%total_step_temp_cali) then
        sum_kin_energy = sum_kin_energy + my_props%kin_energy%val
        if (mod(i,c%adjust_temp_at_interval) == 0) then
          temp_adjust_factor = sqrt(c%adjust_temp_at_interval * 1.5 * c%temperature / &
              sum_kin_energy * (size(my_ps%str%atoms)-1.0) / size(my_ps%str%atoms))
          my_ps%p = my_ps%p * temp_adjust_factor
          sum_kin_energy = 0.0
        end if
      end if
     
      

      ! accumulate the calculated MD property values
        
      call md_accum_properties(my_props)
        
        
      ! print out stuff and interval = average_over_this_many_steps
        
      if (mod(i,c%average_over_this_many_step) == 0) then 
        call md_print_properties(my_props)
        
        if (sum(sum(my_ps%p,1)) > 0.0001) then
          write(*,*) "ERROR:"
          write(*,'(a,3f12.6)') "total momentum ", sum(my_ps%p,1)
          write(*,*) "Serious problem - total momentum different from zero"
          stop
        end if
        
        call md_reset_properties(my_props)
        write(*, '(a,i8,a,f12.4,a)') "MD steps = ", i, " MD run-time = ", time_now, "*10e-13"
      end if
      
      
    end do
    
    
 ! ---------------  finished settling the system -------------- !               
               
               
 ! --------------- calculate rdf and first fom   -------------- !
 
 !   do i = 1, c%average_over_this_many_rdf


 !     do j = 1, c%average_over_this_many_rdf      
      
        ! do one trajectory of length = 1 where pressure_comp and pot_energy is also
        ! calculated
      
        !call trajectory_in_phasespace(my_ps, common_pe_list, 1, c%time_step)
  !      call trajectory_in_phasespace(my_ps, common_pe_list, 1, c%time_step, & 
  !                                  pressure_comp, pot_energy)
      
        !call md_cal_properties(my_ps, my_props, common_pe_list)
  !      call md_cal_properties(my_ps, my_props, common_pe_list, pressure_comp, pot_energy)              
      
               
        ! accumulate the calculated MD property values
        
  !      call md_accum_properties(my_props)
  !    end do
                     
  !    call nn_update_histogram(my_ps%neighb_list)
            
            
  !    my_hist_sum%h = my_hist_sum%h + my_ps%neighb_list%hist%h
      
  !  end do    
          
            
  !my_rdf_sum%g_of_r = my_rdf_sum%g_of_r / c%average_over_this_many_rdf
            
  !density = size(my_ps%str%atoms)/product(my_ps%str%box_edges)
  !          call save_rdf(my_rdf_sum, c%temperature, density)
            
  !my_rdf_sum%g_of_r = 0.0
          
        
        
        
               
               
               
               
               
               
               
               
               
                        

    do i = 1, c%mc_steps
      time_now = c%time_step * i * c%md_steps_per_trajectory  ! perhaps print out 
      
      ! do one trajectory
      
      do i_md = 1, c%md_steps_per_trajectory
      
        !call trajectory_in_phasespace(my_ps, common_pe_list, 1, c%time_step)
        call trajectory_in_phasespace(my_ps, common_pe_list, 1, c%time_step, & 
                                      pressure_comp, pot_energy)
        
        !call md_cal_properties(my_ps, my_props, common_pe_list)
        call md_cal_properties(my_ps, my_props, common_pe_list, pressure_comp, pot_energy)


        ! accumulate the calculated MD property values
          
        call md_accum_properties(my_props)
      
        sum_kin_energy = sum_kin_energy + my_props%kin_energy%val
      end do  
        
      ! print out stuff
        
      call md_print_properties(my_props)
        
      if (sum(sum(my_ps%p,1)) > 0.0001) then
        write(*,*) "ERROR:"
        write(*,'(a,3f12.6)') "total momentum ", sum(my_ps%p,1)
        write(*,*) "Serious problem - total momentum different from zero"
        stop
      end if
      
      call md_reset_properties(my_props)
      write(*, '(a,i8,a,f12.4,a)') "MD steps = ", i, " MD run-time = ", time_now, "*10e-13"
      
      !call nn_update_histogram(my_ps%neighb_list)
      
!      print *, "FOM = ", func_val_nn(my_ps%str, common_fom_list, my_ps%neighb_list)
            
      
      ! adjust the temperature
      temp_adjust_factor = sqrt(c%md_steps_per_trajectory * 1.5 * c%temperature / &
              sum_kin_energy * (size(my_ps%str%atoms)-1.0) / size(my_ps%str%atoms)) 
      my_ps%p = my_ps%p * temp_adjust_factor
      sum_kin_energy = 0.0
      
      
    end do
   
    call save_structure(my_ps%str, "output/argon_structure.xml")
    
    print *, ' '
    print *, 'Job took ', toc(), ' seconds to execute.'
  end subroutine run_mdmc_control

end module mdmc_control_class
