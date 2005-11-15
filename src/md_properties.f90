module md_properties_class
use phasespace_class

implicit none

  public :: md_cal_properties, md_cal_properties_extra
  public :: md_accum_properties
  public :: md_reset_properties
  public :: md_print_properties

  type property
    real(db) :: val   
    real(db) :: sum = 0.0, sum2 = 0.0
    real(db) :: ave, esd
  end type property

  type md_properties
    type (property) :: kin_energy
    type (property) :: tot_energy
    type (property) :: pressure
    
    ! this one in theory should be in type property but
    ! since in the code the below the properties are 
    ! all accumulated in one go I have moved it down for
    ! now
    integer :: n_accum = 0
  end type md_properties

contains

  subroutine md_cal_properties(ps, props, list)
    type (phasespace), intent(in) :: ps
    type (md_properties), intent(out) :: props
    type (pe_list), intent(in) :: list
    
    real(db) :: pot_energy, sum_mass_v2
    integer :: n_atoms
    
    n_atoms = size(ps%str%atoms)
    
    ! here assumes that ps%v2 is calculated at t=t+h/2
    ! and therefore needs to be re-calculated so that
    ! we calculate both the potential and kinetic energy
    ! at t=t+h
    
    sum_mass_v2 = sum(ps%p * ps%p * ps%inv_mass)
    
    ! unit of kinetic energy below is KJ per mol per atom (assuming
    ! mass_unit=amu=g/mol, t_unit=10^-13s, r_unit=AA), so do get the total
    ! energy of the system in 'absolute' units you would multiply this
    ! number by n_atoms and divide by Avogadro's number
    
    props%kin_energy%val = 0.5*sum_mass_v2 / n_atoms
    
    if (ps%neighb_list%ignore_list == .true.) then
      pot_energy = gpe_val(ps%str, list) / n_atoms
    else
      pot_energy = gpe_val_nn(ps%str, list, ps%neighb_list) / n_atoms
    end if
    
    props%tot_energy%val = pot_energy + props%kin_energy%val
    
    
    ! pressure is in units of 16387.72 atm
    ! the expression below is only valid when using non-periodic boundary conditions
    !props%pressure%val = (sum_mass_v2-sum(ps%deriv*ps%r)) / &
    !   (product(ps%str%box_edges)*ndim)
    ! so if periodic boundary condition in use then put for now 
    props%pressure%val = 0.0

  end subroutine md_cal_properties
  
  subroutine md_cal_properties_extra(ps, props, list, pressure_comp, pot_energy)
    type (phasespace), intent(in) :: ps
    type (md_properties), intent(inout) :: props
    type (pe_list), intent(in) :: list
 		real (db), intent(in) :: pressure_comp, pot_energy    
    
    real(db) :: sum_mass_v2
    integer :: n_atoms
    
    n_atoms = size(ps%str%atoms)
    
    ! here assumes that ps%v2 is calculated at t=t+h/2
    ! and therefore needs to be re-calculated so that
    ! we calculate both the potential and kinetic energy
    ! at t=t+h
    
    sum_mass_v2 = sum(ps%p * ps%p * ps%inv_mass)
    
    ! unit of kinetic energy below is KJ per mol per atom (assuming
    ! mass_unit=amu=g/mol, t_unit=10^-13s, r_unit=AA), so do get the total
    ! energy of the system in 'absolute' units you would multiply this
    ! number by n_atoms and divide by Avogadro's number
    
    props%kin_energy%val = 0.5*sum_mass_v2 / n_atoms
    
    props%tot_energy%val = pot_energy / n_atoms + props%kin_energy%val
    
    
    ! pressure is in units of 16387.72 atm
    
    props%pressure%val = (sum_mass_v2+pressure_comp) / &
       (product(ps%str%box_edges)*ndim)

  end subroutine md_cal_properties_extra

  
  subroutine md_accum_properties(props)
    type (md_properties), intent(inout):: props    

    props%kin_energy%sum = props%kin_energy%sum + props%kin_energy%val
    props%tot_energy%sum = props%tot_energy%sum + props%tot_energy%val
    props%pressure%sum = props%pressure%sum + props%pressure%val
    
    props%kin_energy%sum2 = props%kin_energy%sum2 + & 
      props%kin_energy%val*props%kin_energy%val
    props%tot_energy%sum2 = props%tot_energy%sum2 + & 
      props%tot_energy%val*props%tot_energy%val  
    props%pressure%sum2 = props%pressure%sum2 + & 
      props%pressure%val*props%pressure%val
      
    props%n_accum = props%n_accum + 1
      
  end subroutine md_accum_properties


  subroutine md_reset_properties(props)
    type (md_properties), intent(out) :: props
    
    props%kin_energy%sum = 0.0
    props%tot_energy%sum = 0.0
    props%pressure%sum = 0.0
    
    props%kin_energy%sum2 = 0.0
    props%tot_energy%sum2 = 0.0  
    props%pressure%sum2 = 0.0
    
    props%n_accum = 0   

  end subroutine md_reset_properties
  
  
  subroutine md_print_properties(props)
    type (md_properties), intent(inout) :: props   

    props%kin_energy%ave = props%kin_energy%sum / props%n_accum
    props%tot_energy%ave = props%tot_energy%sum / props%n_accum
    props%pressure%ave = props%pressure%sum / props%n_accum
    
    props%kin_energy%esd = sqrt(max(props%kin_energy%sum2/props%n_accum - &
      props%kin_energy%ave*props%kin_energy%ave, 0.0))
    props%tot_energy%esd = sqrt(max(props%tot_energy%sum2/props%n_accum - &
      props%tot_energy%ave*props%tot_energy%ave, 0.0))
    props%pressure%esd = sqrt(max(props%pressure%sum2/props%n_accum - &
      props%pressure%ave*props%pressure%ave, 0.0)) 
    
    write(*, *) " "
    write(*, '(a,i6)') "Average over this many moves = ", props%n_accum
    write(*,'(a,2f12.6)') "Etot = ", props%tot_energy%ave, props%tot_energy%esd    
    write(*,'(a,2f12.6)') "Ekin = ", props%kin_energy%ave, props%kin_energy%esd
    write(*,'(a,2f12.6)') "P = ", props%pressure%ave, props%pressure%esd
    !write(*,'(a,2f12.6)') "Epot = ", props%tot_energy%ave - props%kin_energy%ave
    

  end subroutine md_print_properties
  
end module md_properties_class
