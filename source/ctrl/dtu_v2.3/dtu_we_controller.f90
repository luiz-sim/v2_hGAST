module dtu_we_controller_mod
   !
   ! Main module of the Basic DTU Wind Energy Controller. Interface for HAWC2.
   !
   use dtu_we_controller_fcns
   use turbine_controller_mod
   use safety_system_mod
   implicit none
   integer , save :: CtrlStatus
   real(mk) dump_array(50)
   real(mk), save :: time_old
contains
!**************************************************************************************************
subroutine init_regulation(array1, array2) !j bind(c, name='init_regulation')
   !
   ! Controller parameters initialization.
   !
!j !DEC$ IF .NOT. DEFINED(__MAKEFILE__)
!j !DEC$ ATTRIBUTES DLLEXPORT :: init_regulation
!j !DEC$ END IF
   real(mk), dimension(100), intent(inout) :: array1
   real(mk), dimension(1)  , intent(inout) :: array2
   ! Local vars
   integer i, ifejl
   character(len=32) text32
   real(mk) minimum_pitch_angle
   logical findes
   write(*, *) 'Basic DTU Wind Energy Controller ' //trim(adjustl(vertext32))// ' loaded...'
   ! Input array1 must contain
   ! Overall parameters
   !  constant   1 ; Rated power [kW]
   !  constant   2 ; Minimum rotor (LSS) speed [rad/s]
   !  constant   3 ; Rated rotor (LSS) speed [rad/s]
   !  constant   4 ; Maximum allowable generator torque [Nm]
   !  constant   5 ; Minimum pitch angle, PitchMin [deg],
   !               ; if |PitchMin|>90, then a table of <wsp,PitchMin> is read ;
   !               ; from a file named 'wptable.n', where n=int(PitchMin)
   !  constant   6 ; Maximum pitch angle [deg]
   !  constant   7 ; Maximum pitch velocity operation [deg/s]
   !  constant   8 ; Frequency of generator speed filter [Hz]
   !  constant   9 ; Damping ratio of speed filter [-]
   !  constant  10 ; Frequency of free-free DT torsion mode [Hz], if zero no notch filter used
   ! Partial load control parameters
   !  constant  11 ; Optimal Cp tracking K factor [Nm/(rad/s)^2], ;
   !               ; Qg=K*Omega^2, K=eta*0.5*rho*A*Cp_opt*R^3/lambda_opt^3
   !  constant  12 ; Proportional gain of torque controller [Nm/(rad/s)]
   !  constant  13 ; Integral gain of torque controller [Nm/rad]
   !  constant  14 ; Differential gain of torque controller [Nm/(rad/s^2)]
   ! Full load control parameters
   !  constant  15 ; Generator control switch [1=constant power, 0=constant torque, or interpolation between the two]
   !  constant  16 ; Proportional gain of pitch controller [rad/(rad/s)]
   !  constant  17 ; Integral gain of pitch controller [rad/rad]
   !  constant  18 ; Differential gain of pitch controller [rad/(rad/s^2)]
   !  constant  19 ; Proportional power error gain [rad/W]
   !  constant  20 ; Integral power error gain [rad/(Ws)]
   !  constant  21 ; Coefficient of linear term in aerodynamic gain scheduling, KK1 [deg]
   !  constant  22 ; Coefficient of quadratic term in aerodynamic gain scheduling, KK2 [deg^2] &
   !               ; (if zero, KK1 = pitch angle at double gain)
   !  constant  23 ; Normalized speed where the pitch controller gains are doubled [-]
   ! Cut-in simulation parameters
   !  constant  24 ; Cut-in time [s], no cut-in is simulated if zero or negative
   !  constant  25 ; Time delay for soft start of torque [1/1P]
   ! Cut-out simulation parameters
   !  constant  26 ; Shut-down time [s], no shut-down is simulated if zero or negative
   !  constant  27 ; Time of linear torque cut-out during a generator assisted stop [s]
   !  constant  28 ; Stop type [1=normal, 2=emergency]
   !  constant  29 ; Time delay for pitch stop after shut-down signal [s]
   !  constant  30 ; Maximum pitch velocity during initial period of stop [deg/s]
   !  constant  31 ; Time of initial pitch stop phase [s] (maintains pitch speed specified in constant 30)
   !  constant  32 ; Maximum pitch velocity during final phase of stop [deg/s]
   ! Expert parameters (keep default values unless otherwise given)
   !  constant  33 ; Time for the maximum torque rate = Maximum allowable generator torque/(constant 33 + 0.01s) [s]
   !  constant  34 ; Angle above lowest minimum pitch angle for switch to full load [deg]
   !  constant  35 ; Percentage of the rated speed when the torque limits are fully opened [%]
   !  constant  36 ; Time constant of 1st order filter on wind speed used for minimum pitch [1/1P]
   !  constant  37 ; Time constant of 1st order filter on pitch angle used for gain scheduling [1/1P]
   ! Drivetrain damper
   !  constant  38 ; Proportional gain of active DT damper [Nm/(rad/s)], requires frequency in input 10
   ! Overspeed
   !  constant  39 ; Overspeed percentage before initiating turbine controller alarm (shut-down) [%]
   ! Additional non-linear pitch control term (not used when all zero)
   !  constant  40 ; Rotor speed error scaling factor        Err0       [rad/s]
   !  constant  41 ; Rotor acceleration error scaling factor ErrDot0    [rad/s^2]
   !  constant  42 ; Pitch rate gain                         PitNonLin1 [rad/s]
   ! Storm control command
   !  constant  43 ; Wind speed 'Vstorm' above which derating of rotor speed is used [m/s]
   !  constant  44 ; Cut-out wind speed (only used for derating of rotor speed in storm) [m/s]
   ! Safety system parameters
   !  constant  45 ; Overspeed percentage before initiating safety system alarm (shut-down) [%]
   !  constant  46 ; Max low-pass filtered tower top acceleration level before initiating turbine controller alarm (shut-down) [m/s^2]
   ! Turbine parameter
   !  constant  47 ; Nominal rotor diameter [m]
   ! Parameters for rotor inertia reduction in variable speed region
   !  constant  48 ; Proportional gain on rotor acceleration in variable speed region [Nm/(rad/s^2)] (not used when zero)
   ! Parameters for alternative partial load controller with PI regulated TSR tracking
   !  constant  49 ; Optimal tip speed ratio [-] (only used when K=constant 11 = 0 otherwise  Qg=K*Omega^2 is used)
   ! Parameters for adding aerodynamic drivetrain damping on gain scheduling
   !  constant  50 ; Aerodynamic DT damping coefficient at the operational point of zero pitch angle [Nm/(rad/s)] (not used when zero)
   !  constant  51 ; Coefficient of linear term in aerodynamic DT damping scheduling, KK1 [deg]
   !  constant  52 ; Coefficient of quadratic term in aerodynamic DT damping scheduling, KK2 [deg^2]
   !
   ! Output array2 contains nothing for init
   !
!j---- Added by GL
   PID_gen_var%INIT_OUT= array1( 99) !gentrq
   PID_pit_var%INIT_OUT= array1(100) !collpit
   PitchColRefOld      = array1(100) !collpit !j
   ! Overall parameters
   PeRated             = array1( 1)*1000.0_mk
   GenSpeedRefMin      = array1( 2)
   GenSpeedRefMax      = array1( 3)
   GenTorqueMax        = array1( 4)
   minimum_pitch_angle = array1( 5)*degrad
   PitchStopAng        = array1( 6)*degrad
   PID_pit_var%velmax  = array1( 7)*degrad
   ! Generator speed second order low pass filter
   omega2ordervar%f0   = array1( 8)
   omega2ordervar%zeta = array1( 9)
   power2ordervar%f0   = array1( 8)
   power2ordervar%zeta = array1( 9)
   MoniVar%omega2ordervar%f0           = array1(8)
   MoniVar%omega2ordervar%zeta         = array1(9)
   SafetySystemVar%omega2ordervar%f0   = array1(8)
   SafetySystemVar%omega2ordervar%zeta = array1(9)
   ! Drivetrain mode notch filters for pitch controller
       DT_mode_filt%f0 = array1(10)
   pwr_DT_mode_filt%f0 = array1(10)
   ! Partial load control parameters
   Kopt             = array1(11)
   PID_gen_var%Kpro = array1(12)
   PID_gen_var%Kint = array1(13)
   PID_gen_var%Kdif = array1(14)
   ! Full load control parameters
   TorqueCtrlRatio=min(1.0_mk,max(0.0_mk,array1(15)))
   ! - Gains
   PID_pit_var%kpro(1) = array1(16)
   PID_pit_var%kint(1) = array1(17)
   PID_pit_var%kdif(1) = array1(18)
   PID_pit_var%kpro(2) = array1(19)
   PID_pit_var%kint(2) = array1(20)
   PID_pit_var%kdif(2) = 0.0_mk
   ! - Gain-scheduling
   PitchGSVar%invkk1 = 1.0_mk/(array1(21)*degrad)
   if (array1(22).eq.0.0_mk) then
     PitchGSVar%invkk2 = 0.0_mk
   else
     PitchGSVar%invkk2 = 1.0_mk/(array1(22)*degrad*degrad)
   endif
   rel_limit = array1(23)
   ! Cut-in simulation parameters
   CutinVar%time  = array1(24)
   CutinVar%delay = array1(25)*2.0_mk*pi/GenSpeedRefMax
   ! Cut-out simulation parameters
   CutoutVar%time         = array1(26)
   CutoutVar%torquedelay  = array1(27)
   CutoutVar%stoptype = int(array1(28))
   CutoutVar%pitchdelay   = array1(29)
   CutoutVar%pitchvelmax  = array1(30)*degrad
   CutoutVar%pitchdelay2  = array1(31)
   CutoutVar%pitchvelmax2 = array1(32)*degrad
   if (CutinVar%time .gt. 0.0_mk)then
     CtrlStatus = -2
     generator_cutin=.false.
   endif
   ! Expert parameters (keep default values unless otherwise given)
   PID_gen_var%velmax       = GenTorqueMax/(array1(33)+0.001_mk) !was 0. in ver2.0
   SwitchVar%pitang_upper   = array1(34)*degrad
   SwitchVar%rel_sp_open_Qg = array1(35)*0.01_mk
   wspfirstordervar%tau     = array1(36)*2.0_mk*pi/GenSpeedRefMax
   pitchfirstordervar%tau   = array1(37)*2.0_mk*pi/GenSpeedRefMax
   ! Drivetrain damper
   DT_damper%gain           = array1(38)
   DT_damper%bandpass%f0    = array1(10)
   ! Overspeed
   MoniVar%overspeed = (1.0_mk + array1(39)*0.01_mk)*GenSpeedRefMax
   ! Additional non-linear pitch control term
   Err0       = array1(40)
   ErrDot0    = array1(41)
   PitNonLin1 = array1(42)
  ! Default and derived parameters
   GenTorqueRated = PeRated/GenSpeedRefMax
   MoniVar%rystevagtfirstordervar%tau = 2.0_mk*pi/GenSpeedRefMax
   SafetySystemVar%rystevagtfirstordervar%tau = 2.0_mk*pi/GenSpeedRefMax
   ! Wind speed table
   if (dabs(minimum_pitch_angle).lt.90.0_mk*degrad) then
     OPdatavar%lines=2
     OPdatavar%wpdata(1,1) = 0.0_mk
     OPdatavar%wpdata(2,1) = 99.0_mk
     OPdatavar%wpdata(1,2) = minimum_pitch_angle
     OPdatavar%wpdata(2,2) = minimum_pitch_angle
   else
     write(text32,'(i3)') int(minimum_pitch_angle*raddeg)
     inquire(file='./wpdata.'//trim(adjustl(text32)),exist=findes)
     if (findes) then
       open(88,file='./wpdata.'//trim(adjustl(text32)))
       read(88,*,iostat=ifejl) OPdatavar%lines
       if (ifejl.eq.0) then
         do i=1,OPdatavar%lines
           read(88,*,iostat=ifejl) OPdatavar%wpdata(i,1),OPdatavar%wpdata(i,2)
           if (ifejl.ne.0) then
             write(*,*) ' *** ERROR *** Could not read lines in minimum '&
                      //'pitch table in file wpdata.'//trim(adjustl(text32))
             stop
           endif
           OPdatavar%wpdata(i,2)=OPdatavar%wpdata(i,2)*degrad
         enddo
       else
         write(*,*) ' *** ERROR *** Could not read number of lines '&
                  //'in minimum pitch table in file wpdata.'//trim(adjustl(text32))
         stop
       endif
       close(88)
     else
       write(*,*) ' *** ERROR *** File ''wpdata.'//trim(adjustl(text32))&
                //''' does not exist in the folder'
       stop
     endif
   endif
   ! Storm controller input
   Vstorm  = array1(43) ! [m/s] Vstorm (e.g. 25)
   Vcutout = array1(44) ! [m/s] Vcut-out (e.g. 45)
   if (Vcutout.gt.Vstorm) then
     write (*,'(a,f4.1,a,f4.1,a)') ' Storm control is active above ', Vstorm, &
                                   'm/s until cut-out at ', Vcutout, 'm/s'
   endif
   ! Overspeed monitor
   SafetySystemVar%overspeed = (1.0_mk + array1(45)*0.01_mk)*GenSpeedRefMax
   ! "Rystevagt" monitor
   MoniVar%RysteVagtLevel = array1(46)
   R = 0.5_mk*array1(47)
   ! Alternative partial load controller
   Kopt_dot= array1(48)
   TSR_opt = array1(49)
   if (array1(11) .le. 0.0_mk) then
     PartialLoadControlMode = 2 !Kopt<=0
   else
     PartialLoadControlMode = 1 !Kopt> 0
   endif
   ! Gain scheduling dQdomega
   PitchGSVar%kp_speed = array1(50)
   if (array1(51) .gt. 0.0_mk) then
     PitchGSVar%invkk1_speed = 1.0_mk/(array1(51)*degrad)
   else
     PitchGSVar%invkk1_speed = 0.0_mk
   endif
   if (array1(52) .gt. 0.0_mk) then
     PitchGSVar%invkk2_speed = 1.0_mk/(array1(52)*degrad*degrad)
   else
     PitchGSVar%invkk2_speed = 0.0_mk
   endif
   ! Set parameters that can be modified with advanced options
   ! -Generator torque exclusion zone [disabled]
   ExcluZone%Lwr             = 0.0_mk                             !array1(53)
   ExcluZone%Lwr_Tg          = 0.0_mk                             !array1(54)
   ExcluZone%Hwr             = 0.0_mk                             !array1(55)
   ExcluZone%Hwr_Tg          = 0.0_mk                             !array1(56)
   ExcluZone%time_excl_delay = 0.0_mk                             !array1(57)
   ! -Drive train mode damper [disabled if gain array1(38)=0.]
!j DT_damper%gain            = array1(38)
!j DT_damper%bandpass%f0     = array1(10)
   DT_damper%notch%f0        = 10.0_mk*array1(10)                 !array1(58)
   DT_damper%bandpass%zeta   = 0.02_mk                            !array1(59)
   DT_damper%notch%zeta2     = 0.01_mk                            !array1(60)
   DT_damper%Td              = 0.0_mk                             !array1(61)
   ! -Tower top fore-aft mode damper [disabled if gain=0.]
   TTfa_damper%bandpass%f0   = 10.0_mk                            !array1(62)
   TTfa_damper%notch%f0      = 10.0_mk                            !array1(63)
   TTfa_damper%bandpass%zeta = 0.02_mk                            !array1(64)
   TTfa_damper%notch%zeta2   = 0.01_mk                            !array1(65)
   TTfa_damper%gain          = 0.0_mk                             !array1(66)
   TTfa_damper%Td            = 0.0_mk                             !array1(67)
   TTfa_PWRfirstordervar%tau = 10.0_mk                            !array1(68)
   TTfa_PWR_lower            = 0.0_mk                             !array1(69)
   TTfa_PWR_upper            = 0.0_mk                             !array1(70)
   ! -Tower top side-to- mode filter [disabled if exclusion zone is disabled]
   ExcluZone%notch%f0        = 100.0_mk                           !array1(71)
   ExcluZone%notch%zeta2     = 0.01_mk                            !array1(72)
   ! -"Rystevagt" monitor for Safety System
   SafetySystemVar%RysteVagtLevel = MoniVar%RysteVagtLevel*1.1_mk !array1(73)
   ! Gear Ratio
   GearRatio                 = 1.0_mk
   ! Pitch deviation monitor
   DeltaPitchThreshold       = 0.0_mk
   TAve_Pitch                = 0.0_mk
   ! Initiate the dynamic variables
   stepno                    = 0
   time_old                  = 0.0_mk
!qq--- Added by GL
   time_old                  =-1.0_mk
   AddedPitchRate            = 0.0_mk
   PitchAngles               = 0.0_mk
   AveragedMeanPitchAngles   = 0.0_mk
   AveragedPitchReference    = 0.0_mk
   ! No output
   array2                    = 0.0_mk
   return
end subroutine init_regulation
!**************************************************************************************************
subroutine init_regulation_advanced(array1, array2) !j bind(c,name='init_regulation_advanced')
   !
   ! Controller parameters initialization with additional inputs.
   !
!j !DEC$ IF .NOT. DEFINED(__MAKEFILE__)
!j !DEC$ ATTRIBUTES DLLEXPORT::init_regulation_advanced
!j !DEC$ END IF
   real(mk), dimension(100), intent(inout) ::  array1
   real(mk), dimension(1)  , intent(inout) ::  array2
   ! Torque exclusion zone
   !  constant  53 ; Exclusion zone: Lower speed limit [rad/s]                                                                       (Default 0 used if zero)
   !  constant  54 ; Exclusion zone: Generator torque at lower limit [Nm]                                                            (Default 0 used if zero)
   !  constant  55 ; Exclusion zone: Upper speed limit [rad/s] (if =< 0 then exclusion zone functionality is inactive)
   !  constant  56 ; Exclusion zone: Generator torque at upper limit [Nm]                                                            (Default 0 used if zero) 
   !  constant  57 ; Time constant of reference switching at exclusion zone [s]                                                      (Default 0 used if zero)
   ! DT torsion mode damper
   !  constant  58 ; Frequency of notch filter [Hz]                                                                                  (Default 10 x input 10 used if zero)
   !  constant  59 ; Damping of BP filter [-]                                                                                        (Default 0.02 used if zero) 
   !  constant  60 ; Damping of notch filter [-]                                                                                     (Default 0.01 used if zero) 
   !  constant  61 ; Phase lag of damper [s] =>  max 40*dt                                                                           (Default 0 used if zero) 
   ! Fore-aft Tower mode damper                                                                                                      
   !  constant  62 ; Frequency of BP filter [Hz]                                                                                     (Default 10 used if zero)\\ 
   !  constant  63 ; Frequency of notch fiter [Hz]                                                                                   (Default 10 used if zero)\\ 
   !  constant  64 ; Damping of BP filter [-]                                                                                        (Default 0.02 used if zero)\\
   !  constant  65 ; Damping of notch filter [-]                                                                                     (Default 0.01 used if zero)\\
   !  constant  66 ; Gain of damper [-]                                                                                              (Default 0 used if zero)\\ 
   !  constant  67 ; Phase lag of damper [s] =>  max 40*dt                                                                           (Default 0 used if zero)\\ 
   !  constant  68 ; Time constant of 1st order filter on PWR used for fore-aft Tower mode damper GS [Hz]                            (Default 10 used if zero)
   !  constant  69 ; Lower PWR limit used for fore-aft Tower mode damper GS [-]                                                      (Default 0 used if zero)
   !  constant  70 ; Upper PWR limit used for fore-aft Tower mode damper GS [-]                                                      (Default 0 used if zero) 
   ! Side-to-side Tower mode filter
   !  constant  71 ; Frequency of Tower side-to-sede notch filter [Hz]                                                               (Default 100 used if zero)
   !  constant  72 ; Damping of notch filter [-]                                                                                     (Default 0.01 used if zero)
   ! Safety parameters
   !  constant  73 ; Max low-pass filtered tower top acceleration level before initiating safety system alarm (shut-down) [m/s^2]    (Default 1.1 x input 46 used if zero)
   !  constant  74 ; Time constant of 1st order filter on tower top acceleration [1/1P]                                              (Default 1 used if zero)
   ! Pitch deviation monitor parameters
   !  constant  75 ; Parameters for pitch deviation monitoring. The format is 1,nnn,mmm 
   !               ; where 'nnn' [s] is the period of the moving average and 'mmm' is threshold of the deviation [0.1 deg] (functionality is inactive if value $<$ 1,000,000)
   ! Gear ratio
   !  constant  76 ; Gear ratio used for the calculation of the LSS rotational speeds and the HSS generator torque reference [-] (Default 1 if zero)
   !
   call init_regulation(array1, array2)
   ! Generator torque exclusion zone
   if (array1(53).gt.0.0_mk) ExcluZone%Lwr                      = array1(53)
   if (array1(54).gt.0.0_mk) ExcluZone%Lwr_Tg                   = array1(54)
   if (array1(55).gt.0.0_mk) ExcluZone%Hwr                      = array1(55)
   if (array1(56).gt.0.0_mk) ExcluZone%Hwr_Tg                   = array1(56)
   if (array1(57).gt.0.0_mk) ExcluZone%time_excl_delay          = array1(57)
   ! Drive train mode damper
   if (array1(58).gt.0.0_mk) DT_damper%notch%f0                 = array1(58)
   if (array1(59).gt.0.0_mk) DT_damper%bandpass%zeta            = array1(59)
   if (array1(60).gt.0.0_mk) DT_damper%notch%zeta2              = array1(60)
   if (array1(61).gt.0.0_mk) DT_damper%Td                       = array1(61)
   ! Tower top fore-aft mode damper
   if (array1(62).gt.0.0_mk) TTfa_damper%bandpass%f0            = array1(62)
   if (array1(63).gt.0.0_mk) TTfa_damper%notch%f0               = array1(63)
   if (array1(64).gt.0.0_mk) TTfa_damper%bandpass%zeta          = array1(64)
   if (array1(65).gt.0.0_mk) TTfa_damper%notch%zeta2            = array1(65)
   if (array1(66).gt.0.0_mk) TTfa_damper%gain                   = array1(66)
   if (array1(67).gt.0.0_mk) TTfa_damper%Td                     = array1(67)
   if (array1(68).gt.0.0_mk) TTfa_PWRfirstordervar%tau          = 1.0_mk/(2.0_mk*pi*array1(68))
   if (array1(69).gt.0.0_mk) TTfa_PWR_lower                     = array1(69)
   if (array1(70).gt.0.0_mk) TTfa_PWR_upper                     = array1(70)
   !Tower top side-to- mode filter
   if (array1(71).gt.0.0_mk) ExcluZone%notch%f0                 = array1(71)
   if (array1(72).gt.0.0_mk) ExcluZone%notch%zeta2              = array1(72)
   ! "Rystevagt" monitor for Safety System
   if (array1(73).gt.0.0_mk) SafetySystemVar%RysteVagtLevel     = array1(73)
   if (array1(74).gt.0.0_mk) MoniVar%rystevagtfirstordervar%tau = array1(74)*2.0_mk*pi/GenSpeedRefMax !qq v2.0 array1(75)
   SafetySystemVar%rystevagtfirstordervar%tau                   = array1(74)*2.0_mk*pi/GenSpeedRefMax !qq v2.0 array1(75)
   ! Pitch devaiation monitor
   if (array1(75).gt.1000000.0_mk) then
     DeltaPitchThreshold = (array1(75)-floor(array1(75)/1000.0_mk)*1000.0_mk)*0.1_mk
     TAve_Pitch = (array1(75)-1000000.0_mk-DeltaPitchThreshold*10.0_mk)/1000.0_mk
   endif
   ! Gear ratio
   if (array1(76).gt.0.0_mk) GearRatio = array1(76)
   ! Initialization
   TimerExcl = -0.02_mk
end subroutine init_regulation_advanced
!**************************************************************************************************
subroutine update_regulation(array1, array2) !j bind(c,name='update_regulation')
   !
   ! Controller interface.
   !  - sets DLL inputs/outputs.
   !  - sets controller timers.
   !  - calls the safety system monitor (higher level).
   !
!j !DEC$ IF .NOT. DEFINED(__MAKEFILE__)
!j !DEC$ ATTRIBUTES DLLEXPORT :: update_regulation
!j !DEC$ END IF
   real(mk), dimension(100), intent(inout) :: array1
   real(mk), dimension(100), intent(inout) :: array2
   ! Input array1 must contain
   !
   !    1: general time                            [s]
   !    2: constraint bearing1 shaft_rot 1 only 2  [rad/s] Generator speed (Default LSS, if HSS insert gear ratio in input #76)
   !    3: constraint bearing2 pitch1 1 only 1     [rad]
   !    4: constraint bearing2 pitch2 1 only 1     [rad]
   !    5: constraint bearing2 pitch3 1 only 1     [rad]
   !  6-8: wind free_wind 1 0.0 0.0 hub height     [m/s] global coords at hub height
   !    9: elec. power  ; [W]
   !   10: grid flag  ; [1=no grid,0=grid]
   !   11: Tower top x-acceleration  ; [m/s^2]
   !   12: Tower top y-acceleration  ; [m/s^2]
   !
   ! Output array2 contains
   !
   !    1: Generator torque reference               [Nm]
   !    2: Pitch angle reference of blade 1         [rad]
   !    3: Pitch angle reference of blade 2         [rad]
   !    4: Pitch angle reference of blade 3         [rad]
   !    5: Power reference                          [W]
   !    6: Filtered wind speed                      [m/s]
   !    7: Filtered rotor speed                     [rad/s]
   !    8: Filtered rotor speed error for torque    [rad/s]
   !    9: Bandpass filtered rotor speed            [rad/s]
   !   10: Proportional term of torque contr.       [Nm]
   !   11: Integral term of torque controller       [Nm]
   !   12: Minimum limit of torque                  [Nm]
   !   13: Maximum limit of torque                  [Nm]
   !   14: Torque limit switch based on pitch       [-]
   !   15: Filtered rotor speed error for pitch     [rad/s]
   !   16: Power error for pitch                    [W]
   !   17: Proportional term of pitch controller    [rad]
   !   18: Integral term of pitch controller        [rad]
   !   19: Minimum limit of pitch                   [rad]
   !   20: Maximum limit of pitch                   [rad]
   !   21: Torque reference from DT damper          [Nm]
   !   22: Status signal                            [-]
   !   23: Total added pitch rate                   [rad/s]
   !   24: Filtered pitch angle                     [rad]
   !   25: Flag for mechnical brake                 [0=off/1=on]
   !   26: Flag for emergency pitch stop            [0=off/1=on]
   !   27: LP filtered acceleration level           [m/s^2]
   !   28: Rotor speed exlusion zone region         [-]
   !   29: Filtered tower top acc. for tower damper [m/s^2]
   !   30: Reference pitch from tower damper        [rad]
   !   31: Monitored average of reference pitch     [rad]
   !   32: Monitored ave. of pitch (largest devia.) [rad]
   !
   ! Local variables
   integer GridFlag, EmergPitchStop, ActiveMechBrake
   real(mk) GenSpeed, wsp, PitchVect(3), Pe, TT_acc(2), time
   real(mk) GenTorqueRef, PitchColRef
   EmergPitchStop = 0
   ActiveMechBrake = 0
   ! Time
   time = array1(1)
   !***********************************************************************************************
   ! Increment time step (may actually not be necessary in type2 DLLs)
   !***********************************************************************************************
   if (time .gt. time_old) then
     deltat = time - time_old
!qqv2if (deltat .eq. 0.0_mk) deltat = 0.01_mk
     time_old = time
     stepno = stepno + 1
   endif
   ! Rotor (Generator) speed in LSS
   GenSpeed = array1(2)/GearRatio
   ! Pitch angle
   PitchVect(1) = array1(3)
   PitchVect(2) = array1(4)
   PitchVect(3) = array1(5)
   ! Wind speed as horizontal vector sum
   wsp = dsqrt(array1(6)**2 + array1(7)**2)
   if (stepno.eq.1) then
      Pe = 0.0_mk ! Elec. power
      GridFlag = 0 ! Grid flag
   else
      Pe=array1(9) ! Elec. power
      GridFlag = int(array1(10)) ! Grid flag
   endif
   ! Tower top acceleration
   TT_acc(1) = array1(11)
   TT_acc(2) = array1(12)
   !***********************************************************************************************
   ! Safety system
   !***********************************************************************************************
   if (time .gt. 5.0_mk) then
      call safety_system(stepno, deltat, GenSpeed, TT_acc, EmergPitchStop, ActiveMechBrake, &
                         dump_array)
   endif
   !***********************************************************************************************
   ! Start-up timer monitoring
   !***********************************************************************************************
   if ((CutinVar%time .gt. 0.0_mk).and.(time .gt. CutinVar%time) .and. (CtrlStatus .eq. -2)) then
     CtrlStatus = -1
     TimerStartup = deltat
   endif
   !***********************************************************************************************
   ! Shut-down timer monitoring
   !***********************************************************************************************
   if ((CutoutVar%time .gt. 0.0_mk) .and. (time .gt. CutoutVar%time) .and. (CtrlStatus .eq. 0)) then
     if (stoptype.eq.1) CtrlStatus = 4
     if (stoptype.eq.2) CtrlStatus = 5
     GenSpeed_at_stop = GenSpeed
     GenTorque_at_stop = GenTorqueRefOld
     stoptype = CutoutVar%stoptype
     TimerShutdown = 0.0_mk
     TimerShutdown2 = 0.0_mk
   endif
   !***********************************************************************************************
   ! Wind turbine controller
   !***********************************************************************************************
!qq
   dump_array(50)=time
   call turbine_controller(CtrlStatus, GridFlag, GenSpeed, PitchVect, wsp, Pe, TT_acc, &
                           GenTorqueRef, PitchColRef, dump_array)
   !***********************************************************************************************
   ! Output
   !***********************************************************************************************
                                       !    1: time [sec]
   array2( 1) = GenTorqueRef/GearRatio !    2: Generator torque reference               [Nm]
   array2( 2) = PitchColRef            !    3: Pitch angle reference of blade 1         [rad]
   array2( 3) = PitchColRef            !    4: Pitch angle reference of blade 2         [rad]
   array2( 4) = PitchColRef            !    5: Pitch angle reference of blade 3         [rad]         ! Write into dump array
   array2( 5) = dump_array(1)          !    6: Power reference                          [W]           dump_array(1)  = GenTorqueRef*GenSpeed
   array2( 6) = dump_array(2)          !    7: Filtered wind speed                      [m/s]         dump_array(2)  = WSPfilt
   array2( 7) = dump_array(3)          !    8: Filtered rotor speed                     [rad/s]       dump_array(3)  = GenSpeedFilt
   array2( 8) = dump_array(4)          !    9: Filtered rotor speed error for torque    [rad/s]       dump_array(4)  = GenSpeedFiltErr
   array2( 9) = dump_array(5)          !   10: Bandpass filtered rotor speed            [rad/s]   ### NOT COMMUNICATED! ### 
   array2(10) = dump_array(6)          !   11: Proportional term of torque contr.       [Nm]          dump_array(6)  = PID_gen_var%outpro
   array2(11) = dump_array(7)          !   12: Integral term of torque controller       [Nm]          dump_array(7)  = PID_gen_var%outset
   array2(12) = dump_array(8)          !   13: Minimum limit of torque                  [Nm]          dump_array(8)  = PID_gen_var%outmin
   array2(13) = dump_array(9)          !   14: Maximum limit of torque                  [Nm]          dump_array(9)  = PID_gen_var%outmax
   array2(14) = dump_array(10)         !   15: Torque limit switch based on pitch       [-]           dump_array(10) = switch
   array2(15) = dump_array(11)         !   16: Filtered rotor speed error for pitch     [rad/s]       dump_array(11) = GenSpeedFiltErr
   array2(16) = dump_array(12)         !   17: Power error for pitch                    [W]           dump_array(12) = err_pitch(2)
   array2(17) = dump_array(13)         !   18: Proportional term of pitch controller    [rad]         dump_array(13) = PID_pit_var%outpro
   array2(18) = dump_array(14)         !   19: Integral term of pitch controller        [rad]         dump_array(14) = PID_pit_var%outset
   array2(19) = dump_array(15)         !   20: Minimum limit of pitch                   [rad]         dump_array(15) = PID_pit_var%outmin
   array2(20) = dump_array(16)         !   21: Maximum limit of pitch                   [rad]         dump_array(16) = PID_pit_var%outmax
   array2(21) = dump_array(17)         !   22: Torque reference from DT damper          [Nm]          dump_array(17) = Qdamp_ref
   array2(22) = dump_array(18)         !   23: Status signal                            [-]           dump_array(18) = CtrlStatus
   array2(23) = dump_array(19)         !   24: Total added pitch rate                   [rad/s]       dump_array(19) = AddedPitchRate
   array2(24) = dump_array(20)         !   25: Filtered pitch angle                     [rad]         dump_array(20) = PitchMeanFilt
   array2(25) = ActiveMechBrake        !   26: Flag for mechnical brake                 [0=off/1=on]
   array2(26) = EmergPitchStop         !   27: Flag for emergency pitch stop            [0=off/1=on]
   array2(27) = dump_array(23)         !   28: LP filtered acceleration level           [m/s^2]       dump_array(23) = TTAccFilt
   array2(28) = dump_array(24)         !   29: Rotor speed exlusion zone region         [-]           dump_array(24) = w_region
   array2(29) = dump_array(25)         !   30: Filtered tower top acc. for tower damper [m/s^2]       dump_array(25) = TTfa_acc_filt
   array2(30) = dump_array(26)         !   31: Reference pitch from tower damper        [rad]         dump_array(26) = theta_dam_ref
   array2(31) = dump_array(27)         !   32: Monitored average of reference pitch     [rad]         dump_array(27) = AveragedPitchReference(1)
   array2(32) = dump_array(28)         !   33: Monitored ave. of actual pitch (blade 1) [rad]         dump_array(28) = AveragedMeanPitchAngles(1)
!j  Added by jim
   array2(33) = dump_array(40)         !   34: GenTorque[Min,Max]_full    [Nm]
   array2(34) = dump_array(43)         !   35: GenTorqueMin_partial [Nm]
   array2(35) = dump_array(44)         !   36: GenTorqueMax_partial [Nm]
   array2(36) = dump_array(41)         !   37: switch wmin          [-]
   array2(37) = dump_array(42)         !   38: switch wmax          [-]
   return
end subroutine update_regulation
!**************************************************************************************************
end module dtu_we_controller_mod
