!modifications for built-in controller
!=========================================
!
!1. torque limits
!   1. if me omega oxi me torque
!   2. const power
!2. Speed Demand selector
!
!3. Drive train Damper
!   1. notch filter (6P??)
!   2. 2 frequencies
!
!4. Gain schedule
!
!5. PID different implementation [checks inside for P and I indipendently]
!   proportional [limits rate]
!   integrator [limits value]
!   2nd gain of PI2 is power and not torque!
!   set point variation in PI2 only through integral term so Qtop only in Ki
!
!6. filter spd error instead of speed [maybe lead to different gains?]
!
!
!OMEGA_ref=OMEAGA_max, allthough 2 variables are used in pdf!
!
! cut-in
! cut-out
! blade stuck
! emergency shut down
! brake
! electrical / mechanical losses
!
! min pitch
! derate (storm control) ?reduce omega set?
! stall regulated
!
! library
! DTD
! IPC/IFC
! TTa
! nonlinear gain
! nonlinear added pitch term
!Use Cbeam, only : DT_el, TPERIOD_el, PI2
!Use Cbeam, only : TIME, TREF, R2D, PI2, ICMOD, RAT_GEAR
!
! Init Torque? why zero?
!
!
!----------------------------------------------------------------------
 Subroutine control_init_el
!----------------------------------------------------------------------

 Use Cbeam
 Use Ctrl_mod

   implicit none


!--- Initialization of control variables
   QTC_el  = 0.d0;   QTCP_el = 0.d0;
   QTC1_el = 0.d0;   QTCP1_el= 0.d0;
   QTC2_el = 0.d0;   QTCP2_el= 0.d0;

   if (NBLADE_el == 0) return

   TGenLSS_el        =  0.d0 !TREF * RAT_GEAR !T_HSS = T_LSS/RAT_GEAR
   TGenLSS_D_el      =  0.d0
   TGenLSS_M_el      =  0.d0
   TMLossLSS_el      =  0.d0
   TBrakeLSS_el      =  0.d0
   TIME_BRAKE        = -1.d0

   GenTrqLSS_accel   =  0.d0 !generator  torque      for load reduction (ttop ax)
   BetaF_accel       =  0.d0 !collective flap  angle for load reduction (ttop ax)
   PITCH_accel       =  0.d0 !collective pitch angle for load reduction (ttop ax)

   PITCH_INDIV   (:) = 0.d0;
   PITCH_INDIV_P (:) = 0.d0;
   PITCH_INDIV_PP(:) = 0.d0;
   FLAP_INDIV    (:) = 0.d0;
   FLAP_INDIV_P  (:) = 0.d0;
   FLAP_INDIV_PP (:) = 0.d0;


!------ sellect the controller to be used
      ICMOD_def      = 0       !0: controller is loaded as dll or external subroutine (i.e. oc3)
   if (ICMOD==10) then
      ICMOD_def      = 1       !1: controller is the built-in controller with user-defined subroutines/parameters
      ICMOD          = 1
   endif
!0. OC3
!1. DLL
!2. DTU_v1
!3. DTU_v2

#ifndef CONTROLLER_TYPE
#define CONTROLLER_TYPE 0
#endif

#if   CONTROLLER_TYPE == 0
 write(10,*) 'OC3 Baseline controller is enabled'
#elif CONTROLLER_TYPE == 1
 write(10,*) 'DLL call for controller is enabled'
#elif CONTROLLER_TYPE == 2
 write(10,*) 'DTU_v1 controller is enabled'
#elif CONTROLLER_TYPE == 3
 write(10,*) 'DTU_v2 controller is enabled'
#elif CONTROLLER_TYPE == 33
 write(10,*) 'DTU_v2.3 controller is enabled'
#endif


 END Subroutine control_init_el
!--------------------------------------------------------------------------------
!
!  Subroutine : control   ------------------
!
!  Equation for the controls
!
!--------------------------------------------------------------------------------
 Subroutine control (NTIME, it, ISTEP)
!--------------------------------------------------------------------------------

 Use Ctrl_mod
 Use Cbeam, only : ICMOD, ICMOD_def, OMEGAR, TIME, TREF, UT1_el, NDFBT_el, NQSW !, R2D, PI2, RAT_GEAR

   implicit none

   integer, intent(in) :: NTIME, it, ISTEP

   real(8),save        :: Time_startP=-1.d0, Time_start=-1.d0
   real(8)             :: duration, c_start, c_down

!------------------------------------------------------------------------------------------------------------------------
!---------- Control modes ----------
! 0. no action
! 1. typical variable speed/pitch controller using external controller [from dll or fortran sub-routine] --> ICMOD_def=0
! 2. Parked simulation
! 3. Stall regulated WT [only Generator Torque, and/or tip brake]
! 4. Stall regulated WT [only Generator Torque, and/or tip brake] Airbrake activation
! 5. prescribed omega, pitch from file perform.inp
!10. typical variable speed/pitch controller using the built-in controller                               --> ICMOD_def=1
!------------------------------------------------------------------------------------------------------------------------

!------ Flap Controller (before main call of Control for IPC)
!!!qq call FLAP_Control     (it,ISTEP)


   if     (ICMOD == 0) then
!qq   call pitch_change
!qq   call tip_brake
!------ Controller is deactivated. [constant omega, pitch]
      call pitch_heli
      return
   elseif (ICMOD == 2) then
!------ Parked case. Controller is deactivated. [free-free b.c. for shaft, fixed pitch]
      TGenLSS_el = 0.d0
!qq   call tip_brake
      return
   elseif (ICMOD == 3.or.ICMOD == 4) then
!------ Stall regulated case. Controller is activated. [Generator Torque, constant pitch, tip brake only for specific dlcs]
!     T = B(w-ws)*(ws/ws) = B*slip*ws -> B=T/(slip*ws):: T(w), slip(w)=(w-ws)/ws at rated conditions, Trater=Prated/wrated

!----- update
     if (it==1) &
      Time_startP= Time_start
      Time_start = Time_startP

     if (-UT1_el(NDFBT_el+NQSW) < OMEGAR *0.97d0) then !for start-up
      TGenLSS_el   =  0.d0
      TGenLSS_D_el =  0.d0
     else
        if (Time_start<0.d0) Time_start = TIME

         duration= 0.5d0

         c_start = min(1.d0,     (TIME-Time_start)/duration)
         c_down  = max(0.d0,1.d0-(TIME-Time_start)/duration)

      TGenLSS_el   = c_start*TREF * (-UT1_el(NDFBT_el+NQSW)-OMEGAR)     ! B*(w-ws), w is negative, slip~=0.07 - 0.1
      TGenLSS_D_el = c_start*TREF
     endif
!!                             1      2     3       4              5              6                                7
!!    write(354,'(20e15.5)') TIME,c_start,c_down, Time_start, Time_startP,-UT1_el(NDFBT_el+NQSW)-OMEGAR*0.97d0, TGenLSS_el
      call tip_brake
      return
   elseif (ICMOD == 5) then
      TGenLSS_el = 0.d0
      call set_omega_pitch (NTIME,it)
      return
   endif


!--- Selection of controller [ICMOD_def 0: external controller, 1: built-in contoller]
   if    ( ICMOD_def == 0 ) then
      call controller (NTIME,it,ISTEP)
   elseif( ICMOD_def == 1 ) then
      call ctrl_main  (NTIME,it,ISTEP)
   endif


 END Subroutine control
!----------------------------------------------------------------------
 Subroutine ctrl_main (NTIME,it,ISTEP)
!----------------------------------------------------------------------

 Use Cbeam, only : TIME, TREF, R2D, PI2, ICMOD, RAT_GEAR
 Use Ctrl_mod

   implicit none

   integer, intent(in) :: NTIME, it, ISTEP
!  real(8), save :: array1(100), array2(100)


   goto (1,2,3,4), ISTEP

!--- main part
 1 call ctrl_sensors

!--- Initialization of the generic controller
   if ((NTIME==1).and.(it==1)) &
   call ctrl_init

!--- Main controller call, pass the communication variables
!  array1( 1) = TIME            !  1: general time                            [s]
!  array1( 2) = GenSpeed        !  2: constraint bearing1 shaft_rot 1 only 2  [rad/s] Generator speed (Default LSS, if HSS insert gear ratio in input #76)
!  array1( 3) = BlPitch(1)      !  3: constraint bearing2 pitch1 1 only 1     [rad]
!  array1( 4) = BlPitch(2)      !  4: constraint bearing2 pitch2 1 only 1     [rad]
!  array1( 5) = BlPitch(3)      !  5: constraint bearing2 pitch3 1 only 1     [rad]
!  array1( 6) = WindVel(1)      !  6: wind free_wind 1 0.0 0.0 hub height     [m/s] global coords at hub height
!  array1( 7) = WindVel(2)      !  7: wind free_wind 1 0.0 0.0 hub height     [m/s] global coords at hub height
!  array1( 8) = WindVel(3)      !  8: wind free_wind 1 0.0 0.0 hub height     [m/s] global coords at hub height
!  array1( 9) = array2 (5)      !  9: mech. power  ; [W]
!  array1(10) = 0               ! 10: grid flag  ; [1=no grid,0=grid]
!qq local or global??
!  array1(11) = ttop_acc(1)     ! 11: Tower top x-acceleration  ; [m/s^2]                                                                                   
!  array1(12) = ttop_acc(2)     ! 12: Tower top y-acceleration  ; [m/s^2]                                                                                   


!--- basic power control
   call power_control

!--- IPC
   call IPC_control  (ISTEP)

!--- TTACC (actuator, coll pitch, gen, coll flap)
   call TTACC_control(ISTEP)

!--- Drive train damper, Generator lag and Pitch & Flap Actuators
   call DTD_GenLag_PitFlapAct

   return

!--- Update ctrl vars & output
 2 continue

   QTCP_el (1:NQSC0) = QTC_el (1:NQSC0);
   QTCP1_el(1:NQSC0) = QTC1_el(1:NQSC0);
   QTCP2_el(1:NQSC0) = QTC2_el(1:NQSC0);

   call ctrl_output

   return

!--- Rerun
 3 return

!--- Recall
 4 return


 END Subroutine ctrl_main
!----------------------------------------------------------------------
 Subroutine ctrl_init
!----------------------------------------------------------------------

 Use Ctrl_mod
 Use Cbeam, only : TIME, TREF, R2D, PI2, ICMOD, RAT_GEAR

   implicit none

   real(8) :: array1(100)
   integer :: i, ii, iostat
   real(8) :: kk1, kk2


   open(88, file='controls.inp')
   do i = 1, 74
      ii=(i-1)*(i-9)*(i-13)*(i-22)*(i-34)*(i-41)*(i-52)*(i-62)*(i-63)*(i-70)
      if (ii==0) read(88, *) !** comments...
      read(88, *, iostat=iostat) array1(i)
      if (iostat /= 0) then; write(*,*) ' *** ERROR *** Could not read variable ', i,' in controller input file'; stop; endif;
   enddo !i
   close(88)

!------ Define parameters
   open(23,file='control_init.dat')
   do i=1,74
      write(23,*) array1(i)
   enddo
   close(23)

!------ power controller
      power_rated     = array1( 1) * 1000.d0
      i               = int(array1(2))
      if (i==1) then
         gearctrl     = 1.d0
      else
         gearctrl     = RAT_GEAR
      endif
      omega_min       = array1( 3)
      omega_max       = array1( 4)
      vs_Qmax         = array1( 5)
      vs_Kopt         = array1( 9)
      vs_Kp           = array1(10)
      vs_Ki           = array1(11)
      vs_Kd           = array1(12)
      vp_pitmin       = array1( 6)/R2D
      vp_pitmax       = array1( 7)/R2D
      vp_pitratmax    = array1( 8)/R2D
      vp_pitratmin    =-vp_pitratmax
      vp_ipowtor      = int(array1(13)) !constant  power (1) or torque (2)
      vp_Kp           = array1(14)
      vp_Ki           = array1(15)
      vp_Kd           = array1(16)
      vp_Kp2          = array1(17)
      vp_Ki2          = array1(18)
!------ GS
      kk1             = array1(19)/R2D
      kk2             = array1(20)/R2D
      GS_NL_lim       = array1(21)
!------ Filters
      LP2_freq  (1:2) = array1(22)*PI2  ! [rad/s]
      LP2_damp  (1:2) = array1(23)
      Notch_freq(1)   = array1(24)*PI2  ! [rad/s]
      Notch_damp(1)   = array1(25)
      Notch_freq(2)   = array1(26)*PI2  ! [rad/s]
      Notch_damp(2)   = array1(27)
      Notch_freq(3)   = array1(28)*PI2  ! [rad/s]
      Notch_damp(3)   = array1(29)
      LP1_lag   (1)   = array1(30)      ! [s] omega
      LP1_lag   (2)   = array1(31)      ! [s] pitch
      LP1_lag   (3)   = array1(32)      ! [s] power
      LP1_lag   (4)   = array1(33)      ! [s] torque
!------ Drive Train Damper
      DTD_gain(1)     = array1(34)       !          !400.d0*RAT_GEAR**2
      DTD_freq(1)     = array1(35)*PI2   ! [rad/s]  !  1.870d0 !Hz free-free
      DTD_damp(1)     = array1(36)       ! [-]      !  0.100d0
      DTD_gain(2)     = array1(37)       !          ! 00.d0
      DTD_freq(2)     = array1(38)*PI2   ! [rad/s]  ! 10.000d0 !Hz
      DTD_damp(2)     = array1(39)       ! [-]      !  0.100d0
      MaxGenTrqDTD    = array1(40)       ! [Nm]
!------ inputs for IPC (Blade root out-of-plane momemnts)
      ipc_nhT         = int(array1(41))  ! Number of harmonic (0:no, 1:1P, 2:2P) [-]
      ipc_PitchGain   = array1(42)       ! GAIN_PITCH         (default=1)        [-]
      ipc_FlapGain    = array1(43)       ! GAIN_FLAP          (default=1)        [-]
      GS_IPC0         = array1(44)       ! 0.01d0
      ipc_Ki (1:2)    = array1(45:46)    !
      ipc_Kp (1:2)    = array1(47:48)    ! 2.5d0*Ki(:)
      ipc_MaxPitchAmp = array1(49)/R2D   ! Maximum Pitch Amplitude               [rad]
      ipc_MaxFlapAmp  = array1(50)/R2D   ! Maximum Flap  Amplitude               [rad]
      ipc_ntime_init  = int(array1(51))  ! time step to start IPC ntime_init     [-]
!------ inputs for TTa
      tta_GenGain     = array1(52)       !                                       [-]
      tta_PitchGain   = array1(53)       !                                       [-]
      tta_FlapGain    = array1(54)       !                                       [-]
      GS_TTa0         = array1(55)       ! 0.01d0
      tta_Ki          = array1(56)       !
      tta_Kp          = array1(57)       ! 2.5d0*Ki(:)
      tta_MaxGenAmp   = array1(58)       ! Maximum Gen   Amplitude               [Nm]
      tta_MaxPitchAmp = array1(59)/R2D   ! Maximum Pitch Amplitude               [rad]
      tta_MaxFlapAmp  = array1(60)/R2D   ! Maximum Flap  Amplitude               [rad]
      tta_ntime_init  = int(array1(61))  ! time step to start TTa ntime_init     [-]
!------ inputs for Generator lag
      Gen_lag         = array1(62)       ! [sec]    !  0.02d0         !generator lag                          [sec]
!------ inputs for pitch Actuator
      itypePit_act    = int(array1(63))
      freqPit_act     = array1(64)*PI2   ! [rad/s] #!  1.6d0*PI2      !frequency         of pitch actuator    [rad/s] ### 1  Hz
      dampPit_act     = array1(65)       ! [-]     #!  0.8d0          !damping factor    of pitch actuator    [-]     ### 0.7
      MinPitAng_act   = array1(66)/R2D   ! [rad]    ! -5.d0 / R2D     !lower angle limit of pitch actuator    [rad]
      MaxPitAng_act   = array1(67)/R2D   ! [rad]    ! 90.d0 / R2D     !upper angle limit of pitch actuator    [rad]
      MaxPitRat_act   = array1(68)/R2D   ! [rad/s]  ! 10.d0 / R2D     !speed       limit of pitch actuator    [rad/s]   ### 8
      MaxPitAcc_act   = array1(69)/R2D   ! [rad/s^2]! 15.d0 / R2D     !accelerationlimit of pitch actuator    [rad/s^2] ### 8
!------ inputs for flap  Actuator
      itypeFlap_act   = int(array1(70))
      timeFlap_act    = array1(71)       ! [rad/s] #!  1.6d0*PI2      !frequency         of flap  actuator    [rad/s] ### 1  Hz
      MinFlapAng_act  = array1(72)/R2D   ! [rad]    ! -5.d0 / R2D     !lower angle limit of flap  actuator    [rad]
      MaxFlapAng_act  = array1(73)/R2D   ! [rad]    ! 90.d0 / R2D     !upper angle limit of flap  actuator    [rad]
      MaxFlapRat_act  = array1(74)/R2D   ! [rad/s]  ! 10.d0 / R2D     !speed       limit of flap  actuator    [rad/s]   ### 8

      if (itypePit_act  /=2) then; write(*,*)'only Pitch actuator type 2 is supported atm'; stop; endif
      if (itypeFlap_act /=1) then; write(*,*)'only Flap  actuator type 1 is supported atm'; stop; endif

      Q_rated         = power_rated / omega_max
!------ Gain Schedule kk1_inv,kk2_inv
      if    (kk1>0.d0 .and. kk2>0.d0) then
         kk1_inv = 1.d0/kk1
         kk2_inv = 1.d0/kk2
      elseif(kk1>0.d0               ) then
         kk1_inv = 1.d0/kk1
         kk2_inv = 0.d0
      else
         kk1_inv = 0.d0
         kk2_inv = 0.d0
      endif

!qqqq
!     vp_Kp           = vp_Kp + vp_Kp2 * Q_rated
!     vp_Ki           = vp_Ki + vp_Ki2 * Q_rated
!     vp_Kp2          = 0.d0
!     vp_Ki2          = 0.d0
!--------------------------------------------
!--- Initialization of control variables                                  current configuration
!---    1. DTD1                                               [QTC ]      free-free 1st
!---    2. DTD2                                               [QTC ]      free-free 2nd
!---    3. GEN LAG                                            [QTC1]      ~0.02s
!--- 4- 6. Pitch Actuator  2nd order                          [QTC ]      ~1.6Hz   , z=0.8
!--- 7-10. LOWPASS-1st (genspeed,collpitch,power,Torque)      [QTC1]      Tperiod
!---11-14. LOWPASS-2nd (genspeed,collpitch,power,Torque)      [QTC ]      -//-
!
!---15-18. IPC filters LP2 (Mtilt1P, Myaw1P, Mtilt2P, Myaw2P) [QTC]       2.5P     , z=0.7
!---19-22. IPC Notch       (Mtilt1P, Myaw1P, Mtilt2P, Myaw2P) [QTC]       3.0P     , z=0.2
!---23-26. IPC PI          (Mtilt1P, Myaw1P, Mtilt2P, Myaw2P) [QTC1]
!---   27. Ttop accel  LP2                                    [QTC ]      2.5P     , z=0.7
!---   28. Ttop accel  Notch                                  [QTC ]      3.0P     , z=0.2
!---   29. Ttop accel  PI                                     [QTC1]
!---30-32. Flap Actuator  1st order                           [QTC1]
!
!---   33. PI   Torque                                        [QTC1]
!---   34. PI2  pitchcol                                      [QTC1]
!---   35. LP2  omega err                                     [QTC ]
!---   36. Not2 omega err                                     [QTC ]
!---   37. Not2 omega err                                     [QTC ]
!---   38. Not2 omega err                                     [QTC ]
!---   39. LP2  power err                                     [QTC ]
!--------------------------------------------
   call ctrl_sensors

      GenTrqDem = max (min( vs_Kopt * GenSpeed**2, Q_rated ), 0.d0)
write(*,*)'GenTrqDem - init', GenTrqDem
      QTC_el  = 0.d0;
      QTC1_el = 0.d0;
      QTC2_el = 0.d0;

!DTD  QTC_el  (1:2) = 0.d0                    !BP2
      QTC1_el (3  ) = GenTrqDem               !LP1
      QTC_el  (4:6) = BlPitch(1:3)            !Actuator
      QTC1_el (7  ) = GenSpeed                !LP1
      QTC1_el (8  ) = collpit                 !LP1
      QTC1_el (9  ) = GenTrqDem * GenSpeed    !LP1
      QTC1_el (10 ) = GenTrqDem               !LP1
!!    QTC1_el (10 ) = WindVel(1)              !LP1
      QTC_el  (11 ) = GenSpeed                !LP2
      QTC_el  (12 ) = collpit                 !notch
      QTC_el  (13 ) = GenTrqDem * GenSpeed    !LP2
      QTC_el  (14 ) = GenTrqDem               !notch
      QTC1_el (33 ) = GenTrqDem               !PI   Torque
      QTC1_el (34 ) = collpit                 !PID2 pitch

      QTCP_el       = QTC_el ;   
      QTCP1_el      = QTC1_el;   
      QTCP2_el      = QTC2_el;   


 END Subroutine ctrl_init
!----------------------------------------------------------------------
 Subroutine ctrl_output
!----------------------------------------------------------------------

 Use Ctrl_mod
 Use Cbeam, only : TIME

   implicit none

   integer :: i, ISTEP

   ISTEP = 2 !output ctrl vars


#  if   ASCII == 1
      open(27,file='control.dat',access='append')
       write(27,100)              &
#  elif ASCII == 0
      open(27,file='control.bin',access='append',form='UNFORMATTED')
       write(27    )              &
#  endif
         sngl(TIME                  ),& ! 1
         sngl(vs_PI_Qmin            ),& ! 2
         sngl(vs_PI_Qmax            ),& ! 3
         sngl(vs_Qopt               ),& ! 4
         sngl(vs_Qoptmin*vs_Reg1_lim),& ! 5
         sngl(omega_set             ),& ! 6 
         sngl(vp_PI2_pitmin         ),& ! 7
         sngl(GS                    ),& ! 8
         sngl(GS_NL                 ),& ! 9
         sngl(GS_IPC                ),& !10
         sngl(GS_TTa                ),& !11
         sngl(GenTrqDTDtot*gearctrl ),& !12
         sngl(GenTrqDTD(1)*gearctrl ),& !13
         sngl(GenTrqDTD(2)*gearctrl ),& !14
         sngl(GenSpeed_LP1          ),& !15
         sngl(GenSpeed_LP2          ),& !16
         sngl(collpit_LP1           ),& !17
         sngl(collpit_LP2           ),& !18
         sngl(power_LP1             ),& !19
         sngl(power_LP2             ),& !20
         sngl(GenTrq_LP1            ),& !21
         sngl(GenTrq_LP2            )   !22
      close(27)


#  if   ASCII == 1
      open(27,file='controller_qtc0.dat',access='append')
       write(27,100)              &
#  elif ASCII == 0
      open(27,file='controller_qtc0.bin',access='append',form='UNFORMATTED')
       write(27    )              &
#  endif
         sngl(TIME      )        ,&
        (sngl(QTC_el (i)),i=1,39)
      close(27)

#  if   ASCII == 1
      open(27,file='controller_qtc1.dat',access='append')
       write(27,100)               &
#  elif ASCII == 0
      open(27,file='controller_qtc1.bin',access='append',form='UNFORMATTED')
       write(27    )               &
#  endif
         sngl(TIME       )        ,&
        (sngl(QTC1_el (i)),i=1,39)
      close(27)

#  if   ASCII == 1
      open(27,file='controller_qtc2.dat',access='append')
       write(27,100)               &
#  elif ASCII == 0
      open(27,file='controller_qtc2.bin',access='append',form='UNFORMATTED')
       write(27    )               &
#  endif
         sngl(TIME       )        ,&
        (sngl(QTC2_el (i)),i=1,39)
      close(27)

!--- for Output
   call IPC_control  (ISTEP)
   call TTACC_control(ISTEP)

 100  format(1000e15.5)


 END Subroutine ctrl_output
!----------------------------------------------------------------------
 Subroutine ctrl_sensors
!----------------------------------------------------------------------

 Use Cbeam
 Use Craft, only : HUB_VEL
 Use Ctrl_mod

   implicit none

   integer       :: nbod,nbsub,nb_el,nel, NQACC
   real(8)       :: x(3),dx(3),ddx(3),dxl(3),ddxl(3), AKSI0(3)


!--- generator speed
      GenSpeed = -UT1_el(NDFBT_el+NQSW) * gearctrl
     dGenSpeed = -UT2_el(NDFBT_el+NQSW) * gearctrl

!--- pitch: b_hGAST= -b_COLL -b_INDIV + b_IMB ==> b_COLL= -b_hGAST -b_INDIV + b_IMB
      collpit  = 0.d0
   do nbod     = 1, NBLADE_el
      NQACC    = ACCU%NQSACC_el(nbod-1) 
      BlPitch(nbod) = - (UT_el (NDFBT_el+NQSP+NQACC+5)-PITCH_Imb_el(nbod) ) -PITCH_INDIV(nbod) -PITCH_accel
      collpit  = collpit + BlPitch(nbod) / dble(NBLADE_el)
   enddo

!--- Wind Velocity
   WindVel(1:3) = HUB_VEL (1:3)

!--- tower top acceleration
   if ( NBODBTWT_el == NTOWER_el) then
      nbod     = NBODBT_el
      nbsub    = body(nbod)%NBODSUB_el
      nb_el    = body(nbod)%NBODTGLB_el(nbsub  )
      nel      = subbody(nb_el)%NTEB_el
      AKSI0(:) = 0.d0;
      AKSI0(2) = subbody(nb_el)%HTA_el (nel+1)

      call beam_kinematics (nbod,nbsub,nel,AKSI0,x,dx,ddx,dxl,ddxl)
   else
      ddx  (:) = 0.d0;
   endif
      ttop_acc(1:3) = ddx(1:3) !global (x:fore-aft,y:side-side,z:vertical)

!!!!--- Tower base Acceleration
!!! 1 accel1           = 0.d0
!!!  if (ICASE_el >=3) then
!!!   nbod             = NBLADE_el + 2
!!!   nbsub            = 1
!!!   nb_el            = body(nbod)%NBODTGLB_el(nbsub)
!!!   nel              = 1
!!!   AKSI0(:)         = 0.d0;
!!!
!!!   call beam_kinematics (nbod,nbsub,nel,AKSI0,xfl,dxfl,ddxfl,dxlfl,ddxlfl)
!!!
!!!   accel1           = ddxlfl(1)
!!!  endif
!!!
!!!!--- Tower top Acceleration
!!!   nbod             = NBLADE_el + 2
!!!   nbsub            = body   (nbod)%NBODSUB_el
!!!   nb_el            = body   (nbod)%NBODTGLB_el(nbsub)
!!!   nel              = subbody(nb_el)%NTEB_el
!!!   AKSI0(:)         = 0.d0;
!!!   AKSI0(2)         = subbody(nb_el)%HTA_el  (nel+1)
!!!
!!!   call beam_kinematics (nbod,nbsub,nel,AKSI0,xfl,dxfl,ddxfl,dxlfl,ddxlfl)
!!!
!!!   accel2           = ddxlfl(1)-accel1 !Atop-Abot

!--- Blade root moments (flap / edge) -> pitch (out-of-plane / in-plane) -> psi (tilt / yaw) 1P, 2P
!  call CYCLIC_MOM ( M_NRT_C, psi, NH )
!  Mroot(3,2) = 


 END Subroutine ctrl_sensors
!-----------------------------------------------------------------------
!
! Provides for each blade:
!   M_BLD  : 3 root Moments wrt blade        c.s. (flap        /torsion/edge    ) [N]
!   M_RD   : 3 root Moments wrt rotor disk   c.s. (out-of-plane/torsion/in-plane) [N]
!   M_NR   : 3 root Moments wrt Non-rotating c.s. (tilt        /yaw    /roll    ) [N]
!** psi    : positive azimuth angle                                             [rad]
!   pitch  : positive pitch   angle                                             [rad]
!
!   M_NRT  : the total Moment wrt Non-rotating c.s.                               [N]
!** M_NRT_C: the total Tilting and Yawing moments for the NH harmonics only by
!            considering the Mout-of-plane moment                                 [N]
!
!** only these vars are communicated at the final version
!-----------------------------------------------------------------------
 Subroutine CYCLIC_MOM ( M_NRT_C, psi, NH)
!-----------------------------------------------------------------------

 Use Cbeam

   implicit none

   integer, intent(in   ) :: NH
   real(8), intent(  out) :: M_NRT_C (NH,2)
   real(8), intent(  out) :: psi  (NBLADE_el)

   real(8) :: M_BLD(NBLADE_el,3), M_RD (NBLADE_el,3), M_NR (NBLADE_el,3), M_NRT(3), pitch(NBLADE_el)
   integer :: nbod,nbsub,nb_el,nn_el,nel,nod, I1,I2, h
   real(8) :: AT0G(3,3), AT0L(3,3), AG(3,3), AL(3,3), Apsi(3,3)
   real(8) :: FG  (6)  , FL  (6)


      M_NRT  (    1:3) = 0.d0;
   do nbod             = 1, NBLADE_el
      nbsub            = 1
      nb_el            = body(nbod)%NBODTGLB_el(nbsub)
      nn_el            = body(nbod)%NNODTGLB_el(nbsub)
      nel              = 1
      nod              = 1

!----- wrt blade root after cone and pitch, before pre-curved angle
     !AT0G   (1:3,1:3) = ATP_el   (nbod,1:3,1:3) !Blade c.s. after cone, before pitch 
      AT0G   (1:3,1:3) = ATPWrG_el(nbod,1:3,1:3) !Blade c.s. after cone, before pitch 
      AT0L   (1:3,1:3) = ATPWr_el (nbod,1:3,1:3) !Blade c.s. after cone, and pitch, before pre-curved angles
      AG     (1:3,1:3) = matmul (AT0G (1:3,1:3),transf_mat(nn_el)%A_el(1:3,1:3))
      AL     (1:3,1:3) = matmul (AT0L (1:3,1:3),transf_mat(nn_el)%A_el(1:3,1:3))

      I1               = subbody(nb_el)%NDFPNACC_el ( nod-1      ) + 1
      I2               = subbody(nb_el)%NDFPNACC_el ( nod        )
      FG         (1:6) = subbody(nb_el)%AFLOC_el    ( nel, I1:I2 )
      FL         (1:6) = subbody(nb_el)%AFLOC_el    ( nel, I1:I2 )

!------ Rotate LOADS
      call LOADS_TRANS0_NEW ( FG, AG, 1, 1, 1 )
      call LOADS_TRANS0_NEW ( FL, AL, 1, 1, 1 )

      pitch (nbod    ) = -UT_el(NDFBT_el+NQSP+ACCU%NQSACC_el(nbod-1)+5)
      psi   (nbod    ) = -UT_el(NDFBT_el+NQSW) - UThub_el + PHI0(nbod)

      call ROT_MATRIX0      ( 3, -psi(nbod), Apsi )

      M_BLD (nbod,1:3) = FL (IMDOF_el(4:6))
      M_RD  (nbod,1:3) = FG (IMDOF_el(4:6))
      M_NR  (nbod,1:3) = matmul(Apsi (1:3,1:3), M_RD(nbod,1:3))
      M_NRT (     1:3) = M_NRT       (1:3)    + M_NR(nbod,1:3)
   enddo !nbod

      M_NRT_C(:,:)     = 0.d0;
   do nbod             = 1, NBLADE_el
   do h                = 1, NH
      M_NRT_C(h,1)     =  M_NRT_C(h,1) + M_RD(nbod,1) * dcos(dble(h)*psi(nbod)) !Mtilt
      M_NRT_C(h,2)     =  M_NRT_C(h,2) + M_RD(nbod,1) * dsin(dble(h)*psi(nbod)) !Myaw
   enddo !h
   enddo !nbod


 END Subroutine CYCLIC_MOM
!-----------------------------------------------------------------------------------------------
!--- Power controller ---
!
!--- Filters 1st and 2nd
!--- Gain Schedule         OK --------------------> pending NL #23
!qq- Define pitch min-----------------------------> pending
!--- Speed demand selector OK
!--- Torque limits         OK
!--- PID torque            OK
!--- PID pitch             OK --------------------> pending NL pitch differential term #40-42 additional PI gains if contribution >1
!
!----------------------------------------------------------------------
 Subroutine power_control
!----------------------------------------------------------------------

 Use Cbeam, only : DT_el, TPERIOD_el, PI2, R2D
 Use Ctrl_mod

   implicit none

   integer :: i,ieq,ieq1,ieq2
   real(8) :: VAR , DVAR , DDVAR
   real(8) :: VAR2, DVAR2
   real(8) :: W,W1,W2, D,D1,D2, TVLAG
!--- filters
   real(8) :: varf(4), timef(4), dummy, freqnot(3)


!--- 1st order lowpass filter for genspeed, collpitch, power and UwindX
      varf(1)       = GenSpeed             ; timef(1) = LP1_lag(1) !!TPERIOD_el/2.5 !ieq=7
      varf(2)       = collpit              ; timef(2) = LP1_lag(2) !!TPERIOD_el     !ieq=8
      varf(3)       = GenTrqDem * GenSpeed ; timef(3) = LP1_lag(3) !!TPERIOD_el     !ieq=9
      varf(4)       = GenTrqDem            ; timef(4) = LP1_lag(4) !!TPERIOD_el/2.5 !ieq=10
!!!   varf(4)       = WindVel(1)           ; timef(4) = TPERIOD_el     !ieq=10 !sqrt(WindVel(1)**2+WindVel(2)**2)
   do i             = 1, 4
      ieq           = 7 + i-1 !7-10
      VAR           = varf (i)
      TVLAG         = timef(i)
      
      call filter_lowpass1 ( TVLAG       , DT_el        , VAR          , &
                             QTCP_el(ieq), QTCP1_el(ieq), QTCP2_el(ieq), &
                             QTC_el (ieq), QTC1_el (ieq), QTC2_el (ieq), &
                             dummy                                         )
   enddo !i

!--- filtered vars LP1
   GenSpeed_LP1     = QTC1_el( 7) !vs_PI_Qmax = P/w, GS_NL
   collpit_LP1      = QTC1_el( 8) !GS
   power_LP1        = QTC1_el( 9) !GS_IPC
   GenTrq_LP1       = QTC1_el(10) !omega_set, torque limits
!! Ux_LP1           = QTC1_el(10)

!--- 2nd order lowpass filter for genspeed, collpitch, power and UwindX
   do i             = 1, 4
      ieq           = 11 + i-1 !11-14
      VAR           = varf (i)
      W             = PI2 / timef(i)  ![rad/s] !0.4Hz = 2.5P
      D             = 0.7d0

      call filter_lowpass2 ( W            , D             , DT_el         , &
                             VAR                                          , &
                             QTCP_el (ieq), QTCP1_el (ieq), QTCP2_el (ieq), &
                             QTC_el  (ieq), QTC1_el  (ieq),  QTC2_el (ieq)    )
   enddo !i

!--- filtered vars LP2
   GenSpeed_LP2     = QTC_el(11)
   collpit_LP2      = QTC_el(12)
   power_LP2        = QTC_el(13)
   GenTrq_LP2       = QTC_el(14)
!! Ux_LP2           = QTC_el(14)

!qq- Define pitch min LP Power

!--- Cp tracking
     vs_Qopt        = min (vs_Kopt * GenSpeed**2, vs_Qmax)
     vs_Qoptmin     =      vs_Kopt * omega_min**2

!--- Define omega set for Torque control
     vs_Reg1_lim    = 1.02d0
   if ( GenTrq_LP1 <= vs_Qoptmin * vs_Reg1_lim ) then
     omega_set      = omega_min
   else
     omega_set      = omega_max
   endif

!--- TORQUE_LIMITS, define vs_PI_Qmin, vs_PI_Qmax
   if ( GenTrq_LP1 > vs_Qoptmin * vs_Reg1_lim ) then

     vs_PI_Qmin     = vs_Qopt

    if (vp_ipowtor == 1) then 
!--- constant power
     vs_PI_Qmax     = min (power_rated / GenSpeed_LP1, vs_Qmax)
    elseif (vp_ipowtor == 2) then 
!--- constant Torque
     vs_PI_Qmax     = Q_rated
    endif

!----- Don't reduce torque until fine pitch is reached
!!qq vs_pit_switch  = 0.d0 !0.5d0 / R2D
!!   if (collpit - vp_pitmin > vs_pit_switch ) &
!!   vs_PI_Qmin     = min ( vs_PI_Qmax, GenTrqDem )
   else
     vs_PI_Qmin     = 0.d0
     vs_PI_Qmax     = vs_Qopt
   endif

!--- Gain Schedule
      GS            = 1.d0/(1.d0 + collpit_LP1 * kk1_inv + collpit_LP1**2 * kk2_inv)

!--- GS for IPC and TTa
   if     (power_LP1<=0.8d0*power_rated) then
      GS_IPC        = GS_IPC0
      GS_TTa        = GS_TTa0
   elseif (power_LP1>=      power_rated) then
      GS_IPC        = 1.d0
      GS_TTa        = 1.d0
   else
      GS_IPC        = (power_LP1-0.8d0*power_rated)/(0.2d0*power_rated) * (1.d0-GS_IPC0) + GS_IPC0
      GS_TTa        = (power_LP1-0.8d0*power_rated)/(0.2d0*power_rated) * (1.d0-GS_TTa0) + GS_TTa0
   endif
      GS_IPC        = GS_IPC * GS
      GS_TTa        = GS_TTa * GS

!-- Nonlinear gain to avoid large rotor speed excursion
      GS_NL         = 1.d0
   if(GS_NL_lim > 1.d0 ) &
      GS_NL         = 1.d0 + ( (GenSpeed_LP1-omega_max) / (omega_max*(GS_NL_lim - 1.d0)) )**2
      GS            = GS * GS_NL

!-----------------
!- Torque Control
!-----------------

!------ PI controller for generator torque
   ieq              = 33
   VAR              =  GenSpeed - omega_set
   DVAR             = dGenSpeed

   call PID2_CONTROL( 1.d0        , DT_el                       , & !Gain, DT
                      vs_Kp       , vs_Ki        , vs_Kd        , & !Kp ,Ki ,Kd
                      0.d0        , 0.d0         , 0.d0         , & !Kp2,Ki2,Kd2
                      -1.d30      , 1.d30                       , & !minvar  , maxvar
                      vs_PI_Qmin  , vs_PI_Qmax                  , & !mindvar , maxdvar
                      -1.d30      , 1.d30                       , & !minddvar, maxddvar
                      VAR         , DVAR         , 0.d0         , & !VAR , DVAR, DDVAR
                      0.d0        , 0.d0         , 0.d0         , & !VAR2, DVAR2,DDVAR2
                      QTCP_el(ieq), QTCP1_el(ieq), QTCP2_el(ieq), &
                      QTC_el (ieq), QTC1_el (ieq), QTC2_el (ieq)    )

   GenTrqDem        = QTC1_el(ieq)

!-----------------
!- Pitch Control
!-----------------

!--- omega err: 2nd-order low-pass filter
   ieq              = 35
   VAR              = GenSpeed - omega_max        !omega_set
   W                = LP2_freq(1) !!10.0d0 * (PI2/TPERIOD_el)    ![rad/s]
   D                = LP2_damp(1) !!0.7d0

   call filter_lowpass2 ( W            , D             , DT_el         , &
                          VAR                                          , &
                          QTCP_el (ieq), QTCP1_el (ieq), QTCP2_el (ieq), &
                          QTC_el  (ieq), QTC1_el  (ieq),  QTC2_el (ieq)    )
   
!--- omega err: 3x 2nd-order Notch filters at 1st free-free, 2nd free-free and 3P frequencies
      freqnot(1)    = freq_free1 !DTD_freq(2)
      freqnot(2)    = 3.0d0 * (PI2/TPERIOD_el)
      freqnot(3)    = DTD_freq(1)

   do i             = 1, 3
      ieq           = 36 + i-1 !36-38
      ieq1          = ieq-1
      VAR           = QTC_el (ieq1)
      DVAR          = QTC1_el(ieq1)
      DDVAR         = QTC2_el(ieq1)
      W1            = Notch_freq(i)    ![rad/s]
      W2            = W1
      D1            = 0.001d0
      D2            = Notch_damp(i)  !!0.200d0
      
      call filter_notch2 ( W1   , W2   , D1   , D2   , DT_el             , &
                           VAR  , DVAR , DDVAR                           , &
                           QTCP_el (ieq) , QTCP1_el (ieq), QTCP2_el (ieq), &
                           QTC_el  (ieq) , QTC1_el  (ieq),  QTC2_el (ieq)    )
   enddo !i

!--- Power err: 2nd-order low-pass filter
   ieq              = 39
   VAR              = GenTrqDem * GenSpeed - power_rated
   W                = LP2_freq(2) !!10.0d0 * (PI2/TPERIOD_el)    ![rad/s]
   D                = LP2_damp(2) !!0.7d0
   
   call filter_lowpass2 ( W            , D             , DT_el         , &
                          VAR                                          , &
                          QTCP_el (ieq), QTCP1_el (ieq), QTCP2_el (ieq), &
                          QTC_el  (ieq), QTC1_el  (ieq),  QTC2_el (ieq)    )

!--- pitch: PID2 controller
   ieq              = 34
   ieq1             = 38
   ieq2             = 39
!qq only LP2 and notch at free-free
   VAR              = QTC_el (ieq1) !  omega err
   DVAR             = QTC1_el(ieq1) ! domega err
   DDVAR            = QTC2_el(ieq1) !ddomega err
   VAR2             = QTC_el (ieq2) ! power  err
   DVAR2            = QTC1_el(ieq2) !dpower  err

   call PID2_CONTROL( GS          , DT_el                       , & !Gain, DT
                      vp_Kp       , vp_Ki        , vp_Kd        , & !Kp ,Ki ,Kd
                      vp_Kp2      , vp_Ki2       , 0.d0         , & !Kp2,Ki2,Kd2
                      -1.d30      , 1.d30                       , & !minvar  , maxvar
                      vp_pitmin   , vp_pitmax                   , & !mindvar , maxdvar
                      vp_pitratmin, vp_pitratmax                , & !minddvar, maxddvar
                      VAR         , DVAR         , DDVAR        , & !VAR , DVAR, DDVAR
                      VAR2        , DVAR2        , 0.d0         , & !VAR2, DVAR2,DDVAR2
                      QTCP_el(ieq), QTCP1_el(ieq), QTCP2_el(ieq), &
                      QTC_el (ieq), QTC1_el (ieq), QTC2_el (ieq)    )

   PitAngDem(1:3)   = QTC1_el(ieq);


 END Subroutine power_control
!----------------------------------------------------------------------
 Subroutine DTD_GenLag_PitFlapAct
!----------------------------------------------------------------------
       
 Use Cbeam, only : DT_el, NDFBT_el, NQSP, ACCU, NBLADE_el, UT_el,UT1_el,UT2_el, PITCH_Imb_el
 Use Ctrl_mod

   implicit none

   integer :: i,ieq,nbod
!--- Drive Train Damper
   real(8) :: QVGAIN, QVFREQ, QVDAMP
   real(8) :: Q0M(2)
!--- Gen lag
   real(8) :: Q0MLAG, TVLAG
!--- Actuators
   real(8) :: W,D
   real(8) :: VAR   , DVAR  !, DDVAR
   real(8) :: MinVAR, MaxVAR, MinDVAR, MaxDVAR, MinDDVAR, MaxDDVAR

!------------------------------------------------------------------------------------
!--- DTD and Generator time lag

!--- Drive Train Damper 1st and 2nd free-free modes
      GenTrqDTDtot = 0.d0
      GenTrqDTD(:) = 0.d0;
   do i            = 1, 2
      ieq          = i           !1-2
      QVGAIN       = DTD_gain(i)
      QVFREQ       = DTD_freq(i) ![rad/s]
      QVDAMP       = DTD_damp(i)
      VAR          = dGenSpeed

      call filter_bandpass2 ( QVGAIN      , QVFREQ       , QVDAMP       , &
                              DT_el       , VAR          , Q0M(i)       , &
                              QTCP_el(ieq), QTCP1_el(ieq), QTCP2_el(ieq), &
                              QTC_el (ieq), QTC1_el (ieq), QTC2_el (ieq)    )

      GenTrqDTDtot = GenTrqDTDtot + QTC_el(ieq)
      GenTrqDTD(i) = QTC_el(ieq)
   enddo !i
      GenTrqDTDtot = min ( max( GenTrqDTDtot, -MaxGenTrqDTD  ), MaxGenTrqDTD )


!--- Generator's time lag for torque demand
   ieq          = 3
   VAR          = GenTrqDem + GenTrqDTDtot + GenTrqLSS_accel
   TVLAG        = Gen_lag

   call filter_lowpass1 ( TVLAG       , DT_el        , VAR          , &
                          QTCP_el(ieq), QTCP1_el(ieq), QTCP2_el(ieq), &
                          QTC_el (ieq), QTC1_el (ieq), QTC2_el (ieq), &
                          Q0MLAG                                        )

   TGenLSS_el   =  QTC1_el (ieq)           * gearctrl
   TGenLSS_M_el = (Q0M(1)+Q0M(2)) * Q0MLAG * gearctrl**2
!!!TGenLSS_D_el = GenTrq_D_out             * gearctrl

!------------------------------------------------------------------------------------
!--- Actuators (Pitch and Flap)

!--- Pitch Actuator: 2nd order (Pitch Angle demand)
   do nbod      = 1, NBLADE_el
      i         = NDFBT_el + NQSP + ACCU%NQSACC_el(nbod-1) + 5
      ieq       = 4 + (nbod-1) !4-6
      W         = freqPit_act      ![rad/s]
      D         = dampPit_act
      VAR       = PitAngDem (nbod) + PITCH_INDIV(nbod) + PITCH_accel
      DVAR      = 0.d0 !PitVel_out(nbod)
      MinVAR    = MinPitAng_act
      MaxVAR    = MaxPitAng_act
      MinDVAR   =-MaxPitRat_act
      MaxDVAR   = MaxPitRat_act
      MinDDVAR  =-MaxPitAcc_act
      MaxDDVAR  = MaxPitAcc_act

      call ACTUATOR_2 ( W, D        , DT_el                       , &
                        MinVAR      , MaxVAR                      , &
                        MinDVAR     , MaxDVAR                     , &
                        MinDDVAR    , MaxDDVAR                    , &
                        VAR         , DVAR                        , &
                        QTCP_el(ieq), QTCP1_el(ieq), QTCP2_el(ieq), &
                        QTC_el (ieq), QTC1_el (ieq), QTC2_el (ieq)    )

      UT_el (i) = -QTC_el (ieq) + PITCH_Imb_el(nbod)
      UT1_el(i) = -QTC1_el(ieq)
      UT2_el(i) = -QTC2_el(ieq)
   enddo !nbod


!--- Flap Actuator: 1st order low pass filter (usually 0.1s)
   do nbod      = 1, NBLADE_el
      ieq       = 30 + (nbod-1) !30-32

      TVLAG     = timeFlap_act !sec
      VAR       = FLAP_INDIV (nbod) + BetaF_accel
      MinVAR    = MinFlapAng_act
      MaxVAR    = MaxFlapAng_act
      MinDVAR   =-MaxFlapRat_act 
      MaxDVAR   = MaxFlapRat_act 
      MinDDVAR  =-1.d30
      MaxDDVAR  = 1.d30

      call ACTUATOR_1 ( TVLAG   , DT_el          , &
                        MinVAR  , MaxVAR         , &
                        MinDVAR , MaxDVAR        , &
                        MinDDVAR, MaxDDVAR       , &
                        VAR                      , &
                        QTCP_el(ieq), QTCP1_el(ieq), QTCP2_el(ieq), &
                        QTC_el (ieq), QTC1_el (ieq), QTC2_el (ieq)    )

!!    FlapAng(nbod) =  QTC1_el(ieq)                      * R2D
!!    FlapVel(nbod) =  QTC2_el(ieq)                      * R2D
!!    FlapAcc(nbod) = (QTC2_el(ieq)-QTCP2_el(ieq))/DT_el * R2D
   enddo !nbod


 END Subroutine DTD_GenLag_PitFlapAct
!----------------------------------------------------------------------
 Subroutine IPC_control (ISTEP)
!----------------------------------------------------------------------
       
 Use Cbeam, only : TPERIOD_el, DT_el, PI2, TIME, NBLADE_el, R2D
 Use Ctrl_mod

   implicit none

   integer, intent(in) :: ISTEP
   integer, parameter  :: NH = 2 !Number of harmonics for flap motion
   integer, parameter  :: NB = 3 !Number of blades

   integer :: h, nn, nbod, ieq1, ieq2, ieq3
   real(8) :: GAIN   , CONP  , CONI
   real(8) :: VAR    , DVAR  , DDVAR
   real(8) :: W      , D     , W1,W2,D1,D2

   real(8)       :: theta   (NH,2), Bang(NBLADE_el)
   real(8), save :: M_NRT_C (NH,2)                  !saved for output
   real(8), save :: PVEL (NB), PACC(NB), psi (NB)
   real(8), save :: FVEL (NB), FACC(NB)


   if (NBLADE_el > NB) then; write(*,*)'increase NB param in IPC_control';stop;endif

   goto (1,2), ISTEP

 1 call CYCLIC_MOM ( M_NRT_C, psi, NH )

      theta(:,:)  = 0.d0
   do h           = 1, ipc_nhT
   do nn          = 1, 2 !1: Mtilt, 2: Myaw

!------ 2nd-order low-pass filter at 2.5P
      ieq1        = 15 + (h-1)*2 + nn-1 !15-18  Mtilt1P, Myaw1P, Mtilt2P, Myaw2P
      W           = PI2 / TPERIOD_el *2.5d0   ![rad/s]
      D           = 0.7d0
      VAR         = M_NRT_C(h,nn) !nn:1 Mtilt, nn:2 Myaw
      
      call filter_lowpass2 ( W             , D              , DT_el          , &
                             VAR                                             , &
                             QTCP_el (ieq1), QTCP1_el (ieq1), QTCP2_el (ieq1), &
                             QTC_el  (ieq1), QTC1_el  (ieq1),  QTC2_el (ieq1)    )
      
!------ 2nd-order Notch filter at 3P
      ieq2        = 19 + (h-1)*2 + nn-1 !19-22  Mtilt1P, Myaw1P, Mtilt2P, Myaw2P
      W1          = PI2 / TPERIOD_el *3.d0  ![rad/s]
      W2          = W1
      D1          = 0.001d0
      D2          = 0.200d0
      VAR         = QTC_el (ieq1)
      DVAR        = QTC1_el(ieq1)
      DDVAR       = QTC2_el(ieq1)
   
      call filter_notch2 ( W1   , W2   , D1   , D2   , DT_el                , &
                           VAR  , DVAR , DDVAR                              , &
                           QTCP_el (ieq2) , QTCP1_el (ieq2), QTCP2_el (ieq2), &
                           QTC_el  (ieq2) , QTC1_el  (ieq2),  QTC2_el (ieq2)    )

      if (TIME< dble(ipc_ntime_init)*DT_el) cycle

!------ PI controller for Non-rotating M[tilt,yaw] at 1P & 2P freq.
      ieq3        = 23 + (h-1)*2 + nn-1 !23-26  Mtilt1P, Myaw1P, Mtilt2P, Myaw2P
      GAIN        = GS_IPC
      CONP        = ipc_kp(h)
      CONI        = ipc_ki(h)
      VAR         = QTC_el (ieq2)
      DVAR        = QTC1_el(ieq2)

      call PID2_CONTROL( GAIN,     DT_el                              , & !Gain, DT
                         CONP    , CONI    , 0.d0                     , & !Kp ,Ki ,Kd
                         0.d0    , 0.d0    , 0.d0                     , & !Kp2,Ki2,Kd2
                         -1.d30  , 1.d30                              , & !minvar  , maxvar
                         -1.d30  , 1.d30                              , & !mindvar , maxdvar
                         -1.d30  , 1.d30                              , & !minddvar, maxddvar
                         VAR     , DVAR    , 0.d0                     , & !VAR , DVAR, DDVAR
                         0.d0    , 0.d0    , 0.d0                     , & !VAR2, DVAR2,DDVAR2
                         QTCP_el(ieq3), QTCP1_el(ieq3), QTCP2_el(ieq3), &
                         QTC_el (ieq3), QTC1_el (ieq3), QTC2_el (ieq3)    )

      theta(h,nn) = QTC1_el(ieq3)
   enddo !nn
   enddo !h

!--- Inverse Colleman - Set IPC, IFC angles
      Bang(:)     = 0.d0
   do nbod        = 1, NBLADE_el
   do h           = 1, ipc_nhT       !tilt                         yaw
      Bang(nbod)  = Bang(nbod) + theta(h,1)*dcos(h*psi(nbod)) + theta(h,2)*dsin(h*psi(nbod)) !inv Colleman
   enddo !h

!------ IPC - Cyclic (Individual) Pitch Control
      PITCH_INDIV_PP(nbod) = PITCH_INDIV_P(nbod)
      PITCH_INDIV_P (nbod) = PITCH_INDIV  (nbod)
      PITCH_INDIV   (nbod) = min( max( Bang(nbod)*ipc_PitchGain, -ipc_MaxPitchAmp ), ipc_MaxPitchAmp )
      PVEL          (nbod) = ( 3.d0*PITCH_INDIV(nbod) - 4.d0 * PITCH_INDIV_P(nbod) + PITCH_INDIV_PP(nbod)) / (2*DT_el)
      PACC          (nbod) = (      PITCH_INDIV(nbod) - 2.d0 * PITCH_INDIV_P(nbod) + PITCH_INDIV_PP(nbod)) / (  DT_el**2)

!------ IFC - Cyclic (Individual) Flap Control
      FLAP_INDIV_PP (nbod) = FLAP_INDIV_P(nbod)
      FLAP_INDIV_P  (nbod) = FLAP_INDIV  (nbod)
      FLAP_INDIV    (nbod) = min( max(-Bang(nbod)*ipc_FlapGain, -ipc_MaxFlapAmp ), ipc_MaxFlapAmp )
      FVEL          (nbod) = ( 3.d0*FLAP_INDIV(nbod) - 4.d0 * FLAP_INDIV_P(nbod) + FLAP_INDIV_PP(nbod)) / (2*DT_el)
      FACC          (nbod) = (      FLAP_INDIV(nbod) - 2.d0 * FLAP_INDIV_P(nbod) + FLAP_INDIV_PP(nbod)) / (  DT_el**2)
   enddo !nbod

   return

!--- Individual Pitch Angle, Velocity and Acceleration [deg,deg/s,deg/s^2]
 2 if (ipc_PitchGain > 1.d-10) then
#    if   ASCII == 1
      open(23,file='control_IPC.dat',access='append')
       write(23,'(1000e15.5)')                                 &
#    elif ASCII == 0
      open(23,file='control_IPC.bin',access='append',form='UNFORMATTED')
       write(23              )                                 &
#    endif
         sngl(TIME                 )                          ,& !1
        (sngl(psi        (nbod)*R2D)                          ,& !2, 6, 10
         sngl(PITCH_INDIV(nbod)*R2D)                          ,& !3, 7, 11
         sngl(PVEL       (nbod)*R2D)                          ,& !4, 8, 12
         sngl(PACC       (nbod)*R2D),    nbod = 1, NBLADE_el) ,& !5, 9, 13
       ((sngl(M_NRT_C    (h    ,nn)),    nn=1,2) , h=1,ipc_nhT)  !14: tilt1P, 15:yaw1P, 16: tilt2P, 17:yaw2P
      close(23)
   endif

   if (ipc_FlapGain > 1.d-10) then
#    if   ASCII == 1
      open(23,file='control_IFC.dat',access='append')
       write(23,'(1000e15.5)')                                 &
#    elif ASCII == 0
      open(23,file='control_IFC.bin',access='append',form='UNFORMATTED')
       write(23              )                                 &
#    endif
         sngl(TIME                 )                          ,& !1
        (sngl(psi        (nbod)*R2D)                          ,& !2, 6, 10
         sngl(FLAP_INDIV (nbod)*R2D)                          ,& !3, 7, 11
         sngl(FVEL       (nbod)*R2D)                          ,& !4, 8, 12
         sngl(FACC       (nbod)*R2D),    nbod = 1, NBLADE_el) ,& !5, 9, 13
       ((sngl(M_NRT_C    (h    ,nn)),    nn=1,2) , h=1,ipc_nhT)  !14: tilt1P, 15:yaw1P, 16: tilt2P, 17:yaw2P
      close(23)
   endif


 END Subroutine IPC_control
!----------------------------------------------------------------------
!
!  Tower top acceleration PI control through Gen and Coll Pitch or Flap
!
!----------------------------------------------------------------------
 Subroutine TTACC_control (ISTEP)
!----------------------------------------------------------------------
       
 Use Cbeam, only : TIME, DT_el, R2D, PI2, TPERIOD_el
 Use Cbeam, only : ICASE_el,NBLADE_el,body,subbody
 Use Ctrl_mod

   implicit none

   integer, intent(in) :: ISTEP

   real(8), save :: accel1, accel2, accelF !save for output
   integer :: nbod, nbsub, nb_el, nel, ieq1, ieq2, ieq3

   real (8) :: AKSI0(3), xfl(3), dxfl(3), ddxfl(3), dxlfl(3), ddxlfl(3)

   real(8)  :: W, D, W1,W2, D1,D2
   real(8)  :: VAR, DVAR, DDVAR
   real(8)  :: GAIN, CONP, CONI, COND


   goto (1,2), ISTEP

!--- Tower base Acceleration
 1 accel1           = 0.d0
  if (ICASE_el >=3) then
   nbod             = NBLADE_el + 2
   nbsub            = 1
   nb_el            = body(nbod)%NBODTGLB_el(nbsub)
   nel              = 1
   AKSI0(:)         = 0.d0;

   call beam_kinematics (nbod,nbsub,nel,AKSI0,xfl,dxfl,ddxfl,dxlfl,ddxlfl)

   accel1           = ddxlfl(1)
  endif

!--- Tower top Acceleration
   nbod             = NBLADE_el + 2
   nbsub            = body   (nbod)%NBODSUB_el
   nb_el            = body   (nbod)%NBODTGLB_el(nbsub)
   nel              = subbody(nb_el)%NTEB_el
   AKSI0(:)         = 0.d0;
   AKSI0(2)         = subbody(nb_el)%HTA_el  (nel+1)

   call beam_kinematics (nbod,nbsub,nel,AKSI0,xfl,dxfl,ddxfl,dxlfl,ddxlfl)

   accel2           = ddxlfl(1)-accel1 !Atop-Abot

!--- Tower Acceleration: 2nd-order low-pass filter at 2.5P
   ieq1             = 27
!qqW                = 2.5d0 * (PI2/TPERIOD_el)    ![rad/s]
   W                = 2.0d0 *  PI2   ! 2Hz        ![rad/s]
   D                = 0.7d0
   VAR              = accel2
   
   call filter_lowpass2 ( W             , D              , DT_el          , &
                          VAR                                             , &
                          QTCP_el (ieq1), QTCP1_el (ieq1), QTCP2_el (ieq1), &
                          QTC_el  (ieq1), QTC1_el  (ieq1),  QTC2_el (ieq1)    )
   
!--- Tower Acceleration: 2nd-order Notch filter at 3P
   ieq2             = 28
   W1               = 3.0d0 * (PI2/TPERIOD_el)    ![rad/s]
   W2               = W1
   D1               = 0.001d0
   D2               = 0.200d0
   VAR              = QTC_el (ieq1)
   DVAR             = QTC1_el(ieq1)
   DDVAR            = QTC2_el(ieq1)
   
   call filter_notch2 ( W1   , W2   , D1   , D2   , DT_el                , &
                        VAR  , DVAR , DDVAR                              , &
                        QTCP_el (ieq2) , QTCP1_el (ieq2), QTCP2_el (ieq2), &
                        QTC_el  (ieq2) , QTC1_el  (ieq2),  QTC2_el (ieq2)    )

   accelF           = QTC_el (ieq2) !for output
   GenTrqLSS_accel  = 0.d0
   PITCH_accel      = 0.d0
   BetaF_accel      = 0.d0

   if (TIME< dble(tta_ntime_init)*DT_el) return

!--- Tower Acceleration: PI controller
   ieq3             = 29
   GAIN             = GS_TTa
   CONP             = tta_kp
   CONI             = tta_ki
   COND             = 0.d0
   VAR              = QTC_el (ieq2)
   DVAR             = QTC1_el(ieq2)
   DDVAR            = QTC2_el(ieq2)

   call PID2_CONTROL( GAIN,     DT_el                              , &
                      CONP    , CONI    , COND                     , & !Kp ,Ki ,Kd
                      0.d0    , 0.d0    , 0.d0                     , & !Kp2,Ki2,Kd2
                      -1.d30  , 1.d30                              , & !minvar  , maxvar
                      -1.d30  , 1.d30                              , & !mindvar , maxdvar
                      -1.d30  , 1.d30                              , & !minddvar, maxddvar
                      VAR     , DVAR    , DDVAR                    , & !VAR , DVAR, DDVAR
                      0.d0    , 0.d0    , 0.d0                     , & !VAR2, DVAR2,DDVAR2
                      QTCP_el(ieq3), QTCP1_el(ieq3), QTCP2_el(ieq3), &
                      QTC_el (ieq3), QTC1_el (ieq3), QTC2_el (ieq3)    )

   GenTrqLSS_accel  = min( max(-QTC1_el(ieq3) * tta_GenGain   , -tta_MaxGenAmp   ), tta_MaxGenAmp   )
   PITCH_accel      = min( max(-QTC1_el(ieq3) * tta_PitchGain , -tta_MaxPitchAmp ), tta_MaxPitchAmp )
   BetaF_accel      = min( max(-QTC1_el(ieq3) * tta_FlapGain  , -tta_MaxFlapAmp  ), tta_MaxFlapAmp  )

   return

#    if   ASCII == 1
 2    open(23,file='control_TTa.dat',access='append')
       write(23,'(1000e15.5)')          &
#    elif ASCII == 0
 2    open(23,file='control_TTa.bin',access='append',form='UNFORMATTED')
       write(23              )          &
#    endif
         sngl(TIME                   ) ,& !1
         sngl(GenTrqLSS_accel/1000.d0) ,& !2 [kNm  ]
         sngl(PITCH_accel    *R2D    ) ,& !3 [deg  ]
         sngl(BetaF_accel    *R2D    ) ,& !4 [deg  ]
         sngl(accel1                 ) ,& !5 [m/s^2]
         sngl(accel2                 ) ,& !6 [m/s^2]
         sngl(accelF                 )    !7 [m/s^2]
      close(23)


 END Subroutine TTACC_control
!
!
!
!--------------------------------------------------------------------------------
!
!  Equation for double input PID controller: Gain (Ki + Kp s + Kd s^2) / s^2
!
!--------------------------------------------------------------------------------
 Subroutine PID2_CONTROL ( GAIN,     DT              , &
                           Kp      , Ki      , Kd    , &
                           Kp2     , Ki2     , Kd2   , &
                           min_qc  , max_qc          , &
                           min_dqc , max_dqc         , &
                           min_ddqc, max_ddqc        , &
                           err     , derr    , dderr , &
                           err2    , derr2   , dderr2, &
                           QTCP    , QTCP1   , QTCP2 , &
                           QTC     , QTC1    , QTC2      )
!--------------------------------------------------------------------------------
   implicit none

   real(8), intent(in   ) :: GAIN  , DT
   real(8), intent(in   ) :: Kp    , Ki    , Kd
   real(8), intent(in   ) :: Kp2   , Ki2   , Kd2
   real(8), intent(in   ) :: min_qc, max_qc, min_dqc, max_dqc, min_ddqc, max_ddqc
   real(8), intent(in   ) :: err   , derr  , dderr
   real(8), intent(in   ) :: err2  , derr2 , dderr2
   real(8), intent(in   ) :: QTCP  , QTCP1 , QTCP2
   real(8), intent(  out) :: QTC   , QTC1  , QTC2

   real(8) :: AM, AC, AK, AQ


   AM = 1.d0
   AC = 0.d0
   AK = 0.d0
   AQ = GAIN*( Ki *err  + Kp *derr  +Kd *dderr + &
               Ki2*err2 + Kp2*derr2 +Kd2*dderr2   )

   call Newmark_solve1 ( AM, AC, AK, AQ, DT, &
                         QTCP, QTCP1, QTCP2, &
                         QTC , QTC1 , QTC2     )

   call Newmark_saturate1 ( DT                        , &
                            min_qc  , max_qc          , &
                            min_dqc , max_dqc         , &
                            min_ddqc, max_ddqc        , &
                            QTCP    , QTCP1   , QTCP2 , &
                            QTC     , QTC1    , QTC2      )


 END Subroutine PID2_CONTROL
!--------------------------------------------------------------------------------
!
!  Equation for 2nd-order BandPass Filter : Gain (2zw s) / (s^2 + 2zw s + w^2)
!
!--------------------------------------------------------------------------------
 Subroutine filter_bandpass2 ( Gain  , w     , z     ,&
                               DT    , DVAR  , Q0M   ,&
                               QTCP  , QTCP1 , QTCP2 ,&
                               QTC   , QTC1  , QTC2     )
!--------------------------------------------------------------------------------
   implicit none

   real(8), intent(in   ) :: Gain, w, z, DT, DVAR
   real(8), intent(in   ) :: QTCP, QTCP1, QTCP2
   real(8), intent(  out) :: QTC , QTC1 , QTC2
   real(8), intent(  out) :: Q0M

   real(8) :: AM, AC, AK, AQ
   real(8) :: BITA, GAMA


   AM   = 1.d0
   AC   = 2.d0*z*w
   AK   =        w**2
   AQ   = 2.d0*Gain*z*w * DVAR

   call Newmark_solve1 (AM  , AC   , AK   , AQ, DT,&
                        QTCP, QTCP1, QTCP2,        &
                        QTC , QTC1 , QTC2            )

!--- derivative dQTC/dDVAR, based on Newmark
   BITA = 0.25d0
   GAMA = 0.50d0
   Q0M  = 2.d0*Gain*z*w / (AM/BITA/DT**2 + AC*GAMA/BITA/DT + AK) !dQ0/d(DVAR)
!  Q0M  = 2.d0*Gain*z*w / (AM/(GAMA*DT) + AC + AK*BITA*DT/GAMA)  !dQ1/d(DVAR)


 END Subroutine filter_bandpass2
!--------------------------------------------------------------------------------
!
!  Equation for 2nd-order Notch Filter: (s^2 + 2z1w1 s + w1^2) / (s^2 + 2z2w2 s + w2^2)
!
!--------------------------------------------------------------------------------
 Subroutine filter_notch2 ( W1, W2, Z1, Z2, DT, &
                            VAR , DVAR , DDVAR, &
                            QTCP, QTCP1, QTCP2, &
                            QTC , QTC1 , QTC2     )
!--------------------------------------------------------------------------------
   implicit none

   real(8), intent(in   ) :: W1, W2, Z1, Z2, DT
   real(8), intent(in   ) :: VAR , DVAR , DDVAR
   real(8), intent(in   ) :: QTCP, QTCP1, QTCP2
   real(8), intent(  out) :: QTC , QTC1 , QTC2

   real(8) :: AM, AC, AK, AQ


   AM = 1.d0
   AC = 2.d0*Z2*W2
   AK =         W2**2
   AQ =               DDVAR  &
       +2.d0*Z1*W1   * DVAR  &
       +        W1**2*  VAR

   call Newmark_solve1 (AM, AC, AK, AQ, DT, &
                        QTCP, QTCP1, QTCP2, &
                        QTC , QTC1 , QTC2     )


 END Subroutine filter_notch2
!--------------------------------------------------------------------------------
!
!  Equation for 2nd-order low-pass filter: w^2 / (s^2      + 2z w s + w^2)
!                                         =1   / (s^2 /w^2 + 2z/w s + 1  )
!
!--------------------------------------------------------------------------------
 Subroutine filter_lowpass2 ( W   , Z    , DT   ,&
                              VAR               ,&
                              QTCP, QTCP1, QTCP2,&
                              QTC , QTC1 , QTC2    )
!--------------------------------------------------------------------------------
   implicit none

   real(8), intent(in   ) :: W     , Z     , DT
   real(8), intent(in   ) :: VAR
   real(8), intent(in   ) :: QTCP  , QTCP1 , QTCP2
   real(8), intent(  out) :: QTC   , QTC1  , QTC2

   real(8) :: AM, AC, AK, AQ


   AM = 1.d0
   AC = 2.d0*Z*W
   AK =        W**2
   AQ =        W**2 * VAR

   call Newmark_solve1 (AM, AC, AK, AQ, DT, &
                        QTCP, QTCP1, QTCP2, &
                        QTC , QTC1 , QTC2     )


 END Subroutine filter_lowpass2
!--------------------------------------------------------------------------------
!
!  Equation for 1st-order lowpass filter: 1 / (s (1+s τ)) = 1 / (τ s^2 + s)
!
!--------------------------------------------------------------------------------
 Subroutine filter_lowpass1 ( TVLAG, DT   , VAR  , &
                              QTCP , QTCP1, QTCP2, &
                              QTC  , QTC1 , QTC2 , &
                              Q0M                    )
!--------------------------------------------------------------------------------
   implicit none

   real(8), intent(in   ) :: TVLAG , DT
   real(8), intent(in   ) :: VAR
   real(8), intent(in   ) :: QTCP  , QTCP1 , QTCP2
   real(8), intent(  out) :: QTC   , QTC1  , QTC2
   real(8), intent(  out) :: Q0M

   real(8) :: AM, AC, AK, AQ
   real(8) :: BITA, GAMA


   AM   = TVLAG
   AC   = 1.d0
   AK   = 0.d0
   AQ   = VAR

   call Newmark_solve1 ( AM, AC, AK, AQ, DT, &
                         QTCP, QTCP1, QTCP2, &
                         QTC , QTC1 , QTC2     )

!--- derivative dQTC1/dVAR, based on Newmark
   BITA = 0.25d0
   GAMA = 0.50d0
   Q0M  = 1.d0 / (AM/GAMA/DT + AC + AK*BITA*DT/GAMA) !dQTC1/dVAR


 END Subroutine filter_lowpass1
!----------------------------------------------------------------------
 Subroutine ACTUATOR_1 ( TVLAG   , DT             , &
                         MinVAR  , MaxVAR         , &
                         MinDVAR , MaxDVAR        , &
                         MinDDVAR, MaxDDVAR       , &
                         VAR                      , &
                         QTCP    , QTCP1   , QTCP2, &
                         QTC     , QTC1    , QTC2     )
!-----------------------------------------------------------------------

   implicit none

   real(8), intent(in   ) :: TVLAG , DT
   real(8), intent(in   ) :: MinVAR, MaxVAR, MinDVAR, MaxDVAR, MinDDVAR, MaxDDVAR
   real(8), intent(in   ) :: VAR
   real(8), intent(in   ) :: QTCP  , QTCP1 , QTCP2
   real(8), intent(  out) :: QTC   , QTC1  , QTC2

   real(8) :: AM, AC, AK, AQ


   AM   = TVLAG
   AC   = 1.d0
   AK   = 0.d0
   AQ   = VAR

   call Newmark_solve1 ( AM, AC, AK, AQ, DT, &
                         QTCP, QTCP1, QTCP2, &
                         QTC , QTC1 , QTC2     )

   call Newmark_saturate1 ( DT                        , &
                            MinVAR  , MaxVAR          , &
                            MinDVAR , MaxDVAR         , &
                            MinDDVAR, MaxDDVAR        , &
                            QTCP    , QTCP1   , QTCP2 , &
                            QTC     , QTC1    , QTC2      )

 END Subroutine ACTUATOR_1
!--------------------------------------------------------------------------------
!
!  Equation for 2nd-order pitch actuator: (2z w s + w^2) / (s^2       + 2z w s + w^2) =
!                                         (2z/w s + 1  ) / (s^2 / w^2 + 2z/w s + 1  )
!
!  numinator order < denuminator order
!--------------------------------------------------------------------------------
 Subroutine ACTUATOR_2 ( W, D    , DT             , &
                         MinVAR  , MaxVAR         , &
                         MinDVAR , MaxDVAR        , &
                         MinDDVAR, MaxDDVAR       , &
                         VAR     , DVAR           , &
                         QTCP    , QTCP1   , QTCP2, &
                         QTC     , QTC1    , QTC2     )
!--------------------------------------------------------------------------------
   implicit none

   real(8), intent(in   ) :: W, D  , DT
   real(8), intent(in   ) :: MinVAR, MaxVAR, MinDVAR, MaxDVAR, MinDDVAR, MaxDDVAR
   real(8), intent(in   ) :: VAR   , DVAR
   real(8), intent(in   ) :: QTCP  , QTCP1 , QTCP2
   real(8), intent(  out) :: QTC   , QTC1  , QTC2

   real(8) :: AM, AC, AK, AQ


   AM = 1.d0
   AC = 2.d0*D*W
   AK = W**2
   AQ = W**2 * VAR + 2.d0*D*W * DVAR

   call Newmark_solve1 ( AM, AC, AK, AQ, DT, &
                         QTCP, QTCP1, QTCP2, &
                         QTC , QTC1 , QTC2     )

   call Newmark_saturate1 ( DT                        , &
                            MinVAR  , MaxVAR          , &
                            MinDVAR , MaxDVAR         , &
                            MinDDVAR, MaxDDVAR        , &
                            QTCP    , QTCP1   , QTCP2 , &
                            QTC     , QTC1    , QTC2      )


 END Subroutine ACTUATOR_2
!-----------------------------------------------------------------------
 Subroutine Newmark_solve1 (AM  , AC   , AK   , AQ , DT, &
                            QTCP, QTCP1, QTCP2         , &
                            QTC , QTC1 , QTC2             )
!-----------------------------------------------------------------------

   implicit none

   real(8), intent(in ) :: DT
   real(8), intent(in ) :: AM  , AC   , AK   , AQ
   real(8), intent(in ) :: QTCP, QTCP1, QTCP2
   real(8), intent(out) :: QTC , QTC1 , QTC2

   real(8) :: BITA, GAMA, c1,c2,c3,c4
   real(8) :: AKEQ, UPRE, UTPRE, RQ, U


   BITA  = 0.25d0
   GAMA  = 0.50d0
   c1    =(0.50d0 - BITA)*DT**2
   c2    =(1.00d0 - GAMA)*DT
   c3    = 1.00d0 /(BITA*DT**2)
   c4    = GAMA   /(BITA*DT   )

   AKEQ  = AM*c3     + AC*c4     + AK
   UPRE  = QTCP      + QTCP1*DT  + QTCP2*c1
   UTPRE = QTCP1     +             QTCP2*c2
   RQ    = AQ + ( AM*c3 + AC*c4 )*UPRE - AC * UTPRE
   U     = RQ/AKEQ
   QTC   = U
   QTC1  = UTPRE + c4*(U-UPRE)
   QTC2  =         c3*(U-UPRE)


 END Subroutine Newmark_solve1
!-----------------------------------------------------------------------
 Subroutine Newmark_saturate1 ( DT                        , &
                                MinVAR  , MaxVAR          , &
                                MinDVAR , MaxDVAR         , &
                                MinDDVAR, MaxDDVAR        , &
                                QTCP    , QTCP1   , QTCP2 , &
                                QTC     , QTC1    , QTC2      )
!-----------------------------------------------------------------------

   implicit none

   real(8), intent(in   ) :: DT
   real(8), intent(in   ) :: MinVAR, MaxVAR, MinDVAR, MaxDVAR, MinDDVAR, MaxDDVAR
   real(8), intent(in   ) :: QTCP  , QTCP1 , QTCP2
   real(8), intent(inout) :: QTC   , QTC1  , QTC2

   real(8) :: UPRE, U1PRE, BITA, GAMA, c1,c2


      BITA  = 0.25d0
      GAMA  = 0.50d0
      c1    =(0.50d0 - BITA)*DT**2
      c2    =(1.00d0 - GAMA)*DT
!--- Satate acceleration (DDVAR)
   if (QTC2 < MinDDVAR .or. QTC2 > MaxDDVAR ) then
      UPRE  =  QTCP  + QTCP1*DT + QTCP2*c1
      U1PRE =          QTCP1    + QTCP2*c2
      QTC2  = min ( max( QTC2, MinDDVAR ), MaxDDVAR )
      QTC   = UPRE  + BITA * DT**2 * QTC2
      QTC1  = U1PRE + GAMA * DT    * QTC2
   endif

!--- Saturate velocity (DVAR)
   if (QTC1 < MinDVAR .or. QTC1 > MaxDVAR ) then
      UPRE  =  QTCP  + QTCP1*DT + QTCP2*c1
      U1PRE =          QTCP1    + QTCP2*c2
      QTC1  = min ( max( QTC1, MinDVAR  ), MaxDVAR )
      QTC2  = (QTC1-U1PRE) / (GAMA*DT)
      QTC   = UPRE  + BITA  * DT**2 * QTC2
   endif

!--- Saturate position (VAR)
   if (QTC < MinVAR .or. QTC > MaxVAR ) then
      UPRE  =  QTCP  + QTCP1*DT + QTCP2*c1
      U1PRE =          QTCP1    + QTCP2*c2
      QTC   = min ( max( QTC , MinVAR   ), MaxVAR )
      QTC2  = (QTC-UPRE) / (BITA * DT**2)
      QTC1  = U1PRE + GAMA * DT    * QTC2
   endif


 END Subroutine Newmark_saturate1
!----------------------------------------------------------------------
 Subroutine filterN ( A,B,C,D, QTCP,QTCP1,QTCP2 ,QTC,QTC1,QTC2, DT, N, XX,XX1, YY,YY1 )
!----------------------------------------------------------------------
! [A,B,C,D]=func(n,db,[w1,w2],type,'s')
!  func: cheby1, cheby2, butter, ellip
!  n:order (i.e.2)
!  w1,w2 frequency range [rad/s]
!  type='low' for lowpass, 'high' for highpass,'bandpass' for bandpass and 'stop' for bandstop
!  's' for analog filter
!----------------------------------------------------------------------
!---- from lowpass wind
!   N           = 2
!   damp        = 0.20d0 !0.20 default
!   A (1:2,1:2) = 0.d0;   A (1,1) =-damp;    A (2,1) = 1.00d0
!   B (1:2    ) = 0.d0;   B (1  ) =-damp**2; B (2  ) = damp
!   C (1:2    ) = 0.d0;   C (1  ) = 0.00d0;  C (2  ) = 1.00d0
!   D           = 0.d0
!----------------------------------------------------------------------

   implicit none

   integer, intent(in   ) :: N
   real(8), intent(in   ) :: A (N,N), B  (N), C  (N), D
   real(8), intent(in   ) :: QTCP(N), QTCP1(N), QTCP2(N)
   real(8), intent(in   ) :: DT
   real(8), intent(in   ) :: XX     , XX1    !, XX2
   real(8), intent(  out) :: QTC (N), QTC1 (N), QTC2 (N)
   real(8), intent(  out) :: YY     , YY1    !, YY2

   real(8)  ::  AM(N,N), AC(N,N), AK(N,N), AQ(N)
   integer  ::  i


   AK(1:N,1:N) = 0.d0;
   AC(1:N,1:N) =-A(1:N,1:N)
   AM(1:N,1:N) = 0.d0; forall(i = 1:N) AM(i,i) = 1.d0
   AQ(1:N    ) = B(1:N)*XX

   call NEWMARKN (QTCP,QTCP1,QTCP2, QTC,QTC1,QTC2 ,DT, AM,AC,AK,AQ, N)

   YY  = dot_product (C(1:N), QTC1(1:N)) + D*XX   !value
   YY1 = dot_product (C(1:N), QTC2(1:N)) + D*XX1  !1st time derivative


 END Subroutine filterN
!----------------------------------------------------------------------
 Subroutine NEWMARKN (QTCP,QTCP1,QTCP2, QTC,QTC1,QTC2, DT, AM,AC,AK,AQ, N) 
!----------------------------------------------------------------------

   implicit none

   real(8), intent(in   ) :: DT
   integer, intent(in   ) :: N
   real(8), intent(in   ) :: AM(N,N), AC(N,N) , AK(N,N), AQ(N)
   real(8), intent(in   ) :: QTCP(N), QTCP1(N), QTCP2(N)
   real(8), intent(  out) :: QTC (N), QTC1 (N), QTC2 (N)

   integer :: ierr
   real(8) :: RQ(N), AKEQ(N,N)
   integer :: INDXO(N)
   real(8) :: BITA, GAMA, c1,c2,c3,c4
   real(8) :: UPRE(N), UTPRE(N), d


   BITA           = 0.25d0
   GAMA           = 0.50d0
   c1             =(0.50d0 - BITA)*DT**2
   c2             =(1.00d0 - GAMA)*DT
   c3             = 1.00d0 /(BITA*DT**2)
   c4             = GAMA   /(BITA*DT   )

   UPRE (1:N    ) = QTCP  (1:N    )    + QTCP1(1:N)*DT + QTCP2(1:N)*c1
   UTPRE(1:N    ) = QTCP1 (1:N    )    +                 QTCP2(1:N)*c2
   AKEQ (1:N,1:N) = AM    (1:N,1:N)*c3 +         AC (1:N,1:N)*c4 + AK (1:N,1:N)
   RQ   (1:N    ) = AQ    (1:N    )    + matmul( AM (1:N,1:N)*c3 + AC (1:N,1:N)*c4, UPRE (1:N) ) &
                                       - matmul( AC (1:N,1:N)                     , UTPRE(1:N) )

   call LUDCMP (AKEQ, N, N, INDXO, d , ierr); if (ierr==1) then;write(*,*)'LUDCMP in NEWMARKN error',N;stop;endif
   call LUBKSB (AKEQ, N, N, INDXO, RQ      )

   QTC  (1:N    ) = RQ   (1:N)
   QTC1 (1:N    ) = UTPRE(1:N) + c4*(RQ(1:N)-UPRE(1:N))
   QTC2 (1:N    ) =              c3*(RQ(1:N)-UPRE(1:N))


 END Subroutine NEWMARKN
!----------------------------------------------------------------------------------------
 Subroutine tip_brake
!----------------------------------------------------------------------------------------
!-- This is a special case used for an old a stall regulated WT for the REWIND project
!-- The sub-routines use the pre-curved angles in order to rotate the tip brake.
!-- This means that atm the blade can't be also pre-bended or swept.
!-- This assumption is valid, because the old small blades used to be straight.
!
!   Remarks:
!-- The time at which the brake is applied must be specified by the user [the controller
!-- won't detect any overspeed].
!-- At the moment the tip-brake starts at the last sub-body.
!-- The rotation is a linear ramp (t).
!----------------------------------------------------------------------------------------

 Use Cbeam
 Use Ctrl_mod

   implicit none

   integer       :: nbod, nnsub, nnsubT
   real(8)       :: TIME_fault, TIME_stop, durationBrake, durationGen, Angle_max, OverSpeed, ANGLE, c_down
   real(8), save :: TIME_tipbrake=-1.d0


!---- Settings
   TIME_stop     =  15.0d0                                                    !time for the normal shut-down
   TIME_fault    =9999.0d0                                                    !time for the network loss (zero Tgen)
   durationBrake =   0.5d0                                                    !brake duration to turn linearly from 0 to Angle_max degrees [sec]
   durationGen   =   0.5d0                                                    !Gen   duration to linearky become zero [sec]
   Angle_max     =  85.d0/R2D                                                 !set angle
   OverSpeed     =  1.08d0                                                    !overspeed, if omega>omega_nominal*overspeed --> brake activation

   if     (ICMOD /= 4               ) return
   if     (TIME < TIME_fault .and. &
           TIME < TIME_stop         ) return

   if     (TIME >= TIME_fault) then                                            !time at which generator torque is zero
      TGenLSS_el   = 0.d0
      TGenLSS_D_el = 0.d0
      if ((-UT1_el(NDFBT_el+NQSW)<OMEGAR*OverSpeed).and. &
          (TIME_tipbrake         <0.d0)                    ) return            !still airbrake is not enabled, waiting for overspeed occurance

   elseif (TIME >= TIME_stop ) then                                            !time at which normal shut-down initiates

      c_down       = max(0.d0,1.d0-(TIME-Time_stop)/durationGen)
      TGenLSS_el   = 0.d0 !TGenLSS_el   * c_down
      TGenLSS_D_el = 0.d0 !TGenLSS_D_el * c_down
   endif


   if (TIME_tipbrake<0.d0) TIME_tipbrake = TIME

   IAX_Ben_Swe = 2                                                            !transform the pre-curved angle around pitch axis

   ANGLE = -min(Angle_max, Angle_max * (TIME-TIME_tipbrake)/durationBrake)    !negative as pitch

   if (TIME_stop<0.d0) &
   ANGLE = -Angle_max

   do nbod = 1, NBLADE_el

      body(nbod)%PRECURV(           :) = 0.d0;                                !zero possible pre-bend/sweep angles
                               nnsub   = body(nbod)%NBODSUB_el                !the tip brake starts at the last sub-body
                               nnsubT  = body(nbod)%NBODSUB_el + 1            !Total nodes or each blade
      body(nbod)%PRECURV(nnsub:nnsubT) = ANGLE                                !All nodes after tip brake are rotated
   enddo

   open(122,file='stall_ctrl.dat',access='append')
      write(122,'(10e15.5)') TIME,-ANGLE*R2D,Time_tipbrake, c_down, TGenLSS_el
   close(122)
   

 END Subroutine tip_brake
!----------------------------------------------------------------------------------------
 Subroutine pitch_change
!----------------------------------------------------------------------------------------
!-- This is a special case used for thalis project in order to simulate a simple pitch step.
!
!   Remarks:
!-- Very similar sb with tip_brake used for REWIND project [tip brake]
!-- The rotation is a linear ramp (t).
!----------------------------------------------------------------------------------------

 Use Cbeam
 Use Ctrl_mod

   implicit none

   integer :: nbod,iap,idp0
   real(8) :: TIME_fault, duration, Angle_max, ANGLE


!---- Settings
   TIME_fault = 60.0d0                                                        !time to start the pitch step [sec]
   duration   = 20.d0/360.d0*TPERIOD_el                                       !brake duration to turn linearly from 0 to Anle_max degrees [sec]
   Angle_max  =  2.d0/R2D                                                     !set angle

   if (TIME  <  TIME_fault) return

   if (TIME_BRAKE<0.d0) TIME_BRAKE = TIME

   ANGLE = -min(Angle_max, Angle_max * (TIME-TIME_BRAKE)/duration)            !negative as pitch

!--- Update pitch values based on pitch demand
   do nbod   = 1, NBLADE_el
      iap    = NDFBT_el + NQSP + ACCU%NQSACC_el(nbod-1) + 5
      idp0   = 5

      UT_el  (iap) = ANGLE
      UT1_el (iap) = 0.d0
      UT2_el (iap) = 0.d0
   enddo !nbod
   

 END Subroutine pitch_change
!----------------------------------------------------------------------------------------
 Subroutine pitch_heli
!----------------------------------------------------------------------------------------
!-- This subroutine sets a priori the blade pitch of the helicopter, based on the given 
!    collective and cyclic constants and the azimuth position.
!----------------------------------------------------------------------------------------

 Use Cbeam

   implicit none

   integer :: nbod,iap
   real(8) :: q0,qc,qs,psi,psit,psitt  !wt=azimuth angle


!! nbsub = body(nbod)%NBODSUB_el !last subbody
!! nb_el = body(nbod)%NBODTGLB_el(nbsub)
!! nel   = subbody(nb_el)%NTEB_el
!! nod   = NNPE_el

!! call beam_loads (nbod,nbsub,nel,nod, FGG, FLG, FLL)

   if (IAPPL_el /= 2) return
!--- Update pitch values based on pitch demand
   do nbod         = 1, NBLADE_el
!--- settings (also a PI controller could be added here)
      q0           = PITCH_el  (nbod)
      qc           = PITCHC_el (nbod)
      qs           = PITCHS_el (nbod)

      psi          = UT_el (NDFBT_el+NQSW) + PIhalf + PHI0  (nbod) !- UThub_el
      psit         = UT1_el(NDFBT_el+NQSW)
      psitt        = UT2_el(NDFBT_el+NQSW)

              iap  = NDFBT_el + NQSP + ACCU%NQSACC_el(nbod-1) + 5
      UT_el  (iap) = q0 + qc*        dcos(psi)                        +  qs*        dsin(psi)
      UT1_el (iap) =    - qc*  psit *dsin(psi)                        +  qs*  psit *dcos(psi)
      UT2_el (iap) =    - qc*( psitt*dsin(psi) + psit**2*dcos(psi) )  +  qs*( psitt*dcos(psi) - psit**2*dsin(psi) )
   enddo !nbod


 END Subroutine pitch_heli
!----------------------------------------------------------------------------------------
 Subroutine set_omega_pitch ( NTIME, it )
!----------------------------------------------------------------------------------------
!-- This subroutine sets the blade pitch and the omega based on given prescribed timeseries.
!----------------------------------------------------------------------------------------

 Use Cbeam

   implicit none

   integer, intent(in)          :: NTIME, it
   integer, parameter           :: np=20000
   real(8), dimension(np), save :: time_rd
   real(8), dimension(np), save ::  azim_rd , dazim_rd , ddazim_rd
   real(8), dimension(np), save :: pitch_rd ,dpitch_rd ,ddpitch_rd
   real(8),                save ::  azim_int, dazim_int, ddazim_int
   real(8),                save :: pitch_int,dpitch_int,ddpitch_int
   real(8),                save :: dt_rd    , dt2_rd   , dtsqr_rd
   integer,                save :: nt_rd
   integer                      :: i, nbod


!--- Initialization
   if (NTIME<=1 .and. it<=1) then
     open(1,file="perform.inp")
     do i = 1, np
        read(1,*, end=30) time_rd(i),dazim_rd(i),pitch_rd(i) !sec, rad/s, rad
        dazim_rd(i)=dazim_rd(i)*pi/30.    ! rpm->rad/s
        pitch_rd(i)=pitch_rd(i)*pi/180.   ! deg->rad
     enddo
        write(* ,*) 'end of perform.inp was not reatched - increase np'
        write(10,*) 'end of perform.inp was not reatched - increase np'
        stop
30   close (1)
     nt_rd    = i-1
     dt_rd    = time_rd(2)-time_rd(1)
     dt2_rd   = 2.d0*dt_rd
     dtsqr_rd = dt_rd**2

        write(* ,*)'ICMOD5-nt_rd',nt_rd,dt_rd
        write(10,*)'ICMOD5-nt_rd',nt_rd,dt_rd

                  i  = 1
          azim_rd(i) = 0.d0
        ddazim_rd(i) =                    (dazim_rd(i+1)-dazim_rd(i  )) /    dt_rd
        dpitch_rd(i) =                    (pitch_rd(i+1)-pitch_rd(i  )) /    dt_rd
       ddpitch_rd(i) = (pitch_rd(i+2)-2.d0*pitch_rd(i+1)+pitch_rd(i  )) / dtsqr_rd
     do           i  = 2, nt_rd-1
          azim_rd(i) =   azim_rd(i-1)+     dazim_rd(i  )                *    dt_rd
        ddazim_rd(i) = (dazim_rd(i+1)                   -dazim_rd(i-1)) /   dt2_rd
        dpitch_rd(i) = (pitch_rd(i+1)                   -pitch_rd(i-1)) /   dt2_rd
       ddpitch_rd(i) = (pitch_rd(i+1)-2.d0*pitch_rd(i  )+pitch_rd(i-1)) / dtsqr_rd
     enddo
                  i  = nt_rd
          azim_rd(i) =                      azim_rd(i-1)+dazim_rd(i  )  *    dt_rd
        ddazim_rd(i) = (dazim_rd(i  )-     dazim_rd(i-1)              ) /    dt_rd
        dpitch_rd(i) = (pitch_rd(i  )-     pitch_rd(i-1)              ) /    dt_rd
       ddpitch_rd(i) = (pitch_rd(i  )-2.d0*pitch_rd(i-1)+pitch_rd(i-2)) / dtsqr_rd
   endif


!  if (it>1) return


!--- Set interpolated prescribed azim,pitch and corresponding time derivatives
      call LIN_INT(TIME,   azim_int,time_rd,   azim_rd, nt_rd, np)
      call LIN_INT(TIME,  dazim_int,time_rd,  dazim_rd, nt_rd, np)
      call LIN_INT(TIME, ddazim_int,time_rd, ddazim_rd, nt_rd, np)
      call LIN_INT(TIME,  pitch_int,time_rd,  pitch_rd, nt_rd, np)
      call LIN_INT(TIME, dpitch_int,time_rd, dpitch_rd, nt_rd, np)
      call LIN_INT(TIME,ddpitch_int,time_rd,ddpitch_rd, nt_rd, np)

!--- set omega
   i = NDFBT_el + NQSW
      UT_el  (i) = -  azim_int
      UT1_el (i) = - dazim_int
      UT2_el (i) = -ddazim_int

!--- set pitch
   do nbod       = 1, NBLADE_el
              i  = NDFBT_el + NQSP + ACCU%NQSACC_el(nbod-1) + 5
      UT_el  (i) = -  pitch_int
      UT1_el (i) = - dpitch_int
      UT2_el (i) = -ddpitch_int
   enddo !nbod


 END Subroutine set_omega_pitch
!
!
!
!--- select the controller
!0. OC3
!1. DLL
!2. DTU_v1
!3. DTU_v2
!33.DTU_v2.3

#  ifndef    CONTROLLER_TYPE
#     define CONTROLLER_TYPE 0
#  endif

#  if   CONTROLLER_TYPE == 0
#     include "./ctrl/oc3_control.f90"
#  elif CONTROLLER_TYPE == 1
#     include "./ctrl/dll_control.f90"
#  elif CONTROLLER_TYPE == 2
#     include "./ctrl/10MW_control.f90"
#  elif CONTROLLER_TYPE == 3
#     include "./ctrl/dtu_control_v2.f90"
#  elif CONTROLLER_TYPE == 33
#     include "./ctrl/dtu_control_v2.f90"
#  endif
