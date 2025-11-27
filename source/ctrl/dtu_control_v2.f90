!-------------------------------------
!dtu_we_controller.f90
!          subroutine init_regulation          (array1, array2)
!          subroutine init_regulation_advanced (array1, array2)
!          subroutine update_regulation        (array1, array2)
!
!dtu_we_controller_fcns.f90                    
!          function   switch_spline            (x, x0, x1)
!          function   interpolate              (x, x0, x1, f0, f1)
!          function   GetOptiPitch             (wsp)
!          function   PID                      (stepno,dt,kgain,PIDvar,error)
!          function   PID2                     (stepno,dt,kgain,PIDvar,error,added_term)
!          subroutine damper                   (stepno, dt, x, filters, y, x_filt)      **new
!
!misc_mod.f90                                  
!          function   lowpass1orderfilt        (dt, stepno, filt, x)
!          function   lowpass2orderfilt        (dt, stepno, filt, x)
!          function   notch2orderfilt          (dt,stepno,filt,x)
!          function   bandpassfilt             (dt, stepno, filt, x)
!          function   timedelay                (dt, stepno, filt, Td, x)                **new
!
!turbine_controller.f90
!          subroutine turbine_controller       (CtrlStatus, GridFlag, GenSpeed, PitchVect, wsp, Pe, TTAccVect, &
!          subroutine normal_operation         (GenSpeed, PitchVect, wsp, Pe, TTfa_acc, GenTorqueRef, PitchColRef, &
!          subroutine start_up                 (CtrlStatus, GenSpeed, PitchVect, wsp, GenTorqueRef, PitchColRef, dump_array)
!          subroutine shut_down                (CtrlStatus, GenSpeed, PitchVect, wsp, GenTorqueRef, PitchColRef, dump_array)
!          subroutine monitoring               (CtrlStatus, GridFlag, GenSpeed, TTAcc, dump_array)
!          subroutine torquecontroller         (GenSpeed, GenSpeedFilt, dGenSpeed_dtFilt, PitchMean, WSPfilt, &
!          subroutine pitchcontroller          (GenSpeedFilt, dGenSpeed_dtFilt, PitchMeanFilt, Pe, PitchMin, &
!          subroutine rotorspeedexcl           (GenSpeedFilt, GenTorque, Qg_min_partial, GenTorqueMax_partial, GenSpeedFiltErr, &
!          subroutine drivetraindamper         (GenSpeed, Qdamp_ref, dump_array)
!          subroutine towerdamper              (TTfa_acc, theta_dam_ref, dump_array)
!
!safety_system.f90
!          subroutine safety_system            (stepno, deltat, omega, TTAccVect, EmergPitchStop, ActiveMechBrake, dump_array)

!---------------old v1

!../10MW_controller.f90:2
!          Subroutine init_regulation          (array1,gentrq,collpit)
!          Subroutine update_regulation        (array1,array2)

!10MW_controller_fcns.f90
!          function   switch_spline            (x,x0,x1)
!          function   interpolate              (x,x0,x1,f0,f1)
!          function   GetOptiPitch             (wsp)
!          function   lowpass1orderfilt        (dt,stepno,filt,x)
!          function   lowpass2orderfilt        (dt,stepno,filt,x)
!          function   notch2orderfilt          (dt,stepno,filt,x)
!          function   bandpassfilt             (dt,stepno,filt,x)
!          function   PID                      (stepno,dt,kgain,PIDvar,error)
!          function   PID2                     (stepno,dt,kgain,PIDvar,error)


!-- Interface for the 10MW DTU REFERENCE controller. Sets the communicating variables between gast and discon.
!-- The drive train damper is implimented inside the controller and at moment the gain is ZERO.
!-- The generator lag is introduced here, as well as the pitch actuator.

#  if CONTROLLER_TYPE == 3
 include "./ctrl/dtu_v2/misc_mod.f90"
 include "./ctrl/dtu_v2/dtu_we_controller_fcns.f90"
 include "./ctrl/dtu_v2/safety_system.f90"
 include "./ctrl/dtu_v2/turbine_controller.f90"
 include "./ctrl/dtu_v2/dtu_we_controller.f90"
#  elif CONTROLLER_TYPE == 33
 include "./ctrl/dtu_v2.3/misc_mod.f90"
 include "./ctrl/dtu_v2.3/dtu_we_controller_fcns.f90"
 include "./ctrl/dtu_v2.3/safety_system.f90"
 include "./ctrl/dtu_v2.3/turbine_controller.f90"
 include "./ctrl/dtu_v2.3/dtu_we_controller.f90"
#  endif
!----------------------------------------------------------------------
 Subroutine controller (NTIME,it,ISTEP)
!----------------------------------------------------------------------

 Use Cbeam, only : TIME, TREF, R2D, PI2, ICMOD, RAT_GEAR
 Use Ctrl_mod
 Use dtu_we_controller_mod

   implicit none

   integer, intent(in) :: NTIME, it, ISTEP
   integer       :: i          , idtu, iLSS
   real(8), save :: array1(100), array2(100)
   real(8)       :: array0( 50)
   real(8)       :: kk1, kk2

   idtu=0 !1.dtu-ctrl, 0.built-in

   if (ICMOD == 0) return

   if (ICMOD == 2) then
       TGenLSS_el   = 0.d0
       TMLossLSS_el = 0.d0
       return
   endif

   goto (1,2,3,4), ISTEP

!--- main part
 1 if (it>1.and.idtu==1) return

   call ctrl_sensors

   if ((NTIME==1).and.(it==1)) then
!--- Initialization of the 10Mw controller
      call controller_init (collpit, TREF, array0, array1)

!------ Define parameters
      write(*,*)'array0'
      do i=1,39
      write(*,*)array0(i)
      enddo

!------ inputs
      Gen_lag         = array0( 1)       ! [sec]    !  0.02d0         !generator lag                          [sec]
      itypePit_act    = int(array0(2))
      freqPit_act     = array0( 3)*PI2   ! [rad/s] #!  1.6d0*PI2      !frequency         of pitch actuator    [rad/s] ### 1  Hz
      dampPit_act     = array0( 4)       ! [-]     #!  0.8d0          !damping factor    of pitch actuator    [-]     ### 0.7
      MinPitAng_act   = array0( 5)/R2D   ! [rad]    ! -5.d0 / R2D     !lower angle limit of pitch actuator    [rad]
      MaxPitAng_act   = array0( 6)/R2D   ! [rad]    ! 90.d0 / R2D     !upper angle limit of pitch actuator    [rad]
      MaxPitRat_act   = array0( 7)/R2D   ! [rad/s]  ! 10.d0 / R2D     !speed       limit of pitch actuator    [rad/s]   ### 8
      MaxPitAcc_act   = array0( 8)/R2D   ! [rad/s^2]! 15.d0 / R2D     !accelerationlimit of pitch actuator    [rad/s^2] ### 8
      DTD_gain(1)     = array0( 9)       !          !400.d0*RAT_GEAR**2
      DTD_freq(1)     = array0(10)*PI2   ! [rad/s]  !  1.870d0 !Hz free-free
      DTD_damp(1)     = array0(11)       ! [-]      !  0.100d0
      DTD_gain(2)     = array0(12)       !          ! 00.d0
      DTD_freq(2)     = array0(13)*PI2   ! [rad/s]  ! 10.000d0 !Hz
      DTD_damp(2)     = array0(14)       ! [-]      !  0.100d0
      MaxGenTrqDTD    = array0(15)       ! [Nm]
!------ power controller
      power_rated     = array1( 1) * 1000.d0
      omega_min       = array1( 2)
      omega_max       = array1( 3)
      vs_Qmax         = array1( 4)
      vs_Kopt         = array1(11)
      vs_Kp           = array1(12)
      vs_Ki           = array1(13)
      vs_Kd           = array1(14)
      freq_free1      = array1(10)*PI2
!     Freq_3p         =
      vp_pitmin       = array1( 5)/R2D
      vp_pitmax       = array1( 6)/R2D
      vp_pitratmax    = array1( 7)/R2D
      vp_pitratmin    =-vp_pitratmax
      vp_ipowtor      = int(array1(15)) !constant  power (1) or torque (2)
      vp_Kp           = array1(16)
      vp_Ki           = array1(17)
      vp_Kd           = array1(18)
      vp_Kp2          = array1(19)
      vp_Ki2          = array1(20)
      iLSS            = int(array1(76))
      if (iLSS==1) then
         gearctrl     = 1.d0
      else
         gearctrl     = RAT_GEAR
      endif
!------ GS
      kk1             = array1(21)/R2D
      kk2             = array1(22)/R2D
      GS_NL_lim       = array1(23)
!------ inputs for IPC (Blade root out-of-plane momemnts)
      ipc_nhT         = int(array0(16))  ! Number of harmonic (0:no, 1:1P, 2:2P) [-]
      ipc_PitchGain   = array0(17)       ! GAIN_PITCH         (default=1)        [-]
      ipc_FlapGain    = array0(18)       ! GAIN_FLAP          (default=1)        [-]
      ipc_Ki (1:2)    = array0(19:20)    !
      ipc_Kp (1:2)    = array0(21:22)    ! 2.5d0*Ki(:)
      ipc_MaxPitchAmp = array0(23)/R2D   ! Maximum Pitch Amplitude               [rad]
      ipc_MaxFlapAmp  = array0(24)/R2D   ! Maximum Flap  Amplitude               [rad]
      ipc_ntime_init  = int(array0(25))  ! time step to start IPC ntime_init     [-]
!------ inputs for TTa
      tta_GenGain     = array0(26)       !                                       [-]
      tta_PitchGain   = array0(27)       !                                       [-]
      tta_FlapGain    = array0(28)       !                                       [-]
      tta_Ki          = array0(29)       !
      tta_Kp          = array0(30)       ! 2.5d0*Ki(:)
      tta_MaxGenAmp   = array0(31)       ! Maximum Gen   Amplitude               [Nm]
      tta_MaxPitchAmp = array0(32)/R2D   ! Maximum Pitch Amplitude               [rad]
      tta_MaxFlapAmp  = array0(33)/R2D   ! Maximum Flap  Amplitude               [rad]
      tta_ntime_init  = int(array0(34))  ! time step to start TTa ntime_init     [-]
!------ inputs for flap Actuator
      itypeFlap_act   = int(array0(35))
      timeFlap_act    = array0(36)       ! [rad/s] #!  1.6d0*PI2      !frequency         of flap  actuator    [rad/s] ### 1  Hz
      MinFlapAng_act  = array0(37)/R2D   ! [rad]    ! -5.d0 / R2D     !lower angle limit of flap  actuator    [rad]
      MaxFlapAng_act  = array0(38)/R2D   ! [rad]    ! 90.d0 / R2D     !upper angle limit of flap  actuator    [rad]
      MaxFlapRat_act  = array0(39)/R2D   ! [rad/s]  ! 10.d0 / R2D     !speed       limit of flap  actuator    [rad/s]   ### 8

      if (itypePit_act  /=3) then; write(*,*)'only Pitch actuator type 3 is supported atm'; stop; endif
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
!---   11. LOWPASS-2nd (genspeed)                             [QTC ]      0.4Hz    , z=0.7
!---   12. NOTCH  -2nd (genspeed) free-free                   [QTC ]      free-free, z=0.1
!---   13. LOWPASS-2nd (power   )                             [QTC1]omega 1P       , z=0.7
!---   14. NOTCH  -2nd (power   ) free-free                   [QTC1]omega free-free, z=0.1

!---15-18. IPC filters LP2 (Mtilt1P, Myaw1P, Mtilt2P, Myaw2P) [QTC]       2.5P     , z=0.7
!---19-22. IPC Notch       (Mtilt1P, Myaw1P, Mtilt2P, Myaw2P) [QTC]       3.0P     , z=0.2
!---23-26. IPC PI          (Mtilt1P, Myaw1P, Mtilt2P, Myaw2P) [QTC1]
!---   27. Ttop accel  LP2                                    [QTC ]      2.5P     , z=0.7
!---   28. Ttop accel  Notch                                  [QTC ]      3.0P     , z=0.2
!---   29. Ttop accel  PI                                     [QTC1]
!---30-32. Flap Actuator  1st order                           [QTC1]
!qq
!---   33. PI   Torque                                        [QTC1]
!---   34. PI2  pitchcol                                      [QTC1]
!---   35. LP2  omega err                                     [QTC ]
!---   36. Not2 omega err                                     [QTC ]
!---   37. Not2 omega err                                     [QTC ]
!---   38. Not2 omega err                                     [QTC ]
!---   39. LP2  power err                                     [QTC ]

!--[15-18. Ttop accel  (notch1P 3P,lowpass1P 2.5P,notch2P 3P,lowpass2P 2.5P)]
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

      GS_IPC0       = 0.01d0
      GS_TTa0       = 0.01d0
   endif !initialization


!--- Main controller call, pass the communication variables
   array1( 1) = TIME            !  1: general time                            [s]
   array1( 2) = GenSpeed        !  2: constraint bearing1 shaft_rot 1 only 2  [rad/s] Generator speed (Default LSS, if HSS insert gear ratio in input #76)
   array1( 3) = BlPitch(1)      !  3: constraint bearing2 pitch1 1 only 1     [rad]
   array1( 4) = BlPitch(2)      !  4: constraint bearing2 pitch2 1 only 1     [rad]
   array1( 5) = BlPitch(3)      !  5: constraint bearing2 pitch3 1 only 1     [rad]
   array1( 6) = WindVel(1)      !  6: wind free_wind 1 0.0 0.0 hub height     [m/s] global coords at hub height
   array1( 7) = WindVel(2)      !  7: wind free_wind 1 0.0 0.0 hub height     [m/s] global coords at hub height
   array1( 8) = WindVel(3)      !  8: wind free_wind 1 0.0 0.0 hub height     [m/s] global coords at hub height
   array1( 9) = array2 (5)      !  9: mech. power  ; [W]
   array1(10) = 0               ! 10: grid flag  ; [1=no grid,0=grid]
!qq local or global??
   array1(11) = ttop_acc(1)     ! 11: Tower top x-acceleration  ; [m/s^2]                                                                                   
   array1(12) = ttop_acc(2)     ! 12: Tower top y-acceleration  ; [m/s^2]                                                                                   

  if (idtu==1) then
   call update_regulation(array1, array2)

   GenTrqDem      = array2(1  )
   PitAngDem(1:3) = array2(2:4)

  else

!--- basic power control
   call power_control
  endif

!--- IPC
   call IPC_control  (ISTEP)

!--- TTACC (actuator, coll pitch, gen, coll flap)
   call TTACC_control(ISTEP)

!--- Drive train damper, Generator lag and Pitch & Flap Actuators
   call DTD_GenLag_PitFlapAct

   return

!--- Update for iterations [output results]
 2 continue

   QTCP_el (1:NQSC0) = QTC_el (1:NQSC0);
   QTCP1_el(1:NQSC0) = QTC1_el(1:NQSC0);
   QTCP2_el(1:NQSC0) = QTC2_el(1:NQSC0);

#  if   ASCII == 1
      open(27,file='controller.dat',access='append')
       write(27,100)              &
#  elif ASCII == 0
      open(27,file='controller.bin',access='append',form='UNFORMATTED')
       write(27    )              &
#  endif
         sngl(TIME    )          ,&
        (sngl(array2(i)),i=1,37)
      close(27)


#  if   ASCII == 1
      open(27,file='control.dat',access='append')
       write(27,100)              &
#  elif ASCII == 0
      open(27,file='control.bin',access='append',form='UNFORMATTED')
       write(27    )              &
#  endif
         sngl(TIME                 ),& ! 1
         sngl(vs_PI_Qmin           ),& ! 2
         sngl(vs_PI_Qmax           ),& ! 3
         sngl(vs_Qopt              ),& ! 4
         sngl(vs_Qoptmin * 1.02d0  ),& ! 5
         sngl(omega_set            ),& ! 6
         sngl(vp_PI2_pitmin        ),& ! 7 
         sngl(GS                   ),& ! 8
         sngl(GS_NL                ),& ! 9
         sngl(GS_IPC               ),& !10
         sngl(GS_TTa               ),& !11
         sngl(GenTrqDTDtot*gearctrl),& !12
         sngl(GenTrqDTD(1)*gearctrl),& !13
         sngl(GenTrqDTD(2)*gearctrl),& !14
         sngl(GenSpeed_LP1         ),& !15
         sngl(GenSpeed_LP2         ),& !16
         sngl(collpit_LP1          ),& !17
         sngl(collpit_LP2          ),& !18
         sngl(power_LP1            ),& !19
         sngl(power_LP2            ),& !20
         sngl(GenTrq_LP1           ),& !21
         sngl(GenTrq_LP2           )   !22
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
   return

!--- Rerun
 3 return

!--- Recall
 4 return

 100  format(1000e15.5)


 END Subroutine controller
!
!
!
!----------------------------------------------------------------------
 Subroutine controller_init (collpit, gentrq, array0, array1)
!----------------------------------------------------------------------

 Use dtu_we_controller_mod !subroutine init_regulation

   implicit none

   real(8), intent(in ) :: collpit, gentrq
   real(8), intent(out) :: array0(50), array1(100)

   real(8) ::  array2(1)
   integer :: i, ii, iostat


      array0 = 1.d10;
      array1 = 0.d0;

   open(88, file='controls.inp')
#  if CONTROLLER_TYPE == 33
                 read(88, *) !** title
#  endif
   do i = 1, 76
      ii=(i- 1)*(i-11)*(i-15)*(i-24)*(i-26)*(i-33)*(i-38)*(i-39)*(i-40)*(i-43)*(i-45)* &
#  if CONTROLLER_TYPE == 3
         (i-47)*(i-48)*(i-49)*(i-50)*(i-53)*(i-58)*(i-62)*(i-71)*(i-74)*(i-76)
#  elif CONTROLLER_TYPE == 33
         (i-47)*(i-48)*(i-49)*(i-50)*(i-53)*(i-58)*(i-62)*(i-71)*(i-73)*(i-76)
#  endif
      if (ii==0) read(88, *) !** comments...
      read(88, *, iostat=iostat) array1(i)
      if (iostat /= 0) then; write(*,*) ' *** ERROR *** Could not read variable ', i,' in controller input file'; stop; endif;
   enddo !i

      read(88, *) ![blank line]
   do i = 1, 39
      ii=(i-1)*(i-2)*(i-9)*(i-16)*(i-26)*(i-35)
      if (ii==0) read(88, *) !** comments...
      read(88, *, iostat=iostat) array0(i)
      if (iostat /= 0) then; write(*,*) ' *** ERROR *** Could not read variable ', i,' in controller input file'; stop; endif;
   enddo !i

   close(88)
!qq
      array1( 99) = gentrq
      array1(100) = collpit

#  if CONTROLLER_TYPE == 3
      call init_regulation         (array1, array2)
#  elif CONTROLLER_TYPE == 33
      call init_regulation_advanced(array1, array2)
#  endif


 END Subroutine controller_init
