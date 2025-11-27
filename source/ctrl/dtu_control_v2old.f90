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

 Use Cbeam
 Use Craft, ONLY : HUB_VEL
 Use dtu_we_controller_mod

   implicit none

   integer, intent(in) :: NTIME, it, ISTEP
   integer       :: nbod,NQACC, i
   real(8), save :: GenSpeedLSS, dGenSpeedLSS, BlPitch(3),GenTrqLSS, PitAngDem(3)
   real(8)       :: QTC , QTC1 , QTC2
   real(8)       :: QTCP, QTCP1, QTCP2
   real(8)       :: VAR , DVAR , W    , D, TVLAG ,QVGAIN, QVFREQ, QVDAMP
   real(8)       :: DDVAR , W1, W2, D1, D2
   real(8)       :: MinVAR, MaxVAR, MinDVAR, MaxDVAR, MinDDVAR, MaxDDVAR
   real(8)       :: UPRE, U1PRE
   integer       :: ieq
   real(8), save :: array1(100), array2(100)
   real(8)       :: gentrq
   real(8), save :: collpit !save for output
   real(8)       :: array0( 50)
!--- Generator lag
   real(8), save :: Gen_lag, Q0MLAG
!--- Pitch actuator
   integer, save :: itypePit_act
   real(8), save :: freqPit_act  , dampPit_act
   real(8), save :: MinPitAng_act, MaxPitAng_act
   real(8), save :: MaxPitRat_act, MaxPitAcc_act
!--- Drive Train Damper
   real(8), save, dimension(2) :: DTD_gain, DTD_freq, DTD_damp, Q0M
   real(8), save               :: MaxGenTrqDTD
   real(8)                     :: GenTrqDTD
!--- Tower acceleration
   real(8), dimension(3) :: x,dx,ddx,dxl,ddxl, AKSI0
   integer               :: nbsub,nb_el,nel
!--- filters
   real(8)               :: varf(4), timef(4), dummy
!--- Gain Schedule
   real(8), save         :: GS, GS_NL, GS_IPC, GS_TTa, kk1, kk2, GS_NL_lim, GS_IPC0, kk1_inv, kk2_inv
!--- filtered vars LP1
   real(8)               :: GenSpeed_LP1, collpit_LP1, power_LP1, Ux_LP1, GenSpeed_LP2, GenSpeed_Notch
   real(8), save         :: power_rated
!--- IPC parameters
   integer, save         :: ipc_nhT      , ipc_ntime_init
   real(8), save         :: ipc_PitchGain, ipc_MaxPitchAmp
   real(8), save         :: ipc_FlapGain , ipc_MaxFlapAmp
   real(8), save         :: ipc_Ki(2)    , ipc_Kp(2)
!--- TTa parameters
   integer, save         ::                tta_ntime_init
   real(8), save         :: tta_Gengain  , tta_MaxGenAmp
   real(8), save         :: tta_Pitchgain, tta_MaxPitchAmp
   real(8), save         :: tta_Flapgain , tta_MaxFlapAmp
   real(8), save         :: tta_Ki       , tta_Kp
!--- flap Actuator
   integer, save         :: itypeFlap_act
   real(8), save         :: timeFlap_act , MinFlapAng_act, MaxFlapAng_act, MaxFlapRat_act 


   if (ICMOD == 0) return

   if (ICMOD == 2) then
       TGenLSS_el   = 0.d0
       TMLossLSS_el = 0.d0
       return
   endif

   goto (1,2,3,4), ISTEP

!--- main part
 1 if (it>1) return

      GenSpeedLSS   = -UT1_el(NDFBT_el+NQSW)
     dGenSpeedLSS   = -UT2_el(NDFBT_el+NQSW)

!--- b_hGAST= -b_COLL -b_INDIV + b_IMB ==> b_COLL= -b_hGAST -b_INDIV + b_IMB
      collpit = 0.d0
   do nbod    = 1, NBLADE_el
      NQACC   = ACCU%NQSACC_el(nbod-1) 
      BlPitch(nbod) = - (UT_el (NDFBT_el+NQSP+NQACC+5)-PITCH_Imb_el(nbod) ) -PITCH_INDIV(nbod) -PITCH_accel_el
      collpit = collpit + BlPitch(nbod) / dble(NBLADE_el)
   enddo

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


   if ((NTIME==1).and.(it==1)) then
!--- Initialization of the 10Mw controller
      gentrq  = TREF
      call controller_init (collpit, gentrq, array0, array1)

!------ Define parameters
      write(*,*)'array0'
      do i=1,23
      write(*,*)array0(i)
      enddo

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
!------ inputs
      power_rated     = array1( 1) * 1000.d0
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
   endif


!--- Main controller call, pass the communication variables
   array1( 1) = TIME            !  1: general time                            [s]
   array1( 2) = GenSpeedLSS     !  2: constraint bearing1 shaft_rot 1 only 2  [rad/s] Generator speed (Default LSS, if HSS insert gear ratio in input #76)
   array1( 3) = BlPitch(1)      !  3: constraint bearing2 pitch1 1 only 1     [rad]
   array1( 4) = BlPitch(2)      !  4: constraint bearing2 pitch2 1 only 1     [rad]
   array1( 5) = BlPitch(3)      !  5: constraint bearing2 pitch3 1 only 1     [rad]
   array1( 6) = HUB_VEL(1)      !  6: wind free_wind 1 0.0 0.0 hub height     [m/s] global coords at hub height
   array1( 7) = HUB_VEL(2)      !  7: wind free_wind 1 0.0 0.0 hub height     [m/s] global coords at hub height
   array1( 8) = HUB_VEL(3)      !  8: wind free_wind 1 0.0 0.0 hub height     [m/s] global coords at hub height
   array1( 9) = array2 (5)      !  9: mech. power  ; [W]
   array1(10) = 0               ! 10: grid flag  ; [1=no grid,0=grid]
!qq local or global??
   array1(11) = ddx    (1)      ! 11: Tower top x-acceleration  ; [m/s^2]                                                                                   
   array1(12) = ddx    (2)      ! 12: Tower top y-acceleration  ; [m/s^2]                                                                                   

   call update_regulation(array1, array2)

   GenTrqLSS      = array2(1  )
   PitAngDem(1:3) = array2(2:4)

   if ((NTIME==1).and.(it==1)) then
!--------------------------------------------
!--- Initialization of control variables                                  current configuration
!---    1. DTD1                                               [QTC ]      free-free 1st
!---    2. DTD2                                               [QTC ]      free-free 2nd
!---    3. GEN LAG                                            [QTC1]      ~0.02s
!--- 4- 6. Pitch Actuator  2nd order                          [QTC ]      ~1.6Hz   , z=0.8
!--- 7-10. LOWPASS-1st (genspeed,collpitch,power,UwindX)      [QTC1]      Tperiod
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

!--[15-18. Ttop accel  (notch1P 3P,lowpass1P 2.5P,notch2P 3P,lowpass2P 2.5P)]
!--------------------------------------------
      QTC_el  = 0.d0;   QTCP_el = 0.d0;   QTC0_el  = 0.d0;
      QTC1_el = 0.d0;   QTCP1_el= 0.d0;   QTC01_el = 0.d0;
      QTC2_el = 0.d0;   QTCP2_el= 0.d0;   QTC02_el = 0.d0;

!DTD  QTC_el  (1:2) = GenTrqLSS               ;  QTCP_el  (1:2) = GenTrqLSS               !BP2
      QTC1_el (3  ) = GenTrqLSS               ;  QTCP1_el (3  ) = GenTrqLSS               !LP1
      QTC_el  (4:6) = BlPitch(1:3)            ;  QTCP_el  (4:6) = BlPitch(1:3)            !Actuator
      QTC1_el (7  ) = GenSpeedLSS             ;  QTCP1_el (7  ) = GenSpeedLSS             !LP1
      QTC1_el (8  ) = collpit                 ;  QTCP1_el (8  ) = collpit                 !LP1
      QTC1_el (9  ) = GenTrqLSS * GenSpeedLSS ;  QTCP1_el (9  ) = GenTrqLSS * GenSpeedLSS !LP1
      QTC1_el (10 ) = HUB_VEL(1)              ;  QTCP1_el (10 ) = HUB_VEL(1)              !LP1
      QTC_el  (11 ) = GenSpeedLSS             ;  QTCP_el  (11 ) = GenSpeedLSS             !LP2
      QTC_el  (12 ) = GenSpeedLSS             ;  QTCP_el  (12 ) = GenSpeedLSS             !notch
      QTC_el  (13 ) = GenTrqLSS * GenSpeedLSS ;  QTCP_el  (13 ) = GenTrqLSS * GenSpeedLSS !LP2
      QTC_el  (14 ) = GenTrqLSS * GenSpeedLSS ;  QTCP_el  (14 ) = GenTrqLSS * GenSpeedLSS !notch
   endif

!--- Drive Train Damper 1st and 2nd free-free modes
      GenTrqDTD    = 0.d0
   do i            = 1, 2
      ieq          = i           !1-2
      QVGAIN       = DTD_gain(i)
      QVFREQ       = DTD_freq(i) ![rad/s]
      QVDAMP       = DTD_damp(i)
      VAR          = dGenSpeedLSS
      
      call filter_bandpass2 ( QVGAIN      , QVFREQ       , QVDAMP       , &
                              DT_el       , VAR          , Q0M(i)       , &
                              QTCP_el(ieq), QTCP1_el(ieq), QTCP2_el(ieq), &
                              QTC_el (ieq), QTC1_el (ieq), QTC2_el (ieq)    )

      GenTrqDTD    = GenTrqDTD + QTC_el(ieq) 
   enddo !i
      GenTrqDTD    = min ( max( GenTrqDTD, -MaxGenTrqDTD  ), MaxGenTrqDTD )


!--- 1st order lowpass filter for genspeed, collpitch, power and UwindX
      varf(1)   = GenSpeedLSS             ; timef(1) = TPERIOD_el !ieq=7
      varf(2)   = collpit                 ; timef(2) = TPERIOD_el !ieq=8
      varf(3)   = GenTrqLSS * GenSpeedLSS ; timef(3) = TPERIOD_el !ieq=9
      varf(4)   = HUB_VEL(1)              ; timef(4) = TPERIOD_el !ieq=10 !sqrt(HUB_VEL(1)**2+HUB_VEL(2)**2)
   do i         = 1, 4
      ieq       = 7 + i-1 !7-10
      VAR       = varf (i)
      TVLAG     = timef(i)
      
      call filter_lowpass1 ( TVLAG       , DT_el        , VAR          , &
                             QTCP_el(ieq), QTCP1_el(ieq), QTCP2_el(ieq), &
                             QTC_el (ieq), QTC1_el (ieq), QTC2_el (ieq), &
                             dummy                                         )
   enddo !i


!------  Equation for the low-pass filter: w^2 / (s^2      + 2z w s + w^2)
      ieq       = 11
      W         = 0.4d0  * PI2  ![rad/s]
      D         = 0.7d0
      VAR       = GenSpeedLSS
      
      call filter_lowpass2 ( W             , D             , DT_el         , &
                             VAR                                           , &
                             QTCP_el (ieq) , QTCP1_el (ieq), QTCP2_el (ieq), &
                             QTC_el  (ieq) , QTC1_el  (ieq),  QTC2_el (ieq)    )

!------  Equation for the notch filter: (s^2 + 2z1w1 s + w1^2) / (s^2 + 2z2w2 s + w2^2)
      ieq       = 12
      W1        = 1.80d0 * PI2  ![rad/s]
      W2        = W1
      D1        = 0.001d0
      D2        = 0.100d0
      VAR       = QTC_el (11) !-UT_el (NDFBT_el+NQSW) - (TPERIOD_el/PI2)*TIME  !OMEGA_REF*TIME
      DVAR      = QTC1_el(11) !-UT1_el(NDFBT_el+NQSW) - (TPERIOD_el/PI2)       !OMEGA_REF
      DDVAR     = QTC2_el(11) !-UT2_el(NDFBT_el+NQSW)
   
      call filter_notch2 ( W1   , W2   , D1   , D2   , DT_el             , &
                           VAR  , DVAR , DDVAR                           , &
                           QTCP_el (ieq) , QTCP1_el (ieq), QTCP2_el (ieq), &
                           QTC_el  (ieq) , QTC1_el  (ieq),  QTC2_el (ieq)    )


!------  Equation for the low-pass filter: w^2 / (s^2      + 2z w s + w^2)
      ieq       = 13
      W         = PI2 / TPERIOD_el  ![rad/s]
      D         = 0.7d0
      VAR       = GenTrqLSS * GenSpeedLSS !+ GenTrqDTD !Torque Demand + DTD
      
      call filter_lowpass2 ( W             , D             , DT_el         , &
                             VAR                                           , &
                             QTCP_el (ieq) , QTCP1_el (ieq), QTCP2_el (ieq), &
                             QTC_el  (ieq) , QTC1_el  (ieq),  QTC2_el (ieq)    )

!------  Equation for the notch filter: (s^2 + 2z1w1 s + w1^2) / (s^2 + 2z2w2 s + w2^2)
      ieq       = 14
      W1        = 1.80d0 * PI2  ![rad/s]
      W2        = W1
      D1        = 0.001d0
      D2        = 0.100d0
      VAR       = QTC_el (13) !-UT_el (NDFBT_el+NQSW) - (TPERIOD_el/PI2)*TIME  !OMEGA_REF*TIME
      DVAR      = QTC1_el(13) !-UT1_el(NDFBT_el+NQSW) - (TPERIOD_el/PI2)       !OMEGA_REF
      DDVAR     = QTC2_el(13) !-UT2_el(NDFBT_el+NQSW)
   
      call filter_notch2 ( W1   , W2   , D1   , D2   , DT_el             , &
                           VAR  , DVAR , DDVAR                           , &
                           QTCP_el (ieq) , QTCP1_el (ieq), QTCP2_el (ieq), &
                           QTC_el  (ieq) , QTC1_el  (ieq),  QTC2_el (ieq)    )

!--- filtered vars
   GenSpeed_LP1   = QTC1_el( 7)
   collpit_LP1    = QTC1_el( 8)
   power_LP1      = QTC1_el( 9)
   Ux_LP1         = QTC1_el(10)
   GenSpeed_LP2   = QTC_el (11)
   GenSpeed_Notch = QTC_el (12)
!  power_LP2      = QTC_el (13)
!  power_Notch    = QTC_el (14)


!--- Gain Schedule
   GS    = 1.d0/(1.d0 + collpit_LP1 * kk1_inv + collpit_LP1**2 * kk2_inv)

!-- Nonlinear gain to avoid large rotor speed excursion
   GS_NL = 1.d0
!! if (GS_NL_lim > 1.d0 ) &
!! GS_NL = (1.d0 + GenSpeedFiltErr**2 / (GenSpeedRef_full*(GS_NL_lim - 1.d0))**2)
   GS    = GS * GS_NL

      GS_IPC0 = 0.01d0
   if     (power_LP1<=0.8d0*power_rated) then
      GS_IPC  = GS_IPC0
   elseif (power_LP1>=      power_rated) then
      GS_IPC  = 1.d0
   else
      GS_IPC  = (power_LP1-0.8d0*power_rated)/(0.2d0*power_rated) * (1.d0-GS_IPC0) + GS_IPC0
   endif
      GS_IPC  = GS_IPC * GS
      GS_TTa  = 1.d0

!qq Pending
!----- Power controller ---
!--- Torque limits
!--- PID torque
!--- PID pitch

!--- IPC
   call IPC_control  (GS_IPC, ISTEP)

!--- TTACC (actuator, coll pitch, gen, coll flap)
   call TTACC_control(GS_TTa, ISTEP)


!--- Generator's time lag for torque demand
   ieq          = 3
   VAR          = GenTrqLSS + GenTrqDTD + TGenLSS_accel_el !Torque Demand + DTD
   TVLAG        = Gen_lag

   call filter_lowpass1 ( TVLAG       , DT_el        , VAR          , &
                          QTCP_el(ieq), QTCP1_el(ieq), QTCP2_el(ieq), &
                          QTC_el (ieq), QTC1_el (ieq), QTC2_el (ieq), &
                          Q0MLAG                                        )

   TGenLSS_el   =  QTC1_el (ieq)           ! * RAT_GEAR
   TGenLSS_M_el = (Q0M(1)+Q0M(2)) * Q0MLAG ! * RAT_GEAR**2
!!!TGenLSS_D_el = GenTrq_D_out             ! * RAT_GEAR


!--- Pitch Actuator: 2nd order (Pitch Angle demand)
   do nbod      = 1, NBLADE_el
      i         = NDFBT_el + NQSP + ACCU%NQSACC_el(nbod-1) + 5
      ieq       = 4 + (nbod-1) !4-6
      W         = freqPit_act      ![rad/s]
      D         = dampPit_act
      VAR       = PitAngDem (nbod) + PITCH_INDIV(nbod) + PITCH_accel_el
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
      VAR       = FLAP_INDIV (nbod) + BetaF_accel_el
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


   return

!--- Update for iterations [output results]
#  if   ASCII == 1
 2    open(27,file='controller.dat',access='append')
       write(27,100)              &
#  elif ASCII == 0
 2    open(27,file='controller.bin',access='append',form='UNFORMATTED')
       write(27    )              &
#  endif
         sngl(TIME    )          ,&
        (sngl(array2(i)),i=1,37)
      close(27)

#  if   ASCII == 1
      open(27,file='controller_qtc0.dat',access='append')
       write(27,100)                  &
#  elif ASCII == 0
      open(27,file='controller_qtc0.bin',access='append',form='UNFORMATTED')
       write(27    )                  &
#  endif
         sngl(TIME        )          ,&
        (sngl(QTC_el (ieq)),ieq=1,32),&
         sngl(GS)                    ,&
         sngl(GS_IPC)                ,&
         sngl(collpit)
      close(27)

#  if   ASCII == 1
      open(27,file='controller_qtc1.dat',access='append')
       write(27,100)                  &
#  elif ASCII == 0
      open(27,file='controller_qtc1.bin',access='append',form='UNFORMATTED')
       write(27    )                  &
#  endif
         sngl(TIME         )         ,&
        (sngl(QTC1_el (ieq)),ieq=1,32)
      close(27)

#  if   ASCII == 1
      open(27,file='controller_qtc2.dat',access='append')
       write(27,100)                  &
#  elif ASCII == 0
      open(27,file='controller_qtc2.bin',access='append',form='UNFORMATTED')
       write(27    )                  &
#  endif
         sngl(TIME         )         ,&
        (sngl(QTC2_el (ieq)),ieq=1,32)
      close(27)

!--- for Output
   call IPC_control  (GS_IPC, ISTEP)
   call TTACC_control(GS_TTa, ISTEP)


10 format (50e16.7)
   return

!--- Rerun
 3 return

!--- Recall
 4 return

 100  format(1000e15.5)

 contains
!----------------------------------------------------------------------
 Subroutine IPC_control (GS_IPC, ISTEP)
!----------------------------------------------------------------------
       
 Use Cbeam

   implicit none

!  integer, save       :: ipc_nhT      , ipc_ntime_init
!  real(8), save       :: ipc_PitchGain, ipc_MaxPitchAmp
!  real(8), save       :: ipc_FlapGain , ipc_MaxFlapAmp   !for flap, unused atm
!  real(8), save       :: ipc_Ki(2)    , ipc_Kp(2)

   integer, intent(in) :: ISTEP
   real(8), intent(in) :: GS_IPC
   integer, parameter  :: NH = 2 !Number of harmonics for flap motion
   integer, parameter  :: NB = 3 !Number of blades

   integer :: h, nn, nbod, ieq1, ieq2, ieq3
   real(8) :: GAIN   , CONP  , CONI
   real(8) :: VAR    , DVAR  , DDVAR
   real(8) :: W      , D     , dummy, W1,W2,D1,D2

   real(8)       :: theta   (NH,2), Bang(NBLADE_el)
   real(8), save :: M_NRT_C (NH,2)                  !save for output
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

      call PID2_CONTROL( GAIN,     DT_el                              , &
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
!
!  Tower top acceleration PI control through Gen and Coll Pitch or Flap
!
!----------------------------------------------------------------------
 Subroutine TTACC_control (GS_TTa, ISTEP)
       
   Use Cbeam

   implicit none

!  integer, save      ::                tta_ntime_init
!  real(8), save      :: tta_Gengain  , tta_MaxGenAmp
!  real(8), save      :: tta_Pitchgain, tta_MaxPitchAmp
!  real(8), save      :: tta_Flapgain , tta_MaxFlapAmp
!  real(8), save      :: tta_Ki       , tta_Kp

   integer, intent(in) :: ISTEP
   real(8), intent(in) :: GS_TTa

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
   W                = 2.0d0 *  PI2                ![rad/s]
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
   TGenLSS_accel_el = 0.d0
   PITCH_accel_el   = 0.d0
   BetaF_accel_el   = 0.d0

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

   TGenLSS_accel_el = min( max(-QTC1_el(ieq3) * tta_GenGain   , -tta_MaxGenAmp   ), tta_MaxGenAmp   )
   PITCH_accel_el   = min( max(-QTC1_el(ieq3) * tta_PitchGain , -tta_MaxPitchAmp ), tta_MaxPitchAmp )
   BetaF_accel_el   = min( max(-QTC1_el(ieq3) * tta_FlapGain  , -tta_MaxFlapAmp  ), tta_MaxFlapAmp  )

   return

#    if   ASCII == 1
 2    open(23,file='control_TTa.dat',access='append')
       write(23,'(1000e15.5)')           &
#    elif ASCII == 0
 2    open(23,file='control_TTa.bin',access='append',form='UNFORMATTED')
       write(23              )           &
#    endif
         sngl(TIME                    ) ,& !1
         sngl(TGenLSS_accel_el/1000.d0) ,& !2 [kNm  ] 
         sngl(PITCH_accel_el  *R2D    ) ,& !3 [deg  ]
         sngl(BetaF_accel_el  *R2D    ) ,& !4 [deg  ]
         sngl(accel1                  ) ,& !5 [m/s^2]
         sngl(accel2                  ) ,& !6 [m/s^2]
         sngl(accelF                  )    !7 [m/s^2]
      close(23)


 END Subroutine TTACC_control

 END Subroutine controller
!----------------------------------------------------------------------
 Subroutine controller_init (collpit, gentrq, array0, array1)
!----------------------------------------------------------------------

 Use dtu_we_controller_mod

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
