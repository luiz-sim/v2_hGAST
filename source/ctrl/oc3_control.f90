!-- Interface for the oc3 controller. Sets the communicating variables between gast and discon.
!-- A drive train damper was found to be important for the correct work of the controller.
!-- A pitch actuator may be also introduced, but at the moment the pitch acceleration is calculated
!-- inside the controller (modification) with time derivative of the pitch rating.

!-- controls2.f90 is the original file [1.double precision, 2.pitch acceleration]
!-- controls3.f90 is modiffied for taking into account the iterative process of gast.

!----------------------------------------------------------------------
 Subroutine controller (NTIME,it,ISTEP)
!----------------------------------------------------------------------

 Use Cbeam
 Use Craft

   implicit none

   integer, intent(in) :: NTIME, it, ISTEP
   integer      :: nbod,NQACC, i, istat   , isat
   real(8),save :: GenSpeed, dGenSpeed, GenTrq_out, GenTrq_D_out 
   real(8),save :: BlPitch(3), PitAng_out(3), PitVel_out(3), PitAcc_out(3)
   real(8)      :: QTC , QTC1 , QTC2
   real(8)      :: QTCP, QTCP1, QTCP2
   real(8)      :: VAR , DVAR , W    , D, TVLAG ,QVGAIN, QVFREQ, QVDAMP
   real(8)      :: MinVAR, MaxVAR, MinDVAR, MaxDVAR, MinDDVAR, MaxDDVAR
   real(8)      :: Q0M1, Q0M2 , Q0MLAG
   real(8)      :: UPRE, U1PRE
   integer      :: ieq , ieq2
   real(8),save :: array1(100), array2(100)
!--- Generator lag
   real(8)      :: Gen_lag
!--- Pitch actuator
   real(8)      :: freq_act     , damp_act     , lag_act
   real(8)      :: MinPitAng_act, MaxPitAng_act
   real(8)      :: MaxPitRat_act, MaxPitAcc_act
!--- Drive Train Damper
   real(8)      :: DTD1_gain, DTD1_freq, DTD1_damp
   real(8)      :: DTD2_gain, DTD2_freq, DTD2_damp
   integer      :: Itype_Act  !0: No, 1:1st Angle, 2:1st Rate, 3:2nd Angle, 4:2nd Angle & Rate
   real(8)      :: velF
!--- Notch/Low-pass filters
   real(8)      :: W1, W2, Z1, Z2, DDVAR


!--- Define parameters
   Itype_Act     =  3 !0: No, 1:1st Angle,  2:1st Rate,  3:2nd Angle
   Gen_lag       =  0.02d0         !generator lag                          [sec]
   freq_act      =  1.0d0* PI2     !frequency         of pitch actuator    [rad/s] 1.6d0*PI2
   damp_act      =  0.7d0          !damping factor    of pitch actuator    [-]     0.8d0
   lag_act       =  0.2d0
   MinPitAng_act = -5.d0 / R2D     !lower angle limit of pitch actuator    [rad]
   MaxPitAng_act = 90.d0 / R2D     !upper angle limit of pitch actuator    [rad]
   MaxPitRat_act =  8.d0 / R2D     !speed       limit of pitch actuator    [rad/s]
   MaxPitAcc_act =  8.d0 / R2D     !accelerationlimit of pitch actuator    [rad/s^2]

!from 10MW
!! Gen_lag       =  0.02d0         !generator lag                          [sec]
!! freq_act      =  1.6d0*PI2      !frequency         of pitch actuator    [rad/s] ### 1  Hz
!! damp_act      =  0.8d0          !damping factor    of pitch actuator    [-]     ### 0.7
!! MinPitAng_act = -5.d0 / R2D     !lower angle limit of pitch actuator    [rad]
!! MaxPitAng_act = 90.d0 / R2D     !upper angle limit of pitch actuator    [rad]
!! MaxPitRat_act = 10.d0 / R2D     !speed       limit of pitch actuator    [rad/s]
!! MaxPitAcc_act = 15.d0 / R2D     !accelerationlimit of pitch actuator    [rad/s^2]

!qq
!--- NREL5MW
   DTD1_gain     =  5.000d2
   DTD1_freq     =  1.700d0 !Hz free-free
   DTD1_damp     =  0.100d0
   DTD2_gain     =  0.000d6
   DTD2_freq     = 10.000d0 !Hz
   DTD2_damp     =  0.100d0
!--- V112
!  DTD1_gain     =200.d0
!  DTD1_freq     =  1.800d0 !Hz free-free
!  DTD1_damp     =  0.100d0
!  DTD2_gain     =200.d0
!  DTD2_freq     =  2.550d0 !Hz
!  DTD2_damp     =  0.100d0
!--- E66
!  DTD1_gain     =  5.000d6
!  DTD1_freq     =  2.800d0 !Hz free-free
!  DTD1_damp     =  0.100d0
!  DTD2_gain     =  5.000d6
!  DTD2_freq     = 10.000d0 !Hz
!  DTD2_damp     =  0.100d0

!qq
!! MinPitAng_act = -1.d9 / R2D     !lower angle limit of pitch actuator    [rad]
!! MaxPitAng_act =  1.d9 / R2D     !upper angle limit of pitch actuator    [rad]
!! MaxPitRat_act =  1.d9 / R2D     !speed       limit of pitch actuator    [rad/s]
!! MaxPitAcc_act =  1.d9 / R2D     !accelerationlimit of pitch actuator    [rad/s^2]


   if (ICMOD == 0) return

   if (ICMOD == 2) then
       TGenLSS_el   = 0.d0
       TMLossLSS_el = 0.d0
       return
   endif

   goto (1,2,3,4), ISTEP


!--- main part
 1    GenSpeed      = -RAT_GEAR * UT1_el(NDFBT_el+NQSW)
     dGenSpeed      = -RAT_GEAR * UT2_el(NDFBT_el+NQSW)

!--- b_hGAST= -b_COLL -b_INDIV + b_IMB ==> b_COLL= -b_hGAST -b_INDIV + b_IMB
   do nbod  = 1, NBLADE_el
      NQACC = ACCU%NQSACC_el(nbod-1) 
      BlPitch(nbod) = - (UT_el (NDFBT_el+NQSP+NQACC+5)-PITCH_Imb_el(nbod) ) -PITCH_INDIV(nbod) -PITCH_accel_el
   enddo


!--- Initialization of the controller
      istat   = 1
   if ((NTIME<=1).and.(it==1)) &
      istat   = 0

!--------------------------------------------
!--- Initialization of control variables
!---    1. DTD1
!---    2. DTD2
!---    3. GEN LAG
!--- 4- 6. Pitch Actuator [1st Pitch Angle or 1st Pitch Rate or 2nd Pitch Angle]
!--- 7- 8. Notch filters
!---    9. Low-pass filter
!---10-12. Low pass wind
!--------------------------------------------
   if ((NTIME==1).and.(it==1)) then
      QTC_el  = 0.d0;   QTCP_el = 0.d0;   QTC0_el  = 0.d0;
      QTC1_el = 0.d0;   QTCP1_el= 0.d0;   QTC01_el = 0.d0;
      QTC2_el = 0.d0;   QTCP2_el= 0.d0;   QTC02_el = 0.d0;

     if     (Itype_Act ==  1)                      then !0: No, 1:1st Angle,  2:1st Rate,  3:2nd Angle
      QTCP1_el ( 4: 6) = BlPitch(1:3) !1st order pitch angle
      QTC1_el  ( 4: 6) = BlPitch(1:3)
     elseif (Itype_Act ==  2 .or. Itype_Act ==  3) then !0: No, 1:1st Angle,  2:1st Rate,  3:2nd Angle
      QTCP_el  ( 4: 6) = BlPitch(1:3) !1nd order pitch rate
      QTC_el   ( 4: 6) = BlPitch(1:3)
     endif
      QTC1_el  ( 7: 8) = GenSpeed;  QTCP1_el  ( 7: 8) = GenSpeed
      QTC_el   ( 9   ) = GenSpeed;  QTCP_el   ( 9   ) = GenSpeed
   endif


!  call lowpasswind ( YP_vel, YP1_vel, YP2_vel, Y_vel, Y1_vel, Y2_vel, velA, velF, DT_el )
   call lowpasswind ( QTCP_el(10:11), QTCP1_el(10:11), QTCP2_el(10:11), &
                      QTC_el (10:11), QTC1_el (10:11), QTC2_el (10:11), &
                      HUB_VEL(1)    , velF           , DT_el              )

!--- Notch filter 1 [spd err] 3P
   ieq      = 7
   W1       = 1.11d0 * PI2
   W2       = W1
   Z1       = 0.0d0
   Z2       = 0.2d0
   VAR      = -RAT_GEAR * UT_el (NDFBT_el+NQSW)          !-RAT_GEAR*UT_el (jg) - OMEGA_REF*TIME
   DVAR     =  GenSpeed                                  !-RAT_GEAR*UT1_el(jg) - OMEGA_REF
   DDVAR    = dGenSpeed                                  !-RAT_GEAR*UT2_el(jg)

   call filter_notch2 ( W1   , W2   , Z1   , Z2   , DT_el             , &
                        VAR  , DVAR , DDVAR                           , &
                        QTCP_el (ieq) , QTCP1_el (ieq), QTCP2_el (ieq), &
                        QTC_el  (ieq) , QTC1_el  (ieq),  QTC2_el (ieq)    )
   
!--- Notch filter 2 [spd err] free-free
   ieq      = 8
   ieq2     = 7
   W1       = 2.22d0 * PI2
   W2       = W1
   Z1       = 0.0d0
   Z2       = 0.2d0
   VAR      = QTC_el  (ieq2)
   DVAR     = QTC1_el (ieq2)
   DDVAR    = QTC2_el (ieq2)
   
   call filter_notch2 ( W1   , W2   , Z1   , Z2   , DT_el             , &
                        VAR  , DVAR , DDVAR                           , &
                        QTCP_el (ieq) , QTCP1_el (ieq), QTCP2_el (ieq), &
                        QTC_el  (ieq) , QTC1_el  (ieq),  QTC2_el (ieq)    )

!--- Low-pass filter [spd err]
   ieq      = 9
   ieq2     = 8
   W1       = 1.59d0 * PI2 !10rad/s
   Z1       = 1.0d0
   VAR      = QTC1_el (ieq2)
!qqVAR      = GenSpeed

   call filter_lowpass2 ( W1            , Z1            , DT_el         , &
                          VAR                                           , &
                          QTCP_el (ieq) , QTCP1_el (ieq), QTCP2_el (ieq), &
                          QTC_el  (ieq) , QTC1_el  (ieq),  QTC2_el (ieq)    )

!qqGenSpeed = QTC_el(ieq)
!--- 
   call CONTROLOS3 ( TIME      , GenSpeed  , BlPitch   , velF      , &
                     NBLADE_el , istat     , ICASE_el  , DT_el     , &
                     GenTrq_out, GenTrq_D_out                      , &
                     PitAng_out, PitVel_out, PitAcc_out                )

   if ((NTIME==1).and.(it==1)) then
      QTCP1_el ( 3   ) = GenTrq_out
      QTC1_el  ( 3   ) = GenTrq_out
   endif

!--- Drive Train Damper Freq1
   ieq          = 1
   QVGAIN       = DTD1_gain
   QVFREQ       = DTD1_freq * PI2
   QVDAMP       = DTD1_damp
   DVAR         = dGenSpeed

   call filter_bandpass2 ( QVGAIN      , QVFREQ       , QVDAMP       , &
                           DT_el       , DVAR         , Q0M1         , &
                           QTCP_el(ieq), QTCP1_el(ieq), QTCP2_el(ieq), &
                           QTC_el (ieq), QTC1_el (ieq), QTC2_el (ieq)    )

!--- Drive Train Damper Freq2
   ieq2         = 2
   QVGAIN       = DTD2_gain
   QVFREQ       = DTD2_freq * PI2
   QVDAMP       = DTD2_damp
   DVAR         = dGenSpeed

   call filter_bandpass2 ( QVGAIN       , QVFREQ        , QVDAMP        , &
                           DT_el        , DVAR          , Q0M2          , &
                           QTCP_el(ieq2), QTCP1_el(ieq2), QTCP2_el(ieq2), &
                           QTC_el (ieq2), QTC1_el (ieq2), QTC2_el (ieq2)    )

   TGenLSS_el   = ( GenTrq_out + QTC_el (ieq) + QTC_el (ieq2)) * RAT_GEAR
   TGenLSS_M_el =  (Q0M1+Q0M2)                                 * RAT_GEAR**2
!!!TGenLSS_D_el = GenTrq_D_out                                 * RAT_GEAR

!--- Generator's time lag for torque demand
   VAR          = GenTrq_out + QTC_el (ieq) + QTC_el (ieq2)
   ieq          = 3
   TVLAG        = Gen_lag

   call filter_lowpass1 ( TVLAG       , DT_el        , VAR          , &
                          QTCP_el(ieq), QTCP1_el(ieq), QTCP2_el(ieq), &
                          QTC_el (ieq), QTC1_el (ieq), QTC2_el (ieq), &
                          Q0MLAG                                        )

   TGenLSS_el   =   QTC1_el (ieq)                              * RAT_GEAR
   TGenLSS_M_el =  (Q0M1+Q0M2) * Q0MLAG                        * RAT_GEAR**2


 if     (Itype_Act == 0) then

!--- Set pitch angle, velocity and acceleration
   do nbod      = 1, NBLADE_el
      i         = NDFBT_el + NQSP + ACCU%NQSACC_el(nbod-1) + 5
      UT_el (i) = - PitAng_out (nbod) - PITCH_INDIV(nbod) + PITCH_Imb_el(nbod) - PITCH_accel_el
      UT1_el(i) = - PitVel_out (nbod)
      UT2_el(i) = - PitAcc_out (nbod)
   enddo !nbod

 elseif (Itype_Act == 1) then

!--- Pitch Actuator 1st order Pitch Angle demand
   do nbod      = 1, NBLADE_el
      i         = NDFBT_el + NQSP + ACCU%NQSACC_el(nbod-1) + 5
      ieq       = 4 + (nbod-1) !4-6
      TVLAG     = lag_act
      VAR       = PitAng_out(nbod) + PITCH_INDIV(nbod) + PITCH_accel_el
      MinVAR    =-1.d30
      MaxVAR    = 1.d30
      MinDVAR   = MinPitAng_act
      MaxDVAR   = MaxPitAng_act
      MinDDVAR  =-MaxPitRat_act
      MaxDDVAR  = MaxPitRat_act

      call filter_lowpass1   ( TVLAG       , DT_el        , VAR          , &
                               QTCP_el(ieq), QTCP1_el(ieq), QTCP2_el(ieq), &
                               QTC_el (ieq), QTC1_el (ieq), QTC2_el (ieq), &
                               Q0MLAG                                        )

      call Newmark_saturate1 ( DT_el                                     , &
                               MinVAR      , MaxVAR                      , &
                               MinDVAR     , MaxDVAR                     , &
                               MinDDVAR    , MaxDDVAR                    , &
                               QTCP_el(ieq), QTCP1_el(ieq), QTCP2_el(ieq), &
                               QTC_el (ieq), QTC1_el (ieq), QTC2_el (ieq)    )

      UT_el (i) = - QTC1_el(ieq) + PITCH_Imb_el(nbod)
      UT1_el(i) = - QTC2_el(ieq)
      UT2_el(i) = -(QTC2_el(ieq)-QTCP2_el(ieq))/DT_el
   enddo !nbod

 elseif (Itype_Act == 2) then

!--- Pitch Actuator 1st order Pitch Rate demand
   do nbod      = 1, NBLADE_el
      i         = NDFBT_el + NQSP + ACCU%NQSACC_el(nbod-1) + 5
      ieq       = 4 + (nbod-1) !4-6
      TVLAG     = lag_act
      VAR       = PitVel_out(nbod) !!+ PITCH_INDIV(nbod) + PITCH_accel_el
      MinVAR    = MinPitAng_act
      MaxVAR    = MaxPitAng_act
      MinDVAR   =-MaxPitRat_act
      MaxDVAR   = MaxPitRat_act
      MinDDVAR  =-MaxPitAcc_act
      MaxDDVAR  = MaxPitAcc_act

      call filter_lowpass1   ( TVLAG       , DT_el        , VAR          , &
                               QTCP_el(ieq), QTCP1_el(ieq), QTCP2_el(ieq), &
                               QTC_el (ieq), QTC1_el (ieq), QTC2_el (ieq), &
                               Q0MLAG                                        )

      call Newmark_saturate1 ( DT_el                                     , &
                               MinVAR      , MaxVAR                      , &
                               MinDVAR     , MaxDVAR                     , &
                               MinDDVAR    , MaxDDVAR                    , &
                               QTCP_el(ieq), QTCP1_el(ieq), QTCP2_el(ieq), &
                               QTC_el (ieq), QTC1_el (ieq), QTC2_el (ieq)    )

      UT_el (i) = - QTC_el (ieq) + PITCH_Imb_el(nbod)
      UT1_el(i) = - QTC1_el(ieq)
      UT2_el(i) = - QTC2_el(ieq)
   enddo !nbod

 elseif (Itype_Act == 3) then

!--- Pitch Actuator 2nd order Pitch Angle demand
   do nbod      = 1, NBLADE_el
      i         = NDFBT_el + NQSP + ACCU%NQSACC_el(nbod-1) + 5
      ieq       = 4 + (nbod-1) !4-6
      W         = freq_act     ![rad/s]
      D         = damp_act
      VAR       = PitAng_out(nbod) + PITCH_INDIV(nbod) + PITCH_accel_el
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

 endif !Itype_Act

   return

!--- Update for iterations
 2 call UPDATE_CONTROLOS3

!--- Write control variables for built-in controller
   open(27,file='controller_qtc.dat',access='append')
    write(27,10)TIME,(QTC_el(ieq),QTC1_el(ieq),QTC2_el(ieq),ieq=1,11)
   close(27)

10 format (50e16.7)
   return

!--- Rerun
 3 call CONTROLOS3 ( TIME      , GenSpeed  , BlPitch   , HUB_VEL(1), &
                     NBLADE_el , 3         , ICASE_el  , DT_el     , &
                     GenTrq_out, PitAng_out, PitVel_out, PitAcc_out    )
   return

!--- Recall
 4 call CONTROLOS3 ( TIME      , GenSpeed  , BlPitch   , HUB_VEL(1), &
                     NBLADE_el , 4         , ICASE_el  , DT_el     , &
                     GenTrq_out, PitAng_out, PitVel_out, PitAcc_out    )
   return


 END Subroutine controller
!----------------------------------------------------------------------
 Subroutine UPDATE_CONTROLOS3
!----------------------------------------------------------------------

 Use Cbeam
 Use Craft, ONLY : HUB_VEL

   implicit None

   integer      :: nbod, i 
   Real(8),save :: GenSpeed, BlPitch(3),GenTrq_out, PitAng_out(3), PitVel_out(3), PitAcc_out(3)
   Real(8),save :: Velhub



   do nbod          = 1, NBLADE_el
      i             = NDFBT_el + NQSP + ACCU%NQSACC_el(nbod-1) + 5
      BlPitch(nbod) =-UT_el (i)
   enddo

   call CONTROLOS3 ( TIME      , GenSpeed  , BlPitch   , HUB_VEL(1), &
                     NBLADE_el , 2         , ICASE_el  , DT_el     , &
                     GenTrq_out, PitAng_out, PitVel_out, PitAcc_out    )


 END Subroutine UPDATE_CONTROLOS3

!include './ctrl/controls3_nrel.f90'
 include './ctrl/controls3_nrel2.f90'
!include './ctrl/controls3_v112.f90'
!include './ctrl/controls3_E66.f90'
