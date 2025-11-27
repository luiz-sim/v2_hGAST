!-- Interface for the 10MW DTU REFERENCE controller. Sets the communicating variables between gast and discon.
!-- The drive train damper is implemented inside the controller and at moment the gain is ZERO.
!-- The generator lag is introduced here, as well as
!-- the pitch actuator.

 include "./ctrl/10MW_controller_fcns.f90"
 include "./ctrl/10MW_controller.f90"
!----------------------------------------------------------------------
 Subroutine controller (NTIME,it,ISTEP)
!----------------------------------------------------------------------

 Use Cbeam
 Use Craft

   implicit none

   integer, intent(in) :: NTIME, it, ISTEP
   integer      :: nbod,NQACC, i, ictrl
   real(8),save :: GenSpeedLSS, BlPitch(3),GenTrqLSS, PitAngDem(3)
   real(8)      :: QTC , QTC1 , QTC2
   real(8)      :: QTCP, QTCP1, QTCP2
   real(8)      :: QTC0, QTC01, QTC02
   real(8)      :: VAR , DVAR , W    , D, TVLAG, Q0MLAG
   real(8)      :: MinVAR, MaxVAR, MinDVAR, MaxDVAR, MinDDVAR, MaxDDVAR
   real(8)      :: BITA, GAMMA
   real(8)      :: UPRE, U1PRE
   integer      :: ieq
   real(8),save :: array1(1000), array2(100)
   real(8) :: gentrq, collpit
!--- Generator lag
   real(8) :: Gen_lag
!--- Pitch actuator
   real(8) :: freq_act     , damp_act
   real(8) :: MinPitAng_act, MaxPitAng_act
   real(8) :: MaxPitRat_act, MaxPitAcc_act


!--- Define parameters
   Gen_lag       =  0.02d0         !generator lag                          [sec]
   freq_act      =  1.6d0*PI2      !frequency         of pitch actuator    [rad/s]
   damp_act      =  0.8d0          !damping factor    of pitch actuator    [-]
   MinPitAng_act = -5.d0 / R2D     !lower angle limit of pitch actuator    [rad]
   MaxPitAng_act = 90.d0 / R2D     !upper angle limit of pitch actuator    [rad]
   MaxPitRat_act = 10.d0 / R2D     !speed       limit of pitch actuator    [rad/s]
   MaxPitAcc_act = 15.d0 / R2D     !accelerationlimit of pitch actuator    [rad/s^2]
   BITA          = 0.25d0
   GAMMA         = 0.50d0


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

!--- b_hGAST= -b_COLL -b_INDIV + b_IMB ==> b_COLL= -b_hGAST -b_INDIV + b_IMB
   do nbod  = 1, NBLADE_el
      NQACC = ACCU%NQSACC_el(nbod-1) 
      BlPitch(nbod) = - (UT_el (NDFBT_el+NQSP+NQACC+5)-PITCH_Imb_el(nbod) ) -PITCH_INDIV(nbod) -PITCH_accel_el
   enddo

!--- Initialization of the 10Mw controller
   if ((NTIME<=1).and.(it==1)) then
      gentrq  = TREF
      collpit = (BlPitch(1)+BlPitch(2)+BlPitch(3))/3.d0
     open  (88, file='controls_v1.inp')
!     write(* ,*) 'select ctrl config 1[INN], 2 [AVA]'
      read (88,*) ictrl
     close (88)
      if     (ictrl==1) then; call controller_init_INN (collpit, gentrq) !INNWIND
      elseif (ictrl==2) then; call controller_init_AVA (collpit, gentrq) !AVAVAR
      else                  ; write(*,*)'Bad ctrl input in file controls_v1.inp';stop
      endif
   endif


!--- Main controller call, pass the communication variables
   array1(1) = TIME            ! time
   array1(2) = GenSpeedLSS     ! Generator speed, ref LSS
   array1(3) = BlPitch(1)      ! Pitch 1
   array1(4) = BlPitch(2)      ! Pitch 2
   array1(5) = BlPitch(3)      ! Pitch 3
   array1(6) = HUB_VEL(1)      ! Wind u
   array1(7) = HUB_VEL(2)      ! Wind v
   array1(8) = HUB_VEL(3)      ! Wind w

   call update_regulation(array1, array2)

   GenTrqLSS      = array2(1  )
   PitAngDem(1:3) = array2(2:4)

   if ((NTIME==1).and.(it==1)) then
      QTCP1_el (8   ) = GenTrqLSS
      QTC1_el  (8   ) = GenTrqLSS
      QTCP_el  (9:11) = BlPitch(1:3)
      QTC_el   (9:11) = BlPitch(1:3)
   endif


!--- Generator's time lag for torque demand
   ieq    = 8
   TVLAG  = Gen_lag
   VAR    = GenTrqLSS

   call TIME_LAG  ( TVLAG       , DT_el        , BITA         , GAMMA, &
                    -1.d30      , 1.d30                              , &
                    -1.d30      , 1.d30                              , &
                    -1.d30      , 1.d30                              , &
                    VAR                                              , &
                    QTCP_el(ieq), QTCP1_el(ieq), QTCP2_el(ieq)       , &
                    QTC_el (ieq), QTC1_el (ieq), QTC2_el (ieq)       , &
                    QTC0_el(ieq), QTC01_el(ieq), QTC02_el(ieq)       , &
                    Q0MLAG                                               )

   TGenLSS_el   = QTC1_el (ieq) !! * RAT_GEAR


!--- Pitch Actuator 2nd order Pitch Angle demand
   do nbod      = 1, NBLADE_el
      i         = NDFBT_el + NQSP + ACCU%NQSACC_el(nbod-1) + 5
      ieq       = 9 + (nbod-1) !9-11
      W         = freq_act      ![rad/s]
      D         = damp_act
      VAR       = PitAngDem (nbod) + PITCH_INDIV(nbod) + PITCH_accel_el
      DVAR      = 0.d0 !PitVel_out(nbod)
      MinVAR    = MinPitAng_act
      MaxVAR    = MaxPitAng_act
      MinDVAR   =-MaxPitRat_act
      MaxDVAR   = MaxPitRat_act
      MinDDVAR  =-MaxPitAcc_act
      MaxDDVAR  = MaxPitAcc_act

      call ACTUATOR_2 ( W, D        , DT_el        , BITA         , GAMMA, &
                        MinVAR      , MaxVAR                             , &
                        MinDVAR     , MaxDVAR                            , &
                        MinDDVAR    , MaxDDVAR                           , &
                        VAR         , DVAR                               , &
                        QTCP_el(ieq), QTCP1_el(ieq), QTCP2_el(ieq)       , &
                        QTC_el (ieq), QTC1_el (ieq), QTC2_el (ieq)       , &
                        QTC0_el(ieq), QTC01_el(ieq), QTC02_el(ieq)           )

      UT_el (i) = -QTC_el (ieq) + PITCH_Imb_el(nbod)
      UT1_el(i) = -QTC1_el(ieq)
      UT2_el(i) = -QTC2_el(ieq)
   enddo !nbod

   return

!--- Update for iterations [output results]
 2 open(27,file='controller.dat',access='append')
    write(27,10)TIME,(array2(i),i=1,25)
   close(27)
   open(27,file='controller_qtc.dat',access='append')
    write(27,10)TIME,(QTC_el(ieq),QTC1_el(ieq),QTC2_el(ieq),ieq=8,11)
   close(27)

10 format (50e16.7)
   return

!-- Rerun
 3 return

!-- Recall
 4 return


 END Subroutine controller
!----------------------------------------------------------------------
 Subroutine controller_init_INN (collpit, gentrq)
!----------------------------------------------------------------------
   implicit none

   real(8) :: collpit, gentrq
   real(8) :: array1(1000), rGearRatio, rEff


!--- Initialise demands to current values
!  collpit                           ! Blade pitch angle         [rad]
!  gentrqr                           ! Measured Generator Torque [Nm] LSS
!--- Gear Box
   rGearRatio = 50                   ! Gearbox ratio
   rEff       = 0.94        !0.944   ! Electrical efficiency
!--- prepare call to DTU code
   array1( 1) = 10000./rEff !10000   ! Rated power [kW]
   array1( 2) = 0.628       !0.733   ! Minimum rotor speed [rad/s] --> 6   rpm
   array1( 3) = 1.005                ! Rated   rotor speed [rad/s] --> 9.6 rpm
   array1( 4) = 15.6E6               ! Maximum allowable generator torque [Nm]
   array1( 5) = 0           !100     ! Minimum pitch angle, theta_min [deg],
                                     ! if |theta_min|>90, then a table of <wsp,theta_min> is read from a file named 'wptable.n', where n=int(theta_min)
   array1( 6) = 90                   ! Maximum pitch angle [deg]
   array1( 7) = 10                   ! Maximum pitch velocity operation [deg/s]
   array1( 8) = 0.2                  ! Frequency of generator speed filter [Hz]
   array1( 9) = 0.7                  ! Damping ratio of speed filter [-]
   array1(10) = 3.33 !1.83  !0.64    ! Frequency of free-free DT torsion mode [Hz], if zero no notch filter used
!--- Partial load control parameters
   array1(11) = 1.0013E7    !9.5e6   ! Optimal Cp tracking K factor [Nm/(rad/s)^2],
   array1(12) = 6.83E7      !7.33e7  ! Proportional gain of torque controller [Nm/(rad/s)]
   array1(13) = 1.53E7      !1.32e7  ! Integral gain of torque controller [Nm/rad]
   array1(14) = 0.0                  ! Differential gain of torque controller [Nm/(rad/s^2)]
!--- Full load control parameters
   array1(15) = 2                    ! Generator control switch [1=constant power, 2=constant torque]
   array1(16) = 0.524  !/ 3.d0      !0.592   ! Proportional gain of pitch controller [rad/(rad/s)]
   array1(17) = 0.141  !/10.d0      !0.133   ! Integral gain of pitch controller [rad/rad]
   array1(18) = 0.0                  ! Differential gain of pitch controller [rad/(rad/s^2)]
   array1(19) = 0.4e-8               ! Proportional power error gain [rad/W]
   array1(20) = 0.4e-8               ! Integral power error gain [rad/(Ws)]
   array1(21) = 198         !164.13  ! Coefficient of linear term in aerodynamic gain scheduling, KK1 [deg]
   array1(22) = 693         !702.09  ! Coefficient of quadratic term in aerodynamic gain scheduling, KK2 [deg^2]
                                     ! (if zero, KK1 = pitch angle at double gain)
   array1(23) = 1.3                  ! Relative speed for double nonlinear gain [-]
!--- Cut-in simulation parameters
   array1(24) = 0.0         !0.1     ! Cut-in time [s], if zero no cut-in simulated
   array1(25) = 4.0                  ! Time delay for soft start [1/1P]
!--- Cut-out simulation parameters
   array1(26) = 0.0         !710     ! Cut-out time [s], if zero no cut-out simulated
   array1(27) = 5.0                  ! Time constant for 1st order filter lag of torque cut-out [s]
   array1(28) = 1                    ! Stop type [1=linear two pitch speed stop, 2=exponential pitch speed stop]
   array1(29) = 1.0                  ! Time delay for pitch stop 1 [s]
   array1(30) = 20.0                 ! Maximum pitch velocity during stop 1 [deg/s]
   array1(31) = 1.0                  ! Time delay for pitch stop 2 [s]
   array1(32) = 10.0                 ! Maximum pitch velocity during stop 2 [deg/s]
!--- Expert parameters (keep default values unless otherwise given)
   array1(33) = 0.5                  ! Lower angle above lowest minimum pitch angle for switch [deg]
   array1(34) = 0.5                  ! Upper angle above lowest minimum pitch angle for switch [deg]
   array1(35) = 95.0                 ! Ratio between filtered and reference speed for fully open torque limits [%]
   array1(36) = 5.0                  ! Time constant of 1st order filter on wind speed used for minimum pitch [1/1P]
   array1(37) = 5.0                  ! Time constant of 1st order filter on pitch angle for gain scheduling [1/1P]
!--- Drive train damper
   array1(38) = 1.8e7                ! Proportional gain of DT damper [Nm/(rad/s)], requires frequency in input 10
!----------------------------------------------------------------------------------------------------------------------
!--- Over speed shutdown
!! array1(39) = 10.0    !avatar1500  ! Percent maximum over speed of generator speed before cut-out [%]
!--- Additional non-linear pitch control term
!! array1(40) = 0.0                  ! Err0 [rad/s]
!! array1(41) = 0.0                  ! ErrDot0 [rad/s^2]
!! array1(42) = 0.0                  ! PitNonLin1 [rad/s]

   call init_regulation (array1, gentrq, collpit)


 END Subroutine controller_init_INN
!----------------------------------------------------------------------
 Subroutine controller_init_AVA (collpit, gentrq)
!
!  AVATAR controller according to D4.6 setup (Boorsma e-mail 2016/3/11)
!----------------------------------------------------------------------
   implicit none

   real(8) :: collpit, gentrq
   real(8) :: array1(1000), rGearRatio, rEff


!--- Initialise demands to current values
!  collpit                           ! Blade pitch angle         [rad]
!  gentrqr                           ! Measured Generator Torque [Nm] LSS
!--- Gear Box
   rGearRatio = 50                   ! Gearbox ratio
   rEff       = 0.94                 ! Electrical efficiency
!--- prepare call to DTU code
   array1( 1) = 10000./rEff          ! Rated power [kW]
   array1( 2) = 0.628                ! Minimum rotor speed [rad/s] --> 6   rpm
   array1( 3) = 1.005                ! Rated   rotor speed [rad/s] --> 9.6 rpm
   array1( 4) = 15.6E6               ! Maximum allowable generator torque [Nm]
   array1( 5) = 0                    ! Minimum pitch angle, theta_min [deg],
                                     ! if |theta_min|>90, then a table of <wsp,theta_min> is read from a file named 'wptable.n', where n=int(theta_min)
   array1( 6) = 90                   ! Maximum pitch angle [deg]
   array1( 7) = 10                   ! Maximum pitch velocity operation [deg/s]
   array1( 8) = 0.4                  ! Frequency of generator speed filter [Hz]
   array1( 9) = 0.7                  ! Damping ratio of speed filter [-]
   array1(10) = 1.64                 ! Frequency of free-free DT torsion mode [Hz], if zero no notch filter used
!--- Partial load control parameters
   array1(11) = 9.65E6               ! Optimal Cp tracking K factor [Nm/(rad/s)^2],
   array1(12) = 1.05E8               ! Proportional gain of torque controller [Nm/(rad/s)]
   array1(13) = 1.53E7               ! Integral gain of torque controller [Nm/rad]
   array1(14) = 0.0                  ! Differential gain of torque controller [Nm/(rad/s^2)]
!--- Full load control parameters
   array1(15) = 1                    ! Generator control switch [1=constant power, 2=constant torque]
   array1(16) = 0.762                ! Proportional gain of pitch controller [rad/(rad/s)]
   array1(17) = 0.224                ! Integral gain of pitch controller [rad/rad]
   array1(18) = 0.0                  ! Differential gain of pitch controller [rad/(rad/s^2)]
   array1(19) = 0.4e-8               ! Proportional power error gain [rad/W]
   array1(20) = 0.4e-8               ! Integral power error gain [rad/(Ws)]
   array1(21) = 10.7                 ! Coefficient of linear term in aerodynamic gain scheduling, KK1 [deg]
   array1(22) = 601                  ! Coefficient of quadratic term in aerodynamic gain scheduling, KK2 [deg^2]
                                     ! (if zero, KK1 = pitch angle at double gain)
   array1(23) = 1.3                  ! Relative speed for double nonlinear gain [-]
!--- Cut-in simulation parameters
   array1(24) = 0.0                  ! Cut-in time [s], if zero no cut-in simulated
   array1(25) = 1.0                  ! Time delay for soft start [1/1P]
!--- Cut-out simulation parameters
   array1(26) = 0.0                  ! Cut-out time [s], if zero no cut-out simulated
   array1(27) = 5.0                  ! Time constant for 1st order filter lag of torque cut-out [s]
   array1(28) = 1                    ! Stop type [1=linear two pitch speed stop, 2=exponential pitch speed stop]
   array1(29) = 1.0                  ! Time delay for pitch stop 1 [s]
   array1(30) = 3.0                  ! Maximum pitch velocity during stop 1 [deg/s]
   array1(31) = 3.0                  ! Time delay for pitch stop 2 [s]
   array1(32) = 4.0                  ! Maximum pitch velocity during stop 2 [deg/s]
!--- Expert parameters (keep default values unless otherwise given)
   array1(33) = 2.0                  ! Lower angle above lowest minimum pitch angle for switch [deg]
   array1(34) = 2.0                  ! Upper angle above lowest minimum pitch angle for switch [deg]
   array1(35) = 95.0                 ! Ratio between filtered and reference speed for fully open torque limits [%]
   array1(36) = 2.0                  ! Time constant of 1st order filter on wind speed used for minimum pitch [1/1P]
   array1(37) = 1.0                  ! Time constant of 1st order filter on pitch angle for gain scheduling [1/1P]
!--- Drive train damper
   array1(38) = 0.0                  ! Proportional gain of DT damper [Nm/(rad/s)], requires frequency in input 10
!----------------------------------------------------------------------------------------------------------------------
!--- Over speed shutdown
!! array1(39) = 10.0    !avatar1500  ! Percent maximum over speed of generator speed before cut-out [%]
!--- Additional non-linear pitch control term
!! array1(40) = 0.0                  ! Err0 [rad/s]
!! array1(41) = 0.0                  ! ErrDot0 [rad/s^2]
!! array1(42) = 0.0                  ! PitNonLin1 [rad/s]

   call init_regulation (array1, gentrq, collpit)


 END Subroutine controller_init_AVA
