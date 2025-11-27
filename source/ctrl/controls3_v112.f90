!=======================================================================
Subroutine CONTROLOS3 ( Time      , GenSpeed  , BlPitch   , Velhub    , &
                        NumBl     , iStatus   , icase     , DT        , &
                        GenTrq_out, GenTrq_D_out                      , &
                        PitAng_out, PitVel_out, PitAcc_out                )

   ! This controller is used to implement a variable-speed generator-torque
   ! controller and PI collective blade pitch controller for the NREL 5 MW
   ! Offshore baseline wind turbine.

IMPLICIT NONE

   ! Input Vars
REAL(8), intent (in)    :: Time                                            ! Current simulation time, sec.
REAL(8), intent (in)    :: GenSpeed                                        ! Current  HSS (generator) speed, rad/s.
REAL(8), intent (in)    :: BlPitch   (3)                                   ! Current values of the blade pitch angles, rad.
REAL(8), intent (in)    :: Velhub                                          ! Current velocity at hubheight, m/s.
INTEGER, intent (in)    :: NumBl                                           ! Number of blades, (-).
INTEGER, intent (in)    :: iStatus                                         ! A status flag set by the simulation as follows: 0 if this is the first call, 1 for all subsequent time steps, -1 if this is the final call at the end of the simulation.
INTEGER, intent (in)    :: icase           !jim                            ! 0:onshore, 1:monopile, 2:jacket, 3:floating
REAL(8), intent (in)    :: DT                                              ! Time-step, sec.

   ! Output Vars
REAL(8), intent (out)   :: GenTrq_out                                      ! Communicate Generator's Torque Demand
REAL(8), intent (out)   :: GenTrq_D_out                                    !
REAL(8), intent (out)   :: PitAng_out(3)                                   ! Communicate Pitch angle
REAL(8), intent (out)   :: PitVel_out(3)                                   ! Communicate Pitch vel
REAL(8), intent (out)   :: PitAcc_out(3)                                   ! Communicate Pitch accel

   ! Local Variables:
INTEGER,      SAVE      :: iconst_pow = 1                                  ! flag for region 3 constant power (1) -DEFAULT- or constant torque (0), (-).
REAL(8), PARAMETER      :: R2D        = 180.d0 / dacos(-1.d0) !=57.29578d0 ! Factor to convert radians to degrees.
REAL(8), PARAMETER      :: RPS2RPM    =  30.d0 / dacos(-1.d0) != 9.54930d0 ! Factor to convert radians per second to revolutions per minute.
REAL(8)                 :: Alpha                                           ! Current coefficient in the recursive, single-pole, low-pass filter, (-).
REAL(8), PARAMETER      :: CornerFreq =   1.0d0 !1.570796d0                       !#Corner frequency (-3dB point) in the recursive, single-pole, low-pass filter, rad/s. -- chosen to be 1/4 the blade edgewise natural frequency ( 1/4 of approx. 1Hz = 0.25Hz = 1.570796rad/s)
REAL(8), SAVE           :: GenSpeedF                                       ! Filtered HSS (generator) speed, rad/s.
REAL(8),      SAVE      :: GenTrq, GenTrqD !jim add SAVE                   ! Electrical generator torque, N-m.
REAL(8)                 :: GK                                              ! Current value of the gain correction factor, used in the gain scheduling law of the pitch controller, (-).

REAL(8), SAVE           :: IntSpdErr                                       ! Current integral of speed error w.r.t. time, rad.
REAL(8), SAVE           :: LastGenTrq                                      ! Commanded electrical generator torque the last time the controller was called, N-m.
REAL(8),      SAVE      :: PC_KI!     =   0.008068634d0                    ! Integral gain for pitch controller at rated pitch (zero), (-).      0.0008965149d0 !REDUCE GAINS FOR HYWIND
REAL(8), PARAMETER      :: PC_KK      =   0.1099965d0  !     3.78deg     ??!#Pitch angle were the derivative of the aerodynamic power w.r.t. pitch has increased by a factor of two relative to the derivative at rated pitch (zero), rad.
REAL(8),      SAVE      :: PC_KP!     =   0.01882681d0                     ! Proportional gain for pitch controller at rated pitch (zero), sec.  0.006275604d0  !REDUCE GAINS FOR HYWIND
REAL(8), PARAMETER      :: MIN_Spd    =   5.75d0/RPS2RPM *113.2d0          !#Transitional generator speed (HSS side) between regions 1     and 1 1/2, rad/s.
REAL(8), PARAMETER      :: MAX_Spd    =  13.2d0 /RPS2RPM *113.2d0          !#Desired (reference) HSS speed for pitch controller, rad/s.
REAL(8), PARAMETER      :: PC_MaxPit  =  90.0d0 /R2D                       ! Maximum pitch setting in pitch controller, rad.
REAL(8), PARAMETER      :: PC_MaxRat  =   8.0d0 /R2D                       ! Maximum pitch  rate (in absolute value) in pitch  controller, rad/s.
REAL(8),       SAVE     :: PC_MinPit  =   0.0d0 /R2D                       ! Minimum pitch setting in pitch controller, rad.
REAL(8),       SAVE     :: PC_RefSpd  = MAX_Spd                            ! Desired (reference) HSS speed for pitch controller, rad/s.
REAL(8), SAVE           :: PitCom    (3)                                   ! Commanded pitch of each blade the last time the controller was called, rad.
REAL(8)                 :: PitComI                                         ! Integral term of command pitch, rad.
REAL(8)                 :: PitComP                                         ! Proportional term of command pitch, rad.
REAL(8)                 :: PitComT                                         ! Total command pitch based on the sum of the proportional and integral terms, rad.
REAL(8),       SAVE     :: PitRate   (3)  !jim add SAVE                    ! Pitch rates of each blade based on the current pitch angles and current pitch command, rad/s.
REAL(8)                 :: SpdErr                                          ! Current speed error, rad/s.
REAL(8)                 :: TrqRate                                         ! Torque rate based on the current and last torque commands, N-m/s.
REAL(8), PARAMETER      :: VS_RtPwr   = 3450000.d0/0.94d0                  !#Rated generator power in Region 3, Watts. -- chosen to be 5MW divided by the electrical generator efficiency of 94.4%
REAL(8), PARAMETER      :: VS_Rgn2K   =       0.00335d0                    !#Generator torque constant in Region 2 (HSS side), N-m/(rad/s)^2.
REAL(8), PARAMETER      :: VS_pow     =       3.0d0                        !#Generator torque omega power constant in Region 2 (Default is 2).
REAL(8), PARAMETER      :: VS_MaxTq   = VS_RtPwr/MAX_Spd * 1.1d0           ! 43093.55*1.1 ! Maximum generator torque in Region 3 (HSS side), N-m. -- chosen to be 10% above VS_RtTq = 43.09355kNm
REAL(8), PARAMETER      :: VS_MaxRat  =  VS_MaxTq * 0.3d0  !15000          ! Maximum torque rate (in absolute value) in torque controller, N-m/s.
REAL(8),      SAVE      :: VS_CtInSp  = MIN_Spd                            ! Transitional generator speed (HSS side) between regions 1     and 1 1/2, rad/s.
REAL(8),      SAVE      :: VS_Rgn2Sp  = MIN_Spd   * 1.26d0                 !#Transitional generator speed (HSS side) between regions 1 1/2 and 2    , rad/s.
REAL(8), SAVE           :: VS_TrGnSp                                       ! Transitional generator speed (HSS side) between regions 2     and 2 1/2, rad/s.
REAL(8), PARAMETER      :: VS_Rgn3MP  =   1.0d0 /R2D                       ! Minimum pitch angle at which the torque is computed as if we are in region 3 regardless of the generator speed, rad. -- chosen to be 1.0 degree above PC_MinPit
REAL(8), PARAMETER      :: VS_RtGnSp  = MAX_Spd   * 0.99d0                 ! Rated generator speed (HSS side), rad/s. -- chosen to be 99% of MAX_Spd
REAL(8), SAVE           :: VS_Slope15                                      ! Torque/speed slope of region 1 1/2 cut-in torque ramp , N-m/(rad/s).
REAL(8), SAVE           :: VS_Slope25                                      ! Torque/speed slope of region 2 1/2 induction generator, N-m/(rad/s).
REAL(8), PARAMETER      :: VS_SlPc    =   8.d0                             ! Rated generator slip percentage in Region 2 1/2, %. !jim: Defines VS_SySp = VS_RtGnSp/(1+VS_SlPc/100)
REAL(8), SAVE           :: VS_SySp    = VS_RtGnSp/( 1.d0+VS_SlPc/100.d0 )  ! Synchronous speed of region 2 1/2 induction generator, rad/s.

REAL(8), SAVE           :: LastPitComAve   !jim AVERAGED!!
REAL(8), SAVE           :: LastIntSpdErr   !jim
REAL(8),      SAVE      :: LastGenSpeedF   !jim
REAL(8), SAVE           :: LastPitRate(3)  !jim
REAL(8), SAVE           :: PitAcc     (3)  !jim
REAL(8),      SAVE      :: LastBlPitch(3)  !jim
REAL(8)                 :: Region          !jim output the control region

REAL(8),      SAVE      :: res_vs(9), res_pc(26)                           ! Store resutls for output
REAL(8),      SAVE      :: V(9), P(9)                                      ! velocity-pitch look-up tables, m/s, rad.
INTEGER                 :: I                                               ! Generic index.
INTEGER                 :: K                                               ! Loops through blades.
!--- Cut-out vars
real(8), save :: cutout_time1   =-50.d0 !Time to initiate cut-out (phase 1), sec
real(8), save :: cutout_time2   =  3.d0 !Time to start cut-out phase 2, sec
real(8), save :: cutout_timegen =  5.d0 !Time to perform the Torque linear reduction, sec
real(8), save :: cutout_pitvel1 =  3.d0 !Pitch velocity during cut-out phase 1, deg/s
real(8), save :: cutout_pitvel2 =  4.d0 !Pitch velocity during cut-out phase 2, deg/s
integer, save :: cutout_step    =  0    !cutout phase 0,1,2
real(8), save :: cutout_TrqRate         !cutout Torque rate for linear reduction
!--- Cut-in vars
real(8), save :: cutin_time1    =-30.d0 !Time to initiate cut-in  (phase 1)       , sec
real(8), save :: cutin_time2    =  5.d0 !Time to keep constant pitang1            , sec
real(8), save :: cutin_timer    =  1.d9 !Time at which pitch reaches cutin_pitang1, sec
real(8), save :: cutin_pitang1  = 30.d0 !Pitch angle to stop phase 1              , deg
real(8), save :: cutin_pitvel1  =  4.d0 !Pitch velocity during cut-out phase 1, deg/s
integer, save :: cutin_step     =  0    !cutin phase 0,1,2,3
!--- PID2
integer, save :: IPID2 = 1    !qq       !flag to enable the PID2 option
real(8), save :: error (1)
real(8), save :: errorp(1), errorpp(1)
!--- Newton
real(8)       :: func, dfunc
integer       :: ii
!--- Torque at Region2
real(8)       :: Ftorq, dFtorq, w


    Ftorq(w) =          VS_Rgn2K * w** VS_pow 
   dFtorq(w) = VS_pow * VS_Rgn2K * w**(VS_pow-1.d0)
!   Ftorq(w) =   3d-07*w**6 -   0.0002d0*w**5 +   0.0527d0*w**4 -   7.9161d0*w**3 +   664.59d0*w**2 - 29544d0*w + 543323d0
!  dFtorq(w) = 6*3d-07*w**5 - 5*0.0002d0*w**4 + 4*0.0527d0*w**3 - 3*7.9161d0*w**2 + 2*664.59d0*w    - 29544d0

!--- Initialization, Backup, Output and Update for next timestep
   if (iStatus /=1) then
      call init_back_out (iStatus)
      if (istatus >0) return
   endif


   ! Main control calculations:

!=======================================================================

   ! Filter the HSS (generator) speed measurement:
   ! ---------------------------------------------
   ! NOTE: This is a very simple recursive, single-pole, low-pass filter with
   !       exponential smoothing.

   Alpha     = dexp( -DT*CornerFreq )                                ! a = e^-dt*fc[rad/s]
   GenSpeedF = (1.d0 - Alpha)*GenSpeed + Alpha*LastGenSpeedF         ! Filtered Generator speed
!qqGenSpeedF = GenSpeed
   PC_RefSpd = MAX_Spd

!  call lint1  ( xa, ya, n1, n2, np, x     , y         )
   call lint1  ( V , P , 1 , 5 , 9 , Velhub, PC_MinPit )             ! PC_MinPit: Minimum pitch setting in pitch controller, rad.

!=======================================================================

   ! Controller modes:
   call Cut_in
   call Cut_out
   call Normal_oper

!=======================================================================

   ! Set communication vars:
   GenTrq_out          = GenTrq
   GenTrq_D_out        = GenTrqD
   PitAng_out(1:NumBl) = PitCom     (1:NumBl)
   PitVel_out(1:NumBl) = PitRate    (1:NumBl)
   PitAcc_out(1:NumBl) = PitAcc     (1:NumBl)


 Contains

!---------------------------------
 Subroutine Normal_oper
!---------------------------------

 implicit none

   real(8) :: output


   if (cutout_step > 0) return
   if (cutin_step  > 0) return

   ! Variable-speed torque control:
   ! ------------------------------

   ! Specify region and compute the generator torque:

   if ( GenSpeedF     >= VS_RtGnSp           .OR. &                           ! VS_RtGnSp  = MAX_Spd   * 0.99d0
        LastPitComAve >= PC_MinPit+VS_Rgn3MP .AND. cutin_time1 <= 0.d0 ) then ! We are in region 3 - (power or torque) is constant
      Region  = 3.0d0
     if     (iconst_pow == 0) then
      GenTrq  = VS_RtPwr/MAX_Spd                                              !    Region 3 Constant Torque
      GenTrqD = 0.0d0
     elseif (iconst_pow == 1) then
      GenTrq  = VS_RtPwr/GenSpeedF                                            !    Region 3 Constant Power
!qq
      GenTrqD =-VS_RtPwr/GenSpeedF**2                                         !    Region 3 Constant Power
     endif
   elseif ( GenSpeedF <= VS_CtInSp )  then                                    ! We are in region 1 - torque is zero
      Region  = 1.0d0
      GenTrq  = 0.0d0
      GenTrqD = 0.0d0
   elseif ( GenSpeedF <  VS_Rgn2Sp )  then                                    ! We are in region 1 1/2 - linear ramp in torque from zero to optimal
      Region  = 1.5d0
      GenTrq  = VS_Slope15*( GenSpeedF - VS_CtInSp )
      GenTrqD = VS_Slope15
   elseif ( GenSpeedF <  VS_TrGnSp )  then                                    ! We are in region 2 - optimal torque is proportional to the square of the generator speed
      Region  = 2.0d0
!     GenTrq  =        VS_Rgn2K*GenSpeedF**VS_pow
!     GenTrqD = VS_pow*VS_Rgn2K*GenSpeedF**(VS_pow-1.d0)
      GenTrq  =  Ftorq( GenSpeedF )
      GenTrqD = dFtorq( GenSpeedF )
   else                                                                       ! We are in region 2 1/2 - simple induction generator transition region
      Region  = 2.5d0
      GenTrq  = VS_Slope25*( GenSpeedF - VS_SySp   )
      GenTrqD = VS_Slope25
   endif

   ! Saturate GenTrq using value and rate limits

      GenTrq  = dmin1( GenTrq                      , VS_MaxTq  )  ! Saturate the commanded torque using the maximum torque limit
      if ( iStatus == 0 )  LastGenTrq = GenTrq                    ! Initialize the value of LastGenTrq on the first pass only
      TrqRate = ( GenTrq - LastGenTrq ) / DT                      ! Torque rate (unsaturated)
      TrqRate = dmin1( dmax1( TrqRate, -VS_MaxRat ), VS_MaxRat )  ! Saturate the torque rate using its maximum absolute value
      GenTrq  = LastGenTrq + TrqRate*DT                           ! Saturate the commanded torque using the torque rate limit

!------ VS output vars
      res_vs(1) = Time
      res_vs(2) = GenTrq
      res_vs(3) = TrqRate
      res_vs(4) = LastGenTrq
      res_vs(5) = GenSpeedF
      res_vs(6) = GenSpeed
      res_vs(7) = Region
      res_vs(8) = Alpha
      res_vs(9) = GenTrqD

!=======================================================================

   ! Pitch control:
   ! --------------

   ! Compute the gain scheduling correction factor based on the previously
   !   averaged commanded pitch angle:

         GK         = 1.d0/( 1.d0 + LastPitComAve/PC_KK )

   ! Compute the current speed error and its integral w.r.t. time; saturate the
   !   integral term using the pitch angle limits:

         SpdErr     = GenSpeedF     - PC_RefSpd                             ! Current speed error
         IntSpdErr  = LastIntSpdErr + SpdErr*DT                             ! Current integral of speed error w.r.t. time
         IntSpdErr  = MIN( MAX( IntSpdErr, PC_MinPit/( GK*PC_KI ) ), &
                                           PC_MaxPit/( GK*PC_KI )      )    ! Saturate the integral term using the pitch angle limits, converted to integral speed error limits

   ! Compute the pitch commands associated with the proportional and integral gains:

         PitComP    = GK*PC_KP*   SpdErr                                    ! Proportional term
         PitComI    = GK*PC_KI*IntSpdErr                                    ! Integral term (saturated)

   ! Superimpose the individual commands to get the total pitch command;
   !   saturate the overall command using the pitch angle limits:

         PitComT    = PitComP + PitComI                                     ! Overall command (unsaturated)
         PitComT    = MIN( MAX( PitComT, PC_MinPit ), PC_MaxPit )           ! Saturate the overall command using the pitch angle limits

   ! Saturate the overall commanded pitch using the pitch rate limit:
   ! NOTE: Since the current pitch angle may be different for each blade
   !       (depending on the type of actuator implemented in the structural
   !       dynamics model), this pitch rate limit calculation and the
   !       resulting overall pitch angle command may be different for each
   !       blade.

      do K          = 1, NumBl
         PitRate(K) = ( PitComT - LastBlPitch(K) )/DT                       ! Pitch rate of blade K (unsaturated)
         PitRate(K) = MIN( MAX( PitRate(K), -PC_MaxRat ), PC_MaxRat )       ! Saturate the pitch rate of blade K using its maximum absolute value
!!!      PitCom (K) = LastBlPitch(K) + PitRate(K)*DT                        ! Saturate the overall command of blade K using the pitch rate limit
         PitCom (K) =     BlPitch(K) + PitRate(K)*DT                        ! Saturate the overall command of blade K using the pitch rate limit
         PitAcc (K) = (PitRate(K)-LastPitRate(K))/DT       
      enddo !K

!qq PID2 implementation
 if (IPID2==1) then
     error  = SpdErr
     output = LastBlPitch(1)

     call PID2 ( 1, DT, error, errorp, errorpp, PC_KP, PC_KI, 0.d0, GK, &
                 PC_MinPit, PC_MaxPit, PC_MaxRat, output )

     PitCom (:) = output;
     PitRate(:) = 0.d0;
     PitAcc (:) = 0.d0;
  endif

!------ PC output vars
      res_pc( 1   ) = Time
      res_pc( 2: 4) = PitCom (1:3)
      res_pc( 5: 7) = PitRate(1:3)
      res_pc( 8:10) = PitAcc (1:3)
      res_pc(11   ) = GK
      res_pc(12   ) = SpdErr
      res_pc(13   ) = IntSpdErr
      res_pc(14   ) = PitComP
      res_pc(15   ) = PitComI
      res_pc(16   ) = PitComT
      res_pc(17:19) = BlPitch(1:3)
      res_pc(20:22) = LastBlPitch(1:3)
      res_pc(23   ) = GenSpeed
      res_pc(24   ) = Velhub
      res_pc(25   ) = PC_MinPit
      res_pc(26   ) = PC_RefSpd 


 END Subroutine Normal_oper
!---------------------------------
 Subroutine Cut_in
!---------------------------------

 implicit none


   if      (cutout_step > 0             ) return

   if      (cutin_time1 <= 0.d0 .or. &          ! cut-in will not simulated in this case
            time-cutin_timer>cutin_time2) then  ! goto normal control mode
      cutin_step     = 0
      return
   elseif  (time < cutin_time1          ) then  ! cut-in haven't started yet (still parked)
      cutin_step     = 1
      TrqRate        = 0.d0
      GenTrq         = 0.d0
      PitCom (:)     = LastBlPitch(:)
      PitRate(:)     = 0.d0
      PitAcc (:)     = 0.d0
      return
   elseif (BlPitch(1)> cutin_pitang1/R2D) then  ! cut-in phase I  ( Generator torque is zero
      cutin_step     = 2                        !                   pitch reduction with vel1 until pitang1 is reached )
      TrqRate        = 0.d0
      PitRate(:)     =-cutin_pitvel1/R2D
   elseif (BlPitch(1)<=cutin_pitang1/R2D) then  ! cut-in phase II ( Generator torque linear increase
     if     (cutin_step==2) then                !                   pitch control for low speed      )
      cutin_step     = 3
      cutin_timer    = time                     ! next timestep switch to normal operation
     endif
      TrqRate        = 0.d0
      PitRate(:)     = 0.d0
   endif

   ! Saturate Torque and pitch
      TrqRate        = min( max( TrqRate   ,-VS_MaxRat ), VS_MaxRat )
      GenTrq         = LastGenTrq + TrqRate*DT
      GenTrq         = min( max( GenTrq    , 0.d0      ), VS_MaxTq  )
   do K              = 1, NumBl
      PitRate(K)     = min( max( PitRate(K),-PC_MaxRat ), PC_MaxRat )
!!!   PitCom (K)     = LastBlPitch(K) + PitRate(K)*DT
      PitCom (K)     = LastPitComAve  + PitRate(K)*DT
      PitCom (K)     = min( max( PitCom (K), PC_MinPit ), PC_MaxPit )
      PitAcc (K)     = (PitRate(K)-LastPitRate(K))/DT
   enddo !K

   ! Store results for output
      res_vs(:)      = 0.d0;
      res_vs(1)      = Time
      res_vs(2)      = GenTrq
      res_vs(3)      = TrqRate
      res_vs(4)      = LastGenTrq
      res_vs(6)      = GenSpeed

      res_pc(:    )  = 0.d0;
      res_pc( 1   )  = Time
      res_pc( 2: 4)  = PitCom (1:3)
      res_pc( 5: 7)  = PitRate(1:3)
      res_pc( 8:10)  = PitAcc (1:3)
      res_pc(17:19)  = BlPitch(1:3)
      res_pc(20:22)  = LastBlPitch(1:3)
      res_pc(23   )  = GenSpeed
      res_pc(24   )  = Velhub
      res_pc(25   )  = PC_MinPit


 END Subroutine Cut_in
!---------------------------------
 Subroutine Cut_out
!---------------------------------

 implicit none

   if     (cutin_step  > 0               ) return

   if     (time < cutout_time1  .or. &            ! cut-out haven't started yet
           cutout_time1 <= 0.d0          ) then   ! cut-out will not simulated in this case
      cutout_step    = 0
      return
   elseif (time<cutout_time1+cutout_time2) then   ! cut-out phase I  ( Generator torque is constant     , pitch increase with vel1 )
      cutout_step    = 1
      TrqRate        = 0.d0
      PitRate(:)     = cutout_pitvel1 /R2D
   else                                           ! cut-out phase II ( Generator torque reduces linearly, pitch increase with vel2 )
     if     (cutout_step==1) then
      cutout_step    = 2
      cutout_TrqRate = (0.d0-LastGenTrq)/cutout_timegen
     elseif (BlPitch(1) >= PC_MaxPit) then
      cutout_pitvel2 = 0.d0
     endif
      TrqRate        = cutout_TrqRate
      PitRate(:)     = cutout_pitvel2 /R2D
   endif

   ! Saturate Torque and pitch
      TrqRate        = min( max( TrqRate   ,-VS_MaxRat ), VS_MaxRat )
      GenTrq         = LastGenTrq + TrqRate*DT
      GenTrq         = min( max( GenTrq    , 0.d0      ), VS_MaxTq  )
   do K              = 1, NumBl
      PitRate(K)     = min( max( PitRate(K),-PC_MaxRat ), PC_MaxRat )
      PitCom (K)     = LastBlPitch(K) + PitRate(K)*DT
      PitCom (K)     = min( max( PitCom (K), PC_MinPit ), PC_MaxPit )
      PitAcc (K)     = (PitRate(K)-LastPitRate(K))/DT
   enddo !K


   ! Store results for output
      res_vs(:)      = 0.d0;
      res_vs(1)      = Time
      res_vs(2)      = GenTrq
      res_vs(3)      = TrqRate
      res_vs(4)      = LastGenTrq
      res_vs(6)      = GenSpeed

      res_pc(:    )  = 0.d0;
      res_pc( 1   )  = Time
      res_pc( 2: 4)  = PitCom (1:3)
      res_pc( 5: 7)  = PitRate(1:3)
      res_pc( 8:10)  = PitAcc (1:3)
      res_pc(17:19)  = BlPitch(1:3)
      res_pc(20:22)  = LastBlPitch(1:3)
      res_pc(23   )  = GenSpeed
      res_pc(24   )  = Velhub
      res_pc(25   )  = PC_MinPit


 END Subroutine Cut_out
!---------------------------------
 recursive Subroutine init_back_out (istep)
!---------------------------------

 implicit none

   integer, intent(in) :: istep


 select case (istep)
 case (0)  ! Initialization


   ! Determine some torque control parameters not specified directly:

      iconst_pow = 0 !MAKE REGION 3 CONSTANT TORQUE INSTEAD OF CONSTANT POWER FOR FLOATING WT
   if (icase >= 3) then
      iconst_pow = 0 !MAKE REGION 3 CONSTANT TORQUE INSTEAD OF CONSTANT POWER FOR FLOATING WT
      PC_KI      = 0.0008965149d0 !REDUCE GAINS FOR HYWIND
      PC_KP      = 0.0062756040d0 !REDUCE GAINS FOR HYWIND
   else           !v8      !v7      !v6      !v5      !v4      !v3       !v2      !v1
      PC_KI      = 0.005d0 !0.005d0 !0.004d0 !0.003d0 !0.003d0 !0.0025d0 !0.001d0 !0.0040d0 !0.008068634d0
      PC_KP      = 0.015d0 !0.013d0 !0.016d0 !0.009d0 !0.006d0 !0.0075d0 !0.003d0 !0.0120d0 !0.018826810d0 
!qq
     if (IPID2==1) then
      PC_KI      = 0.005d0 *6 !0.005d0 !0.004d0 !0.003d0 !0.003d0 !0.0025d0 !0.001d0 !0.0040d0 !0.008068634d0
      PC_KP      = 0.015d0 *6 !0.013d0 !0.016d0 !0.009d0 !0.006d0 !0.0075d0 !0.003d0 !0.0120d0 !0.018826810d0 
     endif
   endif
                   !torque at region 1.5/2.0     omega 1.5/2.0 - 1.0/1.5
!     VS_Slope15 = ( VS_Rgn2K * VS_Rgn2Sp**VS_pow ) / ( VS_Rgn2Sp - VS_CtInSp )
      VS_Slope15 = ( Ftorq( VS_Rgn2Sp )           ) / ( VS_Rgn2Sp - VS_CtInSp )
    if     (iconst_pow == 0) then
      VS_Slope25 = ( VS_RtPwr / PC_RefSpd         ) / ( VS_RtGnSp - VS_SySp   )
    elseif (iconst_pow == 1) then
      VS_Slope25 = ( VS_RtPwr / VS_RtGnSp         ) / ( VS_RtGnSp - VS_SySp   )
    endif


   if ( VS_Rgn2K == 0.d0 )  then  ! .TRUE. if the Region 2 torque is flat, and thus, the denominator in the else condition is zero
      VS_TrGnSp  = VS_SySp
   else                           ! .TRUE. if the Region 2 torque is quadratic with speed
!------ based on trionymo T2=T2.5=> Q w^2=B(w-ws)
!     VS_TrGnSp  = ( VS_Slope25 - dsqrt( VS_Slope25*( VS_Slope25 - 4.d0*VS_Rgn2K*VS_SySp ) ) )/( 2.d0*VS_Rgn2K )

!------ find VS_TrGnSp based on Newton iteerations
         VS_TrGnSp = (MIN_Spd+MAX_Spd)/2.d0
      do ii        = 1, 15
!         func     =        VS_TrGnSp** VS_pow      -VS_Slope25/VS_Rgn2K * (VS_TrGnSp-VS_SySp)
!        dfunc     = VS_pow*VS_TrGnSp**(VS_pow-1.d0)-VS_Slope25/VS_Rgn2K
          func     =  Ftorq( VS_TrGnSp) - VS_Slope25 * (VS_TrGnSp-VS_SySp)
         dfunc     = dFtorq( VS_TrGnSp) - VS_Slope25
         write(*,'(i3,5e20.10)') ii, VS_TrGnSp, VS_TrGnSp-func/dfunc, -func/dfunc
         VS_TrGnSp = VS_TrGnSp -func/dfunc
         if (dabs(-func/dfunc)<1.d-10) exit
      enddo
!!!      VS_TrGnSp = 12.95d0/RPS2RPM *113.2d0
   endif

  open  (400,file='gen_control_init.dat')
                !  1          2         3          4         5          6          7           8
   write(400,10) VS_CtInSp, VS_Rgn2Sp, VS_TrGnSp, VS_SySp , VS_RtGnSp, PC_RefSpd, VS_Slope15, VS_Slope25, &
                 VS_RtPwr , VS_MaxTq , VS_MaxRat, VS_Rgn2K, PC_KI    , PC_KP    , VS_Rgn3MP , VS_SlPc 
                !  9         10         11         12        13         14         15          16
  close (400)


   ! Initialize the SAVE variables:
   ! NOTE: LastGenTrq, though SAVE, is initialized in the torque controller, not here.
   LastPitRate(:) = 0.d0
   LastPitComAve  = 0.d0
  do I            = 1, NumBl
   LastPitComAve  = LastPitComAve + BlPitch(I)/NumBl    ! This will ensure that the variable speed controller picks the correct control region and the pitch controller pickes the correct gain on the first call
  enddo
   LastGenSpeedF  = GenSpeed                            ! This will ensure that generator speed filter will use the initial value of the generator speed on the first pass
   PitCom     (:) = BlPitch(:);                         ! This will ensure that the variable speed controller picks the correct control region and the pitch controller pickes the correct gain on the first call
   LastBlPitch(:) = BlPitch(:);        
   GK             = 1.d0/( 1.d0 + LastPitComAve/PC_KK ) ! This will ensure that the pitch angle is unchanged if the initial SpdErr is zero
   LastIntSpdErr  = LastPitComAve /( GK*PC_KI )         ! This will ensure that the pitch angle is unchanged if the initial SpdErr is zero
   IntSpdErr      = 0.d0
!qq
   V(1)=-50.d0 ;  P(1)= 0.00/R2D;
   V(2)=  8.d0 ;  P(2)= 0.00/R2D;
   V(3)= 11.d0 ;  P(3)= 3.40/R2D; !P(3)= 9.40/R2D;
   V(4)= 15.d0 ;  P(4)=12.00/R2D; !P(3)= 9.40/R2D;
   V(5)= 55.d0 ;  P(5)=12.00/R2D; !P(4)= 9.40/R2D;
                  P(:)= 0.00/R2D; !P(4)= 9.40/R2D;

 case (2) !update for next timestep

   LastGenSpeedF  = GenSpeedF
   LastGenTrq     = GenTrq

   LastPitComAve  = 0.d0
  do I            = 1, NumBl
   LastPitComAve  = LastPitComAve + PitCom(I)/NumBl
  enddo

   LastIntSpdErr        = IntSpdErr
   LastBlPitch(1:NumBl) = BlPitch(1:NumBl)
   LastPitRate(1:NumBl) = PitRate(1:NumBl)

   errorpp(:) = errorp(:)
   errorp (:) = error (:)

!----- Output main control vars
!        1      2       3     4          5          6       7     8      9
!      Time,GenTrq,TrqRate,LastGenTrq,GenSpeedF,GenSpeed,Region,Alpha,GenTrqD
     open (400,file='gen_control_VS.dat',access='append')
      write(400,10) (res_vs(I),I=1,9), dble(cutin_step), dble(cutout_step)
     close(400)
!         1      2-4       5-7          8-10       11   12    
!      Time,PitCom(1:3),PitRate(1:3),PitAcc(1:3),GK,SpdErr,&
!      IntSpdErr,PitComP,PitComI,PitComT,BlPitch(1:3),LastBlPitch(1:3),GenSpeed, Velhub, PC_MinPit, PC_RefSpd
!        13       14       15      16     17-19         20-22           23         24       25         26
     open (401,file='gen_control_PC.dat',access='append')
      write(401,10) (res_pc(I),I=1,26)
     close(401)

 case (3) !create backup

   write (2,100) LastGenSpeedF
   write (2,100) LastGenTrq
   write (2,100) LastBlPitch (1:NumBl)
   write (2,100) LastPitComAve
   write (2,100) LastIntSpdErr
   write (2,100) LastPitRate (1:NumBl)

 case (4) !read backup

   call init_back_out (0) !1st perform the initialization

   read (2,100) LastGenSpeedF
   read (2,100) LastGenTrq
   read (2,100) LastBlPitch (1:NumBl)
   read (2,100) LastPitComAve
   read (2,100) LastIntSpdErr
   read (2,100) LastPitRate (1:NumBl)

 end select !istep

  10 format (20000f15.4 )
 100 format (20000e28.17)


 END Subroutine init_back_out

END Subroutine CONTROLOS3
!------------------------------------------------------------------------
 Subroutine PID2 ( N, dt, error, errorp, errorpp, Kp, Ki, Kd, GAIN, &
                   output_min, output_max, doutput_max, output )
!------------------------------------------------------------------------

   implicit none

   integer, intent(in   ) :: N
   real(8), intent(in   ) :: dt
   real(8), intent(in   ) :: error (N)
   real(8), intent(inout) :: errorp(N), errorpp(N)
   real(8), intent(in   ) :: Kp    (N), Ki(N), Kd(N), GAIN(N)
   real(8), intent(in   ) :: output_min, output_max, doutput_max
   real(8), intent(inout) :: output

   integer, save :: initialize =  0
   real(8)       :: derror(N) , dderror(N), doutput
   integer       :: i

 
   if (initialize==0) then
       initialize = 1
       errorp (:) = error(:)
       errorpp(:) = error(:)
   endif

!  error      = measured_value - set_value
   derror(:)  = ( 3.d0*error(:) - 4.d0 * errorp(:) + errorpp(:)) / (2.d0*dt)
  dderror(:)  = (      error(:) - 2.d0 * errorp(:) + errorpp(:)) /       dt**2
   doutput    = 0.d0
  do i        = 1, N
   doutput    = doutput + GAIN(i) * ( Ki(i)*error(i) + Kp(i)*derror(i) + Kd(i)*dderror(i) ) * dt
  enddo

!--- Saturate output
  doutput     = min( max(doutput,-doutput_max),doutput_max)
   output     = output + doutput
   output     = min( max( output,  output_min), output_max)

!--- update, done elsewhere in order to account for iterations
!! errorpp(:) = errorp(:)
!! errorp (:) = error (:)


 END Subroutine PID2
