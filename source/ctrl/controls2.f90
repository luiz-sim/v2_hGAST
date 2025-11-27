!=======================================================================
Subroutine CONTROLOS2 ( Time      , GenSpeed  , BlPitch               , &
                        GenTrq_out, PitAng_out, PitVel_out, PitAcc_out, &
                        NumBl     , iStatus                               )


   ! This controller is used to implement a variable-speed generator-torque
   ! controller and PI collective blade pitch controller for the NREL 5 MW
   ! Offshore baseline wind turbine.

IMPLICIT                        NONE

   ! Local Variables:

REAL(kind=8)                 :: Alpha                                           ! Current coefficient in the recursive, single-pole, low-pass filter, (-).
REAL(kind=8)                 :: BlPitch   (3)                                   ! Current values of the blade pitch angles, rad.
REAL(kind=8)                 :: ElapTime                                        ! Elapsed time since the last call to the controller, sec.
REAL(kind=8), PARAMETER      :: CornerFreq    =       1.570796d0                ! Corner frequency (-3dB point) in the recursive, single-pole, low-pass filter, rad/s. -- chosen to be 1/4 the blade edgewise natural frequency ( 1/4 of approx. 1Hz = 0.25Hz = 1.570796rad/s)
REAL(kind=8)                 :: GenSpeed                                        ! Current  HSS (generator) speed, rad/s.
REAL(kind=8), SAVE           :: GenSpeedF                                       ! Filtered HSS (generator) speed, rad/s.
REAL(kind=8)                 :: GenTrq                                          ! Electrical generator torque, N-m.
REAL(kind=8)                 :: GK                                              ! Current value of the gain correction factor, used in the gain scheduling law of the pitch controller, (-).

REAL(kind=8), SAVE           :: IntSpdErr                                       ! Current integral of speed error w.r.t. time, rad.
REAL(kind=8), SAVE           :: LastGenTrq                                      ! Commanded electrical generator torque the last time the controller was called, N-m.
REAL(kind=8), SAVE           :: LastTime                                        ! Last time this DLL was called, sec.
REAL(kind=8), SAVE           :: LastTimePC                                      ! Last time the pitch  controller was called, sec.
REAL(kind=8), SAVE           :: LastTimeVS                                      ! Last time the torque controller was called, sec.
!REAL(kind=8), PARAMETER      :: OnePlusEps    = 1.0 + EPSILON(OnePlusEps)      ! The number slighty greater than unity in single precision.
REAL(kind=8), PARAMETER      :: PC_DT         =       0.1d0     !! 0.0001       ! Communication interval for pitch  controller, sec.
REAL(kind=8), PARAMETER      :: PC_KI         =       0.006d0   !! 0.008068634d0! Integral gain for pitch controller at rated pitch (zero), (-).

REAL(kind=8), PARAMETER      :: PC_KK         =       0.1099965d0               ! Pitch angle were the derivative of the aerodynamic power w.r.t. pitch has increased by a factor of two relative to the derivative at rated pitch (zero), rad.
REAL(kind=8), PARAMETER      :: PC_KP         =       0.01d0    !! 0.01882681d0 ! Proportional gain for pitch controller at rated pitch (zero), sec.

REAL(kind=8), PARAMETER      :: PC_MaxPit     =       1.570796d0                ! Maximum pitch setting in pitch controller, rad.
REAL(kind=8), PARAMETER      :: PC_MaxRat     =       0.1396263d0               ! Maximum pitch  rate (in absolute value) in pitch  controller, rad/s.
REAL(kind=8), PARAMETER      :: PC_MinPit     =       0.d0                      ! Minimum pitch setting in pitch controller, rad.
REAL(kind=8), PARAMETER      :: PC_RefSpd     =     122.9096d0     ! 1173.7rpm  ! Desired (reference) HSS speed for pitch controller, rad/s.
REAL(kind=8), SAVE           :: PitCom    (3)                                   ! Commanded pitch of each blade the last time the controller was called, rad.
REAL(kind=8)                 :: PitComI                                         ! Integral term of command pitch, rad.
REAL(kind=8)                 :: PitComP                                         ! Proportional term of command pitch, rad.
REAL(kind=8)                 :: PitComT                                         ! Total command pitch based on the sum of the proportional and integral terms, rad.
REAL(kind=8)                 :: PitRate   (3)                                   ! Pitch rates of each blade based on the current pitch angles and current pitch command, rad/s.
REAL(kind=8), SAVE           :: R2D                 !57.29578d0                 ! Factor to convert radians to degrees.
REAL(kind=8), SAVE           :: RPS2RPM             ! 9.5492966d0               ! Factor to convert radians per second to revolutions per minute.
REAL(kind=8)                 :: SpdErr                                          ! Current speed error, rad/s.
REAL(kind=8)                 :: Time                                            ! Current simulation time, sec.
REAL(kind=8)                 :: TrqRate                                         ! Torque rate based on the current and last torque commands, N-m/s.
REAL(kind=8), PARAMETER      :: VS_CtInSp     =      70.16224d0    ! 670rpm     ! Transitional generator speed (HSS side) between regions 1 and 1 1/2, rad/s.
REAL(kind=8), PARAMETER      :: VS_DT         =       0.01d0    !! 0.0001d0     ! Communication interval for torque controller, sec.
REAL(kind=8), PARAMETER      :: VS_MaxRat     =   15000.d0                      ! Maximum torque rate (in absolute value) in torque controller, N-m/s.
REAL(kind=8), PARAMETER      :: VS_MaxTq      =   47402.91d0       !43093.55*1.1! Maximum generator torque in Region 3 (HSS side), N-m. -- chosen to be 10% above VS_RtTq = 43.09355kNm
REAL(kind=8), PARAMETER      :: VS_Rgn2K      =       2.332287d0                ! Generator torque constant in Region 2 (HSS side), N-m/(rad/s)^2.
REAL(kind=8), PARAMETER      :: VS_Rgn2Sp     =      91.21091d0    ! 871rpm     ! Transitional generator speed (HSS side) between regions 1 1/2 and 2, rad/s.
REAL(kind=8), PARAMETER      :: VS_Rgn3MP     =       0.01745329d0 ! 1 deg.     ! Minimum pitch angle at which the torque is computed as if we are in region 3 regardless of the generator speed, rad. -- chosen to be 1.0 degree above PC_MinPit
REAL(kind=8), PARAMETER      :: VS_RtGnSp     =     121.6805d0     ! 1162rpm    ! Rated generator speed (HSS side), rad/s. -- chosen to be 99% of PC_RefSpd
REAL(kind=8), PARAMETER      :: VS_RtPwr      = 5296610.d0         ! 5MW/0.944  ! Rated generator generator power in Region 3, Watts. -- chosen to be 5MW divided by the electrical generator efficiency of 94.4%
REAL(kind=8), SAVE           :: VS_Slope15                                      ! Torque/speed slope of region 1 1/2 cut-in torque ramp , N-m/(rad/s).
REAL(kind=8), SAVE           :: VS_Slope25                                      ! Torque/speed slope of region 2 1/2 induction generator, N-m/(rad/s).
REAL(kind=8), PARAMETER      :: VS_SlPc       =      10.d0                      ! Rated generator slip percentage in Region 2 1/2, %.
REAL(kind=8), SAVE           :: VS_SySp                                         ! Synchronous speed of region 2 1/2 induction generator, rad/s.
REAL(kind=8), SAVE           :: VS_TrGnSp                                       ! Transitional generator speed (HSS side) between regions 2 and 2 1/2, rad/s.

INTEGER                      :: I                                               ! Generic index.
INTEGER                      :: iStatus                                         ! A status flag set by the simulation as follows: 0 if this is the first call, 1 for all subsequent time steps, -1 if this is the final call at the end of the simulation.
INTEGER                      :: K                                               ! Loops through blades.
INTEGER                      :: NumBl                                           ! Number of blades, (-).


REAL(kind=8)                 :: PitAng_out(3)                                   ! Communicate Pitch angle
REAL(kind=8)                 :: PitVel_out(3)                                   ! Communicate Pitch vel
REAL(kind=8)                 :: PitAcc_out(3)                                   ! Communicate Pitch accel
REAL(kind=8)                 :: GenTrq_out                                      ! Communicate Generator's Torque Demand
REAL(kind=8), SAVE           :: PitVel0_out(3)                                  ! Keep Old PITCH VELOCITY

   ! Load variables from calling program

!iStatus      = 
!NumBl        = 

!BlPitch  (1) = 
!BlPitch  (2) = 
!BlPitch  (3) = 
!GenSpeed     = 
!Time         = 
   


   ! Read any External Controller Parameters specified in the User Interface
   !   and initialize variables:

if ( iStatus == 0 )  then  ! .TRUE. if were on the first call to the DLL

   PitVel0_out(1:NumBl) = 0.d0   ! Initialize
   R2D      = 180.d0 / dacos(-1.d0)
   RPS2RPM  =  30.d0 / dacos(-1.d0) 

   ! Determine some torque control parameters not specified directly:

   VS_SySp    = VS_RtGnSp/( 1.d0 +  0.01d0*VS_SlPc )
   VS_Slope15 = ( VS_Rgn2K * VS_Rgn2Sp**2 ) / ( VS_Rgn2Sp - VS_CtInSp )
   VS_Slope25 = ( VS_RtPwr / VS_RtGnSp    ) / ( VS_RtGnSp - VS_SySp   )
!F VS_Slope25 = ( VS_RtPwr / PC_RefSpd    ) / ( VS_RtGnSp - VS_SySp   ) !MAKE REGION 3 CONSTANT TORQUE INSTEAD OF CONSTANT POWER FOR HYWIND

   if ( VS_Rgn2K == 0.d0 )  then  ! .TRUE. if the Region 2 torque is flat, and thus, the denominator in the else condition is zero
      VS_TrGnSp = VS_SySp
   else                           ! .TRUE. if the Region 2 torque is quadratic with speed
      VS_TrGnSp = ( VS_Slope25 - dsqrt( VS_Slope25*( VS_Slope25 - 4.d0*VS_Rgn2K*VS_SySp ) ) )/( 2.d0*VS_Rgn2K )
   endif

  !write(*,*) 'VS_TrGnSp',VS_TrGnSp,VS_SySp,VS_TrGnSp*30d0/3.14159d0

   ! Check validity of input parameters:

   if ( CornerFreq <= 0.d0 )  then
      write(*,*)'CornerFreq must be greater than zero.'
      stop
   endif

   if ( VS_DT     <= 0.d0 )  then
      write(*,*)'VS_DT must be greater than zero.'
      stop
   endif

   if ( VS_CtInSp <  0.d0 )  then
      write(*,*)'VS_CtInSp must not be negative.'
      stop
   endif

   if ( VS_Rgn2Sp <= VS_CtInSp )  then
      write(*,*)'VS_Rgn2Sp must be greater than VS_CtInSp.'
      stop
   endif

   if ( VS_TrGnSp <  VS_Rgn2Sp )  then
      write(*,*)'VS_TrGnSp must not be less than VS_Rgn2Sp.'
      stop
   endif

   if ( VS_SlPc   <= 0.d0 )  then
      write(*,*)'VS_SlPc must be greater than zero.'
      stop
   endif

   if ( VS_MaxRat <= 0.d0 )  then
      write(*,*)'VS_MaxRat must be greater than zero.'
      stop
   endif

   if ( VS_RtPwr  <  0.d0 )  then
      write(*,*)'VS_RtPwr must not be negative.'
      stop
   endif

   if ( VS_Rgn2K  <  0.d0 )  then
      write(*,*)'VS_Rgn2K must not be negative.'
      stop
   endif

   if ( VS_Rgn2K*VS_RtGnSp*VS_RtGnSp > VS_RtPwr/VS_RtGnSp )  then
      write(*,*)'VS_Rgn2K*VS_RtGnSp^2 must not be greater than VS_RtPwr/VS_RtGnSp.'
      stop
   endif

   if ( VS_MaxTq                     < VS_RtPwr/VS_RtGnSp )  then
      write(*,*)'VS_RtPwr/VS_RtGnSp must not be greater than VS_MaxTq.'
      stop
   endif

   if ( PC_DT     <= 0.d0 )  then
      write(*,*)'PC_DT must be greater than zero.'
      stop
   endif

   if ( PC_KI     <= 0.d0 )  then
      write(*,*)'PC_KI must be greater than zero.'
      stop
   endif

   if ( PC_KK     <= 0.d0 )  then
      write(*,*)'PC_KK must be greater than zero.'
      stop
   endif

   if ( PC_RefSpd <= 0.d0 )  then
      write(*,*)'PC_RefSpd must be greater than zero.'
      stop
   endif
   
   if ( PC_MaxRat <= 0.d0 )  then
      write(*,*)'PC_MaxRat must be greater than zero.'
      stop
   endif

   if ( PC_MinPit >= PC_MaxPit )  then
      write(*,*)'PC_MinPit must be less than PC_MaxPit.'
      stop
   endif


   ! Initialize the SAVEd variables:
   ! NOTE: LastGenTrq, though SAVEd, is initialized in the torque controller
   !       below for simplicity, not here.

   GenSpeedF  = GenSpeed                        ! This will ensure that generator speed filter will use the initial value of the generator speed on the first pass
   PitCom     = BlPitch                         ! This will ensure that the variable speed controller picks the correct control region and the pitch controller pickes the correct gain on the first call
   GK         = 1.d0/( 1.d0 + PitCom(1)/PC_KK ) ! This will ensure that the pitch angle is unchanged if the initial SpdErr is zero
   IntSpdErr  = PitCom(1)/( GK*PC_KI )          ! This will ensure that the pitch angle is unchanged if the initial SpdErr is zero

   LastTime   = Time                            ! This will ensure that generator speed filter will use the initial value of the generator speed on the first pass
   LastTimePC = Time - PC_DT - 2.d-8            ! This will ensure that the pitch  controller is called on the first pass 
   LastTimeVS = Time - VS_DT - 2.d-8            ! This will ensure that the torque controller is called on the first pass 


endif !istatus==0



   ! Main control calculations:

if  ( iStatus >= 0 )  then                      ! Only compute control calculations if no error has occured and we are not on the last time step


!=======================================================================


   ! Filter the HSS (generator) speed measurement:
   ! NOTE: This is a very simple recursive, single-pole, low-pass filter with
   !       exponential smoothing.

   ! Update the coefficient in the recursive formula based on the elapsed time
   !   since the last call to the controller:

   Alpha     = dexp( ( LastTime - Time )*CornerFreq )


   ! Apply the filter:

   GenSpeedF = ( 1.d0 - Alpha )*GenSpeed + Alpha*GenSpeedF


!=======================================================================


   ! Variable-speed torque control:

   ! Compute the elapsed time since the last call to the controller:

   ElapTime = Time - LastTimeVS


   ! Only perform the control calculations if the elapsed time is greater than
   !   or equal to the communication interval of the torque controller:
   ! NOTE: Time is scaled by OnePlusEps to ensure that the contoller is called
   !       at every time step when VS_DT = DT, even in the presence of
   !       numerical precision errors.

!! if ( ( Time*OnePlusEps - LastTimeVS ) >= VS_DT )  then
   if (  ElapTime - VS_DT  > 1.d-08 )  then


   ! Compute the generator torque, which depends on which region we are in:

      if ( (   GenSpeedF >= VS_RtGnSp ) .OR. (  PitCom(1) >= VS_Rgn3MP ) )  then ! We are in region 3 - power is constant
         GenTrq = VS_RtPwr/GenSpeedF
!F       GenTrq = VS_RtPwr/PC_RefSpd   !MAKE REGION 3 CONSTANT TORQUE INSTEAD OF CONSTANT POWER FOR HYWIND
         write(400,*)'#3',GenTrq
      elseif ( GenSpeedF <= VS_CtInSp )  then                                    ! We are in region 1 - torque is zero
         GenTrq = 0.d0
         write(400,*)'#1',GenTrq
      elseif ( GenSpeedF <  VS_Rgn2Sp )  then                                    ! We are in region 1 1/2 - linear ramp in torque from zero to optimal
         GenTrq = VS_Slope15*( GenSpeedF - VS_CtInSp )
         write(400,*)'#1.5',GenTrq
      elseif ( GenSpeedF <  VS_TrGnSp )  then                                    ! We are in region 2 - optimal torque is proportional to the square of the generator speed
         GenTrq = VS_Rgn2K*GenSpeedF**2        
         write(400,*)'#2',GenTrq
      else                                                                       ! We are in region 2 1/2 - simple induction generator transition region
         GenTrq = VS_Slope25*( GenSpeedF - VS_SySp   )
         write(400,*)'#2.5',GenTrq
      endif


   ! Saturate the commanded torque using the maximum torque limit:

      GenTrq  = dmin1( GenTrq                      , VS_MaxTq  )  ! Saturate the command using the maximum torque limit


   ! Saturate the commanded torque using the torque rate limit:

      if ( iStatus == 0 )  LastGenTrq = GenTrq                    ! Initialize the value of LastGenTrq on the first pass only
      TrqRate = ( GenTrq - LastGenTrq ) / ElapTime                ! Torque rate (unsaturated)
      TrqRate = dmin1( dmax1( TrqRate, -VS_MaxRat ), VS_MaxRat )  ! Saturate the torque rate using its maximum absolute value
      GenTrq  = LastGenTrq + TrqRate*ElapTime                     ! Saturate the command using the torque rate limit


   ! Reset the values of LastTimeVS and LastGenTrq to the current values:

      LastTimeVS = Time
      LastGenTrq = GenTrq                                      ! Demanded generator torque

      GenTrq_out = GenTrq

      write(400,10)Time,GenTrq,TrqRate,LastGenTrq,GenSpeedF,GenSpeed,Alpha

   endif !Variable Speed (Generator Torque)


!=======================================================================


   ! Pitch control:

   ! Compute the elapsed time since the last call to the controller:

   ElapTime = Time - LastTimePC


   ! Only perform the control calculations if the elapsed time is greater than
   !   or equal to the communication interval of the pitch controller:
   ! NOTE: Time is scaled by OnePlusEps to ensure that the contoller is called
   !       at every time step when PC_DT = DT, even in the presence of
   !       numerical precision errors.

!! if ( ( Time*OnePlusEps - LastTimePC ) >= PC_DT )  then
   if ( ElapTime - PC_DT > 1.d-08 )  then


   ! Compute the gain scheduling correction factor based on the previously
   !   commanded pitch angle for blade 1:

      GK = 1.d0/( 1.d0 + PitCom(1)/PC_KK )


   ! Compute the current speed error and its integral w.r.t. time; saturate the
   !   integral term using the pitch angle limits:

      SpdErr    = GenSpeedF - PC_RefSpd                                 ! Current speed error
      IntSpdErr = IntSpdErr + SpdErr*ElapTime                           ! Current integral of speed error w.r.t. time
      IntSpdErr = MIN( MAX( IntSpdErr, PC_MinPit/( GK*PC_KI ) ), &
                                       PC_MaxPit/( GK*PC_KI )      )    ! Saturate the integral term using the pitch angle limits, converted to integral speed error limits


   ! Compute the pitch commands associated with the proportional and integral
   !   gains:

      PitComP   = GK*PC_KP*   SpdErr                                    ! Proportional term
      PitComI   = GK*PC_KI*IntSpdErr                                    ! Integral term (saturated)


   ! Superimpose the individual commands to get the total pitch command;
   !   saturate the overall command using the pitch angle limits:

      PitComT   = PitComP + PitComI                                     ! Overall command (unsaturated)
      PitComT   = MIN( MAX( PitComT, PC_MinPit ), PC_MaxPit )           ! Saturate the overall command using the pitch angle limits


   ! Saturate the overall commanded pitch using the pitch rate limit:
   ! NOTE: Since the current pitch angle may be different for each blade
   !       (depending on the type of actuator implemented in the structural
   !       dynamics model), this pitch rate limit calculation and the
   !       resulting overall pitch angle command may be different for each
   !       blade.

      do K = 1,NumBl ! Loop through all blades

         PitRate(K) = ( PitComT - BlPitch(K) )/ElapTime                 ! Pitch rate of blade K (unsaturated)
         PitRate(K) = MIN( MAX( PitRate(K), -PC_MaxRat ), PC_MaxRat )   ! Saturate the pitch rate of blade K using its maximum absolute value
         PitCom (K) = BlPitch(K) + PitRate(K)*ElapTime                  ! Saturate the overall command of blade K using the pitch rate limit

      enddo          ! K - all blades


   ! Reset the value of LastTimePC to the current value:

      LastTimePC = Time

      PitAng_out(1:NumBl) =  PitCom     (1:NumBl)
      PitVel_out(1:NumBl) =  PitRate    (1:NumBl)
      PitAcc_out(1:NumBl) = (PitVel_out (1:NumBl)-PitVel0_out (1:NumBl))/ElapTime

      PitVel0_out(1:NumBl) = PitVel_out(1:NumBl) 

   endif  !Pitch control


!=======================================================================


   ! Reset the value of LastTime to the current value:

   LastTime = Time


endif !istatus>=0

10 format (20f15.4)


END Subroutine CONTROLOS2
