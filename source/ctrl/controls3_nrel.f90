!=======================================================================
!Subroutine CONTROLOS3 ( Time      , GenSpeed  , BlPitch, &
!                        GenTrq_out, PitAng_out, PitVel_out, PitAcc_out, &
!                        NumBl     , iStatus   , icase      )
Subroutine CONTROLOS3 ( Time      , GenSpeed  , BlPitch   , Velhub    , &
                        NumBl     , iStatus   , icase     , DT        , &
                        GenTrq_out, PitAng_out, PitVel_out, PitAcc_out    )


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
REAL(8), intent (out)   :: PitAng_out(3)                                   ! Communicate Pitch angle
REAL(8), intent (out)   :: PitVel_out(3)                                   ! Communicate Pitch vel
REAL(8), intent (out)   :: PitAcc_out(3)                                   ! Communicate Pitch accel

   ! Local Variables:
REAL(8)                 :: Alpha                                           ! Current coefficient in the recursive, single-pole, low-pass filter, (-).
REAL(8)                 :: ElapTime                                        ! Elapsed time since the last call to the controller, sec.
REAL(8), PARAMETER      :: CornerFreq    =       1.570796d0   !   1Hz      !#Corner frequency (-3dB point) in the recursive, single-pole, low-pass filter, rad/s. -- chosen to be 1/4 the blade edgewise natural frequency ( 1/4 of approx. 1Hz = 0.25Hz = 1.570796rad/s)
REAL(8), SAVE           :: GenSpeedF                                       ! Filtered HSS (generator) speed, rad/s.
REAL(8),      SAVE      :: GenTrq        !jim add SAVE                     ! Electrical generator torque, N-m.
REAL(8)                 :: GK                                              ! Current value of the gain correction factor, used in the gain scheduling law of the pitch controller, (-).

REAL(8), SAVE           :: IntSpdErr                                       ! Current integral of speed error w.r.t. time, rad.
REAL(8), SAVE           :: LastGenTrq                                      ! Commanded electrical generator torque the last time the controller was called, N-m.
REAL(8), SAVE           :: LastTime                                        ! Last time this DLL was called, sec.
REAL(8), SAVE           :: LastTimePC                                      ! Last time the pitch  controller was called, sec.
REAL(8), SAVE           :: LastTimeVS                                      ! Last time the torque controller was called, sec.
!REAL=8), PARAMETER      :: OnePlusEps    = 1.0 + EPSILON(OnePlusEps)      ! The number slighty greater than unity in single precision.
REAL(8), PARAMETER      :: PC_DT         =       0.0001d0                  ! Communication interval for pitch  controller, sec.
REAL(8),      SAVE      :: PC_KI!        =       0.008068634d0             ! Integral gain for pitch controller at rated pitch (zero), (-).      0.0008965149d0 !REDUCE GAINS FOR HYWIND
REAL(8), PARAMETER      :: PC_KK         =       0.1099965d0               !#Pitch angle were the derivative of the aerodynamic power w.r.t. pitch has increased by a factor of two relative to the derivative at rated pitch (zero), rad.
REAL(8),      SAVE      :: PC_KP!        =       0.01882681d0              ! Proportional gain for pitch controller at rated pitch (zero), sec.  0.006275604d0  !REDUCE GAINS FOR HYWIND
REAL(8), PARAMETER      :: PC_MaxPit     =       1.570796d0   !   90deg    ! Maximum pitch setting in pitch controller, rad.
REAL(8), PARAMETER      :: PC_MaxRat     =       0.1396263d0  !    8deg/s  ! Maximum pitch  rate (in absolute value) in pitch  controller, rad/s.
REAL(8), PARAMETER      :: PC_MinPit     =       0.d0         !    0deg    ! Minimum pitch setting in pitch controller, rad.
REAL(8), PARAMETER      :: PC_RefSpd     =     122.9096d0     ! 1173.7rpm  !#Desired (reference) HSS speed for pitch controller, rad/s.
REAL(8), SAVE           :: PitCom    (3)                                   ! Commanded pitch of each blade the last time the controller was called, rad.
REAL(8)                 :: PitComI                                         ! Integral term of command pitch, rad.
REAL(8)                 :: PitComP                                         ! Proportional term of command pitch, rad.
REAL(8)                 :: PitComT                                         ! Total command pitch based on the sum of the proportional and integral terms, rad.
REAL(8),       SAVE     :: PitRate   (3)  !jim add SAVE                    ! Pitch rates of each blade based on the current pitch angles and current pitch command, rad/s.
REAL(8), SAVE           :: R2D                 !57.29578d0                 ! Factor to convert radians to degrees.
REAL(8), SAVE           :: RPS2RPM             ! 9.5492966d0               ! Factor to convert radians per second to revolutions per minute.
REAL(8)                 :: SpdErr                                          ! Current speed error, rad/s.
REAL(8)                 :: TrqRate                                         ! Torque rate based on the current and last torque commands, N-m/s.
REAL(8), PARAMETER      :: VS_DT         =       0.0001d0                  ! Communication interval for torque controller, sec.
REAL(8), PARAMETER      :: VS_RtPwr      = 5296610.d0         ! 5MW/0.944  !#Rated generator generator power in Region 3, Watts. -- chosen to be 5MW divided by the electrical generator efficiency of 94.4%
REAL(8), PARAMETER      :: VS_MaxTq      = VS_RtPwr/PC_RefSpd * 1.1d0      ! 43093.55*1.1 ! Maximum generator torque in Region 3 (HSS side), N-m. -- chosen to be 10% above VS_RtTq = 43.09355kNm
REAL(8), PARAMETER      :: VS_MaxRat     =  VS_MaxTq * 0.3d0  !15000       ! Maximum torque rate (in absolute value) in torque controller, N-m/s.
REAL(8), PARAMETER      :: VS_Rgn2K      =       2.332287d0                !#Generator torque constant in Region 2 (HSS side), N-m/(rad/s)^2.
REAL(8), PARAMETER      :: VS_CtInSp     =      70.16224d0    ! 670rpm     !#Transitional generator speed (HSS side) between regions 1   and 1 1/2, rad/s.
REAL(8), PARAMETER      :: VS_Rgn2Sp     = VS_CtInSp * 1.3d0  ! 871rpm     ! Transitional generator speed (HSS side) between regions 1/2 and 2    , rad/s.
REAL(8), SAVE           :: VS_TrGnSp                                       ! Transitional generator speed (HSS side) between regions 2   and 2 1/2, rad/s.
REAL(8), PARAMETER      :: VS_Rgn3MP     =       0.01745329d0 ! 1 deg.     ! Minimum pitch angle at which the torque is computed as if we are in region 3 regardless of the generator speed, rad. -- chosen to be 1.0 degree above PC_MinPit
REAL(8), PARAMETER      :: VS_RtGnSp     = PC_RefSpd * 0.99d0 ! 1162rpm    ! Rated generator speed (HSS side), rad/s. -- chosen to be 99% of PC_RefSpd
REAL(8), SAVE           :: VS_Slope15                                      ! Torque/speed slope of region 1 1/2 cut-in torque ramp , N-m/(rad/s).
REAL(8), SAVE           :: VS_Slope25                                      ! Torque/speed slope of region 2 1/2 induction generator, N-m/(rad/s).
REAL(8), PARAMETER      :: VS_SlPc       =      10.d0                      ! Rated generator slip percentage in Region 2 1/2, %. !jim: Defines VS_SySp = VS_RtGnSp/(1+VS_SlPc/100)
REAL(8), SAVE           :: VS_SySp                                         ! Synchronous speed of region 2 1/2 induction generator, rad/s.

INTEGER                 :: I                                               ! Generic index.
INTEGER                 :: K                                               ! Loops through blades.

REAL(8), SAVE           :: LastPitCom      !jim AVERAGED!!
REAL(8), SAVE           :: LastIntSpdErr   !jim
REAL(8),      SAVE      :: LastGenSpeedF   !jim
REAL(8), SAVE           :: LastPitRate(3)  !jim
REAL(8), SAVE           :: PitAcc     (3)  !jim
REAL(8),      SAVE      :: LastBlPitch(3)  !jim
REAL(8)                 :: Region          !jim output the control region


   ! Read any External Controller Parameters specified in the User Interface
   !   and initialize variables:

if ( iStatus == 0 )  then  ! Initialization

   R2D      = 180.d0 / dacos(-1.d0)
   RPS2RPM  =  30.d0 / dacos(-1.d0) 

   ! Determine some torque control parameters not specified directly:

   VS_SySp    = VS_RtGnSp/( 1.d0 + VS_SlPc/100.d0 )
   VS_Slope15 = ( VS_Rgn2K * VS_Rgn2Sp**2 ) / ( VS_Rgn2Sp - VS_CtInSp )

   if (icase >= 3) then
      VS_Slope25 = ( VS_RtPwr / PC_RefSpd    ) / ( VS_RtGnSp - VS_SySp   ) !MAKE REGION 3 CONSTANT TORQUE INSTEAD OF CONSTANT POWER FOR FLOATING WT
      PC_KI      = 0.0008965149d0 !REDUCE GAINS FOR HYWIND
      PC_KP      = 0.006275604d0  !REDUCE GAINS FOR HYWIND
!!    PC_KP      = 0.012367188d0        !gains for innwind
!!    PC_KI      = 0.012367188d0/45d0   !gains for innwind  0.000274826
   else
      VS_Slope25 = ( VS_RtPwr / VS_RtGnSp    ) / ( VS_RtGnSp - VS_SySp   )
      PC_KI      = 0.008068634d0
      PC_KP      = 0.01882681d0 
   endif

   if ( VS_Rgn2K == 0.d0 )  then  ! .TRUE. if the Region 2 torque is flat, and thus, the denominator in the else condition is zero
      VS_TrGnSp = VS_SySp
   else                           ! .TRUE. if the Region 2 torque is quadratic with speed
      VS_TrGnSp = ( VS_Slope25 - dsqrt( VS_Slope25*( VS_Slope25 - 4.d0*VS_Rgn2K*VS_SySp ) ) )/( 2.d0*VS_Rgn2K )
   endif

   open (400,file='oc_control_init.dat')
                !  1          2         3          4         5          6          7           8
   write(400,10) VS_CtInSp, VS_Rgn2Sp, VS_TrGnSp, VS_SySp , VS_RtGnSp, PC_RefSpd, VS_Slope15, VS_Slope25, &
                 VS_RtPwr , VS_MaxTq , VS_MaxRat, VS_Rgn2K, PC_KI    , PC_KP    , VS_Rgn3MP , VS_SlPc 
                !  9         10         11         12        13         14         15          16
   close(400)

   ! Initialize the SAVEd variables:
   ! NOTE: LastGenTrq, though SAVEd, is initialized in the torque controller
   !       below for simplicity, not here.


   LastPitRate (1:NumBl) = 0.d0   ! Initialize

   LastPitCom     = 0.d0
   do I = 1, NumBl
      LastPitCom  = LastPitCom + BlPitch(I)
   enddo
   LastPitCom     = LastPitCom / NumBl               ! This will ensure that the variable speed controller picks the correct control region and the pitch controller pickes the correct gain on the first call


   LastGenSpeedF  = GenSpeed                         ! This will ensure that generator speed filter will use the initial value of the generator speed on the first pass
   PitCom         = BlPitch; !(1:NumBl)              ! This will ensure that the variable speed controller picks the correct control region and the pitch controller pickes the correct gain on the first call
   LastBlPitch    = BlPitch; !(1:NumBl)
   GK             = 1.d0/( 1.d0 + LastPitCom/PC_KK ) ! This will ensure that the pitch angle is unchanged if the initial SpdErr is zero
   LastIntSpdErr  = LastPitCom /( GK*PC_KI )         ! This will ensure that the pitch angle is unchanged if the initial SpdErr is zero

   LastTime       = Time                             ! This will ensure that generator speed filter will use the initial value of the generator speed on the first pass
   LastTimePC     = Time - PC_DT - 2.d-8             ! This will ensure that the pitch  controller is called on the first pass 
   LastTimeVS     = Time - VS_DT - 2.d-8             ! This will ensure that the torque controller is called on the first pass 


elseif ( iStatus == 2 )  then  !update for next timestep

!-global + filter
   LastTime      = Time
   LastGenSpeedF = GenSpeedF

!-var speed
   ElapTime = Time - LastTimeVS
   if ( ElapTime - VS_DT > 1.d-08 ) then
       LastTimeVS = Time
       LastGenTrq = GenTrq
   endif

!-pitch
   ElapTime = Time - LastTimePC
   if ( ElapTime - PC_DT > 1.d-08 ) then

      LastPitCom           = 0.d0
      do I = 1, NumBl
         LastPitCom        = LastPitCom + PitCom(I)
      enddo
      LastPitCom           = LastPitCom / NumBl

      LastIntSpdErr        = IntSpdErr
      LastBlPitch(1:NumBl) = BlPitch(1:NumBl)
      LastPitRate(1:NumBl) = PitRate(1:NumBl)
      LastTimePC           = Time

   endif

  return

elseif ( iStatus == 3 )  then  !create backup

   write (2,100) LastTime
   write (2,100) LastGenSpeedF
   write (2,100) LastTimeVS
   write (2,100) LastGenTrq
   write (2,100) LastBlPitch (1:NumBl)
   write (2,100) LastPitCom
   write (2,100) LastIntSpdErr
   write (2,100) LastPitRate (1:NumBl)
   write (2,100) LastTimePC

  return

elseif ( iStatus == 4 )  then  !read backup

   read (2,100) LastTime
   read (2,100) LastGenSpeedF
   read (2,100) LastTimeVS
   read (2,100) LastGenTrq
   read (2,100) LastBlPitch (1:NumBl)
   read (2,100) LastPitCom
   read (2,100) LastIntSpdErr
   read (2,100) LastPitRate (1:NumBl)
   read (2,100) LastTimePC


   ! Determine some torque control parameters not specified directly:

   R2D        = 180.d0 / dacos(-1.d0)
   RPS2RPM    =  30.d0 / dacos(-1.d0) 

   VS_SySp    = VS_RtGnSp/( 1.d0 + VS_SlPc/100.d0 )
   VS_Slope15 = ( VS_Rgn2K * VS_Rgn2Sp**2 ) / ( VS_Rgn2Sp - VS_CtInSp )

   if (icase >= 3) then
      VS_Slope25 = ( VS_RtPwr / PC_RefSpd    ) / ( VS_RtGnSp - VS_SySp   ) !MAKE REGION 3 CONSTANT TORQUE INSTEAD OF CONSTANT POWER FOR FLOATING WT
      PC_KI      = 0.0008965149d0 !REDUCE GAINS FOR HYWIND
      PC_KP      = 0.006275604d0  !REDUCE GAINS FOR HYWIND
!!    PC_KP      = 0.012367188d0        !gains for innwind
!!    PC_KI      = 0.012367188d0/45d0   !gains for innwind
   else
      VS_Slope25 = ( VS_RtPwr / VS_RtGnSp    ) / ( VS_RtGnSp - VS_SySp   )
      PC_KI      = 0.008068634d0
      PC_KP      = 0.01882681d0 
   endif

   if ( VS_Rgn2K == 0.d0 )  then  ! .TRUE. if the Region 2 torque is flat, and thus, the denominator in the else condition is zero
      VS_TrGnSp = VS_SySp
   else                           ! .TRUE. if the Region 2 torque is quadratic with speed
      VS_TrGnSp = ( VS_Slope25 - dsqrt( VS_Slope25*( VS_Slope25 - 4.d0*VS_Rgn2K*VS_SySp ) ) )/( 2.d0*VS_Rgn2K )
   endif

  return

endif !istatus

 100 format (20000e28.17)



   ! Main control calculations:

!=======================================================================


   ! Filter the HSS (generator) speed measurement:
   ! ---------------------------------------------
   ! NOTE: This is a very simple recursive, single-pole, low-pass filter with
   !       exponential smoothing.

   ! Update the coefficient in the recursive formula based on the elapsed time
   !   since the last call to the controller:

   !a = e^-dt*fc[rad/s]

   Alpha     = dexp( ( LastTime - Time )*CornerFreq )


   ! Apply the filter:

   GenSpeedF = ( 1.d0 - Alpha )*GenSpeed + Alpha*LastGenSpeedF


!=======================================================================


   ! Variable-speed torque control:
   ! ------------------------------

   ! Compute the elapsed time since the last call to the controller:

   ElapTime = Time - LastTimeVS


   ! Only perform the control calculations if the elapsed time is greater than
   !   or equal to the communication interval of the torque controller:
   ! NOTE: Time is scaled by OnePlusEps to ensure that the contoller is called
   !       at every time step when VS_DT = DT, even in the presence of
   !       numerical precision errors.

   if ( ElapTime - VS_DT > 1.d-08 )  then

   ! Specify region and compute the generator torque:

!!!   if ( (   GenSpeedF >= VS_RtGnSp ) .OR. (  PitCom(1) >= VS_Rgn3MP ) )  then ! We are in region 3 - power is constant
      if ( (   GenSpeedF >= VS_RtGnSp ) .OR. ( LastPitCom >= VS_Rgn3MP ) )  then ! We are in region 3 - power is constant
         Region = 3.0d0
        if (icase >= 3) then
         GenTrq = VS_RtPwr/PC_RefSpd   !MAKE REGION 3 CONSTANT TORQUE INSTEAD OF CONSTANT POWER FOR FLOATING WT
        else
         GenTrq = VS_RtPwr/GenSpeedF   !Constant Power
        endif
      elseif ( GenSpeedF <= VS_CtInSp )  then                                    ! We are in region 1 - torque is zero
         Region = 1.0d0
         GenTrq = 0.0d0
      elseif ( GenSpeedF <  VS_Rgn2Sp )  then                                    ! We are in region 1 1/2 - linear ramp in torque from zero to optimal
         Region = 1.5d0
         GenTrq = VS_Slope15*( GenSpeedF - VS_CtInSp )
      elseif ( GenSpeedF <  VS_TrGnSp )  then                                    ! We are in region 2 - optimal torque is proportional to the square of the generator speed
         Region = 2.0d0
         GenTrq = VS_Rgn2K*GenSpeedF**2        
      else                                                                       ! We are in region 2 1/2 - simple induction generator transition region
         Region = 2.5d0
         GenTrq = VS_Slope25*( GenSpeedF - VS_SySp   )
      endif

   ! Saturate GenTrq using value and rate limits

      GenTrq  = dmin1( GenTrq                      , VS_MaxTq  )  ! Saturate the commanded torque using the maximum torque limit
      if ( iStatus == 0 )  LastGenTrq = GenTrq              ! Initialize the value of LastGenTrq on the first pass only
      TrqRate = ( GenTrq - LastGenTrq ) / ElapTime                ! Torque rate (unsaturated)
      TrqRate = dmin1( dmax1( TrqRate, -VS_MaxRat ), VS_MaxRat )  ! Saturate the torque rate using its maximum absolute value
      GenTrq  = LastGenTrq + TrqRate*ElapTime                     ! Saturate the commanded torque using the torque rate limit

      open (400,file='oc_control_VS.dat',access='append')
                !   1      2       3     4          5          6       7     8
      write(400,10)Time,GenTrq,TrqRate,LastGenTrq,GenSpeedF,GenSpeed,Alpha,Region

      close(400)

   endif !Variable Speed (Generator Torque)


!=======================================================================


   ! Pitch control:
   ! --------------

   ! Compute the elapsed time since the last call to the controller:

   ElapTime = Time - LastTimePC

   ! Only perform the control calculations if the elapsed time is greater than
   !   or equal to the communication interval of the pitch controller:
   ! NOTE: Time is scaled by OnePlusEps to ensure that the contoller is called
   !       at every time step when PC_DT = DT, even in the presence of
   !       numerical precision errors.

   if ( ElapTime - PC_DT > 1.d-08 )  then


   ! Compute the gain scheduling correction factor based on the previously
   !   averaged commanded pitch angle:

      GK = 1.d0/( 1.d0 + LastPitCom/PC_KK )


   ! Compute the current speed error and its integral w.r.t. time; saturate the
   !   integral term using the pitch angle limits:

      SpdErr    = GenSpeedF     - PC_RefSpd                                 ! Current speed error
      IntSpdErr = LastIntSpdErr + SpdErr*ElapTime                           ! Current integral of speed error w.r.t. time
      IntSpdErr = MIN( MAX( IntSpdErr, PC_MinPit/( GK*PC_KI ) ), &
                                       PC_MaxPit/( GK*PC_KI )      )        ! Saturate the integral term using the pitch angle limits, converted to integral speed error limits


   ! Compute the pitch commands associated with the proportional and integral
   !   gains:

      PitComP   = GK*PC_KP*   SpdErr                                        ! Proportional term
      PitComI   = GK*PC_KI*IntSpdErr                                        ! Integral term (saturated)


   ! Superimpose the individual commands to get the total pitch command;
   !   saturate the overall command using the pitch angle limits:

      PitComT   = PitComP + PitComI                                         ! Overall command (unsaturated)
      PitComT   = MIN( MAX( PitComT, PC_MinPit ), PC_MaxPit )               ! Saturate the overall command using the pitch angle limits


   ! Saturate the overall commanded pitch using the pitch rate limit:
   ! NOTE: Since the current pitch angle may be different for each blade
   !       (depending on the type of actuator implemented in the structural
   !       dynamics model), this pitch rate limit calculation and the
   !       resulting overall pitch angle command may be different for each
   !       blade.

      do K = 1,NumBl ! Loop through all blades

         PitRate(K) = ( PitComT - LastBlPitch(K) )/ElapTime                 ! Pitch rate of blade K (unsaturated)
         PitRate(K) = MIN( MAX( PitRate(K), -PC_MaxRat ), PC_MaxRat )       ! Saturate the pitch rate of blade K using its maximum absolute value
         PitCom (K) = LastBlPitch(K) + PitRate(K)*ElapTime                  ! Saturate the overall command of blade K using the pitch rate limit

         PitAcc (K) = (PitRate(K)-LastPitRate(K))/ElapTime

      enddo          ! K - all blades


      open (401,file='oc_control_PC.dat',access='append')
!                   1      2-4       5-7          8-10       11   12    
      write(401,10)Time,PitCom(1:3),PitRate(1:3),PitAcc(1:3),GK,SpdErr,&
              IntSpdErr,PitComP,PitComI,PitComT,BlPitch(1:3),LastBlPitch(1:3),GenSpeed
!                  13       14       15      16     17-19         20-22           23

      close(401)

   endif  !Pitch control


!=======================================================================

   ! This will ensure that even if pitch and/or speed  control weren't called 
   !   the previous output will be used!

   GenTrq_out          = GenTrq
   PitAng_out(1:NumBl) = PitCom     (1:NumBl)
   PitVel_out(1:NumBl) = PitRate    (1:NumBl)
   PitAcc_out(1:NumBl) = PitAcc     (1:NumBl)


10 format (50f15.4)


END Subroutine CONTROLOS3

!DLC2.4: overspeed reduce pitch slowly with rate 0.1rad/s-->overspeed
!        if omega_gen>omega_gen * 1.13 ->emergency shut down

! normal    shut down: pitch rate=2deg/s
! emergency shut down: pitch rate=4deg/s, if pitch >80deg -> brake


