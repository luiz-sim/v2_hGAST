!----------------------------------------------------------------------
 Subroutine Init_controller !_g132
!----------------------------------------------------------------------

 Use Cbeam
 Use Craft
 Use Control_1

   implicit None


!---- G132 ---
   write(*,*) 'G132 dll controller'

   I_pit_var     = 1                      ! pitch actuator 0:pitch angle demand, 1:pitch rate demand
   I_pit_kind    = 1                      ! Kind of pitch control [collective(0) or individual(1)]

   if (I_pit_var == 0) then
      write(* ,*)'Error, Pitch angle actuator not supported at the moment'
      write(10,*)'Error, Pitch angle actuator not supported at the moment'
      stop
   endif

   Gain_opt      =     1.66d0             !Optimal mode quadratic speed-torque gain [Nms²/rad²]
   Pitch_VS      =    -0.9d0   / R2D
!** No brake info was given
   TORQUE_BRAKE  =     0.d0   * RAT_GEAR  ! Nm  (transformed to LSS)
   TIME_BRAKE0   =     0.08d0             ! sec
   TIME_BRAKE1   =     1.80d0             ! sec

   MinPit        =     -2.17d0 / R2D      ! rad
   MaxPit        =     90.00d0 / R2D      ! rad
   MinPitRat     =    -10.00d0 / R2D      ! rad/s
   MaxPitRat     =     10.00d0 / R2D      ! rad/s
   MinPit_Act    =     -5.00d0 / R2D      ! rad
   MaxPit_Act    =     90.00d0 / R2D      ! rad
   PitAct_lag    =      0.2d0             ! sec

   Gen_lag       =     0.075d0            ! sec
   MinGen_Speed  =   700.00d0*PI/30d0     ! rad/s
   MaxGen_Speed  =  1383.25d0*PI/30d0     ! rad/s
   NomGen_Speed  =  1120.00d0*PI/30d0     ! rad/s
   NomGen_Torq   = 29700.00d0             ! Nm    
!  max GenTorque:41000

!102.75

!--- Mechanical & Electrical Losses
   N_mech_loss   =  8
   N_pow_loss    = 10


   Allocate ( Pow(N_pow_loss ), Loss (N_pow_loss ),&
              Trq(N_mech_loss), MLoss(N_mech_loss)   )

!----------------  MECHANICAL LOSSES ----------------------
!     INPUT_TRQ [Nm]         LOSSES_TRQ [Nm]
   Trq(1) = -113000.d0;    MLoss(1) =-24000.d0 !jim added the 1st 2 lines
   Trq(2) =       0.d0;    MLoss(2) =  5000.d0
   Trq(3) =  113000.d0;    MLoss(3) = 24000.d0
   Trq(4) =  178000.d0;    MLoss(4) = 22000.d0
   Trq(5) =  288000.d0;    MLoss(5) = 22000.d0
   Trq(6) =  918000.d0;    MLoss(6) = 30500.d0
   Trq(7) = 1350000.d0;    MLoss(7) = 36000.d0
   Trq(8) = 3060000.d0;    MLoss(8) = 77000.d0

!---------------- ELECTRICAL LOSSES -----------------------
! X = Input shaft power [W]
! Y = Loss [W]
!     PWR_INPUT [W]          PWR_LOSS [W]
   Pow( 1) =  406370.d0;    Loss( 1) =  76370.d0
   Pow( 2) =  746030.d0;    Loss( 2) =  86030.d0
   Pow( 3) = 1085240.d0;    Loss( 3) =  95240.d0
   Pow( 4) = 1423750.d0;    Loss( 4) = 103750.d0
   Pow( 5) = 1776950.d0;    Loss( 5) = 126950.d0
   Pow( 6) = 2115590.d0;    Loss( 6) = 135590.d0
   Pow( 7) = 2459370.d0;    Loss( 7) = 149370.d0
   Pow( 8) = 2798380.d0;    Loss( 8) = 158380.d0
   Pow( 9) = 3141100.d0;    Loss( 9) = 171100.d0
   Pow(10) = 3488680.d0;    Loss(10) = 188680.d0


 END Subroutine Init_controller
