!----------------------------------------------------------------------
 Subroutine Init_controller !_g126
!----------------------------------------------------------------------

 Use Cbeam
 Use Craft
 Use Control_1

   implicit None


!---- G126 ---
   write(*,*) 'G126 dll controller'

   I_pit_var     = 1                      ! pitch actuator 0:pitch angle demand, 1:pitch rate demand
   I_pit_kind    = 1                      ! Kind of pitch control [collective(0) or individual(1)]

   if (I_pit_var == 0) then
      write(* ,*)'Error, Pitch angle actuator not supported at the moment'
      write(10,*)'Error, Pitch angle actuator not supported at the moment'
      stop
   endif

   Gain_opt      =     1.29115d0          !Optimal mode quadratic speed-torque gain [Nms²/rad²]
   Pitch_VS      =     0.0d0   / R2D
!** No brake info was given
   TORQUE_BRAKE  = 23000.d0   * RAT_GEAR  ! Nm  (transformed to LSS)
   TIME_BRAKE0   =     0.08d0             ! sec
   TIME_BRAKE1   =     1.66d0             ! sec

   MinPit        =     -1.10d0 / R2D      ! rad
   MaxPit        =     87.00d0 / R2D      ! rad
   MinPitRat     =     -7.00d0 / R2D      ! rad/s
   MaxPitRat     =      7.00d0 / R2D      ! rad/s
   MinPit_Act    =     -3.00d0 / R2D      ! rad
   MaxPit_Act    =     87.00d0 / R2D      ! rad
   PitAct_lag    =      0.2d0             ! sec

   Gen_lag       =     0.075d0            ! sec
   MinGen_Speed  =   680.00d0*PI/30d0     ! rad/s
   MaxGen_Speed  =  1380.00d0*PI/30d0     ! rad/s  !qq should be checked!
   NomGen_Speed  =  1120.00d0*PI/30d0     ! rad/s
   NomGen_Torq   = 22464.00d0             ! Nm    
  !  max GenTorque:25000

!98.214

!--- Mechanical & Electrical Losses
   N_mech_loss   =  7
   N_pow_loss    = 10


   Allocate ( Pow(N_pow_loss ), Loss (N_pow_loss ),&
              Trq(N_mech_loss), MLoss(N_mech_loss)   )

!----------------  MECHANICAL LOSSES ----------------------
!     INPUT_TRQ [Nm]         LOSSES_TRQ [Nm]
   Trq(1) =  -22000.d0;    MLoss(1) =-22000.d0; !jim added the 1st 2 lines
   Trq(2) =       0.d0;    MLoss(2) =     0.d0;
   Trq(3) =   22000.d0;    MLoss(3) = 22000.d0;
   Trq(4) =  463400.d0;    MLoss(4) = 26400.d0;
   Trq(5) =  757900.d0;    MLoss(5) = 30800.d0;
   Trq(6) = 1025400.d0;    MLoss(6) = 37800.d0;
   Trq(7) = 2300000.d0;    MLoss(7) = 80000.d0;


!---------------- ELECTRICAL LOSSES -----------------------
! X = Input shaft power [W]
! Y = Loss [W]
!     PWR_INPUT [W]          PWR_LOSS [W]
   Pow( 1) =  301.43d3;    Loss( 1) =  51.43d3
   Pow( 2) =  559.45d3;    Loss( 2) =  59.45d3
   Pow( 3) =  816.30d3;    Loss( 3) =  66.30d3
   Pow( 4) = 1073.20d3;    Loss( 4) =  73.20d3
   Pow( 5) = 1330.20d3;    Loss( 5) =  80.20d3
   Pow( 6) = 1585.30d3;    Loss( 6) =  85.30d3
   Pow( 7) = 1849.90d3;    Loss( 7) =  99.90d3
   Pow( 8) = 2107.80d3;    Loss( 8) = 107.80d3
   Pow( 9) = 2369.40d3;    Loss( 9) = 119.40d3
   Pow(10) = 2634.70d3;    Loss(10) = 134.70d3


 END Subroutine Init_controller
