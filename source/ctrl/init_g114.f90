!----------------------------------------------------------------------
 Subroutine Init_controller !_g114
!----------------------------------------------------------------------

 Use Cbeam
 Use Craft
 Use Control_1

   implicit None


!---- G114 ---
   write(*,*) 'G114 dll controller'

   I_pit_var     = 1                      ! pitch actuator 0:pitch angle demand, 1:pitch rate demand
   I_pit_kind    = 1                      ! Kind of pitch control [collective(0) or individual(1)]

   if (I_pit_var == 0) then
      write(* ,*)'Error, Pitch angle actuator not supported at the moment'
      write(10,*)'Error, Pitch angle actuator not supported at the moment'
      stop
   endif

   Gain_opt      =     0.27496d0          !Optimal mode quadratic speed-torque gain [Nms²/rad²]
   Pitch_VS      =     0.d0
!** No brake info was given
   TORQUE_BRAKE  =     0.d0   * RAT_GEAR  ! Nm  (transformed to LSS)
   TIME_BRAKE0   =     0.08d0             ! sec
   TIME_BRAKE1   =     1.80d0             ! sec

   MinPit        =     -0.27d0 / R2D      ! rad
   MaxPit        =     87.00d0 / R2D      ! rad
   MinPitRat     =     -9.00d0 / R2D      ! rad/s
   MaxPitRat     =      9.00d0 / R2D      ! rad/s
   MinPit_Act    =     -5.00d0 / R2D      ! rad
   MaxPit_Act    =     87.00d0 / R2D      ! rad
   PitAct_lag    =      0.175d0           ! sec

   Gen_lag       =     0.075d0            ! sec
   MinGen_Speed  =  1050.0d0*PI/30d0      ! rad/s
   MaxGen_Speed  =  1680.0d0*PI/30d0      ! rad/s
   NomGen_Speed  =  1680.0d0*PI/30d0      ! rad/s
   NomGen_Torq   = 11982.d0               ! Nm    
!  max GenTorque:12575
!--- Mechanical & Electrical Losses
   N_mech_loss   =  8
   N_pow_loss    =  6


   Allocate ( Pow(N_pow_loss ), Loss (N_pow_loss ),&
              Trq(N_mech_loss), MLoss(N_mech_loss)   )

!----------------  MECHANICAL LOSSES ----------------------
!     INPUT_TRQ [Nm]         LOSSES_TRQ [Nm]
   Trq(1) = -138000.d0;    MLoss(1) =-24000.d0
   Trq(2) =       0.d0;    MLoss(2) =  5000.d0
   Trq(3) =  138000.d0;    MLoss(3) = 24000.d0
   Trq(4) =  624000.d0;    MLoss(4) = 28000.d0
   Trq(5) = 1030000.d0;    MLoss(5) = 40000.d0
   Trq(6) = 1305000.d0;    MLoss(6) = 48500.d0
   Trq(7) = 1450000.d0;    MLoss(7) = 53000.d0
   Trq(8) = 1594000.d0;    MLoss(8) = 57000.d0

!---------------- ELECTRICAL LOSSES -----------------------
! X = Input shaft power [W]
! Y = Loss [W]
!     PWR_INPUT [W]          PWR_LOSS [W]
   Pow( 1) =   62000.d0;    Loss( 1) =  34000.d0
   Pow( 2) =  334000.d0;    Loss( 2) =  45000.d0
   Pow( 3) =  448000.d0;    Loss( 3) =  62000.d0
   Pow( 4) = 1000000.d0;    Loss( 4) =  76000.d0
   Pow( 5) = 1670000.d0;    Loss( 5) =  94000.d0
   Pow( 6) = 2000000.d0;    Loss( 6) = 108000.d0



 END Subroutine Init_controller
