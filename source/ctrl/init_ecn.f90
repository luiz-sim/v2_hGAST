!----------------------------------------------------------------------
 Subroutine Init_controller !_ecn
!----------------------------------------------------------------------

 Use Cbeam
 Use Craft
 Use Control_1

   implicit None


!---- ECN controller for tri-spar 10MW OFWT for INNWIND 
   write(*,*) 'ECN controller for tri-spar 10MW OFWT for INNWIND'

   I_pit_var     = 1                      ! pitch actuator 0:pitch angle demand, 1:pitch rate demand
   I_pit_kind    = 1                      ! Kind of pitch control [collective(0) or individual(1)]

   if (I_pit_var == 0) then
      write(* ,*)'Error, Pitch angle actuator not supported at the moment'
      write(10,*)'Error, Pitch angle actuator not supported at the moment'
      stop
   endif

   Gain_opt      =     0.0d0              !Optimal mode quadratic speed-torque gain [Nms²/rad²]
   Pitch_VS      =     0.0d0   / R2D      !$$ Both not used
!** No brake info was given
   TORQUE_BRAKE  =     0.d0   * RAT_GEAR  ! Nm  (transformed to LSS)
   TIME_BRAKE0   =     0.08d0             ! sec
   TIME_BRAKE1   =     1.80d0             ! sec

   MinPit        =    -90.00d0 / R2D      ! rad
   MaxPit        =     90.00d0 / R2D      ! rad
   MinPitRat     =    -90.00d0 / R2D      ! rad/s
   MaxPitRat     =     90.00d0 / R2D      ! rad/s
   MinPit_Act    =     -5.00d0 / R2D      ! rad
   MaxPit_Act    =     90.00d0 / R2D      ! rad
   PitAct_lag    =      0.2d0             ! sec

   Gen_lag       =     0.075d0            ! sec
   MinGen_Speed  =     0.00d0             ! rad/s
   MaxGen_Speed  =     1.00d0             ! rad/s
   NomGen_Speed  =     2.00d0             ! rad/s
   NomGen_Torq   =     1.00d7             ! Nm    
!  max GenTorque:41000

!102.75

!--- Mechanical & Electrical Losses
   N_mech_loss   =  2
   N_pow_loss    =  2


   Allocate ( Pow(N_pow_loss ), Loss (N_pow_loss ),&
              Trq(N_mech_loss), MLoss(N_mech_loss)   )

!----------------  MECHANICAL LOSSES ----------------------
!     INPUT_TRQ [Nm]         LOSSES_TRQ [Nm]
   Trq(1) =-1d10;    MLoss(1) = 0.d0 !jim added the 1st 2 lines
   Trq(2) = 1d10;    MLoss(2) = 0.d0

!---------------- ELECTRICAL LOSSES -----------------------
! X = Input shaft power [W]
! Y = Loss [W]
!     PWR_INPUT [W]          PWR_LOSS [W]
   Pow( 1) =-1d12;   Loss( 1) = 0.d0
   Pow( 2) = 1d12;   Loss( 2) = 0.d0


 END Subroutine Init_controller
