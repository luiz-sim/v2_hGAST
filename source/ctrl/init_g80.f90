!----------------------------------------------------------------------
 Subroutine Init_controller !_g80
!----------------------------------------------------------------------

 Use Cbeam
 Use Craft
 Use Control_1

   implicit None


!---- G80 ---
   write(*,*) 'G80 dll controller'

   I_pit_var     = 1                      ! pitch actuator 0:pitch angle demand, 1:pitch rate demand
   I_pit_kind    = 1                      ! Kind of pitch control [collective(0) or individual(1)]

   if (I_pit_var == 0) then
      write(* ,*)'Error, Pitch angle actuator not supported at the moment'
      write(10,*)'Error, Pitch angle actuator not supported at the moment'
      stop
   endif

   Gain_opt      =     0.136d0
   Pitch_VS      =     0.d0

   TORQUE_BRAKE  =  6000.d0   * RAT_GEAR  ! Nm  (transformed to LSS)
   TIME_BRAKE0   =     0.08d0             ! sec
   TIME_BRAKE1   =     1.80d0             ! sec

   MinPit        =     0.00d0             ! rad
   MaxPit        =    87.40d0 / R2D       ! rad
   MinPitRat     =   -17.00d0 / R2D       ! rad/s
   MaxPitRat     =    17.00d0 / R2D       ! rad/s
   MinPit_Act    =     0.00d0 / R2D       ! rad
   MaxPit_Act    =    87.40d0 / R2D       ! rad
   PitAct_lag    =     0.10d0             ! sec

   Gen_lag       =     0.10d0             ! sec
   MinGen_Speed  =    94.2477d0           ! rad/s
   MaxGen_Speed  =   198.9674d0           ! rad/s
   NomGen_Speed  =   175.929d0            ! rad/s
   NomGen_Torq   = 11795.d0               ! Nm

!--- Mechanical & Electrical Losses
   N_mech_loss   =  8
   N_pow_loss    = 13


   Allocate ( Pow(N_pow_loss ), Loss (N_pow_loss ),&
              Trq(N_mech_loss), MLoss(N_mech_loss)   )

!----------------  MECHANICAL LOSSES ----------------------
!     INPUT_TRQ [Nm]         LOSSES_TRQ [Nm]
   Trq(1) = -286000.d0;    MLoss(1) =-27000.d0
   Trq(2) =       0.d0;    MLoss(2) =  5000.d0
   Trq(3) =  286000.d0;    MLoss(3) = 27000.d0
   Trq(4) =  572000.d0;    MLoss(4) = 43000.d0
   Trq(5) =  857000.d0;    MLoss(5) = 47000.d0
   Trq(6) = 1143000.d0;    MLoss(6) = 40000.d0
   Trq(7) = 1300000.d0;    MLoss(7) = 40000.d0
   Trq(8) = 9999999.d0;    MLoss(8) = 40000.d0

!---------------- ELECTRICAL LOSSES -----------------------
! X = Input shaft power [W]
! Y = Loss [W]
!     PWR_INPUT [W]          PWR_LOSS [W]
   Pow( 1) = 0.0000d+00;    Loss( 1) = 0.1780d+05
   Pow( 2) = 0.1178d+06;    Loss( 2) = 0.1780d+05
   Pow( 3) = 0.2194d+06;    Loss( 3) = 0.1940d+05
   Pow( 4) = 0.3220d+06;    Loss( 4) = 0.2200d+05
   Pow( 5) = 0.4240d+06;    Loss( 5) = 0.2400d+05
   Pow( 6) = 0.5255d+06;    Loss( 6) = 0.2550d+05
   Pow( 7) = 0.6270d+06;    Loss( 7) = 0.2700d+05
   Pow( 8) = 0.8332d+06;    Loss( 8) = 0.3320d+05
   Pow( 9) = 0.1039d+07;    Loss( 9) = 0.3920d+05
   Pow(10) = 0.1562d+07;    Loss(10) = 0.6220d+05
   Pow(11) = 0.2080d+07;    Loss(11) = 0.8010d+05
   Pow(12) = 0.2500d+07;    Loss(12) = 0.1000d+06
   Pow(13) = 0.2500d+08;    Loss(13) = 0.1000d+06


 END Subroutine Init_controller
