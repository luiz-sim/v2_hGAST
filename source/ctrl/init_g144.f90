!----------------------------------------------------------------------
 Subroutine Init_controller !_g144
!----------------------------------------------------------------------

 Use Cbeam
 Use Craft
 Use Control_1

   implicit None


!---- G144 ---
   write(*,*) 'G144 dll controller'

   I_pit_var     = 1                      ! pitch actuator 0:pitch angle demand, 1:pitch rate demand
   I_pit_kind    = 1                      ! Kind of pitch control [collective(0) or individual(1)]

   if (I_pit_var == 0) then
      write(* ,*)'Error, Pitch angle actuator not supported at the moment'
      write(10,*)'Error, Pitch angle actuator not supported at the moment'
      stop
   endif

   Gain_opt      =      0.d0
   Pitch_VS      =      0.d0
!--- Brake info
   TORQUE_BRAKE  = 176000.d0               ! Nm LSS
   TIME_BRAKE0   =      0.10d0             ! sec
   TIME_BRAKE1   =      0.40d0             ! sec

   MinPit        =      0.00d0             ! rad
   MaxPit        =     90.00d0 / R2D       ! rad
   MinPitRat     =    -20.00d0 / R2D       ! rad/s
   MaxPitRat     =     20.00d0 / R2D       ! rad/s
   MinPit_Act    =      0.00d0             ! rad
   MaxPit_Act    =     90.00d0 / R2D       ! rad
   PitAct_lag    =      0.35d0             ! sec

   Gen_lag       =      0.10d0             ! sec
   MinGen_Speed  =     75.6d0  * RPM2RAD   ! rad/s
   MaxGen_Speed  =      0.d0               ! rad/s
   NomGen_Speed  =      0.d0               ! rad/s
   NomGen_Torq   =      0.d0               ! Nm

!--- Mechanical & Electrical Losses
   N_mech_loss   =  8
   N_pow_loss    =  7


   Allocate ( Pow(N_pow_loss ), Loss (N_pow_loss ),&
              Trq(N_mech_loss), MLoss(N_mech_loss)   )

!----------------  MECHANICAL LOSSES ----------------------
!	  INPUT_TRQ [Nm]               LOSSES_TRQ [Nm]
   Trq(1) =-6.936007d+006;    MLoss(1) =-1.134830d+005;
   Trq(2) = 0.000000d+000;    MLoss(2) = 1.134830d+005;
   Trq(3) = 6.936007d+006;    MLoss(3) = 1.134830d+005;
   Trq(4) = 6.962225d+006;    MLoss(4) = 8.726480d+004;
   Trq(5) = 6.986849d+006;    MLoss(5) = 6.264060d+004;
   Trq(6) = 7.011553d+006;    MLoss(6) = 3.793660d+004;
   Trq(7) = 7.033437d+006;    MLoss(7) = 1.605300d+004;
   Trq(8) = 7.033437d+007;    MLoss(8) = 1.605300d+004;

!---------------- ELECTRICAL LOSSES -----------------------
! X = Input shaft power [W]
! Y = Loss [W]
!      PWR_INPUT [W]               PWR_LOSS [W]
   Pow(1) = 0.000000d+000;    Loss(1) = 0.000000d+000;
   Pow(2) = 5.379700d+005;    Loss(2) = 1.879700d+005;
   Pow(3) = 2.009400d+006;    Loss(3) = 2.594000d+005;
   Pow(4) = 3.844200d+006;    Loss(4) = 3.442000d+005;
   Pow(5) = 5.728300d+006;    Loss(5) = 4.783000d+005;
   Pow(6) = 7.626000d+006;    Loss(6) = 6.260000d+005;
   Pow(7) = 7.626000d+007;    Loss(7) = 6.260000d+005;


 END Subroutine Init_controller
