!----------------------------------------------------------------------
 Module Control_1
!----------------------------------------------------------------------

 integer             , save :: N_mech_loss
 integer             , save :: N_pow_loss
 integer             , save :: I_pit_var                ! pitch actuator 0:pitch angle demand, 1:pitch rate demand
 integer             , save :: I_pit_kind               ! Kind of pitch control [collective(0) or individual(1)]
 real(8), allocatable, save :: Pow(:), Loss (:)         ! N_pow_loss
 real(8), allocatable, save :: Trq(:), MLoss(:)         ! N_mech_loss
 real(8)             , Save :: TIME_BRAKE0 , TIME_BRAKE1
 real(8)             , Save :: MinPit_Act  , MaxPit_Act
 real(8)             , Save :: MinPit      , MaxPit
 real(8)             , Save :: MinPitRat   , MaxPitRat
 real(8)             , Save :: Gen_lag     , PitAct_lag
 real(8)             , Save :: MinGen_Speed, MaxGen_Speed
 real(8)             , Save :: NomGen_Speed, NomGen_Torq
 real(8)             , Save :: TORQUE_BRAKE, Gain_opt
 real(8)             , Save :: Pitch_VS

 END Module Control_1
!----------------------------------------------------------------------
 Subroutine controller (NTIME, it, ISTEP)
!----------------------------------------------------------------------

 Use Cbeam
 Use Craft
 Use Control_1
 Use Paths

   implicit None

   integer      :: NTIME, it,nbod,NQACC,istat, i, ISTEP
   real(8),save :: GenSpeed, BlPitch(3),GenTrq_out, PitAng_out(3), PitVel_out(3), PitAcc_out(3), PitRate(3), vel0, AFB2B(6)
   real(8)      :: QTC , QTC1 , QTC2
   real(8)      :: QTCP, QTCP1, QTCP2
   real(8)      :: QTC0, QTC01, QTC02
   real(8)      :: VAR , TVLAG
   real(4),save :: Gast2dll(36), dll2Gast(9)
   real(8)      :: Pow1, Loss1
   real(8)      :: Trq1, MLoss1
   integer      :: ieq, jg, nb_el, nn, nod, I1, I2   

!--- Tower acceleration
   real(8), dimension(3) :: x,dx,ddx,dxl,ddxl, AKSI0
   integer               :: nbsub,nel !,nbod,nb_el


   goto (1,2,3,4), ISTEP

!-- Initialization:
!--------------------

 1 if ((NTIME==1).and.(it==1)) &
    call Init_controller


!-- Mechanical Losses, bassed on
!   current shaft Torque at 1st node:
!---------------------------------------

   nb_el =  body(NBLADE_el + 1)%NBODTGLB_el(1)
   Trq1  = -subbody(nb_el)%AFLOC_el ( 1, 6 )

   call lint1 ( Trq, MLoss, 1 , N_mech_loss , N_mech_loss , Trq1, MLoss1 )

   TMLossLSS_el = MLoss1
!  write(*,*)'Mech Loss',Trq1,MLoss1


   if ((ICMOD == 0).or.(it>1)) return ! dll can't run with iterations

   if (ICMOD == 2) then               ! Parked (mechanical losses are applied)
      TGenLSS_el = 0.d0
      return
   endif


!-- Variables passed to dll controller:
!------------------------------------------------

!-- current rotational speed and blade pitch angles:

   jg            =  NDFBT_el + NQSW
   GenSpeed      =-RAT_GEAR * UT1_el(NDFBT_el+NQSW)

!--- b_hGAST= -b_COLL -b_INDIV + b_IMB ==> b_COLL= -b_hGAST -b_INDIV + b_IMB
   do nbod  = 1, NBLADE_el
      NQACC = ACCU%NQSACC_el(nbod-1) 
      BlPitch(nbod) = - (UT_el (NDFBT_el+NQSP+NQACC+5)-PITCH_Imb_el(nbod) ) -PITCH_INDIV(nbod) -PITCH_accel_el
   enddo

      Gast2dll( 1) = 1.                                                            ! istatus        [0:init, 1:run, -1:last]
   if (NTIME==1) then
      Gast2dll( 1) = 0.                                                            ! istatus        [0:init, 1:run, -1:last]
      Gast2dll(13) = real(TREF)                                                    ! gen torque
      Gast2dll(15) = 1.                                                            ! gen contactor  [0:off, 1:HS or varible, 2:LS] -- dll2Gast(6)
      Gast2dll(16) = 0.                                                            ! brake flag     [0:off, 1:on]                  -- dll2Gast(7)
   endif

   Gast2dll( 2) = real(TIME)                                                       ! time
   Gast2dll( 3) = real(GenSpeed)                                                   ! genspeed
   Gast2dll( 4) = real(GenSpeed / RAT_GEAR)                                        ! rotorspeed
   Gast2dll( 5) = real(dmod ( - UT_el(NDFBT_el+NQSW) - UThub_el + PHI0(1) , PI2 )) ! rotor azimuth angle
   Gast2dll( 6) = real(BlPitch(1))                                                 ! bld1 pitch
   Gast2dll( 7) = real(BlPitch(2))                                                 ! bld2 pitch
   Gast2dll( 8) = real(BlPitch(3))                                                 ! bld3 pitch
   Gast2dll( 9) = real(NBLADE_el)                                                  ! numbler of blds
   Gast2dll(10) = real(dsqrt(HUB_VEL(1)**2+HUB_VEL(2)**2))                         ! Horizontal hub vel

!         nb_el = body(NBODBTWT_el)%NBODTGLB_el(body(NBODBTWT_el)%NBODSUB_el)
!         nn    = subbody(nb_el)%NODMB(subbody(nb_el)%NTEB_el,5)!nnod)
!         i     = ACCU%NDFTACC_el(nb_el-1) + subbody(nb_el)%NDFPNBACC_el(nn-1)
!  Gast2dll(11) = real(UT2_el(i+1))                                                ! tower top fore-aft  acceleration
!  Gast2dll(12) =-real(UT2_el(i+4))                                                ! tower top side-side acceleration

!--- tower top acceleration
   if ( NBODBTWT_el == NTOWER_el) then
      nbod     = NBODBT_el
      nbsub    = body(nbod)%NBODSUB_el
      nb_el    = body(nbod)%NBODTGLB_el(nbsub  )
      nel      = subbody(nb_el)%NTEB_el
      AKSI0(:) = 0.d0;
      AKSI0(2) = subbody(nb_el)%HTA_el (nel+1)

      call beam_kinematics (nbod,nbsub,nel,AKSI0,x,dx,ddx,dxl,ddxl)
   else
      ddx  (:) = 0.d0;
   endif
   Gast2dll(11) = real(ddxl(1))                                                    ! tower top fore-aft  acceleration
   Gast2dll(12) =-real(ddxl(3))                                                    ! tower top side-side acceleration
   
          Pow1  = dble(Gast2dll(13)*Gast2dll(3))
          call lint1 ( Pow, Loss, 1 , N_pow_loss, N_pow_loss, Pow1, Loss1 )
         !write(*,*)'Electr Loss',Pow1, Loss1          
   Gast2dll(14) = real(Pow1-Loss1  )                                               ! gen power
!taGast2dll(17) = real(Pow1        )                                               ! measured shaft power
   Gast2dll(17) = real(Trq1)*Gast2dll(4)                                           ! measured shaft power
   
   Gast2dll(18) = real(MinPit      )                                               ! min pitch angle
   Gast2dll(19) = real(MaxPit      )                                               ! max pitch angle
   Gast2dll(20) = real(MinPitRat   )                                               ! min pitch rate
   Gast2dll(21) = real(MaxPitRat   )                                               ! max pitch rate

   Gast2dll(22) = real(DT_el       )                                               ! Communication interval
   Gast2dll(23) = real(MinGen_Speed)                                               ! Minimum Generator Speed                        rad/s
   Gast2dll(24) = real(MaxGen_Speed)                                               ! Nominal Generator Speed                        rad/s  (used for G10X)
   Gast2dll(25) = real(NomGen_Speed)                                               ! Demanded generator speed above rated [Maximum] rad/s  (used for G10X)
   Gast2dll(26) = real(NomGen_Torq )                                               ! Demanded generator torque [Nominal]            Nm     (used for G10X)
   Gast2dll(27) = real(I_pit_var   )                                               ! pitch actuator           [0:pitch angle demand, 1:pitch rate demand]
   Gast2dll(28) = real(I_pit_kind  )                                               ! Kind of pitch control    [0:collective        , 1:individual       ]
   Gast2dll(29) = real(Pitch_VS    )                                               ! below rate pitch angle set point
   Gast2dll(30) = real(Gain_opt    )                                               ! optimal mode gain


!-- Set for each Blade the in-plane and out-of-plane Bending Moments at Root
   do nbod           = 1, NBLADE_el
      nb_el          = body(nbod)%NBODTGLB_el  (1     ) !nbsub=1
      AFB2B  ( 1:6 ) = subbody(nb_el)%AFLOC_el (1, 1:6)

      Gast2dll(30 + nbod) = real(AFB2B (5))                             ! Blade i root flap-wise bending Moment , i=nbod
      Gast2dll(33 + nbod) = real(AFB2B (2))                             ! Blade i root edge-wise bending Moment , i=nbod
   enddo


!-- Call the dll:
!------------------

   call DLL_Interface ( Gast2dll, dll2Gast, file_dllname, file_dllcfg )

!write(*,*) 'dll Interface not called'; stop


!-- Variables from dll controller:
!--------------------------------------

    GenTrq_out = dble(dll2Gast(1))

!--  Uncomment / set the correct time, in order to simulate the short circuit event.
!--  Generator delay integral is important, because it reduces a lot the moment.
!jim
!!! if ((TIME >= 152.46d0 - 1.d-5).and.(TIME <= 152.46d0 + 1.d-5)) then
!!!    GenTrq_out = 123582.d0
!!!    write(*,*)'Short-Circuit',TIME
!!! endif

    if     (I_pit_kind == 1) then
       PitRate(1) = dble(dll2Gast(2))
       PitRate(2) = dble(dll2Gast(3))
       PitRate(3) = dble(dll2Gast(4))
    elseif (I_pit_kind == 0) then
       PitRate(1) = dble(dll2Gast(9))
       PitRate(2) = dble(dll2Gast(9))
       PitRate(3) = dble(dll2Gast(9))
    endif
!                        Time        1
!                        dll2Gast(1) 2  ! gen Torque demand
!                        dll2Gast(2) 3  ! pit1 demand individual
!                        dll2Gast(3) 4  ! pit2 demand individual
!                        dll2Gast(4) 5  ! pit3 demand individual
!                        dll2Gast(5) 6  ! status flag [0:init, 1:run, -1:last]
!                        dll2Gast(6) 7  ! gen contactor  [0:off, 1:HS or varible, 2:LS]
!                        dll2Gast(7) 8  ! brake flag     [0:off, 1:on]
!                        dll2Gast(8) 9  ! demanded collective pitch angle
!                        dll2Gast(9)10  ! demanded collective pitch rate


   if (NTIME==1) then
      ieq             = 8
      QTCP1_el (ieq)  = GenTrq_out
      QTC1_el  (ieq)  = GenTrq_out
      QTCP1_el (9:11) = PitRate(1:3)
      QTC1_el  (9:11) = PitRate(1:3)
   endif



!-- No Drive Train Damper [allready applied in the dll]
!-- Time lag for torque demand:
!----------------------------------------------------------

   ieq    = 8
   TVLAG  = Gen_lag

   VAR    = GenTrq_out     !HSS
   QTCP   = QTCP_el  (ieq)
   QTCP1  = QTCP1_el (ieq)
   QTCP2  = QTCP2_el (ieq)
   QTC    = QTC_el   (ieq)
   QTC1   = QTC1_el  (ieq)
   QTC2   = QTC2_el  (ieq)

   call TIME_LAG  ( TVLAG, DT_el, VAR  ,&
                    QTCP , QTCP1, QTCP2,&
                    QTC  , QTC1 , QTC2 ,&
                    QTC0 , QTC01, QTC02   )

   QTC_el  (ieq) = QTC_el (ieq) + QTC0
   QTC1_el (ieq) = QTC1_el(ieq) + QTC01
   QTC2_el (ieq) = QTC2_el(ieq) + QTC02
   QTC0_el (ieq) = QTC0
   QTC01_el(ieq) = QTC01
   QTC02_el(ieq) = QTC02

   TGenLSS_el    = QTC1_el (ieq) * RAT_GEAR



!-- Pitch actuator --> Pitch rate demand
!-- Time lag for pitch rate demand
!-----------------------------------------

   do nbod   = 1, NBLADE_el
      ieq    = 9 + (nbod-1)
      TVLAG  = PitAct_lag

      VAR    =  PitRate (nbod) + PITCH_INDIV(nbod) + PITCH_accel_el
      QTCP   =  QTCP_el  (ieq)
      QTCP1  =  QTCP1_el (ieq)
      QTCP2  =  QTCP2_el (ieq)
      QTC    =  QTC_el   (ieq)
      QTC1   =  QTC1_el  (ieq)
      QTC2   =  QTC2_el  (ieq)

      call TIME_LAG  ( TVLAG, DT_el, VAR  ,&
                       QTCP , QTCP1, QTCP2,&
                       QTC  , QTC1 , QTC2 ,&
                       QTC0 , QTC01, QTC02   )

      QTC_el  (ieq) = QTC_el (ieq) + QTC0
      QTC1_el (ieq) = QTC1_el(ieq) + QTC01
      QTC2_el (ieq) = QTC2_el(ieq) + QTC02
      QTC0_el (ieq) = QTC0
      QTC01_el(ieq) = QTC01
      QTC02_el(ieq) = QTC02


!-- Calculate pitch angle, rate and acceleration, based on
!      pitch rate demand after the time lag is set

      NQACC            = ACCU%NQSACC_el(nbod-1)
      vel0             =-UTP1_el (NDFBT_el+NQSP+NQACC+5)

!-- Saturate pitch rate
         if (dabs(QTC1_el(ieq)) > MaxPitRat ) QTC1_el(ieq) = dabs(QTC1_el(ieq))/QTC1_el(ieq) * MaxPitRat
      PitAng_out(nbod) = BlPitch(nbod) + QTC1_el(ieq) * DT_el

!-- Saturate pitch Angle
         if ( PitAng_out(nbod) < MinPit_Act ) PitAng_out(nbod) = MinPit_Act
         if ( PitAng_out(nbod) > MaxPit_Act ) PitAng_out(nbod) = MaxPit_Act
      PitVel_out(nbod) = ( PitAng_out(nbod) - BlPitch(nbod) ) / DT_el
      PitAcc_out(nbod) = ( PitVel_out(nbod) - vel0          ) / DT_el


!--Update pitch values based on pitch actuator (pitch rate demands)
      i          = NDFBT_el + NQSP + ACCU%NQSACC_el(nbod-1) + 5

      UT_el  (i) = - PitAng_out (nbod) + PITCH_Imb_el(nbod)
      UT1_el (i) = - PitVel_out (nbod)
      UT2_el (i) = - PitAcc_out (nbod)
   enddo !nbod



!-- Brake Torque:
!----------------

   if ((int(dll2Gast(7))==1).and.(TIME_BRAKE<0.d0)) TIME_BRAKE = TIME
   
   if ((int(dll2Gast(7))==1)) then !.and.(GenSpeed > 75.6d0*RPM2RAD) ) then

     if     (TIME - TIME_BRAKE < TIME_BRAKE0) then 
        TBrakeLSS_el =  0.d0
     elseif (TIME - TIME_BRAKE < TIME_BRAKE1) then
        TBrakeLSS_el = TORQUE_BRAKE * ((TIME-TIME_BRAKE)-TIME_BRAKE0)/(TIME_BRAKE1-TIME_BRAKE0)
     else
        TBrakeLSS_el = TORQUE_BRAKE
     endif

   else
     
     TIME_BRAKE   = -1.d0
     TBrakeLSS_el =  0.d0
     
   endif



!-- for next call
!------------------

    Gast2dll(13) = real(QTC1_el(8))      ! gen torque HSS
    Gast2dll( 1) = dll2Gast(5)           ! status flag    [0:init, 1:run, -1:last]
    Gast2dll(15) = dll2Gast(6)           ! gen contactor  [0:off, 1:HS or varible, 2:LS]
    Gast2dll(16) = dll2Gast(7)           ! brake flag     [0:off, 1:on]



!-- Outputs
!----------------

     open (1,file='control_in.dat' ,access='append')
     open (2,file='control_out.dat',access='append')
     open (3,file='control_var.dat',access='append')
     
      write(1,10) TIME,Gast2dll(1:36)
      write(2,10) TIME,dll2Gast(1:9)
      write(3,10) TIME,QTC1_el(8:11),Pow1,Loss1,Pow1-Loss1, TGenLSS_el, Trq1, TMLossLSS_el, TBrakeLSS_el
                  ! 1          2:5    6     7   8(P_electr)      9       10      11             12
!       8      --> generator torque actual     [TIME_LAG]
!       9:11   --> individual pitch actuator   [TIME_LAG]  Pitch rate  demand

!!      write(1,10) Gast2dll(1:17),dll2Gast(1:9),QTC1_el(8:11),Pow1,Loss1,Pow1-Loss1, TGenLSS_el, Trq1, TMLossLSS_el,TBrakeLSS_el
!!                  !        1:17          18:26        27:30    31   32      33        34         35       36          37

      10 format (130f21.6)

     close (1)
     close (2)
     close (3)

   return

!-- Update for iterations
 2 return


!-- Rerun
 3 return


!-- Recall
 4 call Init_controller
   return


 END Subroutine controller

#define G80  0
#define G114 1
#define G126 2
#define G132 3
#define G144 4
#define ECN  5
#ifndef CONTROLLER_DLL
#define CONTROLLER_DLL G80
#endif

#if   CONTROLLER_DLL == G80
#     include "./ctrl/init_g80.f90"
#elif CONTROLLER_DLL == G114
#     include "./ctrl/init_g114.f90"
#elif CONTROLLER_DLL == G126
#     include "./ctrl/init_g126.f90"
#elif CONTROLLER_DLL == G132
#     include "./ctrl/init_g132.f90"
#elif CONTROLLER_DLL == G144
#     include "./ctrl/init_g144.f90"
#elif CONTROLLER_DLL == ECN
#     include "./ctrl/init_ecn.f90"
#endif

# include "./dll_interface.f90"
!#include "dll_interface_2.f90"
