#include "lib_ctrl.f90"
!----------------------------------------------------------------------
 Subroutine control_init_el
!----------------------------------------------------------------------

 Use Cbeam

   implicit none


   if (NBLADE_el == 0) return

   TGenLSS_el   = 0.d0 !TREF * RAT_GEAR !T_HSS = T_LSS/RAT_GEAR
   TGenLSS_D_el = 0.d0
   TGenLSS_M_el = 0.d0
   TMLossLSS_el = 0.d0
   TBrakeLSS_el = 0.d0

   Allocate (  PITCH_el (NBLADE_el) );   PITCH_el(:) = 0.d0;
   Allocate ( dPITCH_el (NBLADE_el) );  dPITCH_el(:) = 0.d0;
   Allocate (ddPITCH_el (NBLADE_el) ); ddPITCH_el(:) = 0.d0;
   Allocate (   FLAP_el (NBLADE_el) );    FLAP_el(:) = 0.d0;
   Allocate (  dFLAP_el (NBLADE_el) );   dFLAP_el(:) = 0.d0;
   Allocate ( ddFLAP_el (NBLADE_el) );  ddFLAP_el(:) = 0.d0;

!------ sellect the controller to be used
      ICMOD_def = 0       !0: the controller is loaded from external subroutine (i.e. oc3, dtu10mw) or DLL
   if (ICMOD==10) then
      ICMOD_def = 1       !1: the controller is based on the generic built-in controller
      ICMOD     = 1
   endif

!0. OC3
!1. DLL
!2. DTU_v1
!3. DTU_v2

!!#ifndef CONTROLLER_TYPE
!!#define CONTROLLER_TYPE 0
!!#endif
!!
!!#if   CONTROLLER_TYPE == 0
!! write(10,*) 'OC3 Baseline controller is enabled'
!!#elif CONTROLLER_TYPE == 1
!! write(10,*) 'DLL call for controller is enabled'
!!#elif CONTROLLER_TYPE == 2
!! write(10,*) 'DTU_v1 controller is enabled'
!!#elif CONTROLLER_TYPE == 3
!! write(10,*) 'DTU_v2 controller is enabled'
!!#elif CONTROLLER_TYPE == 33
!! write(10,*) 'DTU_v2.3 controller is enabled'
!!#endif


 END Subroutine control_init_el
!--------------------------------------------------------------------------------
!
!  Subroutine : control   ------------------
!
!  Equation for the controls
!
!--------------------------------------------------------------------------------
 Subroutine control (NTIME, it, ISTEP)
!--------------------------------------------------------------------------------

 Use Cbeam, only : ICMOD, ICMOD_def, OMEGAR, TIME, TREF, NDFBT_el, NQSW, NQSP, ACCU, PITCH_Imb_el, NBLADE_el !, R2D, PI2, RAT_GEAR
 Use Cbeam, only : UT_el, UT1_el, UT2_el
 Use Cbeam, only : TGenLSS_el, TGenLSS_D_el, TGenLSS_M_el, TMLossLSS_el, TBrakeLSS_el
 Use Cbeam, only : PITCH_el, dPITCH_el, ddPITCH_el
 Use Cbeam, only :  FLAP_el,  dFLAP_el,  ddFLAP_el
 Use lib_ctrl

   implicit none

   integer, intent(in) :: NTIME, it, ISTEP

   real(8), save       :: Time_startP=-1.d0, Time_start=-1.d0
   real(8)             :: duration, c_start, c_down

   real(8)             :: array_in (50), array_out(50)
   integer             :: i, nbod
!------------------------------------------------------------------------------------------------------------------------
!---------- Control modes ----------
! 0. no action
! 1. typical variable speed/pitch controller using external controller [from dll or fortran sub-routine] --> ICMOD_def=0
! 2. Parked simulation
! 3. Stall regulated WT [only Generator Torque, and/or tip brake]
! 4. Stall regulated WT [only Generator Torque, and/or tip brake] Airbrake activation
! 5. prescribed omega, pitch from file perform.inp
!10. typical variable speed/pitch controller using the built-in controller                               --> ICMOD_def=1
!------------------------------------------------------------------------------------------------------------------------

!------ Flap Controller (before main call of Control for IPC)
!!!qq call FLAP_Control     (it,ISTEP)


   if     (ICMOD == 0) then
!qq   call pitch_change
!qq   call tip_brake
!------ Controller is deactivated. [constant omega, pitch]
      call pitch_heli
      return
   elseif (ICMOD == 2) then
!------ Parked case. Controller is deactivated. [free-free b.c. for shaft, fixed pitch]
      TGenLSS_el = 0.d0
!qq   call tip_brake
      return
   elseif (ICMOD == 3.or.ICMOD == 4) then
!------ Stall regulated case. Controller is activated. [Generator Torque, constant pitch, tip brake only for specific dlcs]
!     T = B(w-ws)*(ws/ws) = B*slip*ws -> B=T/(slip*ws):: T(w), slip(w)=(w-ws)/ws at rated conditions, Trater=Prated/wrated

!----- update
     if (it==1) &
      Time_startP= Time_start
      Time_start = Time_startP

     if (-UT1_el(NDFBT_el+NQSW) < OMEGAR *0.97d0) then !for start-up
      TGenLSS_el   =  0.d0
      TGenLSS_D_el =  0.d0
     else
        if (Time_start<0.d0) Time_start = TIME

         duration= 0.5d0

         c_start = min(1.d0,     (TIME-Time_start)/duration)
         c_down  = max(0.d0,1.d0-(TIME-Time_start)/duration)

      TGenLSS_el   = c_start*TREF * (-UT1_el(NDFBT_el+NQSW)-OMEGAR)     ! B*(w-ws), w is negative, slip~=0.07 - 0.1
      TGenLSS_D_el = c_start*TREF
     endif
!!                             1      2     3       4              5              6                                7
!!    write(354,'(20e15.5)') TIME,c_start,c_down, Time_start, Time_startP,-UT1_el(NDFBT_el+NQSW)-OMEGAR*0.97d0, TGenLSS_el
      call tip_brake
      return
   elseif (ICMOD == 5) then
      TGenLSS_el = 0.d0
      call set_omega_pitch (NTIME,it)
      return
   endif


!--- case ICMOD = 1 
!--- Selection of controller [ICMOD_def 0: external controller, 1: built-in contoller]
   if    ( ICMOD_def == 0 ) then
!     call controller (NTIME,it,ISTEP)
   elseif( ICMOD_def == 1 ) then

    if (ISTEP==1) then

!qq
!     if (it>1) return

      call ctrl_sensors  (array_in)              !- Set sensors

      if ((NTIME==1).and.(it==1)) &
      call lib_ctrl_init (array_in)              !- Initialization of the generic controller

      call lib_ctrl_main (array_in, array_out)   !- Main call

      TGenLSS_el      = array_out(    1) ! = TGenLSS
      TGenLSS_D_el    = array_out(    2) ! = TGenLSS_D
      TGenLSS_M_el    = array_out(    3) ! = TGenLSS_M
        PITCH_el(1:3) = array_out( 4: 6) ! =   pitchOUT(1:3)
       dPITCH_el(1:3) = array_out( 7: 9) ! =  dpitchOUT(1:3)
      ddPITCH_el(1:3) = array_out(10:12) ! = ddpitchOUT(1:3)
         FLAP_el(1:3) = array_out(13:15) ! =    flapOUT(1:3)
        dFLAP_el(1:3) = array_out(16:18) ! =   dflapOUT(1:3)
       ddFLAP_el(1:3) = array_out(19:21) ! =  ddflapOUT(1:3)
      TMLossLSS_el    = array_out(   22) ! = TMLossLSS
      TBrakeLSS_el    = array_out(   23) ! = TBrakeLSS

!------ set blade pitch angle, velocity and acceleration
      do nbod         = 1, NBLADE_el
         i            = NDFBT_el + NQSP + ACCU%NQSACC_el(nbod-1) + 5
         UT_el (i)    = -  PITCH_el(nbod) + PITCH_Imb_el(nbod)
         UT1_el(i)    = - dPITCH_el(nbod)
         UT2_el(i)    = -ddPITCH_el(nbod)
      enddo !nbod

    elseif (ISTEP==2) then

      call lib_ctrl_update                       !- Update and output Control vars

    endif

   endif !ICMOD_def


 END Subroutine control
!----------------------------------------------------------------------
 Subroutine ctrl_sensors (array_in)
!----------------------------------------------------------------------

 Use Cbeam
 Use Craft, only : HUB_VEL

   implicit none

   real(8), intent(   out) :: array_in(50)
   integer                 :: nbod,nbsub,nb_el,nel, NQACC
   real(8)                 :: x(3),dx(3),ddx(3),dxl(3),ddxl(3), AKSI0(3)

!--- Sensors
   real(8) :: GenSpeedLSS, dGenSpeedLSS
   real(8) :: BlPitch (3)
   real(8) :: WindVel (3)
   real(8) :: ttop_acc(3)
   real(8) :: psi     (3)
   real(8) :: M_NRT_C (2,2)

!--- generator speed
      GenSpeedLSS   = -UT1_el(NDFBT_el+NQSW)
     dGenSpeedLSS   = -UT2_el(NDFBT_el+NQSW)

!--- pitch: b_hGAST= -b_COLL -b_INDIV -b_accel + b_IMB ==> b_COLL= -b_hGAST -b_INDIV -b_accel + b_IMB
   do nbod          = 1, NBLADE_el
      NQACC         = ACCU%NQSACC_el(nbod-1) 
      BlPitch(nbod) = -( UT_el(NDFBT_el+NQSP+NQACC+5) - PITCH_Imb_el(nbod) )
   enddo

!--- Wind Velocity
      WindVel(1:3)  = HUB_VEL (1:3)

!--- tower top acceleration
   if ( NBODBTWT_el == NTOWER_el) then
      nbod          = NBLADE_el+2
      nbsub         = body(nbod)%NBODSUB_el
      nb_el         = body(nbod)%NBODTGLB_el(nbsub)
      nel           = subbody(nb_el)%NTEB_el
      AKSI0(:)      = 0.d0;
      AKSI0(2)      = subbody(nb_el)%HTA_el (nel+1)

      call beam_kinematics (nbod,nbsub,nel,AKSI0,x,dx,ddx,dxl,ddxl)
!------ Vasilis INNWIND tri-spar floater
!!!!--- Tower base Acceleration
!!! 1 accel1           = 0.d0
!!!  if (ICASE_el >=3) then
!!!   nbod             = NBLADE_el + 2
!!!   nbsub            = 1
!!!   nb_el            = body(nbod)%NBODTGLB_el(nbsub)
!!!   nel              = 1
!!!   AKSI0(:)         = 0.d0;
!!!
!!!   call beam_kinematics (nbod,nbsub,nel,AKSI0,xfl,dxfl,ddxfl,dxlfl,ddxlfl)
!!!
!!!   accel1           = ddxlfl(1)
!!!  endif
!!!
!!!!--- Tower top Acceleration
!!!   nbod             = NBLADE_el + 2
!!!   nbsub            = body   (nbod)%NBODSUB_el
!!!   nb_el            = body   (nbod)%NBODTGLB_el(nbsub)
!!!   nel              = subbody(nb_el)%NTEB_el
!!!   AKSI0(:)         = 0.d0;
!!!   AKSI0(2)         = subbody(nb_el)%HTA_el  (nel+1)
!!!
!!!   call beam_kinematics (nbod,nbsub,nel,AKSI0,xfl,dxfl,ddxfl,dxlfl,ddxlfl)
!!!
!!!   accel2           = ddxlfl(1)-accel1 !Atop-Abot
   else
      ddx     (:)   = 0.d0;
   endif
      ttop_acc(1:3) = ddx(1:3) !global (x:fore-aft,y:side-side,z:vertical)


!--- Blade root moments (flap / edge) -> pitch (out-of-plane / in-plane) -> psi (tilt / yaw) 1P, 2P
   call CYCLIC_MOM ( M_NRT_C, psi, 2 ) !nh=2

   array_in(    :)  = 999.d10;
   array_in(    1)  = TIME
   array_in(    2)  =  GenSpeedLSS
   array_in(    3)  = dGenSpeedLSS
   array_in( 4: 6)  = BlPitch (1:3  )
   array_in( 7: 9)  = WindVel (1:3  )
   array_in(10:12)  = ttop_acc(1:3  )
   array_in(13:15)  = psi     (1:3  )
   array_in(16:17)  = M_NRT_C (1,1:2)
   array_in(18:19)  = M_NRT_C (2,1:2)

   array_in(   40)  = dble(NBLADE_el)
   array_in(   41)  = DT_el
   array_in(   42)  = TPERIOD_el
   array_in(   43)  = RAT_GEAR


 END Subroutine ctrl_sensors
!-----------------------------------------------------------------------
!
! Provides for each blade:
!   M_BLD  : 3 root Moments wrt blade        c.s. (flap        /torsion/edge    ) [N]
!   M_RD   : 3 root Moments wrt rotor disk   c.s. (out-of-plane/torsion/in-plane) [N]
!   M_NR   : 3 root Moments wrt Non-rotating c.s. (tilt        /yaw    /roll    ) [N]
!** psi    : positive azimuth angle                                             [rad]
!   pitch  : positive pitch   angle                                             [rad]
!
!   M_NRT  : the total Moment wrt Non-rotating c.s.                               [N]
!** M_NRT_C: the total Tilting and Yawing moments for the NH harmonics only by
!            considering the Mout-of-plane moment                                 [N]
!
!** only these vars are communicated at the final version
!-----------------------------------------------------------------------
 Subroutine CYCLIC_MOM ( M_NRT_C, psi, NH)
!-----------------------------------------------------------------------

 Use Cbeam

   implicit none

   integer, intent(in   ) :: NH
   real(8), intent(  out) :: M_NRT_C (NH,2)
   real(8), intent(  out) :: psi  (NBLADE_el)

   real(8) :: M_BLD(NBLADE_el,3), M_RD (NBLADE_el,3), M_NR (NBLADE_el,3), M_NRT(3), pitch(NBLADE_el)
   integer :: nbod,nbsub,nb_el,nn_el,nel,nod, I1,I2, h
   real(8) :: AT0G(3,3), AT0L(3,3), AG(3,3), AL(3,3), Apsi(3,3)
   real(8) :: FG  (6)  , FL  (6)


      M_NRT  (    1:3) = 0.d0;
   do nbod             = 1, NBLADE_el
      nbsub            = 1
      nb_el            = body(nbod)%NBODTGLB_el(nbsub)
      nn_el            = body(nbod)%NNODTGLB_el(nbsub)
      nel              = 1
      nod              = 1

!----- wrt blade root after cone and pitch, before pre-curved angle
     !AT0G   (1:3,1:3) = ATP_el   (nbod,1:3,1:3) !Blade c.s. after cone, before pitch 
      AT0G   (1:3,1:3) = ATPWrG_el(nbod,1:3,1:3) !Blade c.s. after cone, before pitch 
      AT0L   (1:3,1:3) = ATPWr_el (nbod,1:3,1:3) !Blade c.s. after cone, and pitch, before pre-curved angles
      AG     (1:3,1:3) = matmul (AT0G (1:3,1:3),transf_mat(nn_el)%A_el(1:3,1:3))
      AL     (1:3,1:3) = matmul (AT0L (1:3,1:3),transf_mat(nn_el)%A_el(1:3,1:3))

      I1               = subbody(nb_el)%NDFPNACC_el ( nod-1      ) + 1
      I2               = subbody(nb_el)%NDFPNACC_el ( nod        )
      FG         (1:6) = subbody(nb_el)%AFLOC_el    ( nel, I1:I2 )
      FL         (1:6) = subbody(nb_el)%AFLOC_el    ( nel, I1:I2 )

!------ Rotate LOADS
      call LOADS_TRANS0_NEW ( FG, AG, 1, 1, 1 )
      call LOADS_TRANS0_NEW ( FL, AL, 1, 1, 1 )

      pitch (nbod    ) = -UT_el(NDFBT_el+NQSP+ACCU%NQSACC_el(nbod-1)+5)
      psi   (nbod    ) = -UT_el(NDFBT_el+NQSW) - UThub_el + PHI0(nbod)

      call ROT_MATRIX0      ( 3, -psi(nbod), Apsi )

      M_BLD (nbod,1:3) = FL (IMDOF_el(4:6))
      M_RD  (nbod,1:3) = FG (IMDOF_el(4:6))
      M_NR  (nbod,1:3) = matmul(Apsi (1:3,1:3), M_RD(nbod,1:3))
      M_NRT (     1:3) = M_NRT       (1:3)    + M_NR(nbod,1:3)
   enddo !nbod

      M_NRT_C(:,:)     = 0.d0;
   do nbod             = 1, NBLADE_el
   do h                = 1, NH
      M_NRT_C(h,1)     =  M_NRT_C(h,1) + M_RD(nbod,1) * dcos(dble(h)*psi(nbod)) !Mtilt
      M_NRT_C(h,2)     =  M_NRT_C(h,2) + M_RD(nbod,1) * dsin(dble(h)*psi(nbod)) !Myaw
   enddo !h
   enddo !nbod


 END Subroutine CYCLIC_MOM
!----------------------------------------------------------------------------------------
 Subroutine tip_brake
!----------------------------------------------------------------------------------------
!-- This is a special case used for an old a stall regulated WT for the REWIND project
!-- The sub-routines use the pre-curved angles in order to rotate the tip brake.
!-- This means that atm the blade can't be also pre-bended or swept.
!-- This assumption is valid, because the old small blades used to be straight.
!
!   Remarks:
!-- The time at which the brake is applied must be specified by the user [the controller
!-- won't detect any overspeed].
!-- At the moment the tip-brake starts at the last sub-body.
!-- The rotation is a linear ramp (t).
!----------------------------------------------------------------------------------------

 Use Cbeam

   implicit none

   integer       :: nbod, nnsub, nnsubT
   real(8)       :: TIME_fault, TIME_stop, durationBrake, durationGen, Angle_max, OverSpeed, ANGLE, c_down
   real(8), save :: TIME_tipbrake=-1.d0


!---- Settings
   TIME_stop     =  15.0d0                                                    !time for the normal shut-down
   TIME_fault    =9999.0d0                                                    !time for the network loss (zero Tgen)
   durationBrake =   0.5d0                                                    !brake duration to turn linearly from 0 to Angle_max degrees [sec]
   durationGen   =   0.5d0                                                    !Gen   duration to linearky become zero [sec]
   Angle_max     =  85.d0/R2D                                                 !set angle
   OverSpeed     =  1.08d0                                                    !overspeed, if omega>omega_nominal*overspeed --> brake activation

   if     (ICMOD /= 4               ) return
   if     (TIME < TIME_fault .and. &
           TIME < TIME_stop         ) return

   if     (TIME >= TIME_fault) then                                            !time at which generator torque is zero
      TGenLSS_el   = 0.d0
      TGenLSS_D_el = 0.d0
      if ((-UT1_el(NDFBT_el+NQSW)<OMEGAR*OverSpeed).and. &
          (TIME_tipbrake         <0.d0)                    ) return            !still airbrake is not enabled, waiting for overspeed occurance

   elseif (TIME >= TIME_stop ) then                                            !time at which normal shut-down initiates

      c_down       = max(0.d0,1.d0-(TIME-Time_stop)/durationGen)
      TGenLSS_el   = 0.d0 !TGenLSS_el   * c_down
      TGenLSS_D_el = 0.d0 !TGenLSS_D_el * c_down
   endif


   if (TIME_tipbrake<0.d0) TIME_tipbrake = TIME

   IAX_Ben_Swe = 2                                                            !transform the pre-curved angle around pitch axis

   ANGLE = -min(Angle_max, Angle_max * (TIME-TIME_tipbrake)/durationBrake)    !negative as pitch

   if (TIME_stop<0.d0) &
   ANGLE = -Angle_max

   do nbod = 1, NBLADE_el

      body(nbod)%PRECURV(           :) = 0.d0;                                !zero possible pre-bend/sweep angles
                               nnsub   = body(nbod)%NBODSUB_el                !the tip brake starts at the last sub-body
                               nnsubT  = body(nbod)%NBODSUB_el + 1            !Total nodes or each blade
      body(nbod)%PRECURV(nnsub:nnsubT) = ANGLE                                !All nodes after tip brake are rotated
   enddo

   open(122,file='stall_ctrl.dat',access='append')
      write(122,'(10e15.5)') TIME,-ANGLE*R2D,Time_tipbrake, c_down, TGenLSS_el
   close(122)
   

 END Subroutine tip_brake
!----------------------------------------------------------------------------------------
 Subroutine pitch_change
!----------------------------------------------------------------------------------------
!-- This is a special case used for thalis project in order to simulate a simple pitch step.
!
!   Remarks:
!-- Very similar sb with tip_brake used for REWIND project [tip brake]
!-- The rotation is a linear ramp (t).
!----------------------------------------------------------------------------------------

 Use Cbeam

   implicit none

   integer :: nbod,iap,idp0
   real(8) :: TIME_fault, duration, Angle_max, ANGLE
   real(8),save :: TIME_brake = -1.d0


!---- Settings
   TIME_fault = 60.0d0                                                        !time to start the pitch step [sec]
   duration   = 20.d0/360.d0*TPERIOD_el                                       !brake duration to turn linearly from 0 to Anle_max degrees [sec]
   Angle_max  =  2.d0/R2D                                                     !set angle

   if (TIME  <  TIME_fault) return

   if (TIME_brake<0.d0) TIME_brake = TIME

   ANGLE = -min(Angle_max, Angle_max * (TIME-TIME_brake)/duration)            !negative as pitch

!--- Update pitch values based on pitch demand
   do nbod   = 1, NBLADE_el
      iap    = NDFBT_el + NQSP + ACCU%NQSACC_el(nbod-1) + 5
      idp0   = 5

      UT_el  (iap) = ANGLE
      UT1_el (iap) = 0.d0
      UT2_el (iap) = 0.d0
   enddo !nbod
   

 END Subroutine pitch_change
!----------------------------------------------------------------------------------------
 Subroutine pitch_heli
!----------------------------------------------------------------------------------------
!-- This subroutine sets a priori the blade pitch of the helicopter, based on the given 
!    collective and cyclic constants and the azimuth position.
!----------------------------------------------------------------------------------------

 Use Cbeam

   implicit none

   integer :: nbod,iap
   real(8) :: q0,qc,qs,psi,psit,psitt  !wt=azimuth angle


!! nbsub = body(nbod)%NBODSUB_el !last subbody
!! nb_el = body(nbod)%NBODTGLB_el(nbsub)
!! nel   = subbody(nb_el)%NTEB_el
!! nod   = NNPE_el

!! call beam_loads (nbod,nbsub,nel,nod, FGG, FLG, FLL)

   if (IAPPL_el /= 2) return
!--- Update pitch values based on pitch demand
   do nbod         = 1, NBLADE_el
!--- settings (also a PI controller could be added here)
      q0           = PITCH0_el (nbod)
      qc           = PITCHC_el (nbod)
      qs           = PITCHS_el (nbod)

      psi          = UT_el (NDFBT_el+NQSW) + PIhalf + PHI0  (nbod) !- UThub_el
      psit         = UT1_el(NDFBT_el+NQSW)
      psitt        = UT2_el(NDFBT_el+NQSW)

              iap  = NDFBT_el + NQSP + ACCU%NQSACC_el(nbod-1) + 5
      UT_el  (iap) = q0 + qc*        dcos(psi)                        +  qs*        dsin(psi)
      UT1_el (iap) =    - qc*  psit *dsin(psi)                        +  qs*  psit *dcos(psi)
      UT2_el (iap) =    - qc*( psitt*dsin(psi) + psit**2*dcos(psi) )  +  qs*( psitt*dcos(psi) - psit**2*dsin(psi) )
   enddo !nbod


 END Subroutine pitch_heli
!----------------------------------------------------------------------------------------
 Subroutine set_omega_pitch ( NTIME, it )
!----------------------------------------------------------------------------------------
!-- This subroutine sets the blade pitch and the omega based on given prescribed timeseries.
!----------------------------------------------------------------------------------------

 Use Cbeam

   implicit none

   integer, intent(in)          :: NTIME, it
   integer, parameter           :: np=20000
   real(8), dimension(np), save :: time_rd
   real(8), dimension(np), save ::  azim_rd , dazim_rd , ddazim_rd
   real(8), dimension(np), save :: pitch_rd ,dpitch_rd ,ddpitch_rd
   real(8),                save ::  azim_int, dazim_int, ddazim_int
   real(8),                save :: pitch_int,dpitch_int,ddpitch_int
   real(8),                save :: dt_rd    , dt2_rd   , dtsqr_rd
   integer,                save :: nt_rd
   integer                      :: i, nbod


!--- Initialization
   if (NTIME<=1 .and. it<=1) then
     open(1,file="perform.inp")
     do i = 1, np
        read(1,*, end=30) time_rd(i),dazim_rd(i),pitch_rd(i) !sec, rad/s, rad
        dazim_rd(i)=dazim_rd(i)*pi/30.    ! rpm->rad/s
        pitch_rd(i)=pitch_rd(i)*pi/180.   ! deg->rad
     enddo
        write(* ,*) 'end of perform.inp was not reatched - increase np'
        write(10,*) 'end of perform.inp was not reatched - increase np'
        stop
30   close (1)
     nt_rd    = i-1
     dt_rd    = time_rd(2)-time_rd(1)
     dt2_rd   = 2.d0*dt_rd
     dtsqr_rd = dt_rd**2

        write(* ,*)'ICMOD5-nt_rd',nt_rd,dt_rd
        write(10,*)'ICMOD5-nt_rd',nt_rd,dt_rd

                  i  = 1
          azim_rd(i) = 0.d0
        ddazim_rd(i) =                    (dazim_rd(i+1)-dazim_rd(i  )) /    dt_rd
        dpitch_rd(i) =                    (pitch_rd(i+1)-pitch_rd(i  )) /    dt_rd
       ddpitch_rd(i) = (pitch_rd(i+2)-2.d0*pitch_rd(i+1)+pitch_rd(i  )) / dtsqr_rd
     do           i  = 2, nt_rd-1
          azim_rd(i) =   azim_rd(i-1)+     dazim_rd(i  )                *    dt_rd
        ddazim_rd(i) = (dazim_rd(i+1)                   -dazim_rd(i-1)) /   dt2_rd
        dpitch_rd(i) = (pitch_rd(i+1)                   -pitch_rd(i-1)) /   dt2_rd
       ddpitch_rd(i) = (pitch_rd(i+1)-2.d0*pitch_rd(i  )+pitch_rd(i-1)) / dtsqr_rd
     enddo
                  i  = nt_rd
          azim_rd(i) =                      azim_rd(i-1)+dazim_rd(i  )  *    dt_rd
        ddazim_rd(i) = (dazim_rd(i  )-     dazim_rd(i-1)              ) /    dt_rd
        dpitch_rd(i) = (pitch_rd(i  )-     pitch_rd(i-1)              ) /    dt_rd
       ddpitch_rd(i) = (pitch_rd(i  )-2.d0*pitch_rd(i-1)+pitch_rd(i-2)) / dtsqr_rd
   endif


!  if (it>1) return


!--- Set interpolated prescribed azim,pitch and corresponding time derivatives
      call LIN_INT(TIME,   azim_int,time_rd,   azim_rd, nt_rd, np)
      call LIN_INT(TIME,  dazim_int,time_rd,  dazim_rd, nt_rd, np)
      call LIN_INT(TIME, ddazim_int,time_rd, ddazim_rd, nt_rd, np)
      call LIN_INT(TIME,  pitch_int,time_rd,  pitch_rd, nt_rd, np)
      call LIN_INT(TIME, dpitch_int,time_rd, dpitch_rd, nt_rd, np)
      call LIN_INT(TIME,ddpitch_int,time_rd,ddpitch_rd, nt_rd, np)

!--- set omega
   i = NDFBT_el + NQSW
      UT_el  (i) = -  azim_int
      UT1_el (i) = - dazim_int
      UT2_el (i) = -ddazim_int

!--- set pitch
   do nbod       = 1, NBLADE_el
              i  = NDFBT_el + NQSP + ACCU%NQSACC_el(nbod-1) + 5
      UT_el  (i) = -  pitch_int
      UT1_el (i) = - dpitch_int
      UT2_el (i) = -ddpitch_int
   enddo !nbod


 END Subroutine set_omega_pitch
!
!
!
!--- select the controller
!0. OC3
!1. DLL
!2. DTU_v1
!3. DTU_v2
!33.DTU_v2.3

!!#  ifndef    CONTROLLER_TYPE
!!#     define CONTROLLER_TYPE 0
!!#  endif
!!
!!#  if   CONTROLLER_TYPE == 0
!!#     include "./ctrl/oc3_control.f90"
!!#  elif CONTROLLER_TYPE == 1
!!#     include "./ctrl/dll_control.f90"
!!#  elif CONTROLLER_TYPE == 2
!!#     include "./ctrl/10MW_control.f90"
!!#  elif CONTROLLER_TYPE == 3
!!#     include "./ctrl/dtu_control_v2.f90"
!!#  elif CONTROLLER_TYPE == 33
!!#     include "./ctrl/dtu_control_v2.f90"
!!#  endif
