!----------------------------------------------------------------------
! UINFLOW_sea            (XG,UG,AG,pdyn,Icalc)
!   UINFLOW_airy         (XG,UG,AG,pdyn      )
!      ZELEV_airy        (XG,Zita )
!      UA_AIRY_comp      (XG,time, Waveamp,Wavefreq,Wavenum,Wavephase,Waveang,depth,tide,Stretch_Fact,UINF,AINF,pd)
!      UA_AIRY_deep_comp (XG,time, Waveamp,Wavefreq,Wavenum,Wavephase,Waveang,      tide,zita        ,UINF,AINF,pd)
!   UINFLOW_Inter        (XG,UG,AG)
!   VELOCITY_AIRY        (XG,UG   )
!      U_AIRY_comp       (XG,time, Waveamp,Wavefreq,Wavenum,Wavephase,Waveang,depth,tide,UINF)
!      U_AIRY_deep_comp  (XG,time, Waveamp,Wavefreq,Wavenum,Wavephase,Waveang,      tide,UINF)
!   CURRENT_sea          (XG,UG  )
!
! ZELEV_sea              (XG,Zita)
!   ZELEV_airy           (XG,Zita)
!   ZELEV_Inter          (XG,Zita)
!----------------------------------------------------------------------
 Subroutine UINFLOW_sea ( XG, UG, AG, pdyn, Icalc )
!----------------------------------------------------------------------

 Use Hydro

   implicit none

   real(8) :: XG(3), UG(3), AG(3), pdyn
   integer :: Icalc


   UG(1:3) = 0.d0;
   AG(1:3) = 0.d0;
   pdyn    = 0.d0

   if (Icalc==1) then
!------ Calculate both Velocity/Acceleration [Valid for Morison]
      if (WaveMod /= 6) then
         call UINFLOW_AIRY  ( XG, UG, AG, pdyn )    !Airy theory
      else
         call UINFLOW_Inter ( XG, UG, AG       )    !Interpolation of read data
      endif
   else
!------ Calculate only Velocity from Airy theory without Stretching for floater case [Valid for potential approach]
      call VELOCITY_AIRY ( XG, UG )
   endif

!------ Add the Velocity from Currents
   call CURRENT_sea   ( XG, UG )


 END Subroutine UINFLOW_sea
!----------------------------------------------------------------------
 Subroutine CURRENT_sea ( XG, UG )
!----------------------------------------------------------------------

 Use Hydro

   implicit none

   real(8) :: XG(3), UG(3), Uinf


   if ( CurrMod == 1 ) then
      Uinf  = Current_veloc * Start_Fact * &
              ( max(depth_tide + max(XG(3),0.d0), 0.d0) / depth_tide )**Current_polaw  !--XG(3)<=0.d0 no wheeler stretching for the current
      UG(1) = UG(1) + Uinf  * dcos(Current_angle)
      UG(2) = UG(2) + Uinf  * dsin(Current_angle)
   endif


 END Subroutine CURRENT_sea
!----------------------------------------------------------------------
 Subroutine ZELEV_sea ( XG, Zita )
!----------------------------------------------------------------------

 Use Hydro

   implicit none

   real(8) :: XG(3), Zita


   if (WaveMod /= 6) then
      call ZELEV_AIRY  ( XG, Zita )
   else
      call ZELEV_Inter ( XG, Zita )
   endif


 END Subroutine ZELEV_sea
!----------------------------------------------------------------------
 Subroutine UINFLOW_AIRY ( XG, UG, AG, pdyn )
!----------------------------------------------------------------------

! Add velocities/accelerations for all wave components (linear theory)

 Use Hydro

   implicit none

   real(8) :: Waveamp, Wavefreq, Wavenum, Wavephase, Waveang
   real(8) :: XG(3), UG(3), AG(3), pdyn, UINF(3), AINF(3), pd, zita
   real(8) :: UG_sum(3), AG_sum(3), pd_sum, XGc(3) !c:corrected
   real(8) :: Stretch_Factor
   integer :: ifr


   Stretch_Factor = 1.d0
   zita           = 0.d0
!--- Wheeler Stretching method
   if (StretchMod==1) then
      call ZELEV_airy ( XG, zita )
      Stretch_Factor = depth_tide/(depth_tide+zita)
   endif

   XGc(1:3) = XG(1:3)
   if (XG(3)>zita) then
      write(* ,*)'SOS point of calculation is out of the water!'
      write(10,*)'SOS point of calculation is out of the water!'
      write(* ,*)XG(1),XG(2),XG(3),zita
      XGc(3)=zita
   endif


!$omp parallel shared  (XG,XGc,UG,AG,pdyn,time_h,depth_tide,tide,Stretch_Factor,spectral_freq,spectral_ampl,spectral_phase,spectral_wavenum,spectral_angle,no_spectral_comp_deep,no_spectral_comp) &
!$omp          private (UG_sum,AG_sum,pd_sum, UINF,AINF,pd, ifr, Wavefreq, Waveamp, Wavephase, Wavenum, Waveang)
   UG_sum(1:3) = 0.d0;
   AG_sum(1:3) = 0.d0;
   pd_sum      = 0.d0
!$omp do 
   do ifr = 1, no_spectral_comp_deep-1

      Wavefreq    = spectral_freq    (ifr)
      Waveamp     = spectral_ampl    (ifr)
      Wavephase   = spectral_phase   (ifr)
      Wavenum     = spectral_wavenum (ifr)
      Waveang     = spectral_angle

      call UA_AIRY_comp ( XGc, time_h, Waveamp, Wavefreq, Wavenum, Wavephase, Waveang,&
                          depth_tide, tide, Stretch_Factor, UINF, AINF, pd )

      UG_sum(1:3) = UG_sum(1:3) + UINF(1:3)
      AG_sum(1:3) = AG_sum(1:3) + AINF(1:3)
      pd_sum      = pd_sum      + pd
   enddo
!$omp end do

!$omp do 
   do ifr = no_spectral_comp_deep, no_spectral_comp

      Wavefreq    = spectral_freq    (ifr)
      Waveamp     = spectral_ampl    (ifr)
      Wavephase   = spectral_phase   (ifr)
      Wavenum     = spectral_wavenum (ifr)
      Waveang     = spectral_angle

      call UA_AIRY_deep_comp ( XGc, time_h, Waveamp, Wavefreq, Wavenum, Wavephase, &
                               Waveang,tide,zita          ,UINF, AINF, pd     )

      UG_sum(1:3) = UG_sum(1:3) + UINF(1:3)
      AG_sum(1:3) = AG_sum(1:3) + AINF(1:3)
      pd_sum      = pd_sum      + pd
   enddo
!$omp end do
!$omp critical
   UG(1:3) = UG(1:3) + UG_sum(1:3)
   AG(1:3) = AG(1:3) + AG_sum(1:3)
   pdyn    = pdyn    + pd_sum
!$omp end critical
!$omp end parallel

   UG(1:3) = UG(1:3) * Start_Fact;
   AG(1:3) = AG(1:3) * Start_Fact;
   pdyn    = pdyn    * Start_Fact * rhoGrav


 END Subroutine UINFLOW_AIRY
!-----------------------------------------------------------------------
!
!----- Subroutine :UA_AIRY_comp
!
!--    Airy Wave Velocities and Accelerations
!
!      Used in morison's equation
!
!      z=0 at free surface "level" --> so for hydro: -depth <= z <= zita (or 0)
!
!-----------------------------------------------------------------------
 Subroutine UA_AIRY_comp (XG,time, Waveamp,Wavefreq,Wavenum,Wavephase,&
                          Waveang,depth,tide,Stretch_Fact, UINF, AINF, pd )
!-----------------------------------------------------------------------

   implicit none

   real(8) :: time, Waveamp, Wavefreq, Wavenum, Waveang,Wavephase
   real(8) :: depth, XG(3), UINF(3), AINF(3), Stretch_Fact, tide, pd
   real(8) :: th, kzd,c0,c1,c11,c2


   th   = Wavenum * ( XG(1)*dcos(Waveang)+XG(2)*dsin(Waveang) ) - Wavefreq * time + Wavephase
   kzd  = Wavenum * ( XG(3)-tide+depth ) * Stretch_Fact


   c0   = Waveamp*Wavefreq/dsinh(Wavenum*depth)
   c1   = c0 * dcosh(kzd)
   c11  = c1 * dcos (th )
   c2   = c0 * dsinh(kzd)

   UINF(1) = c11 * dcos(Waveang)
   UINF(2) = c11 * dsin(Waveang)
   UINF(3) = c2  * dsin(th)

   c11  = c1 * Wavefreq * dsin(th)
   c2   = c2 * Wavefreq

   AINF(1) = c11 * dcos(Waveang)
   AINF(2) = c11 * dsin(Waveang)
   AINF(3) =-c2  * dcos(th)

   pd      = Waveamp*dcosh(kzd)/dcosh(Wavenum*depth)*dcos(th) !*rho*grav  !pd=-df/dt*rho
!  Zeta    = Waveamp * dcos(th)


 END Subroutine UA_AIRY_comp
!-----------------------------------------------------------------------
 Subroutine UA_AIRY_deep_comp (XG,time, Waveamp,Wavefreq,Wavenum,Wavephase,&
                               Waveang,tide,zita        ,UINF, AINF, pd      )
!-----------------------------------------------------------------------

   implicit none

   real(8) :: time, Waveamp, Wavefreq, Wavenum, Waveang,Wavephase, XG(3), UINF(3), AINF(3), tide
   real(8) :: zita, pd
   real(8) :: th, c0, c1


   th      = Wavenum * ( XG(1)*dcos(Waveang)+XG(2)*dsin(Waveang) ) - Wavefreq * time + Wavephase
   c0      = Waveamp*Wavefreq    * dexp(Wavenum * (XG(3)-tide-zita) )
   c1      = c0 * Wavefreq

   UINF(1) = c0 * dcos(Waveang) * dcos(th)
   UINF(2) = c0 * dsin(Waveang) * dcos(th)
   UINF(3) = c0                 * dsin(th)

   AINF(1) = c1 * dcos(Waveang) * dsin(th)
   AINF(2) = c1 * dsin(Waveang) * dsin(th)
   AINF(3) =-c1                 * dcos(th)

   pd      = Waveamp*dexp(Wavenum * (XG(3)-tide-zita) )*dcos(th) !*rho*grav


 END Subroutine UA_AIRY_deep_comp
!----------------------------------------------------------------------
 Subroutine ZELEV_AIRY ( XG, Zita )
!----------------------------------------------------------------------

! Add contribution of all  wave components (linear theory)

 Use Hydro

   implicit none

   real(8) :: Waveamp, Wavefreq, Wavenum, Waveang,Wavephase, Zita, th, XG(3)
   real(8) :: Zita_sum
   integer :: ifr


   Zita = 0.d0
!$omp parallel shared  (XG,Zita,time_h,spectral_freq,spectral_ampl,spectral_phase,spectral_wavenum,spectral_angle,no_spectral_comp) &
!$omp          private (Zita_sum, ifr, th, Wavefreq, Waveamp, Wavephase, Wavenum, Waveang)
   Zita_sum = 0.d0
!$omp do 
   do ifr = 1, no_spectral_comp

      Wavefreq  = spectral_freq    (ifr)
      Waveamp   = spectral_ampl    (ifr)
      Wavephase = spectral_phase   (ifr)
      Wavenum   = spectral_wavenum (ifr)
      Waveang   = spectral_angle

      th        = Wavenum * ( XG(1)*dcos(Waveang)+XG(2)*dsin(Waveang) ) - Wavefreq * time_h + Wavephase
      Zita_sum  = Zita_sum + Waveamp * dcos(th)
   enddo
!$omp end do
!$omp critical
   Zita = Zita + Zita_sum
!$omp end critical
!$omp end parallel

   Zita = Zita  * Start_Fact


 END Subroutine ZELEV_AIRY 
!
!
!----------------------------------------------------------------------
 Subroutine UINFLOW_Inter ( XG, UG, AG )
!----------------------------------------------------------------------
!
! Find time of current position
! Interpolate 1st: linear time       values
!    -//-     2nd: spline z position values
!
! SOS: Assumes that the read wave has only x, z components and is
!        repeated along the span direction.
!
!----------------------------------------------------------------------

 Use Hydro

   implicit none

   real(8) :: XG (3), UG (3), AG (3), X_time
   real(8) :: Var   (IZMax_wav)
   real(8) :: dVar2 (IZMax_wav)
   real(8) :: var0, tx, RAT1, RAT2, z
   integer :: it1, it2, i


!--- rotate x,y with opposite of wave angle
   X_time  = XG(1)*dcos(spectral_angle) + XG(2)*dsin(spectral_angle) !XG(1)
   z       = XG(3)
   tx      = time_h - X_time / Cvel_wav + T0_wav
   tx      = tx     - dble (int (tx / Tperiod_wav))*Tperiod_wav

   it1     = int(tx / DT_wav) + 1
   it2     = it1 + 1

   RAT2    = (tx-T_wav(it1))/(T_wav(it2)-T_wav(it1))
   RAT1    = 1.d0 - RAT2

!! write(*,*)'tx',time_h, tx , T_wav(it1),T_wav(it2)


   do i=1,3,2
!------ velocity
      Var(1:IZMax_wav) = U_wav(it1,1:IZMax_wav,i) * RAT1 + U_wav(it2,1:IZMax_wav,i) * RAT2

      call  spline1 ( Z_wav, Var, dVar2, IZMax_wav )
      call  splint1 ( Z_wav, Var, dVar2, IZMax_wav, z, var0 )

      UG(i) = var0 * Start_Fact

!------ acceleration
      Var(1:IZMax_wav) = A_wav(it1,1:IZMax_wav,i) * RAT1 + A_wav(it2,1:IZMax_wav,i) * RAT2

      call  spline1 ( Z_wav, Var, dVar2, IZMax_wav )
      call  splint1 ( Z_wav, Var, dVar2, IZMax_wav, z, var0 )

      AG(i) = var0 * Start_Fact
   enddo

!------ Change for wave direction
      UG(2) = UG(1) * dsin(spectral_angle)
      UG(1) = UG(1) * dcos(spectral_angle)
      AG(2) = AG(1) * dsin(spectral_angle)
      AG(1) = AG(1) * dcos(spectral_angle)


 END Subroutine UINFLOW_Inter
!----------------------------------------------------------------------
 Subroutine ZELEV_Inter ( XG, Zita )
!----------------------------------------------------------------------
!
! Find time of current position
! Interpolate 1st: linear time       values
!
!----------------------------------------------------------------------

 Use Hydro

   implicit none

   real(8) :: Zita, XG(3), X_time
   real(8) :: tx, RAT1, RAT2
   integer :: it1, it2


   X_time  = XG(1)*dcos(spectral_angle) + XG(2)*dsin(spectral_angle) !XG(1)
   tx      = time_h - X_time / Cvel_wav + T0_wav
   tx      = tx     - dble (int (tx / Tperiod_wav))*Tperiod_wav

   it1     = tx / DT_wav + 1
   it2     = it1 + 1

   RAT2    = (tx-T_wav(it1))/(T_wav(it2)-T_wav(it1))
   RAT1    = 1.d0 - RAT2

   Zita    = Zelev_wav(it1) * RAT1 + Zelev_wav(it2) * RAT2
   Zita    = Zita  * Start_Fact


 END Subroutine ZELEV_Inter
!
!
!
!----------------------------------------------------------------------
 Subroutine VELOCITY_AIRY ( XG, UG )
!----------------------------------------------------------------------

! Only when ICASE_el = 3
! Add velocities for all wave components (linear theory) for Morison's eq. drag term
! Used only for the floating case, in order to save computational time (no accelerations).
! Also no Stretching.
!

 Use Hydro

   implicit none

   real(8) :: Waveamp, Wavefreq, Wavenum, Wavephase, Waveang
   real(8) :: XG(3), UG(3), UINF(3), UG_sum(3)
   integer :: ifr


!$omp parallel shared  (XG,UG,time_h,depth_tide,tide,spectral_freq,spectral_ampl,spectral_phase,spectral_wavenum,spectral_angle,no_spectral_comp_deep,no_spectral_comp) &
!$omp          private (UG_sum,UINF, ifr, Wavefreq, Waveamp, Wavephase, Wavenum, Waveang)
   UG_sum(1:3) = 0.d0;
!$omp do 
   do ifr = 1, no_spectral_comp_deep-1

      Wavefreq    = spectral_freq    (ifr)
      Waveamp     = spectral_ampl    (ifr)
      Wavephase   = spectral_phase   (ifr)
      Wavenum     = spectral_wavenum (ifr)
      Waveang     = spectral_angle

      call U_AIRY_comp (XG,time_h, Waveamp,Wavefreq,Wavenum,Wavephase,Waveang,depth_tide,tide, UINF)

      UG_sum(1:3) = UG_sum(1:3) + UINF(1:3)
   enddo
!$omp end do

!$omp do 
   do ifr = no_spectral_comp_deep, no_spectral_comp

      Wavefreq    = spectral_freq    (ifr)
      Waveamp     = spectral_ampl    (ifr)
      Wavephase   = spectral_phase   (ifr)
      Wavenum     = spectral_wavenum (ifr)
      Waveang     = spectral_angle

      call U_AIRY_deep_comp (XG,time_h, Waveamp,Wavefreq,Wavenum,Wavephase,Waveang,tide, UINF)

      UG_sum(1:3) = UG_sum(1:3) + UINF(1:3)
   enddo
!$omp end do
!$omp critical
   UG(1:3) = UG(1:3) + UG_sum(1:3)
!$omp end critical
!$omp end parallel

   UG(1:3) = UG(1:3) * Start_Fact


 END Subroutine VELOCITY_AIRY
!-----------------------------------------------------------------------
 Subroutine U_AIRY_comp (XG,time, Waveamp,Wavefreq,Wavenum,Wavephase,Waveang,depth,tide, UINF)
!-----------------------------------------------------------------------

   implicit none

   real(8) :: time, Waveamp, Wavefreq, Wavenum, Waveang,Wavephase
   real(8) :: depth, XG(3), UINF(3), tide
   real(8) :: th, kzd,c0,c1,c2


   th      = Wavenum * ( XG(1)*dcos(Waveang)+XG(2)*dsin(Waveang) ) - Wavefreq * time + Wavephase
   kzd     = Wavenum * ( XG(3)-tide+depth )
   c0      = Waveamp*Wavefreq/dsinh(Wavenum*depth)
   c1      = c0 * dcosh(kzd) * dcos(th)
   c2      = c0 * dsinh(kzd)
   UINF(1) = c1 * dcos(Waveang)
   UINF(2) = c1 * dsin(Waveang)
   UINF(3) = c2 * dsin(th)


 END Subroutine U_AIRY_comp
!-----------------------------------------------------------------------
 Subroutine U_AIRY_deep_comp (XG,time, Waveamp,Wavefreq,Wavenum,Wavephase,Waveang,tide,UINF)
!-----------------------------------------------------------------------

   implicit none

   real(8) :: time, Waveamp, Wavefreq, Wavenum, Waveang,Wavephase, XG(3), UINF(3), tide
   real(8) :: th, c0


   th      = Wavenum * ( XG(1)*dcos(Waveang)+XG(2)*dsin(Waveang) ) - Wavefreq * time + Wavephase
   c0      = Waveamp*Wavefreq * dexp(Wavenum * (XG(3)-tide))
   UINF(1) = c0 * dcos(Waveang) * dcos(th)
   UINF(2) = c0 * dsin(Waveang) * dcos(th)
   UINF(3) = c0                 * dsin(th)


 END Subroutine U_AIRY_deep_comp
