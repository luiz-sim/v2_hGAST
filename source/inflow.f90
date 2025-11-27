!----------------------------------------------------------------------
! IWINDC  0 : uniform (including shear, veer, yaw, inclination and turbulance),
!         1 : read defined scenario from wind.inp                             ,
!         2 : Extreme Operating Gust                      (EOG)               ,
!         3x: Extreme Cohirent gust with Direction change (ECD) +/-(31,32)    ,
!         4x: Extreme Direction Change                    (EDC) +/-(41,42)    ,
!         5x: Extreme Wind Shear                          (EWS) +/-(51,52)ver ,
!                                                               +/-(53,54)hor
!         6 : Normal  Turbulence Model                    (NTM)
!         7 : Extreme Turbulence Model                    (ETM)
!         8 : Extreme Wind       Model                    (EWM)
!         9 : General scaling by providing s1,s2,s3 from sigma_turb.inp
!
! IGRDTYP 0 : no turb    ,
!         1 : disk       ,
!         2 : rectangular
!---------------------------------------------------------------------------------
!
!
!  Pending: 1.    turbulance   contribution to deformation [need preprocessor for deformation calculation]
!           2. OK turbulance scale -> only 1 file /realization for all bins
!           3.    tower shadow contribution to deformation
!           4.    tower shadow for downwind cases
!           5.    tower shadow for yawed    cases [blending of upwind-downwind]
!           6.    tower shadow smoother
!
!  GenUVP: ??Calculate Particle Deformation only the time they are created??
!
!----------------------------------------------------------------------
 Module mod_inflow
!----------------------------------------------------------------------
!--- Main parameters
 integer     , save :: IWINDC
 real(8)     , save :: VELHUB
 real(8)     , save :: VELHUB_sh, VELHUB_t  , VELHUB_eshy, VELHUB_eshz
 real(8)     , save :: WINDYAW  , WINDYAW_t , WINDINC
 real(8)     , save :: VEER     , SHEXP
 real(8)     , save :: TIME_GUST, TIREF_WT  , VREF_WT    , TLAM1
 real(8)     , save :: RTIP
 real(8)     , save :: HUBPOS(3)
 real(8)     , save :: dxVELHUB_t, dxVELHUB_eshy, dxVELHUB_eshz, dxWINDYAW_t   !derivatives wrt x for deformation, d()/dx=d()/dt * 1/Uref
 real(8)     , save :: Uref                                                    !Uref:a constant reference velocity for derivatives calculation
!--- Defined wind scenario
 integer, parameter :: NINP = 1000
 integer     , save :: NUMVEL
 real(8)     , save :: TINP  (NINP), UINP  (NINP)                              !Defined wind conditions [IWINDC= 1]
 real(8)     , save :: YawINP(NINP), IncINP(NINP), SexpINP(NINP)               !Additional parameters   [IWINDC=11]
!--- Tower Shadow
 integer     , save :: ISHADOW
 real(8)     , save :: Rbot_tow, Rtop_tow, H_tow
 real(8)     , save :: RG_tow(3), AT_tow(3,3), DRG_tow(3)
!--- Check points
 integer     , save :: N_check
 real(8)     , save :: X0_check(3,10), UG_check(3,10)                          !check inflow velocity (max dimension 10 points)
!--- Turbulent Inflow
 integer     , save :: IGRDTYP 
 integer     , save :: ITURB
 real(8)     , save :: TIME_turb
!--- Rotation matrices
 real(8)     , save ::  Ainc(3,3)
!--- Partial Wake
 integer     , save              :: IPartWake
 integer     , save              :: NumDeficitP
 real(8)     , save              :: LateralDist                                !lateral distance of the center of the wake wrt the current tower. if the wake is on the right then the position should be negative because positive global y-axis points to the left while looking from upwind (direction of the wind).
 real(8)     , save, allocatable :: deficitR (:), deficitU (:)
 Character*12, save              :: FileDeficit

 END Module mod_inflow
!---------------------------------------------------------------------------------
 Subroutine INFLOW_init ( IWINDC_in   , VELHUB_in  , WINDYAW_in, WINDINC_in, VEER_in     , SHEXP_in  , &
                          TIME_GUST_in, TIREF_WT_in, VREF_WT_in, RTIP_in   , HUBPOS_in   , ISHADOW_in, &
                          Rbot_tow_in , Rtop_tow_in, H_tow_in  , IGRDTYP_in, TIME_turb_in, file_wind     )
!---------------------------------------------------------------------------------

 use mod_inflow

   implicit none

   integer       :: i
   integer       :: IWINDC_in
   real(8)       :: VELHUB_in
   real(8)       :: WINDYAW_in  , WINDINC_in
   real(8)       :: VEER_in     , SHEXP_in
   real(8)       :: TIME_GUST_in, TIREF_WT_in, VREF_WT_in
   real(8)       :: RTIP_in
   real(8)       :: HUBPOS_in(3)
   integer       :: ISHADOW_in
   real(8)       :: Rbot_tow_in, Rtop_tow_in, H_tow_in
   integer       :: IGRDTYP_in
   real(8)       :: TIME_turb_in
   character*256 :: file_wind      !path of turbulent wind file
   real(8)       :: U(3)
   real(8)       :: R2D
!--- Rotation matrices
   real(8)       ::  A   (3,3),  dA   (3,3)
   real(8)       ::  Ayaw(3,3),  dAyaw(3,3)
   real(8)       ::  Atot(3,3)
   real(8)       :: ATtot(3,3)
!--- Partial Wake
   real(8)       :: ReferenceR, ReferenceU
!--- std for turbulence
   real(8)       :: SIG(3)


   IWINDC      = IWINDC_in         !wind case 0:uniform wind, 1:defined wind scenario, 2:EOG, 3x:ECD, 4x:EDC, 5x:EWS
   VELHUB      = VELHUB_in
   WINDYAW     = WINDYAW_in        !wind yaw
   WINDINC     = WINDINC_in        !wind inclination (tilt)
   SHEXP       = SHEXP_in          !shear effect exponental = 0:no shear
   VEER        = VEER_in           !wind veer:: slope of yaw angle / m wrt hub height [rad/m]
   TIME_GUST   = TIME_GUST_in      !needed for 2,3x,4x,5x         **Time to start extreme Gust
   TIREF_WT    = TIREF_WT_in       !needed for 2    4x,5x,6,7,8   **Ti reference (IEC) A: 0.16, B : 0.14, C  : 0.12
   VREF_WT     = VREF_WT_in        !needed for 2          6,7,8   **V  reference (IEC) I:50.00, II:42.50, III:37.50
   RTIP        = RTIP_in
   HUBPOS(1:3) = HUBPOS_in(1:3);
   ISHADOW     = ISHADOW_in
   Rbot_tow    = Rbot_tow_in
   Rtop_tow    = Rtop_tow_in
   H_tow       = H_tow_in
   IGRDTYP     = IGRDTYP_in        !0:no turb, 1:disk, 2:rectangular
   TIME_turb   = TIME_turb_in      !time after which the turbulent wind starts


!--- No turbulence for defined extreme scenarions EOG, ECD, EDC, EWS.
   if (IWINDC== 2.or.              &
       IWINDC==31.or.IWINDC==32.or.&
       IWINDC==41.or.IWINDC==42.or.&
       IWINDC==51.or.IWINDC==52.or.&
       IWINDC==53.or.IWINDC==54     ) IGRDTYP = 0
   if (IGRDTYP == 0)                  ITURB   = 0
   if (IGRDTYP >  0)                  ITURB   = 1
                                      SIG(:)  = 0.d0;

      write(* ,*)
      write(* ,*)'VELHUB',VELHUB
      write(* ,*)

      write(10,*)
      write(10,*)'--INFLOW INIT---' 
      write(10,*)
      write(10,*)'VELHUB    ',VELHUB
      write(10,*)'WINDYAW   ',WINDYAW * 180.d0/dacos(-1.d0)
      write(10,*)'WINDINC   ',WINDINC * 180.d0/dacos(-1.d0)
      write(10,*)'SHEXP     ',SHEXP
      write(10,*)'VEER      ',VEER    * 180.d0/dacos(-1.d0)
      write(10,*)'IWINDC    ',IWINDC
      write(10,*)'TIREF_WT  ',TIREF_WT
      write(10,*)'VREF_WT   ',VREF_WT
      write(10,*)'TIME_GUST ',TIME_GUST
      write(10,*)'ITURB     ',ITURB
      write(10,*)'IGRDTYP   ',IGRDTYP
      write(10,*)'ISHADOW   ',ISHADOW
      write(10,*)'Rbot_tow  ',Rbot_tow 
      write(10,*)'Rtop_tow  ',Rtop_tow 
      write(10,*)'H_tow     ',H_tow 
      write(10,*)


!--- Set TLAM1 needed for IEC scenarios [IWINDC >1]
      TLAM1 = 0.7d0*max(HUBPOS(3),60.d0)


   if     (IWINDC == 1.or.IWINDC == 11) then

!------ defined wind scenario
      open (88,file='wind.inp')
      read (88,*) NUMVEL
      read (88,*)
                  if (NUMVEL>NINP) then
                     write(*,*)'dimension of T,U INP exceed the limit', NUMVEL, NINP; stop
                  endif
      do i = 1, NUMVEL
         if (IWINDC == 1) read (88,*) TINP(i), UINP(i)
         if (IWINDC ==11) read (88,*) TINP(i), UINP(i), YawINP(i), IncINP(i), SexpINP(i)
      enddo
                          R2D       = 180.d0/dacos(-1.d0)
         if (IWINDC ==11) YawINP(:) = YawINP(:) / R2D;
         if (IWINDC ==11) IncINP(:) = IncINP(:) / R2D;

      close (88)

   elseif (IWINDC == 6.or.IWINDC == 7.or.IWINDC == 8) then

!------ Set the std of the turbulent box automatically, based on the wind case

      if (ITURB /= 1) then; write(*,*) 'turbulence should be included for IWINDC 6,7,8'; stop; endif;
!------ NTM Normal  Turbulence Model [IEC 6.3.1.3]
      if (IWINDC==6) SIG(1) = TIREF_WT*(0.75d0*VELHUB+5.6d0)
!------ ETM Extreme Turbulence Model [IEC 6.3.2.3]
      if (IWINDC==7) SIG(1) = 2.d0*TIREF_WT*(0.072d0*(0.1d0*VREF_WT+3.d0)*(VELHUB/2.d0-4.d0)+10.d0)
!------ EWM Extreme Wind speed Model [IEC 6.3.2.1]
      if (IWINDC==8) SIG(1) = 0.11d0*VELHUB
                     SIG(2) = 0.8d0*SIG(1)
                     SIG(3) = 0.5d0*SIG(1)

!qq                  SHEXP  = 0.2d0 (onshore), 0.14 (offshore)     !by default NTM,ETM
!qq                  SHEXP  = 0.11d0                               !by default EWM
      write(10,*)'S - TI based on IEC'
      write(10,*)'S1 - TI1  ',SIG(1), SIG(1)/VELHUB*100.d0
      write(10,*)'S2 - TI2  ',SIG(2), SIG(2)/VELHUB*100.d0
      write(10,*)'S3 - TI3  ',SIG(3), SIG(3)/VELHUB*100.d0

   elseif (IWINDC == 9) then

!------ Read sigma1,2,3 from file sigma_turb.inp
      open (88,file='sigma_turb.inp')
         read(88,*) SIG(1),SIG(2),SIG(3)
      close(88)
      write(10,*)'S - TI from input file'
      write(10,*)'S1 - TI1  ',SIG(1), SIG(1)/VELHUB*100.d0
      write(10,*)'S2 - TI2  ',SIG(2), SIG(2)/VELHUB*100.d0
      write(10,*)'S3 - TI3  ',SIG(3), SIG(3)/VELHUB*100.d0

   elseif (IWINDC == 91) then
!------ Read TI instead of TI_REF 
      SIG(1) = TIREF_WT*VELHUB
      SIG(2) = 0.8d0*SIG(1)
      SIG(3) = 0.5d0*SIG(1)
      write(10,*)'S - TI with given TI'
      write(10,*)'S1 - TI1  ',SIG(1), SIG(1)/VELHUB*100.d0
      write(10,*)'S2 - TI2  ',SIG(2), SIG(2)/VELHUB*100.d0
      write(10,*)'S3 - TI3  ',SIG(3), SIG(3)/VELHUB*100.d0
   endif



!--- Initialize turbulant wind - all proc enter but only one reads
   if (ITURB == 1) call INWINDINT ( file_wind, VELHUB, TIME_turb, IGRDTYP, SIG )


!--- Set check points for output
   N_check             =  1
   if (ITURB == 0) &
   N_check             =  1
   X0_check (1:3,1:10) =  0.d0;
   X0_check (1:3,1   ) =  0.0d0; !hub
   X0_check (2  ,2   ) = 89.9d0; ![0., 89.9,  0. ]
   X0_check (2  ,3   ) =-89.9d0; ![0.,-89.9,  0. ]
   X0_check (3  ,4   ) =-89.9d0; ![0.,  0. ,-89.9]
   X0_check (3  ,5   ) = 89.9d0; ![0.,  0. , 89.9]


!--- set constant Ainc
   call ROT_MATRIX1 ( 2, -WINDINC, A, dA )!, AQ02 )
   Ainc(1:3,1:3) = A(1:3,1:3);

!--- set Uref
   call ROT_MATRIX1 ( 3, -WINDYAW, Ayaw, dAyaw )

       Atot(:,:) = matmul(Ayaw, Ainc);  ATtot(:,:) = transpose( Atot(:,:) );
   U(1  ) = VELHUB
   U(2:3) = 0.d0;
   U(1:3) = matmul(Atot(1:3,1:3),U(1:3))
   Uref   = U(1)
   write(10,*)
   write(10,*) 'Uref',Uref


!--- Partial Wake
   FileDeficit = "partwake.inp"
   open  (77,file=FileDeficit)
    read (77,*,end=1,err=1) NumDeficitP

    goto 2

 1     IPartWake   = 0  !Partial Wake is disabled
       NumDeficitP = 0
#  ifndef    OS
#     define OS 1
#  endif
#  if   OS == 0
       call system ('rm -f partwake.inp')
#  elif OS == 1
       call system ('del partwake.inp')
#  endif

 2  if (NumDeficitP > 0) then
       IPartWake   = 1  !Partial Wake is enabled

       Allocate ( deficitR (NumDeficitP) )
       Allocate ( deficitU (NumDeficitP) )

       read (77,*) LateralDist, ReferenceR, ReferenceU
       read (77,*) !**         radius r/R  deficit1..N


       do i = 1, NumDeficitP
          read (77,*)   deficitR(i),deficitU(i)
          deficitR(i) = deficitR(i)*ReferenceR
          deficitU(i) = deficitU(i)*ReferenceU
       enddo
       if (WINDYAW /= 0.) then; write(*,*) 'Partial wake with wind yaw is not supported'        , WINDYAW; stop; endif;
       if (WINDINC /= 0.) then; write(*,*) 'Partial wake with wind inclination is not supported', WINDINC; stop; endif;
    endif
   close (77)
!AVATAR
    if     (NumDeficitP < 0) then
       IPartWake   = 2
    endif

      write(10,*)'IPartWake  ', IPartWake
   if (IPartWake==1) then
      write(10,*)'LateralDist', LateralDist
      write(10,*)'ReferenceR ', ReferenceR
      write(10,*)'ReferenceU ', ReferenceU
   endif


 END Subroutine INFLOW_init
!-------------------------------------------------------------------------
!
! Subroutine: INFLOW
!
!
! Sets VELHUB_sh
!      VELHUB_t
!      VELHUB_eshy
!      VELHUB_eshz
!      WINDYAW_t
!    dxVELHUB_t
!    dxVELHUB_eshy
!    dxVELHUB_eshz
!    dxWINDYAW_t
!
! Stores RG_tow, AT_tow, DRG_tow for tower shadow calculation in case of floating WT
!
! calcualate UG_check [i.e. HUB_VEL]
!
! Outputs instant hub velocity, wind yaw and inclination [atm is constant]
! needed by BEM.
! 
!------------------------------------------------------------------------
 Subroutine INFLOW (time       , RG_tow_in  , AT_tow_in  , DRG_tow_in, &
                    HUB_VEL_out, VELHUBT_out, WINDYAW_out, WINDINC_out   )
!----------------------------------------------------------------------

 use mod_inflow

   implicit none


   real(8),intent(in ) :: time, RG_tow_in(3), AT_tow_in(3,3), DRG_tow_in(3)
   real(8),intent(out) :: HUB_VEL_out(3), VELHUBT_out, WINDYAW_out, WINDINC_out

   real(8) :: TIREF,VREF,SIG1,Vel,DIAM,Tcg,Vcg,Thcg,TS1,TS2,VGUST, Vsh, dxVsh, adir,Vreference

   real(8) :: X0(3), UG(3), velhub_inst, windyaw_inst
   real(8) :: uturb (3)
   real(8) :: pi, R2D
   integer :: ipoint
!--- Rotation matrices
   real(8) ::  A   (3,3),  dA   (3,3)
   real(8) ::  Ayaw(3,3),  dAyaw(3,3)
   real(8) ::  Atot(3,3)
   real(8) :: ATtot(3,3)
!--- Partial Wake
   real(8) :: RposPWake, U_def
!AVATAR
   real(8) :: RAVATAR,XupWT,YupWT,ZupWT,RrelPW

              
   pi               = dacos(-1.d0)
   R2D              = 180.d0/pi
   RG_tow (1:3    ) = RG_tow_in (1:3    )
   AT_tow (1:3,1:3) = AT_tow_in (1:3,1:3)
   DRG_tow(1:3    ) = DRG_tow_in(1:3    )
   Vreference       = VELHUB

!--------------------------------------------------------------------
!-- Set Hub Velocity based on Wind-Condition Parameter [IWINDC]
!     VELHUB0 = VELHUB_sh *(XAB(3)/HUBPOS(3))**SHEXP + VELHUB_t             &
!                                                    + VELHUB_eshy* XAB(2)  &
!                                                    + VELHUB_eshz*(XAB(3)-HUBpos(3))
!--------------------------------------------------------------------

!--- uniform inflow [initial values] or IWINDC=0
   VELHUB_sh   = VELHUB 
   VELHUB_t    = 0.d0
   VELHUB_eshy = 0.d0
   VELHUB_eshz = 0.d0
   WINDYAW_t   = 0.d0
 dxVELHUB_t    = 0.d0
 dxVELHUB_eshy = 0.d0
 dxVELHUB_eshz = 0.d0
 dxWINDYAW_t   = 0.d0
   adir        = 1.d0               !+:31,41,51,53
   if (mod(IWINDC,2)==0) adir=-1.d0 !-:32,42,52,54
!w write(*,*)'check adir',mod(IWINDC,2), adir


   if     (IWINDC == 1) then

!------ defined wind scenario
      call LIN_INT (TIME, VELHUB_sh, TINP,    UINP, NUMVEL, NINP)

   elseif (IWINDC ==11) then

!------ detailed defined wind scenario
      call LIN_INT (TIME, VELHUB_sh, TINP,    UINP, NUMVEL, NINP)
      call LIN_INT (TIME, WINDYAW_t, TINP,  YawINP, NUMVEL, NINP)
      call LIN_INT (TIME, WINDINC  , TINP,  IncINP, NUMVEL, NINP)
      call LIN_INT (TIME, SHEXP    , TINP, SexpINP, NUMVEL, NINP)

!------ set Ainc (not constant anymore)
      call ROT_MATRIX1 ( 2, -WINDINC, A, dA )
      Ainc(1:3,1:3) = A(1:3,1:3);

!--- set Uref
   call ROT_MATRIX1 ( 3, -WINDYAW, Ayaw, dAyaw )

       Atot(:,:) = matmul(Ayaw, Ainc);  ATtot(:,:) = transpose( Atot(:,:) );

   elseif (IWINDC == 2) then

!------ EOG Extreme operating gust [IEC 6.3.2.2]
      TIREF = TIREF_WT
      VREF  = VREF_WT
      SIG1  = TIREF*(0.75d0*VELHUB+5.6d0) !NTM
      Vel   = 1.12d0*VREF !0.8*1.4=1.12   !EWM steady Ve1
      DIAM  = 2.d0*RTIP
      Tcg   = 10.5d0
      TS1   = TIME_GUST
      TS2   = TS1 + Tcg

      VGUST = dmin1 ( 1.35d0*(Vel-VELHUB), 3.30d0*SIG1/(1.d0 + 0.1d0*(DIAM/TLAM1)) ) !defined at hub-height

      if (TIME >= TS1.and.TIME <= TS2) then

         VELHUB_t  = -0.37d0*VGUST*  dsin(3.d0*pi*(TIME-TS1)/Tcg)*(1.d0-dcos(2.d0*pi*(TIME-TS1)/Tcg))
                                                                                                                
       dxVELHUB_t  = -0.37d0*VGUST*( dcos(3.d0*pi*(TIME-TS1)/Tcg)*(1.d0-dcos(2.d0*pi*(TIME-TS1)/Tcg))*3.d0*pi/Tcg &
                                    +dsin(3.d0*pi*(TIME-TS1)/Tcg)*      dsin(2.d0*pi*(TIME-TS1)/Tcg) *2.d0*pi/Tcg  ) / Vreference
      endif

   elseif ( IWINDC == 31.or.IWINDC == 32 ) then

!------ ECD Extreme coherent gust with direction change positive(31) or negative(32) [IEC 6.3.2.5]
      Tcg   =  10.d0
      Vcg   =  15.d0
      TS1   = TIME_GUST
      TS2   = TS1 + Tcg
      Thcg  = 720.d0/max(VELHUB,4d0)/R2D * adir

      if     (TIME > TS2) then
         VELHUB_t  = Vcg
         WINDYAW_t = Thcg
      elseif (TIME > TS1) then
         VELHUB_t  = 0.5d0*Vcg *(1.d0-dcos(pi*(TIME-TS1)/Tcg))
       dxVELHUB_t  = 0.5d0*Vcg *      dsin(pi*(TIME-TS1)/Tcg) *pi/Tcg / Vreference
         WINDYAW_t = 0.5d0*Thcg*(1.d0-dcos(pi*(TIME-TS1)/Tcg))
       dxWINDYAW_t = 0.5d0*Thcg*      dsin(pi*(TIME-TS1)/Tcg) *pi/Tcg / Vreference
      endif

   elseif ( IWINDC == 41.or.IWINDC == 42 ) then 

!------ EDC+/- Extreme Direction Change [IEC 6.3.2.4]
      TIREF = TIREF_WT
      SIG1  = TIREF*(0.75d0*VELHUB+5.6d0)
      DIAM  = 2.d0*RTIP
      Tcg   =  6.d0
      TS1   = TIME_GUST
      TS2   = TS1 + Tcg
      Thcg  = 4.d0*datan2(SIG1,(VELHUB*(1.d0+0.1d0*(DIAM/TLAM1)))) * adir

      if     (TIME > TS2) then
         WINDYAW_t = Thcg
      elseif (TIME > TS1) then
         WINDYAW_t = 0.5d0*Thcg*(1.d0-dcos(pi*(TIME-TS1)/Tcg))
       dxWINDYAW_t = 0.5d0*Thcg*      dsin(pi*(TIME-TS1)/Tcg) *pi/Tcg / Vreference
      endif

   elseif ( IWINDC == 51.or.IWINDC == 52.or.IWINDC == 53.or.IWINDC == 54 ) then

!------ EWS 51:ver+, 52:ver-, 53:hor+, 54:hor- [IEC 6.2.3.6]
!qq   SHEXP = 0.2d0                              !by default
      TIREF = TIREF_WT
      SIG1  = TIREF*(0.75d0*VELHUB+5.6d0)
      DIAM  = 2.d0*RTIP
      Tcg   = 12.d0
      TS1   = TIME_GUST
      TS2   = TS1 + Tcg
      Vsh   =(2.5d0+1.28d0*SIG1*(DIAM/TLAM1)**0.25d0)*(1.d0-dcos(2.d0*pi*(TIME-TS1)/Tcg))/DIAM * adir
    dxVsh   =(2.5d0+1.28d0*SIG1*(DIAM/TLAM1)**0.25d0)*      dsin(2.d0*pi*(TIME-TS1)/Tcg) /DIAM * adir * 2.d0*pi/Tcg / Vreference

      if ((TIME >= TS1).and.(TIME <= TS2)) then
         if     (IWINDC == 51.or.IWINDC == 52) then  !vertical   +/-
            VELHUB_eshz =   Vsh
          dxVELHUB_eshz = dxVsh
         else!if(IWINDC == 53.or.IWINDC == 54) then  !horizontal +/-
            VELHUB_eshy =   Vsh
          dxVELHUB_eshy = dxVsh
         endif
      endif

   endif !IWINDC


!-- Calculate and Store Hub velocity, and Velocity at defined fixed points for check (the shear effect is not considered here)
!-----------------------------------------------------------------------------------------------------------------------------

!---- check the box, do not apply shear
   do ipoint       = 1, N_check
      velhub_inst  = VELHUB_sh + VELHUB_t
      windyaw_inst = WINDYAW   + WINDYAW_t
      uturb(1:3)   = 0.d0;

!------ set rotation matrices
      call ROT_MATRIX1 ( 3, -windyaw_inst, Ayaw, dAyaw )

       Atot(:,:) = matmul(Ayaw, Ainc);  ATtot(:,:) = transpose( Atot(:,:) );

      if (ITURB.eq.1) then
!--------- X0 point to interpolate wind speed
         X0  (1:3) = matmul( ATtot(1:3,1:3), X0_check(1:3,ipoint) )
     
!qq for partial wake
!!       X0(2)=X0(2)-105.d0
         call VELINWIND (X0,uturb,TIME ,IGRDTYP)
      endif

!qq   uturb(1:3) = 0.d0;

      UG(1  ) = uturb(1  ) + velhub_inst
      UG(2:3) = uturb(2:3);

!------ Partial Wake
      if (IPartWake == 1) then
       
         RposPWake  = dsqrt((X0(2)-LateralDist)**2+X0(3)**2)

         call LIN_INT( RposPWake, U_def,deficitR,deficitU, NumDeficitP, NumDeficitP )

         UG(1)      = UG(1) + U_def
      endif
!AVATAR
      if (IPartWake == 2) then
!------ Control point coordinates with respect to hub (used for turbulent wind)
         X0  (1:3) = X0_check (1:3,ipoint) - HUBPOS(1:3)
         RAVATAR=102.88; XupWT=-10.*RAVATAR; YupWT=RAVATAR; ZupWT=0.
         RrelPW=dsqrt((X0(2)-YupWT)**2+(X0(3)-ZupWT)**2)
         UG(1)     = VELHUB_sh* &
         (0.65+tanh((RrelPW/1.08/RAVATAR)**2.5)-0.65*tanh((RrelPW/0.75/RAVATAR)**1.8))
      endif

      UG_check (1:3,ipoint) = matmul( Atot(1:3,1:3),UG(1:3) )
   enddo!ipoint
 
   HUB_VEL_out(1:3) = UG_check (1:3,1)
   VELHUBT_out      = velhub_inst
   WINDYAW_out      = windyaw_inst
   WINDINC_out      = WINDINC


 END Subroutine INFLOW
!------------------------------------------------------------------------
!---- Subroutine : UINFLOW
!
!--   Determine the wind inflow velocity UG at point XG defined wrt the 
!     global c.s.
!
!     velhub_inst =  VELHUB_sh *(XAB(3)/HUBPOS(3))**SHEXP   ! ==velhub_log
!                  + VELHUB_t                            &  ! defined based on IEC for EOG or ECD
!                  + VELHUB_eshy* XAB(2)                 &  !   //      //      //     EWS
!                  + VELHUB_eshz*(XAB(3)-HUBpos(3))         !   //      //      //     EWS       
!
!     windyaw_inst = WINDYAW + WINDYAW_t + (z-z_hub)*Veer   ! WINDYAW_t is defined based on IEC for EDC or ECD
!     Yaw          = Yaw0    + Yawt(t)   + (z-z_hub)*Veer       
!
!     WINDYAW       :initial mean wind yaw
!     WINDYAW_t     :time dependant term, based on IWINDC
!     windyaw_inst  :instant wind yaw including VEER
!     Veer          :slope of yaw angle / m wrt hub height [rad/m]
!
!
!     Attention: Any modification should be also done in DINFLOW!!
!
!------------------------------------------------------------------------
 Subroutine UINFLOW (time, XG, UG, include_turb, include_tshadow)
!------------------------------------------------------------------------

 use mod_inflow

   implicit none

   integer,intent(in ) :: include_turb     !flag to include turbulance contribution if =1
   integer,intent(in ) :: include_tshadow  !flag to include tower shadow effect     if =1
   real(8),intent(in ) :: XG(3), time      !Control point coordinates wrt to the absolute frame
   real(8),intent(out) :: UG(3)

   real(8) :: X_tow(3), U_tow(3)
   real(8) :: X0   (3)
   real(8) :: uturb(3)
   real(8) :: velhub_log, velhub_inst, windyaw_inst
!--- Rotation matrices
   real(8) ::  Ayaw(3,3),  dAyaw(3,3)
   real(8) ::  Atot(3,3)
   real(8) :: ATtot(3,3)
!--- Partial Wake
   real(8) :: RposPWake, U_def
!AVATAR
   real(8) :: RAVATAR,XupWT,YupWT,ZupWT,RrelPW


!--- comes from GenUVP qq
!P if (XG(3) > 0.d0.and.SHEXP > 0.d0) then

!     write (*,'(a,9f15.3)')'uinf', VELHUB_sh, XG,HUBPOS,SHEXP 
      velhub_log = VELHUB_sh *( XG (3)/max(1.d-05,HUBPOS(3)) )**SHEXP 
!P else
!P    velhub_log = 0.1d0
!P endif

!--- Include Shear Effect
   velhub_inst  =   velhub_log                       &
                  + VELHUB_t                         &
                  + VELHUB_eshy* XG (2)              &
                  + VELHUB_eshz*(XG (3)-HUBPOS(3))        
                          
   windyaw_inst =   WINDYAW                          &
                  + WINDYAW_t                        &
                  + (XG (3)-HUBPOS(3))*VEER

   uturb(1:3)   = 0.d0;

!--- set rotation matrices
   call ROT_MATRIX1 ( 3, -windyaw_inst, Ayaw, dAyaw )

    Atot(:,:) = matmul(Ayaw, Ainc);  ATtot(:,:) = transpose( Atot(:,:) );
 
!--- Include Turbulence
   if (ITURB == 1 .and. include_turb == 1) then
!------ Control point coordinates with respect to hub (used for turbulent wind)
      X0  (1:3) = XG (1:3) - HUBPOS(1:3)

!------ Rotate hub centered coordinates for yaw and inclination
      X0  (1:3) = matmul( ATtot(1:3,1:3), X0 (1:3) )

!qq for partial wake
!!    X0(2)=X0(2)-105.d0
      call VELINWIND (X0, uturb, time, IGRDTYP)
   endif !ITURB

   UG(1  ) = uturb(1  ) + velhub_inst
   UG(2:3) = uturb(2:3);
!--- Rotate Wind for inclination and wind yaw
   UG(1:3) = matmul( Atot(1:3,1:3),UG(1:3) )


!--- Partial Wake
   if (IPartWake == 1) then
!------ Control point coordinates with respect to hub (used for turbulent wind)
      X0  (1:3) = XG (1:3) - HUBPOS(1:3)

      RposPWake = dsqrt((X0(2)-LateralDist)**2+X0(3)**2)

      call LIN_INT( RposPWake, U_def,deficitR,deficitU, NumDeficitP, NumDeficitP )

      UG(1)     = UG(1) + U_def
   endif
!AVATAR
   if (IPartWake == 2) then
!------ Control point coordinates with respect to hub (used for turbulent wind)
      X0  (1:3) = XG (1:3) - HUBPOS(1:3)
      RAVATAR=102.88; XupWT=-10.*RAVATAR; YupWT=RAVATAR; ZupWT=0.
      RrelPW=dsqrt((X0(2)-YupWT)**2+(X0(3)-ZupWT)**2)
      UG(1)     = VELHUB_sh* &
      (0.65+tanh((RrelPW/1.08/RAVATAR)**2.5)-0.65*tanh((RrelPW/0.75/RAVATAR)**1.8))
   endif


!------ Control point coordinates with respect to the tower frame (used for tower shadow)
   X_tow (1:3) = matmul(AT_tow(1:3,1:3),XG(1:3)-RG_tow(1:3))
   U_tow (1:3) = DRG_tow(1:3)

!--- Include tower shadow effect (modify UG inside)
   if (include_tshadow==1) &
   call TOWER_SHADOW (X_tow, U_tow, UG)


 END Subroutine UINFLOW
!------------------------------------------------------------------------
 Subroutine DINFLOW_dp (time, XG, UG, DG, include_turb, include_tshadow)
!----------------------------------------------------------------------

 use mod_inflow

   implicit none

   integer,intent(in ) :: include_turb     !flag to include turbulance contribution if =1
   integer,intent(in ) :: include_tshadow  !flag to include tower shadow effect     if =1
   real(8),intent(in ) :: XG(3), time      !Control point coordinates wrt to the absolute frame
   real(8),intent(out) :: UG(3), DG(3,3)

   real(8) :: X_tow(3), U_tow(3)
   real(8) :: X0   (3), Uo   (3)  , dUo  (3)
   real(8) :: uturb(3), dturb(3,3)
   real(8) :: velhub_log, velhub_inst, windyaw_inst, dzvelhub_log
!--- Rotation matrices
   real(8) ::  Ayaw(3,3),  dAyaw(3,3)
   real(8) ::  Atot(3,3),  dAtot(3,3)
   real(8) :: ATtot(3,3), dATtot(3,3)
!--- Partial Wake
   real(8) :: RposPWake, U_def
!AVATAR
   real(8) :: RAVATAR,XupWT,YupWT,ZupWT,RrelPW


!--- Partial Wake
!  The partial wake is added in GenUVP for ISOPE paper, but 2 simplifications are made.
!  1. Deformation terms dU/dy = dU/dz = 0 (because i.e. dU/dy=dU/dr * dr/dy=dU/dr * 1/cosfi which becomes singular at +-pi/2)
!  2. The deficit is assumed constant along the wake equal to the deficit at the rotor.
!qqif (IPartWake == 1) then; write(*,*) 'GenUVP with prescribed partial wake is not supported'; stop; endif;
      

!P if (XG(3) > 0.d0.and.SHEXP > 0.d0)  then
      velhub_log = VELHUB_sh *(XG (3)/max(1.d-05,HUBPOS(3)))** SHEXP 
    dzvelhub_log = VELHUB_sh *(XG (3)/max(1.d-05,HUBPOS(3)))**(SHEXP-1.d0) * SHEXP/max(1.d-05,HUBPOS(3))
!P    velhub_log = VELHUB_sh *(XG (3)/HUBPOS(3))** SHEXP 
!P  dzvelhub_log = VELHUB_sh *(XG (3)/HUBPOS(3))**(SHEXP-1.d0) * SHEXP/HUBPOS(3)
!P else
!P    velhub_log = 0.1d0
!P  dzvelhub_log = 0.0d0
!P endif

!--- instant velocity and yaw angle
   velhub_inst  =   velhub_log                       &
                  + VELHUB_t                         &
                  + VELHUB_eshy* XG (2)              &
                  + VELHUB_eshz*(XG (3)-HUBPOS(3))        
                          
   windyaw_inst =   WINDYAW                          &
                  + WINDYAW_t                        &
                  + (XG (3)-HUBPOS(3))*VEER

   uturb(1:3    ) = 0.d0;
   dturb(1:3,1:3) = 0.d0;

!--- set rotation matrices
   call ROT_MATRIX1 ( 3, -windyaw_inst, Ayaw, dAyaw )

    Atot(:,:) = matmul( Ayaw, Ainc);  ATtot(:,:) = transpose(  Atot(:,:) );
   dAtot(:,:) = matmul(dAyaw, Ainc); dATtot(:,:) = transpose( dAtot(:,:) );
 
!--- Include Turbulence
   if (ITURB == 1 .and. include_turb == 1) then
!------ Control point coordinates with respect to hub (used for turbulent wind)
      X0  (1:3) = XG (1:3) - HUBPOS(1:3)

!!----- Rotate hub centered coordinates for yaw and inclination
      X0  (1:3) = matmul( ATtot(1:3,1:3), X0 (1:3) )

!qq for partial wake
!!    X0(2)=X0(2)-105.d0
      call VELINWIND (X0, uturb, time, IGRDTYP)
!     call DEFINWIND (X0, dturb, time, IGRDTYP)
   endif !ITURB

   UG(1  ) = uturb(1  ) + velhub_inst
   UG(2:3) = uturb(2:3);
   Uo(1:3) = UG(1:3)
!--- Rotate Wind for inclination and wind yaw
   UG(1:3) = matmul( Atot(1:3,1:3),UG(1:3) )


!--- Partial Wake
   if (IPartWake == 1) then
!------ Control point coordinates with respect to hub (used for turbulent wind)
      X0  (1:3) = XG (1:3) - HUBPOS(1:3)

      RposPWake = dsqrt((X0(2)-LateralDist)**2+X0(3)**2)

      call LIN_INT( RposPWake, U_def,deficitR,deficitU, NumDeficitP, NumDeficitP )

      UG(1)     = UG(1) + U_def

!------ deformation
!!    fi = datan2(X0(3),X0(2)-LateralDist)
!!    call LIN_INT( RposPWake, Ur_def,deficitR,deficitUr, NumDeficitP, NumDeficitP )
!!    dturb(1,2) = Ur_def * dcos(fi)
!!    dturb(1,3) = Ur_def * dsin(fi)
   endif
!AVATAR
   if (IPartWake == 2) then
!------ Control point coordinates with respect to hub (used for turbulent wind)
      X0  (1:3) = XG (1:3) - HUBPOS(1:3)
      RAVATAR=102.88; XupWT=-10.*RAVATAR; YupWT=RAVATAR; ZupWT=0.
      RrelPW=dsqrt((X0(2)-YupWT)**2+(X0(3)-ZupWT)**2)
      UG(1)     = VELHUB_sh* &
      (0.65+tanh((RrelPW/1.08/RAVATAR)**2.5)-0.65*tanh((RrelPW/0.75/RAVATAR)**1.8))
   endif


!------ Control point coordinates with respect to the tower frame (used for tower shadow)
   X_tow (1:3) = matmul(AT_tow(1:3,1:3),XG(1:3)-RG_tow(1:3))
   U_tow (1:3) = DRG_tow(1:3)

!--- Include tower shadow effect (modify UG inside)
!  if (include_tshadow==1) &
!  call TOWER_SHADOW  (X_tow, U_tow, UG)
!  if (include_tshadow==1) &
!  call TOWER_SHADOWD (X_tow, U_tow, DG)

!--- calc Dix
   dUo(1    ) = dturb (1  ,1) + dxVELHUB_t                      &
                              + dxVELHUB_eshy* XG (2)           &
                              + dxVELHUB_eshz*(XG (3)-HUBPOS(3))
   dUo(2:3  ) = dturb (2:3,1)
    DG(1:3,1) = matmul ( Atot(1:3,1:3),dUo(1:3))              + &
                matmul (dAtot(1:3,1:3), Uo(1:3)) * dxWINDYAW_t 

!--- calc Diy
   dUo(1    ) = dturb (1  ,2) +   VELHUB_eshy
   dUo(2:3  ) = dturb (2:3,2)
    DG(1:3,2) = matmul ( Atot(1:3,1:3),dUo(1:3))


!--- calc Diz
   dUo(1    ) = dturb (1  ,3) + dzvelhub_log                    &
                              +   VELHUB_eshz
   dUo(2:3  ) = dturb (2:3,3)
    DG(1:3,3) = matmul ( Atot(1:3,1:3),dUo(1:3))              + &
                matmul (dAtot(1:3,1:3), Uo(1:3)) * VEER


 END Subroutine DINFLOW_dp
!------------------------------------------------------------------
!
! ----Subroutine TOWER_SHADOW -------------------------------------
!
!------------------------------------------------------------------
 Subroutine TOWER_SHADOW (XAB, UTOW, UVEL)
!------------------------------------------------------------------

 use mod_inflow !ISHADOW, Rbot_tow, Rtop_tow,H_tow

   implicit none

   real(8) :: XAB(3), UVEL(3), RDIST2, Tap_tow, R_tow, UTSH(2), UTOW(3)


   if ( ISHADOW == 0    ) return
   if ( XAB(3)  >  H_tow) return

   RDIST2    = (XAB(1)**2+XAB(2)**2)**2 + 1.d-15
   Tap_tow   = (Rbot_tow - Rtop_tow)/H_tow
   R_tow     =  Rbot_tow -  Tap_tow*XAB(3)
   UTSH(1)   = ( (UVEL(1)-UTOW(1))*(XAB(2)**2-XAB(1)**2) &
                -(UVEL(2)-UTOW(2))* 2.d0*XAB(1)*XAB(2)   ) /RDIST2*R_tow**2
   UTSH(2)   = (-(UVEL(2)-UTOW(2))*(XAB(2)**2-XAB(1)**2) &
                -(UVEL(1)-UTOW(1))* 2.d0*XAB(1)*XAB(2)   ) /RDIST2*R_tow**2
   UVEL(1:2) = UVEL(1:2) + UTSH(1:2)


 END Subroutine TOWER_SHADOW
!----------------------------------------------------------------------
 Subroutine INFLOW_out ( time )
!----------------------------------------------------------------------

 use mod_inflow

   implicit none

   real(8),intent(in ) :: time
   integer             :: i, j


#ifndef ASCII
#define ASCII 1
#endif

#if   ASCII == 0
   open (1,file='velhub.bin', access='append',form='UNFORMATTED')
    write (1)                         &
#elif ASCII == 1
   open (1,file='velhub.dat', access='append')
    write (1,100)                     &
#endif
     sngl(time)                      ,&
    (sngl(UG_check      (i,1)),i=1,3),&
     sngl(dsqrt(UG_check(1,1)**2 + &
                UG_check(2,1)**2 + &
                UG_check(3,1)**2    ))
   close (1)


#if   ASCII == 0
   open (1,file='velcheck.bin', access='append',form='UNFORMATTED')
    write (1)                                 &
#elif ASCII == 1
   open (1,file='velcheck.dat', access='append')
    write (1,100)                             &
#endif
      sngl(time)                             ,&
    ((sngl(UG_check   (i,j)),i=1,3),j=2,N_check)
   close (1)

 100  format (150en16.5e3)

 END Subroutine INFLOW_out
!----------------------------------------------------------------------
 Subroutine ROT_MATRIX1 ( IAX, Q0, AQ0, AQ01 )
!----------------------------------------------------------------------

   implicit none

   integer :: IAX
   real(8) :: AQ0(3,3), AQ01(3,3), CQ0, SQ0, Q0


   CQ0  = dcos(Q0)
   SQ0  = dsin(Q0)
   AQ0  =  0.d0;
   AQ01 =  0.d0;

   goto (1,2,3), IAX

 1 AQ0  (1,1) =  1.d0;
   AQ0  (2,2) =   CQ0;   AQ01 (2,2) = - SQ0;
   AQ0  (2,3) = - SQ0;   AQ01 (2,3) = - CQ0;
   AQ0  (3,2) =   SQ0;   AQ01 (3,2) =   CQ0;
   AQ0  (3,3) =   CQ0;   AQ01 (3,3) = - SQ0;
  return
 
 
 2 AQ0  (1,1) =   CQ0;   AQ01 (1,1) = - SQ0;
   AQ0  (1,3) =   SQ0;   AQ01 (1,3) =   CQ0;
   AQ0  (2,2) =  1.d0;
   AQ0  (3,1) = - SQ0;   AQ01 (3,1) = - CQ0;
   AQ0  (3,3) =   CQ0;   AQ01 (3,3) = - SQ0;
  return

 
 3 AQ0  (1,1) =   CQ0;   AQ01 (1,1) = - SQ0;
   AQ0  (1,2) = - SQ0;   AQ01 (1,2) = - CQ0;
   AQ0  (2,1) =   SQ0;   AQ01 (2,1) =   CQ0;
   AQ0  (2,2) =   CQ0;   AQ01 (2,2) = - SQ0;
   AQ0  (3,3) =  1.d0;
  return
 

 END Subroutine ROT_MATRIX1


#include "../../hGAST/source/anemos.f90"
