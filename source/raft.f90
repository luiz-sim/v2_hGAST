!----------------------------------------------------------------------
 Module Craft
!----------------------------------------------------------------------
 integer, parameter :: MBLADE      =   3
 integer, parameter :: NSTRIPM     = 300
 integer, parameter :: NSPANM      = 300
 integer, parameter :: NANGM       = 360

 integer, parameter :: AngFLAPSmax = 30   ! Maximum flap angles in Curves
 integer, parameter :: NNODEM      = 200  ! Maximum airfoil points
 integer, parameter :: NFOURM      = 20   ! Maximum Fourier coefficients
 integer, parameter :: NWAKE       = 10   ! WAKE 
 real(8), parameter :: XCLPI       = 0.25 ! XC Pitching Centre

!----------------------------------------------------------------------
 type :: type_node       ! <NBLADE, NSTRIP+1>
!----------------------------------------------------------------------
    real(8) :: xBld (3)
    real(8) :: xRD  (3)
 END type type_node
!----------------------------------------------------------------------
 type :: type_strip      ! <NBLADE, NSTRIP>
!----------------------------------------------------------------------
    real(8) :: rBld, drBld
    real(8) :: rRD , drRD
    real(8) :: chord
    real(8) :: xae, zae
    real(8) :: fa
    real(8) :: floss
    real(8) :: AINDTB    , AINDTB_PRE
    real(8) :: AINDPTB   , AINDPTB_PRE
    real(8) :: AINDTBmean, AINDTB_PREmean
    real(8) :: ALPHATB
    real(8) :: PHITB
    real(8) :: twist     !old TWISTTB
    real(8) :: pitch     !old THY0TB
    real(8) :: THYTB
    real(8) :: THY1TB
    real(8) :: THY2TB
    real(8) :: WEFFXTB, WEFFZTB !WEFFTB(3)
    real(8) :: Mach
    real(8) :: Reynolds
    real(8) :: CLIFTTB
    real(8) :: CDRAGTB
    real(8) :: CNORMTB
    real(8) :: CTANGTB
    real(8) :: CMOM4TB
    real(8) :: xi_sk
    real(8) :: gama_sk
    real(8) :: psi0_sk
    real(8) :: Coef_sk
!----- loading
    real(8) :: FCP  (6)
    real(8) :: FCPL (6)
!---- velocities
    real(8) :: UBODY   (3) !RD     c.s.
    real(8) :: UWINDTB (3) !RD     c.s.
    real(8) :: UWINDGTB(3) !global c.s.
    real(8) :: UWAKIZTB
    real(8) :: WEFF_RD (3) !RD     c.s.
    real(8) :: Ui_ax 
    real(8) :: Ui_per
    real(8) :: A_curv(3,3), A(3,3)
!----- Beam correspondance
    real(8) :: HTAe2a
    integer :: NSUBe2a
    integer :: NELe2a
 END type type_strip
!--------------------------------------------------------------------------
 type (type_strip), allocatable, save :: strip    (:,:  )                    !NBLADE_ae, NSTRIP
 type (type_node ), allocatable, save :: node     (:,:  )                    !NBLADE_ae, NSTRIP+1
 real(8)          , allocatable, save :: AP_ae    (:,:,:), ATP_ae    (:,:,:) !(NBLADE_ae,3,3)
 real(8)          , allocatable, save :: A_pit_ae (:,:,:,:)                  !(NBLADE_ae,NSTRIPM,3,3)
 real(8)          , allocatable, save :: PHI0_ae  (:    )                    !(NBLADE_ae    )
 integer          ,              save :: NBLADE_ae
 real(8)          ,              save :: DT_ae, TIME_ae, AZIM_ae, OMEGA_ae, CONE_ae, Hhub_ae
 real(8)          ,              save :: TILT_ae, YAW_ae
 real(8)          ,              save :: PI_ae, PI2_ae , R2D_ae
 real(8)          ,              save :: A_cone_ae(3,3), AT_cone_ae(3,3)
 real(8)          ,              save :: HUBPOS_ae(  3)
!--------------------------------------------------------------------------
!- main parameters
 integer, save      :: NSTRIP
 integer, save      :: INOINDUCT, IDYNWAKE, IDYNSTALL, IDYNSTRIP
 integer, save      :: ITIPLOS  , IHUBLOS
 integer, save      :: IPARKED
 real(8), save      :: HUB_VEL (3)
 real(8), save      :: VELHUBT  , WINDYAW , WINDINC  , SHEXP
 real(8), save      :: AIRDEN   , SSPEED
 real(8), save      :: RHUB     , RTIP
 real(8), save      :: RHCONE   , RTCONE
!- read profilb, airfoils [lookup tables]
 integer, save      :: IspanANGmax
 integer, save      :: NSPANB2   (MBLADE                         )
 real(8), save      :: RCB2      (MBLADE,NSPANM                  )
 integer, save      :: NANGB     (MBLADE,NSPANM                  )
 real(8), save      :: AATB      (MBLADE,NSPANM,            NANGM)
 real(8), save      :: CLTB      (MBLADE,NSPANM,AngFLAPSmax,NANGM)
 real(8), save      :: CDTB      (MBLADE,NSPANM,AngFLAPSmax,NANGM)
 real(8), save      :: CMTB      (MBLADE,NSPANM,AngFLAPSmax,NANGM)
!- OUTPUT [Integrated loads/power]
 real(8), save      :: TTHRUST, TTORQUE, POWER
 real(8), save      :: TORQUE    (MBLADE)
 real(8), save      :: THRUST    (MBLADE)
!- Aero2Elast Vectors
 real(8), save      :: FSTRIP_el (6,MBLADE,NSTRIPM)
 real(8), save      :: XCAER_el  (3,MBLADE,NSTRIPM)
 real(8), save      :: RSTRIP_el (  MBLADE,NSTRIPM)
!-----------------------------------------------------------------------------------------------
!/FLAPS/
 integer, save              :: INDXFlap, Nflaps
 integer, save              :: IFLAPLOC   (MBLADE,NSTRIPM            )
 integer, save              :: NFLAPANGS  (MBLADE,NSPANM             )
 real(8), save              :: FLAPANGS   (MBLADE,NSPANM ,AngFLAPSmax)
 real(8), save              :: FLCONTR    (MBLADE,NSTRIPM            )
 real(8), save, allocatable :: Rflaps (:,:)
 integer, save              :: IFLAP_typ    !jim - Flap type 1:collective, 2:cyclic, 10:ctrl
 real(8), save              :: FLAP_ampl    !jim - Amplitude of motion [deg]
 real(8), save              :: NTIME_INI

 END Module Craft
!----------------------------------------------------------------------
 Module Caero
!----------------------------------------------------------------------

 Use Craft, ONLY : MBLADE,NSTRIPM,NANGM

!- AEROPAR
 real(8), save :: AZER_ae
 real(8), save :: AZER0_ae
 real(8), save :: DCLLIN_ae, DCMLIN_ae
 real(8), save :: CDLIN_ae,  CMLIN_ae
 real(8), save :: CMO_ae
 real(8), save :: DCL_ae , DCD_ae , DCM_ae
 real(8), save :: DACL_ae, DACD_ae, DACM_ae
 real(8), save :: CL_ae  , CD_ae  , CM_ae
 real(8), save :: AS_L_ae
 real(8), save :: AK_L_ae
 real(8), save :: AL_L_ae
 real(8), save :: AAL_L_ae
 real(8), save :: ASIG_L_ae
 real(8), save :: AD_L_ae
 real(8), save :: AA0_L_ae, AR0_L_ae
 real(8), save :: AA2_L_ae, AR2_L_ae, AE2_L_ae
 real(8), save :: ASIG_D_ae
 real(8), save :: AA0_D_ae, AR0_D_ae
 real(8), save :: AA2_D_ae, AR2_D_ae, AE2_D_ae
 real(8), save :: AS_M_ae
 real(8), save :: ASIGP_M_ae
 real(8), save :: ASIG_M_ae
 real(8), save :: AD_M_ae
 real(8), save :: AA0_M_ae, AR0_M_ae
 real(8), save :: AA2_M_ae, AR2_M_ae, AE2_M_ae

!- ONERA VARIABLES
 real(8), save :: CLOLDA  (3,MBLADE*NSTRIPM)
 real(8), save :: CLOLDB  (3,MBLADE*NSTRIPM)
 real(8), save :: CDOLDB  (3,MBLADE*NSTRIPM)
 real(8), save :: CMOLDB  (3,MBLADE*NSTRIPM)
 real(8), save :: CLIFTA  (3,MBLADE*NSTRIPM)
 real(8), save :: CLIFTB  (3,MBLADE*NSTRIPM)
 real(8), save :: CDRAGB  (3,MBLADE*NSTRIPM)
 real(8), save :: CMOMB   (3,MBLADE*NSTRIPM)
 real(8), save :: VELEFF  (  MBLADE*NSTRIPM)
 real(8), save :: VELEFFP (3,MBLADE*NSTRIPM)
 real(8), save :: ANGEFF  (  MBLADE*NSTRIPM)
 real(8), save :: ANGEFFP (3,MBLADE*NSTRIPM)                   !!DW0(MBLADE*NSTRIPM)
!- BEDDOES-LEISHMAN VARS
 real(8), save :: BX1     (3,MBLADE*NSTRIPM)
 real(8), save :: BX2     (3,MBLADE*NSTRIPM)
 real(8), save :: BX3     (3,MBLADE*NSTRIPM)
 real(8), save :: BX4     (3,MBLADE*NSTRIPM)
 real(8), save :: BXOLD1  (3,MBLADE*NSTRIPM)
 real(8), save :: BXOLD2  (3,MBLADE*NSTRIPM)
 real(8), save :: BXOLD3  (3,MBLADE*NSTRIPM)
 real(8), save :: BXOLD4  (3,MBLADE*NSTRIPM)

 integer, save :: NTIME_un

!- POLAR
 real(8), save :: AGON(NANGM)
 real(8), save :: CLST(NANGM)
 real(8), save :: CDST(NANGM)
 real(8), save :: CMST(NANGM)

 integer, save :: NPP

 END Module Caero
!----------------------------------------------------------------------
 Module Foilfs_mod
!----------------------------------------------------------------------

 Use Craft, ONLY : MBLADE,NSTRIPM,NSPANM,NNODEM,NFOURM,NWAKE,XCLPI,AngFLAPSmax
 Use Craft, ONLY : AIRDEN
 Use Cbeam, ONLY : DT_el, PI, NBLADE_el, TIME

!/FOILFS/
 integer, save              :: NFOUR
 integer, save              :: NNODE 
 real(8), save              :: XNODE     ( 2,NNODEM)
 real(8), save              :: XNODEL    ( 2,NNODEM)
 real(8), save              :: XCP       (   NNODEM)
 real(8), save              :: YNODEOLD  (   NNODEM) 
 real(8), save              :: YNODEOLDN (MBLADE, NSTRIPM,NNODEM)
 real(8), save              :: YTEMP     (MBLADE, NSTRIPM,NNODEM)
 real(8), save              :: XTEMP     (MBLADE, NSTRIPM,NNODEM)
 real(8), save              :: XC4N      (MBLADE, NSTRIPM,NNODEM)
 real(8), save              :: YCP       (   NNODEM)
 real(8), save              :: XCPL      (   NNODEM)
 real(8), save              :: YCPL      (   NNODEM)
 real(8), save              :: DYXCP     (   NNODEM)
 real(8), save              :: DYT       (   NNODEM)
 real(8), save              :: DX        (   NNODEM)
 real(8), save              :: WWASH     (   NNODEM)
 real(8), save              :: WWASH_1   (   NNODEM)
 real(8), save              :: AFOUR     (   NFOURM)
 real(8), save              :: XWAKE       (MBLADE, NSTRIPM, NWAKE)
 real(8), save              :: YWAKE       (MBLADE, NSTRIPM, NWAKE)
 real(8), save              :: DXW         (MBLADE, NSTRIPM, NWAKE)
 real(8), save              :: ALPHAFOIL   (MBLADE, NSTRIPM)
 real(8), save              :: AZER0pot_ae (MBLADE, NSTRIPM)
 real(8), save              :: AZER_ae0    (MBLADE, NSTRIPM) 
 real(8), save, allocatable :: GWAKE    (:,:,:)
 real(8), save, allocatable :: circ     (:,:,:)
 real(8), save, allocatable :: circ_int (:,:,:)
 real(8), save, allocatable :: circ_old (:,:,:)
 real(8), save, allocatable :: cmpot    (:,:  )
 real(8), save, allocatable :: clf      (:,:,:)
 real(8), save, allocatable :: cdf      (:,:,:)
 real(8), save, allocatable :: cmf      (:,:,:)
!/GEOMETRY/
 real(8), save, allocatable :: RRLOC    (:,:  )
 real(8), save, allocatable :: XAIRNODE (:,:,:)
 real(8), save, allocatable :: YAIRNODE (:,:,:)
 real(8), save, allocatable :: XNSTR    (:,:,:)
 real(8), save, allocatable :: YNSTR    (:,:,:)
 real(8), save, allocatable :: NSTRIPL1 (:,:  )
 real(8), save, allocatable :: NSTRIPL2 (:,:  )
 real(8), save, allocatable :: MH       (:    )
 integer, save              :: NPOINTS  ( 5,30)

!/FLAPS/
 real(8), save, allocatable :: FLAPEX   (:)   !Flap ext.
 real(8), save, allocatable :: FLAPMODE (:)   !Flap modeshape
 real(8), save, allocatable :: XFLAP    (:)   !Flap ext. in strips

 END Module Foilfs_mod


!----------------------------------------------------------------------
 Subroutine RAFT_init
!----------------------------------------------------------------------

 Use Cbeam   !R2D,NBLADE_el, Hhub, modelpath,windfile ... TIME,PHI0,OMEGA
 Use Craft
 Use Paths
 Use Caero
 Use Foilfs_mod

   implicit none

   real(8)            :: ang, RR1, RR2, RRT, RR, ARA, PAR, DRREF, RBLADE, TIME_turb
   integer            :: nbod, NSTRIP1,NSTRIP2, i, istrip, ispan, iang, NANGmax

   integer            :: IWINDC   , ISHADOW , IGRDTYP
   real(8)            :: VELHUB   ,                    VEER
   real(8)            :: TIME_GUST, TIREF_WT, VREF_WT, RTIP_in  , HUBPOS(3)
   real(8)            :: Rbot_tow , Rtop_tow, H_tow

   character*256      :: PROFILTB (MBLADE       )
   character*256      :: GEOMFIL  (MBLADE       )
   character*256      :: CLCDCMB  (MBLADE,NSPANM)

   real(8)            :: RCB(NSPANM), VAR(NSPANM)
   integer, parameter :: nnp=301
   real(8)            :: xBld(nnp,3), Xnel(nnp), Ynel(nnp), Znel(nnp)
   integer            :: nbsub, nb_el, nel, n, n1, n2 , ip  !nbod,i
   real(8)            :: r, chd, twi, xae, zae, fa

   real(8)            :: DRC0   (NSTRIPM       )
!--- read geomb.inp
   integer            :: NSPANB1(MBLADE        )
   real(8)            :: RCB1   (MBLADE, NSPANM)
   real(8)            :: CHORDTB(MBLADE, NSPANM)
   real(8)            :: TWTB   (MBLADE, NSPANM)
   real(8)            :: XAERTB (MBLADE, NSPANM)
   real(8)            :: ZAERTB (MBLADE, NSPANM)
!--- transformation matrices
   real(8)            :: Ex(3),Ey(3),Ez(3),met, A(3,3)
!--- flaps
   real(8)            :: RLD, X1, X2, Y1, Y2  !,RLOC
   integer            :: ib, ifl, fts, lstrip !, iblade
   integer            :: L, K, ist, size, M, F, x, ga, gb, D, G !, point_num 
   character*256      :: AIRGEOM  (MBLADE,NSPANM)


   PI_ae   = PI; PI2_ae = PI2; R2D_ae = R2D; CONE_ae = CONE; DT_ae = DT_el; NBLADE_ae = NBLADE_el
   Hhub_ae = Hhub;

   Allocate ( AP_ae    (NBLADE_ae,        3,3), ATP_ae    (NBLADE_ae,3,3) );
   Allocate ( A_pit_ae (NBLADE_ae,NSTRIPM,3,3), PHI0_ae   (NBLADE_ae    ) ); PHI0_ae = PHI0;


   open (12, file=trim(file_rotorin)) !'rotorin.inp'

!---- read wind parameters
    read (12,*)                    !title
    read (12,*) VELHUB             !UINFLOW
    read (12,*) WINDYAW            !wind yaw
    read (12,*) WINDINC            !wind inclination (tilt)
    read (12,*) AIRDEN             !density
    read (12,*) SSPEED             !Speed of sound 340 m/s [define Mach]
    read (12,*) SHEXP              !shear effect exponental = 0:no shear
    read (12,*) VEER               !wind veer:: slope of yaw angle / m wrt hub height [deg/m]
    read (12,*) IWINDC             !wind case 0:uniform wind, 1:defined wind scenario, 2:EOG, 3x:ECD, 4x:EDC, 5x:EWS, 6:NTM, 7:ETM, 8:EWM, 9:read sigma
    read (12,*) VREF_WT            !reference velocity [m/s]   I:50.00, II:42.50, III:37.50
    read (12,*) TIREF_WT           !reference TI       [ - ]   A: 0.16, B : 0.14, C  : 0.12
    read (12,*) TIME_GUST          !Time to start extreme Gust [sec]
    read (12,*) IGRDTYP            !0:no turb, 1:disk, 2:rectangular
    read (12,*) TIME_turb          !time after which the turbulent wind starts
    read (12,*)

!---- read general parameters
    read (12,*) IPARKED            ![0:N, 1:Y]
    read (12,*) INOINDUCT          ![0:N, 1:Y]
    read (12,*) IDYNWAKE           ![0:N, 1:Y]
    read (12,*) IDYNSTALL          ![0:N, 1:ONERA, 2:Boddes-Leishman]
    read (12,*) IDYNSTRIP          ![after which strip to apply DS]
    read (12,*)

!---- read blade geometry parameters
    read (12,*) RHUB               !Rhub
    read (12,*) RTIP               !Rtip
    read (12,*) NSTRIP             !number of strips

!---- read correction parameters
    read (12,*) ITIPLOS            !calculate tip loss          when ITIPLOS = 1
    read (12,*) IHUBLOS            !calculate hub loss          when IHUBLOS = 1
    read (12,*) ISHADOW            !include tower shadow effect when ISHADOW = 1
    read (12,*) Rbot_tow           !R tower bottom
    read (12,*) Rtop_tow           !R tower top
    read (12,*,end=1,err=1)
    read (12,*,end=1,err=1) INDXFlap !0: NO IPC - NO IFC, 1: IFC+IPC, 2:IPC
    goto 2
 1  write(10,*) 'No flap is defined in BEM'
    INDXFlap = 0

 2  if     (INDXFlap == 1) then
       read(12,*) Nflaps           !Number of flaps
       read(12,*) NFOUR            !Number of Fouries coeff
       read(12,*) NTIME_INI        !Time to start flap action
       read(12,*)
       read(12,*) 
       Allocate (Rflaps   (Nflaps,2) )
       Allocate (FLAPEX   (Nflaps)   )
       Allocate (FLAPMODE (Nflaps)   )
       Allocate (XFLAP    (NSTRIP)   )
      do L = 1, Nflaps
       read(12,*) Rflaps(L,1), Rflaps(L,2), FLAPEX(L), FLAPMODE(L) !Rstart, Rend, flap extent
      enddo
       read(12,*) IFLAP_typ        !Flap type 1:collective, 2:cyclic, 20:IFCS tuner, 21: IFCS, 22:IFCS+IPC, 30: IFC+IPC (PID)
       read(12,*) FLAP_ampl        !Amplitude of motion [deg]
    elseif (INDXFlap == 2) then
       IFLAP_typ = 22 !for IPC
       read(12,*)                  !Number of flaps
       read(12,*)                  !Number of Fouries coeff
       read(12,*) NTIME_INI        !Time to start flap action
    endif !INDXFlap
   close (12)
    FLCONTR (:,:) = 0.d0; !zero flap angle

   WINDYAW  = WINDYAW / R2D_ae
   WINDINC  = WINDINC / R2D_ae
   VEER     = VEER    / R2D_ae

!--- no Dynamic Stall for WindYaw greater than 60deg.
   ang = 60.1d0
   if  ((dabs(WINDYAW)>ang/R2D_ae).and.(dabs(WINDYAW)<(360.d0-ang)/R2D_ae)) IDYNSTALL = 0
!--- Parked for EWM
   if  (IWINDC== 8)                    IPARKED=1


   do nbod           = 1, NBLADE_ae
      GEOMFIL (nbod) = file_geomb
      PROFILTB(nbod) = file_profilb
   enddo

   NTIME_un = 5
   nbod     = 1
   nbsub    = body(nbod)%NBODSUB_el
   nb_el    = body(nbod)%NBODTGLB_el(nbsub)
   nel      = subbody(nb_el)%NTEB_el
   n2       = subbody(nb_el)%INODNEL_el(nel,2)

   if (dabs(RTIP-beam_timo(n2)%XL_el(2)-Hhub_ae)>1.d-3) then
      write(*,*)'check RTIP',RTIP,beam_timo(n2)%XL_el(2)+Hhub_ae
      stop
   endif

   RTIP         = beam_timo(n2)%XL_el(2)+Hhub_ae
   RTCONE       = RTIP * dcos(CONE_ae)
   RHCONE       = RHUB * dcos(CONE_ae) !or Hhub_ae * dcos(CONE_ae)
   RTIP_in      = RTCONE
   HUBPOS (1:3) = HUBpos_el (1:3)
   H_tow        = Htow

!--- Initialize the inflow part
   call INFLOW_init ( IWINDC   , VELHUB  , WINDYAW, WINDINC, VEER     , SHEXP    ,&
                      TIME_GUST, TIREF_WT, VREF_WT, RTIP_in, HUBPOS   , ISHADOW  ,&
                      Rbot_tow , Rtop_tow, H_tow  , IGRDTYP, TIME_turb, file_wind   )

   write(10,*)
   write(10,*)'--RAFT INIT---'
   write(10,*)
   write(10,*)'AIRDEN    ',AIRDEN
   write(10,*)'SSPEED    ',SSPEED
   write(10,*)'IPARKED   ',IPARKED
   write(10,*)'INOINDUCT ',INOINDUCT
   write(10,*)'IDYNWAKE  ',IDYNWAKE
   write(10,*)'IDYNSTALL ',IDYNSTALL
   write(10,*)'IDYNSTRIP ',IDYNSTRIP
   write(10,*)'RHUB      ',RHUB
   write(10,*)'RTIP      ',RTIP
   write(10,*)'NSTRIP    ',NSTRIP
   write(10,*)'ITIPLOS   ',ITIPLOS
   write(10,*)'IHUBLOS   ',IHUBLOS
   write(10,*)
   write(10,*)'INDXFlap  ',INDXFlap 

  if (INDXFlap > 0) then
   write(10,*)'Nflaps    ',Nflaps   
   write(10,*)'NFOUR     ',NFOUR    
   write(10,*)'NTIME_INI ',NTIME_INI
   write(10,*)'IFLAP_typ ',IFLAP_typ 
   write(10,*)'FLAP_ampl ',FLAP_ampl 
  endif


   Allocate (strip(NBLADE_ae,NSTRIP  ))
   Allocate (node (NBLADE_ae,NSTRIP+1))

!--- Define Blade geometrical parameters
   RR1              = RHUB - Hhub_ae
   RR2              = RTIP - RHUB
   RRT              = RTIP - Hhub_ae
   NSTRIP1          = (RR1/RRT)*NSTRIP
   NSTRIP1          = max(NSTRIP1,1)
   NSTRIP2          = NSTRIP - NSTRIP1

  if (INDXFlap/=1) then
   do istrip        = 1, NSTRIP1
      DRC0(istrip)  =  RR1/RTIP/NSTRIP1
   enddo

   ARA = 0.95d0
   PAR = 1.0d0
   do i=1,NSTRIP2-1
      PAR = PAR + ARA**i
   enddo

   DRREF            = RR2/PAR

   do istrip = 1, NSTRIP2
      DRC0(NSTRIP1+istrip)  =  DRREF*ARA**(istrip-1)/RTIP
!!    DR0 (NSTRIP1+istrip)  =  DRC0(NSTRIP1+istrip) * dcos(CONE_ae)
   enddo

  else!if (INDXFlap==1) then
      D=1
  if  ((Rflaps(1,1)==Rhub).and.(Rflaps(Nflaps,2)==Rtip)) then
      size           = 2*Nflaps-1
      allocate(RRLOC(size,4))
      ist            = 1
   do K              = 1, Nflaps
      RRLOC(ist  ,1) = Rflaps(K,1)
      RRLOC(ist  ,2) = Rflaps(K,2)
      RRLOC(ist  ,3) = RRLOC(ist,2)-RRLOC(ist,1)
      RRLOC(ist  ,4) = 1
     if (K/=Nflaps) then
      RRLOC(ist+1,1) = Rflaps(K,2)
      RRLOC(ist+1,2) = Rflaps(K+1,1)
      RRLOC(ist+1,3) = RRLOC(ist+1,2)-RRLOC(ist+1,1)
      RRLOC(ist+1,4) = 0
     endif
      ist            = ist+2
   enddo !K
  elseif (Rflaps(1,1)==Rhub) then
      size           = 2*Nflaps
      allocate(RRLOC(size,4))
      ist            = 1
   do K              = 1, Nflaps
      RRLOC(ist  ,1) = Rflaps(K,1)
      RRLOC(ist  ,2) = Rflaps(K,2)
      RRLOC(ist  ,3) = RRLOC(ist,2)-RRLOC(ist,1)
      RRLOC(ist  ,4) = 1
      RRLOC(ist+1,1) = Rflaps(K,2)
     if (K==NFLAPS) then
      RRLOC(ist+1,2) = Rtip
     else
      RRLOC(ist+1,2) = Rflaps(K+1,1)
     endif
      RRLOC(ist+1,3) = RRLOC(ist+1,2)-RRLOC(ist+1,1)
      RRLOC(ist+1,4) = 0
      ist=ist+2
   enddo !K
  elseif (Rflaps(Nflaps,2)==Rtip) then
      size           = 2*Nflaps
      allocate(RRLOC(size,4))
      ist            = 2*Nflaps
   do K              = Nflaps,1,-1
      RRLOC(ist  ,2) = Rflaps(K,2)
      RRLOC(ist  ,1) = Rflaps(K,1)
      RRLOC(ist  ,3) = RRLOC(ist,2)-RRLOC(ist,1)
      RRLOC(ist  ,4) = 1
      RRLOC(ist-1,2) = Rflaps(K,1)
     if (K==1) then
      RRLOC(ist-1,1) = Rhub
     else
      RRLOC(ist-1,1) = Rflaps(K-1,2)
     endif
      RRLOC(ist-1,3) = RRLOC(ist-1,2)-RRLOC(ist-1,1)
      RRLOC(ist+1,4) = 0
      ist            = ist-2
   enddo
  else
      size=2*Nflaps+1
      allocate(RRLOC(size,4))
      ist=1
   do K              = 1, Nflaps+1      
     if (K==1) then
      RRLOC(ist  ,1) = Rhub
      RRLOC(ist  ,2) = Rflaps(K,1)
      RRLOC(ist  ,3) = RRLOC(ist,2)-RRLOC(ist,1)
      RRLOC(ist  ,4) = 0
      RRLOC(ist+1,1) = Rflaps(K,1)
      RRLOC(ist+1,2) = Rflaps(K,2)
      RRLOC(ist+1,3) = RRLOC(ist+1,2)-RRLOC(ist+1,1)
      RRLOC(ist+1,4) = 1
     elseif (K==Nflaps+1) then
      RRLOC(ist  ,2) = Rtip
      RRLOC(ist  ,1) = Rflaps(K-1,2)
      RRLOC(ist  ,3) = RRLOC(ist,2)-RRLOC(ist,1)
      RRLOC(ist  ,4) = 0
     else
      RRLOC(ist  ,1) = Rflaps(K-1,2)
      RRLOC(ist  ,2) = Rflaps(K,1)
      RRLOC(ist  ,3) = RRLOC(ist,2)-RRLOC(ist,1)
      RRLOC(ist  ,4) = 0
      RRLOC(ist+1,1) = Rflaps(K,1)
      RRLOC(ist+1,2) = Rflaps(K,2)
      RRLOC(ist+1,3) = RRLOC(ist+1,2)-RRLOC(ist+1,1)
      RRLOC(ist+1,4) = 1 
     endif
      ist            = ist+2
   enddo !K
  endif

!-- Strip creation
   ga = max(nint(RR1/RRT*Nstrip),1)  ! Strips in Hhub-Rhub region
   gb = Nstrip-ga                    ! Strips in flap region
   allocate(NSTRIPL1(ga,5))
   allocate(NSTRIPL2(gb,5))
      
   do    F             = 1, size-1
         x             = max(nint((RRLOC(F,3)/(RRT))*gb),1)
      do M             = 1, x
         NSTRIPL2(D,1) = RRLOC(F,1)+(M-1)*(RRLOC(F,3)/x)
         NSTRIPL2(D,2) = RRLOC(F,1)+M*(RRLOC(F,3)/x)
         NSTRIPL2(D,3) = NSTRIPL2(D,2)-NSTRIPL2(D,1)
         NSTRIPL2(D,4) = NSTRIPL2(D,3)/2+NSTRIPL2(D,1)
        if (RRLOC(F,4)==1) then
         NSTRIPL2(D,5) = 1
        else
         NSTRIPL2(D,5) = 0
        endif
         D             = D+1
      enddo !M
   enddo !F
      
      G             = D-1
   do M             = 1, (gb-G)
      NSTRIPL2(D,1) = RRLOC(size,1)+(M-1)*RRLOC(size,3)/(gb-G)
      NSTRIPL2(D,2) = RRLOC(size,1)+M*RRLOC(size,3)/(gb-G)
      NSTRIPL2(D,3) = NSTRIPL2(D,2)-NSTRIPL2(D,1)
      NSTRIPL2(D,4) = NSTRIPL2(D,3)/2+NSTRIPL2(D,1)
     if (RRLOC(size,4)==1) then
      NSTRIPL2(D,5) = 1
     else
      NSTRIPL2(D,5) = 0
     endif
      D             = D+1
   enddo !M
      
   do M             = 1, ga
      NSTRIPL1(M,1) = Hhub+(M-1)*RR1/ga
      NSTRIPL1(M,2) = Hhub+M*RR1/ga
      NSTRIPL1(M,3) = NSTRIPL1(M,2)-NSTRIPL1(M,1)
      NSTRIPL1(M,4) = NSTRIPL1(M,3)/2+NSTRIPL1(M,1)
      NSTRIPL1(M,5) = 0
   enddo
     
   do    ib                 = 1, NBLADE_el
      do M                  = 1, ga
         DRC0    (       M) = NSTRIPL1(M,3)/RTIP
         IFLAPLOC(ib,    M) = NSTRIPL1(M,5)
      enddo !M
      
      do M                  = 1, gb
         DRC0    (    ga+M) = NSTRIPL2(M,3)/RTIP
         IFLAPLOC(ib, ga+M) = NSTRIPL2(M,5)
      enddo !M
   enddo !ib     
  endif !INDXFlap
!qq CHECK
!! INDXFlap=0;  IFLAPLOC(:, :) =0;
   RR               = (Hhub_ae/RTIP)

!! open (38,file=trim(modelpath_inp)//'aero.inp')

   write(10,*)
   write(10,*)'**Aero strips**'
   write(10,*)'NSTRIP ',NSTRIP
   write(10,*)'NSTRIP1',NSTRIP1
   write(10,*)'NSTRIP2',NSTRIP2
   write(10,*)'RHUB   ',RHUB
   write(10,*)'RTIP   ',RTIP

   do istrip            = 1, NSTRIP
       RBLADE           = RR + DRC0(istrip)/2d0
!jim new
!!     read (38,*,end=33) RBLADE,DRC0(istrip)
!!     RBLADE           = RBLADE/RTIP
!!     DRC0(istrip)     = DRC0(istrip)/RTIP
!jim new
       RR = RBLADE + DRC0(istrip)/2d0

       write(10,'(i4,4f16.3,i3)')      &
        istrip                        ,&
        RBLADE                  *RTIP ,&
       (RBLADE-DRC0(istrip)/2d0)*RTIP ,&
       (RBLADE+DRC0(istrip)/2d0)*RTIP ,&
               DRC0(istrip)     *RTIP ,&
        IFLAPLOC(1, istrip) 

!      write(10,*) istrip,DRC0(istrip)*RTIP,RBLADE*RTIP,IFLAPLOC(1, istrip) 

       xBld (istrip  ,2) = (RBLADE-DRC0(istrip)/2d0)*RTIP !blade the cone has not been considered yet
   enddo
       xBld (NSTRIP+1,2) = (RBLADE+DRC0(NSTRIP)/2d0)*RTIP

!! close(38)
   write(10,*)

!!   goto 34
!!33 write(*,*)'error read aero.inp'
!!   stop
34 continue


!--- Set bend-sweep coordinates
   call ROT_MATRIX0 ( 1, CONE_ae, A_cone_ae )
   AT_cone_ae  = transpose (A_cone_ae)

   do    nbod  = 1, NBLADE_ae
         nbsub = 1
         nb_el = body(nbod)%NBODTGLB_el(nbsub)
         nel   = 1
         n1    = subbody(nb_el)%INODNEL_el(nel,1)
         nbsub = body(nbod)%NBODSUB_el
         nb_el = body(nbod)%NBODTGLB_el(nbsub)
         nel   = subbody(nb_el)%NTEB_el
         n2    = subbody(nb_el)%INODNEL_el(nel,2)
         i     = 0
      do n     = n1, n2
         if (mod(n-n1,2)/=0) cycle
         i     = i+1
         Xnel (i) = beam_timo(n)%XL_el(1)
         Ynel (i) = beam_timo(n)%XL_el(2)+Hhub_ae
         Znel (i) = beam_timo(n)%XL_el(3)
!w       write(200,'(5f15.5,2i3)') Ynel(i),Xnel(i),Znel(i),n ,i
      enddo !i
         i     = i+1
         Xnel (i) = beam_timo(n2)%XL_el(1)
         Ynel (i) = beam_timo(n2)%XL_el(2)+Hhub_ae
         Znel (i) = beam_timo(n2)%XL_el(3)
!w       write(200,'(3f15.5,2i3)') Ynel(i),Xnel(i),Znel(i),n2,i
         n  = i
      do i  = 1, NSTRIP+1

         call lint1   ( Ynel, Xnel, 1, n, nnp, xBld(i,2), xBld(i,1) )
         call lint1   ( Ynel, Znel, 1, n, nnp, xBld(i,2), xBld(i,3) )

         node(nbod,i)%xBld(1:3) =                           xBld(i,1:3)
         node(nbod,i)%xRD (1:3) = matmul(A_cone_ae(1:3,1:3),xBld(i,1:3))
!w       write(201,'(3f15.5,i3)') xBld(i,2),xBld(i,1),xBld(i,3),i
      enddo !in
!w       write(200,*); write(200,*)
!w       write(201,*); write(201,*)
      do i          = 1, NSTRIP
         ip         = i+1
         strip(nbod,i)%rRD   = (node(nbod,ip)%xRD (2)+node(nbod,i)%xRD (2)) / 2.d0
         strip(nbod,i)%rBld  = (node(nbod,ip)%xBld(2)+node(nbod,i)%xBld(2)) / 2.d0
         strip(nbod,i)%drRD  =  node(nbod,ip)%xRD (2)-node(nbod,i)%xRD(2)
         strip(nbod,i)%drBld = dsqrt( (node(nbod,ip)%xRD(1)-node(nbod,i)%xRD(1))**2 +&
                                      (node(nbod,ip)%xRD(2)-node(nbod,i)%xRD(2))**2 +&
                                      (node(nbod,ip)%xRD(3)-node(nbod,i)%xRD(3))**2   )
      enddo !i
   enddo !nb


!--- Set transformation matrices
   do    nbod     = 1, NBLADE_ae
      do i        = 1, NSTRIP
         ip       = i+1
!--------- Set Ez=i x Ey_bend
         Ex(1)    = 1.d0;  Ex(2:3) = 0.d0;
         Ey(1)    = 0.d0;  Ey(2:3) = node(nbod,ip)%xBld(2:3) - node(nbod,i)%xBld(2:3);
         met      = dsqrt(dot_product (Ey,Ey))
         Ey(1:3)  = Ey(1:3)/met
         call EXTEPR ( Ex, Ey, Ez )
         met      = dsqrt(dot_product (Ez,Ez))
         Ez(1:3)  = Ez(1:3)/met

!--------- Set Ey
         Ey(1:3)  = node(nbod,ip)%xBld(1:3) - node(nbod,i)%xBld(1:3);
         met      = dsqrt(dot_product (Ey,Ey))
         Ey(1:3)  = Ey(1:3)/met

!--------- Set Ex=Ey x Ez
         call EXTEPR ( Ey, Ez, Ex )
         met      = dsqrt(dot_product (Ex,Ex))
         Ex(1:3)  = Ex(1:3)/met

         A(1:3,1) = Ex(1:3)
         A(1:3,2) = Ey(1:3)
         A(1:3,3) = Ez(1:3)

         strip(nbod,i)%A_curv(1:3,1:3) = A(1:3,1:3)
      enddo !i
   enddo !nb


!--- Read Input Data for the BLADES
   if (INDXFlap==1) then
      allocate(XAIRNODE(NBLADE_el,Nstrip,NNODEM))  ! Allocate matrices of airfoil nodes
      allocate(YAIRNODE(NBLADE_el,Nstrip,NNODEM))  ! Maximum nodes allowed in airfoil geom=100
      allocate(XNSTR   (NBLADE_el,Nstrip,NNODEM))
      allocate(YNSTR   (NBLADE_el,Nstrip,NNODEM))
      allocate(MH(NNODEM))
   endif !INDXFlap

   NANGmax     = 0

   do nbod =1, NBLADE_ae

      open (13,file=trim(GEOMFIL(nbod)) ) !'geomb.inp'

      read (13,*)
      read (13,*) NSPANB1(nbod)
      read (13,*)

      do ispan =1, NSPANB1(nbod)

         read (13,*) RCB1    (nbod,ispan),&   !R
                     CHORDTB (nbod,ispan),&   !chord
                     TWTB    (nbod,ispan),&   !twist
                     XAERTB  (nbod,ispan),&   !distance between elastic
                     ZAERTB  (nbod,ispan)     !and aerodynamic axis
      enddo !ispan

      close(13)

     open (14,file=trim(PROFILTB(nbod)) ) !'profilb.inp'

      read (14,*)
      read (14,*) NSPANB2(nbod)
      read (14,*)

      do ispan = 1, NSPANB2(nbod)

         if (INDXFlap==1) then; read (14,*) RCB2(nbod,ispan), CLCDCMB (nbod,ispan), AIRGEOM (nbod,ispan)
         else                 ; read (14,*) RCB2(nbod,ispan), CLCDCMB (nbod,ispan)
         endif

        open  (15 , file=trim(dir_airfoils)//trim(CLCDCMB(nbod,ispan)) )        !airfoil polars

         read (15 ,*) NFLAPANGS  (nbod,ispan)                                   !Flap Angle (Num) total given
         read (15 ,*) (FLAPANGS  (nbod,ispan,ifl), ifl=1, NFLAPANGS(nbod,ispan))!Flap Num
         read (15 ,*)
         read (15 ,*) NANGB(nbod,ispan)                                         !Angles total given

!--------- Search which span position has the most Angles given
         if ( NANGmax<NANGB(nbod,ispan) ) then
              NANGmax=NANGB(nbod,ispan)
              IspanANGmax = ispan
         endif

         do iang =1, NANGB(nbod,ispan)

            read (15 ,*)   AATB(nbod,ispan,iang),    &                              !alphaf
                         ( CLTB(nbod,ispan,ifl,iang),&                              !CL
                           CDTB(nbod,ispan,ifl,iang),&                              !CD
                           CMTB(nbod,ispan,ifl,iang), ifl=1,NFLAPANGS(nbod,ispan) ) !CM
         enddo !iang
        close (15)

       if (INDXFlap==1) then
            open (16 , file=trim(dir_airfoils)//trim(AIRGEOM(nbod,ispan)) )         !Airfoil geometry  

            read (16 ,*)  NPOINTS(nbod,ispan)                                       !Points of airfoil     
            read (16 ,*)
 
            do fts=1, NPOINTS(nbod,ispan)
               read (16,*) XAIRNODE(nbod,ispan,fts), YAIRNODE(nbod,ispan,fts)      ! X,Ynode of airfoil at ispan at nbod.
            enddo  !fts

           close (16)
       endif !INDXFlap

      enddo !ispan
     close(14)


      if (INDXFlap==1) then
!--------- Express coordinates XAIRNODE, YAIRONDE on common the common basis
!--------- of the first strip (CYLINDER) for every body
         do       ispan = 2, NSPANB2(nbod)
            do    I     = 1, NPOINTS(nbod,1)
               do fts   = 1, NPOINTS(nbod, ispan)
                  if (  XAIRNODE(nbod,1,I) <  XAIRNODE(nbod,ispan,fts) ) then
                    Y1                       = YAIRNODE(nbod, ispan, fts-1)
                    X1                       = XAIRNODE(nbod, ispan, fts-1)
                    Y2                       = YAIRNODE(nbod, ispan, fts)
                    X2                       = XAIRNODE(nbod, ispan, fts) 
                    MH(I) = XAIRNODE(nbod, 1, I)* (Y2-Y1)/(X2-X1) + (Y1*X2-Y2*X1)/(X2-X1)                     
                    exit
                  elseif (XAIRNODE(nbod, 1, I) == XAIRNODE(nbod, ispan, fts)) then
                    MH(I) = YAIRNODE(nbod, ispan, fts)
                    exit
                  endif
               enddo ! fts
            enddo ! I

            do I = 1, NPOINTS(nbod, ispan)
               XAIRNODE(nbod, ispan, I) = XAIRNODE(nbod, 1, I)
               YAIRNODE(nbod, ispan, I) = MH(I)
            enddo ! I 
         enddo !ispan
      endif !INDXFlap


!qq
                          !CLTB(:,:,:,:)=0.d0;
                          !CDTB(:,:,:,:)=0.d0;

!w    do I             = 1, NSPANB1(nbod)
!w       write(202,'(5f15.5)') RCB1(nbod,I), CHORDTB(nbod,I), TWTB(nbod,I), XAERTB(nbod,I), ZAERTB(nbod,I)
!w    enddo !I

      do i  = 1, NSTRIP
         r  = strip(nbod,i)%rBld
         if (r<RCB1(nbod,1).or.r>RCB1(nbod,NSPANB1(nbod))) then
            write(* ,*)'radius not found in geomb.inp', r, RCB1(nbod,1)
            write(10,*)'radius not found in geomb.inp'
            stop
         endif
         n      = NSPANB1(nbod  );
         RCB(:) = RCB1   (nbod,:);  !    lint1( xa   , ya , n1, n2, np    , x, y   )
         VAR(:) = CHORDTB(nbod,:);  call lint1( RCB  , VAR, 1 , n , NSPANM, r, chd );  strip(nbod,i)%chord = chd           !chord
         VAR(:) = TWTB   (nbod,:);  call lint1( RCB  , VAR, 1 , n , NSPANM, r, twi );  strip(nbod,i)%twist = twi / R2D_ae  !twist
         VAR(:) = XAERTB (nbod,:);  call lint1( RCB  , VAR, 1 , n , NSPANM, r, xae );  strip(nbod,i)%xae   = xae           !distance between elastic
         VAR(:) = ZAERTB (nbod,:);  call lint1( RCB  , VAR, 1 , n , NSPANM, r, zae );  strip(nbod,i)%zae   = zae           !and aerodynamic axis
!w       write(203,'(5f15.5)') r,chd,twi,xae,zae
      enddo !i
!w       write(202,*); write(202,*)
!w       write(203,*); write(203,*)
   enddo !nbod

!jim
!!!IspanANGmax = 1
   write(*,*) 'IspanANGmax',IspanANGmax


!--- Initialize induction factors and set fa (dynamic inflow term)
   do nbod   = 1, NBLADE_ae
   do istrip = 1, NSTRIP
      r      =  strip(nbod,istrip)%rRD/RTCONE
      call FATINT (r, fa)!Fa(r/R)

      strip(nbod, istrip)%fa          = fa
      strip(nbod, istrip)%AINDTB      = 0.3d0  !a
      strip(nbod, istrip)%AINDTB_PRE  = 0.3d0  !a
      strip(nbod, istrip)%AINDPTB     = 0.d0   !a'
      strip(nbod, istrip)%AINDPTB_PRE = 0.d0   !a'
   enddo
   enddo

   call Elast2Aero_bem_init


   if (INDXFlap==1) then
!------ Find Geometry of airfoil at radious RLD
      do    nbod   = 1, NBLADE_el
         do lstrip = 1, NSTRIP

!!          RLD=RADSTRIP(lstrip)*RTIP
            RLD=strip(nbod,lstrip)%rRD
            call AIRGEOMETRY (nbod, RLD, lstrip, RCB2, XAIRNODE, YAIRNODE, XNSTR, YNSTR, NSPANB2, NPOINTS(nbod,1), NSTRIP)

            if (IFLAPLOC(nbod, lstrip) == 1) then
               do i=1, Nflaps
                  if ( (RLD.ge.Rflaps(i,1)) .and. (RLD.le.Rflaps(i,2)) ) then
                   XFLAP(lstrip)=FLAPEX(i)
                  endif
               enddo !Nflaps
             endif ! IFLAPLOC
         enddo ! lstrip
     enddo ! nbod

!--- GWAKE wake vorticity
     allocate (GWAKE(NBLADE_el, NSTRIP, NWAKE))
     allocate (circ(NBLADE_el, NSTRIP, 0: NTIMEM))
     allocate (circ_int(NBLADE_el, NSTRIP, NNODEM))
     allocate (circ_old(NBLADE_el, NSTRIP, NNODEM))
     allocate (cmpot(NBLADE_el, NSTRIP))
     allocate (clf(NBLADE_el, NSTRIP, NTIMEM))
     allocate (cmf(NBLADE_el, NSTRIP, NTIMEM))
     allocate (cdf(NBLADE_el, NSTRIP, NTIMEM))
   endif !INDXFlap


 END Subroutine RAFT_init
!--------------------------------------------------------------------
 Subroutine Elast2Aero_bem_init
!--------------------------------------------------------------------

 Use Cbeam
 Use Craft

   implicit none

   real(8) :: RB(3), AB(3,3)
   real(8) :: HLOC,ALLOC,HTA,HLOC1,HLOC2
   integer :: nbod,ibsub,ib_el,nn0_el,nn_el,iel,lstrip ,n1,n2


   do    nbod   = 1, NBLADE_ae
   do    lstrip = 1, NSTRIP
         HLOC   = strip(nbod,lstrip)%rBld
      do ibsub  = 1, body(nbod)%NBODSUB_el
         ib_el  = body(nbod)%NBODTGLB_el(ibsub)
         nn0_el = body(nbod)%NNODTGLB_el (1)
         nn_el  = body(nbod)%NNODTGLB_el (ibsub)
         AB (1:3,1:3) = matmul ( transf_mat(nn0_el)%AT_el(1:3,1:3),  transf_mat(nn_el)%A_el(1:3,1:3)                               )
         RB (1:3)     = matmul ( transf_mat(nn0_el)%AT_el(1:3,1:3), (transf_mat(nn_el)%R_el(    1:3)-transf_mat(nn0_el)%R_el(1:3)) )
         RB (2)       = RB(2) + Hhub

!--------- Search which element HLOC lies in
         do iel   = 1, subbody(ib_el)%NTEB_el
            n1    = subbody(ib_el)%INODNEL_el(iel,1)
            n2    = subbody(ib_el)%INODNEL_el(iel,2)
            HLOC1 = beam_timo(n1)%XL_el(2) + Hhub
            HLOC2 = beam_timo(n2)%XL_el(2) + Hhub

!           if ((nbod==1).and.(lstrip==NSTRIP)) write(10,*) 'hloc1,hloc2', HLOC1,HLOC2
!           if ((nbod==1).and.(lstrip==NSTRIP)) write(* ,*) 'hloc1,hloc2', HLOC1,HLOC2

            if ( (HLOC >= HLOC1).and.(HLOC < HLOC2) ) then
               ALLOC    =  subbody(ib_el)%ALENG_el (iel)
               HTA      = (HLOC-HLOC1)/(HLOC2-HLOC1)*ALLOC

               strip(nbod,lstrip)%NSUBe2a = ibsub
               strip(nbod,lstrip)%NELe2a  = iel
               strip(nbod,lstrip)%HTAe2a  = HTA

               write(10,'(a,4i3,4f15.3)') 'e2a', nbod,lstrip, ibsub,iel,HTA,HLOC1,HLOC,HLOC2
!w             write(* ,'(a,4i3,4f15.3)') 'e2a', nbod,lstrip, ibsub,iel,HTA,HLOC1,HLOC,HLOC2

               goto 99
            endif
         enddo !iel
      enddo !ibsub

      write(* ,*)'HLOC not found',nbod,HLOC,HLOC1,HLOC2
      write(10,*)'HLOC not found',nbod,HLOC,HLOC1,HLOC2
      stop

99 enddo !lstrip
   enddo !nbod


 END Subroutine Elast2Aero_bem_init
!----------------------------------------------------------------------
 Subroutine Elast2Aero_bem
!----------------------------------------------------------------------

 Use Cbeam !FLAP_el,dFLAP_el,ddFLAP_el,TIME, UT_el,UT1_el,TILT,YAW
 Use Craft

   implicit none

   integer :: nbod, nbsub, nn_el, i, lstrip
   real(8) :: XG(3), V_RD(3), THY, THY1, THY2, THYP, pitch_ae
   real(8) :: A(3,3)
   real(8) :: FlapAng(MBLADE), FlapVel(MBLADE), FlapAcc(MBLADE), psiout(MBLADE), psi, ramp
   real(8), save :: TIME_wr=0.d0


!--- Vars communication from the structural/dynamic part
   call Elast2Aero_bem0        !- Calculation of Blade deformation angles, velocities, acclerations [UT[ ,1,2]els2aer]


      TIME_ae            = TIME
      AZIM_ae            =-UT_el (NDFBT_el+NQSW) - UThub_el
      OMEGA_ae           =-UT1_el(NDFBT_el+NQSW)
      TILT_ae            = TILT
      YAW_ae             = YAW
   do nbod               = 1, NBLADE_ae
      ATP_ae(nbod,:,:)   =           ATP_el(nbod,:,:)
      AP_ae (nbod,:,:)   = transpose(ATP_el(nbod,:,:))

      do lstrip          = 1, NSTRIP

!--------- Call subroutine that defines rigid body motion and elastic velocities
         call Elast2Aero_bem1 ( nbod, lstrip, XG, V_RD, THY, THY1, THY2, THYP, pitch_ae )

         call ROT_MATRIX0 ( 2, pitch_ae, A )

         A_pit_ae (nbod,lstrip,:,:)    =         A        (:,:)
         A        (            :,:)    = matmul (A_cone_ae(:,:), A_pit_ae(nbod,lstrip,:,:))

         strip(nbod,lstrip)%UBODY(1:3) =  V_RD(1:3)
         strip(nbod,lstrip)%THYTB      = -THY
         strip(nbod,lstrip)%THY1TB     = -THY1
         strip(nbod,lstrip)%THY2TB     = -THY2
         strip(nbod,lstrip)%pitch      = -THYP
         strip(nbod,lstrip)%A(1:3,1:3) =  matmul( A(1:3,1:3), strip(nbod,lstrip)%A_curv(1:3,1:3) ) !from RD to blade local without twist and torsion [but with pitch]
      enddo
   enddo


!--- Set instantaneous hub position
   if (NBODBT_el > NBLADE_el) then
      nbsub          = body(NSHAFT_el)%NBODSUB_el
      nn_el          = body(NSHAFT_el)%NNODTGLB_el(nbsub+1)
      HUBPOS_ae(1:3) = transf_mat(nn_el)%R_el(1:3)
   else
      HUBPOS_ae(1:3) = HUBpos_el(1:3)
   endif


!--- Set flap angle   !case 1:1P, 2:3P, 3:6P, 10:ctrl IFC (110:tune-spinner,111:spinner)
   if ( INDXFlap == 0 ) return
   if (TIME<dble(NTIME_INI)*DT_el) return

         ramp          = 1.d0
   if (TIME<dble(NTIME_INI)*DT_el+TPERIOD_el) &
         ramp          = (TIME-dble(NTIME_INI)*DT_el) /TPERIOD_el
!--- 1p, 3p and 6p collective flap motion without delay, or motion defined by the controller
   do    nbod          = 1, NBLADE_el
         psi           = -UT_el(NDFBT_el+NQSW)            !-UThub_el       !collective
!!       psi           = -UT_el(NDFBT_el+NQSW) +PHI0(nbod)!-UThub_el       !cyclic
         psiout (nbod) = dmod (psi,PI2)
         FlapAng(nbod) = 0.d0
      if     (IFLAP_typ == 1) then
         FlapAng(nbod) = Flap_ampl*cos( 1.d0*psi ) * ramp
         FlapVel(nbod) =-Flap_ampl*sin( 1.d0*psi ) * ramp *  1.d0 * UT1_el(NDFBT_el+NQSW)
         FlapAcc(nbod) =-Flap_ampl*cos( 1.d0*psi ) * ramp * (1.d0 * UT1_el(NDFBT_el+NQSW))**2
      elseif (IFLAP_typ == 2) then
         FlapAng(nbod) = Flap_ampl*cos( 3.d0*psi ) * ramp
         FlapVel(nbod) =-Flap_ampl*sin( 3.d0*psi ) * ramp *  3.d0 * UT1_el(NDFBT_el+NQSW)
         FlapAcc(nbod) =-Flap_ampl*cos( 3.d0*psi ) * ramp * (3.d0 * UT1_el(NDFBT_el+NQSW))**2
      elseif (IFLAP_typ == 3) then
         FlapAng(nbod) = Flap_ampl*cos( 6.d0*psi ) * ramp
         FlapVel(nbod) =-Flap_ampl*sin( 6.d0*psi ) * ramp *  6.d0 * UT1_el(NDFBT_el+NQSW)
         FlapAcc(nbod) =-Flap_ampl*cos( 6.d0*psi ) * ramp * (6.d0 * UT1_el(NDFBT_el+NQSW))**2
      elseif (IFLAP_typ ==10) then !set by the controller
         write(*,*)'no controller for flap';stop
         FlapAng(nbod) =   FLAP_el(nbod) * R2D
         FlapVel(nbod) =  dFLAP_el(nbod) * R2D
         FlapAcc(nbod) = ddFLAP_el(nbod) * R2D
      endif

      do lstrip = 1, NSTRIP
         if (IFLAPLOC(nbod,lstrip) == 1) FLCONTR (nbod,lstrip) = FlapAng(nbod) !degrees
      enddo
   enddo !nbod

!qq move to WRITEOUT_ae for iterations
   if(TIME>TIME_wr) then
      TIME_wr = TIME
#     if   ASCII == 1
         open(23,file='flap.dat',access='append')
          write(23,'(50f15.5)')                     &
#     elif ASCII == 0
         open(23,file='flap.bin',access='append',form='UNFORMATTED')
          write(23            )                     &
#     endif
            sngl(TIME)                            , & ! 1
           (sngl(psiout (i)*R2D )                 , & ! 2,6,10
            sngl(FlapAng(i)     )                 , & ! 3,7,11
            sngl(FlapVel(i)     )                 , & ! 4,8,12
            sngl(FlapAcc(i)     ), i=1, NBLADE_el)    ! 5,9,13
         close(23)
   endif


 END Subroutine Elast2Aero_bem
!----------------------------------------------------------------------
 Subroutine Elast2Aero_bem0 !old ELST2AER2
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   real(8) :: A   (3,3), T0 (3), T0p (3)
   real(8) :: AT0 (3,3), T1 (3), T1p (3)
   real(8) :: ATPA(3,3), T2 (3), T2p (3)
   integer :: nbod, nbsub, nel, nb_el, nn_el, nnod, i, j !!, nn0_el!, nn, ii
   real(8), dimension(NDFPEM_el) :: UT, UT1, UT2


   do nbod          = 1, NBLADE_el
!!    nn0_el        = body(nbod)%NNODTGLB_el(1)
!!    AT0(1:3,1:3)  = transf_mat(nn0_el)%AT_el(1:3,1:3) !after cone, after pitch, afte   precurve rotation
!!    AT0(1:3,1:3)  = ATPWr_el(nbod,1:3,1:3)            !after cone, after pitch, before precurve rotation
      AT0(1:3,1:3)  = ATP_el  (nbod,1:3,1:3)            !before cone [wrt RD]

      T0 (1:3) = 0.d0;  T1 (1:3) = 0.d0;  T2 (1:3) = 0.d0;
      T0p(1:3) = 0.d0;  T1p(1:3) = 0.d0;  T2p(1:3) = 0.d0;

      nb_el         = body(nbod)%NBODTGLB_el(1)

      UTels2aer  (nb_el,1:3) = 0.d0;
      UT1els2aer (nb_el,1:3) = 0.d0;
      UT2els2aer (nb_el,1:3) = 0.d0;

      do nbsub      = 1, body(nbod)%NBODSUB_el - 1
         nb_el      = body(nbod)%NBODTGLB_el(nbsub)
         nn_el      = body(nbod)%NNODTGLB_el(nbsub)

         A(1:3,1:3) = transf_mat(nn_el)%A_el(1:3,1:3)
         ATPA       = matmul( AT0,A );
!        ATPA       =         transf_mat(nn_el)%ATPA_el;    !from rotor-plane before pitch to blade c.s.

         nel        = subbody(nb_el)%NTEB_el
         nnod       = NNPE_el
         i          = subbody(nb_el)%NDFPNACC_el (nnod-1)

         call LOCAL_UT_1 ( nb_el, nel, UT, UT1, UT2 )

         do j       = 1, 3
            T0(j)   = UT (i+IMDOF_el(j+3))
            T1(j)   = UT1(i+IMDOF_el(j+3))
            T2(j)   = UT2(i+IMDOF_el(j+3))
         enddo

         T0         = matmul( ATPA, T0 )
         T1         = matmul( ATPA, T1 )
         T2         = matmul( ATPA, T2 )

         T0p(1:3)   = T0p(1:3) + T0(1:3)
         T1p(1:3)   = T1p(1:3) + T1(1:3)
         T2p(1:3)   = T2p(1:3) + T2(1:3)

         UTels2aer  (nb_el+1,1:3) = T0p(1:3)
         UT1els2aer (nb_el+1,1:3) = T1p(1:3)
         UT2els2aer (nb_el+1,1:3) = T2p(1:3)
      enddo !nbsub
   enddo !nbod


 END Subroutine Elast2Aero_bem0 !old ELST2AER2
!----------------------------------------------------------------------
! Provides
!   XG      : the global position vector of the strip
!   V_RD    : the total body velocity wrt the RD
!   THY     : the blade torsion deflection wrt the blade local c.s.
!   THY[1,2]: the blade torsion & pitch [velocity,acceleration]
!             wrt the blade local c.s.
!   THYP    : the blade pitch angle projected on the blade local c.s.
!----------------------------------------------------------------------
 Subroutine Elast2Aero_bem1 ( nbod, lstrip, XG, V_RD, THY, THY1, THY2, THYP, pitch_ae ) !old LOCAL_UT_2
!----------------------------------------------------------------------

 Use Cbeam
 Use Craft

   implicit none

   integer, intent(in ) :: nbod, lstrip
   real(8), intent(out) :: XG(3), V_RD(3), THY, THY1, THY2, THYP, pitch_ae

   real(8), dimension(NEQPEM_el,NDFPEM_el) :: SHAPE
   real(8), dimension(NDFPEM_el          ) ::    UT,    UT1,    UT2
   real(8), dimension(NEQPEM_el          ) ::   UTT,   UTT1,   UTT2
!  real(8), dimension(NEQPEM_el,NDFPEM_el) :: DSHAPE                 !for Euler-Bernoulli
!  real(8), dimension(NEQPEM_el          ) ::  eUTT,  eUTT1,  eUTT2  !for Euler-Bernoulli
!  real(8), dimension(NEQPEM_el          ) :: eDUTT, eDUTT1, eDUTT2  !for Euler-Bernoulli
!  real(8), dimension(3                  ) ::    A0,     A1,     A2  !for Euler-Bernoulli
   real(8), dimension(3                  ) :: VG, AKSI0, TH_Bld, TH1_Bld, TH2_Bld, THP,THP1,THP2

   integer :: nbsub, nb_el, nn_el, nel, NEQPE, NDFPE, i
   real(8) :: ALLOC, HTA0, HTA, HTAL, PHPX, PHPZ,  SS(3,6), AT(3,3), ATPA (3,3)
!  integer :: nn0_el
!  real(8) :: ATPA0 (3,3)


   nbsub     = strip(nbod,lstrip)%NSUBe2a
   nel       = strip(nbod,lstrip)%NELe2a
   HTA       = strip(nbod,lstrip)%HTAe2a
   nb_el     = body(nbod)%NBODTGLB_el(nbsub)
   nn_el     = body(nbod)%NNODTGLB_el(nbsub)
   HTA0      = subbody(nb_el)%HTA_el (nel)
   NEQPE     = subbody(nb_el)%NEQPE_el
   NDFPE     = subbody(nb_el)%NDFPE_el
   PHPX      = subbody(nb_el)%PHIX_el (nel)
   PHPZ      = subbody(nb_el)%PHIZ_el (nel)
   ALLOC     = subbody(nb_el)%ALENG_el(nel)
   HTAL      = HTA/ALLOC
   AKSI0 (:) = 0.d0;
!qqAKSI0 (1) = strip(nbod,lstrip)%xae
   AKSI0 (2) = HTA0+HTA
!qqAKSI0 (3) = strip(nbod,lstrip)%zae

   if     (body(nbod)%IDOF_el==0) then

      UTT   = 0.d0;
      UTT1  = 0.d0;
      UTT2  = 0.d0;

   elseif (IORDbl_el==1) then

      call SSHAPEFUNC15 ( HTAL, ALLOC, PHPX, PHPZ, SHAPE, NDFPE, NEQPE )
      call LOCAL_UT_1   ( nb_el, nel, UT, UT1, UT2 )

      UTT   = matmul (SHAPE, UT )
      UTT1  = matmul (SHAPE, UT1)
      UTT2  = matmul (SHAPE, UT2)

!!!elseif (IORDbl_el==2) then

!!!   call DSHAPEFUNC15eu ( HTAL, ALLOC, SHAPE, DSHAPE, NDFPE ,NEQPE, NDFPEM_el,NEQPEM_el )
!!!   call LOCAL_UT_1     (nb_el, nel, UT, UT1, UT2)

!!!   eUTT  (1:NEQPE) = matmul ( SHAPE (1:NEQPE,1:NDFPE), UT (1:NDFPE) )
!!!   eUTT1 (1:NEQPE) = matmul ( SHAPE (1:NEQPE,1:NDFPE), UT1(1:NDFPE) )
!!!   eUTT2 (1:NEQPE) = matmul ( SHAPE (1:NEQPE,1:NDFPE), UT2(1:NDFPE) )
!!!   eDUTT (1:NEQPE) = matmul ( DSHAPE(1:NEQPE,1:NDFPE), UT (1:NDFPE) )
!!!   eDUTT1(1:NEQPE) = matmul ( DSHAPE(1:NEQPE,1:NDFPE), UT1(1:NDFPE) )
!!!   eDUTT2(1:NEQPE) = matmul ( DSHAPE(1:NEQPE,1:NDFPE), UT2(1:NDFPE) )

!!!   A0 (1) = eDUTT(3);   A1 (1) = eDUTT1(3);   A2 (1) = eDUTT2(3);
!!!   A0 (2) = eUTT (4);   A1 (2) = eUTT1 (4);   A2 (2) = eUTT2 (4);
!!!   A0 (3) = eDUTT(1);   A1 (3) = eDUTT1(1);   A2 (3) = eDUTT2(1);

!!!   UTT (IMDOF_el(1:3)) = eUTT (1:3);   UTT (IMDOF_el(4:6)) = A0 (1:3);
!!!   UTT1(IMDOF_el(1:3)) = eUTT1(1:3);   UTT1(IMDOF_el(4:6)) = A1 (1:3);
!!!   UTT2(IMDOF_el(1:3)) = eUTT2(1:3);   UTT2(IMDOF_el(4:6)) = A2 (1:3);
   endif

!--- Global position vector and body velocity
   call Calc_SS ( SS, AKSI0 )
   XG     (1:3) = transf_mat(nn_el)% R_el(1:3) + matmul ( transf_mat(nn_el)% A_el(1:3,1:3), matmul( SS(1:3,1:6), UTT (IMDOF_el(1:6)) ) + AKSI0(1:3) )
   VG     (1:3) = transf_mat(nn_el)%DR_el(1:3) + matmul ( transf_mat(nn_el)%DA_el(1:3,1:3), matmul( SS(1:3,1:6), UTT (IMDOF_el(1:6)) ) + AKSI0(1:3) ) &
                                               + matmul ( transf_mat(nn_el)% A_el(1:3,1:3), matmul( SS(1:3,1:6), UTT1(IMDOF_el(1:6)) )              )

!--- Body velocity wrt rotor plane
   V_RD   (1:3) = matmul ( ATP_el (nbod,1:3,1:3), VG (1:3) )


!--- (Blade only) Elastic Rotational deformation, Velocity and Acceleration wrt blade local c.s.
   ATPA (1:3,1:3) = transf_mat(nn_el )%ATPA_el(1:3,1:3)           !from Bld(elast) to RD
   AT   (1:3,1:3) = transpose( strip(nbod,lstrip)%A(1:3,1:3) )    !from RD to Bld(aero)
!  AT   (1:3,1:3) = matmul ( transf_mat(nn_el)% AT_el(1:3,1:3), transpose(ATP_el (nbod,1:3,1:3)) )  !from RD to Bld
!----                         wrt RD                                                           wrt blade local c.s.
   TH_Bld   (1:3) = UTels2aer  (nb_el,1:3) + matmul( ATPA(1:3,1:3), UTT (IMDOF_el(4:6)) );     TH_Bld (1:3) = matmul ( AT(1:3,1:3), TH_Bld (1:3) )
   TH1_Bld  (1:3) = UT1els2aer (nb_el,1:3) + matmul( ATPA(1:3,1:3), UTT1(IMDOF_el(4:6)) );     TH1_Bld(1:3) = matmul ( AT(1:3,1:3), TH1_Bld(1:3) )
   TH2_Bld  (1:3) = UT2els2aer (nb_el,1:3) + matmul( ATPA(1:3,1:3), UTT2(IMDOF_el(4:6)) );     TH2_Bld(1:3) = matmul ( AT(1:3,1:3), TH2_Bld(1:3) )


!--- Pitch Position, Velocity and Acceleration (as set by the controller)
   i               = NDFBT_el + NQSP + ACCU%NQSACC_el(nbod-1) + 5
   pitch_ae        = UT_el  (i)
!--- handle air-brake through modification of the pitch angle
!!if ((ICMOD==3.or.ICMOD==4).and.IAX_Ben_Swe == 2) &   !transform the pre-curved angle around pitch axis
  if (IAX_Ben_Swe == 2) &   !transform the pre-curved angle around pitch axis
   pitch_ae        = UT_el  (i) + body(nbod)%PRECURV(nbsub)
!                                                                          !wrt RD             wrt blade local c.s.
   THP  = 0.d0;   THP (2) = pitch_ae  ;   THP (1:3) = matmul(A_cone_ae(1:3,1:3), THP (1:3));   THP (1:3) = matmul ( AT(1:3,1:3), THP (1:3) )
   THP1 = 0.d0;   THP1(2) = UT1_el (i);   THP1(1:3) = matmul(A_cone_ae(1:3,1:3), THP1(1:3));   THP1(1:3) = matmul ( AT(1:3,1:3), THP1(1:3) )
   THP2 = 0.d0;   THP2(2) = UT2_el (i);   THP2(1:3) = matmul(A_cone_ae(1:3,1:3), THP2(1:3));   THP2(1:3) = matmul ( AT(1:3,1:3), THP2(1:3) )


!--- output the local y component
   THYP =              THP (2)
   THY  = TH_Bld (2)
   THY1 = TH1_Bld(2) + THP1(2)
   THY2 = TH2_Bld(2) + THP2(2)


 END Subroutine Elast2Aero_bem1 !old LOCAL_UT_2
!
!
!
!----------------------------------------------------------------------------------
!--- Performs Blade Element Momentum Theory Calculations
!----------------------------------------------------------------------
 Subroutine RAFT ( NTIME )
!----------------------------------------------------------------------

 Use Craft

   implicit none

   integer, intent(in)  :: NTIME

   real(8), allocatable :: x(:) ![2 o r2*NBLADE_ae]
   real(8) :: xRD
   logical :: check
   integer :: nb, nbo, i0, BEM_CPL, BEM_eq, nb_1, nb_2
!--- skewed wake modeling
   real(8) :: gama_sk, xi_sk, psi0_sk
   real(8) :: vec1(3),vec2(3),vec_ep(3), ayaw, atilt, prod_ext, prod_int, AMEAN, AMEAN_PRE
   external:: BEMT

   integer ::     NTIMEn,lstrip,ibladen
   common /status/NTIMEn,lstrip,ibladen
   save   /status/


   call Elast2Aero_bem         !- Vars communication from structural/dynamic part
   call UINFLOW_bem            !- Inflow velocity estimation
   call foilfs_init (NTIME)


!--- Solve BEM equations with Newton-Rapson: define a, a' and the aerodynamic loading
   do    lstrip       = 1, NSTRIP
!--------- skewed inflow angle gama_sk (i.e yaw+tilt)
         ayaw         = WINDYAW + YAW_ae
         atilt        = 0.d0 !WINDINC + TILT_ae
         vec1(1  )    = 1.d0;
         vec1(2:3)    = 0.d0;

         call ROTATEP(3,-ayaw ,vec1,vec2  )
         call ROTATEP(2,-atilt,vec2,vec2  )
         call EXTEPR (   vec1 ,vec2,vec_ep)

         prod_ext     = dsqrt(dot_product(vec_ep,vec_ep))
         prod_int     = dot_product(vec1(1:3),vec2(1:3))  
         gama_sk      = datan2(prod_ext,prod_int)
         psi0_sk      = datan2(vec2(3) ,-vec2(2))
!qq      psi0_sk      = 0.d0
!        gama_sk      = WINDYAW !+YAW, WINDINC + TILT

!w       write(*,*) 'psi0',psi0_sk*R2D_ae
!w       write(*,*) 'gama',gama_sk*R2D_ae, WINDYAW*R2D_ae!;stop

!--------- Wake skewed angle xi_sk
         AMEAN        = 0.d0
         AMEAN_PRE    = 0.d0
        do nb         = 1, NBLADE_ae
         AMEAN        = AMEAN     + strip(nb, lstrip)%AINDTB     / dble(NBLADE_ae)
         AMEAN_PRE    = AMEAN_PRE + strip(nb, lstrip)%AINDTB_PRE / dble(NBLADE_ae)
        enddo !nb
         xi_sk        = datan2(dsin(gama_sk), dcos(gama_sk)*(1.d0-AMEAN))

         strip(1:NBLADE_ae,lstrip)%xi_sk          = xi_sk
         strip(1:NBLADE_ae,lstrip)%gama_sk        = gama_sk
         strip(1:NBLADE_ae,lstrip)%psi0_sk        = psi0_sk
         strip(1:NBLADE_ae,lstrip)%AINDTB_PREmean = AMEAN_PRE


!---------- Decoupled or Coupled BEM
          NTIMEn  = NTIME
          BEM_CPL = max(IDYNWAKE-1,0)
          if (NTIME<2.or.lstrip<0) &
          BEM_CPL = 0
10    if (BEM_CPL==0) then
          nb_1    = 1
          nb_2    = NBLADE_ae
          BEM_eq  = 2
      else
          nb_1    = 0
          nb_2    = 0
          BEM_eq  = 2*NBLADE_ae
      endif

          Allocate (x(BEM_eq))

!------ call Newton
      do nb       = nb_1, nb_2
         ibladen  = nb
      do nbo      = 1, NBLADE_ae
         i0       = (nbo-1)*2
         if (BEM_CPL==0) &
         i0       = 0
         x(i0+1)  = 0.01d0   !a
         x(i0+2)  = 0.00d0   !a'
      enddo !nbo

!w    write(*,*)'RAFT',lstrip
      call newt (x,BEM_eq,check) !calls BEMT uses NTIMEn,lstrip,ibladen

      if (check.eqv..true.) then
         if(BEM_CPL==1) then
            BEM_CPL= 0
            write(* ,'(a,2i4)')'Attention, Newton local solution found in RAFT Coupled',lstrip,nb
            write(10,'(a,2i4)')'Attention, Newton local solution found in RAFT Coupled',lstrip,nb
            Deallocate(x)
            goto 10
         else
            write(* ,'(a,2i4)')'Attention, Newton local solution found in RAFT Decoupled',lstrip,nb
            write(10,'(a,2i4)')'Attention, Newton local solution found in RAFT Decoupled',lstrip,nb
         endif
      endif

!     call aero_loa_bem (iblade, lstrip)
!qq   call VIV
      enddo !nb
      Deallocate(x)

      if(BEM_CPL==0) then
         AMEAN    = 0.d0
        do nb     = 1, NBLADE_ae
         AMEAN    = AMEAN + strip(nb, lstrip)%AINDTB / dble(NBLADE_ae)
        enddo !nb
        strip(1:NBLADE_ae,lstrip)%AINDTBmean = AMEAN
      endif
   enddo !lstrip


!--------- Calculate Blade/Rotor Integrated Loads
         TTHRUST    =0.d0; TTORQUE   =0.d0; POWER=0.d0;
          THRUST(:) =0.d0;  TORQUE(:)=0.d0;
   do    nb         = 1, NBLADE_ae
      do lstrip     = 1, NSTRIP
         xRD        = (node(nb,lstrip+1)%xRD (1)+node(nb,lstrip)%xRD (1)) / 2.d0
         THRUST(nb) = THRUST(nb) + strip(nb,lstrip)%FCP(3)*strip(nb,lstrip)%drBLD
         TORQUE(nb) = TORQUE(nb) + strip(nb,lstrip)%FCP(1)*strip(nb,lstrip)%drBLD*strip(nb,lstrip)%rRD &
                                 - strip(nb,lstrip)%FCP(2)*strip(nb,lstrip)%drBLD*                 xRD
      enddo !lstrip
!--------- Total Aerodynamic Thrust, Torque and Power calculation
         TTHRUST    = TTHRUST + THRUST(nb)
         TTORQUE    = TTORQUE + TORQUE(nb)
         POWER      = POWER   + TORQUE(nb) * OMEGA_ae !qq elastic velocity is not included in OMEGA herein (otherwise power for each strip should be calculated using its elastic velocity)
   enddo !nb

   call Aero2Elast_bem         !- Communicate the aerodynamic loads to the structural/dynamic part


 END Subroutine RAFT
!-------------------------------------------------------------------------
 Subroutine BEMT (n,x,ff) !called from newton
!-------------------------------------------------------------------------

 Use Craft
 Use Caero !C[L,D,M]_ae
 Use Foilfs_mod

   implicit none

   integer, intent(in ) :: n
   real(8), intent(in ) :: x (n)
   real(8), intent(out) :: ff(n)

   integer :: INOBEM,IDYNWA,istrip,IhighA, IYC, iblade, nb, i0, nb_1, nb_2
   real(8) :: UINDZ,torsion,dtorsion,ddtorsion,pitch,chord,twist
   real(8) :: F,fe,APREV,APPREV,A,AP,WEFF,PHI,ALPHA, PHI_RD
   real(8) :: AMACH,CLIFTDS, CDRAGDS, CMOMDS
   real(8) :: Ct,Cn,Cm,SIG
   real(8) :: CTrst, bita
   real(8) :: visc, Re, dr, ds
   real(8) :: rRD, rBld
   real(8) :: lamda_r, WdiaV_2
   real(8) :: AA(3,3), AT(3,3)
   real(8) :: WEFF_RD(3), WEFFP_RD(3), WEFFP_Bld(3), Ct_RD, Cn_RD
   real(8) :: ctad_dy, ctad, cmad, ctbe, cmbe, tsr
   real(8) :: CONST
   real(8) :: FCP(3), FCPL(3), MCP(3), MCPL(3)
   real(8) :: Amean, APREVmean
!--- skewed wake modeling
   real(8) :: gama_sk, xi_sk, psi0_sk, mi_sk, F_sk, K_sk, Coef_sk, psi
   integer :: NTIME,lstrip,ibladen
!--- for FOILFS
   real(8) :: ANGMOVE, AFLAP, ALPHA0

   common /status/NTIME,lstrip,ibladen
   save   /status/


      IYC          = 5
      visc         = 1.511d-05

      Amean        = 0.d0
   if (ibladen==0) then
      nb_1         = 1
      nb_2         = NBLADE_ae
   do nb           = 1, NBLADE_ae
      Amean        = Amean + x((nb-1)*2+1) / dble(NBLADE_ae)
   enddo !nb
   else
      nb_1         = ibladen
      nb_2         = ibladen
   endif
      
   do iblade       = nb_1, nb_2
      i0           = (iblade-1)*2
      if (ibladen/=0) &
      i0           = 0
      AA (1:3,1:3) = strip(iblade,lstrip)%A(1:3,1:3)
      AT (1:3,1:3) = transpose(AA(1:3,1:3))
      rRD          = strip(iblade,lstrip)%rRD
      rBld         = strip(iblade,lstrip)%rBld
      chord        = strip(iblade,lstrip)%chord
      ds           = strip(iblade,lstrip)%drBld
      dr           = strip(iblade,lstrip)%drRD
      twist        = strip(iblade,lstrip)%twist
      pitch        = strip(iblade,lstrip)%pitch
      torsion      = strip(iblade,lstrip)%THYTB
     dtorsion      = strip(iblade,lstrip)%THY1TB
    ddtorsion      = strip(iblade,lstrip)%THY2TB
      UINDZ        = strip(iblade,lstrip)%UWAKIZTB
      WEFF_RD(1:3) = strip(iblade,lstrip)%UWINDTB(1:3) - strip(iblade,lstrip)%UBODY(1:3)
      tsr          = maX(OMEGA_ae*RTCONE/HUB_VEL(1), 0.1d0)

      INOBEM       = INOINDUCT !0
      IDYNWA       = IDYNWAKE  !1
!     if (lstrip<10) &
      if (ibladen/=0)&
      IDYNWA       = min(IDYNWA,1)
      if (NTIME<2) &
      IDYNWA       = 0
      istrip       = (iblade-1)*NSTRIP + lstrip
      APREV        = strip(iblade,lstrip)%AINDTB_PRE     !a prev time step
      APPREV       = strip(iblade,lstrip)%AINDPTB_PRE    !a'prev time step
      APREVmean    = strip(iblade,lstrip)%AINDTB_PREmean !a prev time step mean value for dynamic inflow
      A            = x(i0+1)                             !a
      AP           = x(i0+2)                             !a'
      IhighA       = 0

      if ( (rBld < RHUB) .or. (IPARKED == 1) .or. (INOBEM == 1) ) then
         A         = 0.d0
         AP        = 0.d0
         INOBEM    = 1
         F         = 1.d0
      endif

!--------------------------------------------------------------
!------ Cylinder skewed wake model
      xi_sk        = strip(iblade,lstrip)%xi_sk
      gama_sk      = strip(iblade,lstrip)%gama_sk
      psi0_sk      = strip(iblade,lstrip)%psi0_sk

!------ Flow expansion function F_sk(mi_sk) and K-function K(xi_sk)
      mi_sk        = rRD/RTCONE
     if    (IYC==0) then
      F_sk         = mi_sk
      K_sk         = 0.d0
     elseif(IYC==1) then
      F_sk         = mi_sk + 0.4d0*mi_sk**3 + 0.4d0*mi_sk**5                         !1.Oye         , 1992
      K_sk         = 1.d0                                         ! dtan(xi_sk/2.d0) !1.Colleman    , 1945
     elseif(IYC==2) then
      F_sk         = mi_sk                                                           !2.Glauert     , 1926
      K_sk         = 4.d0/3.d0*(1.d0-1.8d0*(dsin(gama_sk)/tsr)**2)!*dtan(xi_sk/2.d0) !2.Meijer-Drees, 1949
     elseif(IYC==3) then
      F_sk         = mi_sk                                                           !3.Glauert     , 1926
      K_sk         = 1.d0                                         ! dtan(xi_sk/2.d0) !3.Colleman    , 1945
     elseif(IYC==4) then
      F_sk         = mi_sk                                                           !3.Glauert     , 1926
      K_sk         = 0.736d0                                      ! dtan(xi_sk/2.d0) !
     elseif(IYC==5) then
      F_sk         = mi_sk                                                           !3.Glauert     , 1926
      K_sk         = 0.5d0                                        ! dtan(xi_sk/2.d0) !
     endif

!     K_sk         = 0.736d0                                      !*dtan(xi_sk/2.d0) !  Pitt-Peters , 1981 [pi*15/64]

      psi          = AZIM_ae + PHI0_ae(iblade) - psi0_sk
      Coef_sk      = (1.d0 + K_sk*F_sk*dtan(xi_sk/2.d0)*dsin(psi))
      if (NTIME<2) &
      Coef_sk      = 1.d0
      A            = A*Coef_sk
!--------------------------------------------------------------
!--------- limit1 a, a'
      if (dabs(AP)>0.25) then
!!       write(* ,'(a,2i4)')"circ induction factor a'> 0.25",iblade,lstrip
!w       write(10,'(a,2i4)')"circ induction factor a'> 0.25",iblade,lstrip
         AP        = max(min(AP,0.25d0),-0.25d0)
      endif
      if (dabs(A )>0.99d0) then
!!       write(* ,'(a,2i4)')'axial induction factor a > 0.99',iblade,lstrip
!w       write(10,'(a,2i4)')'axial induction factor a > 0.99',iblade,lstrip
         A         = 0.99d0 !APREV
         AP        = APPREV
         IhighA    = 1
      endif
!     A            = max(min(A ,0.95d0),-0.95d0)
!     AP           = max(min(AP,0.25d0),-0.25d0)
!------ limit2 WEFF_RD(3)
      if (IPARKED==0.and.INOINDUCT==0)&
      WEFF_RD  (3  ) = max(WEFF_RD(3),0.1d0)
      WEFFP_RD (1  ) = WEFF_RD(1)*(1.d0+AP)
      WEFFP_RD (2  ) = WEFF_RD(2)
      if (A<=1.d0) &
      WEFFP_RD (3  ) = WEFF_RD(3)*(1.d0-A )
      if (A >1.d0) &
      WEFFP_RD (3  ) =-WEFF_RD(3)*(1.d0-A )
      WEFFP_Bld(1:3) = matmul(AT(1:3,1:3), WEFFP_RD(1:3))
!------ limit3 WEFFP_Bld(3)
!     if (IPARKED==0.and.INOINDUCT==0)&
!     WEFFP_Bld(3  ) = dmax1(WEFFP_Bld(3),1.0d0)
!     WEFFP_RD (1:3) = matmul(AA(1:3,1:3), WEFFP_Bld(1:3))


!------ find AoA and Weff wrt chord c.s.
      WEFF   = dsqrt (WEFFP_Bld(1)**2+WEFFP_Bld(3)**2)
      PHI    = datan2(WEFFP_Bld(3),   WEFFP_Bld(1)   )
      PHI_RD = datan2(WEFFP_RD (3),   WEFFP_RD (1)   )  !for tip/hub loss correction
      AMACH  = dmin1 (WEFF/SSPEED ,   0.8d0          )
      Re     = WEFF*chord /visc
      ALPHA  = PHI-twist-torsion


10    if     (ALPHA  >  PI_ae) then
         ALPHA = -PI2_ae + ALPHA
         goto 10
      elseif (ALPHA  < -PI_ae) then
         ALPHA =  PI2_ae + ALPHA
         goto 10
      endif


!------ Calculation of the local Cl,Cd,Cm [CL_ae,CD_ae,CM_ae]
      call LOCAL_AERO_PARAM (iblade,rBld,ALPHA,AMACH,lstrip)


!------ Calculate TIP & HUB losses
         F     = 1.d0
!------ Tiploss Calculation
      if (ITIPLOS /= 0 .and. INOBEM == 0) then
         fe    = max(dble(NBLADE_ae)*(RTCONE/rRD-1.d0) / (2.d0 *dabs(dsin(PHI_RD)+1.d-06)), -50.d0) !Prandtl
!        fe    =  dble(NBLADE_ae)*(1.d0-rRD/RTCONE) / (2.d0 *(1.d0-A)) * tsr                        !Takis
         F     = (2.d0/PI_ae)*dacos(dexp(-fe))
      endif
!------ Hubloss calculation
      if (IHUBLOS /= 0 .and. INOBEM == 0) then
         fe    = max(dble(NBLADE_ae)*(rRD/RHCONE-1.d0) / (2.d0 *dabs(dsin(PHI_RD)+1.d-06)), -50.d0)
         F     = F * (2.d0/PI_ae)*dacos(dexp(-fe))
      endif


!---- FOILFS
    if ( IFLAPLOC(iblade, lstrip) == 1 )  then          
         ANGMOVE =-torsion  !- (THY-THY0)  
         AFLAP   = FLCONTR (iblade, lstrip)
         ALPHA0  = PHI-twist
!        ALPHA0  = PHI-THET-THY0 !phi-twi-pitch
         
         call FOILFS (iblade, lstrip, NTIME, ANGMOVE, AFLAP, ALPHA0, WEFF, chord)
         
         CLIFTDS = clf (iblade, lstrip, NTIME)
         CDRAGDS = cdf (iblade, lstrip, NTIME)
         CMOMDS  = cmf (iblade, lstrip, NTIME)

      if     ( (IDYNSTALL == 1).and.(lstrip>IDYNSTRIP) ) then
         call ONERALDM (rRD, chord, WEFF, ALPHA , torsion+pitch, dtorsion, ddtorsion, &
                        CLIFTDS, CDRAGDS, CMOMDS, istrip, lstrip,iblade, NTIME          )
      endif
   else
!------ Dynamic Stall model
      if     ( (IDYNSTALL == 1).and.(lstrip>IDYNSTRIP) ) then
         call ONERALDM (rRD, chord, WEFF, ALPHA , torsion+pitch, dtorsion, ddtorsion, &
                        CLIFTDS, CDRAGDS, CMOMDS, istrip, lstrip,iblade, NTIME          )
      elseif ( (IDYNSTALL == 2).and.(lstrip>IDYNSTRIP) ) then
         call BEDDLEISH(rRD, chord, WEFF, ALPHA , torsion+pitch, &
                        CLIFTDS, CDRAGDS, CMOMDS, istrip, NTIME   )
      else
         CLIFTDS = CL_ae
         CDRAGDS = CD_ae
         CMOMDS  = CM_ae
      endif
   endif !FLAPLOC

!------ Cn = normal, Ct = tangential at BLD
      Ct    = -( CDRAGDS*dcos(PHI) - CLIFTDS*dsin(PHI) )   !rotate by -PHI
      Cn    =    CDRAGDS*dsin(PHI) + CLIFTDS*dcos(PHI)
      Cm    =    CMOMDS
!------ Cn = normal, Ct = tangential at RD
      Ct_RD = -(AA(1,1) * (-Ct)  +  AA(1,3) * Cn)
      Cn_RD =   AA(3,1) * (-Ct)  +  AA(3,3) * Cn


!---- Glaurt's correction CT= (0.425+1.39a)F                ,a>0.33
!---- Takis 's correction CT= (0.64 +0.8 a)F                ,a>0.40
!---- Buhl's   correction CT= 8/9+(4F-40/9)a+(50/9-4F)a**2  ,a>0.40
!--------- BEMT equations
      if (INOBEM == 0.and.IhighA==0) then
         SIG      = dble(NBLADE_ae)*chord/(PI2_ae*rRD) * (ds/dr)  !- modified local solidity (including dr_Bld/dr_RD)
         lamda_r  = WEFF_RD(1)/WEFF_RD(3)
         WdiaV_2  = WEFF**2/WEFF_RD(3)**2
         bita     = 0.33d0
        if     (dabs(A)<=bita) then
                                   if    (     A <0.d0 ) then
!!                                   write(* ,'(a,2i4,4f15.5)')'axial induction factor < 0',iblade,lstrip,PHI_RD,PHI,A,AP
!w                                   write(10,'(a,2i4,4f15.5)')'axial induction factor < 0',iblade,lstrip,PHI_RD,PHI,A,AP
                                   endif
         CTrst    = 4.d0*A*(1.d0-A)*F                                                                      !Glauert's theory default w/o skewed wake
!        CTrst    = 4.d0*A* dsqrt( 1.d0-A*(2.d0*dcos(gama_sk)-A) )*F                                       !Glauert's theory
!        CTrst    = 4.d0*A*F* dcos(gama_sk) + dtan(xi_sk/2.d0)*dsin(gama_sk) - A*dcos(xi_sk/2.d0)**(-2.d0) !Vortex    theory
        else!if(dabs(A)<=1.d0) then
         CTrst    = (0.425d0+1.39d0*A)*F !*dcos(gama_sk)
!        CTrst    = (Ct1-4.d0*(sqrt(Ct1)-1.d0)*(1.d0-A))*F, bita     = 1.d0-sqrt(Ct1)/2.d0 Wind Energy HandBook for Ct1=1.816 = Katopis
!        CTrst    = (0.64d0 +0.8d0 *A)*F                                  ! a>0.40
!        CTrst    = (8.d0 + (36.d0*F-40.d0)*A+(50.d0-36.d0*F)*A**2)/9d0   ! a>0.40
       !else !if (     A >1.d0 ) then 
       !                             write(* ,'(a,2i4,4f15.5)')'axial induction factor > 1',iblade,lstrip,PHI_RD,PHI,A,AP
       !                             write(10,'(a,2i4,4f15.5)')'axial induction factor > 1',iblade,lstrip,PHI_RD,PHI,A,AP
       !                           ! write(* ,'(a,2i4,4f15.5)')'propeller-brake region    ',iblade,lstrip,PHI_RD,PHI,A,AP;
       !                           ! write(10,'(a,2i4,4f15.5)')'propeller-brake region    ',iblade,lstrip,PHI_RD,PHI,A,AP;
       ! CTrst    =-4.d0*A*(1.d0-A)*F    !*dcos(gama_sk)                                                   !propeller-brake
!      ! CTrst    = (0.425d0+1.39d0*A)*F !*dcos(gama_sk)
        endif !A<=bita

        if     (IDYNWA == 0) then
         ctad_dy  = 0.d0
        elseif (IDYNWA == 1) then
         ctad_dy  = 4.d0*RTCONE*strip(iblade,lstrip)%fa / WEFF_RD(3) * (A    /Coef_sk-APREV    )/DT_ae
        elseif (IDYNWA == 2) then
         ctad_dy  = 4.d0*RTCONE*strip(iblade,lstrip)%fa / WEFF_RD(3) * (Amean/Coef_sk-APREVmean)/DT_ae
        endif
         ctad     = ctad_dy + CTrst                                                                        !thrust equ. actuator disk
         cmad     = 4.d0*AP*(1.d0-A)*lamda_r*F                                                             !torque equ. actuator disk
         ctbe     = SIG * Cn_RD * WdiaV_2                                                                  !thrust equ. blade element
         cmbe     = SIG * Ct_RD * WdiaV_2                                                                  !torque equ. blade element
         ff(i0+1) = ctad - ctbe
         ff(i0+2) = cmad - cmbe
      else !INOBEM==1.or.IhighA==1
         ff(i0+1) = 0.d0
         ff(i0+2) = 0.d0
      endif

!------ Calculate loads
      CONST       = 0.5d0 * AIRDEN * WEFF**2 * chord
!------ The loads relative to the Blade local c.s.
      FCPL (:)    = 0.d0;
      FCPL (1)    =-CONST*Ct        !--- Cyrcumferential
      FCPL (3)    = CONST*Cn        !--- Axial
      MCPL (:)    = 0.d0;
      MCPL (2)    = CONST*Cm*chord  !--- Pitching Moment

!------ The loads on the strip relative to the rotor plane, before pitch rotation
      FCP = matmul(AA,FCPL);
      MCP = matmul(AA,MCPL);

!------ Store strip results
      strip(iblade,lstrip)%AINDTB    = A /Coef_sk
      strip(iblade,lstrip)%AINDPTB   = AP
      strip(iblade,lstrip)%ALPHATB   = ALPHA
      strip(iblade,lstrip)%PHITB     = PHI
      strip(iblade,lstrip)%WEFFXTB   = WEFFP_Bld(1)
      strip(iblade,lstrip)%WEFFZTB   = WEFFP_Bld(3)
      strip(iblade,lstrip)%Mach      = AMACH
      strip(iblade,lstrip)%Reynolds  = Re
      strip(iblade,lstrip)%CLIFTTB   = CLIFTDS
      strip(iblade,lstrip)%CDRAGTB   = CDRAGDS
      strip(iblade,lstrip)%CNORMTB   = Cn_RD
      strip(iblade,lstrip)%CTANGTB   = Ct_RD
      strip(iblade,lstrip)%CMOM4TB   = Cm
      strip(iblade,lstrip)%floss     = F
      strip(iblade,lstrip)%xi_sk     = xi_sk
      strip(iblade,lstrip)%gama_sk   = gama_sk
      strip(iblade,lstrip)%Coef_sk   = Coef_sk
      strip(iblade,lstrip)%UWAKIZTB  =-A *WEFF_RD(3)*(Coef_sk-1.d0)
      strip(iblade,lstrip)%Ui_ax     =-A *WEFF_RD(3)
      strip(iblade,lstrip)%Ui_per    = AP*WEFF_RD(1)
      strip(iblade,lstrip)%WEFF_RD(:)= WEFF_RD(:)
      strip(iblade,lstrip)%AINDTBmean= Amean

      strip(iblade,lstrip)%FCP (1:3) = FCP (1:3)
      strip(iblade,lstrip)%FCPL(1:3) = FCPL(1:3)
      strip(iblade,lstrip)%FCP (4:6) = MCP (1:3)
      strip(iblade,lstrip)%FCPL(4:6) = MCPL(1:3)
   enddo !iblade


 END Subroutine BEMT
!-----------------------------------------------------------
 Subroutine aero_loa_bem (iblade, lstrip)
!-----------------------------------------------------------

 Use Craft

   implicit none

   integer, intent(in ) :: iblade, lstrip

   real(8) :: chord , WEFF, PHI, Cl,Cd,Cm, Cn,Ct, CONST, AA(3,3)
   real(8) :: FCP(3), FCPL(3), MCP(3), MCPL(3)


   chord       = strip(iblade,lstrip)%chord
   PHI         = strip(iblade,lstrip)%PHITB
   Cl          = strip(iblade,lstrip)%CLIFTTB
   Cd          = strip(iblade,lstrip)%CDRAGTB
   Cm          = strip(iblade,lstrip)%CMOM4TB
   AA(1:3,1:3) = strip(iblade,lstrip)%A(1:3,1:3)
   WEFF        = dsqrt (strip(iblade,lstrip)%WEFFXTB**2 + &
                        strip(iblade,lstrip)%WEFFZTB**2    )

!--- Cn = normal, Ct = tangential at BLD
   Ct          = -( Cd*dcos(PHI) - Cl*dsin(PHI) )   !rotate by -PHI
   Cn          =    Cd*dsin(PHI) + Cl*dcos(PHI)
!--- Calculate loads
   CONST       = 0.5d0 * AIRDEN * WEFF**2 * chord
!--- The loads relative to the Blade local c.s.
   FCPL (:)    = 0.d0;
   FCPL (1)    =-CONST*Ct !tangent
   FCPL (3)    = CONST*Cn !normal
   MCPL (:)    = 0.d0;
   MCPL (2)    = CONST*Cm*chord  !--- Pitching Moment

!--- The loads on the strip relative to the rotor plane, before pitch rotation
   FCP         = matmul(AA,FCPL);
   MCP         = matmul(AA,MCPL);

   strip(iblade,lstrip)%FCP (1:3) = FCP (1:3)
   strip(iblade,lstrip)%FCPL(1:3) = FCPL(1:3)
   strip(iblade,lstrip)%FCP (4:6) = MCP (1:3)
   strip(iblade,lstrip)%FCPL(4:6) = MCPL(1:3)


 END Subroutine aero_loa_bem
!-----------------------------------------------------------
!
! ----Subroutine FATINT-------------------------------------
!
!     Calculates the function Fa(r/R) for the dynamic term
!
!-----------------------------------------------------------
 Subroutine FATINT (RRAT, Fa) !,Ft
!-----------------------------------------------------------

   implicit none

   real(8), intent(in ) :: RRAT      !r/R wrt Rotor Disk
   real(8), intent(out) :: Fa !,Ft

   real(8) :: pi, dphi,phi1,phi2,Fa1,Fa2
   integer :: N, j

      pi   = dacos(-1.d0)
      N    = 360
      Fa   = 0.0d0
      dphi = 2.d0*pi/dble(N)

   do j    = 1, N
      phi1 = dble(j-1)*dphi
      phi2 = dble(j  )*dphi
      Fa1  = (1.d0-RRAT*dcos(phi1))/(1.d0+RRAT**2-2.d0*RRAT*dcos(phi1))**1.5d0
      Fa2  = (1.d0-RRAT*dcos(phi2))/(1.d0+RRAT**2-2.d0*RRAT*dcos(phi2))**1.5d0
      Fa   = Fa + (Fa1+Fa2)*dphi/2.d0
   enddo
      Fa   = 2.d0*pi/Fa


 END Subroutine FATINT
!--------------------------------------------------------------------
 Subroutine WRITEOUTA
!--------------------------------------------------------------------

 Use Craft

   implicit none

   integer      :: nb  , lstrip, j
   real(8)      :: PSIB, WEFF
   character    :: CNUM1(2), CNUM2(3)
   character*80 :: outfil
 

# ifndef    ASCII
#    define ASCII 1
# endif

      PSIB   = (AZIM_ae+PHI0_ae(1))*R2D_ae
      PSIB   = dmod (PSIB,360.d0)
# if   ASCII == 0
   open  (1, file='LOADS_aer.bin',access='append',form='UNFORMATTED')
   write (1)                                                     &
# elif ASCII == 1
   open  (1, file='LOADS_aer.dat',access='append')
   write (1,100)                                                 &
# endif
     sngl( TIME_ae)                                             ,& !    1: Time   [sec]
     sngl( PSIB)                                                ,& !    2: Azimuth[deg]
     sngl( TTHRUST  /(0.5d0*AIRDEN*VELHUBT**2*PI_ae*RTCONE**2)) ,& !    3: Ct     [ - ]
!    sngl(-TTORQUE  /(0.5d0*AIRDEN*VELHUBT**2*PI_ae*RTCONE**3)) ,& !    4: Cm=Cp/λ[ - ] 
     sngl(-POWER    /(0.5d0*AIRDEN*VELHUBT**3*PI_ae*RTCONE**2)) ,& !    4: Cp     [ - ]
     sngl( OMEGA_ae*RTCONE/VELHUBT)                             ,& !    5: lamda  [ - ]
     sngl( TTHRUST  /1000.d0)                                   ,& !    6: Thrust [kN ]
     sngl(-TTORQUE  /1000.d0)                                   ,& !    7: Torque [kNm]
     sngl(-POWER    /1000.d0)                                   ,& !    8: Power  [kW ]
     sngl( POWER/TTORQUE*30.d0/PI_ae)                           ,& !    9: Omega  [rpm]
    (sngl( THRUST(j)/1000.d0),         j=1,NBLADE_ae)           ,& !10-13: Thrust [kN ] of each blade
    (sngl(-TORQUE(j)/1000.d0),         j=1,NBLADE_ae)              !14-16: Torque [kNm] of each blade
   close(1)


   do nb     = 1, NBLADE_ae
      call INT_2_CHAR ( 2, CNUM1, nb     )
   do lstrip = 1, NSTRIP
      call INT_2_CHAR ( 3, CNUM2, lstrip )

      PSIB   = (AZIM_ae+PHI0_ae(nb))*R2D_ae
      PSIB   = dmod (PSIB,360.d0)
      WEFF   = dsqrt(strip(nb,lstrip)%WEFFXTB**2+strip(nb,lstrip)%WEFFZTB**2)

#    if   ASCII == 0
      outfil = 'stripl'//CNUM1(1)//CNUM1(2)//'_'//CNUM2(1)//CNUM2(2)//CNUM2(3)//'.bin'
      open (2,file=outfil, access='append',form='UNFORMATTED')
      write(2)                                   &
#    elif ASCII == 1
      outfil = 'stripl'//CNUM1(1)//CNUM1(2)//'_'//CNUM2(1)//CNUM2(2)//CNUM2(3)//'.dat'
      open (2,file=outfil, access='append')
      write(2,100)                               &
#    endif
       sngl(TIME_ae)                            ,&  !    1: Time                  [sec]
       sngl(PSIB)                               ,&  !    2: Azimuth               [deg]
       sngl(strip(nb,lstrip)%rRD)               ,&  !    3: Radius_RD             [m  ]
       sngl(strip(nb,lstrip)%ALPHATB    *R2D_ae),&  !    4: AoA                   [deg]
       sngl(strip(nb,lstrip)%PHITB      *R2D_ae),&  !    5: Phi_BLD               [deg]
       sngl(WEFF)                               ,&  !    6: Weff_BLD              [m/s]
       sngl(strip(nb,lstrip)%AINDTB     *        &  !    7: a                     [ - ]
            strip(nb,lstrip)%Coef_sk    )       ,&
       sngl(strip(nb,lstrip)%AINDPTB    )       ,&  !    8: a'                    [ - ]
       sngl(strip(nb,lstrip)%CLIFTTB    )       ,&  !    9: Cl                    [ - ]
       sngl(strip(nb,lstrip)%CDRAGTB    )       ,&  !   10: Cd                    [ - ]
       sngl(strip(nb,lstrip)%CNORMTB    )       ,&  !   11: Cn                    [ - ]
       sngl(strip(nb,lstrip)%CTANGTB    )       ,&  !   12: Ct                    [ - ]
       sngl(strip(nb,lstrip)%CMOM4TB    )       ,&  !   13: Cm                    [ - ]
       sngl(strip(nb,lstrip)%FCPL(1))/1000.     ,&  !   14: Fx_BLD (Cyrcumfer.)   [kN/m]
       sngl(strip(nb,lstrip)%FCPL(3))/1000.     ,&  !   15: Fz_BLD (axial)        [kN/m]
       sngl(strip(nb,lstrip)%FCPL(5))/1000.     ,&  !   16: My_BLD                [kN ]
       sngl(strip(nb,lstrip)%twist      *R2D_ae),&  !   17: Twist                 [deg]
       sngl(strip(nb,lstrip)%pitch      *R2D_ae),&  !   18:    Pitch_BLD          [deg]
       sngl(strip(nb,lstrip)%THYTB      *R2D_ae),&  !   19:    Torsion_BLD        [deg]
       sngl(strip(nb,lstrip)%THY1TB     *R2D_ae),&  !   20:  d(Torsion+Pitch)_BLD [deg/s]
       sngl(strip(nb,lstrip)%THY2TB     *R2D_ae),&  !   21: dd(Torsion+Pitch)_BLD [deg/s^2]

       sngl(strip(nb,lstrip)%Mach       )       ,&  !   22: Mach                  [ - ]
       sngl(strip(nb,lstrip)%Reynolds   )       ,&  !   23: Reynolds              [ - ]
       sngl(strip(nb,lstrip)%floss      )       ,&  !   24: Loss Coeff            [ - ]
       sngl(strip(nb,lstrip)%gama_sk    *R2D_ae),&  !   25: Flow skewed angle     [deg]
       sngl(strip(nb,lstrip)%xi_sk      *R2D_ae),&  !   26: Wake skewed angle     [deg]
       sngl(strip(nb,lstrip)%Coef_sk    )       ,&  !   27: Yaw correct Coeff     [ - ]
       sngl(strip(nb,lstrip)%Ui_ax      )       ,&  !   28: Axial  ind. vel.      [m/s]
       sngl(strip(nb,lstrip)%Ui_per     )       ,&  !   29: Periph ind. vel.      [m/s]

      (sngl(strip(nb,lstrip)%UWINDGTB(j)),j=1,3),&  !30-32: Uwind_Glob            [m/s]
      (sngl(strip(nb,lstrip)%UBODY(j))   ,j=1,3),&  !33-35: Ubody_RD              [m/s]
       sngl(strip(nb,lstrip)%WEFFXTB    )       ,&  !   36: Weffx_BLD             [m/s]
       sngl(strip(nb,lstrip)%WEFFZTB    )       ,&  !   37: Weffz_BLD             [m/s]
       sngl(strip(nb,lstrip)%UWAKIZTB   )       ,&  !   38: Axial ind. vel dwe to skewed wake [m/s]

      (sngl(strip(nb,lstrip)%FCP(j))/1000.,j=1,6),& !39-44: Loads RD              [kN/m,kN]
      (sngl(strip(nb,lstrip)%WEFF_RD(j) ) ,j=1,3),& !45-47: WEFF_RD no induction  [m/s]

       sngl(strip(nb,lstrip)%AINDTBmean )           !   48: Mean a for dyn. wake  [ - ]
      close (2)
      enddo
   enddo

   call INFLOW_out ( TIME_ae )

 100  format (150f15.5)


 END Subroutine WRITEOUTA
!------------------------------------------------------------------------
!---- Subroutine : UINFLOW_bem
!
!--   Determine the wind inflow velocity for each aero-strip of all the
!     blades wrt the deformed rotor plane (before pitch) c.s.
!
!     1. Turbulence is also added
!     2. Wind scenarios are set here
!     3. HUB_VEL(1:3), BLD_VEL_loc(1:3), BLD_VEL_glob(1:3) are set also [outputs]
!
!     The wind will be "ok" if the cone and the pre-sweep/bend angles are zero.
!     Otherwise the wind has to be rotated in order to take into account the
!     two angles, because atm only the pitch and the twist rotations are performed.
!
!     WINDYAW0      :initial mean wind yaw [it might chage based on IWINDC]
!     WINDTAW       :actual  mean wind yaw [mean: at hub height because VEER=0]
!     windyaw_strip :actual  mean wind yaw including VEER
!
!------------------------------------------------------------------------
 Subroutine UINFLOW_bem
!----------------------------------------------------------------------

 Use Craft

   implicit none

   real(8) :: RG_tow_in(3), AT_tow_in(3,3), DRG_tow_in(3)
   integer :: iblade, lstrip
   real(8) :: XG(3), UG(3), UL(3), A(3,3)
!!!real(8) :: V_RD(3), THY, THY1, THY2, THYP


   call get_tower_ref  ( RG_tow_in, AT_tow_in, DRG_tow_in )

   call INFLOW (TIME_ae, RG_tow_in, AT_tow_in, DRG_tow_in,&      !inputs
                HUB_VEL, VELHUBT  , WINDYAW  , WINDINC      )    !outputs [Craft module]


!--- Loop for all blades
   do iblade = 1, NBLADE_ae
!------ For all strips
      do lstrip = 1, NSTRIP

!!!      call Elast2Aero_bem1 ( iblade, lstrip, XG, V_RD, THY, THY1, THY2, THYP, pitch_ae )

!--------- XG: Control point coordinates with respect to the absolute frame (used for shear exp)
         A  (1:3,1:3) = matmul( matmul(AP_ae(iblade,:,:),A_cone_ae(:,:)),  A_pit_ae(iblade,lstrip,:,:) )
         XG (1:3) = HUBPOS_ae(1:3) + matmul (  A(1:3,1:3),     (node(iblade,lstrip+1)%xBld(1:3)+node(iblade,lstrip)%xBld(1:3))/2.d0  )
         call UINFLOW (TIME_ae, XG, UG, 1, 1) !include_turb, include_tshadow

!--------- limit4
         if (IPARKED==0.and.INOINDUCT==0)&
         UG(1)=max(UG(1),1.0d0)

!--------- Transform the Velocity wrt the deformed rotor plane (before Cone, pitch) c.s.
         UL(1:3) = matmul ( ATP_ae(iblade,1:3,1:3),UG(1:3) )

         strip(iblade,lstrip)%UWINDTB (1:3) = UL(1:3)
         strip(iblade,lstrip)%UWINDGTB(1:3) = UG(1:3)
      enddo !lstrip
   enddo !iblade


 END Subroutine UINFLOW_bem
!
!
!
!---------------------------------------------------------------------------
!   Subroutine : LOCAL_AERO_PARAM
!
!   Set local Mach and incidence dependent parameters
!   We consider the A0 for every flap angle, but the CLslope at zero flap angle.
!
!---------------------------------------------------------------------------
 Subroutine LOCAL_AERO_PARAM ( NB, HLOC, AEFF, AMACH, lstrip )
!--------------------------------------------------------------------------

 Use Craft
 Use Caero
 Use Foilfs_mod

   implicit none

   real(8) :: AAT0  (NANGM ), AAT1  (NANGM )
   real(8) :: CLT0  (NANGM ), CLT1  (NANGM )
   real(8) :: CDT0  (NANGM ), CDT1  (NANGM )
   real(8) :: CMT0  (NANGM ), CMT1  (NANGM )
   real(8) :: HLOC, AEFF, AMACH, HLOC0, A1, A2
   real(8) :: CL1,CD1,CM1, CL2,CD2,CM2, AZER, AZER0
   real(8) :: CL0,CD0,CM0, CLSLO, A0, A0M,A0P, CLSLO0
   real(8) :: CL0M,CD0M,CM0M,CL0P,CD0P,CM0P
   real(8) :: CMSLO1, CMSLO2,CMSLO, CMO, CDSLO,CDO
   real(8) :: ACDMIN, AEFF0,AEFFM,AEFFP,CL,CD,CM
   real(8) :: CLP,CDP,CMP,CLM,CDM,CMM, CLLIN,CDLIN,CMLIN
   real(8) :: DCL,DCD,DCM,BETA2,BETA0,PAR1,PAR2
   integer :: NB, i,imin, lstrip


   HLOC0 = HLOC

!--- Interpolation for Flap angles
   call  CLCDCM_CURV (NB, HLOC0, FLCONTR(NB, lstrip),IspanANGmax, &
                      RCB2   , FLAPANGS, AATB, CLTB, CDTB, CMTB,  &
                      MBLADE , NSPANB2 , NSPANM   ,               &
                      NANGB  , NANGM   , NFLAPANGS,               &
                      CLT0 , CDT0 , CMT0 , AAT0 , AngFLAPSmax)

!--- Interpolation for zero Flap angle
   call  CLCDCM_CURV (NB, HLOC0, 0.d0               ,IspanANGmax, &
                      RCB2   , FLAPANGS, AATB, CLTB, CDTB, CMTB,  &
                      MBLADE , NSPANB2 , NSPANM   ,               &
                      NANGB  , NANGM   , NFLAPANGS,               &
                      CLT1 , CDT1 , CMT1 , AAT1 , AngFLAPSmax)

   NPP   = NANGB(NB,IspanANGmax)

   do i=1,NPP
      AGON(i) = AAT0(i)
      CLST(i) = CLT0(i)
      CDST(i) = CDT0(i)
      CMST(i) = CMT0(i)
   enddo

   A1    = 1.d0 !  0.d0  (BO105)
   A2    =-4.d0 ! -2.d0

   call FOIL_CLCDCM (A1   ,NPP,NANGM,AAT0,CLT0,CDT0,CMT0,CL1,CD1,CM1)
   call FOIL_CLCDCM (A2   ,NPP,NANGM,AAT0,CLT0,CDT0,CMT0,CL2,CD2,CM2)

   A1    = A1/R2D_ae
   A2    = A2/R2D_ae

!--- Zero lift angle [AZER] at current flap angle
   AZER  =  A1 - (A1-A2)/(CL1 -CL2 )*CL1

   A1    =  1.d0 !  0.d0  (BO105)
   A2    = -4.d0 ! -2.d0

   call FOIL_CLCDCM (A1   ,NPP,NANGM,AAT1,CLT1,CDT1,CMT1,CL1,CD1,CM1)
   call FOIL_CLCDCM (A2   ,NPP,NANGM,AAT1,CLT1,CDT1,CMT1,CL2,CD2,CM2)

   A1    = A1/R2D_ae
   A2    = A2/R2D_ae

!--- Zero lift angle [AZER0] at zero flap angle
   AZER0 =  A1 - (A1-A2)/(CL1 -CL2)*CL1

   A1    =  0.d0 ! AZER*R2D_ae - 1.0
   A2    =  4.d0 ! AZER*R2D_ae + 1.0

   call FOIL_CLCDCM (A1   ,NPP,NANGM,AAT0,CLT0,CDT0,CMT0,CL1,CD1,CM1)
   call FOIL_CLCDCM (A2   ,NPP,NANGM,AAT0,CLT0,CDT0,CMT0,CL2,CD2,CM2)

   A1    = A1/R2D_ae
   A2    = A2/R2D_ae

!--- Slope of Cl around AZER [CLSLO]
   CLSLO = (CL2-CL1)/(A2-A1)

   A1    =  0.d0 ! AZER*180./PI - 1.0 !  -2.d0  (BO105)
   A2    =  4.d0 ! AZER*180./PI + 1.0 !   0.d0

   call FOIL_CLCDCM (A1   ,NPP,NANGM,AAT1,CLT1,CDT1,CMT1,CL1,CD1,CM1)
   call FOIL_CLCDCM (A2   ,NPP,NANGM,AAT1,CLT1,CDT1,CMT1,CL2,CD2,CM2)

   A1    = A1/R2D_ae
   A2    = A2/R2D_ae
   
!--- Slope of Cl around AZER0 [CLSLO0]
   CLSLO0= (CL2-CL1)/(A2-A1)

   A0    = 0.d0
   A0M   =-0.1d0
   A0P   = 0.1d0

   call FOIL_CLCDCM (A0   ,NPP,NANGM,AAT0,CLT0,CDT0,CMT0,CL0 ,CD0 ,CM0 )
   call FOIL_CLCDCM (A0M  ,NPP,NANGM,AAT0,CLT0,CDT0,CMT0,CL0M,CD0M,CM0M)
   call FOIL_CLCDCM (A0P  ,NPP,NANGM,AAT0,CLT0,CDT0,CMT0,CL0P,CD0P,CM0P)

   A0    = A0 / R2D_ae
   A0M   = A0M/ R2D_ae
   A0P   = A0P/ R2D_ae

!--- CM at aoa 0deg and corresponding slope [CMO, CMSLO]
   CMSLO1= (CM0P-CM0 )/(A0P-A0 )
   CMSLO2= (CM0 -CM0M)/(A0 -A0M)
   CMSLO = 0.5d0*(CMSLO1+CMSLO2)
   CMO   = CM0

   CD0 = 1000.d0
   do i = 1, NPP
      if (CD0 > CDT0(i)) then
          CD0 = CDT0(i)
          imin=i
      endif
   enddo

!--- minimum CD and corresponding slope at aoa [CDO, CDSLO, ACDMIN]
   CDSLO      = 0.d0
   CDO        = CD0
   ACDMIN     = AAT0(imin)/ R2D_ae

   AEFF0      = AEFF * R2D_ae
   AEFFP      = AEFF * R2D_ae + 0.05d0
   AEFFM      = AEFF * R2D_ae - 0.05d0

   call FOIL_CLCDCM (AEFF0,NPP,NANGM,AAT0,CLT0,CDT0,CMT0,CL ,CD ,CM )
   call FOIL_CLCDCM (AEFFP,NPP,NANGM,AAT0,CLT0,CDT0,CMT0,CLP,CDP,CMP)
   call FOIL_CLCDCM (AEFFM,NPP,NANGM,AAT0,CLT0,CDT0,CMT0,CLM,CDM,CMM)


   AEFFP      = AEFF + 0.05d0 / R2D_ae
   AEFFM      = AEFF - 0.05d0 / R2D_ae


  if ( IFLAPLOC(NB, lstrip) ==  1 ) then 
   CLLIN      =        CLSLO0*           (AEFF - AZER)
  else
   CLLIN      =  0.5d0*CLSLO0*dsin(2.0d0*(AEFF - AZER))
  endif
   CDLIN      =  CDO
   CMLIN      =  CMO + CMSLO*AEFF

   DCL        =  CLLIN - CL
   DCD        =  CDLIN - CD
   DCM        =  CMLIN - CM

   AZER_ae    = AZER
   AZER0_ae   = AZER0

   DCLLIN_ae  = CLSLO0
   DCMLIN_ae  = CMSLO

   CDLIN_ae   = CDLIN
   CMLIN_ae   = CMLIN
   CMO_ae     = CMO

   DCL_ae     = DCL
   DCD_ae     = DCD
   DCM_ae     = DCM

   DACL_ae    = (CLP-CLM)/(AEFFP-AEFFM)
   DACD_ae    = (CDP-CDM)/(AEFFP-AEFFM)
   DACM_ae    = (CMP-CMM)/(AEFFP-AEFFM)

   CL_ae      = CL
   CD_ae      = CD
   CM_ae      = CM


!--- ONERA variables (AMACH, DCL)
   BETA2      = 1.d0-AMACH**2
   BETA0      = dsqrt(BETA2)

   AS_L_ae    =       PI_ae+ 5.d0  *PI_ae*(BETA2**0.285d0-1d0)
   AK_L_ae    = 0.5d0*PI_ae+ 1.96d0*PI_ae*(BETA0         -1d0)
   AL_L_ae    = 0.17d0     - 0.13d0* AMACH
   AAL_L_ae   = 0.53d0     + 0.25d0*(BETA0  -1d0)
   ASIG_L_ae  = 2.d0*PI_ae/BETA0
   AD_L_ae    = 0.d0*dabs(DCL) !check is ZERO
   AA0_L_ae   = 0.30d0
   AR0_L_ae   = 0.18d0
   AA2_L_ae   = 0.20d0
   AR2_L_ae   = 0.18d0
   AE2_L_ae   =-1.50d0

   ASIG_D_ae  = 0.d0
   AA0_D_ae   = 0.25d0
   AR0_D_ae   = 0.22d0
   AA2_D_ae   = 0.d0
   AR2_D_ae   = 0.20d0
   AE2_D_ae   =-1.15d0

   PAR1       = -3.d0*PI_ae/16.d0*(-1.26d0-1.53d0*datan(15d0*(AMACH-0.7d0)))
   PAR2       = -     PI_ae/2.d0 * (1.d0+1.4d0*AMACH**2)

   AS_M_ae    = PAR1
   ASIGP_M_ae = 0.5d0*PAR2
   ASIG_M_ae  = PAR2-PAR1
   AD_M_ae    = 0.d0*dabs(DCL) !check is ZERO
   AA0_M_ae   = 0.25d0
   AR0_M_ae   = 0.22d0
   AA2_M_ae   = 0.05d0
   AR2_M_ae   = 0.20d0
   AE2_M_ae   = 1.425d0 ! CVas1 limit value 0.50


 END Subroutine LOCAL_AERO_PARAM
!------------------------------------------------------------------
 Subroutine CLCDCM_CURV (nb, y, AFL,  Ismax            ,      &
                         RC, FLAPANGS  , AAT, CLT, CDT, CMT,  &
                         MBLADE, NSPAN, NSPANM,               &
                         NANG,NANGM,NFLAPANGS ,               &
                         CLT0,CDT0,CMT0,AAT0, AngFLAPSmax)
!------------------------------------------------------------------

   implicit none

!! integer, intent(in ) :: nb, Ismax, MBLADE, NSPANM, NMACHM, NANGM
!! real(8), intent(in ) :: y, AMACH
!! real(8), intent(in ) :: AAT   (MBLADE,NSPANM,       NANGM)
!! real(8), intent(in ) :: CLT   (MBLADE,NSPANM,NMACHM,NANGM)
!! real(8), intent(in ) :: CDT   (MBLADE,NSPANM,NMACHM,NANGM)
!! real(8), intent(in ) :: CMT   (MBLADE,NSPANM,NMACHM,NANGM)
!! real(8), intent(out) :: AAT0  (                     NANGM)
!! real(8), intent(out) :: CLT0  (                     NANGM)
!! real(8), intent(out) :: CDT0  (                     NANGM)
!! real(8), intent(out) :: CMT0  (                     NANGM)

!! real(8) :: CL1   (                     NANGM)
!! real(8) :: CL2   (                     NANGM)
!! real(8) :: CD1   (                     NANGM)
!! real(8) :: CD2   (                     NANGM)
!! real(8) :: CM1   (                     NANGM)
!! real(8) :: CM2   (                     NANGM)
!! real(8) :: AMACHT(MBLADE,       NMACHM      )
!! real(8) :: RC    (MBLADE,       NSPANM      )
!! integer :: NANG  (MBLADE,       NSPANM      )
!! integer :: NMACH (MBLADE                    )
!! integer :: NSPAN (MBLADE                    )
!! integer :: i, ispan, imach, ii,jj, iang,iang1
!! real(8) :: CL01,CD01,CM01,CL02,CD02,CM02,A, HTA



   integer :: iflap, nb, Ismax, MBLADE
   integer :: NSPANM, NANGM, AngFLAPSmax, i, ispan
   integer :: ii,jj, iang, iang1
   real(8) :: AFL, y, A, CL01, CL02, CD01, CD02, CM01, CM02,  HTA
   real(8) :: AAT (MBLADE,NSPANM,            NANGM)
   real(8) :: CLT (MBLADE,NSPANM,AngFLAPSmax,NANGM)
   real(8) :: CDT (MBLADE,NSPANM,AngFLAPSmax,NANGM)
   real(8) :: CMT (MBLADE,NSPANM,AngFLAPSmax,NANGM)
   real(8) :: AAT0(NANGM)
   real(8) :: CLT0(NANGM)
   real(8) :: CDT0(NANGM)
   real(8) :: CMT0(NANGM)
   real(8) :: CL1 (NANGM), CL2(NANGM)
   real(8) :: CD1 (NANGM), CD2(NANGM)
   real(8) :: CM1 (NANGM), CM2(NANGM)
   real(8) :: FLAPANGS (MBLADE, NSPANM   ,AngFLAPSmax)
   integer :: NANG     (MBLADE, NSPANM)
   integer :: NSPAN    (MBLADE)
   integer :: NFLAPANGS(MBLADE, NSPANM)
   real (8):: RC       (MBLADE, NSPANM)


!--- Find span position 
   do i=1,NSPAN(nb)-1
      if (y.ge.RC(nb,i).and.y.le.RC(nb,i+1)) then
         ispan = i
         goto 11
      endif
   enddo !i

   write(* ,*) 'ispan not found',nb,y,AFL
   write(10,*) 'ispan not found',nb,y,AFL
   stop
11 continue

!--- Find flap interval
     do i=1,NFLAPANGS(nb,ispan)-1
         if ((AFL.ge.FLAPANGS(nb,ispan,i)).and.(AFL.le.FLAPANGS(nb,ispan,i+1) ) ) then
            iflap = i
            goto 12
         endif
     enddo !i

   write(* ,*) 'iflap not found',nb,y,AFL
   write(10,*) 'iflap not found',nb,y,AFL
   stop
12 continue

   ii = ispan
   jj = iflap

!--- Interpolation accoring to flap angle START
!--- CL,CD,CM  for 1st position (ispan) for current FLAP by means of linear interpolation
!   for all given angles(ispan)
   do iang=1,NANG(nb,ii)            ! All angles for body nb at span ii (first)

      CL1(iang) =   CLT(nb,ii,jj,iang)                       &
                 + (CLT(nb,ii,jj+1,iang)-CLT(nb,ii,jj,iang)) &
      /(FLAPANGS(nb,ii,jj+1)-FLAPANGS(nb,ii,jj))*(AFL-FLAPANGS(nb,ii,jj))
      CD1(iang) =   CDT(nb,ii,jj,iang)                       &
                 + (CDT(nb,ii,jj+1,iang)-CDT(nb,ii,jj,iang)) &
      /(FLAPANGS(nb,ii,jj+1)-FLAPANGS(nb,ii,jj))*(AFL-FLAPANGS(nb,ii,jj))
      CM1(iang) =   CMT(nb,ii,jj,iang)                       &
                 + (CMT(nb,ii,jj+1,iang)-CMT(nb,ii,jj,iang)) &
      /(FLAPANGS(nb,ii,jj+1)-FLAPANGS(nb,ii,jj))*(AFL-FLAPANGS(nb,ii,jj))
   enddo !iang

!--- CL,CD,CM  for 2nd position (ispan) for current FLAP by means of linear interpolation
!   for all given angles
   do iang=1,NANG(nb,ii+1)           ! All angles for body nb at span ii+1 (second)

      CL2(iang) =  CLT(nb,ii+1,jj,iang)                         &
                + (CLT(nb,ii+1,jj+1,iang)-CLT(nb,ii+1,jj,iang)) &
      /(FLAPANGS(nb,ii+1,jj+1)-FLAPANGS(nb,ii+1,jj))*(AFL-FLAPANGS(nb,ii+1,jj))
      CD2(iang) =  CDT(nb,ii+1,jj,iang)                         &
                + (CDT(nb,ii+1,jj+1,iang)-CDT(nb,ii+1,jj,iang)) &
      /(FLAPANGS(nb,ii+1,jj+1)-FLAPANGS(nb,ii+1,jj))*(AFL-FLAPANGS(nb,ii+1,jj))
      CM2(iang) =  CMT(nb,ii+1,jj,iang)                         &
                + (CMT(nb,ii+1,jj+1,iang)-CMT(nb,ii+1,jj,iang)) &
      /(FLAPANGS(nb,ii+1,jj+1)-FLAPANGS(nb,ii+1,jj))*(AFL-FLAPANGS(nb,ii+1,jj))
   enddo !iang

   do iang=1,NANG(nb,Ismax)          ! Do loop with span with most angles

      A = AAT(nb,Ismax,iang)
!--- Interpolation according to flap angle END

!--- Angle interpolation of the span

!--- CL,CD,CM  for 1st position (ispan) for current FLAP by means of linear interpolation 
!   for all given angles(ispan) of the 1st ispan
      do iang1=1,NANG(nb,ii)-1

         if(A.ge.AAT(nb,ii,iang1).and.A.le.AAT(nb,ii,iang1+1)) then

            CL01 = CL1(iang1)+(CL1(iang1+1)-CL1(iang1))/             &
                              (AAT(nb,ii,iang1+1)-AAT(nb,ii,iang1))  &
                             *(A-AAT(nb,ii,iang1))
            CD01 = CD1(iang1)+(CD1(iang1+1)-CD1(iang1))/             &
                              (AAT(nb,ii,iang1+1)-AAT(nb,ii,iang1))  &
                             *(A-AAT(nb,ii,iang1))
            CM01 = CM1(iang1)+(CM1(iang1+1)-CM1(iang1))/             &
                              (AAT(nb,ii,iang1+1)-AAT(nb,ii,iang1))  &
                             *(A-AAT(nb,ii,iang1))
            goto 13

         endif

      enddo!iang1

      write(* ,*) 'C[L,D,M]01 not found',iang,A
      write(10,*) 'C[L,D,M]01 not found',iang,A
      stop
13    continue

!--- CL,CD,CM  for 2nd position (ispan) for current Mach by means of linear interpolation
!   for all given angles(ispan) of the 1st ispan
      do iang1=1,NANG(nb,ii+1)-1

        if(A.ge.AAT(nb,ii+1,iang1).and.A.le.AAT(nb,ii+1,iang1+1)) then

          CL02 = CL2(iang1)+(CL2(iang1+1)-CL2(iang1))/                 &
                            (AAT(nb,ii+1,iang1+1)-AAT(nb,ii+1,iang1))  &
                           *(A-AAT(nb,ii+1,iang1))
          CD02 = CD2(iang1)+(CD2(iang1+1)-CD2(iang1))/                 &
                            (AAT(nb,ii+1,iang1+1)-AAT(nb,ii+1,iang1))  &
                           *(A-AAT(nb,ii+1,iang1))
          CM02 = CM2(iang1)+(CM2(iang1+1)-CM2(iang1))/                 &
                            (AAT(nb,ii+1,iang1+1)-AAT(nb,ii+1,iang1))  &
                           *(A-AAT(nb,ii+1,iang1))
          goto 14

        endif

      enddo!iang1

      write(* ,*) 'C[L,D,M]02 not found',iang,A
      write(10,*) 'C[L,D,M]02 not found',iang,A
      stop
14    continue

!-- Outouts the CL,CD,CM  for all the angles of the 1st ispan for the current Y and Flap Angle
      HTA         = (y-RC(nb,ii))/(RC(nb,ii+1)-RC(nb,ii))
      CLT0 (iang) = (1.d0-HTA)*CL01 + HTA*CL02
      CDT0 (iang) = (1.d0-HTA)*CD01 + HTA*CD02
      CMT0 (iang) = (1.d0-HTA)*CM01 + HTA*CM02
      AAT0 (iang) =  A
   enddo !iang


 END Subroutine CLCDCM_CURV
!----------------------------------------------------------------------
 Subroutine FOIL_CLCDCM (ALPHA,NPANG,NPANGM,ANGT,CLT,CDT,CMT,CL,CD,CM)
!----------------------------------------------------------------------

   implicit none

   integer, intent(in ) :: NPANG,NPANGM
   real(8), intent(in ) :: ALPHA
   real(8), intent(in ) :: ANGT(NPANGM),CLT(NPANGM),CDT(NPANGM), CMT(NPANGM)
   real(8), intent(out) :: CL, CD, CM

   real(8) :: A, PER
   integer :: I

   A = ALPHA

   if (A >  180.d0) A = -360.d0 + A
   if (A < -180.d0) A =  360.d0 + A

   if ( (A > 180.d0).or.(A <-180.d0) ) then
      write(10,*) 'error in FOIL_CLCDCM, A out of range!',A
      stop
   endif

   do I = 1, NPANG-1
      if ((A >= ANGT(I)).and.(A < ANGT(I+1))) then

         PER = (A-ANGT(I))/(ANGT(I+1)-ANGT(I))
         CL  = PER*(CLT(I+1)-CLT(I))+CLT(I)
         CD  = PER*(CDT(I+1)-CDT(I))+CDT(I)
         CM  = PER*(CMT(I+1)-CMT(I))+CMT(I)

         return
      endif
   enddo

   write(10,*) 'Angle not found in FOIL_CLCDCM',A; stop


 END Subroutine FOIL_CLCDCM
!--------------------------------------------------------------------------
!
! --  Subroutine :ONERALDM
! --  ONERA stall model for the lift drag and moment ----------------------
!
!--------------------------------------------------------------------------
 Subroutine ONERALDM ( RLD, CHD, WEFF, ALPHA, THY, THY1, THY2, &
                       CLIFTDS,CDRAGDS,CMOMDS,istrip,lstrip,iblade,NTIME )
!--------------------------------------------------------------------------

 Use Craft
 Use Caero
 Use Foilfs_mod

   implicit none

   real(8) :: RLD, CHD, WEFF, ALPHA, THY, THY1, THY2
   real(8) :: CLIFTDS,CDRAGDS,CMOMDS
   integer :: istrip,NTIME, lstrip,iblade
!--- local
   real(8) :: ACHORD,AEFF,AEFFM,ARSPAN,CLLIN_ae,DCL3D,DCD3D,DCLIFT,DCDRAG
   real(8) :: ALFA1,ALFA2,ALFA1M,ALFA2M,CPAR,COEF,COEFP,COEF1,COEFP1,WEFFT,AEFFT
   real(8) :: THYT, THYTT,W0T,W1,W1T,TEFF,QTOP1,QTOP2,QTOP3,ABIG
   real(8) :: A11,A22,A33,CON,CON1,CON2,CON3
   real(8) :: AA_L_ae, AR_L_ae, AR_L0_ae, AE_L_ae, AL_L0_ae
   real(8) :: AA_D_ae, AR_D_ae, AR_D0_ae, AE_D_ae
   real(8) :: AA_M_ae, AR_M_ae, AR_M0_ae, AE_M_ae
   real(8) :: G1L,G2L,G2D,G2M,ASIG_D
   real(8) :: CLLIN_pot


   ACHORD = CHD
   ARSPAN = RLD

   AEFF   = ALPHA
   AEFFM  = AEFF*180.d0/PI_ae

   ANGEFF (istrip) = AEFF
   VELEFF (istrip) = WEFF

!--- Definition of the DCL = CLlinear - CL steady  ------------------------------
  if ( IFLAPLOC(iblade, lstrip) == 1 ) then
   CLLIN_ae  = DCLLIN_ae * (AEFF-AZER_ae)
   CLLIN_pot = 2*PI_ae*(AEFF-AZER_ae0(iblade,lstrip))  
  else
   CLLIN_ae = 0.5d0*DCLLIN_ae*dsin(2.0d0*(AEFF-AZER_ae))
  endif


!--- 3D effects
   DCL3D  =  CLLIN_ae - CL_ae
   DCD3D  =-(CDLIN_ae - CD_ae)

   if (AEFFM.gt.25.d0.and.AEFFM.lt.40d0) then
    DCL3D = (1d0-(AEFFM-25d0)/15d0)*DCL3D
    DCD3D = (1d0-(AEFFM-25d0)/15d0)*DCD3D
   else if (AEFFM.ge.40d0) then
    DCL3D = 0.d0
    DCD3D = 0.d0
   endif

   DCLIFT = 2.2d0*(ACHORD/ARSPAN)**1.5d0*(dcos(THY))**4*DCL3D
   DCDRAG = 2.2d0*(ACHORD/ARSPAN)**1.5d0*(dcos(THY))**4*DCD3D

!--- 3D effects not active
   DCLIFT = 0.d0
   DCDRAG = 0.d0

   CL_ae = CL_ae + DCLIFT
   CD_ae = CD_ae + DCDRAG

!--- Parameters of the equation
   if     (NTIME.le. NTIME_un  ) then

     CLIFTDS = CL_ae
     CDRAGDS = CD_ae
     CMOMDS  = CM_ae
     return

   elseif (NTIME.eq. NTIME_un+1) then

     ANGEFFP (3,istrip) = ANGEFF (istrip)
     VELEFFP (3,istrip) = VELEFF (istrip)
     CLOLDA  (3,istrip) = CLLIN_ae*WEFF
     CLOLDB  (3,istrip) = 0.d0
     CDOLDB  (3,istrip) = 0.d0
     CMOLDB  (3,istrip) = 0.d0

     CLIFTDS = CL_ae
     CDRAGDS = CD_ae
     CMOMDS  = CM_ae
     return

   elseif (NTIME.eq. NTIME_un+2) then

     ANGEFFP (2,istrip) = ANGEFF (istrip)
     VELEFFP (2,istrip) = VELEFF (istrip)
     CLOLDA  (2,istrip) = CLLIN_ae*WEFF
     CLOLDB  (2,istrip) = 0.d0
     CDOLDB  (2,istrip) = 0.d0
     CMOLDB  (2,istrip) = 0.d0

     CLIFTDS = CL_ae
     CDRAGDS = CD_ae
     CMOMDS  = CM_ae
     return

   elseif (NTIME.eq. NTIME_un+3) then

     ANGEFFP (1,istrip) = ANGEFF (istrip)
     VELEFFP (1,istrip) = VELEFF (istrip)
     CLOLDA  (1,istrip) = CLLIN_ae*WEFF
     CLOLDB  (1,istrip) = 0.d0
     CDOLDB  (1,istrip) = 0.d0
     CMOLDB  (1,istrip) = 0.d0

     CLIFTDS = CL_ae
     CDRAGDS = CD_ae
     CMOMDS  = CM_ae
     return

   endif


!--- For high values of the incidence
   ALFA1  =  60.d0
   ALFA2  =  80.d0
   ALFA1m = -60.d0
   ALFA2m = -80.d0
   ABIG   = 100.d0

   if      ((AEFFM.gt.ALFA1).and.(AEFFM.lt.ALFA2)) then    !a1<a<a2
     CPAR  = (AEFFM-ALFA1)/(ALFA2-ALFA1)
     COEF  =  1.d0*(1.d0-CPAR) + ABIG*CPAR
     COEFP =  1.d0*(1.d0-CPAR) + 5.d0*CPAR
   else if (AEFFM.gt.ALFA2) then                           !   a>a2
     COEF  = ABIG
     COEFP = 5.d0
   else if ((AEFFM.lt.ALFA1m).and.(AEFFM.gt.ALFA2m)) then  !   a1m<a<a2m
     CPAR  = (AEFFM-ALFA1m)/(ALFA2m-ALFA1m)
     COEF  =  1.d0*(1.d0-CPAR) + ABIG*CPAR
     COEFP =  1.d0*(1.d0-CPAR) + 5.d0*CPAR
   else if (AEFFM.lt.ALFA2m) then                          !       a<a2m
     COEF  = ABIG
     COEFP = 5.d0
   else                                                    !  a1m<a<a1
     COEF  =  1.d0
     COEFP =  1.d0
   endif
   COEF1 = 1.d0/COEF    !0.01
   COEFP1= 1.d0/COEFP   !0.2


!--- Second order derivatives of W0 and THY
   WEFFT= ( 3.d0*VELEFF (  istrip) - 4.d0*VELEFFP(1,istrip) + VELEFFP(2,istrip) )/(2.d0*DT_ae)
   AEFFT= ( 3.d0*ANGEFF (  istrip) - 4.d0*ANGEFFP(1,istrip) + ANGEFFP(2,istrip) )/(2.d0*DT_ae)

   THYT =-THY1
   THYTT=-THY2

   W0T  = (WEFFT*dsin(AEFF) + WEFF*dcos(AEFF)*AEFFT) *COEF1
   W1   = 0.5d0*ACHORD* THYT                         *COEF1
   W1T  = 0.5d0*ACHORD* THYTT                        *COEF1
   TEFF = 0.5d0*ACHORD/WEFF


!--- Lift Equation for the attached flow
   QTOP1    = CLOLDA(1,istrip)
   QTOP2    = CLOLDA(2,istrip)

   AL_L0_ae = AL_L_ae*COEFP

   A22      =  -AL_L0_ae/TEFF
   A33      =   0.5d0*AL_L0_ae/TEFF*WEFF*DCLLIN_ae*dsin(2.0d0*(AEFF-AZER_ae)) &
               +      AL_L_ae /TEFF*ASIG_L_ae    *W1                          &
               +     (AAL_L_ae*DCLLIN_ae+AD_L_ae)*W0T                         &
               +      AAL_L_ae*ASIG_L_ae         *W1T

   CON      = 1.d0-2.d0/3.d0*DT_ae*A22
   CON1     =      2.d0/3.d0*DT_ae*A33
   CON2     =      4.d0/3.d0
   CON3     =    - 1.d0/3.d0

   CLIFTA(1,istrip)= 1.d0/CON*(CON1+CON2*QTOP1+CON3*QTOP2)


!--- Lift Equation for the separated flow
   QTOP1    = CLOLDB(1,istrip)
   QTOP2    = CLOLDB(2,istrip)
   QTOP3    = CLOLDB(3,istrip)

   AA_L_ae  =  AA0_L_ae  +  AA2_L_ae*DCL_ae**2
   AR_L_ae  = (AR0_L_ae  +  AR2_L_ae*DCL_ae**2)**2
   AE_L_ae  =               AE2_L_ae*DCL_ae**2

   AR_L0_ae =  AR_L_ae*COEFP

   A11      = - AA_L_ae /TEFF
   A22      = - AR_L0_ae/TEFF**2
   A33      = - AR_L0_ae/TEFF**2*WEFF*DCL_ae - AE_L_ae /TEFF   *W0T

   CON      = 2.d0/DT_ae**2-3.d0*A11/(2.d0*DT_ae)-A22
   CON1     = 5.d0/DT_ae**2-2.d0*A11/DT_ae
   CON2     = A11/(2.d0*DT_ae)-4.d0/(DT_ae**2)
   CON3     = 1.d0/DT_ae**2

   CLIFTB(1,istrip) =  CON1/CON*QTOP1 + CON2/CON*QTOP2 + CON3/CON*QTOP3 + A33/CON


!--- Drag Equation for the separated flow
   QTOP1    = CDOLDB(1,istrip)
   QTOP2    = CDOLDB(2,istrip)
   QTOP3    = CDOLDB(3,istrip)

   AA_D_ae  =  AA0_D_ae  +  AA2_D_ae*DCL_ae**2
   AR_D_ae  = (AR0_D_ae  +  AR2_D_ae*DCL_ae**2)**2
   AE_D_ae  =               AE2_D_ae*DCL_ae**2

   AR_D0_ae =  AR_D_ae*COEFP

   A11      = - AA_D_ae /TEFF
   A22      = - AR_D0_ae/TEFF**2
   A33      = - AR_D0_ae/TEFF**2*WEFF*DCD_ae - AE_D_ae /TEFF   *W0T

   CON      = 2.d0/DT_ae**2-3.d0*A11/(2.d0*DT_ae)-A22
   CON1     = 5.d0/DT_ae**2-2.d0*A11/DT_ae
   CON2     = A11/(2.d0*DT_ae)-4.d0/(DT_ae**2)
   CON3     = 1.d0/DT_ae**2

   CDRAGB(1,istrip) =  CON1/CON*QTOP1 + CON2/CON*QTOP2 + CON3/CON*QTOP3 + A33/CON


!--- Moment Equation for the separated flow
   QTOP1    = CMOLDB(1,istrip)
   QTOP2    = CMOLDB(2,istrip)
   QTOP3    = CMOLDB(3,istrip)

   AA_M_ae  =  AA0_M_ae  +  AA2_M_ae*DCL_ae**2
   AR_M_ae  = (AR0_M_ae  +  AR2_M_ae*DCL_ae**2)**2
   AE_M_ae  =               AE2_M_ae*DCL_ae**2

   AR_M0_ae =  AR_M_ae*COEFP

   A11      = - AA_M_ae /TEFF
   A22      = - AR_M0_ae/TEFF**2
   A33      = - AR_M0_ae/TEFF**2*WEFF*DCM_ae - AE_M_ae /TEFF   *W0T

   CON      = 2.d0/DT_ae**2-3.d0*A11/(2.d0*DT_ae)-A22
   CON1     = 5.d0/DT_ae**2-2.d0*A11/DT_ae
   CON2     = A11/(2.d0*DT_ae)-4.d0/(DT_ae**2)
   CON3     = 1.d0/DT_ae**2

   CMOMB (1,istrip) =  CON1/CON*QTOP1 + CON2/CON*QTOP2 + CON3/CON*QTOP3 + A33/CON


!--- Final value of CL-CD-CM
   G1L      = CLIFTA(1,istrip)
   G2L      = CLIFTB(1,istrip)
   CLIFTDS  = 1.d0/WEFF**2*( WEFF *(G1L+G2L) + 0.5d0*AS_L_ae*ACHORD*W0T + 0.5d0*AK_L_ae*ACHORD*W1T )

   G2D      = CDRAGB(1,istrip)
   ASIG_D   = ASIG_D_ae*AEFF
   CDRAGDS  = 1.d0/WEFF**2*( WEFF *G2D + WEFF**2*CDLIN_ae + 0.5d0*ASIG_D*ACHORD*W0T )

   G2M      = CMOMB (1,istrip)
   CMOMDS   = 1.d0/WEFF**2*( WEFF*G2M +WEFF**2*CMLIN_ae + 0.5d0* (ASIGP_M_ae+AD_M_ae)*ACHORD*W0T &
                            +   ASIG_M_ae*WEFF*W1       + 0.5d0*             AS_M_ae *ACHORD*W1T   )

!--- FOILFS Additional
   if ( IFLAPLOC(iblade,lstrip)==1 ) then
     CLIFTDS  = clf (iblade, lstrip, NTIME) + 1/WEFF*G2L  - (CLLIN_pot-CLLIN_ae)
     CDRAGDS  = 1/WEFF*G2D + cdf (iblade, lstrip, NTIME)  + CDLIN_ae
     CMOMDS   = 1/WEFF*G2M + cmf (iblade, lstrip, NTIME)  + CMLIN_ae - cmpot(iblade,lstrip)
   endif


 END Subroutine ONERALDM
!--------------------------------------------------------------------------
! --  Subroutine :BEDDLEISH
! --  Beddoes-Leishman stall model for the lift drag and moment -----------
!--------------------------------------------------------------------------
 Subroutine BEDDLEISH ( RLD    , CHD    , WEFF  , ALPHA , THY  ,&
                        CLIFTDS, CDRAGDS, CMOMDS, istrip, NTIME   )
 Use Craft
 Use Caero

   implicit none

   real(8) :: RLD, CHD, WEFF, ALPHA, THY
   real(8) :: CLIFTDS,CDRAGDS,CMOMDS
   integer :: istrip,NTIME
!--- local
   real(8) :: ASTP(4)
   real(8) :: ACHORD,AEFF,AEFFM,ARSPAN,CLLIN_ae,DCL3D,DCD3D,DCLIFT,DCDRAG
   real(8) :: WEFFT,AEFFT, TEFF,COEF1
   real(8) :: TP,TF,DFALFL,BA1,BA2,B1,B2,AA1,AA2,AS34,X1,X2,X3,X4,COEF2,AE,AEM
   real(8) :: CLSTAE,CDSTAE,CMSTAE,XSST,XSST0,DT1,ACOEF,CLP,CLPP,AREF,AREFM
   real(8) :: CLSTAF,CDSTAF,CMSTAF,XSUN1,CNV,CLFS,DCDIND,DCDFPP,DCMFPP
   real(8) :: ASTFPP,ASTFAE


   ACHORD = CHD
   ARSPAN = dabs(RLD)

   AEFF   = ALPHA
   AEFFM  = AEFF * R2D_ae

   ANGEFF (istrip) = AEFF
   VELEFF (istrip) = WEFF

!--- Definition the DCL = CLlinear - CL steady
   CLLIN_ae = 0.5d0*DCLLIN_ae*dsin(2.0d0*(AEFF -AZER_ae))

!--- 3D effects
   DCL3D =  CLLIN_ae - CL_ae
   DCD3D =-(CDLIN_ae - CD_ae)

   if (AEFFM > 25.d0.and.AEFFM < 40.d0) then
      DCL3D = (1.d0-(AEFFM-25.d0)/15.d0)*DCL3D
      DCD3D = (1.d0-(AEFFM-25.d0)/15.d0)*DCD3D
   elseif (AEFFM >= 40.d0) then
      DCL3D = 0.d0
      DCD3D = 0.d0
   endif

      DCLIFT = 2.2d0*(ACHORD/ARSPAN)**1.5d0*(dcos(THY))**4*DCL3D
      DCDRAG = 2.2d0*(ACHORD/ARSPAN)**1.5d0*(dcos(THY))**4*DCD3D

!------ 3D effects not active
      DCLIFT = 0.d0
      DCDRAG = 0.d0

      CL_ae  = CL_ae + DCLIFT
      CD_ae  = CD_ae + DCDRAG

!--- Definition of ast
   call FOIL_AST (AGON,CLST,CMST,DCLLIN_ae,CMO_ae,AZER_ae,ASTP,NPP,NANGM)


!--- Beddoes-Leishman hysteresis parameters
   TP = 1.7d0
   TF = 3.0d0


!--- Parameters of the equation
   if     (NTIME <= NTIME_un  ) then

      CLIFTDS = CL_ae
      CDRAGDS = CD_ae
      CMOMDS  = CM_ae
      return

   elseif (NTIME == NTIME_un+1) then

      ANGEFFP (3,istrip) = ANGEFF (istrip)
      VELEFFP (3,istrip) = VELEFF (istrip)

      BXOLD1(3,istrip) = 0.d0
      BXOLD2(3,istrip) = 0.d0
      BXOLD3(3,istrip) = 0.d0
      BXOLD4(3,istrip) = 0.d0

      CLIFTDS = CL_ae
      CDRAGDS = CD_ae
      CMOMDS  = CM_ae

      return

   elseif (NTIME == NTIME_un+2) then

      ANGEFFP (2,istrip) = ANGEFF (istrip)
      VELEFFP (2,istrip) = VELEFF (istrip)

      BXOLD1(2,istrip) = 0.d0
      BXOLD2(2,istrip) = 0.d0
      BXOLD3(2,istrip) = 0.d0
      BXOLD4(2,istrip) = 0.d0

      CLIFTDS = CL_ae
      CDRAGDS = CD_ae
      CMOMDS  = CM_ae

      return

   elseif (NTIME == NTIME_un+3) then

      ANGEFFP (1,istrip) = ANGEFF (istrip)
      VELEFFP (1,istrip) = VELEFF (istrip)

      BXOLD1(1,istrip) = 0.d0
      BXOLD2(1,istrip) = 0.d0
      BXOLD3(1,istrip) = 0.d0
      BXOLD4(1,istrip) = 0.d0

      CLIFTDS = CL_ae
      CDRAGDS = CD_ae
      CMOMDS  = CM_ae

      return

   endif


!--- First order derivatives of WEFF, AEFF
   WEFFT  = ( 3.d0*VELEFF (  istrip) &
             -4.d0*VELEFFP(1,istrip) &
             +     VELEFFP(2,istrip) )/(2.d0*DT_ae)

   AEFFT  = ( 3.d0*ANGEFF (  istrip) &
             -4.d0*ANGEFFP(1,istrip) &
             +     ANGEFFP(2,istrip) )/(2.d0*DT_ae)

   DFALFL = AEFFT*ACHORD/WEFF

   TEFF   = 0.5d0*ACHORD/WEFF


!--- Lift Equation for the attached flow
   BA1=0.2940d0
   BA2=0.3310d0
   B1 =0.0664d0
   B2 =0.3266d0

   AA1=-2.d0*WEFF/ACHORD*(B1+ACHORD*WEFFT/(2.d0*WEFF**2))
   AA2=-2.d0*WEFF/ACHORD*(B2+ACHORD*WEFFT/(2.d0*WEFF**2))

   AS34=AEFF + 0.5d0*DFALFL


   COEF1 = 1.d0/(3.d0/(2.d0*DT_ae)-AA1)
   X1    = 2.d0/DT_ae*       COEF1*BXOLD1(1,istrip) &
          -1.d0/(2.d0*DT_ae)*COEF1*BXOLD1(2,istrip) &
          +                  COEF1*AS34*2.d0*WEFF/ACHORD*B1*BA1
   COEF2 = 1.d0/(3.d0/(2.d0*DT_ae)-AA2)
   X2    = 2.d0/DT_ae       *COEF2*BXOLD2(1,istrip) &
          -1.d0/(2.d0*DT_ae)*COEF2*BXOLD2(2,istrip) &
          +                  COEF2*AS34*2.d0*WEFF/ACHORD*B2*BA2

   BX1(1,istrip) = X1
   BX2(1,istrip) = X2

!--- The effective angle of attack
   AE  = AS34*(1.d0-BA1-BA2) + X1 + X2
   AEM = AE * R2D_ae

   call FOIL_CLCDCM (AEM,NPP,NANGM,AGON, &
                     CLST,CDST,CMST,CLSTAE,CDSTAE,CMSTAE)


   if (dabs((AE-AZER_ae)* R2D_ae) < 1.d0) then
      XSST0 =  1.d0
      goto 14
   endif

   if (dabs(AEM) > 30.d0) then
      XSST0 =  0.d0
      goto 14
   endif

   XSST0 = (2.d0*dsqrt(CLSTAE/(DCLLIN_ae*(AE-AZER_ae))) -1.d0)**2


   if (XSST0 > 1.d0) XSST0 = 1.d0
   if (XSST0 < 0.d0) XSST0 = 0.d0

14 continue

!--- Total Lift Force Coefficient
   CLP= DCLLIN_ae*(AE-AZER_ae) + PI_ae/2.d0*DFALFL


!--- Lift Equation for the separated flow
   DT1  = DT_ae * WEFF/ACHORD
   ACOEF=1.d0/(3.d0/(2.d0*DT1)+1.d0/TP)
   X3 = 2.d0/DT1       *ACOEF*BXOLD3(1,istrip) &
       -1.d0/(2.d0*DT1)*ACOEF*BXOLD3(2,istrip) &
       +            ACOEF*CLP/TP

   CLPP=X3

   BX3(1,istrip) = X3

!---
   AREF=CLPP/DCLLIN_ae + AZER_ae
   AREFM=AREF* R2D_ae

   call FOIL_CLCDCM (AREFM,NPP,NANGM,AGON, &
                     CLST,CDST,CMST,CLSTAF,CDSTAF,CMSTAF)

   if (dabs((AREF-AZER_ae)* R2D_ae).lt.1.d0) then
   XSST =  1.d0
   goto 15
   endif

   if (dabs(AREFM).gt.30.d0) then
   XSST =  0.d0
   goto 15
   endif

   XSST = (2.d0*dsqrt(CLSTAF/(DCLLIN_ae*(AREF-AZER_ae))) -1.d0)**2


   if (XSST.gt.1.d0) XSST = 1.d0
   if (XSST.lt.0.d0) XSST = 0.d0

15 continue

   ACOEF= 1.d0/(3.d0/(2.d0*DT1)+1.d0/TF)
   X4   = 2.d0/DT1       *ACOEF*BXOLD4(1,istrip) &
         -1.d0/(2.d0*DT1)*ACOEF*BXOLD4(2,istrip) &
         +                ACOEF*XSST/TF

   XSUN1=X4

   BX4(1,istrip) = X4

   if (XSUN1 > 1.d0) XSUN1 = 1.d0
   if (XSUN1 < 0.d0) XSUN1 = 0.d0

!--- dynamic stall term not included
   CNV = 0.d0


!--- dynamic lift
   if (XSST0 > 0.999999d0) then
      CLFS = CLSTAE/2.d0
   else
      CLFS =(CLSTAE - DCLLIN_ae*(AE-AZER_ae)*XSST0)/(1.d0-XSST0)
   endif
      CLIFTDS=DCLLIN_ae*(AE-AZER_ae)*XSUN1+PI_ae/2.d0*DFALFL+CLFS*(1.d0-XSUN1)


!--- Drag Equation
   DCDIND = (AEFF-AE)*CLIFTDS
   DCDFPP = (CDSTAE-CDLIN_ae)*( (0.5d0*(1.d0-dsqrt(XSUN1)))**2 &
                               -(0.5d0*(1.d0-dsqrt(XSST0)))**2 )


   CDRAGDS = CDSTAE + DCDIND + DCDFPP


!--- Moment Equation
   ASTFPP = ASTP(1)          &
           +ASTP(2)*XSUN1    &
           +ASTP(3)*XSUN1**2 &
           +ASTP(4)*XSUN1**3
   ASTFAE = ASTP(1)          &
           +ASTP(2)*XSST0    &
           +ASTP(3)*XSST0**2 &
           +ASTP(4)*XSST0**3

   DCMFPP =-CLIFTDS*(ASTFPP - ASTFAE)

   CMOMDS = CMSTAE + DCMFPP - PI_ae/4.d0*DFALFL


 END Subroutine BEDDLEISH
!---------------------------------------------------------------------
!     Subroutine defining ast parameter
!
!---------------------------------------------------------------------
 Subroutine FOIL_AST (AGON,CLST,CMST,CLSLO,CMO,AZER,BSYS,NPP,NMAX)
!---------------------------------------------------------------------

   implicit none

   integer :: NPP, NMAX
   real(8) :: AGON(NMAX),CLST(NMAX),CMST(NMAX)
   real(8) :: CLSLO,CMO,AZER
   real(8) :: BSYS(4)
!--- local
   real(8) :: XSEP(NMAX),AST (NMAX),ASYS(4,4)
   integer :: INDX(4)
   real(8) :: pi,R2D,ALFA,ALFAM
   integer :: i,k,NPP1, ierr


   pi  = dacos(-1.d0)
   R2D = 180d0/pi

   k = 0
   do i=1,NPP
      ALFA  = AGON(i)/R2D
      ALFAM = AGON(i)
      if (ALFAM.lt.0.d0.or.ALFAM.gt.30.d0) cycle
      k = k + 1
      if (ALFA.lt.AZER+0.01d0.and.ALFA.gt.AZER-0.01d0) then
         XSEP(k)  = 1.d0
      else
         XSEP(k)  = (2.d0*sqrt(CLST(i)/(CLSLO*(ALFA-AZER))) -1.d0)**2
         if (XSEP(k).gt.1.d0) XSEP(k) = 1.d0
         if (XSEP(k).lt.0.d0) XSEP(k) = 0.d0
      endif
      AST (k)  = (CMST(i)-CMO)/CLST(i)
   enddo

   NPP1 = k

   ASYS = 0.d0;
   BSYS = 0.d0;

!  do    k = 1, NNP
!     do j = 1, 4
!     do i = 1, 4
!        ASYS(i,j) = ASYS(i,j) + XSEP(k)**(i-1)*XSEP(k)**(j-1)
!     enddo
!        BSYS(  j) = BSYS(  j) + AST(k)        *XSEP(k)**(j-1)
!     enddo
!  enddo

   do i         = 1, NPP1
      ASYS(1,1) = ASYS(1,1) + 1.d0
      ASYS(1,2) = ASYS(1,2) + XSEP(i)
      ASYS(1,3) = ASYS(1,3) + XSEP(i)**2
      ASYS(1,4) = ASYS(1,4) + XSEP(i)**3

      ASYS(2,1) = ASYS(2,1) + 1.d0       *XSEP(i)
      ASYS(2,2) = ASYS(2,2) + XSEP(i)    *XSEP(i)
      ASYS(2,3) = ASYS(2,3) + XSEP(i)**2 *XSEP(i)
      ASYS(2,4) = ASYS(2,4) + XSEP(i)**3 *XSEP(i)

      ASYS(3,1) = ASYS(3,1) + 1.d0       *XSEP(i)**2
      ASYS(3,2) = ASYS(3,2) + XSEP(i)    *XSEP(i)**2
      ASYS(3,3) = ASYS(3,3) + XSEP(i)**2 *XSEP(i)**2
      ASYS(3,4) = ASYS(3,4) + XSEP(i)**3 *XSEP(i)**2

      ASYS(4,1) = ASYS(4,1) + 1.d0       *XSEP(i)**3
      ASYS(4,2) = ASYS(4,2) + XSEP(i)    *XSEP(i)**3
      ASYS(4,3) = ASYS(4,3) + XSEP(i)**2 *XSEP(i)**3
      ASYS(4,4) = ASYS(4,4) + XSEP(i)**3 *XSEP(i)**3

      BSYS(1)   = BSYS(1) + AST(i)
      BSYS(2)   = BSYS(2) + AST(i)*XSEP(i)
      BSYS(3)   = BSYS(3) + AST(i)*XSEP(i)**2
      BSYS(4)   = BSYS(4) + AST(i)*XSEP(i)**3
   enddo


   CALL DGESV (4, 1, ASYS, 4, INDX, BSYS, 4, ierr)
!  call LUDCMP (ASYS , 4, 4, INDX, d   )
!  call LUBKSB (ASYS , 4, 4, INDX, BSYS)


 END Subroutine FOIL_AST
!---------------------------------------------------------------------
!
!-- Subroutine : Aero2Elast_bem  -------------------------------------
!
!----------------------------------------------------------------------
 Subroutine Aero2Elast_bem
!----------------------------------------------------------------------

 Use Cbeam
 Use Craft

   implicit none

   integer :: lstrip, nboda


   do nboda  = 1, NBLADE_ae
   do lstrip = 1, NSTRIP
!------ Aerodynamic Loads before pitch angle
      FSTRIP_el (1:6,nboda,lstrip) = strip(nboda, lstrip)%FCP (1:6)
      XCAER_el  (1  ,nboda,lstrip) = strip(nboda, lstrip)%xae
      XCAER_el  (3  ,nboda,lstrip) = strip(nboda, lstrip)%zae
   enddo
   enddo


 END Subroutine Aero2Elast_bem
!
!
!
!-----------------------------------------------------------------------
!
! Subroutine : UPDATE_AERO
!
! Updates previous solutions for a and a'
! Updates ONERA or BL variables
!
!----------------------------------------------------------------------
 Subroutine UPDATE_AERO (NTIME)
!----------------------------------------------------------------------

 Use Craft

   implicit none

   integer, intent(in) :: NTIME
   integer :: lstrip, iblad


!--- Save solution for induction factors
   do iblad  = 1, NBLADE_ae
   do lstrip = 1, NSTRIP
      strip(iblad,lstrip)%AINDTB_PRE  = strip(iblad,lstrip)%AINDTB
      strip(iblad,lstrip)%AINDPTB_PRE = strip(iblad,lstrip)%AINDPTB
   enddo
   enddo

   call UPDATE_ONERA (NTIME)
   call UPDATE_BL    (NTIME)
   call UPDATE_FOILFS


 END Subroutine UPDATE_AERO
!----------------------------------------------------------------------
 Subroutine UPDATE_ONERA (NTIME)
!----------------------------------------------------------------------

 Use Craft
 Use Caero

   implicit none

   integer, intent(in) :: NTIME
   integer :: lstrip, istrip, iblad


   if (IDYNSTALL /= 1     ) return
   if (NTIME <= NTIME_un+3) return

      istrip = 0
   do iblad  = 1, NBLADE_ae
   do lstrip = 1, NSTRIP
      istrip = istrip + 1

      VELEFFP (3,istrip) = VELEFFP (2,istrip)
      VELEFFP (2,istrip) = VELEFFP (1,istrip)
      VELEFFP (1,istrip) = VELEFF  (  istrip)
      ANGEFFP (3,istrip) = ANGEFFP (2,istrip)
      ANGEFFP (2,istrip) = ANGEFFP (1,istrip)
      ANGEFFP (1,istrip) = ANGEFF  (  istrip)
      CLOLDA  (3,istrip) = CLOLDA  (2,istrip)
      CLOLDA  (2,istrip) = CLOLDA  (1,istrip)
      CLOLDA  (1,istrip) = CLIFTA  (1,istrip)
      CLOLDB  (3,istrip) = CLOLDB  (2,istrip)
      CLOLDB  (2,istrip) = CLOLDB  (1,istrip)
      CLOLDB  (1,istrip) = CLIFTB  (1,istrip)
      CDOLDB  (3,istrip) = CDOLDB  (2,istrip)
      CDOLDB  (2,istrip) = CDOLDB  (1,istrip)
      CDOLDB  (1,istrip) = CDRAGB  (1,istrip)
      CMOLDB  (3,istrip) = CMOLDB  (2,istrip)
      CMOLDB  (2,istrip) = CMOLDB  (1,istrip)
      CMOLDB  (1,istrip) = CMOMB   (1,istrip)
   enddo
   enddo


 END Subroutine UPDATE_ONERA
!-----------------------------------------------------------------------
 SUBROUTINE UPDATE_BL (NTIME)
!-----------------------------------------------------------------------

 Use Craft
 Use Caero

   implicit none

   integer, intent(in) :: NTIME
   integer :: lstrip, istrip, iblad


   if (IDYNSTALL /= 2     ) return
   if (NTIME <= NTIME_un+3) return

      istrip = 0
   do iblad  = 1, NBLADE_ae
   do lstrip = 1, NSTRIP
      istrip = istrip + 1

      VELEFFP (3,istrip) = VELEFFP (2,istrip)
      VELEFFP (2,istrip) = VELEFFP (1,istrip)
      VELEFFP (1,istrip) = VELEFF  (  istrip)
      ANGEFFP (3,istrip) = ANGEFFP (2,istrip)
      ANGEFFP (2,istrip) = ANGEFFP (1,istrip)
      ANGEFFP (1,istrip) = ANGEFF  (  istrip)
      BXOLD1  (3,istrip) = BXOLD1  (2,istrip)
      BXOLD1  (2,istrip) = BXOLD1  (1,istrip)
      BXOLD1  (1,istrip) = BX1     (1,istrip)
      BXOLD2  (3,istrip) = BXOLD2  (2,istrip)
      BXOLD2  (2,istrip) = BXOLD2  (1,istrip)
      BXOLD2  (1,istrip) = BX2     (1,istrip)
      BXOLD3  (3,istrip) = BXOLD3  (2,istrip)
      BXOLD3  (2,istrip) = BXOLD3  (1,istrip)
      BXOLD3  (1,istrip) = BX3     (1,istrip)
      BXOLD4  (3,istrip) = BXOLD4  (2,istrip)
      BXOLD4  (2,istrip) = BXOLD4  (1,istrip)
      BXOLD4  (1,istrip) = BX4     (1,istrip)
   enddo
   enddo
 

 END SUBROUTINE UPDATE_BL
!
!
!
!--- Newton routines
!----------------------------------------------------------------------
!
!  Given an initial guess x(1:n) for a root in n dimensions, find the root by a globally
!  convergent Newton’s method. The vector of functions to be zeroed, called fvec(1:n)
!  in the routine below, is returned by a user-supplied subroutine that must be called
!  and have the declaration subroutine BEMT (n,x,fvec).
!  The output quantity check is false on a normal return and true if the routine has
!  converged to a local minimum of the function fmin defined below. In this case try
!  restarting from a different initial guess.
!
!  Parameters: 
!      maxits  is the maximum number of iterations;
!      tolf    sets the convergence criterion on function values;
!      tolmin  sets the criterion for deciding whether spurious convergence
!              to a minimum of fmin has occurred;
!      tolx    is  the convergence criterion on δx;
!      stpmx   is the scaled maximum step length allowed in line searches.

!----------------------------------------------------------------------
 Subroutine newt (x,n,check)
!----------------------------------------------------------------------

   implicit none

   integer, intent(in   ) :: n
   real(8), intent(inout) :: x(n)
   logical, intent(  out) :: check
!--- USES fdjac,fmin,lnsrch,lubksb,ludcmp
   integer, parameter :: maxits=200
   real(8), parameter :: tolf=1.d-4,tolmin=1.d-6,tolx=1.d-7,stpmx=100.d0
   real(8) :: fvec(n)
   integer :: i,its,indx(n), ierr
   real(8) :: d,den,f,fold,stpmax,sum,temp,test,fjac(n,n),g(n)
   real(8) :: p(n),xold(n),fmin
   external fmin

   integer ::     NTIME,lstrip,ibladen
   common /status/NTIME,lstrip,ibladen
   save   /status/


      f    = fmin(n,x,fvec)
      test = 0.d0
   do i    = 1, n
     if (dabs(fvec(i)).gt.test) test=dabs(fvec(i))
   enddo !i

   if (test.lt..01d0*tolf) then; check=.false.; return; endif

      sum       = dot_product(x(1:n),x(1:n))
      stpmax    = stpmx*max(dsqrt(sum),dble(n))

   do its       = 1, maxits
      call fdjac(n,x,fvec,fjac)

      g   (1:n) = matmul(transpose(fjac(1:n,1:n)),fvec(1:n))
      xold(1:n) = x(1:n)
      fold      = f
      p   (1:n) =-fvec(1:n)
      
      call ludcmp(fjac,n,n,indx,d,ierr);
      if (ierr==1) then; write(*,*)'newton1a',lstrip; check=.true.; return; endif
      call lubksb(fjac,n,n,indx,p)
      call lnsrch(n,xold,fold,g,p,x,f,fvec,stpmax,check,fmin,ierr)
      if (ierr==1) then; write(*,*)'newton1b',lstrip; check=.true.; return; endif
      
         test = 0.d0
      do i    = 1, n
        if (dabs(fvec(i)).gt.test) test=dabs(fvec(i))
      enddo !i

!w    write(*,*) its,test
      if (test.lt.tolf) then; check=.false.; return; endif

      if (check) then
            test = 0.d0
            den  = max(f,.5d0*n)
         do i    = 1, n
            temp = dabs(g(i))*max(dabs(x(i)),1.d0)/den
           if (temp.gt.test) test=temp
         enddo !i
         if (test.lt.tolmin) then
           check=.true.; write(*,*)'newton2',lstrip
         else
           check=.false.
         endif
         return
      endif

         test = 0.d0
      do i    = 1, n
         temp = (dabs(x(i)-xold(i)))/max(dabs(x(i)),1.d0)
         if (temp.gt.test) test=temp
      enddo !i
      if (test.lt.tolx) return
   enddo !its

   write(* ,*) 'maxits exceeded in newt',lstrip
   write(10,*) 'maxits exceeded in newt',lstrip


 END Subroutine newt
!----------------------------------------------------------------------
!
!  Computes forward-difference approximation to Jacobian. On input, x(1:n) is the point
!  at which the Jacobian is to be evaluated, fvec(1:n) is the vector of function values at
!  the point, and n is the physical dimension of the Jacobian array df(1:n,1:n) which is
!  output. subroutine BEMT (n,x,f) is a fixed-name, user-supplied routine that returns
!  the vector of functions at x.
!
!  Parameters: 
!  EPS is the approximate square root of the machine precision.
!
!----------------------------------------------------------------------
 Subroutine fdjac (n,x,fvec,df)
!----------------------------------------------------------------------

   implicit none

!--- USES BEMT 
   integer, intent(in   ) :: n
   real(8), intent(in   ) :: fvec(n)
   real(8), intent(inout) :: x(n)
   real(8), intent(  out) :: df(n,n)

   real(8), parameter :: eps=1.d-4
   real(8)            :: h,temp,f(n)
   integer            :: i,j


   do j    = 1, n
      temp = x(j)
      h    = eps*dabs(temp)
      if (h.eq.0.d0) h=eps
      x(j) = temp+h
      h    = x(j)-temp
      call BEMT (n,x,f)
      x(j) = temp
      do i = 1, n
         df(i,j)=(f(i)-fvec(i))/h
     enddo !i
   enddo !j


 END Subroutine fdjac
!----------------------------------------------------------------------
!
!  Returns f = 1 F · F at x. subroutine BEMT (n,x,f) is a fixed-name, user-supplied
!  routine that returns the vector of functions at x.
!
!----------------------------------------------------------------------
 Function fmin (n,x,fvec)
!----------------------------------------------------------------------
   implicit none

!--- USES BEMT 
   integer, intent(in)  :: n
   real(8), intent(in)  :: x(n)
   real(8), intent(out) :: fvec(n)
   real(8) :: fmin


   call BEMT (n,x,fvec)

   fmin = 0.5d0*dot_product(fvec(1:n),fvec(1:n))


 END Function fmin
!----------------------------------------------------------------------
!
!  Given an n-dimensional point xold(1:n), the value of the function and gradient there,
!  fold and g(1:n), and a direction p(1:n), finds a new point x(1:n) along the direction
!  p from xold where the function func has decreased “sufficiently.” The new function value
!  is returned in f. stpmax is an input quantity that limits the length of the steps so that you
!  do not try to evaluate the function in regions where it is undefined or subject to overflow.
!  p is usually the Newton direction. The output quantity check is false on a normal exit.
!  It is true when x is too close to xold. In a minimization algorithm, this usually signals
!  convergence and can be ignored. However, in a zero-finding algorithm the calling program
!  should check whether the convergence is spurious.
!  Parameters: ALF ensures sufficient decrease in function value; TOLX is the convergence
!  criterion on ∆x.
!
!----------------------------------------------------------------------
   Subroutine lnsrch (n,xold,fold,g,p,x,f,fvec,stpmax,check,func,ierr) !fmin
!----------------------------------------------------------------------

   implicit none

   integer :: n, ierr
   logical :: check
   real(8), parameter :: alf=1.d-4, tolx=1.d-7
   real(8) :: f,fold,stpmax,g(n),p(n),x(n),xold(n),func,fvec(n)
   external func !fmin
!--- USES func   !fmin
   integer :: i
   real(8) :: a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope
   real(8) :: sum,temp,test,tmplam


   ierr  = 0
   check = .false.
   sum   = dsqrt(dot_product(p(1:n),p(1:n)))
   if(sum.gt.stpmax)then
      p(1:n)=p(1:n)*stpmax/sum
   endif
   slope=dot_product(g(1:n),p(1:n))
   test=0.d0
   do i=1,n
     temp=dabs(p(i))/max(dabs(xold(i)),1.d0)
     if(temp.gt.test)test=temp
   enddo
   alamin=tolx/test
   alam=1.d0

1     x(1:n)=xold(1:n)+alam*p(1:n)
      f=func(n,x,fvec)
      if(alam.lt.alamin)then
        x(1:n)=xold(1:n)
        check=.true.
        return
      elseif (f.le.fold+alf*alam*slope) then
        return
      else
        if (alam.eq.1.d0) then
          tmplam=-slope/(2.d0*(f-fold-slope))
        else
          rhs1=f-fold-alam*slope
          rhs2=f2-fold2-alam2*slope
          a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
          b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
          if(a.eq.0.d0)then
            tmplam=-slope/(2.d0*b)
          else
            disc=b*b-3.d0*a*slope
            if(disc.lt.0.d0) then; write(*,*) 'roundoff problem in lnsrch'; ierr=1; return; endif
            tmplam=(-b+dsqrt(disc))/(3.d0*a)
          endif
          if(tmplam.gt.0.5d0*alam)tmplam=0.5d0*alam
        endif
      endif
      alam2 = alam
      f2    = f
      fold2 = fold
      alam  = max(tmplam,0.1d0*alam)
   goto 1


 END Subroutine lnsrch
!----------------------------------------------------------------------
!--- LU decomposition routines
!----------------------------------------------------------------------
 Subroutine lubksb (a,n,np,indx,b)
!----------------------------------------------------------------------

   implicit none

   integer, intent(in   ) :: n,np,indx(n)
   real(8), intent(in   ) :: a(np,np)
   real(8), intent(inout) :: b(n)

   integer :: i,ii,j,ll
   real(8) :: sum


   ii=0
   do i=1,n
     ll=indx(i)
     sum=b(ll)
     b(ll)=b(i)
     if (ii.ne.0)then
       do j=ii,i-1
         sum=sum-a(i,j)*b(j)
       enddo !j
     else if (sum.ne.0.d0) then
       ii=i
     endif
     b(i)=sum
   enddo !i
   do i=n,1,-1
     sum=b(i)
     do j=i+1,n
       sum=sum-a(i,j)*b(j)
     enddo !j
     b(i)=sum/a(i,i)
   enddo !i


 END Subroutine lubksb
!----------------------------------------------------------------------
 Subroutine ludcmp (a,n,np,indx,d,ierr)
!----------------------------------------------------------------------

   implicit none

   integer, intent(in   ) :: n,np
   integer, intent(  out) :: indx(n), ierr
   real(8), intent(inout) :: a(np,np)
   real(8), intent(  out) :: d

   integer, parameter :: nmax=500
   real(8), parameter :: tiny=1.0d-20
   integer :: i,imax,j,k
   real(8) :: aamax,dum,sum,vv(nmax)


   ierr=1
   d=1.d0
   do i=1,n
     aamax=0.d0
     do j=1,n
       if (dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
     enddo
     if (aamax.eq.0.d0) then; write(*,*)'singular matrix in ludcmp'; return; endif
     vv(i)=1.d0/aamax
   enddo
   do j=1,n
     do i=1,j-1
       sum=a(i,j)
       do k=1,i-1
         sum=sum-a(i,k)*a(k,j)
       enddo !k
       a(i,j)=sum
     enddo !i
     aamax=0.d0
     do i=j,n
       sum=a(i,j)
       do k=1,j-1
         sum=sum-a(i,k)*a(k,j)
       enddo !k
       a(i,j)=sum
       dum=vv(i)*dabs(sum)
       if (dum.ge.aamax) then
         imax=i
         aamax=dum
       endif
     enddo !i
     if (j.ne.imax)then
       do k=1,n
         dum=a(imax,k)
         a(imax,k)=a(j,k)
         a(j,k)=dum
       enddo !k
       d=-d
       vv(imax)=vv(j)
     endif
     indx(j)=imax
     if(a(j,j).eq.0.d0)a(j,j)=TINY
     if(j.ne.n)then
       dum=1.d0/a(j,j)
       do i=j+1,n
         a(i,j)=a(i,j)*dum
       enddo !i
     endif
   enddo !j

   ierr=0


 END Subroutine ludcmp

include "foilfs.f90"
!include "minpack.f90"
