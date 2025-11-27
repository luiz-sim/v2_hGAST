!----------------------------------------------------------------------
 Module mod_types_el
!----------------------------------------------------------------------
 type :: type_boundco    ! NBODT_el  !total sub-bodies
!----------------------------------------------------------------------
!---- Boundary Conditions could be either zero or 
!---- depended [upon other sub-body/element/node]
    integer              :: nbcpb           !<                          ! number of b.c per sub-body
    integer, allocatable :: nel     (:    ) !<(nbcpb    )               ! nel  at which b.c. are applied
    integer, allocatable :: nod     (:    ) !<(nbcpb    )               ! nod  at which b.c. are applied [1 or NNPE_el]
    integer, allocatable :: nbcon   (:    ) !<(nbcpb    )               ! body which sets the b.c. (in turn get loads), when zero: fixed conditions (nelcon,nodcon dummy)
    integer, allocatable :: nelcon  (:    ) !<(nbcpb    )               ! nel  which sets the b.c. (in turn get loads)
    integer, allocatable :: nodcon  (:    ) !<(nbcpb    )               ! nod  which sets the b.c. (in turn get loads) [1 or NNPE_el]
    integer, allocatable :: indx    (:,:  ) !<(nbcpb,6  )               ! 0:free, 1:b.c. for each dofs
    real(8), allocatable :: Atrans1 (:,:,:) !<(nbcpb,6,6)               ! transformation matrix for B.C
    real(8), allocatable :: Atrans2 (:,:,:) !<(nbcpb,6,6)               ! transformation matrix for B2B loads

 END type type_boundco
!----------------------------------------------------------------------
 type :: type_foundation ! NBODTFND_el  !total sub-bodies being part of the foundation
!----------------------------------------------------------------------
!---- Foundation (Springs)
    integer              :: nb                                          ! sub-body           at which foundation acts
    integer              :: nteb_fnd                                    ! number of elements at which foundation acts
    integer, allocatable :: nel     (:    ) !(nteb_fnd            )     ! nel of sub-body    at which foundation acts
    real(8), allocatable :: Stiff   (:,:,:) !(nteb_fnd,NPFPE,NDFPE)     ! Local element Stiffness matrix from foundation
    real(8), allocatable :: Damp    (:,:,:) !(nteb_fnd,NDFPE,NDFPE)     ! Local element Damping   matrix from foundation)

 END type type_foundation
!----------------------------------------------------------------------
 type :: type_beam_timo  ! INODNELPROP_el, used in <beam_timo>
!----------------------------------------------------------------------
!---- Mass variables
    real(8)              :: DENS_el                                     ! Density/Length
    real(8)              :: XCMD_el                                     ! Mass (Gravity) Centre X
    real(8)              :: ZCMD_el                                     ! Mass (Gravity) Centre Z
    real(8)              :: AMOMX_el                                    ! Mass Moment around z direction
    real(8)              :: AMOMZ_el                                    ! Mass Moment around x direction
    real(8)              :: RIXX_el                                     ! Mass inertia along z direction
    real(8)              :: RIZZ_el                                     ! Mass inertia along x direction
    real(8)              :: RIXZ_el                                     ! Cross mass inertia
    real(8)              :: POLI_el                                     ! Polar moment of Inertia
!---- Timoshenko Full Stiffness Matrix
    real(8)              :: K      (21)                                 ! Upper symmetric stiffness matrix in the order of hGAST
!---- Bend-Sweep offsets
    real(8)              :: XL_el(3)                                    ! body local coordinates

 END type type_beam_timo
!----------------------------------------------------------------------
 type :: type_substr_mor  ! INODNELPROP_el, used in <beam_prop1>
!----------------------------------------------------------------------

    real(8)              :: Cm_el                                       ! Morison Inertia Coefficient, Cm = 1+Ca (Ca added-mass coeff.)
    real(8)              :: Cd_el                                       ! Morison Drag    Coefficient
    real(8)              :: DIAMET_el                                   ! Outer Diameter
    real(8)              :: INDIAMET_el                                 ! Inner Diameter

 END type type_substr_mor
!----------------------------------------------------------------------
 type :: type_conc_mass   ! NCOMA_el, used in <conc_mass>
!----------------------------------------------------------------------

    integer              :: NBCOMA_el                                   ! Body     at which mass is located
    integer              :: NSBCOMA_el                                  ! Sub-body at which mass is located
    integer              :: NELCOMA_el                                  ! Element  at which mass is located
    real(8)              :: Hta                                         ! local position along elastic axis at which mass is attached
    real(8)              :: Xoff                                        ! X   local offset wrt elastic axis
    real(8)              :: Yoff                                        ! Y   local offset wrt HCOMA_el
    real(8)              :: Zoff                                        ! Z   local offset wrt elastic axis
    real(8)              :: Mass                                        ! concentrated  mass
    real(8)              :: Ixx                                         ! Ixx wrt local mass c.s.
    real(8)              :: Iyy                                         ! Iyy wrt local mass c.s.
    real(8)              :: Izz                                         ! Izz wrt local mass c.s.
    real(8)              :: Ixy                                         ! Ixy
    real(8)              :: Ixz                                         ! Ixz
    real(8)              :: Iyz                                         ! Iyz
    real(8)              :: IPx                                         ! polar inertia around x axis wrt local mass c.s.
    real(8)              :: IPy                                         ! polar inertia around y axis wrt local mass c.s.
    real(8)              :: IPz                                         ! polar inertia around z axis wrt local mass c.s.

 END type type_conc_mass
!----------------------------------------------------------------------
 type :: type_transf_mat   ! NNODT_el, used in <transf_mat>
!----------------------------------------------------------------------

!---- Transformation Matrices for every sub-body node
    real(8)              :: A_el       (  3,3)                          ! global base [A_el . Xloc  = Xglob]
    real(8)              :: AT_el      (  3,3)                          ! local  base [AT_el. Xglob = Xloc ]
    real(8)              :: DA_el      (  3,3)                          ! 1st time derivative of A_el
    real(8)              :: R_el       (    3)                          ! Position Vector
    real(8)              :: DR_el      (    3)                          ! 1st time derivative of R_el
    real(8)              :: ATDA_el    (  3,3)                          !
    real(8)              :: ATDDA_el   (  3,3)                          !
    real(8)              :: ATDR_el    (    3)                          !
    real(8)              :: ATDDR_el   (    3)                          !
!---- Qs
    real(8), allocatable :: A0_el      (:,:,:)     !<( NQS, 3, 3 )      !
    real(8), allocatable :: AT0_el     (:,:,:)     !<( NQS, 3, 3 )      !
    real(8), allocatable :: ATDA0_el   (:,:,:)     !<( NQS, 3, 3 )      !
    real(8), allocatable :: ATDA1_el   (:,:,:)     !<( NQS, 3, 3 )      !
    real(8), allocatable :: ATDDA0_el  (:,:,:)     !<( NQS, 3, 3 )      !
    real(8), allocatable :: ATDDA1_el  (:,:,:)     !<( NQS, 3, 3 )      !
    real(8), allocatable :: ATDDA2_el  (:,:,:)     !<( NQS, 3, 3 )      !
    real(8), allocatable :: R0_el      (  :,:)     !<( NQS,    3 )      !
    real(8), allocatable :: ATDDR0_el  (  :,:)     !<( NQS,    3 )      !
    real(8), allocatable :: ATDDR1_el  (  :,:)     !<( NQS,    3 )      !
    real(8), allocatable :: ATDDR2_el  (  :,:)     !<( NQS,    3 )      !
    real(8), allocatable :: ATDR0_el   (  :,:)     !<( NQS,    3 )      !
    real(8), allocatable :: ATDR1_el   (  :,:)     !<( NQS,    3 )      !
!---- only for blades
    real(8), allocatable :: ATPA_el    (  :,:)     !
!   real(8), allocatable :: DATPA_el   (  :,:)     !
!   real(8), allocatable :: DDATPA_el  (  :,:)     !
!---- only for jacket
    real(8), allocatable :: Aja_el     (  :,:)     !
    real(8), allocatable :: Rja_el     (    :)     !

 END type type_transf_mat
!----------------------------------------------------------------------
 type :: type_QSnode       ! NNODT_el, used in 
!----------------------------------------------------------------------

    integer              :: Tot_el    (2)       !
    integer, allocatable :: IORD_el   (:)       ! NQS1
   
 END type type_QSnode
!----------------------------------------------------------------------
 type :: type_ACCU        ! 1, used in 
!----------------------------------------------------------------------

    integer, allocatable :: NDFTBACC_el  (:)       !<(0:NBODBT_el)                    ! Accumulative number of dofs per     body
    integer, allocatable :: NQSACC_el    (:)       !<(0:NBODBT_el)
    integer, allocatable :: NQSBACC_el   (:,:)     !<(  NBODBT_el,0:NSUBperBM_el+1)
!Subbody
    integer, allocatable :: NDFTACC_el   (:)       !<(0:NBODT_el)                     ! Accumulative number of dofs per sub-body

 END type type_ACCU
!----------------------------------------------------------------------
 type :: type_body        ! NBOBDT_el, used in 
!----------------------------------------------------------------------
!Body
    integer              :: NBODSUB_el             !
    integer              :: NTEBB_el               !
    integer              :: IDOF_el                !
!Body, Subbody                                                                                                                                      
    integer, allocatable :: NBODTGLB_el  (:)       !<(  NBODSUB_el      )
    integer, allocatable :: NNODTGLB_el  (:)       !<(  NBODSUB_el+1    )
!-- Precurve angle
    real(8), allocatable :: PRECURV      (:)       !<(0:NBODSUB_el+1    )
    real(8), allocatable :: PRECURV_mat  (:,:,:)   !<(0:NBODSUB_el+1,3,3)

 END type type_body
!----------------------------------------------------------------------
 type :: type_subbody      ! NBODT_el, used in 
!----------------------------------------------------------------------
!-- Subbody
    integer              :: IBTYP_el              !
    integer              :: IAENUM                !
    integer              :: NTEB_el               !                     ! Number of elements  per sub-body
    integer              :: NNTB_el               !                     ! Number of nodes     per sub-body
    integer              :: NDFTB_el              !                     ! Number of dofs      per sub-body
!!! integer              :: NNPE_el               !                     ! Number of nodes     per element
    integer              :: NDFPE_el              !                     ! Number of dofs      per element
    integer              :: NEQPE_el              !                     ! Number of equations per element
    integer, allocatable :: NDFPN_el      (:  )   !<(  NNPE_el)         ! Number of dofs      per node
    integer, allocatable :: NDFPNACC_el   (:  )   !<(0:NNPE_el)         ! Accumulative number of dofs per element node
!-- Subbody,element,node
    integer, allocatable :: NODMB         (:,:)   !<(  NTEB_el, NNPE_el)! Local element to glabal body nodes correspondance
    integer, allocatable :: NODCB         (:  )   !<(  NNTB_el         )
    integer, allocatable :: NDFPNBACC_el  (:  )   !<(0:NNTB_el         )! Accumulative number of dofs per sub-body node
    integer, allocatable :: INODNEL_el    (:,:)   !<(  NTEB_el, 2      )
!-- Local Element Matrices
    real(8), allocatable :: AMLOC_el      (:,:,:) !<(NTEB_el, NDFPE_el, NDFPE_el)
    real(8), allocatable :: ACLOC_el      (:,:,:) !<(NTEB_el, NDFPE_el, NDFPE_el)
    real(8), allocatable :: AKLOC_el      (:,:,:) !<(NTEB_el, NDFPE_el, NDFPE_el)
    real(8), allocatable :: AFLOC_el      (:,:  ) !<(NTEB_el, NDFPE_el          )
    real(8), allocatable :: AMLOCNQ_el    (:,:,:) !<(NTEB_el, NDFPE_el, NQS     )
    real(8), allocatable :: ACLOCNQ_el    (:,:,:) !<(NTEB_el, NDFPE_el, NQS     )
    real(8), allocatable :: AKLOCNQ_el    (:,:,:) !<(NTEB_el, NDFPE_el, NQS     )
!-- Other Variables
    real(8)              :: ALENGB_el             !<                    ! Length of sub-body
    real(8), allocatable :: HTA_el       (:  )    !<(NTEB_el+1   )      ! Position along elastic axis of the nodes of each element
    real(8), allocatable :: ALENG_el     (:  )    !<(NTEB_el     )      ! Length of element
    real(8), allocatable :: PHIX_el      (:  )    !<(NTEB_el     )      ! 
    real(8), allocatable :: PHIZ_el      (:  )    !<(NTEB_el     )      ! 
    real(8), allocatable :: FCPL_el      (:,:)    !<(NTEB_el  , 3)      ! External Force  at the local system of the body
    real(8), allocatable :: AMPL_el      (:,:)    !<(NTEB_el  , 3)      ! External Moment at the local system of the body
    real(8), allocatable :: Rinit        (:,:)    !<(NTEB_el+1, 6)
!#ELEMENT HTA_el is replaced by XL(2)
!!  real(8), allocatable :: XL            (:,:  ) !(nel+1,3  )
!!  real(8), allocatable :: AE            (:,:,:) !(nel  ,3,3)

!-- Damping Coefficients
!   real(8), allocatable :: COEFM_el                                    ! Rayleigh   Damping mass      coefficient [NBODBT_el]
!   real(8), allocatable :: COEFK_el                                    ! Rayleigh   Damping stiffness coefficient [NBODBT_el]
!   real(8), allocatable :: CCRIT_el                                    ! Structural Damping coefficients          [0:NBODBTWT_el+1]
!   real(8),             :: FREQ0_el                                    !

 END type type_subbody

!----------------------------------------------------------------------
 type :: type_buoyancy_tap ! NNBuoyaTap_el, used in 
!----------------------------------------------------------------------

    integer              :: nbod                                        ! at which body [given only counting the jacket bodies and modified inside the code]
    integer              :: nel                                         ! at which nel
    integer              :: nod                                         ! 1:start, 2:end wrt local body c.s. [start or end of the element]
    integer              :: itype                                       ! 1: outside[normal case], 2: inside [flooded]
    real(8)              :: area                                        ! area of the wet surface [i.e. pi*(Rout**2-Rin**2)]
    real(8)              :: Ca                                          ! Inertia  coefficient
    real(8)              :: Cd                                          ! Drag     coefficient
    integer              :: ifroude                                     ! flag for including Froude-Krylov Force    in Morison's equ. [0:no, 1:include]
    integer              :: idyn_pres                                   ! flag for including dynamic pressure  term in Bernoulli equ. [0:no, 1:include]
    integer              :: irel                                        ! flag for including relative system's term in Bernoulli equ. [0:no, 1:include]
    real(8)              :: Output_res(7)                               ! results

 END type type_buoyancy_tap


 END Module mod_types_el
!
!
!
!----------------------------------------------------------------------
 Module mod_matrix_vars_el
!----------------------------------------------------------------------

    real(8), allocatable, save :: AM_el    (:,:)                        ! Total MASS          Matrix                 !<(NDFT_el,NDFT_el)
    real(8), allocatable, save :: AC_el    (:,:)                        ! Total DAMPING       Matrix                 !<(NDFT_el,NDFT_el)
    real(8), allocatable, save :: AK_el    (:,:)                        ! Total STIFFNESS     Matrix                 !<(NDFT_el,NDFT_el)
    real(8), allocatable, save :: AQ_el    (:  )                        ! Total RHS (LOAD)    Matrix                 !<(NDFT_el        )
    real(8), allocatable, save :: CDAMPSTR (:,:)                        ! Total Modal Damping Matrix                 !<(NDFT_el,NDFT_el)

                                                                        ! .:deformation, 1:velocity, 2:acceleration
    real(8), allocatable, save :: UT_el    (:), UT1_el (:), UT2_el (:)  ! Current  time/iteration solution           !<(NDFT_el)
    real(8), allocatable, save :: UTP_el   (:), UTP1_el(:), UTP2_el(:)  ! Previous time-step      solution           !<(NDFT_el)
    real(8), allocatable, save :: UT0_el   (:), UT01_el(:), UT02_el(:)  ! Current  iteration      perturbation       !<(NDFT_el)

    integer, allocatable, save :: INDSYSQ  (:)                          ! Flag for Q's reduction [0:no red. ,-1:don't include q at all, >0, reduction (== to the set elastic dof) !<(NQS    )
    integer, allocatable, save :: INDSYSB  (:)                          ! switch solution parameters (UT_el, UT1_el, UT2_el) to the reduced order matrices [NDFT0_el,NDFT0_el]    !<(NDFT_el)

 END Module mod_matrix_vars_el
!----------------------------------------------------------------------
 Module mod_transf_mat
!----------------------------------------------------------------------

 Use mod_types_el

!---- Transformation Matrices
    type (type_transf_mat ), allocatable, save :: transf_mat  (:)      ! [NNODT_el]
    real(8)                , allocatable, save :: ATP_el      (:,:,:)  ! [NBLADE_el, 3, 3] Blade c.s. after cone, before pitch, before pre-curved angles
!!  real(8)                , allocatable, save :: DATP_el     (:,:,:)  ! [NBLADE_el, 3, 3]
!!  real(8)                , allocatable, save :: DDATP_el    (:,:,:)  ! [NBLADE_el, 3, 3]
    real(8)                , allocatable, save :: ATPWr_el    (:,:,:)  ! [NBLADE_el, 3, 3] Blade c.s. after cone, after  pitch, before pre-curved angles - only for writeout
    real(8)                , allocatable, save :: ATPWrG_el   (:,:,:)  ! [NBLADE_el, 3, 3] Blade c.s. after cone, before pitch, before pre-curved angles - only for writeout
    real(8)                ,              save :: HUBpos_el   (    3)
!---- for the nacelle drag
    real(8)                ,              save :: Anacelle_el (  3,3)
!---- for the yaw actuator
    real(8)                ,              save :: ATyaw_el    (  3,3)
    real(8)                , allocatable, save :: AT0yaw_el   (:,:,:)  ! [NQS      , 3, 3]
!!  real(8)                ,              save :: Ayaw_el     (  3,3)
!!  real(8)                , allocatable, save :: A0yaw_el    (:,:,:)  ! [NQS      , 3, 3]
!---- for the rotating shaft
    real(8)                ,              save :: ATsftNoR_el (  3,3)
!---- for the floater
    real(8)                ,              save ::   AT_float  (3,3  )
    real(8)                ,              save ::    A_float  (3,3  )
    real(8)                ,              save ::   DA_float  (3,3  )
    real(8)                ,              save ::  DDA_float  (3,3  )
    real(8)                ,              save ::    R_float  (3    )
    real(8)                ,              save ::   DR_float  (3    )
    real(8)                ,              save ::  DDR_float  (3    )
    real(8)                ,              save ::   A0_float  (3,3,3)
    real(8)                ,              save ::  DA0_float  (3,3,3)
    real(8)                ,              save :: DDA0_float  (3,3,3)

 END Module mod_transf_mat
!
!
!
!----------------------------------------------------------------------
 Module Cbeam
!----------------------------------------------------------------------

 Use mod_types_el
 Use mod_matrix_vars_el
 Use mod_transf_mat


!--- Accumalative Variables
 type (type_ACCU)      ,              save :: ACCU

!--- Body Variables
 type (type_body   )   , allocatable, save :: body     (:)     ! NBODBT_el

!--- Sub-body Variables
 type (type_subbody)   , allocatable, save :: subbody  (:)     ! NBODT_el
 type (type_QSnode )   , allocatable, save :: QSnode   (:)     ! NNODT_el

!--- Parameters
 integer               , parameter         :: NNPE_el     =  2     !     number of nodes     per element
 integer               , parameter         :: NDFPEM_el   = 12     ! Max number of d.o.f.    per element
 integer               , parameter         :: NEQPEM_el   =  6     ! Max number of equations per element
 integer               , parameter         :: NBCPBODM_el =  2     ! Max number of Boundary  Conditions per Body (jacket)

!--- General Variables
 real(8)               ,              save :: PI
 real(8)               ,              save :: GRAV
 real(8)               ,              save :: R2D
 real(8)               ,              save :: PI2
 real(8)               ,              save :: PIhalf
 real(8)               ,              save :: RPM2RAD
 real(8)               ,              save :: RLXs
 real(8)               ,              save :: BIG

!--- Consentrated Masses
 integer               ,              save :: NCOMA_el
 type (type_conc_mass ), allocatable, save :: conc_mass (:)        ! [NCOMA_el]

!--- Beam Structural Properties
 integer               ,              save :: INODNELPROP_el
 integer               ,              save :: INODNELPROPbl_el
 integer               ,              save :: INODNELPROPWT_el
 integer               ,              save :: INODNELPROPJA_el
 type (type_beam_timo ), allocatable, save :: beam_timo  (:)       ! [INODNELPROP_el]
 type (type_substr_mor), allocatable, save :: substr_mor (:)       ! [INODNELPROP_el] -->  INODNELPROPsubstr_el

!--- Matrix integration using Gauss Integration method
 integer               , parameter         :: IORDER = 6
 real(8)               , allocatable, save :: XGAUSS(:), WGAUSS(:) ! [IORDER]

!--- Boundary Conditions
 type (type_boundco   ), allocatable, save :: boundco   (:)         ! [NBODT_el]

!--- Foundation
 integer               ,              save :: NBODTFND_el
 type (type_foundation), allocatable, save :: foundation(:)         ! [NBODTFND_el]

!--- Application Parameters
 integer               ,              save :: IAPPL_el         ! Application: 0: HAWT, 1: VAWT, 2:Helicopter
 integer               ,              save :: ICASE_el         ! WT type: 0: onshore, 1: monopile, 2: bottom based tripod/jacket, 3: rigid floater, 4: flexible floater
 integer               ,              save :: IYNmodal         ! solution 0: FEM    , 1: modal
 integer               ,              save :: IYNmodal_ut      ! switch to FEM solution[0] or modal solution[1] for elastic deflection, velocity, acceleration in sb LOCAL_UT_1
 integer               ,              save :: IREDUCT_el       ! reduction        flag: 0: no reduction, 1: matrix reduction enabled
 integer               ,              save :: IREDUCTJA_el     ! reduction jacket flag: 0: no reduction, 1: system equations reduction enabled
 integer               ,              save :: ITRUSS_el        ! truss element    flag: 0: not used, 1: extrernal force, 2: coupled with the structural part
 integer               ,              save :: IEIG             ! eigen-value      flag: 0: no (time domain), 1:eigen value without gravity, 2:eigen value including gravity
 integer               ,              save :: IMODAL           ! modal damping    flag: 0: no, 1:calculate modal damping, 2:read modal damping from file
 integer               ,              save :: IAERO_el         ! aerodynamic      flag: 0: no, 1: BEM, 2: GENUVP
 integer               ,              save :: IORDbl_el        ! beam order       flag: 1: 1st order Timoshenko, 2: 2nd order Euller-Bernoulli
 integer               ,              save :: ITYPE_FOUND_el   ! foundation       flag: 0: Apparently fixed AP, 1: Concentrated springs CS, 2: Distributed springs DS
 integer               ,              save :: ICMOD            ! controller       flag: 0: fixed omega and pitch, 1: variable speed-pitch control, 2: free rotor with fixed pitch
                                                               !                                              10: variable speed-pitch control built-in controller 
 integer               ,              save :: ICMOD_def        ! which controller  : 0: external conrtoller, 1: default built-in controller

 integer               , allocatable, save :: IWRITE_el (:,:)  ! flag for writing global/local signals (loads,deform) [NBODBT_el+1,2] 1:global, 2:local

!--- Control Vars
 real(8)               ,              save :: TGenLSS_el        ! Generator actual Torque at LSS
 real(8)               ,              save :: TGenLSS_D_el      ! Generator Damping term through linearization basically for stall regulated WT
 real(8)               ,              save :: TGenLSS_M_el      ! Generator Inertia term through linearization (atm DTD and Gen Lag)
 real(8)               ,              save :: TMLossLSS_el      ! Torque due to mechanichal loss
 real(8)               ,              save :: TBrakeLSS_el      ! Actual Braque torque
 real(8)               , allocatable, save :: PITCH_el (:), dPITCH_el (:), ddPITCH_el (:) ! Actual Pitch angle, velocitiy and acceleration [NBLADE_el]
 real(8)               , allocatable, save ::  FLAP_el (:),  dFLAP_el (:),  ddFLAP_el (:) ! Actual Pitch angle, velocitiy and acceleration [NBLADE_el]

!--- Initial Conditions
 integer               ,              save :: ICALCINIT        ! calculate initial state: 0: all deflections are set zero, 1: initial static calculation
 real(8)               , allocatable, save :: PHI0        (:)  ! Initial Azimuth   angle  for each Blade          [NBLADE_el]
 real(8)               , allocatable, save :: PITCH0_el   (:)  ! Initial Pitch     angle  for each Blade          [NBLADE_el]
 real(8)               , allocatable, save :: PITCHC_el   (:)  ! Initial Pitch Cosine Cyc for each Blade for Heli [NBLADE_el]
 real(8)               , allocatable, save :: PITCHS_el   (:)  ! Initial Pitch   Sine Cyc for each Blade for Heli [NBLADE_el]
 real(8)               , allocatable, save :: PITCH_Imb_el(:)  ! Pitch   Imbalance angle  for each Blade          [NBLADE_el]
 real(8)               ,              save :: OMEGA            ! Initial rotational speed       [LSS]
 real(8)               ,              save :: OMEGAR           ! Rated/Nomilan rotational speed [LSS] --> [transformed to rad/s]
 real(8)               ,              save :: TREF             ! Initial controller torque or B_GEN

!--- WT Characteristics
 real(8)               ,              save :: Htow             !
 real(8)               ,              save :: Htow0            !
 real(8)               ,              save :: Hsh              !
 real(8)               ,              save :: Hshof            !
 real(8)               ,              save :: Hhub             !
 integer               ,              save :: IAX_Ben_Swe      !
 real(8)               ,              save :: YAW              !
 real(8)               ,              save :: TILT             !
 real(8)               ,              save :: CONE             !
 real(8)               ,              save :: TPERIOD_el       !
 real(8)               ,              save :: RAT_GEAR         !
 integer               ,              save :: IACT_YAW         !
 real(8)               ,              save :: CSTIF_YAW        !
 real(8)               ,              save :: CDAMP_YAW        !

!--- Simulation Parameters
 real(8)               ,              save :: TIME             ! total time [current]
 real(8)               ,              save :: DT_el            ! time-step
 real(8)               ,              save :: BITA_el          ! parameter of NEWMARK integration shceme [β>=γ/2: conditionally stable]
 real(8)               ,              save :: GAMMA_el         ! parameter of NEWMARK integration shceme [γ=0.5: 2nd order, γ>0.5: numerical damping]
 integer               ,              save :: IRERUN           !
 integer               ,              save :: NTIMEBACK        !
 real(8)               ,              save :: MAXERR           !
 integer               ,              save :: ITERMAX          !
 real(8)               ,              save :: GRAVITY_el       !

!--- Damping Coefficients
 real(8)               , allocatable, save :: COEFM_el   (:)   ! Rayleigh   Damping mass      coefficient [NBODBT_el]
 real(8)               , allocatable, save :: COEFK_el   (:)   ! Rayleigh   Damping stiffness coefficient [NBODBT_el]
 real(8)               , allocatable, save :: CCRIT_el   (:)   ! Structural Damping coefficients          [0:NBODBTWT_el+1]
 real(8)               ,              save :: FREQ0_el         !

!--- Nacelle/Tower Drag Variables
 integer               ,              save :: IDRAG_TOW                                ! Index for activation (1) or not (0) the tower drag
 real(8)               ,              save :: CDRAG_TOW_el     , Z_DRAG_LIM            ! Tower Drag Coeff and Height of application (above it)
 integer               ,              save :: IDRAG_NAC                                ! Index for activation (1) or not (0) the nacelle drag
 real(8)               ,              save :: NAC_LX_el, NAC_LY_el, NAC_LZ_el          ! Nacelle length wrt global c.s.
 real(8)               ,              save :: NAC_RX_el, NAC_RY_el, NAC_RZ_el          ! Center of nacelle wrt tower top last beam point
 integer               ,              save :: NANG_NAC                                 ! number of angles for Cl,Cd
 real(8)               ,              save :: Alpha_NAC(900), CL_NAC(900), CD_NAC(900) ! angle, cl, cd

!--- Bodies/Sub-bodies/Elements
 integer               ,              save :: NTIMEM           !
 integer               ,              save :: NBODBT_el        ! Total Number of bodies
 integer               ,              save :: NBLADE_el        ! Total Number of blades
 integer               ,              save :: NSHAFT_el        ! Body  number of shaft
 integer               ,              save :: NTOWER_el        ! Body  number of tower
 integer               ,              save :: NJACKET_el       ! Body  number of jacket
 integer               ,              save :: NFLOATER_el      ! Body  number of floaters [dofs are 6*NFLOATER_el]

 integer               ,              save :: NBODT_el         ! Total Number of sub-bodies
 integer               ,              save :: NNODT_el         ! Total Number of nodes per sub-bodies
 integer               ,              save :: NDFT_el          ! Total Number of dofs before reduction
 integer               ,              save :: NDFT0_el         ! Total Number of dofs after  reduction
 integer               ,              save :: NDFBT_el         ! Total Number of elastic dofs before reduction
 integer               ,              save :: NDFBT0_el        ! Total Number of elastic dofs after  reduction
 integer               ,              save :: NDFBTWT0_el      ! Total Number of WT (blade, shaft, tower)     bodies 
 integer               ,              save :: NDFBTWT_el       ! Total Number of WT (blade, shaft, tower)     bodies 
 integer               ,              save :: NQS              ! Total Number of qs
 integer               ,              save :: NQSP             ! Number of primary qs [floater, omega, yaw]
 integer               ,              save :: NQSW             !
 integer               ,              save :: NQSY             !
 integer               ,              save :: IQfl
 integer               ,              save :: IQfl_rot
 integer               ,              save :: IQfl_tr
 integer               ,              save :: IMDOF_el(6)      !

 integer               ,              save :: NBODBTWT_el      ! Total Number of WT (blade, shaft, tower)     bodies 
 integer               ,              save :: NBODTWT_el       ! Total Number of WT (blade, shaft, tower) sub-bodies 
 integer               ,              save :: NNODTWT_el       !
 integer               ,              save :: NBODTJA_el       !
 integer               ,              save :: NBODBTJA_el      !
 integer               ,              save :: NSUBperBM_el     !
 integer               ,              save :: NODperSUBM_el    !
 integer               ,              save :: NELperSUBM_el    !
 integer               ,              save :: NELperBODYM_el   !

 real(8)               ,              save :: UThub_el         !
 real(8)               ,              save :: UThub1_el        !
 real(8)               ,              save :: UThub2_el        !
 real(8)               ,              save :: tilt_ttop        !
 real(8)               ,              save :: yaw_ttop

!--- Add the angular elastic deformations
 real(8)               , allocatable, save :: UTels2aer  (:,:) !(nb_el rotor, 3 )
 real(8)               , allocatable, save :: UT1els2aer (:,:) !(nb_el rotor, 3 )
 real(8)               , allocatable, save :: UT2els2aer (:,:) !(nb_el rotor, 3 )

!--- Jacket Buoyancy tapers
 integer                 ,              save :: NBuoyaTap_el
 type (type_buoyancy_tap), allocatable, save :: buoyancy_tap    (:) ! [NBuoyaTap_el]
 integer                 , allocatable, save :: IFloodbod_el    (:) ! [ NBODBT_el  ] indicates if the body is free flooded (1) or not (0).
 real(8)                 ,              save :: BuoyancyT_el    (6)
 real(8)                 ,              save :: BuoyancyTtap_el (6)
 real(8)                 ,              save :: BuoyancyTside_el(6)

!--- Dummy Variables for paper
 real(8)               ,              save :: Var_Paper(20)


 END Module Cbeam 

!----------------------------------------------------------------------
 Module coupling_cfd_0
!----------------------------------------------------------------------
 integer             , save :: NBODSUB_el_cfd, NTEB_el_cfd, NNTB_el_cfd
 integer, allocatable, save :: body_NBODSUB_el(:), body_NBODTGLB_el(:,:), & 
                                                   body_NNODTGLB_el(:,:) 
 integer, allocatable, save :: subbody_NTEB_el(:),subbody_NEQPE_el(:), &
                               subbody_NDFPE_el(:), &
                               subbody_NODMB(:,:,:),subbody_NDFPNBACC_el(:,:)
 integer             , save :: IMDOF_elcfd(6)
 real(8), allocatable, save :: subbody_HTA_el(:,:), subbody_ALENG_el(:,:), &
                               subbody_PHIX_el(:,:), subbody_PHIZ_el(:,:)

 integer, allocatable, save :: ACCU_NDFTACC_el(:) 
 real(8), allocatable, save :: transf_mat_R_el(:,:),   transf_mat_DR_el(:,:), &
                               transf_mat_A_el(:,:,:), transf_mat_AT_el(:,:,:), &
                               transf_mat_DA_el(:,:,:)
 END Module coupling_cfd_0
!----------------------------------------------------------------------
 Module coupling_cfd_f
!----------------------------------------------------------------------
 real(8), allocatable, save :: F_cfd(:,:,:)
 END Module coupling_cfd_f

!----------------------------------------------------------------------
 Module Paths
!----------------------------------------------------------------------
!-- Open files
 character*256, save :: file_dfile_el   !
 character*256, save :: file_rotorin    !
 character*256, save :: file_hydro      !
 character*256, save :: file_jacket     !
 character*256, save :: file_machi      !
 character*256, save :: file_geomb      !
 character*256, save :: file_profilb    !
 character*256, save :: dir_airfoils    !
 character*256, save :: file_drag       !
 character*256, save :: file_floater    !
 character*256, save :: file_diffract   !
 character*256, save :: file_retard     !
 character*256, save :: file_morison    !
 character*256, save :: file_truss      !
 character*256, save :: file_wind       !
 character*256, save :: dir_wave        !
 character*256, save :: file_foundation !
 character*256, save :: file_cstr       !
 character*256, save :: file_dllname    !
 character*256, save :: file_dllcfg     !

 END Module Paths
!----------------------------------------------------------------------
 Program hGAST
!----------------------------------------------------------------------

 Use Cbeam
#ifdef HAVE_MPI
 Use MPI
#endif
 
   implicit None

   integer :: ntime, ntime0
   integer :: my_rank,ierr

 ! call omp_set_num_threads(1)
 ! call mkl_set_num_threads(1)
 ! call mkl_domain_set_num_threads(1)

#ifdef HAVE_MPI
 call MPI_INIT     (                       ierr)
 call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)
#else
 ierr    = 0
 my_rank = 0
#endif

!--- Initialization
                    call INITIA_el

#ifdef HAVE_MPI
 call MPI_BCAST(IYNmodal,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

   if (IYNmodal==0) call GELAST0
 if (my_rank==0) then
#ifdef HAVE_MODAL
                    call INITIA_md
   if (IYNmodal==1) call GELAST0_md
#endif
                    call RECALL  (NTIME0)
 endif !my_rank

#ifdef HAVE_MPI
 call MPI_BCAST(ntime0  ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
 call MPI_BCAST(NTIMEM  ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
 call MPI_BCAST(TIME    ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
 call MPI_BCAST(DT_el   ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
 call MPI_BCAST(ITERMAX ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
#endif

!--- for all time steps
   do ntime = ntime0, NTIMEM
      TIME  = TIME + DT_el 

#ifdef HAVE_MPI
 call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif

      if (IYNmodal==0) call GELAST    (ntime) 
#ifdef HAVE_MODAL
      if (IYNmodal==1) call GELAST_md (ntime) 
#endif
      if (my_rank ==0 &
     .and.IYNmodal==0) call RERUN     (ntime)
   enddo

!--- Finalization
 if (my_rank==0) &
   call FINIT_el

#ifdef HAVE_MPI
 call MPI_BARRIER (MPI_COMM_WORLD,ierr)
 call MPI_FINALIZE(ierr)
#endif


 END Program hGAST
