!----------------------------------------------------------------------
 module param_tr
!----------------------------------------------------------------------

 integer, parameter         :: NNPE_tr      = 2
 integer, parameter         :: NDFPE_tr     = 6
 integer, parameter         :: NDFPN_tr     = 3
 integer, parameter         :: NEQPE_tr     = 3

 END module param_tr
!----------------------------------------------------------------------
 module mod_types_tr
!----------------------------------------------------------------------
 Use param_tr

!----------------------------------------------------------------------
 type :: type_body_tr     ! NBODT_tr, used in <body_tr>
!----------------------------------------------------------------------

!---- Bodies (truss elements)
    integer :: INODE_tr  (2)        !
    real(8) :: ALENG_tr             !
    real(8) :: ADIAM_tr             !
    real(8) :: DENS_tr              !
    real(8) :: WaterFz_tr           !
    real(8) :: EAS1_tr              !
    real(8) :: EAS2_tr              !
    real(8) :: CDAMP_tr             !
    real(8) :: Ca_tr                !
    real(8) :: Cdnorm_tr            !
    real(8) :: Cdtang_tr            !
    real(8) :: Cdfric_tr            !

    real(8) :: AMLOC_tr  (NDFPE_tr,NDFPE_tr)
    real(8) :: ACLOC_tr  (NDFPE_tr,NDFPE_tr)
    real(8) :: AKLOC_tr  (NDFPE_tr,NDFPE_tr)
    real(8) :: AQLOC_tr  (         NDFPE_tr)

 END type type_body_tr
!----------------------------------------------------------------------
 type :: type_constrain_tr  ! NCONSTR_tr, used in <constrain_tr>
!----------------------------------------------------------------------

!---- Constrains per node
    integer :: INODC_tr 
    integer :: INDXC_tr  (NDFPN_tr)

 END type type_constrain_tr
!----------------------------------------------------------------------
 type :: type_connect_tr  ! NCONNECT_tr, used in <connect_tr>
!----------------------------------------------------------------------
                               
!---- Connections with the floater
    integer :: IBODCONNECT_tr             ! which truss element body
    integer :: INODCONNECT_tr             ! which truss element node
    real(8) :: ALOADB         (3)         ! Global  loads of the connection
    real(8) :: RLOC_CONNECT_tr(3)         ! Global   distance            of the connection wrt reference point when all qs are zero
    real(8) ::   UT_CONNECT_tr(3)         ! current  Global position     of the connection node
    real(8) ::  UT1_CONNECT_tr(3)         ! current  Global velocity     of the connection node
    real(8) ::  UT2_CONNECT_tr(3)         ! current  Global acceleration of the connection node
    real(8) ::  UTP_CONNECT_tr(3)         ! previous Global position     of the connection node
    real(8) :: UTP1_CONNECT_tr(3)         ! previous Global velocity     of the connection node
    real(8) :: UTP2_CONNECT_tr(3)         ! previous Global acceleration of the connection node
!---- for ICASE_el = 4
    integer :: nbod_ja                    ! which elastic(beam) body
    integer :: nel_ja                     ! which elastic(beam) nel
    real(8) :: RLOC_ja(3)                 ! local distance between the connection and the elastic(beam) body - wrt the start of the element

 END type type_connect_tr
!----------------------------------------------------------------------
 type :: type_conc_mass_tr ! NCON_tr, used in <conc_mass_tr>
!----------------------------------------------------------------------
                               
!---- Concentrated masses (Buoys)
    integer :: IBCON_tr
    integer :: INCON_tr
    real(8) :: AMASSCON_tr
    real(8) :: DIAMCON_tr

 END type type_conc_mass_tr
!
!
!
 END module mod_types_tr
!----------------------------------------------------------------------
 module truss
!----------------------------------------------------------------------
 
 Use param_tr
 Use mod_types_tr
                               
!----------------------------------------------------------------------
!-- Body Variables
 integer                              , save :: NBODT_tr
 type (type_body_tr )    , allocatable, save :: body_tr      (:)
!----------------------------------------------------------------------
!-- Constrain Variables
 integer                              , save :: NCONSTR_tr
 type (type_constrain_tr), allocatable, save :: constrain_tr (:)
!----------------------------------------------------------------------
!-- Connection Variables
 integer                              , save :: NCONNECT_tr
 type (type_connect_tr)  , allocatable, save :: connect_tr   (:)
!----------------------------------------------------------------------
!-- Concentrated Mass Variables
 integer                              , save :: NCON_tr
 type (type_conc_mass_tr), allocatable, save :: conc_mass_tr (:)
!----------------------------------------------------------------------
!-- General Variables/Parameters
 real(8), save              :: pi_tr
 real(8), save              :: GRAV_tr
 real(8), save              :: AROW_tr
 real(8), save              :: AMVISCW_tr
 real(8), save              :: Depth_tr
 real(8), save              :: DT_tr
 real(8), save              :: TIME_tr
 real(8), save              :: BITA_tr
 real(8), save              :: GAMMA_tr
 real(8), save              :: MAXERR_tr
 real(8), save              :: COEFM_tr
 real(8), save              :: COEFK_tr

 integer, save              :: NDFT_tr
 integer, save              :: NDFT0_tr
 integer, save              :: NNODE_tr
 integer, save              :: NSubTime_tr
 integer, save              :: ITERMAX_tr
 integer, save              :: IYN_Morison_tr
 integer, save              :: IYN_SeabedInt_tr
 integer, save              :: IREDUCT_tr
 integer, save              :: IWRITE_truss         !-- 1 output elem
!----------------------------------------------------------------------
!-- Global truss matrices
 real(8), allocatable, save :: AM_tr          (:,:  )      !(NDFT_tr,NDFT_tr)
 real(8), allocatable, save :: AC_tr          (:,:  )      !(NDFT_tr,NDFT_tr)
 real(8), allocatable, save :: AK_tr          (:,:  )      !(NDFT_tr,NDFT_tr)
 real(8), allocatable, save :: AQ_tr          (:    )      !(NDFT_tr        )
 real(8), allocatable, save :: UT_tr          (:    )      !(NDFT_tr        )
 real(8), allocatable, save :: UT1_tr         (:    )      !(NDFT_tr        )
 real(8), allocatable, save :: UT2_tr         (:    )      !(NDFT_tr        )
 real(8), allocatable, save :: UTP_tr         (:    )      !(NDFT_tr        )
 real(8), allocatable, save :: UTP1_tr        (:    )      !(NDFT_tr        )
 real(8), allocatable, save :: UTP2_tr        (:    )      !(NDFT_tr        )
 real(8), allocatable, save :: UTPP_tr        (:    )      !(NDFT_tr        )
 real(8), allocatable, save :: UTPP1_tr       (:    )      !(NDFT_tr        )
 real(8), allocatable, save :: UTPP2_tr       (:    )      !(NDFT_tr        )
 integer, allocatable, save :: INDSYSB_tr     (:    )      !(NDFT_tr        )
!----------------------------------------------------------------------
!-- Dof accumulative and Initial Global position of each node 
 integer, allocatable, save :: NDFACC_tr      (:    )      !(0:NNODE_tr     )
 real(8), allocatable, save :: XG_tr          (:,:  )      !(3,NNODE_tr     )
!----------------------------------------------------------------------


END module truss

!1. check initial solution                                --->ok
!2. check the acceleration when substeps are introdused   --->ok
!3. lower yaw damping                                     --->ok [fixed with morison length]
!4. matrix reduction??                                    --->ok
!5. structs                                               --->ok 
!6. eigen-value coupled
!7. loads at jacket members when ICASE_el = 4             --->ok

!-- ask about DT and damping ratio for eigen value
