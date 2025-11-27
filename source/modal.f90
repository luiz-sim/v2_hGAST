! 1. subbodies - prebend :: [modal anal + sb_loads_md {expressions without RG}]
! 2. modal damping
! 3. bod2bod            OK

! read modes            OK
! geload_md             OK
! floatloads_md         OK
!
!
!
! REMEMBER
!==========
! 1. re-enable bod2bod  loads GAST [eigen]
! 2. re-enable traction force GAST [math] / stifloc_tract [vgale to return]
!----------------------------------------------------------------------
! TO DO
!========
! 1. modes                                                                 OK [1.subbodies consideration in modal anal new. 2.modal damping]
! 2. prepro                                                                OK [can't incorporate Fy]
! 3. matrix_md/matrloc_md                                                  OK
! 4. initia                                                                OK
! 5. qs_equate                                                             OK [only the modal part of the q assignment]
! 6. matrix_reduction                                                      OK [for static cases                       ]
! 7. concentrated mass/inertia                                             OK
! 8. write out [modes consideration]                                       OK [still writeout loads]
! 9. keep old U for writing purpose  [so writout is streight forward]      OK

! 1. loads application: Int{F^t N^t IIae Nae}                              OK
! 2. bod2bod loads [calc loads for transfer, Fy]                           OK [reaction force from FEM]
! 3. GENLOAD_md                                                            OK
! 4. check with voutsi F^t * F                                             OK
! 5. read kept modes / enable coupling / identification etc...             OK [only chose kept modes]
! 6. modal damping
!
!----------------------------------------------------------------------
 Module modal_mod
!----------------------------------------------------------------------
 type :: type_subbody_md   ! NBODT_el, used in 
!----------------------------------------------------------------------

!---- Subbody_md
    integer              :: Nmode_TR_md                        ! total modes
    integer              :: itor_coupl                         ! switch for torsion   coupling  [0:disable, 1:enable]
    integer              :: iext_coupl                         ! switch for extension coupling  [0:disable, 1:enable]
    integer, allocatable :: mode                  (:      )    !                          (Nmode          )
    real(8), allocatable :: FTT_TR_md             (:,:    )    ! truncated modal matrix   (ndftb ,    Nmode)
    real(8), allocatable :: FTT_TR_T_md           (:,:    )    ! transpose  >>     >>     (Nmode ,    ndftb)
!---- Integrals from Prepro
    real(8), allocatable :: FT_Int_K_F_bod        (:,    :)    !                          (Nmode ,    Nmode)
    real(8), allocatable :: FT_Int_NTII33SN_F_bod (:,:,:,:)    !                          (Nmode ,3,3,Nmode)
    real(8), allocatable :: FT_Int_NTII3_bod      (:,:    )    !                          (Nmode ,3   [1]  )
    real(8), allocatable :: FT_Int_NTII33r0_bod   (:,:,:  )    !                          (Nmode ,3,3,[1]  )
!   real(8), allocatable :: FT_Int_NTIIaeNae_bod  (:,    :)    !                          (Nmode ,   ,6xNstrips)

!   real(8), allocatable :: Fy                    (:,:  )      !                          (nelb  ,gauss)
!   real(8), allocatable :: Load                  (:,:  )      !                          (nelb+1      )

    real(8), allocatable :: LoadIntQ_L            (:,:  )      !  (6, neltb+1  )
    real(8), allocatable :: LoadIntQ              (:    )      !  (6           )
    real(8), allocatable :: LoadIntM              (:,:  )      !  (6, ndftb    )
    real(8), allocatable :: LoadIntC              (:,:  )      !  (6, ndftb    )
    real(8), allocatable :: LoadIntK              (:,:  )      !  (6, ndftb    )
    real(8), allocatable :: LoadIntMq             (:,:  )      !  (6, NQS      )
    real(8), allocatable :: LoadIntCq             (:,:  )      !  (6, NQS      )
    real(8), allocatable :: LoadIntKq             (:,:  )      !  (6, NQS      )


 END type type_subbody_md

 type (type_subbody_md), allocatable, save :: subbody_md (:)   ! NBODT_el

 integer,                parameter         :: NMODEM_md   = 300 ! subbody Max number of d.o.f.[modes]
 integer,                allocatable, save :: NDFTACC_md  (:)   ! subbody accumulated dof for modal solution                      ! [0:NBODT_el]
 integer,                allocatable, save :: INDSYSQ_bod (:)   ! subbody correspondance of a dependent Q, need for reduction     ! [  NQS     ]
 real(8),                allocatable, save :: UT_md       (:)   ! reconstruction of FEM solution vector, based on modal solution  ! [  NDFT_el ] : FEM
 real(8),                allocatable, save :: UT1_md      (:)
 real(8),                allocatable, save :: UT2_md      (:)
 integer,                             save :: Integr_md         !Loads calculation flag 0:reaction, 1:integration

 END Module modal_mod
!----------------------------------------------------------------------
! 
!     Subroutine :INITIA_md
!     
!     Open input FILES, perform initial Calculations
!
!----------------------------------------------------------------------
 Subroutine INITIA_md
!----------------------------------------------------------------------

 use Cbeam
 use modal_mod

   implicit none

!--- dfile vars
   integer           :: Nmode, NDFTB
   integer           :: nbod , nbsub, nb_el
   integer           :: mode(NMODEM_md), ir, i
   integer           :: neltb
   integer           :: itor, iext


   if (IYNmodal==0) return


   open (80,file='modal.inp')
   read (80,*)
   read (80,*) Integr_md         !Loads calculation flag 0:reaction, 1:integration
   write(*,*)
   write(*,*)
   write(*,*) 'Initialization of the modal part'


      Allocate ( subbody_md     (NBODT_el ))
      Allocate ( NDFTACC_md   (0:NBODT_el ))
                 NDFTACC_md (0) = 0

   do nbod  = 1, NBODBTWT_el
   do nbsub = 1, body(nbod)%NBODSUB_el
      nb_el = body(nbod)%NBODTGLB_el(nbsub)
!!    NELtb = 1 !subbody   (nb_el)%NTEB_el
      ndftb = subbody   (nb_el)%NDFTB_el
!     NDFPE = subbody   (nb_el)%NDFPE_el
!     Nmode = subbody_md(nb_el)%Nmode_TR_md
      neltb = subbody   (nb_el)%NTEB_el

      read (80,*) itor, iext, Nmode, ir, (mode(i),i=1,Nmode*ir)

      if (ir/=0.and.ir/=1) then
         write(*,*)'wrong ir input in modal.inp, ir must be 0 or 1'
         stop
      endif
      if ( (itor/=0.and.itor/=1).or.&
           (iext/=0.and.iext/=1)      ) then
         write(*,*)'wrong flag for torsion or extension coupling consideration in modal.inp, must be 0 or 1'
         stop
      endif
      if (nbod==NBLADE_el+1.and.itor==0) then
         write(*,*)'Error. Torsion dof is disabled from shaft in modal.inp'
         stop
      endif

      Allocate ( subbody_md(nb_el)%mode                 (Nmode          ) )

                    subbody_md(nb_el)%Nmode_TR_md = Nmode
                    subbody_md(nb_el)%itor_coupl  = itor   ! switch for torsion   coupling  [0:disable, 1:enable]
                    subbody_md(nb_el)%iext_coupl  = iext   ! switch for extension coupling  [0:disable, 1:enable]
      do i=1,Nmode
         if (ir==0) subbody_md(nb_el)%mode (i)=i
         if (ir==1) subbody_md(nb_el)%mode (i)=mode(i)
      enddo

      Allocate ( subbody_md(nb_el)%FT_Int_K_F_bod       (Nmode    ,Nmode) )
      Allocate ( subbody_md(nb_el)%FT_Int_NTII33SN_F_bod(Nmode,3,3,Nmode) )
      Allocate ( subbody_md(nb_el)%FT_Int_NTII3_bod     (Nmode,3        ) )
      Allocate ( subbody_md(nb_el)%FT_Int_NTII33r0_bod  (Nmode,3,3      ) )
      Allocate ( subbody_md(nb_el)%FTT_TR_md            (Ndftb,    Nmode) )
      Allocate ( subbody_md(nb_el)%FTT_TR_T_md          (Nmode,    ndftb) )

     if (Integr_md==1) then
      Allocate ( subbody_md(nb_el)%LoadIntQ_L           (6, neltb+1 ) )
      Allocate ( subbody_md(nb_el)%LoadIntQ             (6          ) )
      Allocate ( subbody_md(nb_el)%LoadIntM             (6, ndftb   ) )
      Allocate ( subbody_md(nb_el)%LoadIntC             (6, ndftb   ) )
      Allocate ( subbody_md(nb_el)%LoadIntK             (6, ndftb   ) )
      Allocate ( subbody_md(nb_el)%LoadIntMq            (6, NQS     ) )
      Allocate ( subbody_md(nb_el)%LoadIntCq            (6, NQS     ) )
      Allocate ( subbody_md(nb_el)%LoadIntKq            (6, NQS     ) )
     endif
                NDFTACC_md(nb_el) = NDFTACC_md(nb_el-1) + subbody_md(nb_el)%Nmode_TR_md

      write(*,'(10i5)') nbod,nbsub,nb_el, Nmode, NDFTACC_md (nb_el)

      if ( Nmode > NMODEM_md ) then
         write(*,*) 'Increase the Modes per body parameter'
         stop
      endif
   enddo !nbsub
   enddo !nbod

!--- find maximum dimensions for main allocations
   Nmode = maxval ( subbody_md (1:NBODT_el)%Nmode_TR_md )
   NDFTB = maxval ( subbody    (1:NBODT_el)%NDFTB_el    )

   write(*,*)'max modes',Nmode
   write(*,*)'max ndftb',NDFTB

!--- Prepare the modal damping matrix, which in calculated inside GELAST0_MODAL [MODAL_ANAL_md]
   Deallocate ( CDAMPSTR )
   Allocate   ( CDAMPSTR (NDFTACC_md(NBODT_el)+NQS, NDFTACC_md(NBODT_el)+NQS) )

!--- Form the modal matrix and it's transpose, also calculate modal damping
   call GELAST0_MODAL

!--- Calculate the integrals FtNtII33SNF......
!qqcall MATRIX_prepro_md


!--- Set vars based on modal dimensions
   Deallocate ( AM_el , AC_el  , AK_el  , AQ_el   , &
                UT_el , UTP_el , UT0_el , INDSYSB , &
                UT1_el, UTP1_el, UT01_el, INDSYSQ , &
                UT2_el, UTP2_el, UT02_el              )

   Allocate   (  UT_md    (NDFT_el) ) ! Current iteration deformation  from modal to fem dofs
   Allocate   ( UT1_md    (NDFT_el) ) ! Current iteration velocity     from modal to fem dofs
   Allocate   ( UT2_md    (NDFT_el) ) ! Current iteration acceleration from modal to fem dofs

!--- change dimensions
   NDFBT_el    = NDFTACC_md  (NBODT_el)
   NDFT_el     = NDFBT_el + NQS
   IYNmodal_ut = 1           ! switch to FEM solution[0] or modal solution[1] for elastic deflection, velocity, acceleration in sb LOCAL_UT_1

   write(*,*) 'NDFBT_el',NDFBT_el
   write(*,*) 'NDFT_el ',NDFT_el


!--- Allocate main matrices
   Allocate ( AM_el       (NDFT_el,NDFT_el) )
   Allocate ( AC_el       (NDFT_el,NDFT_el) )
   Allocate ( AK_el       (NDFT_el,NDFT_el) )
   Allocate ( AQ_el       (NDFT_el        ) )
   Allocate ( INDSYSB     (NDFT_el        ) )
   Allocate ( INDSYSQ     (NQS            ) )
   Allocate ( INDSYSQ_bod (NQS            ) ) !new
   Allocate ( UT_el       (NDFT_el), UT1_el (NDFT_el), UT2_el (NDFT_el) ) ! Current iteration deformation, velocity, acceleration
   Allocate ( UTP_el      (NDFT_el), UTP1_el(NDFT_el), UTP2_el(NDFT_el) ) ! Previous Time Step
   Allocate ( UT0_el      (NDFT_el), UT01_el(NDFT_el), UT02_el(NDFT_el) ) ! Perturbation


 END Subroutine INITIA_md
!-----------------------------------------------------------------------------------------------------
! -- Subroutine :MATRIX_md
!
!   AM_el , AK_el, AC_el, AQ_el 
!   mass, stiffness damping and forcing matrices for the complete structure
!----------------------------------------------------------------------
 Subroutine MATRIX_md
!----------------------------------------------------------------------

 use Cbeam
 use modal_mod

   implicit none

   real(8) :: AMLOC    (NMODEM_md,NMODEM_md)
   real(8) :: ACLOC    (NMODEM_md,NMODEM_md)
   real(8) :: AKLOC    (NMODEM_md,NMODEM_md)
   real(8) :: AFLOC    (NMODEM_md          )
   real(8) :: AMLOCNQ  (NMODEM_md,NQS      )
   real(8) :: ACLOCNQ  (NMODEM_md,NQS      )
   real(8) :: AKLOCNQ  (NMODEM_md,NQS      )
   real(8) :: AMLOC0   (NDFPEM_el,NDFPEM_el)!--- new
   real(8) :: ACLOC0   (NDFPEM_el,NDFPEM_el)!--- new
   real(8) :: AKLOC0   (NDFPEM_el,NDFPEM_el)
   real(8) :: AFLOC0   (NDFPEM_el          )
   real(8) :: AMLOCNQ0 (NDFPEM_el,NQS      )!--- new
   real(8) :: ACLOCNQ0 (NDFPEM_el,NQS      )!--- new
   real(8) :: AKLOCNQ0 (NDFPEM_el,NQS      )
   integer :: nbod, nbsub, nb_el, nn_el, nel, i, j, IQRayl, IMAX, JMAX, nn, iq,nq
   integer :: NDFPE, NDFTB, Nmode
   real(8) :: COEFM, COEFK
!! real(8), allocatable :: Int_NtIIaeFae_bod (:), Int_NtKfy_N_bod  (:,:)
   real(8), allocatable :: AM_sb   (:,:)
   real(8), allocatable :: AC_sb   (:,:)
   real(8), allocatable :: AK_sb   (:,:)
   real(8), allocatable :: AF_sb   (:  )
   real(8), allocatable :: AMNQ_sb (:,:)
   real(8), allocatable :: ACNQ_sb (:,:)
   real(8), allocatable :: AKNQ_sb (:,:)
   real(8), allocatable :: ftt_t   (:,:)
   real(8), allocatable :: ftt     (:,:)
   integer :: ncm, nbod_cm, nbsub_cm, nel_cm


!--- Zero the global matrices
   AM_el  = 0.d0;
   AC_el  = 0.d0;
   AK_el  = 0.d0;
   AQ_el  = 0.d0;

!--- Zero local elements matrices
   do nb_el = 1, NBODT_el
      subbody   (nb_el)%AMLOC_el   ( :, :, : ) = 0.d0; !nel, ndfpe,ndfpe
      subbody   (nb_el)%ACLOC_el   ( :, :, : ) = 0.d0; !nel, ndfpe,ndfpe
      subbody   (nb_el)%AKLOC_el   ( :, :, : ) = 0.d0; !nel, ndfpe,ndfpe
      subbody   (nb_el)%AFLOC_el   ( :, :    ) = 0.d0; !nel, ndfpe
      if (nb_el>NBODTWT_el.and.ICASE_el==2) cycle
      subbody   (nb_el)%AMLOCNQ_el ( :, :, : ) = 0.d0; !nel, ndfpe,nqs
      subbody   (nb_el)%ACLOCNQ_el ( :, :, : ) = 0.d0; !nel, ndfpe,nqs
      subbody   (nb_el)%AKLOCNQ_el ( :, :, : ) = 0.d0; !nel, ndfpe,nqs
   enddo

   IQRayl    = 1 !0:all qs Rayleigh, 1:body qs with Rayleigh, 2:no Rayleigh to qs
   do nbod   = 1, NBODBT_el
   do nbsub  = body(nbod)%NBODSUB_el, 1, -1
      nb_el  = body(nbod)%NBODTGLB_el(nbsub)
      nn_el  = body(nbod)%NNODTGLB_el(nbsub)
      NDFPE  = subbody   (nb_el)%NDFPE_el
      NDFTB  = subbody   (nb_el)%NDFTB_el
      Nmode  = subbody_md(nb_el)%Nmode_TR_md

      AMLOC   = 0.d0; ACLOC   = 0.d0; AKLOC   = 0.d0; AFLOC = 0.d0;
      AMLOCNQ = 0.d0; ACLOCNQ = 0.d0; AKLOCNQ = 0.d0;

!!    call MATRLOC_md ( AMLOC  , ACLOC  , AKLOC  , AFLOC, &
!!                      AMLOCNQ, ACLOCNQ, AKLOCNQ,        &
!!                      nb_el  , nn_el  , Nmode             )

!!    Allocate ( Int_NtIIaeFae_bod (NDFTB), Int_NtKfy_N_bod (NDFTB,NDFTB) )
!!               Int_NtIIaeFae_bod = 0.d0;  Int_NtKfy_N_bod = 0.d0;
      Allocate ( AM_sb   (NDFTB,NDFTB) );  AM_sb   = 0.d0;
      Allocate ( AC_sb   (NDFTB,NDFTB) );  AC_sb   = 0.d0;
      Allocate ( AK_sb   (NDFTB,NDFTB) );  AK_sb   = 0.d0;
      Allocate ( AF_sb   (NDFTB      ) );  AF_sb   = 0.d0;
      Allocate ( AMNQ_sb (NDFTB,NQS  ) );  AMNQ_sb = 0.d0;
      Allocate ( ACNQ_sb (NDFTB,NQS  ) );  ACNQ_sb = 0.d0;
      Allocate ( AKNQ_sb (NDFTB,NQS  ) );  AKNQ_sb = 0.d0;

!------ Virtual works which depend upon time
      do nel = subbody(nb_el)%NTEB_el,1,-1

         call MATRLOC          ( AMLOC0  , ACLOC0  , AKLOC0  , AFLOC0     ,&
                                 AMLOCNQ0, ACLOCNQ0, AKLOCNQ0             ,&
                                                     nb_el   , nn_el , nel   )
!!       AKLOC0 = 0.d0; AFLOC0 = 0.d0; AKLOCNQ0 = 0.d0;
         call STIFLOC_tract_md ( AKLOC0  , AFLOC0  , nb_el   ,         nel   )

         call FORCLOC          ( AFLOC0  , AKLOCNQ0,                       &
                                 nbod    , nbsub   , nb_el   , nn_el , nel   )


         do ncm      = 1, NCOMA_el
            nbod_cm  = conc_mass(ncm)%NBCOMA_el
            nbsub_cm = conc_mass(ncm)%NSBCOMA_el
            nel_cm   = conc_mass(ncm)%NELCOMA_el

            if (nbod_cm==nbod.and.nbsub_cm==nbsub.and.nel_cm==nel)  &

            call CON_MATRLOC ( AMLOC0  , ACLOC0  , AKLOC0  , AFLOC0,&
                               AMLOCNQ0, ACLOCNQ0, AKLOCNQ0,        &
                               ncm     , nb_el   , nn_el   , nel      )
         enddo !ncm


         call STORE_LOC_MATRIX_el ( AMLOC0  , ACLOC0  , AKLOC0  , AFLOC0,&
                                    AMLOCNQ0, ACLOCNQ0, AKLOCNQ0,        &
                                    nb_el   , nn_el   , nel                )

!--------- assembly FEM
         nn       = subbody(nb_el)%NODMB       (nel,1)!nnod)
         i        = subbody(nb_el)%NDFPNBACC_el(nn-1 )
         j        = i
!!       Int_NtKfy_N_bod   (i+1:i+NDFPE,j+1:j+NDFPE) = Int_NtKfy_N_bod   (i+1:i+NDFPE,j+1:j+NDFPE) + AKLOC0 (1:NDFPE,1:NDFPE)
!!       Int_NtIIaeFae_bod (i+1:i+NDFPE            ) = Int_NtIIaeFae_bod (i+1:i+NDFPE            ) + AFLOC0 (1:NDFPE        )

         AM_sb   (i+1:i+NDFPE,j+1:j+NDFPE) = AM_sb   (i+1:i+NDFPE,j+1:j+NDFPE) + AMLOC0   (1:NDFPE,1:NDFPE)
         AC_sb   (i+1:i+NDFPE,j+1:j+NDFPE) = AC_sb   (i+1:i+NDFPE,j+1:j+NDFPE) + ACLOC0   (1:NDFPE,1:NDFPE)
         AK_sb   (i+1:i+NDFPE,j+1:j+NDFPE) = AK_sb   (i+1:i+NDFPE,j+1:j+NDFPE) + AKLOC0   (1:NDFPE,1:NDFPE)
         AF_sb   (i+1:i+NDFPE            ) = AF_sb   (i+1:i+NDFPE            ) + AFLOC0   (1:NDFPE        )
        do iq = 1, QSnode(nn_el)%Tot_el(2)
           nq = QSnode(nn_el)%IORD_el(iq)
         AMNQ_sb (i+1:i+NDFPE,         nq) = AMNQ_sb (i+1:i+NDFPE,         nq) + AMLOCNQ0 (1:NDFPE,     nq)
         ACNQ_sb (i+1:i+NDFPE,         nq) = ACNQ_sb (i+1:i+NDFPE,         nq) + ACLOCNQ0 (1:NDFPE,     nq)
         AKNQ_sb (i+1:i+NDFPE,         nq) = AKNQ_sb (i+1:i+NDFPE,         nq) + AKLOCNQ0 (1:NDFPE,     nq)
        enddo !iq
      enddo !nel


!------ modal trunctation
      Allocate ( ftt      (NDFTB,Nmode) )
      Allocate ( ftt_t    (Nmode,NDFTB) )

      ftt  (1:NDFTB,1:Nmode) = subbody_md(nb_el)%FTT_TR_md   (1:NDFTB,1:Nmode)
      ftt_t(1:Nmode,1:NDFTB) = subbody_md(nb_el)%FTT_TR_T_md (1:Nmode,1:NDFTB)

!!    AKLOC(1:Nmode,1:Nmode) = AKLOC(1:Nmode,1:Nmode) + matmul( matmul( ftt_t (1:Nmode,1:NDFTB), Int_NtKfy_N_bod   (1:NDFTB,1:NDFTB) ), ftt (1:NDFTB,1:Nmode) )
!!    AFLOC(1:Nmode        ) = AFLOC(1:Nmode        ) +         matmul( ftt_t (1:Nmode,1:NDFTB), Int_NtIIaeFae_bod (1:NDFTB        ) )

      AMLOC  (1:Nmode,1:Nmode) = matmul( ftt_t (1:Nmode,1:NDFTB), matmul( AM_sb  (1:NDFTB,1:NDFTB) , ftt (1:NDFTB,1:Nmode) ) )
      ACLOC  (1:Nmode,1:Nmode) = matmul( ftt_t (1:Nmode,1:NDFTB), matmul( AC_sb  (1:NDFTB,1:NDFTB) , ftt (1:NDFTB,1:Nmode) ) )
      AKLOC  (1:Nmode,1:Nmode) = matmul( ftt_t (1:Nmode,1:NDFTB), matmul( AK_sb  (1:NDFTB,1:NDFTB) , ftt (1:NDFTB,1:Nmode) ) )
      AFLOC  (1:Nmode        ) = matmul( ftt_t (1:Nmode,1:NDFTB),         AF_sb  (1:NDFTB        ) )
     do iq = 1, QSnode(nn_el)%Tot_el(2)
        nq = QSnode(nn_el)%IORD_el(iq)
      AMLOCNQ(1:Nmode,     nq) = matmul( ftt_t (1:Nmode,1:NDFTB),         AMNQ_sb(1:NDFTB,     nq) )
      ACLOCNQ(1:Nmode,     nq) = matmul( ftt_t (1:Nmode,1:NDFTB),         ACNQ_sb(1:NDFTB,     nq) )
      AKLOCNQ(1:Nmode,     nq) = matmul( ftt_t (1:Nmode,1:NDFTB),         AKNQ_sb(1:NDFTB,     nq) )
     enddo


!------ assembly modal
      IMAX  = Nmode
      JMAX  = Nmode
      i     = NDFTACC_md(nb_el-1)
      j     = i
      COEFM = COEFM_el(nbod)
      COEFK = COEFK_el(nbod)

      call ASSEMBLY_MATRIX_el ( AMLOC  , ACLOC  , AKLOC    , AFLOC         ,&
                                AMLOCNQ, ACLOCNQ, AKLOCNQ  , COEFM , COEFK ,&
                                         nn_el  , i        , j     , IQRayl,&
                                IMAX   , JMAX   , NMODEM_md, NMODEM_md        )

!------ loads communication between bodies
      call B_LOADS_md ( nbod, nbsub )


!!    Deallocate ( Int_NtIIaeFae_bod, Int_NtKfy_N_bod )
      Deallocate ( AM_sb  , AC_sb  , AK_sb  , AF_sb )
      Deallocate ( AMNQ_sb, ACNQ_sb, AKNQ_sb        )
      Deallocate ( ftt    , ftt_t                   )
   enddo !nbsub
   enddo !nbod

!----- change local element force vectors for output
!  if (Integr_md==1) then
!     do nbod   = 1, NBODBT_el
!     do nbsub  = 1, body(nbod)%NBODSUB_el
!        nb_el  = body(nbod)%NBODTGLB_el(nbsub)
!     do nel    = 1, subbody(nb_el)%NTEB_el
!        subbody(nb_el)%AFLOC_el ( nel, 1:6  ) = subbody_md(nb_el)%LoadIntQ_L ( 1:6, nel   )
!qq12x12
!        subbody(nb_el)%AFLOC_el ( nel,10:15 ) = subbody_md(nb_el)%LoadIntQ_L ( 1:6, nel+1 )
!     enddo !nel
!     enddo !nbsub
!     enddo !nbod
!  endif


 END Subroutine MATRIX_md
!----------------------------------------------------------------------
 Subroutine modal2fem 
!----------------------------------------------------------------------

 use Cbeam
 use modal_mod

   implicit none

   integer :: nbod,nbsub,nb_el,imod,ifemT,imodT, NDFTB, Nmode


   UT_md (:)   = 0.d0
   UT1_md(:)   = 0.d0
   UT2_md(:)   = 0.d0

   do nbod     = 1, NBODBT_el
   do nbsub    = 1, body(nbod)%NBODSUB_el
      nb_el    = body(nbod)%NBODTGLB_el(nbsub)

      NDFTB    = subbody(nb_el)%NDFTB_el
      Nmode    = subbody_md(nb_el)%Nmode_TR_md

      do imod  = 1, Nmode
         ifemT = ACCU%NDFTACC_el(nb_el-1)
         imodT =      NDFTACC_md(nb_el-1) + imod

         UT_md (ifemT+1:ifemT+NDFTB) = UT_md (ifemT+1:ifemT+NDFTB) + subbody_md(nb_el)%FTT_TR_md (1:NDFTB,imod) * UT_el (imodT) !UT_el :Acoef
         UT1_md(ifemT+1:ifemT+NDFTB) = UT1_md(ifemT+1:ifemT+NDFTB) + subbody_md(nb_el)%FTT_TR_md (1:NDFTB,imod) * UT1_el(imodT) !UT1_el:dAcoef /dt
         UT2_md(ifemT+1:ifemT+NDFTB) = UT2_md(ifemT+1:ifemT+NDFTB) + subbody_md(nb_el)%FTT_TR_md (1:NDFTB,imod) * UT2_el(imodT) !UT2_el:d2Acoef/dt2
      enddo !imod

   enddo !nbsub
   enddo !nbod

!--------- Q's
         ifemT = ACCU%NDFTACC_el(NBODT_el)
         imodT =      NDFTACC_md(NBODT_el)
         UT_md (ifemT+1:ifemT+NQS) = UT_el (imodT+1:imodT+NQS)
         UT1_md(ifemT+1:ifemT+NQS) = UT1_el(imodT+1:imodT+NQS)
         UT2_md(ifemT+1:ifemT+NQS) = UT2_el(imodT+1:imodT+NQS)


 END Subroutine modal2fem 
!----------------------------------------------------------------------
 Subroutine Get_Acoef (nb_el, Nmode, Acoef0, Acoef1, Acoef2)
!----------------------------------------------------------------------

 use Cbeam
 use modal_mod

   implicit none

   integer :: nb_el
   integer :: Nmode
   integer :: i
   real(8) :: Acoef0 (NMODEM_md)
   real(8) :: Acoef1 (NMODEM_md)
   real(8) :: Acoef2 (NMODEM_md)


   i               = NDFTACC_md(nb_el-1)
   Acoef0(1:Nmode) = UT_el (i+1:i+Nmode)
   Acoef1(1:Nmode) = UT1_el(i+1:i+Nmode)
   Acoef2(1:Nmode) = UT2_el(i+1:i+Nmode)


 END Subroutine Get_Acoef
!-------------------------------------------------------------------------
!   Subroutine :MATRLOC_md
! 
!-------------------------------------------------------------------------
 Subroutine MATRLOC_md ( AMLOC, ACLOC, AKLOC, AFLOC, AMLOCNQ, ACLOCNQ, AKLOCNQ,&
                         nb_el, nn_el, Nmode                                     )
!-------------------------------------------------------------------------

 use Cbeam
 use modal_mod

   implicit none

   real(8) :: ATDAx2  (3,3)
   real(8) :: ATDDA   (3,3)
   real(8) :: ATDDR   (  3)
   real(8) :: AT      (3,3)

   real(8) :: ATDAx20 (3,3)
   real(8) :: ATDAx21 (3,3)
   real(8) :: ATDDA0  (3,3)
   real(8) :: ATDDA1  (3,3)
   real(8) :: ATDDA2  (3,3)
   real(8) :: AT0     (3,3)
   real(8) :: ATDDR0  (  3)
   real(8) :: ATDDR1  (  3)
   real(8) :: ATDDR2  (  3)

   real(8) :: AKLOC   (NMODEM_md,NMODEM_md)
   real(8) :: ACLOC   (NMODEM_md,NMODEM_md)
   real(8) :: AMLOC   (NMODEM_md,NMODEM_md)
   real(8) :: AFLOC   (NMODEM_md          )
   real(8) :: AMLOCNQ (NMODEM_md,NQS      )
   real(8) :: ACLOCNQ (NMODEM_md,NQS      )
   real(8) :: AKLOCNQ (NMODEM_md,NQS      )
   real(8) :: Acoef0  (NMODEM_md          )
   real(8) :: Acoef1  (NMODEM_md          )
   real(8) :: Acoef2  (NMODEM_md          )
   integer :: nb_el, nn_el
   integer :: Nmode
   integer :: n, iq, mm,nn


!--- Initialize local matrices
   AMLOC   = 0.d0;
   ACLOC   = 0.d0;
   AKLOC   = 0.d0;
   AFLOC   = 0.d0;
   AMLOCNQ = 0.d0;
   ACLOCNQ = 0.d0;
   AKLOCNQ = 0.d0;

!--- Pass Transformation matrices to local arrays
   ATDAx2(1:3,1:3) = transf_mat(nn_el)%ATDA_el  (1:3,1:3) * 2.d0
   ATDDA (1:3,1:3) = transf_mat(nn_el)%ATDDA_el (1:3,1:3)
   AT    (1:3,1:3) = transf_mat(nn_el)%AT_el    (1:3,1:3)
   ATDDR (1:3    ) = transf_mat(nn_el)%ATDDR_el (1:3    )


   call Get_Acoef (nb_el, Nmode, Acoef0, Acoef1, Acoef2)


      AKLOC(1:Nmode,1:Nmode) =                            subbody_md(nb_el)%FT_Int_K_F_bod        (1:Nmode,      1:Nmode)
   do mm = 1, 3
      AMLOC(1:Nmode,1:Nmode) = AMLOC(1:Nmode,1:Nmode)  +  subbody_md(nb_el)%FT_Int_NTII33SN_F_bod (1:Nmode,mm,mm,1:Nmode)                   ![I]: diagonal terms only
   do nn = 1, 3
      ACLOC(1:Nmode,1:Nmode) = ACLOC(1:Nmode,1:Nmode)  +  subbody_md(nb_el)%FT_Int_NTII33SN_F_bod (1:Nmode,mm,nn,1:Nmode) * ATDAx2(mm,nn)
      AKLOC(1:Nmode,1:Nmode) = AKLOC(1:Nmode,1:Nmode)  +  subbody_md(nb_el)%FT_Int_NTII33SN_F_bod (1:Nmode,mm,nn,1:Nmode) * ATDDA (mm,nn)
      AFLOC(1:Nmode        ) = AFLOC(1:Nmode        )  -  subbody_md(nb_el)%FT_Int_NTII33r0_bod   (1:Nmode,mm,nn        ) * ATDDA (mm,nn)
   enddo !nn
      AFLOC(1:Nmode        ) = AFLOC(1:Nmode        )  +  subbody_md(nb_el)%FT_Int_NTII3_bod      (1:Nmode,mm           ) * AT    (mm,3 ) * (-GRAV) &
                                                       -  subbody_md(nb_el)%FT_Int_NTII3_bod      (1:Nmode,mm           ) * ATDDR (mm   )
   enddo !mm
      AFLOC(1:Nmode        ) = AFLOC(1:Nmode        )  -  matmul (  AKLOC (1:Nmode,1:Nmode), Acoef0 ( 1:Nmode )) &
                                                       -  matmul (  ACLOC (1:Nmode,1:Nmode), Acoef1 ( 1:Nmode )) &
                                                       -  matmul (  AMLOC (1:Nmode,1:Nmode), Acoef2 ( 1:Nmode ))

!--- [M,C,K] matrices qs terms
   do iq = 1, QSnode(nn_el)%Tot_el(2)
      n  =    QSnode(nn_el)%IORD_el(iq)

      ATDAx20(1:3,1:3) = transf_mat(nn_el)%ATDA0_el  (iq,1:3,1:3) * 2.d0
      ATDAx21(1:3,1:3) = transf_mat(nn_el)%ATDA1_el  (iq,1:3,1:3) * 2.d0
      ATDDA0 (1:3,1:3) = transf_mat(nn_el)%ATDDA0_el (iq,1:3,1:3)
      ATDDA1 (1:3,1:3) = transf_mat(nn_el)%ATDDA1_el (iq,1:3,1:3)
      ATDDA2 (1:3,1:3) = transf_mat(nn_el)%ATDDA2_el (iq,1:3,1:3)
      AT0    (1:3,1:3) = transf_mat(nn_el)%AT0_el    (iq,1:3,1:3)
      ATDDR0 (1:3    ) = transf_mat(nn_el)%ATDDR0_el (iq,1:3    )
      ATDDR1 (1:3    ) = transf_mat(nn_el)%ATDDR1_el (iq,1:3    )
      ATDDR2 (1:3    ) = transf_mat(nn_el)%ATDDR2_el (iq,1:3    )


      do mm = 1, 3
         AMLOCNQ(1:Nmode,n) = AMLOCNQ(1:Nmode,n)  +          subbody_md(nb_el)%FT_Int_NTII3_bod      (1:Nmode,mm                               ) * ATDDR2  (mm   )
         ACLOCNQ(1:Nmode,n) = ACLOCNQ(1:Nmode,n)  +          subbody_md(nb_el)%FT_Int_NTII3_bod      (1:Nmode,mm                               ) * ATDDR1  (mm   )
         AKLOCNQ(1:Nmode,n) = AKLOCNQ(1:Nmode,n)  +          subbody_md(nb_el)%FT_Int_NTII3_bod      (1:Nmode,mm                               ) * ATDDR0  (mm   ) &
                                                  -          subbody_md(nb_el)%FT_Int_NTII3_bod      (1:Nmode,mm                               ) * AT0     (mm,3 ) * (-GRAV)
      do nn = 1, 3
         AMLOCNQ(1:Nmode,n) = AMLOCNQ(1:Nmode,n)  +          subbody_md(nb_el)%FT_Int_NTII33r0_bod   (1:Nmode,mm,nn        )                     * ATDDA2  (mm,nn) &
                                                  +  matmul( subbody_md(nb_el)%FT_Int_NTII33SN_F_bod (1:Nmode,mm,nn,1:Nmode), Acoef0 (1:Nmode )) * ATDDA2  (mm,nn)

         ACLOCNQ(1:Nmode,n) = ACLOCNQ(1:Nmode,n)  +          subbody_md(nb_el)%FT_Int_NTII33r0_bod   (1:Nmode,mm,nn        )                     * ATDDA1  (mm,nn) &
                                                  +  matmul( subbody_md(nb_el)%FT_Int_NTII33SN_F_bod (1:Nmode,mm,nn,1:Nmode), Acoef0 (1:Nmode )) * ATDDA1  (mm,nn) &
                                                  +  matmul( subbody_md(nb_el)%FT_Int_NTII33SN_F_bod (1:Nmode,mm,nn,1:Nmode), Acoef1 (1:Nmode )) * ATDAx21 (mm,nn)

         AKLOCNQ(1:Nmode,n) = AKLOCNQ(1:Nmode,n)  +          subbody_md(nb_el)%FT_Int_NTII33r0_bod   (1:Nmode,mm,nn        )                     * ATDDA0  (mm,nn) &
                                                  +  matmul( subbody_md(nb_el)%FT_Int_NTII33SN_F_bod (1:Nmode,mm,nn,1:Nmode), Acoef0 (1:Nmode )) * ATDDA0  (mm,nn) &
                                                  +  matmul( subbody_md(nb_el)%FT_Int_NTII33SN_F_bod (1:Nmode,mm,nn,1:Nmode), Acoef1 (1:Nmode )) * ATDAx20 (mm,nn)
      enddo !nn
      enddo !mm
   enddo !iq


 END Subroutine MATRLOC_md
!-----------------------------------------------------------------------------------------------------
!
! -- Subroutine :MATRIX_prepro_md
!
!   AM_el , AK_el, AC_el, AQ_el  dependent only to structure and not of time
!   mass, stifness damping and forcing matrices for the complete structure
!
!----------------------------------------------------------------------
 Subroutine MATRIX_prepro_md
!----------------------------------------------------------------------

 use Cbeam
 use mod_types_el
 use modal_mod

   implicit none
   
   integer :: nbod, nbsub, nb_el, nn_el, nel, nn
   integer :: i,j, m ,n
   integer :: NDFTB, Nmode
   integer :: NDFPE
   integer :: nbod_cm, nbsub_cm, nel_cm, ncm

!--- local matrices
   real(8) :: NtKN     (NDFPEM_el,NDFPEM_el)
   real(8) :: NtII33SN (NDFPEM_el,NDFPEM_el)
   real(8) :: NtII3    (NDFPEM_el          )
   real(8) :: NtII33r0 (NDFPEM_el          )
!--- local matrices all mn terms
   real(8) :: Int_K         (NDFPEM_el,    NDFPEM_el)       !local Int {(N't*K1*N')+(N't*K2*N)+(Nt*K3*N')+(Nt*K4*N)}   structural terms
   real(8) :: Int_NTII33SN  (NDFPEM_el,3,3,NDFPEM_el)       !local Int {Nt*II*@@*S*N},  @@(t)mn [3x3]
   real(8) :: Int_NTII3     (NDFPEM_el,3            )       !local Int {Nt*II*@     },   @(t)m  [3x1]
   real(8) :: Int_NTII33r0  (NDFPEM_el,3,3          )       !local Int {Nt*II*@@r0  },  @@(t)mn [3x3]
!--- body matrices
   real(8),allocatable :: Int_K_bod             (:,    :)   !(NDFTB,    NDFTB)
   real(8),allocatable :: Int_NTII33SN_bod      (:,:,:,:)   !(NDFTB,3,3,NDFTB)
   real(8),allocatable :: Int_NTII3_bod         (:,:    )   !(NDFTB,3,       )
   real(8),allocatable :: Int_NTII33r0_bod      (:,:,:  )   !(NDFTB,3,3      )


   do nbod  = 1, NBODBT_el
   do nbsub = 1,body(nbod)%NBODSUB_el
      nb_el = body(nbod)%NBODTGLB_el(nbsub)
      nn_el = body(nbod)%NNODTGLB_el(nbsub)
      NDFPE = subbody   (nb_el)%NDFPE_el
      NDFTB = subbody   (nb_el)%NDFTB_el
      Nmode = subbody_md(nb_el)%Nmode_TR_md

!------ Allocate body matrices
      Allocate (Int_K_bod          (NDFTB,    NDFTB))
      Allocate (Int_NTII33SN_bod   (NDFTB,3,3,NDFTB))
      Allocate (Int_NTII3_bod      (NDFTB,3        ))
      Allocate (Int_NTII33r0_bod   (NDFTB,3,3      ))

!------ initialize body matrices
      Int_K_bod        = 0.d0;
      Int_NTII33SN_bod = 0.d0;
      Int_NTII3_bod    = 0.d0;
      Int_NTII33r0_bod = 0.d0;

      do nel = subbody(nb_el)%NTEB_el,1,-1
!-------------------------------------------
! A. gel the local matrices [for all nels]
!-------------------------------------------
         do m = 1, 3
         do n = 1, 3
              call MATRLOC_prepro_md ( NtKN , NtII33SN, NtII3, NtII33r0, &
                                      nb_el, nel     , m    , n           )

           do ncm      = 1, NCOMA_el
              nbod_cm  = conc_mass(ncm)%NBCOMA_el
              nbsub_cm = conc_mass(ncm)%NSBCOMA_el
              nel_cm   = conc_mass(ncm)%NELCOMA_el

              if (nbod_cm==nbod.and.nbsub_cm==nbsub.and.nel_cm==nel)  &
              call CON_MATRLOC_prepro_md ( NtII33SN, NtII3, NtII33r0, &
                                      ncm, nb_el, nel     , m    , n   )
           enddo !ncm

!             call aer_MATRLOC_prepro_md ( NtIIaeNae, nb_el, nel    )

           Int_K       (1:NDFPEM_el,     1:NDFPEM_el) = NtKN     (1:NDFPEM_el,1:NDFPEM_el); !local Int {(N't*K1*N')+(N't*K2*N)+(Nt*K3*N')+(Nt*K4*N)}   structural terms
           Int_NTII33SN(1:NDFPEM_el,m,n ,1:NDFPEM_el) = NtII33SN (1:NDFPEM_el,1:NDFPEM_el); !local Int {Nt*II*@@*S*N},  @@(t)mn [3x3]
           Int_NTII3   (1:NDFPEM_el,m               ) = NtII3    (1:NDFPEM_el            ); !local Int {Nt*II*@     },   @(t)m  [3x1]
           Int_NTII33r0(1:NDFPEM_el,m,n             ) = NtII33r0 (1:NDFPEM_el            ); !local Int {Nt*II*@@r0  },  @@(t)mn [3x3]
         enddo !n
         enddo !m
!-------------------------------------------------------------------
! B. assembly the local matrices to  the global matrix of the body
!-------------------------------------------------------------------
         nn       = subbody(nb_el)%NODMB       (nel,1)!nnod)
         i        = subbody(nb_el)%NDFPNBACC_el(nn-1 )
         j        = i
         do m = 1, 3
         do n = 1, 3  
            Int_NTII33SN_bod (i+1:i+NDFPE,m,n,j+1:j+NDFPE) = Int_NTII33SN_bod (i+1:i+NDFPE,m,n,j+1:j+NDFPE)  +  Int_NTII33SN  (1:NDFPE,m,n,1:NDFPE)
            Int_NTII33r0_bod (i+1:i+NDFPE,m,n            ) = Int_NTII33r0_bod (i+1:i+NDFPE,m,n            )  +  Int_NTII33r0  (1:NDFPE,m,n        )
         enddo !n
            Int_NTII3_bod    (i+1:i+NDFPE,m              ) = Int_NTII3_bod    (i+1:i+NDFPE,m              )  +  Int_NTII3     (1:NDFPE,m          )
         enddo !m
            Int_K_bod        (i+1:i+NDFPE,    j+1:j+NDFPE) = Int_K_bod        (i+1:i+NDFPE,    j+1:j+NDFPE)  +  Int_K         (1:NDFPE,    1:NDFPE)
      enddo !nel
!-------------------------------------------------------------------------------------------------
! C. reduction of the global matrix, based on the trunction Fi and Fi^t matrices from modal anal
!-------------------------------------------------------------------------------------------------
         do m = 1, 3
         do n = 1, 3                    !(Nmode,3,3,Nmode)
            subbody_md(nb_el)%FT_Int_NTII33SN_F_bod (1:Nmode,m,n,1:Nmode) = matmul( matmul( subbody_md(nb_el)%FTT_TR_T_md (1:Nmode,    1:NDFTB), Int_NTII33SN_bod (1:NDFTB,m,n,1:NDFTB) ), subbody_md(nb_el)%FTT_TR_md (1:NDFTB,1:Nmode) )
            subbody_md(nb_el)%FT_Int_NTII33r0_bod   (1:Nmode,m,n        ) =         matmul( subbody_md(nb_el)%FTT_TR_T_md (1:Nmode,    1:NDFTB), Int_NTII33r0_bod (1:NDFTB,m,n        ) )
         enddo !n
            subbody_md(nb_el)%FT_Int_NTII3_bod      (1:Nmode,m          ) =         matmul( subbody_md(nb_el)%FTT_TR_T_md (1:Nmode,    1:NDFTB), Int_NTII3_bod    (1:NDFTB,m          ) )
         enddo !m
            subbody_md(nb_el)%FT_Int_K_F_bod        (1:Nmode,    1:Nmode) = matmul( matmul( subbody_md(nb_el)%FTT_TR_T_md (1:Nmode,    1:NDFTB), Int_K_bod        (1:NDFTB,    1:NDFTB) ), subbody_md(nb_el)%FTT_TR_md (1:NDFTB,1:Nmode) )

      Deallocate (Int_K_bod          )
      Deallocate (Int_NTII33SN_bod   )
      Deallocate (Int_NTII3_bod      )
      Deallocate (Int_NTII33r0_bod   )

   enddo !nbsub
   enddo !nbod


 END Subroutine MATRIX_prepro_md
!-------------------------------------------------------------------------
!   Subroutine :MATRLOC_preprop_md
! 
!   Perform the spartial intergration along each subbody and store the
!   following integrals, splitted based on the elements of A or R into
!   3x3=9 or 3x1=3 sub_matrices:
!
!   local Int {Nt*II*@@*S*N},  @@(t)mn [3x3]
!   local Int {Nt*II*@     },   @(t)m  [3x1]
!   local Int {Nt*II*@@r0  },  @@(t)mn [3x3]
!   local Int {(N't*K1*N')+(N't*K2*N)+(Nt*K3*N')+(Nt*K4*N)} [pure structural terms]
!-------------------------------------------------------------------------
 Subroutine MATRLOC_prepro_md ( NtKN , NtII33SN, NtII3, NtII33r0, &
                                nb_el, nel     , i    , j           )
!-------------------------------------------------------------------------

 use Cbeam

   implicit none

   integer :: i, j
   real(8) :: AZERO   (3,3)
   real(8) :: AUNIT   (3,3)
   real(8) :: RZERO   (3  )
   real(8) :: RUNIT   (3  )
   real(8) :: NtKN0    (NDFPEM_el,NDFPEM_el)
   real(8) :: NtII33SN0(NDFPEM_el,NDFPEM_el)
   real(8) :: NtII30   (NDFPEM_el          )
   real(8) :: NtII33r00(NDFPEM_el          )
   real(8) :: NtKN     (NDFPEM_el,NDFPEM_el)
   real(8) :: NtII33SN (NDFPEM_el,NDFPEM_el)
   real(8) :: NtII3    (NDFPEM_el          )
   real(8) :: NtII33r0 (NDFPEM_el          )
   real(8), dimension(NEQPEM_el,NEQPEM_el) :: AM, AK1, AK2, AK3, AK4
   real(8) :: AQ       (NEQPEM_el          )
   real(8) :: SHAPE    (NEQPEM_el,NDFPEM_el)
   real(8) :: DSHAPE   (NEQPEM_el,NDFPEM_el)
   real(8) :: SHAPET   (NDFPEM_el,NEQPEM_el)
   real(8) :: DSHAPET  (NDFPEM_el,NEQPEM_el)
   real(8) :: SHAPETAM (NDFPEM_el,NEQPEM_el)
   integer :: nb_el, nel, IB, NEQPE, NDFPE, k, n1, n2
   integer :: ii, jj, IND(6)
   real(8) :: ALLOC, PA01 , PA02  , allocwgauss
   real(8) :: DENS , AMOMX, AMOMZ , RIXX  , RIZZ  , RIXZ  , POLI  , HTA   , HTAL
   real(8) :: PHPX , PHPZ , HTA1  , HTA0
   real(8) :: Kf(21)


!--- Initialize local matrices
   NtKN      = 0.d0;
   NtII33SN  = 0.d0;
   NtII3     = 0.d0;
   NtII33r0  = 0.d0;

!--- Set A, R based on i,j [inputs]
   AZERO (1:3,1:3) = 0.d0
   RZERO (1:3    ) = 0.d0
   AUNIT (1:3,1:3) = 0.d0
   RUNIT (1:3    ) = 0.d0
   RUNIT (i      ) = 1.d0
   AUNIT (i  ,j  ) = 1.d0


   IB      = subbody(nb_el)%IBTYP_el
   NEQPE   = subbody(nb_el)%NEQPE_el
   NDFPE   = subbody(nb_el)%NDFPE_el
   n1      = subbody(nb_el)%INODNEL_el(nel,1)
   n2      = subbody(nb_el)%INODNEL_el(nel,2)
   HTA0    = subbody(nb_el)%HTA_el    (nel)
   ALLOC   = subbody(nb_el)%ALENG_el  (nel)
   PHPX    = subbody(nb_el)%PHIX_el   (nel)
   PHPZ    = subbody(nb_el)%PHIZ_el   (nel)
   IND(:)  = IMDOF_el(:)


   do k        = 1, IORDER
      HTAL     = ( XGAUSS(k) + 1.d0 )/2.d0   ! 0 < HTAL   < 1 HTAL=HTA/L
      HTA      = HTAL*ALLOC                  ! 0 < HTA    < L
      HTA1     = HTA0 + HTA
      PA01     = 1.d0 - HTAL
      PA02     = HTAL

      call SHAPEFUNC15 ( HTAL, ALLOC, PHPX, PHPZ, SHAPE, DSHAPE, NDFPE, NEQPE )

      SHAPET   = transpose (SHAPE )
      DSHAPET  = transpose (DSHAPE)
!------ mass
      DENS     = PA01 * beam_timo(n1)%DENS_el   +  PA02 * beam_timo(n2)%DENS_el
      AMOMX    = PA01 * beam_timo(n1)%AMOMX_el  +  PA02 * beam_timo(n2)%AMOMX_el
      AMOMZ    = PA01 * beam_timo(n1)%AMOMZ_el  +  PA02 * beam_timo(n2)%AMOMZ_el
      RIXX     = PA01 * beam_timo(n1)%RIXX_el   +  PA02 * beam_timo(n2)%RIXX_el
      RIZZ     = PA01 * beam_timo(n1)%RIZZ_el   +  PA02 * beam_timo(n2)%RIZZ_el
      RIXZ     = PA01 * beam_timo(n1)%RIXZ_el   +  PA02 * beam_timo(n2)%RIXZ_el
      POLI     = PA01 * beam_timo(n1)%POLI_el   +  PA02 * beam_timo(n2)%POLI_el
!------ stiffness
      Kf(1:21) = PA01 * beam_timo(n1)%K(1:21)   +  PA02 * beam_timo(n2)%K(1:21)


!------ local Int {Nt*II*@@*S*N},  @@(t)mn [3x3]
      call IIxAxS (AM, AUNIT, DENS, AMOMX, AMOMZ, RIXX, RIXZ, RIZZ, 1.d0)
      SHAPETAM  = matmul( SHAPET  , AM    )
      NtII33SN0 = matmul( SHAPETAM, SHAPE )

!------ local Int {Nt*II*@     },   @(t)m  [3x1]
      AQ(1:6)   = 0.d0
      call IIxAtDDR_AtDDAr0 ( AQ, AZERO, RUNIT, DENS, AMOMX, AMOMZ, RIXX, RIXZ, RIZZ, HTA1, 1.d0 )
      NtII30    =  matmul( SHAPET, AQ )

!------ local Int {Nt*II*@@r0  },  @@(t)mn [3x3]
      AQ(1:6)   = 0.d0
      call IIxAtDDR_AtDDAr0 ( AQ, AUNIT, RZERO, DENS, AMOMX, AMOMZ, RIXX, RIXZ, RIZZ, HTA1, 1.d0 )
      NtII33r00 =  matmul( SHAPET, AQ )


!------ local Int {(N't*K1*N')+(N't*K2*N)+(Nt*K3*N')+(Nt*K4*N)}   [pure structural terms]
!K1
      AK1(1:6   ,1:6     ) = 0.d0;
      AK1(IND(1),IND(1:6)) = Kf( 1: 6)
      AK1(IND(2),IND(2:6)) = Kf( 7:11)
      AK1(IND(3),IND(3:6)) = Kf(12:15)
      AK1(IND(4),IND(4:6)) = Kf(16:18)
      AK1(IND(5),IND(5:6)) = Kf(19:20)
      AK1(IND(6),IND(  6)) = Kf(   21)
      do ii=2,6
      do jj=1,ii-1
      AK1(IND(ii),IND(jj)) = AK1(IND(jj),IND(ii))
      enddo
      enddo
!K2
      AK2(1:6   ,1:6   ) = 0.d0;
      AK2(1:6   ,IND(4)) =-AK1(1:6   ,IND(3))
      AK2(1:6   ,IND(6)) = AK1(1:6   ,IND(1))
!K3 [-]
      AK3(1:6   ,1:6   ) = 0.d0;
      AK3(IND(4),1:6   ) = AK1(IND(3),1:6   )
      AK3(IND(6),1:6   ) =-AK1(IND(1),1:6   )
!K4 [-]
      AK4                = 0.d0;
      AK4(1:6   ,IND(4)) =-AK3(1:6   ,IND(3))
      AK4(1:6   ,IND(6)) = AK3(1:6   ,IND(1))

      SHAPETAM           =          matmul( DSHAPET  ,  AK1   )
      NtKN0              =          matmul(  SHAPETAM, DSHAPE )
      SHAPETAM           =          matmul( DSHAPET  ,  AK2   )
      NtKN0              = NtKN0  + matmul(  SHAPETAM,  SHAPE )
      SHAPETAM           =          matmul(  SHAPET  ,  AK3   )
      NtKN0              = NtKN0  - matmul(  SHAPETAM, DSHAPE )
      SHAPETAM           =          matmul(  SHAPET  ,  AK4   )
      NtKN0              = NtKN0  - matmul(  SHAPETAM,  SHAPE )

!------ Replace
      allocwgauss = 0.5d0*ALLOC*WGAUSS(k)
      NtKN        = NtKN     + allocwgauss*NtKN0      ! Nt'K1N+Nt ... K2 + K3  + K4   
      NtII33SN    = NtII33SN + allocwgauss*NtII33SN0  ! Nt. II. @@ .S .N
      NtII3       = NtII3    + allocwgauss*NtII30     ! Nt. II. @
      NtII33r0    = NtII33r0 + allocwgauss*NtII33r00  ! Nt. II. @@ .ro
   enddo !k


 END Subroutine MATRLOC_prepro_md
!-------------------------------------------------------------------------
!   Subroutine :CON_MATRLOC_prepro_md
! 
!  All Local Matrices [MASS, DAMPING, STIFFNESS, FORCES]
!  with concentrated masses calculated 
!  with respect to the body FEM local system.
!
!-------------------------------------------------------------------------
 Subroutine CON_MATRLOC_prepro_md ( NtII33SN, NtII3, NtII33r0,      &
                                    ncm     , nb_el, nel     , i, j   )
!-------------------------------------------------------------------------

 use Cbeam
 use mod_types_el

   implicit none

   real(8) :: AZERO   (3,3)
   real(8) :: AUNIT   (3,3)
   real(8) :: RZERO   (3  )
   real(8) :: RUNIT   (3  )

   real(8) :: NtII33SN (NDFPEM_el,NDFPEM_el)
   real(8) :: NtII3    (NDFPEM_el          )
   real(8) :: NtII33r0 (NDFPEM_el          )

   real(8) :: NtII33SN0(NDFPEM_el,NDFPEM_el)
   real(8) :: NtII30   (NDFPEM_el          )
   real(8) :: NtII33r00(NDFPEM_el          )

   real(8) :: AM       (NEQPEM_el,NEQPEM_el)
   real(8) :: AQ       (NEQPEM_el          )

   real(8) :: SHAPE    (NEQPEM_el,NDFPEM_el)
   real(8) :: DSHAPE   (NEQPEM_el,NDFPEM_el)
   real(8) :: SHAPET   (NDFPEM_el,NEQPEM_el)
   real(8) :: SHAPETAM (NDFPEM_el,NEQPEM_el)

   integer :: i, j

   real(8) :: ALLOC, y0  , HTAL, Hta, PHPX, PHPZ, RAT
   real(8) :: Mass  , Sxoff , Syoff , Szoff
   real(8) :: IxxTot, IyyTot, IzzTot, IxyTot   , IxzTot, Iyztot
   real(8) :: Sy0   , Ixy0  , Iyz0  , Iyy0
   real(8) :: Xoff  , Yoff  , Zoff

   integer :: nb_el, nel, NEQPE, NDFPE, ncm, nbsha_el


!--- Set A, R based on i,j [inputs]
   AZERO (1:3,1:3) = 0.d0
   RZERO (1:3    ) = 0.d0
   AUNIT (1:3,1:3) = 0.d0
   RUNIT (1:3    ) = 0.d0
   RUNIT (i      ) = 1.d0
   AUNIT (i  ,j  ) = 1.d0

   NEQPE   = subbody(nb_el)%NEQPE_el
   NDFPE   = subbody(nb_el)%NDFPE_el
   y0      = subbody(nb_el)%HTA_el   (nel)
   ALLOC   = subbody(nb_el)%ALENG_el (nel)
   PHPX    = subbody(nb_el)%PHIX_el  (nel)
   PHPZ    = subbody(nb_el)%PHIZ_el  (nel)
   Hta     = conc_mass(ncm)%Hta
   HTAL    = Hta/ALLOC


   call SHAPEFUNC15 ( HTAL, ALLOC, PHPX, PHPZ, SHAPE, DSHAPE, NDFPE, NEQPE )

   SHAPET  = transpose (SHAPE )


!--- Vars
   Mass    = conc_mass(ncm)%Mass
!  Hta     = conc_mass(ncm)%Hta
   Xoff    = conc_mass(ncm)%Xoff
   Yoff    = conc_mass(ncm)%Yoff
   Zoff    = conc_mass(ncm)%Zoff
   Sxoff   = Mass*Xoff
   Syoff   = Mass*Yoff
   Szoff   = Mass*Zoff
   IxxTot  = conc_mass(ncm)%Ixx + Mass*Xoff**2
   IyyTot  = conc_mass(ncm)%Iyy + Mass*Yoff**2
   IzzTot  = conc_mass(ncm)%Izz + Mass*Zoff**2
   IxyTot  = conc_mass(ncm)%Ixy + Mass*Xoff*Yoff
   IxzTot  = conc_mass(ncm)%Ixz + Mass*Xoff*Zoff
   IyzTot  = conc_mass(ncm)%Iyz + Mass*Yoff*Zoff

   Sy0     = Mass     *(y0+Hta)
   Ixy0    = Mass*Xoff*(y0+hta)
   Iyy0    = Mass*Yoff*(y0+hta)
   Iyz0    = Mass*Zoff*(y0+hta)

   nbsha_el= body(NBLADE_el+1)%NBODTGLB_el (1)

   RAT     = 1.d0
   if ( (nb_el == nbsha_el).and.(y0 < 0.01d0) ) &
   RAT     = RAT_GEAR**2

!--- local Int {Nt*II*@@*S*N},  @@(t)mn [3x3]
   call Con_IIxAxS ( AM    , AUNIT , RAT          ,&
                     Mass  , Sxoff , Syoff , Szoff,&
                     IxxTot, IyyTot, IzzTot       ,&
                     IxyTot, IxzTot, Iyztot       ,&
                     IMDOF_el                        )
   SHAPETAM  = matmul( SHAPET  , AM    )
   NtII33SN0 = matmul( SHAPETAM, SHAPE )

!--- local Int {Nt*II*@     },   @(t)m  [3x1]
   AQ(1:6)   = 0.d0
   call Con_IIxAtDDR_AtDDAr0 ( AQ    , AZERO , RUNIT , RAT   ,&
                               Mass  , Sxoff , Syoff , Szoff ,&
                               IxxTot, IyyTot, IzzTot        ,&
                               IxyTot, IxzTot, Iyztot        ,&
                               Sy0   , Ixy0  , Iyz0  , Iyy0  ,&
                               IMDOF_el                         )
   NtII30    =  matmul( SHAPET, AQ )

!--- local Int {Nt*II*@@r0  },  @@(t)mn [3x3]
   AQ(1:6)   = 0.d0
   call Con_IIxAtDDR_AtDDAr0 ( AQ    , AUNIT , RZERO , RAT   ,&
                               Mass  , Sxoff , Syoff , Szoff ,&
                               IxxTot, IyyTot, IzzTot        ,&
                               IxyTot, IxzTot, Iyztot        ,&
                               Sy0   , Ixy0  , Iyz0  , Iyy0  ,&
                               IMDOF_el                         )
   NtII33r00 =  matmul( SHAPET, AQ )

!--- Replace
   NtII33SN    = NtII33SN + NtII33SN0  ! Nt. II. @@ .S .N
   NtII3       = NtII3    + NtII30     ! Nt. II. @
   NtII33r0    = NtII33r0 + NtII33r00  ! Nt. II. @@ .ro


 END Subroutine CON_MATRLOC_prepro_md

#include "gast2modal_sb.f90"
#include "loads_integr.f90"
