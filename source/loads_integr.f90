!----------------------------------------------------------------------
!
!  Subrouine : BOD2BOD_LOADS_int  ----------------
!
!  Loads Communication between WT's bodies (1st Subbody with last one of the previous body)
!
!----------------------------------------------------------------------
 Subroutine BOD2BOD_LOADS_int ( nbod1, nbod2 )
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   real(8)             :: Aba      (3,3    )
   real(8)             :: Aba0     (3,3,NQS)
   real(8)             :: Rba      (3      )
   real(8)             :: Rba0     (3  ,NQS)
   real(8)             :: RGba     (3      )
   real(8),allocatable :: AMB2B    (:  ,:  )  !(6,NDFPEM_el)
   real(8),allocatable :: ACB2B    (:  ,:  )  !(6,NDFPEM_el)
   real(8),allocatable :: AKB2B    (:  ,:  )  !(6,NDFPEM_el)
   real(8)             :: AMB2BNQ  (6  ,NQS)
   real(8)             :: ACB2BNQ  (6  ,NQS)
   real(8)             :: AKB2BNQ  (6  ,NQS)
   real(8)             :: AFB2B    (6      )
   integer             :: nbod1, nbod2, i, j, NDFPE, iq, nq
   integer             :: nbodB, nnsubB, nbsubB, nnB_el,          nbB_el, nelB, nodB
   integer             :: nbodA, nn    ,         nnA_el, nnA2_el, nbA_el, nelA, nodA
   integer             :: IMAX , JMAX  , IQRayl, IcalcR
   real(8)             :: COEFM, COEFK
!--- integration
   integer             :: integr


   integr     = 0 !Loads calculation flag. 0:reaction, 1:integration


   if (NBODBTWT_el == NBLADE_el) return !only blades, so no loading transfer
   if (NBODBTWT_el == nbod2    ) return !only communicate blades to shaft


   do nbodB   = nbod1, nbod2
!------ Body (B)
      nnsubB  = 1
      nbsubB  = 1
      nelB    = 1
      nodB    = 1
      nnB_el  = body(nbodB)%NNODTGLB_el(nnsubB)
      nbB_el  = body(nbodB)%NBODTGLB_el(nbsubB)
      NDFPE   = subbody(nbB_el)%NDFPE_el

!------ Body (A)
      nbodA   = nbod2+1
      nnA_el  = body(nbodA)%NNODTGLB_el (body(nbodA)%NBODSUB_el   )
      nnA2_el = body(nbodA)%NNODTGLB_el (body(nbodA)%NBODSUB_el+1 )
      nbA_el  = body(nbodA)%NBODTGLB_el (body(nbodA)%NBODSUB_el   )
      nelA    = subbody(nbA_el)%NTEB_el
      nodA    = NNPE_el
      nn      = subbody(nbA_el)%NODMB(nelA,nodA)
      i       = ACCU%NDFTACC_el(nbA_el-1)+subbody(nbA_el)%NDFPNBACC_el(nn-1     )!+1:6
      j       = ACCU%NDFTACC_el(nbB_el-1)!1:NDFPE !+subbody(nbB_el)%NDFPNBACC_el(nnB-1)+mdof

!------ Aba, Aba0, Rba, Rba0
      RGba(1:3    ) = transf_mat(nnB_el)%R_el(1:3) - transf_mat(nnA2_el)%R_el(1:3)
      Aba (1:3,1:3) = matmul ( transf_mat(nnA_el)%AT_el(1:3,1:3), transf_mat(nnB_el)%A_el(1:3,1:3) )
      Rba (1:3    ) = matmul ( transf_mat(nnA_el)%AT_el(1:3,1:3), RGba(1:3)                        )

      do iq   = 1, QSnode(nnA_el)%Tot_el(2)
         nq   = QSnode(nnB_el)%IORD_el(iq)
         Aba0(1:3,1:3,nq) = matmul( transf_mat(nnA_el)%AT_el (   1:3,1:3), transf_mat(nnB_el)%A0_el(iq,1:3,1:3) ) + &
                            matmul( transf_mat(nnA_el)%AT0_el(iq,1:3,1:3), transf_mat(nnB_el)%A_el (   1:3,1:3) )

         Rba0(1:3,    nq) = matmul( (transf_mat(nnB_el)%R0_el(iq,1:3) - transf_mat(nnA2_el)%R0_el(iq,1:3) ), transf_mat(nnA_el)%AT_el (   1:3,1:3) ) + &
                            matmul(  RGba                    (   1:3)                                      , transf_mat(nnA_el)%AT0_el(iq,1:3,1:3) )
      enddo!iq 

      do iq   = QSnode(nnA_el)%Tot_el(2)+1, QSnode(nnA2_el)%Tot_el(2)
         nq   = QSnode(nnB_el)%IORD_el(iq)
         Aba0(1:3,1:3,nq) = matmul( transf_mat(nnA_el)%AT_el (   1:3,1:3), transf_mat(nnB_el)%A0_el(iq,1:3,1:3) )

         Rba0(1:3,    nq) = matmul( (transf_mat(nnB_el)%R0_el(iq,1:3) - transf_mat(nnA2_el)%R0_el(iq,1:3) ), transf_mat(nnA_el)%AT_el (   1:3,1:3) )
      enddo!iq 

      do iq   = QSnode(nnA2_el)%Tot_el(2)+1, QSnode(nnB_el)%Tot_el(2)
         nq   = QSnode(nnB_el)%IORD_el(iq)
         Aba0(1:3,1:3,nq) = matmul( transf_mat(nnA_el)%AT_el (   1:3,1:3), transf_mat(nnB_el)%A0_el(iq,1:3,1:3) )

         Rba0(1:3,    nq) = matmul( (transf_mat(nnB_el)%R0_el(iq,1:3)                                     ), transf_mat(nnA_el)%AT_el (   1:3,1:3) )
      enddo!iq 

      IcalcR  = 1   !0: no R calculations, 1: R calculations
      IQRayl  = 2   !0: all qs Rayleigh  , 1: body qs with Rayleigh, 2: no Rayleigh to qs
      IMAX    = 6
!     JMAX    = subbody(nbB_el)%NDFPE_el
      COEFM   = 0.d0
      COEFK   = 0.d0

      if (integr==0) then
!--------- Loads from Reaction, based of FEM

      JMAX    = subbody(nbB_el)%NDFPE_el
         Allocate ( AMB2B (6,JMAX) )
         Allocate ( ACB2B (6,JMAX) )
         Allocate ( AKB2B (6,JMAX) )

         call LOC_LOADS_el       ( AMB2B  , ACB2B  , AKB2B  , AFB2B,&
                                   AMB2BNQ, ACB2BNQ, AKB2BNQ       ,&
                                   nbB_el  ,nnB_el  ,nelB   , nodB    )

      else
!--------- Loads from Integration

      JMAX    = subbody   (nbB_el)%NDFTB_el
         Allocate ( AMB2B (6,JMAX) )
         Allocate ( ACB2B (6,JMAX) )
         Allocate ( AKB2B (6,JMAX) )

         call LOC_LOADS_int_el   ( AMB2B  , ACB2B  , AKB2B  , AFB2B,&
                                   AMB2BNQ, ACB2BNQ, AKB2BNQ       ,&
                                   nbodB  , nbsubB , JMAX             )


      endif
         call LOADS_TRANS_el     ( AMB2B  , ACB2B  , AKB2B  , AFB2B,&
                                   AMB2BNQ, ACB2BNQ, AKB2BNQ       ,&
                                   Aba0   , Aba    , Rba0   , Rba  ,&
                                            nnB_el , IcalcR        ,&
                                   JMAX   , JMAX                      )


         call ASSEMBLY_MATRIX_el ( AMB2B  , ACB2B  , AKB2B  , AFB2B         ,&
                                   AMB2BNQ, ACB2BNQ, AKB2BNQ, COEFM , COEFK ,&
                                            nnB_el , i      , j     , IQRayl,&
                                   IMAX   , JMAX   , 6      , JMAX            )

         DeAllocate ( AMB2B, ACB2B, AKB2B )
   enddo !nbod


 END Subroutine BOD2BOD_LOADS_int
!----------------------------------------------------------------------
 Subroutine  Prevbody_LOADS_md ( AFB2B_L, nbod, nbsub )
!----------------------------------------------------------------------

 Use Cbeam
 Use modal_mod

   implicit none

   integer :: nbod, nbsub
   real(8) :: Aba      (3,3        )
   real(8) :: Rba      (3          )
   real(8) :: RGba     (3          )
   real(8) :: AFB2B_L  (6          )
   real(8) :: AFB2B    (6          )
   integer :: nbodB, nbsubB, nnB_el, nbB_el , nbodB1, nbodB2
   integer :: nbodA, nbsubA, nnA_el, nnA2_el


   AFB2B_L  (1:6) = 0.d0;

   if (nbsub==body(nbod)%NBODSUB_el) then
!--------- body loads
      if     (nbod<=NBLADE_el  ) then !blades
         return
      elseif (nbod==NBLADE_el+1) then !shaft
         nbodB1 = 1
         nbodB2 = NBLADE_el
         nbsubB = 1
      elseif (nbod==NBLADE_el+2) then !tower
         nbodB1 = NBLADE_el+1
         nbodB2 = NBLADE_el+1
!!       if (Integr_md==1) &            !Loads calculation flag 0:reaction, 1:integration
!!       nbodB1 = 1
         nbsubB = 1
      endif
   else
!--------- subbody loads
!qq ------ give expresions for subbodies without RGab
         write(*,*)'stopped in B_LOADS_md, not ready for subbodies'
         stop
         return
         nbodB1 = nbod
         nbodB2 = nbod
         nbsubB = nbsub+1
   endif


   do nbodB   = nbodB1, nbodB2
!------ SubBody (B)
      nnB_el  = body      (nbodB )%NNODTGLB_el(nbsubB)
      nbB_el  = body      (nbodB )%NBODTGLB_el(nbsubB)

!------ SubBody (A)
      nbodA   = nbod
      nbsubA  = nbsub
      nnA2_el = body      (nbodA )%NNODTGLB_el(nbsubA+1 ) ! for R during body loads with offset
      nnA_el  = body      (nbodA )%NNODTGLB_el(nbsubA   )

!------ Aba, Rba
      RGba(1:3    ) = transf_mat(nnB_el)%R_el(1:3) - transf_mat(nnA2_el)%R_el(1:3)
      Aba (1:3,1:3) = matmul ( transf_mat(nnA_el)%AT_el(1:3,1:3), transf_mat(nnB_el)%A_el(1:3,1:3) )
      Rba (1:3    ) = matmul ( transf_mat(nnA_el)%AT_el(1:3,1:3), RGba(1:3)                        )

!------ load end load [loads transfer]
      AFB2B ( 1:6 ) = subbody_md(nbB_el)%LoadIntQ_L ( 1:6, 1 )

!------ transfer from SubBody (B) to SubBody (A) c.s.
      call LOADS_TRANS0_NEW ( AFB2B, Aba, 1, 1   , 1   )! Aba.Q
      call LOADS_CALC_M0    ( AFB2B, Rba, 1, 1   , 1   )! Rba x Aba.Q

      AFB2B_L(1:6) = AFB2B_L(1:6) + AFB2B(1:6)
   enddo !nbod


 END Subroutine  Prevbody_LOADS_md
!----------------------------------------------------------------------
 Subroutine LOC_LOADS_int_store_el
!----------------------------------------------------------------------

 Use Cbeam
 Use modal_mod

   implicit none

   real(8),allocatable :: AMB2B (:,:) !   (6,NDFTB    )
   real(8),allocatable :: ACB2B (:,:) !   (6,NDFTB    )
   real(8),allocatable :: AKB2B (:,:) !   (6,NDFTB    )
   real(8)             :: AFB2B    (6          )
   real(8)             :: AMB2BNQ  (6,NQS      )
   real(8)             :: ACB2BNQ  (6,NQS      )
   real(8)             :: AKB2BNQ  (6,NQS      )
   real(8)             :: AFB2B_L  (6          )

   real(8)             :: AMB2B0   (6,NDFPEM_el)
   real(8)             :: ACB2B0   (6,NDFPEM_el)
   real(8)             :: AKB2B0   (6,NDFPEM_el)
   real(8)             :: AFB2B0   (6          )
   real(8)             :: AMB2BNQ0 (6,NQS      )
   real(8)             :: ACB2BNQ0 (6,NQS      )
   real(8)             :: AKB2BNQ0 (6,NQS      )

   integer             :: nbod, nbsub, nb_el, nn_el, nel, nn, j
   integer             :: JMAX, ndftb
   integer             :: iq, nq
!  real(8)             :: COEFM, COEFK
   integer             :: ncm, nbod_cm, nbsub_cm, nel_cm
   real(8)             :: HTAf, Fy


   if (Integr_md/=1) return

!--- Zero local elements matrices
   do nb_el = 1, NBODT_el
      subbody_md(nb_el)%LoadIntQ_L (:, :) = 0.d0; 
      subbody_md(nb_el)%LoadIntQ   (:   ) = 0.d0; 
      subbody_md(nb_el)%LoadIntM   (:, :) = 0.d0; 
      subbody_md(nb_el)%LoadIntC   (:, :) = 0.d0; 
      subbody_md(nb_el)%LoadIntK   (:, :) = 0.d0; 
      subbody_md(nb_el)%LoadIntMq  (:, :) = 0.d0; 
      subbody_md(nb_el)%LoadIntCq  (:, :) = 0.d0; 
      subbody_md(nb_el)%LoadIntKq  (:, :) = 0.d0; 
   enddo
!--- Zero Buoyancy Forces
!  BuoyancyT_el     = 0.d0
!  BuoyancyTtap_el  = 0.d0
!  BuoyancyTside_el = 0.d0


   do nbod  = 1, NBODBT_el
   do nbsub = 1, body(nbod)%NBODSUB_el
      nb_el =    body(nbod )%NBODTGLB_el(nbsub)
      nn_el =    body(nbod )%NNODTGLB_el(nbsub)
      JMAX  = subbody(nb_el)%NDFPE_el
      ndftb = subbody(nb_el)%NDFTB_el

      Allocate ( AMB2B(6,ndftb) )
      Allocate ( ACB2B(6,ndftb) )
      Allocate ( AKB2B(6,ndftb) )

!------ previous bodies to AFB2B_L      
      call Prevbody_LOADS_md ( AFB2B_L, nbod, nbsub )

                                         nel   = subbody(nb_el)%NTEB_el + 1
      subbody_md(nb_el)%LoadIntQ_L (1:6, nel ) = AFB2B_L( 1:6 )

      AMB2B   = 0.d0;   ACB2B   = 0.d0;   AKB2B   = 0.d0;   AFB2B   = 0.d0;
      AMB2BNQ = 0.d0;   ACB2BNQ = 0.d0;   AKB2BNQ = 0.d0;

   do nel   = subbody(nb_el)%NTEB_el,1,-1
      nn    = subbody(nb_el)%NODMB       (nel,1)!nnod)
      j     = subbody(nb_el)%NDFPNBACC_el(nn-1 )

      call MATRLOC_Fint       ( AMB2B0  , ACB2B0  , AKB2B0  , AFB2B0     ,&
                                AMB2BNQ0, ACB2BNQ0, AKB2BNQ0             ,&
                                                    nb_el   , nn_el , nel   )
      call FORCLOC_Fint       ( AFB2B0  , AKB2BNQ0,                       &
                                nbod    , nbsub   , nb_el   , nn_el , nel   )

      Fy = AFB2B_L( IMDOF_el(2) ) + AFB2B0  ( IMDOF_el(2) )

      call STIFLOC_tract_Fint ( AKB2B0  , AFB2B0  , Fy      , nb_el , nel   )

!     call Buoyancy_hyd       ( AFB2B0  , AKB2BNQ0                       ,&
!                               nbod    , nbsub   , nb_el   , nn_el , nel   )

!     call Morison_hyd        ( AMB2B0  , ACB2B0  , AKB2B0  , AFB2B0     ,&
!                               AMB2BNQ0, ACB2BNQ0, AKB2BNQ0             ,&
!                               nbod    , nbsub   , nb_el   , nn_el , nel   )

!!!   call FLOODMATRLOC       ( AMB2B0  , ACB2B0  , AKB2B0  , AFB2B0     ,&
!!!                             AMB2BNQ0, ACB2BNQ0, AKB2BNQ0             ,&
!!!                             nbod    , nbsub   , nb_el   , nn_el  , nel  )

      do ncm      = 1, NCOMA_el
         nbod_cm  = conc_mass(ncm)%NBCOMA_el
         nbsub_cm = conc_mass(ncm)%NSBCOMA_el
         nel_cm   = conc_mass(ncm)%NELCOMA_el

         if (nbod_cm==nbod.and.nbsub_cm==nbsub.and.nel_cm==nel) then

         call CON_MATRLOC_Fint ( AMB2B0  , ACB2B0  , AKB2B0  , AFB2B0,&
                                 AMB2BNQ0, ACB2BNQ0, AKB2BNQ0,        &
                                 ncm     , nb_el   , nn_el   , nel      )
         endif
      enddo !ncm


!------ transfer moments from previous nels
      HTAf        = subbody(nb_el)%ALENG_el  (nel)

      AFB2B_L( IMDOF_el(4)            ) = AFB2B_L( IMDOF_el(4)            ) + AFB2B_L( IMDOF_el(3)            ) *HTAf
      AFB2B_L( IMDOF_el(6)            ) = AFB2B_L( IMDOF_el(6)            ) - AFB2B_L( IMDOF_el(1)            ) *HTAf
      AMB2B  ( IMDOF_el(4), j+1:ndftb ) = AMB2B  ( IMDOF_el(4), j+1:ndftb ) + AMB2B  ( IMDOF_el(3), j+1:ndftb ) *HTAf
      AMB2B  ( IMDOF_el(6), j+1:ndftb ) = AMB2B  ( IMDOF_el(6), j+1:ndftb ) - AMB2B  ( IMDOF_el(1), j+1:ndftb ) *HTAf
      ACB2B  ( IMDOF_el(4), j+1:ndftb ) = ACB2B  ( IMDOF_el(4), j+1:ndftb ) + ACB2B  ( IMDOF_el(3), j+1:ndftb ) *HTAf
      ACB2B  ( IMDOF_el(6), j+1:ndftb ) = ACB2B  ( IMDOF_el(6), j+1:ndftb ) - ACB2B  ( IMDOF_el(1), j+1:ndftb ) *HTAf
      AKB2B  ( IMDOF_el(4), j+1:ndftb ) = AKB2B  ( IMDOF_el(4), j+1:ndftb ) + AKB2B  ( IMDOF_el(3), j+1:ndftb ) *HTAf
      AKB2B  ( IMDOF_el(6), j+1:ndftb ) = AKB2B  ( IMDOF_el(6), j+1:ndftb ) - AKB2B  ( IMDOF_el(1), j+1:ndftb ) *HTAf
      AFB2B  ( IMDOF_el(4)            ) = AFB2B  ( IMDOF_el(4)            ) + AFB2B  ( IMDOF_el(3)            ) *HTAf
      AFB2B  ( IMDOF_el(6)            ) = AFB2B  ( IMDOF_el(6)            ) - AFB2B  ( IMDOF_el(1)            ) *HTAf

      do iq = 1, QSnode(nn_el)%Tot_el(2)
         nq = QSnode(nn_el)%IORD_el(iq)

      AMB2BNQ( IMDOF_el(4), nq        ) = AMB2BNQ( IMDOF_el(4), nq        ) + AMB2BNQ( IMDOF_el(3), nq        ) *HTAf
      AMB2BNQ( IMDOF_el(6), nq        ) = AMB2BNQ( IMDOF_el(6), nq        ) - AMB2BNQ( IMDOF_el(1), nq        ) *HTAf
      ACB2BNQ( IMDOF_el(4), nq        ) = ACB2BNQ( IMDOF_el(4), nq        ) + ACB2BNQ( IMDOF_el(3), nq        ) *HTAf
      ACB2BNQ( IMDOF_el(6), nq        ) = ACB2BNQ( IMDOF_el(6), nq        ) - ACB2BNQ( IMDOF_el(1), nq        ) *HTAf
      AKB2BNQ( IMDOF_el(4), nq        ) = AKB2BNQ( IMDOF_el(4), nq        ) + AKB2BNQ( IMDOF_el(3), nq        ) *HTAf
      AKB2BNQ( IMDOF_el(6), nq        ) = AKB2BNQ( IMDOF_el(6), nq        ) - AKB2BNQ( IMDOF_el(1), nq        ) *HTAf
      enddo !iq


!------ add current element contribution
      AFB2B_L( 1:6             ) = AFB2B_L( 1:6             ) + AFB2B0  ( 1:6         )
      AMB2B  ( 1:6, j+1:j+JMAX ) = AMB2B  ( 1:6, j+1:j+JMAX ) + AMB2B0  ( 1:6, 1:JMAX )
      ACB2B  ( 1:6, j+1:j+JMAX ) = ACB2B  ( 1:6, j+1:j+JMAX ) + ACB2B0  ( 1:6, 1:JMAX )
      AKB2B  ( 1:6, j+1:j+JMAX ) = AKB2B  ( 1:6, j+1:j+JMAX ) + AKB2B0  ( 1:6, 1:JMAX )
      AFB2B  ( 1:6             ) = AFB2B  ( 1:6             ) + AFB2B0  ( 1:6         )
     
      do iq = 1, QSnode(nn_el)%Tot_el(2)
         nq = QSnode(nn_el)%IORD_el(iq)

      AMB2BNQ( 1:6, nq         ) = AMB2BNQ( 1:6, nq         ) + AMB2BNQ0 ( 1:6, nq    )
      ACB2BNQ( 1:6, nq         ) = ACB2BNQ( 1:6, nq         ) + ACB2BNQ0 ( 1:6, nq    )
      AKB2BNQ( 1:6, nq         ) = AKB2BNQ( 1:6, nq         ) + AKB2BNQ0 ( 1:6, nq    )
      enddo !iq

!------ store load along sub-body length [Fy, writeout]
      subbody_md(nb_el)%LoadIntQ_L (1:6,    nel   ) = AFB2B_L( 1:6          )
   enddo !nel

!------ store end load [loads transfer]
      subbody_md(nb_el)%LoadIntM   ( 1:6, 1:ndftb ) = AMB2B  ( 1:6, 1:ndftb )
      subbody_md(nb_el)%LoadIntC   ( 1:6, 1:ndftb ) = ACB2B  ( 1:6, 1:ndftb )
      subbody_md(nb_el)%LoadIntK   ( 1:6, 1:ndftb ) = AKB2B  ( 1:6, 1:ndftb )
      subbody_md(nb_el)%LoadIntQ   ( 1:6          ) = AFB2B  ( 1:6          )
      do iq = 1, QSnode(nn_el)%Tot_el(2)
         nq = QSnode(nn_el)%IORD_el(iq)
      subbody_md(nb_el)%LoadIntMq  ( 1:6, nq      ) = AMB2BNQ( 1:6, nq      )
      subbody_md(nb_el)%LoadIntCq  ( 1:6, nq      ) = ACB2BNQ( 1:6, nq      )
      subbody_md(nb_el)%LoadIntKq  ( 1:6, nq      ) = AKB2BNQ( 1:6, nq      )
      enddo !iq

      Deallocate ( AMB2B, ACB2B, AKB2B )
   enddo !nbsub
   enddo !nbod


 END Subroutine LOC_LOADS_int_store_el
!----------------------------------------------------------------------
 Subroutine LOC_LOADS_int_el   ( AMB2B  , ACB2B  , AKB2B  , AFB2B,&
                                 AMB2BNQ, ACB2BNQ, AKB2BNQ       ,&
                                 nbod   , nbsub  , NDFTB )
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   integer :: NDFTB
   real(8) :: AMB2B    (6,NDFTB    )
   real(8) :: ACB2B    (6,NDFTB    )
   real(8) :: AKB2B    (6,NDFTB    )
   real(8) :: AFB2B    (6          )
   real(8) :: AMB2BNQ  (6,NQS      )
   real(8) :: ACB2BNQ  (6,NQS      )
   real(8) :: AKB2BNQ  (6,NQS      )

   real(8) :: AMB2B0   (6,NDFPEM_el)
   real(8) :: ACB2B0   (6,NDFPEM_el)
   real(8) :: AKB2B0   (6,NDFPEM_el)
   real(8) :: AFB2B0   (6          )
   real(8) :: AMB2BNQ0 (6,NQS      )
   real(8) :: ACB2BNQ0 (6,NQS      )
   real(8) :: AKB2BNQ0 (6,NQS      )

   integer :: nbod, nbsub, nb_el, nn_el, nel, nn, j
   integer :: JMAX
   integer :: iq, nq
!  real(8) :: COEFM, COEFK
   integer :: ncm, nbod_cm, nbsub_cm, nel_cm


   write(*,*) 'at the moment LOC_LOADS_int_el is not used'
   stop
!--- Zero the global matrices
   AMB2B   = 0.d0;   ACB2B   = 0.d0;   AKB2B   = 0.d0;   AFB2B   = 0.d0;
   AMB2BNQ = 0.d0;   ACB2BNQ = 0.d0;   AKB2BNQ = 0.d0;

!--- Zero Buoyancy Forces
!  BuoyancyT_el     = 0.d0
!  BuoyancyTtap_el  = 0.d0
!  BuoyancyTside_el = 0.d0


!  do nbod = 1, NBODBT_el
!  do nbsub=1,body(nbod)%NBODSUB_el

!  COEFM  = COEFM_el(nbod)
!  COEFK  = COEFK_el(nbod)
   nb_el  = body(nbod)%NBODTGLB_el(nbsub)
   nn_el  = body(nbod)%NNODTGLB_el(nbsub)

   JMAX   = subbody(nb_el)%NDFPE_el

   do nel = subbody(nb_el)%NTEB_el,1,-1
      nn  = subbody(nb_el)%NODMB       (nel,1)!nnod)
      j   = subbody(nb_el)%NDFPNBACC_el(nn-1 )

      call MATRLOC_Fint       ( AMB2B0  , ACB2B0  , AKB2B0  , AFB2B0     ,&
                                AMB2BNQ0, ACB2BNQ0, AKB2BNQ0             ,&
                                                    nb_el   , nn_el , nel   )
      call FORCLOC_Fint       ( AFB2B0  , AKB2BNQ0,                       &
                                nbod    , nbsub   , nb_el   , nn_el , nel   )

!!    call STIFLOC_tract_Fint ( AKB2B0  , AFB2B0  , nb_el   ,         nel   )

!     call Buoyancy_hyd       ( AFB2B0  , AKB2BNQ0                       ,&
!                               nbod    , nbsub   , nb_el   , nn_el , nel   )

!     call Morison_hyd        ( AMB2B0  , ACB2B0  , AKB2B0  , AFB2B0     ,&
!                               AMB2BNQ0, ACB2BNQ0, AKB2BNQ0             ,&
!                               nbod    , nbsub   , nb_el   , nn_el , nel   )

!!!   call FLOODMATRLOC       ( AMB2B0  , ACB2B0  , AKB2B0  , AFB2B0     ,&
!!!                             AMB2BNQ0, ACB2BNQ0, AKB2BNQ0             ,&
!!!                             nbod    , nbsub   , nb_el   , nn_el  , nel  )

      do ncm      = 1, NCOMA_el
         nbod_cm  = conc_mass(ncm)%NBCOMA_el
         nbsub_cm = conc_mass(ncm)%NSBCOMA_el
         nel_cm   = conc_mass(ncm)%NELCOMA_el

         if (nbod_cm==nbod.and.nbsub_cm==nbsub.and.nel_cm==nel) then

         call CON_MATRLOC_Fint ( AMB2B0  , ACB2B0  , AKB2B0  , AFB2B0,&
                                 AMB2BNQ0, ACB2BNQ0, AKB2BNQ0,        &
                                 ncm     , nb_el   , nn_el   , nel      )
         endif
      enddo !ncm




!     ACLOC0R( 1:6,   1:  JMAX ) =                      COEFM * AMLOC0  ( 1:6, 1:JMAX ) + &
!                                                       COEFK * AKLOC0  ( 1:6, 1:JMAX )
      AMB2B  ( 1:6, j+1:j+JMAX ) = AMB2B  ( 1:6, j+1:j+JMAX ) + AMB2B0  ( 1:6, 1:JMAX )
      ACB2B  ( 1:6, j+1:j+JMAX ) = ACB2B  ( 1:6, j+1:j+JMAX ) + ACB2B0  ( 1:6, 1:JMAX )!+ ACLOC0R( 1:IMAX, 1:JMAX )
      AKB2B  ( 1:6, j+1:j+JMAX ) = AKB2B  ( 1:6, j+1:j+JMAX ) + AKB2B0  ( 1:6, 1:JMAX )
      AFB2B  ( 1:6             ) = AFB2B  ( 1:6             ) + AFB2B0  ( 1:6         )
!                                 - matmul ( ACLOC0R( 1:6, 1:JMAX ), UT1_el( j+1:j+JMAX ) ) 
     
!------ Qs without Rayleigh Damping
      do iq = 1, QSnode(nn_el)%Tot_el(2)
         nq = QSnode(nn_el)%IORD_el(iq)

!     ACLOC0RQ(1:6             ) =                      COEFM * AMLOCNQ0 ( 1:6, nq    ) + &
!                                                       COEFK * AKLOCNQ0 ( 1:6, nq    )
      AMB2BNQ( 1:6, nq         ) = AMB2BNQ( 1:6, nq         ) + AMB2BNQ0 ( 1:6, nq    )
      ACB2BNQ( 1:6, nq         ) = ACB2BNQ( 1:6, nq         ) + ACB2BNQ0 ( 1:6, nq    )
      AKB2BNQ( 1:6, nq         ) = AKB2BNQ( 1:6, nq         ) + AKB2BNQ0 ( 1:6, nq    )
!     AFB2B  ( 1:6             ) = AFB2B  ( 1:6             ) - ACLOC0RQ ( 1:6        ) * UT1_el(j)
      enddo !iq

   enddo !nel


 END Subroutine LOC_LOADS_int_el
!-------------------------------------------------------------------------
!   Subroutine :MATRLOC
! 
!  All Local Matrices [MASS, DAMPING, STIFFNESS, FORCES]
!  with linearly varying properties calculated 
!  with respect to the body FEM local system.
!
!-------------------------------------------------------------------------
 Subroutine MATRLOC_Fint (  AMLOC, ACLOC, AKLOC, AFLOC, AMLOCNQ, ACLOCNQ, AKLOCNQ,&
                            nb_el, nn_el, nel                                       )
!-------------------------------------------------------------------------

 Use Cbeam

   implicit none

   real(8) :: ATDAx2  (3,3)
   real(8) :: ATDDA   (3,3)
   real(8) :: ATDDR   (  3)
   real(8) :: AT      (3,3)
   real(8) :: AKLOC   (6,NDFPEM_el)
   real(8) :: ACLOC   (6,NDFPEM_el)
   real(8) :: AMLOC   (6,NDFPEM_el)
   real(8) :: AFLOC   (6          )
   real(8) :: AMLOCNQ (6,NQS      )
   real(8) :: ACLOCNQ (6,NQS      )
   real(8) :: AKLOCNQ (6,NQS      )
   real(8) :: AKLOC0  (6,NDFPEM_el)
   real(8) :: ACLOC0  (6,NDFPEM_el)
   real(8) :: AMLOC0  (6,NDFPEM_el)
   real(8) :: AFLOC0  (6          )
   real(8) :: UT      (NDFPEM_el )
   real(8) :: UT1     (NDFPEM_el )
   real(8) :: UT2     (NDFPEM_el )
   real(8) :: UTT1    (NEQPEM_el )
   real(8) :: UTT0    (NEQPEM_el )
!f real(8) :: DUT0    (NEQPEM_el )
   real(8) :: AM      (NEQPEM_el,NEQPEM_el)
   real(8) :: AQ      (NEQPEM_el          )
   real(8) :: SHAPE   (NEQPEM_el,NDFPEM_el)
   real(8) :: DSHAPE  (NEQPEM_el,NDFPEM_el)
!f real(8) :: SHAPET  (NDFPEM_el,NEQPEM_el)
!f real(8) :: DSHAPET (NDFPEM_el,NEQPEM_el)
!f real(8) :: SHAPETAM(NDFPEM_el,NEQPEM_el)
   integer :: nb_el, nn_el, nel, IB, NEQPE, NDFPE, k, n1, n2
   real(8) :: ALLOC, PA01 , PA02  , allocwgauss   ,                 F_Grav
   real(8) :: DENS , AMOMX, AMOMZ , RIXX  , RIZZ  , RIXZ  , POLI  , HTA   , HTAL
!f real(8) :: EIXX , EIZZ , EIXZ  , GIT   , EA    , EAX   , EAZ   , TRACT , XCM
!f real(8) :: GAX  , GAZ  , GAX_X , GAZ_Z , PHPX  , PHPZ  , HTA1  , HTA0  , ZCM
   real(8) :: PHPX , PHPZ , HTA1  , HTA0  , XCM   , ZCM   , HTAf


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


   IB      = subbody(nb_el)%IBTYP_el
   NEQPE   = subbody(nb_el)%NEQPE_el
   NDFPE   = subbody(nb_el)%NDFPE_el
   n1      = subbody(nb_el)%INODNEL_el(nel,1)
   n2      = subbody(nb_el)%INODNEL_el(nel,2)
   HTA0    = subbody(nb_el)%HTA_el    (nel)
   ALLOC   = subbody(nb_el)%ALENG_el  (nel)
   PHPX    = subbody(nb_el)%PHIX_el   (nel)
   PHPZ    = subbody(nb_el)%PHIZ_el   (nel)


   call LOCAL_UT_1 (nb_el, nel, UT, UT1, UT2)

!----------------------------!
! IIxAxS           : set AM  !
! IIxAtxF_Grav     : set AQ  !
! IIxAtDDR_AtDDAr0 : sum AQ  !
!----------------------------!

   do k        = 1, IORDER
      HTAL     = ( XGAUSS(k) + 1.d0 )/2.d0   ! 0 < HTAL   < 1 HTAL=HTA/L
      HTA      = HTAL*ALLOC                  ! 0 < HTA    < L
      HTA1     = HTA0 + HTA
      PA01     = 1.d0 - HTAL
      PA02     = HTAL
!aaa
      HTAf     = HTA1
      HTAf     = HTA

      call SHAPEFUNC15 ( HTAL, ALLOC, PHPX, PHPZ, SHAPE, DSHAPE, NDFPE, NEQPE )

      UTT1     = matmul ( SHAPE , UT1 )
      UTT0     = matmul ( SHAPE , UT  )
!f    DUT0     = matmul ( DSHAPE, UT  ) !for TRACT force

!f    SHAPET   = transpose (SHAPE )
!f    DSHAPET  = transpose (DSHAPE)
!------ mass
      DENS     = PA01 * beam_timo(n1)%DENS_el   +  PA02 * beam_timo(n2)%DENS_el
      AMOMX    = PA01 * beam_timo(n1)%AMOMX_el  +  PA02 * beam_timo(n2)%AMOMX_el
      AMOMZ    = PA01 * beam_timo(n1)%AMOMZ_el  +  PA02 * beam_timo(n2)%AMOMZ_el
      RIXX     = PA01 * beam_timo(n1)%RIXX_el   +  PA02 * beam_timo(n2)%RIXX_el
      RIZZ     = PA01 * beam_timo(n1)%RIZZ_el   +  PA02 * beam_timo(n2)%RIZZ_el
      RIXZ     = PA01 * beam_timo(n1)%RIXZ_el   +  PA02 * beam_timo(n2)%RIXZ_el
      POLI     = PA01 * beam_timo(n1)%POLI_el   +  PA02 * beam_timo(n2)%POLI_el
!------ stiffness
!f    EIXX     = PA01 * beam_timo(n1)%EIXX_el   +  PA02 * beam_timo(n2)%EIXX_el
!f    EIZZ     = PA01 * beam_timo(n1)%EIZZ_el   +  PA02 * beam_timo(n2)%EIZZ_el
!f    EIXZ     = PA01 * beam_timo(n1)%EIXZ_el   +  PA02 * beam_timo(n2)%EIXZ_el
!f    GIT      = PA01 * beam_timo(n1)%GIT_el    +  PA02 * beam_timo(n2)%GIT_el
!f    EA       = PA01 * beam_timo(n1)%EA_el     +  PA02 * beam_timo(n2)%EA_el
!f    EAX      = PA01 * beam_timo(n1)%EAX_el    +  PA02 * beam_timo(n2)%EAX_el
!f    EAZ      = PA01 * beam_timo(n1)%EAZ_el    +  PA02 * beam_timo(n2)%EAZ_el
!f    GAX      = PA01 * beam_timo(n1)%GAX_el    +  PA02 * beam_timo(n2)%GAX_el
!f    GAZ      = PA01 * beam_timo(n1)%GAZ_el    +  PA02 * beam_timo(n2)%GAZ_el
!f    GAX_X    = PA01 * beam_timo(n1)%GAX_X_el  +  PA02 * beam_timo(n2)%GAX_X_el
!f    GAZ_Z    = PA01 * beam_timo(n1)%GAZ_Z_el  +  PA02 * beam_timo(n2)%GAZ_Z_el
!------ gravity
      XCM      = PA01 * beam_timo(n1)%XCMD_el   +  PA02 * beam_timo(n2)%XCMD_el
      ZCM      = PA01 * beam_timo(n1)%ZCMD_el   +  PA02 * beam_timo(n2)%ZCMD_el
      F_Grav   =-DENS*GRAV

!-------------------------------------------
!  calculation of the local MASS matrix [M]
!
!
!  Nt*(rdA)*II*S*N*u
!-------------------------------------------
      AM       = 0.d0;
      AM(1,1)  = DENS
      AM(6,1)  = AMOMX
      AM(2,2)  = RIZZ
      AM(3,2)  = AMOMZ
      AM(5,2)  =-RIXZ
      AM(2,3)  = AMOMZ
      AM(3,3)  = DENS
      AM(5,3)  =-AMOMX
      AM(4,4)  = DENS
      AM(6,4)  =-AMOMZ
      AM(2,5)  =-RIXZ
      AM(3,5)  =-AMOMX
      AM(5,5)  = RIXX
      AM(1,6)  = AMOMX
      AM(4,6)  =-AMOMZ
      AM(6,6)  = POLI
      AM(IMDOF_el(4),IMDOF_el(1:6)) = AM(IMDOF_el(4),IMDOF_el(1:6))+AM(IMDOF_el(3),IMDOF_el(1:6))*HTAf
      AM(IMDOF_el(6),IMDOF_el(1:6)) = AM(IMDOF_el(6),IMDOF_el(1:6))-AM(IMDOF_el(1),IMDOF_el(1:6))*HTAf

!f    SHAPETAM = matmul( SHAPET  , AM    )
!f    AMLOC0   = matmul( SHAPETAM, SHAPE )
      AMLOC0   = matmul( AM      , SHAPE )
      AFLOC0   = matmul( AMLOC0  , UT2   )

!----------------------------------------------
!
!  calculation of the local DAMPING matrix [C]
!
!----------------------------------------------
      call IIxAxS (AM, ATDAx2, DENS, AMOMX, AMOMZ, RIXX, RIXZ, RIZZ, 1.d0)

      AM(IMDOF_el(4),IMDOF_el(1:6)) = AM(IMDOF_el(4),IMDOF_el(1:6))+AM(IMDOF_el(3),IMDOF_el(1:6))*HTAf
      AM(IMDOF_el(6),IMDOF_el(1:6)) = AM(IMDOF_el(6),IMDOF_el(1:6))-AM(IMDOF_el(1),IMDOF_el(1:6))*HTAf
!f    SHAPETAM = matmul( SHAPET  , AM    )
!f    ACLOC0   = matmul( SHAPETAM, SHAPE )
      ACLOC0   = matmul( AM      , SHAPE )
      AFLOC0   = AFLOC0 + matmul( ACLOC0  , UT1 )

!----------------------------------------------
!
!  calculation of the local STIFFNESS matrix [K]
!
!----------------------------------------------
      call IIxAxS (AM, ATDDA, DENS, AMOMX, AMOMZ, RIXX, RIXZ, RIZZ, 1.d0)

      AM(IMDOF_el(4),IMDOF_el(1:6)) = AM(IMDOF_el(4),IMDOF_el(1:6))+AM(IMDOF_el(3),IMDOF_el(1:6))*HTAf
      AM(IMDOF_el(6),IMDOF_el(1:6)) = AM(IMDOF_el(6),IMDOF_el(1:6))-AM(IMDOF_el(1),IMDOF_el(1:6))*HTAf
!f    SHAPETAM = matmul( SHAPET , AM    )
!f    AKLOC0   = AKLOC0 + matmul( SHAPETAM, SHAPE )
      AKLOC0   =          matmul( AM      , SHAPE )
      AFLOC0   = AFLOC0 + matmul( AKLOC0  , UT )

!----------------------------------------------
!
!  calculation of the local FORCING matrix [Q]
!
!----------------------------------------------

!-- Forces
! Gravity ct
! A--> AT
   AQ(1:6)  = 0.d0
   
   if (int(Var_Paper(8))==1) &
      call IIxAtxF_Grav ( AQ, AT, F_Grav, XCM, ZCM )               !AQ= -II.At.F_Grav

!------ F:: R --> ATDDR, A --> ATDDA
      call IIxAtDDR_AtDDAr0 ( AQ, ATDDA, ATDDR, DENS, AMOMX, AMOMZ, RIXX, RIXZ, RIZZ, HTA1, 1.d0 )

      AQ(IMDOF_el(4)) = AQ(IMDOF_el(4)) + AQ(IMDOF_el(3))*HTAf
      AQ(IMDOF_el(6)) = AQ(IMDOF_el(6)) - AQ(IMDOF_el(1))*HTAf
!f    AFLOC0  = AFLOC0 + matmul( SHAPET, AQ )
      AFLOC0  = AFLOC0 +                 AQ

!------ Replace
      allocwgauss = 0.5d0*ALLOC*WGAUSS(k)

      AKLOC = AKLOC + allocwgauss*AKLOC0
      ACLOC = ACLOC + allocwgauss*ACLOC0
      AMLOC = AMLOC + allocwgauss*AMLOC0
      AFLOC = AFLOC - allocwgauss*AFLOC0


!------ For all q's [Replace inside]
      call MATRLOC_Q_Fint ( AMLOCNQ    , ACLOCNQ, AKLOCNQ,         UTT0 , UTT1 , F_Grav, &
                            DENS       , AMOMX  , AMOMZ  , RIXX  , RIXZ , RIZZ , HTA1  , &
                            allocwgauss, XCM    , ZCM    , nn_el ,        1.d0 , HTAf      )
   enddo !k


 END Subroutine MATRLOC_Fint
!-------------------------------------------------------------------------------
!Subroutine MATRLOC_Q_Fint ( AMLOCNQ    , ACLOCNQ, AKLOCNQ, SHAPET, UTT0 , UTT1 , F_Grav, &
 Subroutine MATRLOC_Q_Fint ( AMLOCNQ    , ACLOCNQ, AKLOCNQ,         UTT0 , UTT1 , F_Grav, &
                             DENS       , AMOMX  , AMOMZ  , RIXX  , RIXZ , RIZZ , HTA1  , &
                             allocwgauss, XCM    , ZCM    , nn_el ,        RAT  , HTAf      )
!f                           allocwgauss, XCM    , ZCM    , nn_el , NDFPE, RAT              )
!-------------------------------------------------------------------------------

 Use Cbeam

   implicit none

   real(8) :: ATDAx20 (3,3)
   real(8) :: ATDAx21 (3,3)
   real(8) :: ATDDA0  (3,3)
   real(8) :: ATDDA1  (3,3)
   real(8) :: ATDDA2  (3,3)
   real(8) :: AT0     (3,3)
   real(8) :: ATDDR0  (  3)
   real(8) :: ATDDR1  (  3)
   real(8) :: ATDDR2  (  3)
   real(8) :: AMLOCNQ (6,NQS      )
   real(8) :: ACLOCNQ (6,NQS      )
   real(8) :: AKLOCNQ (6,NQS      )
   real(8) :: AMLOCNQ0(6          )
   real(8) :: ACLOCNQ0(6          )
   real(8) :: AKLOCNQ0(6          )
   real(8) :: UTT0    (NEQPEM_el          )
   real(8) :: UTT1    (NEQPEM_el          )
   real(8) :: AM      (NEQPEM_el,NEQPEM_el)
   real(8) :: AQ      (NEQPEM_el          )
!f real(8) :: SHAPET  (NDFPEM_el,NEQPEM_el)
!f integer :: nn_el, n, NDFPE, iq
   integer :: nn_el, n,        iq
   real(8) :: HTA1, allocwgauss, DENS, AMOMX, AMOMZ, RIXX, RIZZ, RIXZ, F_Grav, XCM, ZCM, RAT, HTAf


   do iq = 1, QSnode(nn_el)%Tot_el(2)
      n  = QSnode(nn_el)%IORD_el(iq)

      ATDAx20(1:3,1:3) = transf_mat(nn_el)%ATDA0_el  (iq,1:3,1:3) * 2.d0
      ATDAx21(1:3,1:3) = transf_mat(nn_el)%ATDA1_el  (iq,1:3,1:3) * 2.d0
      ATDDA0 (1:3,1:3) = transf_mat(nn_el)%ATDDA0_el (iq,1:3,1:3)
      ATDDA1 (1:3,1:3) = transf_mat(nn_el)%ATDDA1_el (iq,1:3,1:3)
      ATDDA2 (1:3,1:3) = transf_mat(nn_el)%ATDDA2_el (iq,1:3,1:3)
      AT0    (1:3,1:3) = transf_mat(nn_el)%AT0_el    (iq,1:3,1:3)
      ATDDR0 (1:3    ) = transf_mat(nn_el)%ATDDR0_el (iq,1:3    )
      ATDDR1 (1:3    ) = transf_mat(nn_el)%ATDDR1_el (iq,1:3    )
      ATDDR2 (1:3    ) = transf_mat(nn_el)%ATDDR2_el (iq,1:3    )

!---------------
!-- AKLOCNQ --
!---------------

!----- GRAV AT0
!jim checks for the paper
   AQ(1:6)  = 0.d0
   
   if (int(Var_Paper(8))==1) &
     call IIxAtxF_Grav ( AQ, AT0, F_Grav, XCM, ZCM )               !AQ= -II.At.F_Grav

!------ ATDDA0
      call IIxAxS (AM, ATDDA0, DENS, AMOMX, AMOMZ, RIXX, RIXZ, RIZZ, RAT)

      AQ = AQ + matmul ( AM, UTT0 )

!------ 2xATDA0
      call IIxAxS (AM, ATDAx20, DENS, AMOMX, AMOMZ, RIXX, RIXZ, RIZZ, RAT)

      AQ = AQ + matmul ( AM, UTT1 )

!------ ATDDR0, ATDDA0
      call IIxAtDDR_AtDDAr0 ( AQ, ATDDA0, ATDDR0, DENS, AMOMX, AMOMZ, RIXX, RIXZ, RIZZ, HTA1, RAT )

      AQ(IMDOF_el(4)) = AQ(IMDOF_el(4)) + AQ(IMDOF_el(3))*HTAf
      AQ(IMDOF_el(6)) = AQ(IMDOF_el(6)) - AQ(IMDOF_el(1))*HTAf
!f    AKLOCNQ0 =  matmul( SHAPET  , AQ )
      AKLOCNQ0 =                    AQ

!---------------
!-- ACLOCNQ --
!---------------

!------ ATDDA1
      call IIxAxS (AM, ATDDA1, DENS, AMOMX, AMOMZ, RIXX, RIXZ, RIZZ, RAT)

      AQ = matmul ( AM, UTT0 )

!------ 2xATDA1
      call IIxAxS (AM, ATDAx21, DENS, AMOMX, AMOMZ, RIXX, RIXZ, RIZZ, RAT)

      AQ = AQ + matmul ( AM, UTT1 )

!------ ATDDR1, ATDDA1
      call IIxAtDDR_AtDDAr0 ( AQ, ATDDA1, ATDDR1, DENS, AMOMX, AMOMZ, RIXX, RIXZ, RIZZ, HTA1, RAT )

      AQ(IMDOF_el(4)) = AQ(IMDOF_el(4)) + AQ(IMDOF_el(3))*HTAf
      AQ(IMDOF_el(6)) = AQ(IMDOF_el(6)) - AQ(IMDOF_el(1))*HTAf
!f    ACLOCNQ0 =  matmul( SHAPET  , AQ )
      ACLOCNQ0 =                    AQ

!---------------
!-- AMLOCNQ --
!---------------

!------ ATDDA2
      call IIxAxS (AM, ATDDA2, DENS, AMOMX, AMOMZ, RIXX, RIXZ, RIZZ, RAT)

      AQ = matmul ( AM, UTT0 )

!------ ATDDR2, ATDDA2
      call IIxAtDDR_AtDDAr0 ( AQ, ATDDA2, ATDDR2, DENS, AMOMX, AMOMZ, RIXX, RIXZ, RIZZ, HTA1, RAT )

      AQ(IMDOF_el(4)) = AQ(IMDOF_el(4)) + AQ(IMDOF_el(3))*HTAf
      AQ(IMDOF_el(6)) = AQ(IMDOF_el(6)) - AQ(IMDOF_el(1))*HTAf
!f    AMLOCNQ0 =  matmul( SHAPET  , AQ )
      AMLOCNQ0 =                    AQ

!------ Replace
!f    AKLOCNQ(1:NDFPE,n) = AKLOCNQ(1:NDFPE,n) + allocwgauss * AKLOCNQ0(1:NDFPE)
!f    ACLOCNQ(1:NDFPE,n) = ACLOCNQ(1:NDFPE,n) + allocwgauss * ACLOCNQ0(1:NDFPE)
!f    AMLOCNQ(1:NDFPE,n) = AMLOCNQ(1:NDFPE,n) + allocwgauss * AMLOCNQ0(1:NDFPE)
      AKLOCNQ(1:6    ,n) = AKLOCNQ(1:6    ,n) + allocwgauss * AKLOCNQ0(1:6    )
      ACLOCNQ(1:6    ,n) = ACLOCNQ(1:6    ,n) + allocwgauss * ACLOCNQ0(1:6    )
      AMLOCNQ(1:6    ,n) = AMLOCNQ(1:6    ,n) + allocwgauss * AMLOCNQ0(1:6    )
   enddo !iq


 END Subroutine MATRLOC_Q_Fint
!-------------------------------------------------------------------------
!   Subroutine :FORCLOC
!  
!  Application of external forces, except gravity and hydrodynamic, which are:
!   1. blade Aerodynamics
!   2. tower Drag force
!   3. External local
!   4. external global [like gravity]
!
!-------------------------------------------------------------------------
 Subroutine FORCLOC_Fint (  AFLOC, AKLOCNQ,&
                            nbod , nbsub, nb_el, nn_el, nel )
!-------------------------------------------------------------------------

 Use Cbeam

   implicit none

   real(8) :: AT      (3,3)
   real(8) :: AT0     (3,3)
   real(8) :: AFLOC   (6         )
   real(8) :: AKLOCNQ (6,NQS     )
   real(8) :: AKLOCNQ0(6         )
   real(8) :: AFLOC0  (6         )
   real(8) :: UT      (NDFPEM_el )
   real(8) :: UT1     (NDFPEM_el )
   real(8) :: UT2     (NDFPEM_el )
   real(8) :: UTT1    (NEQPEM_el )
   real(8) :: UTT0    (NEQPEM_el )
!f real(8) :: DUT0    (NEQPEM_el )
   real(8) :: AQ      (NEQPEM_el  )
   real(8) :: AIIA    (NEQPEM_el,6)
   real(8) :: FORC    (          6)
   real(8) :: F_ae    (          6)
   real(8) :: SHAPE   (NEQPEM_el,NDFPEM_el)
   real(8) :: DSHAPE  (NEQPEM_el,NDFPEM_el)
!f real(8) :: SHAPET  (NDFPEM_el,NEQPEM_el)
   integer :: nbod, nbsub, nb_el, nn_el, nel, IB, NEQPE, NDFPE, k, n1, n2
   real(8) :: ALLOC, PA01 , PA02  , allocwgauss   , XC    , ZC    , F_Grav
   real(8) :: DENS ,                                                HTA   , HTAL
   real(8) ::                                                               XCM
   real(8) ::                               PHPX  , PHPZ  , HTA1  , HTA0  , ZCM, HTAf, HTA1b
   real(8) :: DIAM !, ds   , Jacob
   real(8) :: FX, XLoa, ZLoa
   integer :: n, iq


!--- Pass Transformation matrices to local arrays
   AT    (1:3,1:3) = transf_mat(nn_el)%AT_el    (1:3,1:3)


   IB      = subbody(nb_el)%IBTYP_el
   NEQPE   = subbody(nb_el)%NEQPE_el
   NDFPE   = subbody(nb_el)%NDFPE_el
   n1      = subbody(nb_el)%INODNEL_el(nel,1)
   n2      = subbody(nb_el)%INODNEL_el(nel,2)
   HTA0    = subbody(nb_el)%HTA_el    (nel)
   ALLOC   = subbody(nb_el)%ALENG_el  (nel)
   PHPX    = subbody(nb_el)%PHIX_el   (nel)
   PHPZ    = subbody(nb_el)%PHIZ_el   (nel)


   call LOCAL_UT_1 (nb_el, nel, UT, UT1, UT2)


   do k        = 1, IORDER
      HTAL     = ( XGAUSS(k) + 1.d0 )/2.d0   ! 0 < HTAL   < 1 HTAL=HTA/L
      HTA      = HTAL*ALLOC                  ! 0 < HTA    < L
      HTA1     = HTA0 + HTA
      PA01     = 1.d0 - HTAL
      PA02     = HTAL
!aaa
      HTAf     = HTA1
      HTAf     = HTA

      call SHAPEFUNC15 ( HTAL, ALLOC, PHPX, PHPZ, SHAPE, DSHAPE, NDFPE, NEQPE )

      UTT1     = matmul ( SHAPE , UT1 )
      UTT0     = matmul ( SHAPE , UT  )
!f    DUT0     = matmul ( DSHAPE, UT  ) !for TRACT force

!f    SHAPET   = transpose (SHAPE )
!------ gravity
      DENS     = PA01 * beam_timo(n1)%DENS_el   +  PA02 * beam_timo(n2)%DENS_el
      XCM      = PA01 * beam_timo(n1)%XCMD_el   +  PA02 * beam_timo(n2)%XCMD_el
      ZCM      = PA01 * beam_timo(n1)%ZCMD_el   +  PA02 * beam_timo(n2)%ZCMD_el
      F_Grav   =-DENS*GRAV
!------ BEMT
      HTA1b    = PA01 * beam_timo(n1)%XL_el(2)  +  PA02 * beam_timo(n2)%XL_el(2)

!----------------------------------------------
!
!  calculation of the local FORCING matrix [Q]
!
!----------------------------------------------
! Here AQ is set, so no need to zero it!
!jim checks for the paper
   AQ(1:6)  = 0.d0
   
!!!if (int(Var_Paper(8))==1) &
!!!  call IIxAtxF_Grav ( AQ, AT, F_Grav, XCM, ZCM )               !AQ= -II.At.F_Grav

     XLoa = 0.d0 !XCM
     ZLoa = 0.d0 !ZCM

!--- paper: call external loads
   if (Var_Paper(1)/=0.d0) then
      FX       =-   GRAV/GRAVITY_el *(Var_Paper(1) + Var_Paper(2) * dcos(Var_Paper(3)*TIME))
!------ global Fx [flap for stand-still zero pitch]
      if (int(Var_Paper(10))==1) call IIxAtxF_GravX( AQ, AT, FX, XLoa, ZLoa, 1 )          !AQ is added
!------ global Fy [edge for stand-still zero pitch]
      if (int(Var_Paper(10))==2) call IIxAtxF_GravX( AQ, AT, FX, XLoa, ZLoa, 2 )          !AQ is added
!------ global Fz [extention for stand-still zero pitch]
      if (int(Var_Paper(10))==3) call IIxAtxF_GravX( AQ, AT, FX, XLoa, ZLoa, 3 )          !AQ is added
!------ global combined Fx/Fy [flap/ege for stand-still zero pitch]
      if (int(Var_Paper(10))==10) then
         call IIxAtxF_GravX( AQ, AT, FX, XLoa, ZLoa, 1 )          !AQ is added
         call IIxAtxF_GravX( AQ, AT, FX, XLoa, ZLoa, 2 )          !AQ is added
      endif
      if (int(Var_Paper(10))==11) then
!        call IIxAtxF_GravX( AQ, AT, FX, XLoa, ZLoa, 1 )          !AQ is added
!------ global Fy [edge for stand-still zero pitch]
         call IIxAtxF_GravX( AQ, AT, FX, XLoa, ZLoa, 2 )          !AQ is added
         FX=20000.d0
         call IIxAtxF_GravX( AQ, AT, FX, XLoa, ZLoa, 3 )          !AQ is added
      endif
   endif


!------ Aerodynamic forces
      if (IB == 1.and.IAERO_el>0) then
        if     (IAERO_el == 1) then

           call LOCAL_FORC0_ae   ( nb_el, nbod, nbsub, HTA1, HTA1b, F_ae, XC, ZC )

        elseif (IAERO_el >= 1) then

           write(*,*) 'genuvp/cfd is disabled in modal version'; stop


        endif

         call DIAGO(6,AIIA)

         AIIA(2,3) =  XC       ! AIIA(IMDOF_el(4),IMDOF_el(2)) = -ZC
         AIIA(5,3) = -ZC       ! AIIA(IMDOF_el(5),IMDOF_el(1)) =  ZC
         AIIA(6,1) =  ZC       ! AIIA(IMDOF_el(5),IMDOF_el(3)) = -XC
         AIIA(6,4) = -XC       ! AIIA(IMDOF_el(6),IMDOF_el(2)) =  XC

         FORC(IMDOF_el(1:6)) = F_ae(1:6)

         AQ        = AQ - matmul (AIIA,FORC)
      endif

      if ((nbod == NTOWER_el).and.(IAERO_el /= 0)) then
         DIAM  = PA01 * substr_mor(n1)%DIAMET_el  +  PA02 * substr_mor(n2)%DIAMET_el
         call DRAG_tow ( AQ, HTA1, DIAM, nn_el )
      endif

!------ External forces
      XC  = XCM  !0.d0
      ZC  = ZCM  !0.d0

      call DIAGO(6,AIIA)

      AIIA (2,3) =  XC !--new
      AIIA (5,3) = -ZC !--new
      AIIA (6,1) =  ZC
      AIIA (6,4) = -XC

      FORC(IMDOF_el(1:3)) = subbody(nb_el)%FCPL_el(nel,1:3)
      FORC(IMDOF_el(4:6)) = subbody(nb_el)%AMPL_el(nel,1:3)

      AQ      = AQ - matmul( AIIA, FORC )

      AQ(IMDOF_el(4)) = AQ(IMDOF_el(4)) + AQ(IMDOF_el(3))*HTAf
      AQ(IMDOF_el(6)) = AQ(IMDOF_el(6)) - AQ(IMDOF_el(1))*HTAf
!f    AFLOC0  =    + matmul( SHAPET, AQ )
      AFLOC0  =                      AQ

!------ Replace
      allocwgauss = 0.5d0*ALLOC*WGAUSS(k)

      AFLOC = AFLOC - allocwgauss*AFLOC0

      do iq = 1, QSnode(nn_el)%Tot_el(2)
   
         n  = QSnode(nn_el)%IORD_el(iq)
   
         AT0    (1:3,1:3) = transf_mat(nn_el)%AT0_el    (iq,1:3,1:3)
   
!---------------
!-- AKLOCNQ --
!---------------
!-------- GRAV AT0
!jim checks for the paper
         AQ(1:6)  = 0.d0

!!!      if (int(Var_Paper(8))==1) &
!!!        call IIxAtxF_Grav ( AQ, AT0, F_Grav, XCM, ZCM )               !AQ= -II.At.F_Grav

           XLoa = 0.d0 !XCM
           ZLoa = 0.d0 !ZCM

!--------- paper: call external loads
         if (Var_Paper(1)/=0.d0) then
            FX       =-   GRAV/GRAVITY_el *(Var_Paper(1) + Var_Paper(2) * dcos(Var_Paper(3)*TIME))
!------------ global Fx [flap for stand-still zero pitch]
            if (int(Var_Paper(10))==1) call IIxAtxF_GravX( AQ, AT0, FX, XLoa, ZLoa, 1 )          !AQ is added
!------------ global Fy [edge for stand-still zero pitch]
            if (int(Var_Paper(10))==2) call IIxAtxF_GravX( AQ, AT0, FX, XLoa, ZLoa, 2 )          !AQ is added
!------------ global Fz [extention for stand-still zero pitch]
            if (int(Var_Paper(10))==3) call IIxAtxF_GravX( AQ, AT0, FX, XLoa, ZLoa, 3 )          !AQ is added
!------------ global combined Fx/Fy [flap/ege for stand-still zero pitch]
            if (int(Var_Paper(10))==10) then
               call IIxAtxF_GravX( AQ, AT0, FX, XLoa, ZLoa, 1 )          !AQ is added
               call IIxAtxF_GravX( AQ, AT0, FX, XLoa, ZLoa, 2 )          !AQ is added
            endif
            if (int(Var_Paper(10))==11) then
!              call IIxAtxF_GravX( AQ, AT0, FX, XLoa, ZLoa, 1 )          !AQ is added
!------------ global Fy [edge for stand-still zero pitch]
               call IIxAtxF_GravX( AQ, AT0, FX, XLoa, ZLoa, 2 )          !AQ is added
               FX=20000.d0
               call IIxAtxF_GravX( AQ, AT0, FX, XLoa, ZLoa, 3 )          !AQ is added
            endif
         endif

         AQ(IMDOF_el(4)) = AQ(IMDOF_el(4)) + AQ(IMDOF_el(3))*HTAf
         AQ(IMDOF_el(6)) = AQ(IMDOF_el(6)) - AQ(IMDOF_el(1))*HTAf
!f       AKLOCNQ0 =  matmul( SHAPET  , AQ )
         AKLOCNQ0 =                    AQ
!--------- Replace
!        AKLOCNQ(1:NDFPE,n) = AKLOCNQ(1:NDFPE,n) + allocwgauss * AKLOCNQ0(1:NDFPE)
         AKLOCNQ(1:6    ,n) = AKLOCNQ(1:6    ,n) + allocwgauss * AKLOCNQ0(1:6    )
      enddo !iq
   enddo !k


 END Subroutine FORCLOC_Fint
!-------------------------------------------------------------------------
 Subroutine CON_MATRLOC_Fint( AMLOC, ACLOC, AKLOC, AFLOC, AMLOCNQ, ACLOCNQ, AKLOCNQ,&
                              ncal , nb_el, nn_el, nel                                )
!-------------------------------------------------------------------------

 Use Cbeam

   implicit none

   real(8) :: ATDAx2  (3,3)
   real(8) :: ATDDA   (3,3)
   real(8) :: ATDDR   (  3)
   real(8) :: AT      (3,3)
   real(8) :: AKLOC   (6,NDFPEM_el)!NDFPEM_el
   real(8) :: ACLOC   (6,NDFPEM_el)!NDFPEM_el
   real(8) :: AMLOC   (6,NDFPEM_el)!NDFPEM_el
   real(8) :: AFLOC   (6          )!NDFPEM_el
   real(8) :: AMLOCNQ (6,NQS      )!NDFPEM_el
   real(8) :: ACLOCNQ (6,NQS      )!NDFPEM_el
   real(8) :: AKLOCNQ (6,NQS      )!NDFPEM_el
   real(8) :: AKLOC0  (6,NDFPEM_el)!NDFPEM_el
   real(8) :: ACLOC0  (6,NDFPEM_el)!NDFPEM_el
   real(8) :: AMLOC0  (6,NDFPEM_el)!NDFPEM_el
   real(8) :: AFLOC0  (6          )!NDFPEM_el
   real(8) :: UT      (NDFPEM_el )         
   real(8) :: UT1     (NDFPEM_el )         
   real(8) :: UT2     (NDFPEM_el )         
   real(8) :: UTT1    (NEQPEM_el )         
   real(8) :: UTT0    (NEQPEM_el )         
   real(8) :: AM      (NEQPEM_el,NEQPEM_el)
   real(8) :: AQ      (NEQPEM_el          )
   real(8) :: SHAPE   (NEQPEM_el,NDFPEM_el)
   real(8) :: DSHAPE  (NEQPEM_el,NDFPEM_el)
!f real(8) :: SHAPET  (NDFPEM_el,NEQPEM_el)
!f real(8) :: DSHAPET (NDFPEM_el,NEQPEM_el)
!f real(8) :: SHAPETAM(NDFPEM_el,NEQPEM_el)
   real(8) :: ALLOC, y0  , HTAL, Hta, PHPX, PHPZ, RAT
   real(8) :: Diag(3,3)
   real(8) :: Mass  , Sxoff , Syoff , Szoff
   real(8) :: IxxTot, IyyTot, IzzTot, IxyTot   , IxzTot, Iyztot
   real(8) :: Sy0   , Ixy0  , Iyz0  , Iyy0
   real(8) :: Xoff , Yoff, Zoff
   integer :: nb_el, nn_el, nel, NEQPE, NDFPE, ncal, nbsha_el
   real(8) :: HTA1, HTAf


!--- Pass Transformation matrices to local arrays
   ATDAx2(1:3,1:3) = transf_mat(nn_el)%ATDA_el  (1:3,1:3) * 2.d0
   ATDDA (1:3,1:3) = transf_mat(nn_el)%ATDDA_el (1:3,1:3)
   AT    (1:3,1:3) = transf_mat(nn_el)%AT_el    (1:3,1:3)
   ATDDR (1:3    ) = transf_mat(nn_el)%ATDDR_el (1:3    )


   NEQPE   = subbody(nb_el)%NEQPE_el
   NDFPE   = subbody(nb_el)%NDFPE_el

   y0      = subbody(nb_el)%HTA_el   (nel)
   ALLOC   = subbody(nb_el)%ALENG_el (nel)

   PHPX    = subbody(nb_el)%PHIX_el  (nel)
   PHPZ    = subbody(nb_el)%PHIZ_el  (nel)


   call LOCAL_UT_1 (nb_el, nel, UT, UT1, UT2)

   Hta     = conc_mass(ncal)%Hta
   HTAL    = Hta/ALLOC
   HTA1    = y0 + Hta
!aaa
   HTAf    = HTA1
   HTAf    = HTA


   call SHAPEFUNC15 ( HTAL, ALLOC, PHPX, PHPZ, SHAPE, DSHAPE, NDFPE, NEQPE )

   UTT1    = matmul ( SHAPE , UT1 )
   UTT0    = matmul ( SHAPE , UT  )

!f SHAPET  = transpose (SHAPE )
!f DSHAPET = transpose (DSHAPE)


!--- Vars
   Mass    = conc_mass(ncal)%Mass
!  Hta     = conc_mass(ncal)%Hta
   Xoff    = conc_mass(ncal)%Xoff
   Yoff    = conc_mass(ncal)%Yoff
   Zoff    = conc_mass(ncal)%Zoff
   Sxoff   = Mass*Xoff
   Syoff   = Mass*Yoff
   Szoff   = Mass*Zoff
   IxxTot  = conc_mass(ncal)%Ixx + Mass*Xoff**2
   IyyTot  = conc_mass(ncal)%Iyy + Mass*Yoff**2
   IzzTot  = conc_mass(ncal)%Izz + Mass*Zoff**2
   IxyTot  = conc_mass(ncal)%Ixy + Mass*Xoff*Yoff
   IxzTot  = conc_mass(ncal)%Ixz + Mass*Xoff*Zoff
   IyzTot  = conc_mass(ncal)%Iyz + Mass*Yoff*Zoff

   Sy0     = Mass     *(y0+Hta)
   Ixy0    = Mass*Xoff*(y0+Hta)
   Iyy0    = Mass*Yoff*(y0+Hta)
   Iyz0    = Mass*Zoff*(y0+Hta)

   nbsha_el= body(min(NBLADE_el+1,NBODBTWT_el))%NBODTGLB_el (1)

   RAT     = 1.d0
   if ( (nb_el == nbsha_el).and.(y0 < 0.01d0) ) &
   RAT     = RAT_GEAR**2
!-------------------------------------------
!
!  calculation of the local MASS matrix [M]
!
!  Nt*(rdA)*II*S*N*u
!-------------------------------------------
   call DIAGO (3,Diag)

   call Con_IIxAxS ( AM    , Diag  , RAT          ,&
                     Mass  , Sxoff , Syoff , Szoff,&
                     IxxTot, IyyTot, IzzTot       ,&
                     IxyTot, IxzTot, Iyztot       ,&
                     IMDOF_el                        )

   AM(IMDOF_el(4),IMDOF_el(1:6)) = AM(IMDOF_el(4),IMDOF_el(1:6))+AM(IMDOF_el(3),IMDOF_el(1:6))*HTAf
   AM(IMDOF_el(6),IMDOF_el(1:6)) = AM(IMDOF_el(6),IMDOF_el(1:6))-AM(IMDOF_el(1),IMDOF_el(1:6))*HTAf
!f SHAPETAM = matmul( SHAPET  , AM    )
!f AMLOC0   = matmul( SHAPETAM, SHAPE )
   AMLOC0   = matmul(       AM, SHAPE )
   AFLOC0   = matmul( AMLOC0  , UT2   )
!----------------------------------------------
!
!  calculation of the local DAMPING matrix [C]
!
!----------------------------------------------
   call Con_IIxAxS ( AM    , ATDAx2, RAT          ,&
                     Mass  , Sxoff , Syoff , Szoff,&
                     IxxTot, IyyTot, IzzTot       ,&
                     IxyTot, IxzTot, Iyztot       ,&
                     IMDOF_el                        )

   AM(IMDOF_el(4),IMDOF_el(1:6)) = AM(IMDOF_el(4),IMDOF_el(1:6))+AM(IMDOF_el(3),IMDOF_el(1:6))*HTAf
   AM(IMDOF_el(6),IMDOF_el(1:6)) = AM(IMDOF_el(6),IMDOF_el(1:6))-AM(IMDOF_el(1),IMDOF_el(1:6))*HTAf
!f SHAPETAM = matmul( SHAPET  , AM    )
!f ACLOC0   = matmul( SHAPETAM, SHAPE )
   ACLOC0   = matmul(       AM, SHAPE )
   AFLOC0   = AFLOC0 + matmul( ACLOC0  , UT1 )
!----------------------------------------------
!
!  calculation of the local STIFFNESS matrix [K]
!
!----------------------------------------------
   call Con_IIxAxS ( AM    , ATDDA , RAT          ,&
                     Mass  , Sxoff , Syoff , Szoff,&
                     IxxTot, IyyTot, IzzTot       ,&
                     IxyTot, IxzTot, Iyztot       ,&
                     IMDOF_el                        )

   AM(IMDOF_el(4),IMDOF_el(1:6)) = AM(IMDOF_el(4),IMDOF_el(1:6))+AM(IMDOF_el(3),IMDOF_el(1:6))*HTAf
   AM(IMDOF_el(6),IMDOF_el(1:6)) = AM(IMDOF_el(6),IMDOF_el(1:6))-AM(IMDOF_el(1),IMDOF_el(1:6))*HTAf
!f SHAPETAM = matmul( SHAPET , AM    )
!f AKLOC0   = matmul( SHAPETAM, SHAPE )
   AKLOC0   = matmul(       AM, SHAPE )
   AFLOC0   = AFLOC0 + matmul( AKLOC0  , UT )
!----------------------------------------------
!
!  calculation of the local FORCING matrix [Q]
!
!----------------------------------------------
   call Con_IIxAtxF_Grav     ( AQ, AT, Mass, Sxoff, Syoff, Szoff, GRAV, IMDOF_el )   !-- AQ= -II.At.F_Grav

   call Con_IIxAtDDR_AtDDAr0 ( AQ    , ATDDA , ATDDR , RAT   ,&
                               Mass  , Sxoff , Syoff , Szoff ,&
                               IxxTot, IyyTot, IzzTot        ,&
                               IxyTot, IxzTot, Iyztot        ,&
                               Sy0   , Ixy0  , Iyz0  , Iyy0  ,&
                               IMDOF_el                         )

   AQ(IMDOF_el(4)) = AQ(IMDOF_el(4)) + AQ(IMDOF_el(3))*HTAf
   AQ(IMDOF_el(6)) = AQ(IMDOF_el(6)) - AQ(IMDOF_el(1))*HTAf
!f AFLOC0  = AFLOC0 + matmul( SHAPET, AQ )
   AFLOC0  = AFLOC0 +                 AQ

!--- Replace
   AKLOC = AKLOC + AKLOC0 
   ACLOC = ACLOC + ACLOC0 
   AMLOC = AMLOC + AMLOC0 
   AFLOC = AFLOC - AFLOC0 

!----------------------------------------------
! QS
!----------------------------------------------

!---- For all q's
   call CON_MATRLOC_Q_Fint ( AMLOCNQ, ACLOCNQ, AKLOCNQ       ,&
!f                           SHAPET , UTT0   , UTT1          ,&
                                      UTT0   , UTT1          ,&
                             Mass   , Sxoff  , Syoff , Szoff ,&
                             IxxTot , IyyTot , IzzTot        ,&
                             IxyTot , IxzTot , Iyztot        ,&
                             Sy0    , Ixy0   , Iyz0  , Iyy0  ,&
                             nn_el  ,          RAT   , HTAf     )
!f                           nn_el  , NDFPE  , RAT              )


 END Subroutine CON_MATRLOC_Fint
!-------------------------------------------------------------------------------
 Subroutine CON_MATRLOC_Q_Fint ( AMLOCNQ, ACLOCNQ, AKLOCNQ       ,&
!f                               SHAPET , UTT0   , UTT1          ,&
                                          UTT0   , UTT1          ,&
                                 Mass   , Sxoff  , Syoff , Szoff ,&
                                 IxxTot , IyyTot , IzzTot        ,&
                                 IxyTot , IxzTot , Iyztot        ,&
                                 Sy0    , Ixy0   , Iyz0  , Iyy0  ,&
                                 nn_el  ,          RAT   , HTAf     )
!f                               nn_el  , NDFPE  , RAT              )
!-------------------------------------------------------------------------------

 Use Cbeam

   implicit none

   real(8) :: ATDAx20 (3,3)
   real(8) :: ATDAx21 (3,3)
   real(8) :: ATDDA0  (3,3)
   real(8) :: ATDDA1  (3,3)
   real(8) :: ATDDA2  (3,3)
   real(8) :: AT0     (3,3)
   real(8) :: ATDDR0  (  3)
   real(8) :: ATDDR1  (  3)
   real(8) :: ATDDR2  (  3)
   real(8) :: AMLOCNQ (6,NQS) !NDFPEM_el
   real(8) :: ACLOCNQ (6,NQS) !NDFPEM_el
   real(8) :: AKLOCNQ (6,NQS) !NDFPEM_el
   real(8) :: AMLOCNQ0(6    ) !NDFPEM_el
   real(8) :: ACLOCNQ0(6    ) !NDFPEM_el
   real(8) :: AKLOCNQ0(6    ) !NDFPEM_el
   real(8) :: UTT0    (NEQPEM_el          )
   real(8) :: UTT1    (NEQPEM_el          )
   real(8) :: AM      (NEQPEM_el,NEQPEM_el)
   real(8) :: AQ      (NEQPEM_el          )
!f real(8) :: SHAPET  (NDFPEM_el,NEQPEM_el)
!f integer :: nn_el, n, NDFPE, iq
   integer :: nn_el, n,        iq
   real(8) :: RAT
   real(8) :: Mass  , Sxoff , Syoff , Szoff
   real(8) :: IxxTot, IyyTot, IzzTot, IxyTot   , IxzTot, Iyztot
   real(8) :: Sy0   , Ixy0  , Iyz0  , Iyy0
   real(8) :: HTAf


   do iq = 1, QSnode(nn_el)%Tot_el(2)

      n  = QSnode(nn_el)%IORD_el(iq)

      ATDAx20(1:3,1:3) = transf_mat(nn_el)%ATDA0_el  (iq,1:3,1:3) * 2.d0
      ATDAx21(1:3,1:3) = transf_mat(nn_el)%ATDA1_el  (iq,1:3,1:3) * 2.d0
      ATDDA0 (1:3,1:3) = transf_mat(nn_el)%ATDDA0_el (iq,1:3,1:3)
      ATDDA1 (1:3,1:3) = transf_mat(nn_el)%ATDDA1_el (iq,1:3,1:3)
      ATDDA2 (1:3,1:3) = transf_mat(nn_el)%ATDDA2_el (iq,1:3,1:3)
      AT0    (1:3,1:3) = transf_mat(nn_el)%AT0_el    (iq,1:3,1:3)
      ATDDR0 (1:3    ) = transf_mat(nn_el)%ATDDR0_el (iq,1:3    )
      ATDDR1 (1:3    ) = transf_mat(nn_el)%ATDDR1_el (iq,1:3    )
      ATDDR2 (1:3    ) = transf_mat(nn_el)%ATDDR2_el (iq,1:3    )
!---------------
!-- AKLOCNQ --
!---------------
     call Con_IIxAtxF_Grav     ( AQ, AT0, Mass, Sxoff, Syoff, Szoff, GRAV, IMDOF_el )   !-- AQ= -II.At.F_Grav


!------ ATDDA0
      call Con_IIxAxS          ( AM    , ATDDA0 , RAT          ,&
                                 Mass  , Sxoff  , Syoff , Szoff,&
                                 IxxTot, IyyTot , IzzTot       ,&
                                 IxyTot, IxzTot , Iyztot       ,&
                                 IMDOF_el                         )

      AQ = AQ + matmul ( AM, UTT0 )

!------ 2xATDA0
      call Con_IIxAxS           ( AM    , ATDAx20, RAT          ,&
                                  Mass  , Sxoff  , Syoff , Szoff,&
                                  IxxTot, IyyTot , IzzTot       ,&
                                  IxyTot, IxzTot , Iyztot       ,&
                                  IMDOF_el                         )

      AQ = AQ + matmul ( AM, UTT1 )

!------ ATDDR0, ATDDA0
      call Con_IIxAtDDR_AtDDAr0 ( AQ    , ATDDA0, ATDDR0, RAT   ,&
                                  Mass  , Sxoff , Syoff , Szoff ,&
                                  IxxTot, IyyTot, IzzTot        ,&
                                  IxyTot, IxzTot, Iyztot        ,&
                                  Sy0   , Ixy0  , Iyz0  , Iyy0  ,&
                                  IMDOF_el                         )

      AQ(IMDOF_el(4)) = AQ(IMDOF_el(4)) + AQ(IMDOF_el(3))*HTAf
      AQ(IMDOF_el(6)) = AQ(IMDOF_el(6)) - AQ(IMDOF_el(1))*HTAf
!f    AKLOCNQ0 = matmul( SHAPET  , AQ )
      AKLOCNQ0 =                   AQ
!---------------
!-- ACLOCNQ --
!---------------
!------ ATDDA1
      call Con_IIxAxS          ( AM    , ATDDA1 , RAT          ,&
                                 Mass  , Sxoff  , Syoff , Szoff,&
                                 IxxTot, IyyTot , IzzTot       ,&
                                 IxyTot, IxzTot , Iyztot       ,&
                                 IMDOF_el                         )

      AQ = matmul ( AM, UTT0 )

!------ 2xATDA1
      call Con_IIxAxS           ( AM    , ATDAx21, RAT          ,&
                                  Mass  , Sxoff  , Syoff , Szoff,&
                                  IxxTot, IyyTot , IzzTot       ,&
                                  IxyTot, IxzTot , Iyztot       ,&
                                  IMDOF_el                         )

      AQ = AQ + matmul ( AM, UTT1 )

!------ ATDDR1, ATDDA1
      call Con_IIxAtDDR_AtDDAr0 ( AQ    , ATDDA1, ATDDR1, RAT   ,&
                                  Mass  , Sxoff , Syoff , Szoff ,&
                                  IxxTot, IyyTot, IzzTot        ,&
                                  IxyTot, IxzTot, Iyztot        ,&
                                  Sy0   , Ixy0  , Iyz0  , Iyy0  ,&
                                  IMDOF_el                         )

      AQ(IMDOF_el(4)) = AQ(IMDOF_el(4)) + AQ(IMDOF_el(3))*HTAf
      AQ(IMDOF_el(6)) = AQ(IMDOF_el(6)) - AQ(IMDOF_el(1))*HTAf
!f    ACLOCNQ0 = matmul( SHAPET  , AQ )
      ACLOCNQ0 =                   AQ
!---------------
!-- AMLOCNQ --
!---------------
!------ ATDDA2
      call Con_IIxAxS          ( AM    , ATDDA2 , RAT          ,&
                                 Mass  , Sxoff  , Syoff , Szoff,&
                                 IxxTot, IyyTot , IzzTot       ,&
                                 IxyTot, IxzTot , Iyztot       ,&
                                 IMDOF_el                         )

      AQ = matmul ( AM, UTT0 )

!------ ATDDR2, ATDDA2
      call Con_IIxAtDDR_AtDDAr0 ( AQ    , ATDDA2, ATDDR2, RAT   ,&
                                  Mass  , Sxoff , Syoff , Szoff ,&
                                  IxxTot, IyyTot, IzzTot        ,&
                                  IxyTot, IxzTot, Iyztot        ,&
                                  Sy0   , Ixy0  , Iyz0  , Iyy0  ,&
                                  IMDOF_el                         )

      AQ(IMDOF_el(4)) = AQ(IMDOF_el(4)) + AQ(IMDOF_el(3))*HTAf
      AQ(IMDOF_el(6)) = AQ(IMDOF_el(6)) - AQ(IMDOF_el(1))*HTAf
!f    AMLOCNQ0 = matmul( SHAPET  , AQ )
      AMLOCNQ0 =                   AQ

!------ Replace
      AKLOCNQ(1:6,n) = AKLOCNQ(1:6,n) + AKLOCNQ0(1:6) !1:NDFPE
      ACLOCNQ(1:6,n) = ACLOCNQ(1:6,n) + ACLOCNQ0(1:6) !1:NDFPE
      AMLOCNQ(1:6,n) = AMLOCNQ(1:6,n) + AMLOCNQ0(1:6) !1:NDFPE
   enddo !iq


 END Subroutine CON_MATRLOC_Q_Fint
!-------------------------------------------------------------------------
!   Subroutine :STIFLOC_tract_Fint
! 
!
!-------------------------------------------------------------------------
 Subroutine STIFLOC_tract_Fint ( AKLOC, AFLOC, Fy, nb_el, nel )
!-------------------------------------------------------------------------

 Use Cbeam
 Use modal_mod

   implicit none

   real(8) :: AKLOC   (6,NDFPEM_el)
   real(8) :: AFLOC   (6          )
   real(8) :: AKLOC0  (6,NDFPEM_el)
   real(8) :: AFLOC0  (6          )
   real(8) :: UT      (NDFPEM_el )
   real(8) :: UT1     (NDFPEM_el )
   real(8) :: UT2     (NDFPEM_el )
   real(8) :: DUT0    (NEQPEM_el )
   real(8) :: AM      (NEQPEM_el,NEQPEM_el)
   real(8) :: AQ      (NEQPEM_el          )
   real(8) :: SHAPE   (NEQPEM_el,NDFPEM_el)
   real(8) :: DSHAPE  (NEQPEM_el,NDFPEM_el)
!f real(8) :: SHAPET  (NDFPEM_el,NEQPEM_el)
!f real(8) :: SHAPETAM(NDFPEM_el,NEQPEM_el)
   integer :: nb_el, nel, NEQPE, NDFPE, k, n1, n2
   real(8) :: ALLOC, PA01 , PA02  , allocwgauss
   real(8) ::                                                       HTA   , HTAL
   real(8) ::                               EA    , EAX   , EAZ   , TRACT
   real(8) ::                               PHPX  , PHPZ  , HTA1  , HTA0  , HTAf, Fy

!qq
!! return

   NEQPE   = subbody(nb_el)%NEQPE_el
   NDFPE   = subbody(nb_el)%NDFPE_el
   n1      = subbody(nb_el)%INODNEL_el(nel,1)
   n2      = subbody(nb_el)%INODNEL_el(nel,2)
   HTA0    = subbody(nb_el)%HTA_el    (nel)
   ALLOC   = subbody(nb_el)%ALENG_el  (nel)
   PHPX    = subbody(nb_el)%PHIX_el   (nel)
   PHPZ    = subbody(nb_el)%PHIZ_el   (nel)


   call LOCAL_UT_1 (nb_el, nel, UT, UT1, UT2)


   do k        = 1, IORDER
      HTAL     = ( XGAUSS(k) + 1.d0 )/2.d0   ! 0 < HTAL   < 1 HTAL=HTA/L
      HTA      = HTAL*ALLOC                  ! 0 < HTA    < L
      HTA1     = HTA0 + HTA
      PA01     = 1.d0 - HTAL
      PA02     = HTAL
!aaa
      HTAf     = HTA1
      HTAf     = HTA

      call SHAPEFUNC15 ( HTAL, ALLOC, PHPX, PHPZ, SHAPE, DSHAPE, NDFPE, NEQPE )

      DUT0     = matmul ( DSHAPE, UT  ) !for TRACT force

!f    SHAPET   = transpose (SHAPE )
!------ stiffness
!qq
write(*,*) 'this part should be fixed for the modal version of hGAST';stop
      EA       = 0.d0!PA01 * beam_timo(n1)%EA_el     +  PA02 * beam_timo(n2)%EA_el
      EAX      = 0.d0!PA01 * beam_timo(n1)%EAX_el    +  PA02 * beam_timo(n2)%EAX_el
      EAZ      = 0.d0!PA01 * beam_timo(n1)%EAZ_el    +  PA02 * beam_timo(n2)%EAZ_el

!K3 TRACT
      if (Integr_md==1) then
         TRACT    = PA01 * Fy + PA02 * subbody_md(nb_el)%LoadIntQ_L (IMDOF_el(2), nel+1)

         AM       = 0.d0;
         AM(2,1)  =  -TRACT
         AM(5,4)  =   TRACT
      else
!--------- Evaluate Axial Force
         TRACT    = EA*DUT0(3) + EAZ*DUT0(2) - EAX*DUT0(5)

         AM       = 0.d0;
         AM(2,1)  =  -TRACT
         AM(2,2)  =  -EAZ*DUT0(1)
         AM(2,3)  =  -EA *DUT0(1)
         AM(2,5)  =   EAX*DUT0(1)

         AM(5,2)  =   EAZ*DUT0(4)
         AM(5,3)  =   EA *DUT0(4)
         AM(5,4)  =   TRACT
         AM(5,5)  =  -EAX*DUT0(4)
      endif

      AM(IMDOF_el(4),IMDOF_el(1:6)) = AM(IMDOF_el(4),IMDOF_el(1:6))+AM(IMDOF_el(3),IMDOF_el(1:6))*HTAf
      AM(IMDOF_el(6),IMDOF_el(1:6)) = AM(IMDOF_el(6),IMDOF_el(1:6))-AM(IMDOF_el(1),IMDOF_el(1:6))*HTAf

!f    SHAPETAM = matmul( SHAPET  , AM     )
!f    AKLOC0   =+matmul( SHAPETAM, DSHAPE )
      AKLOC0   =+matmul(       AM, DSHAPE )

      AQ       = 0.d0;
      AQ(2)    = -TRACT * DUT0(1)
      AQ(5)    =  TRACT * DUT0(4)
      AQ(IMDOF_el(4)) = AQ(IMDOF_el(4)) + AQ(IMDOF_el(3))*HTAf
      AQ(IMDOF_el(6)) = AQ(IMDOF_el(6)) - AQ(IMDOF_el(1))*HTAf
      
!f    AFLOC0   =+matmul( SHAPET, AQ )
      AFLOC0   =+                AQ

!------ Replace
      allocwgauss = 0.5d0*ALLOC*WGAUSS(k)
      AKLOC = AKLOC + allocwgauss*AKLOC0
      AFLOC = AFLOC - allocwgauss*AFLOC0
   enddo !k


 END Subroutine STIFLOC_tract_Fint
