include 'jacket.f90'
include 'stifcon.f90'
!-----------------------------------------------------------------------------------------------------
!
! -- Subroutine :MATRIX_el
!
!   AM_el , AK_el, AC_el, AQ_el 
!   mass, stifness damping and forcing matrices for the complete structure
!
!----------------------------------------------------------------------
 Subroutine MATRIX_el
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   real(8) :: AMLOC0   (NDFPEM_el,NDFPEM_el)
   real(8) :: ACLOC0   (NDFPEM_el,NDFPEM_el)
   real(8) :: AKLOC0   (NDFPEM_el,NDFPEM_el)
   real(8) :: AFLOC0   (NDFPEM_el          )
   real(8) :: AMLOCNQ0 (NDFPEM_el,NQS      )
   real(8) :: ACLOCNQ0 (NDFPEM_el,NQS      )
   real(8) :: AKLOCNQ0 (NDFPEM_el,NQS      )
   integer :: nbod, nbsub, nb_el, nn_el, nel, nn, i,j
   integer :: IMAX, JMAX, IQRayl
   real(8) :: COEFM, COEFK


!--- Zero the global matrices
   AM_el  = 0.d0;
   AC_el  = 0.d0;
   AK_el  = 0.d0;
   AQ_el  = 0.d0;
!--- Zero local elements matrices
   do nb_el = 1, NBODT_el
      subbody(nb_el)%AMLOC_el   ( :, :, : ) = 0.d0; !nel, ndfpe,ndfpe
      subbody(nb_el)%ACLOC_el   ( :, :, : ) = 0.d0; !nel, ndfpe,ndfpe
      subbody(nb_el)%AKLOC_el   ( :, :, : ) = 0.d0; !nel, ndfpe,ndfpe
      subbody(nb_el)%AFLOC_el   ( :, :    ) = 0.d0; !nel, ndfpe
      if (nb_el>NBODTWT_el.and.ICASE_el==2) cycle
      subbody(nb_el)%AMLOCNQ_el ( :, :, : ) = 0.d0; !nel, ndfpe,nqs
      subbody(nb_el)%ACLOCNQ_el ( :, :, : ) = 0.d0; !nel, ndfpe,nqs
      subbody(nb_el)%AKLOCNQ_el ( :, :, : ) = 0.d0; !nel, ndfpe,nqs
   enddo
!--- Zero Buoyancy Forces
   BuoyancyT_el     = 0.d0
   BuoyancyTtap_el  = 0.d0
   BuoyancyTside_el = 0.d0


            IQRayl= 1 !0:all qs Rayleigh, 1:body qs with Rayleigh, 2:no Rayleigh to qs
   do       nbod  = 1, NBODBT_el
            COEFM = COEFM_el(nbod)
            COEFK = COEFK_el(nbod)
      do    nbsub = 1, body   (nbod)%NBODSUB_el
            nb_el =    body   (nbod)%NBODTGLB_el(nbsub)
            nn_el =    body   (nbod)%NNODTGLB_el(nbsub)
            IMAX  =    subbody(nb_el)%NDFPE_el
            JMAX  =    subbody(nb_el)%NDFPE_el
         do nel   = 1, subbody(nb_el)%NTEB_el
            nn    =                                subbody(nb_el)%NODMB       (nel,1)!nnod)
            i     = ACCU%NDFTACC_el(nb_el-1)   +   subbody(nb_el)%NDFPNBACC_el(nn-1 )
            j     = i

!------------ 1st order Timoshenko beam
            call MATRLOC             ( AMLOC0  , ACLOC0  , AKLOC0  , AFLOC0     ,&
                                       AMLOCNQ0, ACLOCNQ0, AKLOCNQ0             ,&
                                                           nb_el   , nn_el , nel   )
            call FORCLOC             ( AFLOC0  , AKLOCNQ0,                       &
                                       nbod    , nbsub   , nb_el   , nn_el , nel   )

            call Buoyancy_hyd        ( AFLOC0  , AKLOCNQ0                       ,&
                                       nbod    , nbsub   , nb_el   , nn_el , nel   )

            call Morison_hyd         ( AMLOC0  , ACLOC0  , AKLOC0  , AFLOC0     ,&
                                       AMLOCNQ0, ACLOCNQ0, AKLOCNQ0             ,&
                                       nbod    , nbsub   , nb_el   , nn_el , nel   )

            call FLOODMATRLOC        ( AMLOC0  , ACLOC0  , AKLOC0  , AFLOC0     ,&
                                       AMLOCNQ0, ACLOCNQ0, AKLOCNQ0             ,&
                                       nbod    , nbsub   , nb_el   , nn_el , nel   )

            call STORE_LOC_MATRIX_el ( AMLOC0  , ACLOC0  , AKLOC0  , AFLOC0     ,&
                                       AMLOCNQ0, ACLOCNQ0, AKLOCNQ0,             &
                                       nb_el   , nn_el   , nel                     )

            call ASSEMBLY_MATRIX_el  ( AMLOC0  , ACLOC0  , AKLOC0  , AFLOC0        ,&
                                       AMLOCNQ0, ACLOCNQ0, AKLOCNQ0, COEFM , COEFK ,&
                                                 nn_el   , i       , j     , IQRayl,&
                                       IMAX    , JMAX    , NDFPEM_el, NDFPEM_el       )
         enddo !nel
      enddo !nbsub
   enddo !nbod


!--- Hydrodynamic Force [Buoyancy, Morison] for end surfaces [and free flooded members ends]
   call BuoyancyTap

!--- Add the drag on the nacelle
   call DRAG_nac


 END Subroutine MATRIX_el
!-----------------------------------------------------------------------------------------------------
!
! -- Subroutine :CON_MATRIX_el
!
!----------------------------------------------------------------------
 Subroutine CON_MATRIX_el
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   real(8) :: AMLOC0   (NDFPEM_el,NDFPEM_el)
   real(8) :: ACLOC0   (NDFPEM_el,NDFPEM_el)
   real(8) :: AKLOC0   (NDFPEM_el,NDFPEM_el)
   real(8) :: AFLOC0   (NDFPEM_el          )
   real(8) :: AMLOCNQ0 (NDFPEM_el,NQS      )
   real(8) :: ACLOCNQ0 (NDFPEM_el,NQS      )
   real(8) :: AKLOCNQ0 (NDFPEM_el,NQS      )
   integer :: ncm  , nbod, nbsub, nb_el, nn_el, nel, nn, i,j
   integer :: IMAX , JMAX, IQRayl
   real(8) :: COEFM, COEFK


!--- AM_el , AK_el, AC_el        
!--- mass, stifness and damping matrices of concentrated masses
!--- for the complete structure

   IQRayl   = 1 !0:all qs Rayleigh, 1:body qs with Rayleigh, 2:no Rayleigh to qs
   do ncm   = 1, NCOMA_el
      nbod  = conc_mass(ncm)%NBCOMA_el
      nbsub = conc_mass(ncm)%NSBCOMA_el
      nel   = conc_mass(ncm)%NELCOMA_el
      nb_el = body(nbod)%NBODTGLB_el (nbsub)
      nn_el = body(nbod)%NNODTGLB_el (nbsub)
      nn    =                                subbody(nb_el)%NODMB       (nel,1)!nnod)
      i     = ACCU%NDFTACC_el(nb_el-1)   +   subbody(nb_el)%NDFPNBACC_el(nn-1 )
      j     = i
      COEFM = COEFM_el(nbod)
      COEFK = COEFK_el(nbod)
      IMAX  = subbody(nb_el)%NDFPE_el
      JMAX  = subbody(nb_el)%NDFPE_el

      AMLOC0   = 0.d0; ACLOC0   = 0.d0; AKLOC0   = 0.d0; AFLOC0  = 0.d0;
      AMLOCNQ0 = 0.d0; ACLOCNQ0 = 0.d0; AKLOCNQ0 = 0.d0;

      call CON_MATRLOC         ( AMLOC0  , ACLOC0  , AKLOC0  , AFLOC0,&
                                 AMLOCNQ0, ACLOCNQ0, AKLOCNQ0,        &
                                 ncm     , nb_el   , nn_el   , nel      )

      call STORE_LOC_MATRIX_el ( AMLOC0  , ACLOC0  , AKLOC0  , AFLOC0,&
                                 AMLOCNQ0, ACLOCNQ0, AKLOCNQ0,        &
                                 nb_el   , nn_el   , nel                )

      call ASSEMBLY_MATRIX_el  ( AMLOC0  , ACLOC0  , AKLOC0  , AFLOC0        ,&
                                 AMLOCNQ0, ACLOCNQ0, AKLOCNQ0, COEFM , COEFK ,&
                                           nn_el   , i       , j     , IQRayl,&
                                 IMAX    , JMAX    , NDFPEM_el, NDFPEM_el       )
   enddo !ncm


 END Subroutine CON_MATRIX_el
!----------------------------------------------------------------------
!
! -- Subroutine :STORE_LOC_MATRIX_el -----
!
!    Store Local Matrices for Communicating / Writing Loads (and other use)
!
!----------------------------------------------------------------------
 Subroutine STORE_LOC_MATRIX_el ( AMLOC0  , ACLOC0  , AKLOC0  , AFLOC0,&
                                  AMLOCNQ0, ACLOCNQ0, AKLOCNQ0,        &
                                  nb_el   , nn_el   , nel                )
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   real(8) :: AMLOC0   (NDFPEM_el,NDFPEM_el)
   real(8) :: ACLOC0   (NDFPEM_el,NDFPEM_el)
   real(8) :: AKLOC0   (NDFPEM_el,NDFPEM_el)
   real(8) :: AFLOC0   (NDFPEM_el          )
   real(8) :: AMLOCNQ0 (NDFPEM_el,NQS      )
   real(8) :: ACLOCNQ0 (NDFPEM_el,NQS      )
   real(8) :: AKLOCNQ0 (NDFPEM_el,NQS      )
   integer :: nb_el, nn_el, nel, iq, nq


      subbody(nb_el)%AMLOC_el   ( nel, 1:NDFPEM_el, 1:NDFPEM_el ) = subbody(nb_el)%AMLOC_el   ( nel, 1:NDFPEM_el, 1:NDFPEM_el ) + AMLOC0   ( 1:NDFPEM_el, 1:NDFPEM_el )
      subbody(nb_el)%ACLOC_el   ( nel, 1:NDFPEM_el, 1:NDFPEM_el ) = subbody(nb_el)%ACLOC_el   ( nel, 1:NDFPEM_el, 1:NDFPEM_el ) + ACLOC0   ( 1:NDFPEM_el, 1:NDFPEM_el )
      subbody(nb_el)%AKLOC_el   ( nel, 1:NDFPEM_el, 1:NDFPEM_el ) = subbody(nb_el)%AKLOC_el   ( nel, 1:NDFPEM_el, 1:NDFPEM_el ) + AKLOC0   ( 1:NDFPEM_el, 1:NDFPEM_el )
      subbody(nb_el)%AFLOC_el   ( nel, 1:NDFPEM_el              ) = subbody(nb_el)%AFLOC_el   ( nel, 1:NDFPEM_el              ) + AFLOC0   ( 1:NDFPEM_el              )

   do iq = 1, QSnode(nn_el)%Tot_el(2)
      nq = QSnode(nn_el)%IORD_el(iq)

      subbody(nb_el)%AMLOCNQ_el ( nel, 1:NDFPEM_el,          iq ) = subbody(nb_el)%AMLOCNQ_el ( nel, 1:NDFPEM_el,          iq ) + AMLOCNQ0 ( 1:NDFPEM_el,          nq )
      subbody(nb_el)%ACLOCNQ_el ( nel, 1:NDFPEM_el,          iq ) = subbody(nb_el)%ACLOCNQ_el ( nel, 1:NDFPEM_el,          iq ) + ACLOCNQ0 ( 1:NDFPEM_el,          nq )
      subbody(nb_el)%AKLOCNQ_el ( nel, 1:NDFPEM_el,          iq ) = subbody(nb_el)%AKLOCNQ_el ( nel, 1:NDFPEM_el,          iq ) + AKLOCNQ0 ( 1:NDFPEM_el,          nq )
   enddo !iq


 END Subroutine STORE_LOC_MATRIX_el
!----------------------------------------------------------------------
!
! -- Subroutine :ASSEMBLY_MATRIX_el -----
!
!    Assembly Local Matrices to Global (M, C, K, Q).
!    Called from BodyLoads and Matrix Subroutines
!
!    IQRayl: 0: all qs Rayleigh, 1: body qs with Rayleigh, 2: no Rayleigh to qs
!
!----------------------------------------------------------------------
 Subroutine ASSEMBLY_MATRIX_el ( AMLOC0  , ACLOC0  , AKLOC0  , AFLOC0        , &
                                 AMLOCNQ0, ACLOCNQ0, AKLOCNQ0, COEFM , COEFK , &
                                           nn_el   , i       , j     , IQRayl, &
                                           IMAX    , JMAX    , Idi   , Jdi       )
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   integer :: IMAX,JMAX, Idi,Jdi, nn_el, i, j, IQRayl
   integer :: iq, nq, IQr
   real(8) :: AMLOC0   (Idi ,Jdi ) !  real(8) :: AMLOC0   (IMAX,JMAX) !NDFPEM_el)
   real(8) :: ACLOC0   (Idi ,Jdi ) !  real(8) :: ACLOC0   (IMAX,JMAX) !NDFPEM_el)
   real(8) :: AKLOC0   (Idi ,Jdi ) !  real(8) :: AKLOC0   (IMAX,JMAX) !NDFPEM_el)
   real(8) :: AFLOC0   (Idi      ) !  real(8) :: AFLOC0   (IMAX     )
   real(8) :: AMLOCNQ0 (Idi ,NQS ) !  real(8) :: AMLOCNQ0 (IMAX,NQS )
   real(8) :: ACLOCNQ0 (Idi ,NQS ) !  real(8) :: ACLOCNQ0 (IMAX,NQS )
   real(8) :: AKLOCNQ0 (Idi ,NQS ) !  real(8) :: AKLOCNQ0 (IMAX,NQS )
   real(8) :: ACLOC0R  (IMAX,JMAX)                                    !NDFPEM_el)
   real(8) :: ACLOC0RQ (IMAX     )
   real(8) :: COEFM, COEFK


                  IQr = 0
   if (IQRayl/=0) IQr = QSnode(nn_el)%Tot_el(IQRayl)


   ACLOC0R( 1:  IMAX,   1:  JMAX ) =                           COEFM * AMLOC0  ( 1:IMAX, 1:JMAX)  + &
                                                               COEFK * AKLOC0  ( 1:IMAX, 1:JMAX)
   AM_el( i+1:i+IMAX, j+1:j+JMAX ) = AM_el( i+1:i+IMAX, j+1:j+JMAX ) + AMLOC0  ( 1:IMAX, 1:JMAX )
   AC_el( i+1:i+IMAX, j+1:j+JMAX ) = AC_el( i+1:i+IMAX, j+1:j+JMAX ) + ACLOC0  ( 1:IMAX, 1:JMAX ) + ACLOC0R( 1:IMAX, 1:JMAX )
   AK_el( i+1:i+IMAX, j+1:j+JMAX ) = AK_el( i+1:i+IMAX, j+1:j+JMAX ) + AKLOC0  ( 1:IMAX, 1:JMAX )
   AQ_el( i+1:i+IMAX             ) = AQ_el( i+1:i+IMAX             ) + AFLOC0  ( 1:IMAX         ) - &
                                     matmul ( ACLOC0R( 1:IMAX, 1:JMAX ), UT1_el( j+1:j+JMAX ) ) 

!--- Qs without Rayleigh Damping
   do iq = 1, IQr
      nq = QSnode(nn_el)%IORD_el(iq)
      j  = NDFBT_el + nq

      AM_el( i+1:i+IMAX, j ) = AM_el( i+1:i+IMAX, j ) + AMLOCNQ0 ( 1:IMAX, nq )
      AC_el( i+1:i+IMAX, j ) = AC_el( i+1:i+IMAX, j ) + ACLOCNQ0 ( 1:IMAX, nq )
      AK_el( i+1:i+IMAX, j ) = AK_el( i+1:i+IMAX, j ) + AKLOCNQ0 ( 1:IMAX, nq )
   enddo !iq


!--- Body's Qs with Rayleigh Damping
   do iq = IQr + 1, QSnode(nn_el)%Tot_el(2)
      nq = QSnode(nn_el)%IORD_el(iq)
      j  = NDFBT_el + nq

      ACLOC0RQ(1:  IMAX    ) =                  COEFM * AMLOCNQ0 ( 1:IMAX, nq ) + &
                                                COEFK * AKLOCNQ0 ( 1:IMAX, nq )
      AM_el( i+1:i+IMAX, j ) = AM_el( i+1:i+IMAX, j ) + AMLOCNQ0 ( 1:IMAX, nq )
      AC_el( i+1:i+IMAX, j ) = AC_el( i+1:i+IMAX, j ) + ACLOCNQ0 ( 1:IMAX, nq ) + ACLOC0RQ ( 1:IMAX )
      AK_el( i+1:i+IMAX, j ) = AK_el( i+1:i+IMAX, j ) + AKLOCNQ0 ( 1:IMAX, nq )
      AQ_el( i+1:i+IMAX    ) = AQ_el( i+1:i+IMAX    ) - ACLOC0RQ ( 1:IMAX     ) * UT1_el(j)
   enddo !iq


 END Subroutine ASSEMBLY_MATRIX_el
!-----------------------------------------------------------------------
!
!-- Subroutine :BOUNDCO_el
!
!   Apply BOUNDary COnditions
!
!
!----------------------------------------------------------------------
 Subroutine BOUNDCO_el(IP)
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   integer, intent (in) :: IP
   integer              :: nb_el, NDFB0


!--- Boundary Conditions for Movements and Rotations of each subbody
!--- Zero the 1st node of 1st element 6 dofs
   if (ITYPE_FOUND_el>0) then
      write(*,*)'error in BOUNDC0_el, because of the foundation'
      stop
   endif

   do nb_el=1,NBODTWT_el
     !iel    = 1 
     !inodel = 1
     !inode  =                            subbody(nb_el)%NODMB       (iel,inodel)
     !NDFB0  = ACCU%NDFTACC_el(nb_el-1) + subbody(nb_el)%NDFPNBACC_el(inode-1   )
      NDFB0  = ACCU%NDFTACC_el(nb_el-1)

      call ZERO_DOF ( IP, NDFB0, 1, 6 )
   enddo! nb_el


 END Subroutine BOUNDCO_el
!----------------------------------------------------------------------
 Subroutine ZERO_DOF ( IP, ndof0, ndof1, ndof2 )
!----------------------------------------------------------------------

 Use Cbeam

   Implicit none

   integer, Intent (In) :: IP, ndof0, ndof1, ndof2
   integer              :: ndof, ii


   AM_el ( ndof0+ndof1:ndof0+ndof2 , 1:NDFT_el ) = 0.d0
   AC_el ( ndof0+ndof1:ndof0+ndof2 , 1:NDFT_el ) = 0.d0
   AK_el ( ndof0+ndof1:ndof0+ndof2 , 1:NDFT_el ) = 0.d0

   if (IP == 0) then
      do ndof = ndof1, ndof2
         ii =  ndof0 + ndof
         AK_el(ii,ii) = 1.d0
         AQ_el(ii   ) = -UT_el (ii)
      enddo !ndof
   else
      do ndof = ndof1, ndof2
         ii =  ndof0 + ndof
         AM_el(ii,ii) = 1.d0
         AQ_el(ii   ) = -UT2_el(ii)
      enddo !ndof
   endif


 END Subroutine ZERO_DOF

!----------------------------------------------------------------------
!
!--- Boundary Conditions for Movements and Rotations of the last node
!    of the blade in case of VAWT.
!
!----------------------------------------------------------------------
 Subroutine BOUNDCO_VAWT_el(IP)
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   integer, intent (in) :: IP
   integer              :: nbod,nbsub,nb_el,nn_el,nbc,iel,inodel,inode,ndof0
   integer              :: nn_elcon !,nb_elcon,nelcon,nodcon
   integer              :: iq,nq,i
   real(8)              :: Ainit(3,3) !, Alhs(3,3),Arhs(3,3)


   if (IAPPL_el /= 1) return!0:HAWT, 1:VAWT

   do nbod     = 1, NBLADE_el
      nbsub    = body(nbod)%NBODSUB_el                                  !last subbody of each blade
      nb_el    = body(nbod)%NBODTGLB_el(nbsub  )
      nn_el    = body(nbod)%NNODTGLB_el(nbsub+1)                        !last node    of subbody
      nbc      = boundco(nb_el)%nbcpb                                   !last b.c.
      iel      = subbody(nb_el)%NTEB_el !boundco (nb_el)%nel   (nbc)
      inodel   = NNPE_el                !boundco (nb_el)%nod   (nbc)
      inode    =                            subbody(nb_el)%NODMB       (iel,inodel)
      ndof0    = ACCU%NDFTACC_el(nb_el-1) + subbody(nb_el)%NDFPNBACC_el(inode-1   )

      Ainit(1:3,1:3) = boundco (nb_el)%Atrans1(nbc,1:3,1:3) !matmul (transf_mat(nnT2_el)%AT_el(1:3,1:3), transf_mat(nn_el)%A_el(1:3,1:3))

!     nb_elcon = body(NTOWER_el)%NBODTGLB_el ( body(NTOWER_el)%NBODSUB_el  )               !boundco (nb_el)%nbcon (nbc)
      nn_elcon = body(NTOWER_el)%NNODTGLB_el ( body(NTOWER_el)%NBODSUB_el+1)                                                      !last node of tower: ttop                                     
!     nelcon   = subbody(body(NTOWER_el)%NBODTGLB_el(body(NTOWER_el)%NBODSUB_el))%NTEB_el  !boundco (nb_el)%nelcon(nbc)
!     nodcon   = NNPE_el                                                                   !boundco (nb_el)%nodcon(nbc)
                !boundco (nb_el)%indx  (nbc,1:6) = 1


!--- 1. set deflections: RG_bld_end = RG_ttop
      do i=1,3
            AM_el ( ndof0+IMDOF_el(i), 1:NDFT_el ) = 0.d0;
            AC_el ( ndof0+IMDOF_el(i), 1:NDFT_el ) = 0.d0;
            AK_el ( ndof0+IMDOF_el(i), 1:NDFT_el ) = 0.d0;
      enddo

      do iq   = 1, QSnode(nn_elcon)%Tot_el(2)
         nq   = QSnode(nn_el)%IORD_el(iq)

         do i=1,3
            AK_el ( ndof0+IMDOF_el(i) , NDFBT_el+nq ) =  transf_mat(nn_el)%R0_el(iq,i) - transf_mat(nn_elcon)%R0_el(iq,i)
         enddo
      enddo!iq 

      do iq   = QSnode(nn_elcon)%Tot_el(2)+1, QSnode(nn_el)%Tot_el(2) -3 !3 last rotations do not contribute to the R0
         nq   = QSnode(nn_el)%IORD_el(iq)

         do i=1,3
            AK_el ( ndof0+IMDOF_el(i), NDFBT_el+nq ) =  transf_mat(nn_el)%R0_el(iq,i)
         enddo
      enddo!iq 

         do i=1,3
            AQ_el ( ndof0+IMDOF_el(i)              ) = -transf_mat(nn_el)%R_el (   i) + transf_mat(nn_elcon)%R_el (   i ) 
         enddo


!--- 2. set roations: {AT_t.A_b}_init = AT_t.A_b = AT_t . A_b0 . A_th_last_linear
!                   ==> [A_th_last_linear] = AT_b . A_t . {AT_t.A_b}_init
      do i=4,6
            AM_el ( ndof0+IMDOF_el(i), 1:NDFT_el ) = 0.d0;
            AC_el ( ndof0+IMDOF_el(i), 1:NDFT_el ) = 0.d0;
            AK_el ( ndof0+IMDOF_el(i), 1:NDFT_el ) = 0.d0;
      enddo

   !!!nn_el   = body(nbod)%NNODTGLB_el(nbsub)  !1st node of last subbody
   !!!do iq   = 1, QSnode(nn_elcon)%Tot_el(2)
   !!!   nq   = QSnode(nn_el)%IORD_el(iq)

   !!!   Alhs(1:3,1:3) = - matmul ( matmul( transf_mat(nn_el)%AT_el (   1:3,1:3), transf_mat(nn_elcon)%A0_el(iq,1:3,1:3) ) + &
   !!!                              matmul( transf_mat(nn_el)%AT0_el(iq,1:3,1:3), transf_mat(nn_elcon)%A_el (   1:3,1:3) )    , Ainit(1:3,1:3) )

!! !!!   AK_el ( ndof0+IMDOF_el(4) , NDFBT_el+nq ) = Alhs(3,2)
   !!!   AK_el ( ndof0+IMDOF_el(4) , NDFBT_el+nq ) = (Alhs(3,2)-Alhs(2,3)) / 2.d0
   !!!   AK_el ( ndof0+IMDOF_el(5) , NDFBT_el+nq ) = (Alhs(1,3)-Alhs(3,1)) / 2.d0
   !!!   AK_el ( ndof0+IMDOF_el(6) , NDFBT_el+nq ) = (Alhs(2,1)-Alhs(1,2)) / 2.d0
   !!!enddo!iq 

   !!!do iq   = QSnode(nn_elcon)%Tot_el(2)+1, QSnode(nn_el)%Tot_el(2)
   !!!   nq   = QSnode(nn_el)%IORD_el(iq)

   !!!   Alhs(1:3,1:3) = - matmul ( matmul( transf_mat(nn_el)%AT0_el(iq,1:3,1:3), transf_mat(nn_elcon)%A_el (   1:3,1:3) )    , Ainit(1:3,1:3) )

!! !!!   AK_el ( ndof0+IMDOF_el(4) , NDFBT_el+nq ) = Alhs(3,2)
   !!!   AK_el ( ndof0+IMDOF_el(4) , NDFBT_el+nq ) = (Alhs(3,2)-Alhs(2,3)) / 2.d0
   !!!   AK_el ( ndof0+IMDOF_el(5) , NDFBT_el+nq ) = (Alhs(1,3)-Alhs(3,1)) / 2.d0
   !!!   AK_el ( ndof0+IMDOF_el(6) , NDFBT_el+nq ) = (Alhs(2,1)-Alhs(1,2)) / 2.d0
   !!!enddo!iq 
   !!!do i=4,6
!! !!!do i=4,4
   !!!   AK_el ( ndof0+IMDOF_el(i), ndof0+IMDOF_el(i) ) = 1.d0 
   !!!enddo!i

   !!!   Arhs(1:3,1:3) =   matmul ( matmul( transf_mat(nn_el)%AT_el (   1:3,1:3), transf_mat(nn_elcon)%A_el (   1:3,1:3) )    , Ainit(1:3,1:3) )
!! !!!   AQ_el ( ndof0+IMDOF_el(4)              ) =  Arhs(3,2)                   - UT_el(ndof0+IMDOF_el(4))
   !!!   AQ_el ( ndof0+IMDOF_el(4)              ) = (Arhs(3,2)-Arhs(2,3)) / 2.d0 - UT_el(ndof0+IMDOF_el(4))
   !!!   AQ_el ( ndof0+IMDOF_el(5)              ) = (Arhs(1,3)-Arhs(3,1)) / 2.d0 - UT_el(ndof0+IMDOF_el(5))
   !!!   AQ_el ( ndof0+IMDOF_el(6)              ) = (Arhs(2,1)-Arhs(1,2)) / 2.d0 - UT_el(ndof0+IMDOF_el(6))

!------ zero rotations
      do i = 4, 6
         if (IP==0) then
            AK_el ( ndof0+IMDOF_el(i), ndof0+IMDOF_el(i) ) = 1.d0
            AQ_el ( ndof0+IMDOF_el(i)                    ) = - UT_el (ndof0+IMDOF_el(i))
         else!if (IP==1)&
            AM_el ( ndof0+IMDOF_el(i), ndof0+IMDOF_el(i) ) = 1.d0
            AQ_el ( ndof0+IMDOF_el(i)                    ) = - UT2_el(ndof0+IMDOF_el(i))
         endif
      enddo

!!    write(*,*         ) 'Arhs',nbod
!!    write(*,'(3f15.7)') Arhs(1  ,1:3)
!!    write(*,'(3f15.7)') Arhs(2  ,1:3)
!!    write(*,'(3f15.7)') Arhs(3  ,1:3)
   enddo !nbod


 END Subroutine BOUNDCO_VAWT_el
!----------------------------------------------------------------------
 Subroutine BOUNDCO_VAWT2_el(IP)
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   integer, intent (in) :: IP
   integer              :: nbod,nbsub,nb_el,nn_el,nbc,iel,inodel,inode,ndof0
   integer              :: nn_elcon !,nb_elcon,nelcon,nodcon
   integer              :: iq,nq,j,i
   real(8)              :: Hlen, Y0_U(3)


   if (IAPPL_el /= 1) return!0:HAWT, 1:VAWT

!--- Boundary Conditions for Movements and Rotations of each subbody
!--- Zero the 1st node of 1st element 6 dofs

!--- 1. set deflections RG_bld_end = RG_ttop
   do nbod     = 1, NBLADE_el
      nbsub    = body(nbod)%NBODSUB_el                                  !last subbody of each blade
      nb_el    = body(nbod)%NBODTGLB_el(nbsub)
      nn_el    = body(nbod)%NNODTGLB_el(nbsub)                        !1st  node    of last subbody
      nbc      = boundco(nb_el)%nbcpb                                   !last b.c.
      iel      = subbody(nb_el)%NTEB_el !boundco (nb_el)%nel   (nbc)
      inodel   = NNPE_el                !boundco (nb_el)%nod   (nbc)
      inode    =                            subbody(nb_el)%NODMB       (iel,inodel)
      ndof0    = ACCU%NDFTACC_el(nb_el-1) + subbody(nb_el)%NDFPNBACC_el(inode-1   )
      Hlen     = subbody(nb_el)%ALENGB_el
     do i=1,3
      Y0_U(i)  = UT_el(ndof0+IMDOF_el(i))
     enddo
      Y0_U(2)  = Y0_U(2)+Hlen

!     nb_elcon = body(NTOWER_el)%NBODTGLB_el ( body(NTOWER_el)%NBODSUB_el  )               !boundco (nb_el)%nbcon (nbc)
      nn_elcon = body(NTOWER_el)%NNODTGLB_el ( body(NTOWER_el)%NBODSUB_el+1)                                                      !last node of tower: ttop                                     
!     nelcon   = subbody(body(NTOWER_el)%NBODTGLB_el(body(NTOWER_el)%NBODSUB_el))%NTEB_el  !boundco (nb_el)%nelcon(nbc)
!     nodcon   = NNPE_el                                                                   !boundco (nb_el)%nodcon(nbc)
                !boundco (nb_el)%indx  (nbc,1:6) = 1
      do i     = 1, 3
            AM_el ( ndof0+IMDOF_el(i), 1:NDFT_el ) = 0.d0;
            AC_el ( ndof0+IMDOF_el(i), 1:NDFT_el ) = 0.d0;
            AK_el ( ndof0+IMDOF_el(i), 1:NDFT_el ) = 0.d0;
      enddo


      do iq    = 1, QSnode(nn_elcon)%Tot_el(2)                         !Qs up to tower
         nq    = QSnode(nn_el)%IORD_el(iq)

         do i  = 1, 3
            AK_el ( ndof0+IMDOF_el(i), NDFBT_el+nq ) =  transf_mat(nn_el)%R0_el(iq,i) - transf_mat(nn_elcon)%R0_el(iq,i)   &
                                                       +dot_product(transf_mat(nn_el)%A0_el(iq,i,1:3), Y0_U(1:3))
         enddo
      enddo!iq 

      do iq    = QSnode(nn_elcon)%Tot_el(2)+1, QSnode(nn_el)%Tot_el(2) !Qs from tower to Blade
         nq    = QSnode(nn_el)%IORD_el(iq)

         do i  = 1, 3
            AK_el ( ndof0+IMDOF_el(i), NDFBT_el+nq ) =  transf_mat(nn_el)%R0_el(iq,i)                                      &
                                                       +dot_product(transf_mat(nn_el)%A0_el(iq,i,1:3), Y0_U(1:3))
         enddo
      enddo!iq 

      do i     = 1, 3
         do j  = 1, 3
            AK_el ( ndof0+IMDOF_el(i), ndof0+IMDOF_el(j) ) =  transf_mat(nn_el)%A_el(i,j)
         enddo
      enddo

         do i  = 1, 3
            AQ_el ( ndof0+IMDOF_el(i)             ) = -transf_mat(nn_el)%R_el (   i) + transf_mat(nn_elcon)%R_el (   i )  &
                                                      -dot_product(transf_mat(nn_el)%A_el(i,1:3), Y0_U(1:3))
         enddo
!stop
   enddo !nbod


 END Subroutine BOUNDCO_VAWT2_el
!----------------------------------------------------------------------
!
!  Subrouine : SBOD2SBOD_LOADS    ------------------
!
!  Loads Communication between subbodies of same body
!
!----------------------------------------------------------------------
 Subroutine SBOD2SBOD_LOADS
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   real(8) :: Aba      (3,3        )
   real(8) :: Aba0     (3,3,NQS    )
   real(8) :: Rba      (3          )
   real(8) :: Rba0     (3  ,NQS    )
   real(8) :: AMB2B    (6,NDFPEM_el)
   real(8) :: ACB2B    (6,NDFPEM_el)
   real(8) :: AKB2B    (6,NDFPEM_el)
   real(8) :: AMB2BNQ  (6,NQS      )
   real(8) :: ACB2BNQ  (6,NQS      )
   real(8) :: AKB2BNQ  (6,NQS      )
   real(8) :: AFB2B    (6          )
   integer :: nbod, i, j, NDFPE, iq, nq
   integer :: nnsubB, nbsubB, nnB_el, nbB_el, nelB, nodB
   integer :: nn    ,         nnA_el, nbA_el, nelA, nodA    
   integer :: IMAX  , JMAX  , IQRayl, IcalcR
   real(8) :: COEFM , COEFK


   do nbod    = 1, NBODBTWT_el
   do nbsubB  = 2, body(nbod)%NBODSUB_el
 
!------ subBody (B)
      nnsubB  = nbsubB
      nelB    = 1
      nodB    = 1
      nnB_el  = body(nbod)%NNODTGLB_el(nnsubB  )
      nbB_el  = body(nbod)%NBODTGLB_el(nbsubB  )
      NDFPE   = subbody(nbB_el)%NDFPE_el

!------ subBody (A)
      nnA_el  = body(nbod)%NNODTGLB_el(nnsubB-1)
      nbA_el  = body(nbod)%NBODTGLB_el(nbsubB-1)
      nelA    = subbody(nbA_el)%NTEB_el
      nodA    = NNPE_el
      nn      = subbody(nbA_el)%NODMB(nelA,nodA)
      i       = ACCU%NDFTACC_el(nbA_el-1)+subbody(nbA_el)%NDFPNBACC_el(nn-1     )!+1:6
      j       = ACCU%NDFTACC_el(nbB_el-1)!1:NDFPE !+subbody(nbB_el)%NDFPNBACC_el(mn-1)+mdof

!------ Aba, Aba0
      Aba(1:3,1:3) = matmul (transf_mat(nnA_el)%AT_el(1:3,1:3), transf_mat(nnB_el)%A_el(1:3,1:3))

      do iq   = 1, QSnode(nnA_el)%Tot_el(2)
         nq   = QSnode(nnB_el)%IORD_el(iq)

         Aba0( 1:3, 1:3, nq ) = matmul( transf_mat(nnA_el)%AT_el (   1:3,1:3), transf_mat(nnB_el)%A0_el(iq,1:3,1:3) ) + &
                                matmul( transf_mat(nnA_el)%AT0_el(iq,1:3,1:3), transf_mat(nnB_el)%A_el (   1:3,1:3) )
      enddo!iq 

      do iq   = QSnode(nnA_el)%Tot_el(2)+1, QSnode(nnB_el)%Tot_el(2)
         nq   = QSnode(nnB_el)%IORD_el(iq)

         Aba0( 1:3, 1:3, nq ) = matmul( transf_mat(nnA_el)%AT_el (   1:3,1:3), transf_mat(nnB_el)%A0_el(iq,1:3,1:3) )
      enddo!iq 

      IcalcR  = 0   !0: no R calculations, 1: R calculations
      IQRayl  = 1   !0: all qs Rayleigh  , 1: body qs with Rayleigh, 2: no Rayleigh to qs
      IMAX    = 6
      JMAX    = subbody(nbB_el)%NDFPE_el
      COEFM   = COEFM_el(nbod)
      COEFK   = COEFK_el(nbod)

      call LOC_LOADS_el       ( AMB2B  , ACB2B  , AKB2B  , AFB2B,&
                                AMB2BNQ, ACB2BNQ, AKB2BNQ       ,&
                                nbB_el  ,nnB_el  ,nelB   , nodB    )


      call LOADS_TRANS_el     ( AMB2B  , ACB2B  , AKB2B  , AFB2B,&
                                AMB2BNQ, ACB2BNQ, AKB2BNQ       ,&
                                Aba0   , Aba    , Rba0   , Rba  ,&
                                         nnB_el , IcalcR        ,&
                                JMAX   , NDFPEM_el                 )


      call ASSEMBLY_MATRIX_el ( AMB2B  , ACB2B  , AKB2B  , AFB2B         ,&
                                AMB2BNQ, ACB2BNQ, AKB2BNQ, COEFM , COEFK ,&
                                         nnB_el , i      , j     , IQRayl,&
                                IMAX   , JMAX   , 6      , NDFPEM_el       )
   enddo !nbsub
   enddo !nbod


 END Subroutine SBOD2SBOD_LOADS
!----------------------------------------------------------------------
!
!  Subrouine : BOD2BOD_LOADS  ----------------
!
!  Loads Communication between WT's bodies (1st Subbody with last one of the previous body)
!
!----------------------------------------------------------------------
 Subroutine BOD2BOD_LOADS ( nbod1, nbod2 )
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   real(8) :: Aba      (3,3        )
   real(8) :: Aba0     (3,3,NQS    )
   real(8) :: Rba      (3          )
   real(8) :: Rba0     (3  ,NQS    )
   real(8) :: RGba     (3          )
   real(8) :: AMB2B    (6,NDFPEM_el)
   real(8) :: ACB2B    (6,NDFPEM_el)
   real(8) :: AKB2B    (6,NDFPEM_el)
   real(8) :: AMB2BNQ  (6,NQS      )
   real(8) :: ACB2BNQ  (6,NQS      )
   real(8) :: AKB2BNQ  (6,NQS      )
   real(8) :: AFB2B    (6          )
   integer :: nbod1, nbod2, i, j, NDFPE, iq, nq
   integer :: nbodB, nnsubB, nbsubB, nnB_el,          nbB_el, nelB, nodB
   integer :: nbodA, nn    ,         nnA_el, nnA2_el, nbA_el, nelA, nodA
   integer :: IMAX , JMAX  , IQRayl, IcalcR
   real(8) :: COEFM, COEFK


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
      JMAX    = subbody(nbB_el)%NDFPE_el
      COEFM   = 0.d0
      COEFK   = 0.d0

      call LOC_LOADS_el       ( AMB2B  , ACB2B  , AKB2B  , AFB2B,&
                                AMB2BNQ, ACB2BNQ, AKB2BNQ       ,&
                                nbB_el  ,nnB_el  ,nelB   , nodB    )


      call LOADS_TRANS_el     ( AMB2B  , ACB2B  , AKB2B  , AFB2B,&
                                AMB2BNQ, ACB2BNQ, AKB2BNQ       ,&
                                Aba0   , Aba    , Rba0   , Rba  ,&
                                         nnB_el , IcalcR        ,&
                                JMAX   , NDFPEM_el                 )


      call ASSEMBLY_MATRIX_el ( AMB2B  , ACB2B  , AKB2B  , AFB2B         ,&
                                AMB2BNQ, ACB2BNQ, AKB2BNQ, COEFM , COEFK ,&
                                         nnB_el , i      , j     , IQRayl,&
                                IMAX   , JMAX   , 6      , NDFPEM_el       )
   enddo !nbod


 END Subroutine BOD2BOD_LOADS
!----------------------------------------------------------------------
!
!  Subrouine : BOD2BOD_LOADS  ----------------
!
!  Loads Communication between WT's bodies (1st Subbody with last one of the previous body)
!
!----------------------------------------------------------------------
 Subroutine BOD2BOD_LOADS_VAWT
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   real(8) :: Aba      (3,3        )
   real(8) :: Aba0     (3,3,NQS    )
   real(8) :: Rba      (3          )
   real(8) :: Rba0     (3  ,NQS    )
   real(8) :: RGba     (3          )
   real(8) :: AMB2B    (6,NDFPEM_el)
   real(8) :: ACB2B    (6,NDFPEM_el)
   real(8) :: AKB2B    (6,NDFPEM_el)
   real(8) :: AMB2BNQ  (6,NQS      )
   real(8) :: ACB2BNQ  (6,NQS      )
   real(8) :: AKB2BNQ  (6,NQS      )
   real(8) :: AFB2B    (6          )
   integer :: i, j, NDFPE, iq, nq
   integer :: nbodB,         nnB_el, nnB2_el, nbB_el, nelB, nodB
   integer :: nbodA, nn    , nnA_el, nnA2_el, nbA_el, nelA, nodA
   integer :: IMAX , JMAX  , IQRayl, IcalcR
   real(8) :: COEFM, COEFK


   if (IAPPL_el /= 1) return!0:HAWT, 1:VAWT
   if (NBODBTWT_el /= NTOWER_el) return !only blades and shaft, so no loading transfer


   do nbodB   = 1, NBLADE_el
!------ Body (B)
      nnB_el  = body(nbodB)%NNODTGLB_el (body(nbodB)%NBODSUB_el   )
      nnB2_el = body(nbodB)%NNODTGLB_el (body(nbodB)%NBODSUB_el+1 )
      nbB_el  = body(nbodB)%NBODTGLB_el (body(nbodB)%NBODSUB_el   )
      nelB    = subbody(nbB_el)%NTEB_el
      nodB    = NNPE_el
      NDFPE   = subbody(nbB_el)%NDFPE_el

!------ Body (A)
      nbodA   = NTOWER_el
      nnA_el  = body(nbodA)%NNODTGLB_el (body(nbodA)%NBODSUB_el   )
      nnA2_el = body(nbodA)%NNODTGLB_el (body(nbodA)%NBODSUB_el+1 )
      nbA_el  = body(nbodA)%NBODTGLB_el (body(nbodA)%NBODSUB_el   )
      nelA    = subbody(nbA_el)%NTEB_el
      nodA    = NNPE_el
      nn      = subbody(nbA_el)%NODMB(nelA,nodA)
      i       = ACCU%NDFTACC_el(nbA_el-1)+subbody(nbA_el)%NDFPNBACC_el(nn-1     )!+1:6
      j       = ACCU%NDFTACC_el(nbB_el-1)!1:NDFPE !+subbody(nbB_el)%NDFPNBACC_el(nnB-1)+mdof

!------ Aba, Aba0, Rba, Rba0
      RGba(1:3    ) = transf_mat(nnB2_el)%R_el(1:3) - transf_mat(nnA2_el)%R_el(1:3)
      Aba (1:3,1:3) = matmul ( transf_mat(nnA_el)%AT_el(1:3,1:3), transf_mat(nnB_el)%A_el(1:3,1:3) )
      Rba (1:3    ) = matmul ( transf_mat(nnA_el)%AT_el(1:3,1:3), RGba(1:3)                        )
      Aba0(:,:,:  ) = 0.d0;
      Rba0(:,  :  ) = 0.d0;

      do iq   = 1, QSnode(nnA_el)%Tot_el(2)
         nq   = QSnode(nnA_el)%IORD_el(iq)
         Aba0(1:3,1:3,nq) = matmul(  transf_mat(nnA_el )%AT_el (   1:3,1:3),  transf_mat(nnB_el )%A0_el(iq,1:3,1:3) ) + &
                            matmul(  transf_mat(nnA_el )%AT0_el(iq,1:3,1:3),  transf_mat(nnB_el )%A_el (   1:3,1:3) )

         Rba0(1:3,    nq) = matmul( (transf_mat(nnB2_el)%R0_el (iq,1:3    ) - transf_mat(nnA2_el)%R0_el(iq,1:3    ) ), transf_mat(nnA_el)%AT_el (   1:3,1:3) ) + &
                            matmul(  RGba                      (   1:3    )                                          , transf_mat(nnA_el)%AT0_el(iq,1:3,1:3) )
      enddo!iq 

      do iq   = QSnode(nnA_el)%Tot_el(2)+1, QSnode(nnA2_el)%Tot_el(2)
         nq   = QSnode(nnA2_el)%IORD_el(iq)
         Aba0(1:3,1:3,nq) = matmul(  transf_mat(nnA_el)%AT_el  (   1:3,1:3),  transf_mat(nnB_el )%A0_el(iq,1:3,1:3) )

         Rba0(1:3,    nq) = matmul( (transf_mat(nnB2_el)%R0_el (iq,1:3    ) - transf_mat(nnA2_el)%R0_el(iq,1:3    ) ), transf_mat(nnA_el)%AT_el (   1:3,1:3) )
      enddo!iq 

      do iq   = QSnode(nnA2_el)%Tot_el(2)+1, QSnode(nnB_el)%Tot_el(2)
         nq   = QSnode(nnB_el)%IORD_el(iq)
         Aba0(1:3,1:3,nq) = matmul(  transf_mat(nnA_el )%AT_el (   1:3,1:3),  transf_mat(nnB_el )%A0_el(iq,1:3,1:3) )

         Rba0(1:3,    nq) = matmul( (transf_mat(nnB2_el)%R0_el (iq,1:3    )                                         ), transf_mat(nnA_el)%AT_el (   1:3,1:3) )
      enddo!iq 

      do iq   = QSnode(nnB_el)%Tot_el(2)+1, QSnode(nnB2_el)%Tot_el(2)
         nq   = QSnode(nnB2_el)%IORD_el(iq)
         Aba0(1:3,1:3,nq) = 0.d0

         Rba0(1:3,    nq) = matmul( (transf_mat(nnB2_el)%R0_el (iq,1:3    )                                         ), transf_mat(nnA_el)%AT_el (   1:3,1:3) )
      enddo!iq 

      IcalcR  = 1   !0: no R calculations, 1: R calculations
      IQRayl  = 2   !0: all qs Rayleigh  , 1: body qs with Rayleigh, 2: no Rayleigh to qs
      IMAX    = 6
      JMAX    = subbody(nbB_el)%NDFPE_el
      COEFM   = 0.d0
      COEFK   = 0.d0

      call LOC_LOADS_el       ( AMB2B  , ACB2B  , AKB2B  , AFB2B,&
                                AMB2BNQ, ACB2BNQ, AKB2BNQ       ,&
                                nbB_el  ,nnB_el  ,nelB   , nodB    )


      call LOADS_TRANS_el     ( AMB2B  , ACB2B  , AKB2B  , AFB2B,&
                                AMB2BNQ, ACB2BNQ, AKB2BNQ       ,&
                                Aba0   , Aba    , Rba0   , Rba  ,&
                                         nnB2_el, IcalcR        ,&
                                JMAX   , NDFPEM_el                 )


      call ASSEMBLY_MATRIX_el ( AMB2B  , ACB2B  , AKB2B  , AFB2B         ,&
                                AMB2BNQ, ACB2BNQ, AKB2BNQ, COEFM , COEFK ,&
                                         nnB2_el, i      , j     , IQRayl,&
                                IMAX   , JMAX   , 6      , NDFPEM_el       )
   enddo !nbod


 END Subroutine BOD2BOD_LOADS_VAWT
!----------------------------------------------------------------------
!
!  Subrouine : LOADS_TRANS_el  ----------------
!
!  Transfer the Local loads to another (local) c.o. system
!
!----------------------------------------------------------------------
 Subroutine LOADS_TRANS_el ( AMB2B  , ACB2B  , AKB2B  , AFB2B,&
                             AMB2BNQ, ACB2BNQ, AKB2BNQ       ,&
                             Aba0   , Aba    , Rba0   , Rba  ,&
                                      nn_el  ,IcalcR         ,&
                             JMAX   , Jdi                       )
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   integer :: JMAX, Jdi
   real(8) :: Aba      (3,3    )
   real(8) :: Aba0     (3,3,NQS)
   real(8) :: Rba      (3      )
   real(8) :: Rba0     (3  ,NQS)
   real(8) :: Aba0q    (3,3    )
   real(8) :: Rba0q    (3      )
   real(8) :: AMB2B    (6,Jdi  ) !(6,NDFPEM_el)
   real(8) :: ACB2B    (6,Jdi  ) !(6,NDFPEM_el)
   real(8) :: AKB2B    (6,Jdi  ) !(6,NDFPEM_el)
   real(8) :: AMB2BNQ  (6,NQS  )
   real(8) :: ACB2BNQ  (6,NQS  )
   real(8) :: AKB2BNQ  (6,NQS  )
   real(8) :: AFB2B    (6      )
   real(8) :: FLOC     (6      )
   real(8) :: AbaQ     (6      )
   real(8) :: Q        (6      )
   integer :: nn_el, IcalcR, iq, nq


!--- Move LOADS to another C. S. (A) using Aba,Aba0, Rba,Rba0
   Q          = AFB2B;

   call LOADS_TRANS0_NEW ( AMB2B, Aba, 1, JMAX, Jdi )! Aba.M
   call LOADS_TRANS0_NEW ( ACB2B, Aba, 1, JMAX, Jdi )! Aba.C
   call LOADS_TRANS0_NEW ( AKB2B, Aba, 1, JMAX, Jdi )! Aba.K
   call LOADS_TRANS0_NEW ( AFB2B, Aba, 1, 1   , 1   )! Aba.Q
                                                          

   if (IcalcR == 1) then

      AbaQ    = AFB2B; ! Aba.Q

      call LOADS_CALC_M0    ( AMB2B, Rba, 1, JMAX, Jdi )! Rba x Aba.M
      call LOADS_CALC_M0    ( ACB2B, Rba, 1, JMAX, Jdi )! Rba x Aba.C
      call LOADS_CALC_M0    ( AKB2B, Rba, 1, JMAX, Jdi )! Rba x Aba.K
      call LOADS_CALC_M0    ( AFB2B, Rba, 1, 1   , 1   )! Rba x Aba.Q


!------ For all q's
      do iq   = 1, QSnode(nn_el)%Tot_el(2)
         nq   = QSnode(nn_el)%IORD_el(iq)
 
         call LOADS_TRANS0_NEW ( AMB2BNQ, Aba, nq, nq, NQS )! Aba.Mq
         call LOADS_TRANS0_NEW ( ACB2BNQ, Aba, nq, nq, NQS )! Aba.Cq
         call LOADS_TRANS0_NEW ( AKB2BNQ, Aba, nq, nq, NQS )! Aba.Kq

         call LOADS_CALC_M0    ( AMB2BNQ, Rba, nq, nq, NQS )! Rba x Aba.Mq
         call LOADS_CALC_M0    ( ACB2BNQ, Rba, nq, nq, NQS )! Rba x Aba.Cq
         call LOADS_CALC_M0    ( AKB2BNQ, Rba, nq, nq, NQS )! Rba x Aba.Kq

         FLOC = Q;

         Aba0q(1:3,1:3) = Aba0(1:3,1:3,nq)
         Rba0q(1:3    ) = Rba0(1:3,    nq)

         call LOADS_TRANS0_NEW ( FLOC, Aba0q, 1, 1, 1 )! Aba0i.Q
         call LOADS_CALC_M0    ( FLOC, Rba  , 1, 1, 1 )! Rba   x Aba0i.Q
!check   call LOADS_CALC_M0    ( FLOC, Rba0q, 1, 1, 1 )! Rba0i x Aba0i.Q

         AKB2BNQ(1:6,nq) = AKB2BNQ(1:6,nq) - FLOC(1:6)


         FLOC = AbaQ;

!--------- no add to FLOC
         call LOADS_CALC_M00   ( FLOC, Rba0q , 1, 1, 1 )! Rba0i x Aba.Q

         AKB2BNQ(2  ,nq) = AKB2BNQ(2  ,nq) - FLOC(2  )
         AKB2BNQ(5:6,nq) = AKB2BNQ(5:6,nq) - FLOC(5:6)
      enddo!iq 

   else !IcalcR==0

!------ For all q's
      do iq   = 1, QSnode(nn_el)%Tot_el(2)
         nq   = QSnode(nn_el)%IORD_el(iq)

         call LOADS_TRANS0_NEW ( AMB2BNQ, Aba, nq, nq, NQS )! Aba.Mq
         call LOADS_TRANS0_NEW ( ACB2BNQ, Aba, nq, nq, NQS )! Aba.Cq
         call LOADS_TRANS0_NEW ( AKB2BNQ, Aba, nq, nq, NQS )! Aba.Kq


         FLOC = Q;

         Aba0q(1:3,1:3) = Aba0(1:3,1:3,nq)

         call LOADS_TRANS0_NEW ( FLOC, Aba0q, 1, 1, 1 )! Aba0i.Q

         AKB2BNQ(1:6,nq) = AKB2BNQ(1:6,nq) - FLOC(1:6)
      enddo!iq 

   endif


 END Subroutine LOADS_TRANS_el
!----------------------------------------------------------------------
!
!  Subrouine : LOC_LOADS_el  ----------------
!
!  Set the Local loads of nod 1 or NNPE_el

!  IQRayl = 0: all qs Rayleigh  , 1: body qs with Rayleigh, 2: no Rayleigh to qs
!
!----------------------------------------------------------------------
 Subroutine LOC_LOADS_el ( AMB2B  , ACB2B  , AKB2B  , AFB2B ,&
                           AMB2BNQ, ACB2BNQ, AKB2BNQ        ,&
                           nb_el  , nn_el  ,nel    , nod       )
!----------------------------------------------------------------------

 Use Cbeam

   implicit none
  
   real(8) :: AMB2B    (6,NDFPEM_el)
   real(8) :: ACB2B    (6,NDFPEM_el)
   real(8) :: AKB2B    (6,NDFPEM_el)
   real(8) :: AMB2BNQ  (6,NQS      )
   real(8) :: ACB2BNQ  (6,NQS      )
   real(8) :: AKB2BNQ  (6,NQS      )
   real(8) :: AFB2B    (6          )
   integer :: nb_el, nn_el, nel, nod, NDFPE, iq, nq, j, mn, I1, I2


!--- Set the loads of Body
   I1     = subbody(nb_el)%NDFPNACC_el  (nod-1   ) + 1
   I2     = subbody(nb_el)%NDFPNACC_el  (nod     )
   mn     = subbody(nb_el)%NODMB        (nel  , 1)      !for all element nodes in j direction
   j      =           ACCU%NDFTACC_el   (nb_el-1 ) + &
            subbody(nb_el)%NDFPNBACC_el (mn-1    )
   NDFPE  = subbody(nb_el)%NDFPE_el


   AMB2B  ( 1:6, 1:NDFPE ) = subbody(nb_el)%AMLOC_el ( nel, I1:I2, 1:NDFPE )
   ACB2B  ( 1:6, 1:NDFPE ) = subbody(nb_el)%ACLOC_el ( nel, I1:I2, 1:NDFPE )
   AKB2B  ( 1:6, 1:NDFPE ) = subbody(nb_el)%AKLOC_el ( nel, I1:I2, 1:NDFPE )
   AFB2B  ( 1:6          ) = subbody(nb_el)%AFLOC_el ( nel, I1:I2          )

   AMB2BNQ (:,:) = 0.d0; ACB2BNQ (:,:) = 0.d0; AKB2BNQ (:,:) = 0.d0;

!--- For all q's
   do iq = 1, QSnode(nn_el)%Tot_el(2)
      nq = QSnode(nn_el)%IORD_el(iq)

      AMB2BNQ ( 1:6, nq ) = subbody(nb_el)%AMLOCNQ_el ( nel, I1:I2, iq )
      ACB2BNQ ( 1:6, nq ) = subbody(nb_el)%ACLOCNQ_el ( nel, I1:I2, iq )
      AKB2BNQ ( 1:6, nq ) = subbody(nb_el)%AKLOCNQ_el ( nel, I1:I2, iq )
   enddo


 END Subroutine LOC_LOADS_el
!--------------------------------------------------------------------------------
!
!  Subroutine : QS_EQUAT   ------------------
!
!  Equations for the q d.o.f
!
!-- Jimmys check the zero of the qs Q=Qo+dq and when ZER0_DOF is necesary.
!----------------------------------------------------------------------
 Subroutine QS_EQUAT (IP)
!----------------------------------------------------------------------

 Use Cbeam
 Use Hydro

   implicit none

   integer :: jj, j, ii, i, NQ_el, nbod, nnsub
   integer :: nel, mnod, mn, IQ, ju, ju1, IP  , nbpre_el, nq, nf


   if (NQS  == 0) return !case of jacket alone

   if (IAPPL_el == 2) then
!------ Zero the disabled heli qs
         nf = 1
      do iq = 1, IQfl

         if ( floater(nf)% IQfl_active_el(iq) == 1 ) cycle

         i           = NDFBT_el + (nf-1)*IQfl + iq
         INDSYSQ(iq) = -1

         call ZERO_DOF ( IP,  i, 0, 0 )
      enddo!iq
      goto 1
   endif

!--- NQ = 1:IQfl
!--- Floater's 6 dof equations
   if (ICASE_el < 3) goto 1

!----------------- FLOATER's EQUATIONS ------------------------------

!--- Communicate Forces from Wind Turbine     to 6 Qs
   call FLOATLOADS_WT26Qs
!--- Communicate Forces from flexible Floater to 6 Qs
   call FLOATLOADS_JA26Qs


! WT transfers loads to 1st floater by default atm
! 2nd floater is freely floating
! the only coupling between the 2 floaters is via radiation added damping coefficients


!--- Hydrodynamic M, C, K, Q and moorings Q_truss
   do nf = 1, NFLOATER_el

!------ Assembly floater Matrices to global matrices
      do    ii = 1, IQfl
            i  = NDFBT_el + (nf-1)*IQfl + ii
         do jj = 1, IQfl
            j  = NDFBT_el + (nf-1)*IQfl + jj

            AC_el(i,j) = AC_el(i,j) + floater(nf)% AC_morison (ii,jj)
            AK_el(i,j) = AK_el(i,j) + floater(nf)% AK_total   (ii,jj)
!                                   + floater(nf)% AK_morison (ii,jj)
            AQ_el(i  ) = AQ_el(i  ) - floater(nf)% AK_total   (ii,jj)*UT_el (j)
         enddo !jj

         do jj = 1, IQfl*NFLOATER_el
            j  = NDFBT_el + jj

            AM_el(i,j) = AM_el(i,j) + floater(nf)% AM_total   (ii,jj)
            AC_el(i,j) = AC_el(i,j) + floater(nf)% AC_total   (ii,jj)
            AQ_el(i  ) = AQ_el(i  ) - floater(nf)% AM_total   (ii,jj)*UT2_el(j) &
                                    - floater(nf)% AC_total   (ii,jj)*UT1_el(j)
         enddo !jj

            AQ_el(i  ) = AQ_el(i  ) + floater(nf)% Q_total    (ii   ) &
                                    + floater(nf)% Q_truss    (ii   ) &
                                    + floater(nf)% AQ_morison (ii   )
         INDSYSQ (ii ) = 0
      enddo !ii


!------ Zero the disabled floater's qs
      do iq = 1, IQfl

         if ( floater(nf)% IQfl_active_el(iq) == 1 ) cycle

         i           = NDFBT_el + (nf-1)*IQfl + iq
         INDSYSQ(iq) = -1

         call ZERO_DOF ( IP,  i, 0, 0 )
      enddo!iq

   enddo !nf


!--- NQ = IQfl+1
!--- rotor azimuth rotation equation
 1 nq = NQSW
   iq = NDFBT_el + nq

   if ((IEIG > 0).or.(IMODAL > 0)) then
!------ Eigen-value calculations ---------------------------------
      if ( (IEIG == 1).and.(IMODAL<=0).and.(ICMOD >= 1) )  then
!--------- drive train: free-free
         INDSYSQ(nq) = 0
         call GENLOAD
      else
!--------- drive train: free-fixed
         INDSYSQ(nq) = -1
         call ZERO_DOF ( IP, iq, 0, 0 )
      endif
!--- time domain or static analysis ------------------------------
   elseif (ICMOD == 0) then
!------ drive train: free-fixed
      INDSYSQ(nq) = 0

      call ZERO_DOF ( IP, iq, 0, 0 )
   elseif (ICMOD == 5) then
!------ drive train: free-fixed, set prescribed values
      INDSYSQ(nq) =-1

      call ZERO_DOF ( IP, iq, 0, 0 )
   else !ICMOD 1,2,3,4,[10-->1]
!------ drive train: free-free
      INDSYSQ(nq) = 0

      call GENLOAD
   endif


!--- NQ = IQfl+2
!--- yaw actuator equation
   nq = NQSY
   iq = NDFBT_el + nq

   if (IACT_YAW == 0) then
      INDSYSQ(nq) = -1

      call ZERO_DOF ( IP, iq, 0, 0 )
   else
      INDSYSQ(nq) = 0

      call ACTYAWEQ
   endif


!--- NQ = NQSP+1:NQS
!--- all WT SubBodies q's
   do nbod = 1, NBODBTWT_el

      if ( (ICASE_el == 2).and.(nbod == NBODBTWT_el) ) then
!--------- Communicate Jacket's qs to WT. By default communicate the last body [NBODBT_el]
         call QS_TOWER_WT_JA (IP,NBODBT_el)
      else
!--------- Fixed conditions at the beginning of the body (nnsub=1)
        !nnsub = 1
         nq    = NQSP + ACCU%NQSACC_el (nbod-1) !+ ACCU%NQSBACC_el(nbod,nnsub-1)

         INDSYSQ(nq+1:nq+6) = -1

         iq    = nq   + NDFBT_el
         call ZERO_DOF ( IP, iq, 1, 6 )
      endif

!------ Assign the first 6 qs to the last 6 elastic dof
      do    nnsub    = 2, body(nbod)%NBODSUB_el+1
            nbpre_el = body(nbod)%NBODTGLB_el(nnsub-1)
            nel      = subbody(nbpre_el)%NTEB_el
            mnod     = NNPE_el
            mn       =                               subbody(nbpre_el)%NODMB(nel,mnod)
               ju    = ACCU%NDFTACC_el(nbpre_el-1) + subbody(nbpre_el)%NDFPNBACC_el(mn-1)
            NQ_el    = NQSP + ACCU%NQSACC_el (nbod-1) + ACCU%NQSBACC_el(nbod, nnsub-1)
         do IQ       = 1, 6
            NQ_el    = NQ_el + 1
            i        = NDFBT_el + NQ_el
            ju1      = ju + IMDOF_el(IQ)

            INDSYSQ(NQ_el)  = ju1

            if (IP == 0) then
               AK_el(i,i  ) = 1.d0
               AK_el(i,ju1) =-1.d0
               AQ_el(i    ) = UT_el (ju1)-UT_el (i)
            else
               AM_el(i,i  ) = 1.d0
               AM_el(i,ju1) =-1.d0
               AQ_el(i    ) = UT2_el(ju1)-UT2_el(i)
            endif
         enddo !IQ

      enddo !nnsub
   enddo !nbod


 END Subroutine QS_EQUAT
!----------------------------------------------------------------------
!
!  Subroutine : MATRIX_REDUCT   ------------------
!
!  Reduce the Global Matrices for:
!
!       1. fixed conditions (bodies or q's)
!       2. dependent qs
!       3. dependent deflections (common node (jacket case))
!
!----------------------------------------------------------------------
 Subroutine MATRIX_REDUCT (IP)
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   integer :: INDQ(NQS), nq, i, j, ju, jq, ii, jj, nb_el, ndof, NQS0
   integer ::            n_rem, nbc, IP, indx1(6),indx2(6)
   integer :: nbod,nbsub,iel,inodel,inode,iq,iq0
   integer :: iecon,incon,inodc,ibcon,isbcon,ibsubcon
 

   NDFT0_el    = NDFT_el
   NQS0        = NQS
   NDFBT0_el   = NDFBT_el
   NDFBTWT0_el = NDFBTWT_el

!----------------------
!-- No reduction
!----------------------
   if (IREDUCT_el == 0) then

!------ Boundary Conditions
      call BOUNDCO_el(IP)
      call BOUNDCOJack_el(IP)
      call BOUNDCO_VAWT_el(IP)

      INDSYSQ (1:NQS) = 0
      do i = 1, NDFT0_el
         INDSYSB(i)   = i
      enddo

      return
   endif
      call BOUNDCO_VAWT_el(IP) !atm no reduction for the blade tip for the VAWT

!----------------------------------------------------------
!-- WT sub-bodies: Add dependent Qs to their contribution
!----------------------------------------------------------
      NQS0 = 0
   do nq   = 1, NQS
      ju   = INDSYSQ (nq)

      if (ju == 0) then
!--------- Kept Qs [not reduced]
         NQS0 = NQS0 + 1
         INDQ(NQS0) = nq
         cycle
      elseif (ju < 0) then
!--------- Disabled Qs
         cycle
      endif
!--------- Reduced Qs
      jq   = NDFBT_el + nq
      do i = 1, NDFT_el
         AM_el(i,ju) = AM_el(i,ju) + AM_el(i,jq)
         AC_el(i,ju) = AC_el(i,ju) + AC_el(i,jq)
         AK_el(i,ju) = AK_el(i,ju) + AK_el(i,jq)
      enddo !i
   enddo !nq

!----------------------------------------------------------
!-- JACKET: Add dependent dofs, to dofs which will be kept
!----------------------------------------------------------
   do    nbod     = 1+NBODBTWT_el, NBODBT_el
   do    nbsub    = 1, body(nbod)%NBODSUB_el
         nb_el    =    body(nbod)%NBODTGLB_el(nbsub)
      do nbc      = 1, boundco (nb_el)%nbcpb
         ibcon    =    boundco (nb_el)%nbcon(nbc)
         if (ibcon == 0) cycle

         iel      = boundco (nb_el)%nel (nbc)
         inodel   = boundco (nb_el)%nod (nbc)
         inode    =                            subbody(nb_el)%NODMB       (iel,inodel)
         iq0      = ACCU%NDFTACC_el(nb_el-1) + subbody(nb_el)%NDFPNBACC_el(inode-1   )

         ibsubcon = 1
         isbcon   = body    (ibcon)%NBODTGLB_el(ibsubcon)
         incon    = boundco (nb_el)%nodcon(nbc)
         iecon    = boundco (nb_el)%nelcon(nbc)
         inodc    =                               subbody(isbcon)%NODMB       (iecon,incon)
         ju       = ACCU%NDFTACC_el  (isbcon-1) + subbody(isbcon)%NDFPNBACC_el(inodc-1    )

         do ndof  = 1, 6
            if (boundco (nb_el)%indx(nbc,ndof) == 0) cycle
            do iq = iq0+1, iq0+6
               AM_el ( 1:NDFT_el, ju+ndof ) = AM_el ( 1:NDFT_el, ju+ndof ) + AM_el ( 1:NDFT_el, iq ) * boundco (nb_el)%Atrans1 (nbc,iq-iq0,ndof)
               AC_el ( 1:NDFT_el, ju+ndof ) = AC_el ( 1:NDFT_el, ju+ndof ) + AC_el ( 1:NDFT_el, iq ) * boundco (nb_el)%Atrans1 (nbc,iq-iq0,ndof)
               AK_el ( 1:NDFT_el, ju+ndof ) = AK_el ( 1:NDFT_el, ju+ndof ) + AK_el ( 1:NDFT_el, iq ) * boundco (nb_el)%Atrans1 (nbc,iq-iq0,ndof)
            enddo
         enddo

      enddo !nbc

   enddo !nbsub
   enddo !nbod

!-------------------------------------------------
!- Reduce matrix for boundary conditions
!-------------------------------------------------
!
!-   1. WT zero condition (take care of foundation)
!-   2. JA dependent dofs
!
!----------------
!- Keep LINES
!----------------

!--- Blades and shaft zero dofs
         ii      = 0
   do    nbod    = 1, min( NSHAFT_el, NBODBTWT_el ) !!NBODBTWT_el - 1
      do nbsub   = 1, body(nbod)%NBODSUB_el
         nb_el   =    body(nbod)%NBODTGLB_el(nbsub)

!--------- For all non-zero dofs
         do ndof = 7, subbody(nb_el)%NDFTB_el
            ii   = ii + 1
            i    = ACCU%NDFTACC_el(nb_el-1) + ndof

            INDSYSB(ii) = i

            AM_el ( ii, 1:NDFT_el ) = AM_el ( i, 1:NDFT_el )
            AC_el ( ii, 1:NDFT_el ) = AC_el ( i, 1:NDFT_el )
            AK_el ( ii, 1:NDFT_el ) = AK_el ( i, 1:NDFT_el )
            AQ_el ( ii            ) = AQ_el ( i            )
         enddo !ndof
         NDFBT0_el   = NDFBT0_el   - 6
         NDFBTWT0_el = NDFBTWT0_el - 6
      enddo !nbsub
   enddo !nbod


!--- Tower Foundation [only for tower's 6dofs of 1st subbody:
   if (NBODBTWT_el==NTOWER_el) then

         nbod  = NBODBTWT_el
         nbsub = 1
         nb_el = body(nbod)%NBODTGLB_el(nbsub)

!------ For all non-zero dofs
      if (ITYPE_FOUND_el > 0) then
            n_rem = 0
            nbc   = 1
         do ndof  = 1, 6
            if ((boundco (nb_el)%nbcon(nbc     ) == 0).and. &
                (boundco (nb_el)%indx (nbc,ndof) == 1)        ) then

               n_rem = n_rem + 1
               cycle
            endif

            ii    = ii + 1
            i     = ACCU%NDFTACC_el(nb_el-1) + ndof

            INDSYSB(ii) = i

            AM_el ( ii, 1:NDFT_el ) = AM_el ( i, 1:NDFT_el )
            AC_el ( ii, 1:NDFT_el ) = AC_el ( i, 1:NDFT_el )
            AK_el ( ii, 1:NDFT_el ) = AK_el ( i, 1:NDFT_el )
            AQ_el ( ii            ) = AQ_el ( i            )
         enddo !ndof
      else !if (ITYPE_FOUND_el == 0) then
         n_rem = 6
      endif !foundation

      do ndof = 7, subbody(nb_el)%NDFTB_el
         ii   = ii + 1
         i    = ACCU%NDFTACC_el(nb_el-1) + ndof

         INDSYSB(ii) = i

         AM_el ( ii, 1:NDFT_el ) = AM_el ( i, 1:NDFT_el )
         AC_el ( ii, 1:NDFT_el ) = AC_el ( i, 1:NDFT_el )
         AK_el ( ii, 1:NDFT_el ) = AK_el ( i, 1:NDFT_el )
         AQ_el ( ii            ) = AQ_el ( i            )
      enddo !ndof
      NDFBT0_el   = NDFBT0_el   - n_rem
      NDFBTWT0_el = NDFBTWT0_el - n_rem
   endif !NBODBTWT_el=NTOWER_el


!--- Jacket dependent dofs
   do nb_el=1+NBODTWT_el,NBODT_el

!------ default values if no bc is applied
      n_rem      = 0
      indx1(1:6) = 0
      indx2(1:6) = 0

      do nbc = 1, boundco (nb_el)%nbcpb
         if      (boundco (nb_el)%nod (nbc) == 1) then
            do ndof = 1, 6
               if (boundco (nb_el)%indx(nbc,ndof) == 0) cycle
               n_rem = n_rem+1
               indx1(ndof) = 1
            enddo
         else !if(boundco (nb_el)%nod (nbc) == NNPE_el) then
            do ndof = 1, 6
               if (boundco (nb_el)%indx(nbc,ndof) == 0) cycle
               n_rem = n_rem+1
               indx2(ndof) = 1
            enddo
         endif
      enddo !nbc

      do ndof = 1, 6
         if (indx1(ndof) == 1) cycle  !keep the non dependent dofs
         ii   = ii + 1
         i    = ACCU%NDFTACC_el(nb_el-1) + ndof
         
         INDSYSB(ii) = i
         
         AM_el ( ii, 1:NDFT_el ) = AM_el ( i, 1:NDFT_el )
         AC_el ( ii, 1:NDFT_el ) = AC_el ( i, 1:NDFT_el )
         AK_el ( ii, 1:NDFT_el ) = AK_el ( i, 1:NDFT_el )
         AQ_el ( ii            ) = AQ_el ( i            )
      enddo
      do ndof =  7, subbody(nb_el)%NDFTB_el - 6
         ii   = ii + 1
         i    = ACCU%NDFTACC_el(nb_el-1) + ndof
         
         INDSYSB(ii) = i
         
         AM_el ( ii, 1:NDFT_el ) = AM_el ( i, 1:NDFT_el )
         AC_el ( ii, 1:NDFT_el ) = AC_el ( i, 1:NDFT_el )
         AK_el ( ii, 1:NDFT_el ) = AK_el ( i, 1:NDFT_el )
         AQ_el ( ii            ) = AQ_el ( i            )
      enddo
      do ndof =  subbody(nb_el)%NDFTB_el - 5, subbody(nb_el)%NDFTB_el
         i    = ndof-subbody(nb_el)%NDFTB_el + 6
         if (indx2(i)    == 1) cycle  !keep the non dependent dofs
         ii   = ii + 1
         i    = ACCU%NDFTACC_el(nb_el-1) + ndof
         
         INDSYSB(ii) = i
         
         AM_el ( ii, 1:NDFT_el ) = AM_el ( i, 1:NDFT_el )
         AC_el ( ii, 1:NDFT_el ) = AC_el ( i, 1:NDFT_el )
         AK_el ( ii, 1:NDFT_el ) = AK_el ( i, 1:NDFT_el )
         AQ_el ( ii            ) = AQ_el ( i            )
      enddo

      NDFBT0_el = NDFBT0_el - n_rem
   enddo !nb_el


!--- TRUSS COUPLED ---------
   call MOORINGS_COUPLED_tr ( 3, ii )

   NDFBT0_el = ii


!--- Dependent Q's----------
   do nq = 1, NQS0
      ii = ii + 1
      i  = NDFBT_el   + INDQ(nq)

      INDSYSB(ii) = i
 
      AM_el ( ii, 1:NDFT_el ) = AM_el ( i, 1:NDFT_el )
      AC_el ( ii, 1:NDFT_el ) = AC_el ( i, 1:NDFT_el )
      AK_el ( ii, 1:NDFT_el ) = AK_el ( i, 1:NDFT_el )
      AQ_el ( ii            ) = AQ_el ( i            )
   enddo !nq

   NDFT0_el  = ii

!------------------------
!- Keep COLUMNS
!------------------------

!--- Blades and shaft zero dofs
         jj      = 0
   do    nbod    = 1, min( NSHAFT_el, NBODBTWT_el ) !!NBODBTWT_el - 1
      do nbsub   = 1, body(nbod)%NBODSUB_el
         nb_el   =    body(nbod)%NBODTGLB_el(nbsub)

!--------- For all non-zero dofs
         do ndof = 7, subbody(nb_el)%NDFTB_el
            jj   = jj + 1
            j    = ACCU%NDFTACC_el(nb_el-1) + ndof
     
            AM_el ( 1:NDFT0_el, jj ) = AM_el ( 1:NDFT0_el, j )
            AC_el ( 1:NDFT0_el, jj ) = AC_el ( 1:NDFT0_el, j )
            AK_el ( 1:NDFT0_el, jj ) = AK_el ( 1:NDFT0_el, j )
         enddo !ndof
      enddo !nbsub
   enddo !nbod


!--- Tower Foundation [only for tower's 6dofs of 1st subbody]
   if (NBODBTWT_el==NTOWER_el) then

         nbod    = NBODBTWT_el
         nbsub   = 1
         nb_el   = body(nbod)%NBODTGLB_el(nbsub)

!------ For all non-zero dofs
      if (ITYPE_FOUND_el > 0) then
            nbc  = 1
         do ndof = 1, 6
            if ((boundco (nb_el)%nbcon(nbc     ) == 0).and. &
                (boundco (nb_el)%indx (nbc,ndof) == 1)        ) cycle

            jj   = jj + 1
            j    = ACCU%NDFTACC_el(nb_el-1) + ndof
     
            AM_el ( 1:NDFT0_el, jj ) = AM_el ( 1:NDFT0_el, j )
            AC_el ( 1:NDFT0_el, jj ) = AC_el ( 1:NDFT0_el, j )
            AK_el ( 1:NDFT0_el, jj ) = AK_el ( 1:NDFT0_el, j )
         enddo
      endif !foundation

         do ndof = 7, subbody(nb_el)%NDFTB_el
            jj   = jj + 1
            j    = ACCU%NDFTACC_el(nb_el-1) + ndof
         
            AM_el ( 1:NDFT0_el, jj ) = AM_el ( 1:NDFT0_el, j )
            AC_el ( 1:NDFT0_el, jj ) = AC_el ( 1:NDFT0_el, j )
            AK_el ( 1:NDFT0_el, jj ) = AK_el ( 1:NDFT0_el, j )
         enddo
   endif !NBODBTWT_el=NTOWER_el


!--- Jacket dependent dofs
   do nb_el=1+NBODTWT_el,NBODT_el

!------ default values if no bc is applied
      n_rem      = 0
      indx1(1:6) = 0
      indx2(1:6) = 0

      do nbc = 1, boundco (nb_el)%nbcpb
         if      (boundco (nb_el)%nod (nbc) == 1) then
            do ndof = 1, 6
               if (boundco (nb_el)%indx(nbc,ndof) == 0) cycle
               n_rem = n_rem+1
               indx1(ndof) = 1
            enddo
         else !if(boundco (nb_el)%nod (nbc) == NNPE_el) then
            do ndof = 1, 6
               if (boundco (nb_el)%indx(nbc,ndof) == 0) cycle
               n_rem = n_rem+1
               indx2(ndof) = 1
            enddo
         endif
      enddo !nbc

      do ndof = 1, 6
         if (indx1(ndof) == 1) cycle  !keep the non dependent dofs
         jj   = jj + 1
         j    = ACCU%NDFTACC_el(nb_el-1) + ndof

         AM_el ( 1:NDFT0_el, jj ) = AM_el ( 1:NDFT0_el, j )
         AC_el ( 1:NDFT0_el, jj ) = AC_el ( 1:NDFT0_el, j )
         AK_el ( 1:NDFT0_el, jj ) = AK_el ( 1:NDFT0_el, j )
      enddo
      do ndof =  7, subbody(nb_el)%NDFTB_el - 6
         jj   = jj + 1
         j    = ACCU%NDFTACC_el(nb_el-1) + ndof

         AM_el ( 1:NDFT0_el, jj ) = AM_el ( 1:NDFT0_el, j )
         AC_el ( 1:NDFT0_el, jj ) = AC_el ( 1:NDFT0_el, j )
         AK_el ( 1:NDFT0_el, jj ) = AK_el ( 1:NDFT0_el, j )
      enddo
      do ndof =  subbody(nb_el)%NDFTB_el - 5, subbody(nb_el)%NDFTB_el
         i    = ndof-subbody(nb_el)%NDFTB_el + 6
         if (indx2(i)    == 1) cycle  !keep the non dependent dofs
         jj   = jj + 1
         j    = ACCU%NDFTACC_el(nb_el-1) + ndof

         AM_el ( 1:NDFT0_el, jj ) = AM_el ( 1:NDFT0_el, j )
         AC_el ( 1:NDFT0_el, jj ) = AC_el ( 1:NDFT0_el, j )
         AK_el ( 1:NDFT0_el, jj ) = AK_el ( 1:NDFT0_el, j )
      enddo
   enddo !nb_el


!--- TRUSS COUPLED ---------
   call MOORINGS_COUPLED_tr ( 4, jj )


!--- Dependent Q's----------
   do nq = 1, NQS0
      jj = jj + 1
      j  = NDFBT_el + INDQ(nq)

      AM_el ( 1:NDFT0_el, jj ) = AM_el ( 1:NDFT0_el, j )
      AC_el ( 1:NDFT0_el, jj ) = AC_el ( 1:NDFT0_el, j )
      AK_el ( 1:NDFT0_el, jj ) = AK_el ( 1:NDFT0_el, j )
   enddo !nq


 END Subroutine MATRIX_REDUCT
!----------------------------------------------------------------------
!
!     Performes time integration, solving either the
!     static or the dynamic system of equations.
!       itype=1: the STATIC problem is solved K.u = Q
!       itype=2: the NEWMARK method is used for the 
!                dynamic system               M.ddu + C.du + K.u = Q
!
!    Newmark b-method:
!   ==================
!    qt=dq/dt, qtt=d^2/dt^2
!
!    qt[t+dt] =            qt[t] + dt   { (1  -gama) qtt[t] + gama qtt[t+dt] }
!    q [t+dt] = q [t] + dt qt[t] + dt^2 { (1/2-bita) qtt[t] + bita qtt[t+dt] }
!----------------------------------------------------------------------
 Subroutine Time_Integrate (it, itype, ICONV, ERR)
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   integer, intent(in ) :: it, itype
   real(8), intent(in ) :: ERR
   integer, intent(out) :: ICONV

   integer, allocatable       :: Ipivot  (:)
   real(8)                    :: c1, c2, c3, c4, c5, UPRE(NDFT0_el), UTPRE(NDFT0_el), RESDT,RESDTM, UTR1,UTR2
   integer                    :: i , is, ierr(2), N, imax
   integer                    :: nt
!--- jacket reduction
   integer,              save :: N_W1 ,N_W2 ,N_W 
   integer,              save :: N_J1 ,N_J2 ,N_J 
   integer,              save :: N_Q1 ,N_Q2 ,N_Q 
   integer,              save :: N_WQ1,N_WQ2,N_WQ
   integer                    :: LWORK
   real(8), allocatable       :: WORK    (:  )  !(LWORK)
   integer, allocatable, save :: Ipivot_J(:  )  !(N_J       )
   real(8), allocatable, save :: AK_J    (:,:)  !(N_J  ,N_J )
   real(8), allocatable, save :: AK_JINV (:,:)  !(N_J  ,N_J )
   real(8), allocatable       :: X_J0    (:  )  !(N_J       )
   real(8), allocatable       :: X_JW    (:,:)  !(N_J  ,N_W )
   real(8), allocatable       :: X_JQ    (:,:)  !(N_J  ,N_Q )
   real(8), allocatable       :: AK_WQ   (:,:)  !(N_WQ ,N_WQ)
   real(8), allocatable       :: AQ_WQ   (:  )  !(N_WQ      )


!--- Form the system's matrices, depending on the solution type.
   if     (itype == 2) then
!--- Newmark-b method
                   c1 = (0.5d0 - BITA_el )  * DT_el**2
                   c2 = (1.0d0 - GAMMA_el)  * DT_el
                   c3 =  1.d0    / (BITA_el * DT_el**2)
                   c4 = GAMMA_el / (BITA_el * DT_el   )
                   c5 = GAMMA_el            * DT_el
                   N  = NDFT0_el

      UPRE  (1:N    ) = UTP_el (INDSYSB(1:N))  +  UTP1_el(INDSYSB(1:N))*DT_el  +  UTP2_el(INDSYSB(1:N)) * c1
      UTPRE (1:N    ) = UTP1_el(INDSYSB(1:N))                                  +  UTP2_el(INDSYSB(1:N)) * c2
      AK_el (1:N,1:N) = AK_el(1:N,1:N)       + &
                        AC_el(1:N,1:N) * c4  + &
                        AM_el(1:N,1:N) * c3
      AQ_el (1:N    ) = AQ_el(1:N) + matmul( AM_el(1:N,1:N),              c3 * (UPRE(1:N)-UT_el(INDSYSB(1:N))) + UT2_el(INDSYSB(1:N)) ) &
                                   - matmul( AC_el(1:N,1:N), UTPRE(1:N) - c4 * (UPRE(1:N)-UT_el(INDSYSB(1:N))) - UT1_el(INDSYSB(1:N)) )
!--- AK_el check
!if(it==1)then
! nt=int(TIME/DT_el)
! write(*,*)'nt in Time_integrate',nt
! do i=NDFBTWT0_el+1, NDFBT0_el
!   write(10000+nt,5) (AK_el(i,is),is=NDFBTWT0_el+1, NDFBT0_el)
! enddo!i
!5 format(10000e28.17)
!endif
   elseif (itype == 1 ) then
!--- Static solution
                   c3 = 0.d0
                   c5 = 0.d0
                   N  = NDFT0_el
      UPRE  (1:N    ) = 0.d0;
      UTPRE (1:N    ) = 0.d0;
   endif


   if (IREDUCTJa_el ==1 .and.IMODAL<=0.and.IEIG<=0.and.ICALCINIT<=0) then
      nt=int(TIME/DT_el)
    if(it==1.and.nt==1)then

      write (* ,*)'initialize the reduced solution for jacket'
      write (10,*)'initialize the reduced solution for jacket'

      N_W1  =        1;  N_W2  = NDFBTWT0_el;  N_W  = N_W2  - N_W1  + 1 !WT
      N_J1  = N_W2 + 1;  N_J2  = NDFBT0_el  ;  N_J  = N_J2  - N_J1  + 1 !JA
      N_Q1  = N_J2 + 1;  N_Q2  = NDFT0_el   ;  N_Q  = N_Q2  - N_Q1  + 1 !Qs
      N_WQ1 = N_W2 + 1;  N_WQ2 = N_W2 + N_Q ;  N_WQ = N_WQ2 - N_W1  + 1 !WT+Qs

      write(*,'(a,3i5)')'WT',N_W1 , N_W2 , N_W
      write(*,'(a,3i5)')'JA',N_J1 , N_J2 , N_J
      write(*,'(a,3i5)')'Qs',N_Q1 , N_Q2 , N_Q
      write(*,'(a,3i5)')'WQ',N_W1 , N_WQ2, N_WQ

      LWORK = 8*N_J
     
      Allocate ( WORK    (LWORK     ) )
      Allocate ( Ipivot_J(N_J       ) )
      Allocate ( AK_J    (N_J  ,N_J ) )
      Allocate ( AK_JINV (N_J  ,N_J ) )
     
      AK_J(1:N_J,1:N_J) = AK_el (N_J1:N_J2,N_J1:N_J2)

!------ LU decomposition and Inverse only once for jacket dofs
      call DGETRF (     N_J, N_J, AK_J   , N_J, Ipivot_J,              ierr(1))
      AK_JINV = AK_J;
      call DGETRI (     N_J,      AK_JINV, N_J, Ipivot_J, WORK, LWORK, ierr(2))
     
      if (ierr(1)>0) write(*,*) 'LU Factorization     error in Ja',ierr(1)
      if (ierr(2)>0) write(*,*) 'Inverse calculation  error in Ja',ierr(2)

      Deallocate ( WORK )
     endif!it=1,nt=1
      Allocate ( X_J0    (N_J       ) )
      Allocate ( X_JW    (N_J  ,N_W ) )
      Allocate ( X_JQ    (N_J  ,N_Q ) )
      Allocate ( AK_WQ   (N_WQ ,N_WQ) )
      Allocate ( AQ_WQ   (N_WQ      ) )

      X_J0(1:N_J      ) = AQ_el (N_J1:N_J2)
      X_JW(1:N_J,1:N_W) = -matmul( AK_JINV(1:N_J,1:N_J),AK_el(N_J1:N_J2,N_W1:N_W2) )
      X_JQ(1:N_J,1:N_Q) = -matmul( AK_JINV(1:N_J,1:N_J),AK_el(N_J1:N_J2,N_Q1:N_Q2) )

      call DGETRS ('N', N_J,  1 , AK_J   , N_J, Ipivot_J, X_J0, N_J  , ierr(2))

      if (ierr(2)>0) write(*,*) 'LU Back Substitution error in Ja',ierr(2)

!------ Set AK_WQ, AQ_WQ
      AK_WQ(N_W1 :N_W2 ,N_W1 :N_W2 ) = AK_el(N_W1:N_W2,N_W1:N_W2) + matmul( AK_el(N_W1:N_W2,N_J1:N_J2), X_JW(1:N_J,1:N_W) )  !A_WW
      AK_WQ(N_WQ1:N_WQ2,N_W1 :N_W2 ) = AK_el(N_Q1:N_Q2,N_W1:N_W2) + matmul( AK_el(N_Q1:N_Q2,N_J1:N_J2), X_JW(1:N_J,1:N_W) )  !A_QW

      AK_WQ(N_W1 :N_W2 ,N_WQ1:N_WQ2) = AK_el(N_W1:N_W2,N_Q1:N_Q2) + matmul( AK_el(N_W1:N_W2,N_J1:N_J2), X_JQ(1:N_J,1:N_Q) )  !A_WQ
      AK_WQ(N_WQ1:N_WQ2,N_WQ1:N_WQ2) = AK_el(N_Q1:N_Q2,N_Q1:N_Q2) + matmul( AK_el(N_Q1:N_Q2,N_J1:N_J2), X_JQ(1:N_J,1:N_Q) )  !A_QQ

      AQ_WQ(N_W1 :N_W2             ) = AQ_el(N_W1:N_W2          ) - matmul( AK_el(N_W1:N_W2,N_J1:N_J2), X_J0(1:N_J      ) )  !B_W
      AQ_WQ(N_WQ1:N_WQ2            ) = AQ_el(N_Q1:N_Q2          ) - matmul( AK_el(N_Q1:N_Q2,N_J1:N_J2), X_J0(1:N_J      ) )  !B_Q

      Allocate   ( Ipivot (N_WQ) )

!------ Solve the reduced linear system
      call DGETRF (     N_WQ, N_WQ, AK_WQ, N_WQ, Ipivot,              ierr(1))
      call DGETRS ('N', N_WQ,  1  , AK_WQ, N_WQ, Ipivot, AQ_WQ, N_WQ, ierr(2))
     
      if (ierr(1)>0) write(*,*) 'LU Factorization     error in WQ',ierr(1)
      if (ierr(2)>0) write(*,*) 'LU Back Substitution error in WQ',ierr(2)
     
      Deallocate ( Ipivot )

!------ Substitute AQ_el
      AQ_el(N_W1:N_W2) = AQ_WQ(N_W1 :N_W2 )
      AQ_el(N_Q1:N_Q2) = AQ_WQ(N_WQ1:N_WQ2)
      AQ_el(N_J1:N_J2) = X_J0 (    1:N_J  ) + matmul( X_JW(1:N_J,1:N_W),AQ_WQ(N_W1 :N_W2 ) ) &
                                            + matmul( X_JQ(1:N_J,1:N_Q),AQ_WQ(N_WQ1:N_WQ2) )

      Deallocate ( X_J0, X_JW, X_JQ, AK_WQ, AQ_WQ )
   else
!--- Solve the linear system of equations using LU decomposition
      Allocate ( Ipivot (NDFT0_el)   )
     
      call DGETRF (     NDFT0_el, NDFT0_el, AK_el, NDFT_el, Ipivot,                 ierr(1))
      call DGETRS ('N', NDFT0_el,   1     , AK_el, NDFT_el, Ipivot, AQ_el, NDFT_el, ierr(2))
     
      if (ierr(1)>0) write(*,*) 'LU Factorization    ',ierr(1)
      if (ierr(2)>0) write(*,*) 'LU Back Substitution',ierr(2)
     
      Deallocate ( Ipivot )
   endif !IREDUCTJa_el=1


!--- calculate the maximum and the mean errors, perform convergence check
   RESDT  = 0.d0
   RESDTM = 0.d0
   imax   = 1

   do i = 1, NDFT0_el
      RESDT  = RESDT + dabs(AQ_el(i))
      if (RESDTM<dabs(AQ_el(i))) then
         if (ICMOD==0.and.i==NDFBT0_el+NQSW) cycle !eliminate omega q dof
         imax   = i
         RESDTM = dabs(AQ_el(i))
      endif
   enddo

   if (((RESDTM >= 0.d0).or.(RESDTM <= 0.d0)).eqv.(.FALSE.)) then
      write (10,*) '** NaN Error **'
      stop
   endif

                  RESDT = RESDT/dble(NDFT0_el  )
   if (ICMOD==0)  RESDT = RESDT/dble(NDFT0_el-1) !eliminate omega q dof

   if (RESDTM < ERR) ICONV = 1

   write (* ,1) it, RESDT, RESDTM, imax, UT_el(INDSYSB(imax)) + AQ_el(imax)
   write (10,1) it, RESDT, RESDTM, imax, UT_el(INDSYSB(imax)) + AQ_el(imax)

1  format (2x,'it = ', i3,4x,'  RESDT =',2(e23.16,2x),i5,  4x,e23.16)


!--- Substitute the solution for all active dofs
   do i           = 1, NDFT0_el
      is          = INDSYSB (i)
   
      UT0_el (is) = AQ_el(i)
      UT_el  (is) = UT_el (is) + AQ_el(i)*RLXs

      UTR2        = c3 * ( UT_el(is)-UPRE(i) )
      UT02_el(is) = UTR2     - UT2_el(is)
      UT2_el (is) = UTR2

      UTR1        = UTPRE(i) + c5 * UTR2
      UT01_el(is) = UTR1     - UT1_el(is)
      UT1_el (is) = UTR1    
   enddo

   call GET_REDUSED_DOFS


 END Subroutine Time_Integrate
!----------------------------------------------------------------------
 Subroutine GET_REDUSED_DOFS
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   integer :: i, j, ju, iq, iq0, ju0, nb_el, nbc, is
   integer :: nbod,nbsub,iel,inodel,inode,iecon,incon,inodc
   integer :: ibcon,isbcon,ibsubcon


   if (IREDUCT_el == 0) return ! No reduction

!--- WT: Dependent Qs
   do iq = 1, NQS
      if (INDSYSQ (iq) < 1) cycle

      i  = NDFBT_el + iq
      is = INDSYSQ (iq)

      UT0_el (i) = UT0_el (is)
      UT01_el(i) = UT01_el(is)
      UT02_el(i) = UT02_el(is)

      UT_el  (i) = UT_el  (is)
      UT1_el (i) = UT1_el (is)
      UT2_el (i) = UT2_el (is)
   enddo

!--- JACKET: Dependent elastic dofs
   do nbod  = 1+NBODBTWT_el, NBODBT_el
   do nbsub = 1, body(nbod)%NBODSUB_el

      nb_el = body(nbod)%NBODTGLB_el(nbsub)

      do nbc = 1, boundco (nb_el)%nbcpb

         ibcon    = boundco (nb_el)%nbcon (nbc)
         if (ibcon == 0) cycle

         iel      = boundco (nb_el)%nel (nbc)
         inodel   = boundco (nb_el)%nod (nbc)
         inode    =                             subbody(nb_el )%NODMB       (iel,inodel)
         iq0      = ACCU%NDFTACC_el(nb_el-1 ) + subbody(nb_el )%NDFPNBACC_el(inode-1   )

         ibsubcon = 1
         isbcon   = body    (ibcon)%NBODTGLB_el(ibsubcon)
         incon    = boundco (nb_el)%nodcon(nbc)
         iecon    = boundco (nb_el)%nelcon(nbc)
         inodc    =                             subbody(isbcon)%NODMB       (iecon,incon)
         ju0      = ACCU%NDFTACC_el(isbcon-1) + subbody(isbcon)%NDFPNBACC_el(inodc-1    )

         do i  = 1, 6
            if (boundco (nb_el)%indx(nbc,i) == 0) cycle
            iq = iq0+i

            UT_el  (iq) = 0.d0
            UT1_el (iq) = 0.d0
            UT2_el (iq) = 0.d0
            UT0_el (iq) = 0.d0
            UT01_el(iq) = 0.d0
            UT02_el(iq) = 0.d0

            do j  = 1, 6
               ju = ju0+j

               UT_el  (iq) = UT_el  (iq) + UT_el  (ju) * boundco (nb_el)%Atrans1 (nbc,i,j)
               UT1_el (iq) = UT1_el (iq) + UT1_el (ju) * boundco (nb_el)%Atrans1 (nbc,i,j)
               UT2_el (iq) = UT2_el (iq) + UT2_el (ju) * boundco (nb_el)%Atrans1 (nbc,i,j)

               UT0_el (iq) = UT0_el (iq) + UT0_el (ju) * boundco (nb_el)%Atrans1 (nbc,i,j)
               UT01_el(iq) = UT01_el(iq) + UT01_el(ju) * boundco (nb_el)%Atrans1 (nbc,i,j)
               UT02_el(iq) = UT02_el(iq) + UT02_el(ju) * boundco (nb_el)%Atrans1 (nbc,i,j)
            enddo !j
         enddo !i

      enddo !nbc
   enddo !nbsub
   enddo !nbod


 END Subroutine GET_REDUSED_DOFS
!--------------------------------------------------------------------------------
!
!  Subrouine : GENLOAD  ----------------
!
!  Equation for rotational speed
!
!----------------------------------------------------------------------
 Subroutine GENLOAD
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   real(8) :: AMB2B    (1,NDFPEM_el)
   real(8) :: ACB2B    (1,NDFPEM_el)
   real(8) :: AKB2B    (1,NDFPEM_el)
   real(8) :: AMB2BNQ  (1,NQS      )
   real(8) :: ACB2BNQ  (1,NQS      )
   real(8) :: AKB2BNQ  (1,NQS      )
   real(8) :: AFB2B    (1          )
   integer :: i, j, NDFPE, iq, nq
   integer :: nbod, nbsub, nb_el, nn_el, nel, IMAX, JMAX, IQRayl
   real(8) :: COEFM, COEFK


   if (NBODBTWT_el == NBLADE_el) return


   nbod  = NBLADE_el + 1
   nbsub = 1
   nel   = 1
   nb_el = body(nbod)%NBODTGLB_el(nbsub)
   nn_el = body(nbod)%NNODTGLB_el(nbsub)
   NDFPE = subbody(nb_el)%NDFPE_el
   i     = IMDOF_el(5)
   j     = ACCU%NDFTACC_el   ( nb_el-1 )

      
!--- Loads of Wind Turbine
   AMB2B  ( 1, 1:NDFPE ) = subbody(nb_el)%AMLOC_el ( nel, i, 1:NDFPE )
   ACB2B  ( 1, 1:NDFPE ) = subbody(nb_el)%ACLOC_el ( nel, i, 1:NDFPE )
   AKB2B  ( 1, 1:NDFPE ) = subbody(nb_el)%AKLOC_el ( nel, i, 1:NDFPE )
   AFB2B  ( 1          ) = subbody(nb_el)%AFLOC_el ( nel, i          ) &
                            + TGenLSS_el + TMLossLSS_el + TBrakeLSS_el ! QTC1_el(6)*RAT_GEAR

!--- For all Q's [without Rayleigh Damping]
   do iq = 1, QSnode(nn_el)%Tot_el(2)
      nq = QSnode(nn_el)%IORD_el(iq)

      AMB2BNQ ( 1, nq ) = subbody(nb_el)%AMLOCNQ_el ( nel, i, iq )
      ACB2BNQ ( 1, nq ) = subbody(nb_el)%ACLOCNQ_el ( nel, i, iq )
      AKB2BNQ ( 1, nq ) = subbody(nb_el)%AKLOCNQ_el ( nel, i, iq )
   enddo

!------------------- Damping term for stall-regulated WT
                   nq   = NQSW
      ACB2BNQ ( 1, nq ) = ACB2BNQ ( 1, nq ) + TGenLSS_D_el
      AMB2BNQ ( 1, nq ) = AMB2BNQ ( 1, nq ) + TGenLSS_M_el


!--- LOADS Communicated to Generator (A)
   i    = NDFBT_el + NQSW-1
   j    = ACCU%NDFTACC_el(nb_el-1)

   IQRayl  = 2   !0: all qs Rayleigh  , 1: body qs with Rayleigh, 2: no Rayleigh to qs
   COEFM   = 0.d0
   COEFK   = 0.d0
   IMAX    = 1
   JMAX    = subbody(nb_el)%NDFPE_el

   call ASSEMBLY_MATRIX_el  ( AMB2B  , ACB2B  , AKB2B  , AFB2B         ,&
                              AMB2BNQ, ACB2BNQ, AKB2BNQ, COEFM , COEFK ,&
                                       nn_el  , i      , j     , IQRayl,&
                              IMAX   , JMAX   , 1      , NDFPEM_el       )


 END Subroutine GENLOAD
!----------------------------------------------------------------------
!
!  Subrouine : ACTYAWEQ ----------------
!
!  Equation for free yaw movement
!
!----------------------------------------------------------------------
 Subroutine ACTYAWEQ
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   real(8) :: Aba      (3,3        )
   real(8) :: Aba0     (3,3,NQS    )
   real(8) :: Rba      (3          )
   real(8) :: Rba0     (3  ,NQS    )
   real(8) :: RGba     (3          )
   real(8) :: AMB2B    (6,NDFPEM_el)
   real(8) :: ACB2B    (6,NDFPEM_el)
   real(8) :: AKB2B    (6,NDFPEM_el)
   real(8) :: AMB2BNQ  (6,NQS      )
   real(8) :: ACB2BNQ  (6,NQS      )
   real(8) :: AKB2BNQ  (6,NQS      )
   real(8) :: AFB2B    (6          )
   integer :: i, j, NDFPE, iq, nq, i0
   integer :: nbodB, nnsubB, nbsubB, nnB_el, nbB_el, nelB, nodB
   integer :: nbodA,                 nnA2_el
   integer :: IcalcR


   if (NBODBTWT_el /= NTOWER_el) return

!--- Body (B)
   nbodB     = NBLADE_el + 1
   nnsubB    = 1
   nbsubB    = 1
   nelB      = 1
   nodB      = 1
   nnB_el    = body(nbodB)%NNODTGLB_el(nnsubB)
   nbB_el    = body(nbodB)%NBODTGLB_el(nbsubB)
   NDFPE     = subbody(nbB_el)%NDFPE_el
   j         = ACCU%NDFTACC_el (nbB_el-1)

!--- Body (A)
   nbodA     = NBLADE_el + 2
   nnA2_el   = body(nbodA)%NNODTGLB_el(body(nbodA)%NBODSUB_el+1 )
   i         = NDFBT_el + NQSY
   i0        = IMDOF_el(4) !x-component, ATyaw is defined wrt the shaft c.s.


!--- Aba, Aba0, Rba, Rba0
!--- Yaw Axis, Tower top
   RGba(1:3    ) = transf_mat(nnB_el)%R_el(1:3) - transf_mat(nnA2_el)%R_el(1:3)
   Aba (1:3,1:3) = matmul ( ATyaw_el(1:3,1:3), transf_mat(nnB_el)%A_el(1:3,1:3) )
   Rba (1:3    ) = matmul ( ATyaw_el(1:3,1:3), RGba(1:3)                        )

   do iq = 1, QSnode(nnA2_el)%Tot_el(2)
      nq = QSnode(nnB_el)%IORD_el(iq)

      Aba0(1:3,1:3,nq) = matmul( ATyaw_el (   1:3,1:3), transf_mat(nnB_el)%A0_el(iq,1:3,1:3) ) + &
                         matmul( AT0yaw_el(nq,1:3,1:3), transf_mat(nnB_el)%A_el (   1:3,1:3) )

      Rba0(1:3,    nq) = matmul( (transf_mat(nnB_el)%R0_el(iq,1:3) - transf_mat(nnA2_el)%R0_el(iq,1:3) ), ATyaw_el (   1:3,1:3) ) + &
                         matmul(  RGba                    (   1:3)                                      , AT0yaw_el(nq,1:3,1:3) )
   enddo!iq 

   do iq = QSnode(nnA2_el)%Tot_el(2)+1, QSnode(nnB_el)%Tot_el(2)
      nq = QSnode(nnB_el)%IORD_el(iq)

      Aba0(1:3,1:3,nq) = matmul( ATyaw_el (   1:3,1:3), transf_mat(nnB_el)%A0_el(iq,1:3,1:3) ) + &
                         matmul( AT0yaw_el(nq,1:3,1:3), transf_mat(nnB_el)%A_el (   1:3,1:3) )

      Rba0(1:3,    nq) = matmul( (transf_mat(nnB_el)%R0_el(iq,1:3)                                     ), ATyaw_el (   1:3,1:3) ) + &
                         matmul(  RGba                    (   1:3)                                      , AT0yaw_el(nq,1:3,1:3) )
   enddo!iq 

   IcalcR  = 1   !0: no R calculations, 1: R calculations

   call LOC_LOADS_el        ( AMB2B  , ACB2B  , AKB2B  , AFB2B,&
                              AMB2BNQ, ACB2BNQ, AKB2BNQ       ,&
                              nbB_el  ,nnB_el  ,nelB   , nodB    )

   call LOADS_TRANS_el      ( AMB2B  , ACB2B  , AKB2B  , AFB2B,&
                              AMB2BNQ, ACB2BNQ, AKB2BNQ       ,&
                              Aba0   , Aba    , Rba0   , Rba  ,&
                                       nnB_el , IcalcR        ,&
                              NDFPE  , NDFPEM_el                 )

   AM_el( i, j+1:j+NDFPE ) = AM_el( i, j+1:j+NDFPE ) + AMB2B ( i0, 1:NDFPE )
   AC_el( i, j+1:j+NDFPE ) = AC_el( i, j+1:j+NDFPE ) + ACB2B ( i0, 1:NDFPE )
   AK_el( i, j+1:j+NDFPE ) = AK_el( i, j+1:j+NDFPE ) + AKB2B ( i0, 1:NDFPE )
   AQ_el( i              ) = AQ_el( i              ) + AFB2B ( i0          )

   do iq = 1, QSnode(nnB_el)%Tot_el(2)
      nq = QSnode(nnB_el)%IORD_el(iq)
      j  = NDFBT_el + nq

      AM_el( i, j ) = AM_el( i, j ) + AMB2BNQ ( i0, nq )
      AC_el( i, j ) = AC_el( i, j ) + ACB2BNQ ( i0, nq )
      AK_el( i, j ) = AK_el( i, j ) + AKB2BNQ ( i0, nq )
   enddo !iq

!--- Stiffness and Damping of the yaw mechanism
   AK_el ( i, i ) = AK_el ( i, i ) + CSTIF_YAW
   AC_el ( i, i ) = AC_el ( i, i ) + CDAMP_YAW
   AQ_el ( i    ) = AQ_el ( i    ) - CDAMP_YAW*UT1_el(i) &
                                   - CSTIF_YAW*UT_el (i)


 END Subroutine ACTYAWEQ
!----------------------------------------------------------------------
!
!  NAC_LX_el           [m] Longitudinal length 
!  NAC_LY_el           [m] Side         length
!  NAC_LZ_el           [m] Height
!
!  NAC_FROND_X         [m] nacelle offset wrt tower top x global axis [fore-aft]
!  NAC_FROND_Y         [m] nacelle offset wrt tower top y global axis [side    ]
!  NAC_FROND_Z         [m] nacelle offset wrt tower top z global axis [up-down ]
!
!  CDRAG_NAC_side_el   [-] Cd coeff of the side  surface
!  CDRAG_NAC_fore_el   [-] Cd coeff of the front surface
!
!----------------------------------------------------------------------
 Subroutine DRAG_nac
!----------------------------------------------------------------------

 Use Cbeam
 Use Craft

   implicit none

   real(8) :: F(3), M(3), A(3,3)
   real(8) :: HTAL, ALLOC, PHIX, PHIZ
   real(8) :: AQ      (NEQPEM_el          )
   real(8) :: SHAPET  (NDFPEM_el,NEQPEM_el)
   real(8) :: SHAPE   (NEQPEM_el,NDFPEM_el)
   real(8) :: AFLOC0  (NDFPEM_el          )

   integer :: nbod, nbsub, nb_el, nn_el, nel, nn, i, j, NDFPE, NEQPE
   real(8), dimension(3)   :: x,dx,ddx,dxl,ddxl, AKSI0                 !--- Tower acceleration
   real(8), dimension(3)   :: Rcalc, RcalcL, Vnac_wind, Vnac_rigid, Vnac_tot

   real(8)                 :: alpha, cl,cd
!--- Rotation matrices
   real(8), dimension(3,3) :: Ayaw,  dAyaw, ATnacelle
   real(8)                 :: tang(3), norm(3), area, unitv(3)
   real(8)                 :: Fnac_lift, Fnac_drag


   if (IDRAG_NAC == 0 ) return

   ATnacelle   = Transpose(Anacelle_el)
!--- set point of calculation wrt tower top
   Rcalc (1  ) = NAC_RX_el
   Rcalc (2  ) = NAC_RY_el
   Rcalc (3  ) = NAC_RZ_el !- Hoff_drag
   RcalcL(1:3) = matmul(ATnacelle, Rcalc)

!--- get nacelle kinematics (qq atm tower top)
   nbod     = NTOWER_el
   nbsub    = body(nbod)%NBODSUB_el
   nb_el    = body(nbod)%NBODTGLB_el(nbsub)
   nel      = subbody(nb_el)%NTEB_el
   AKSI0(:) = 0.d0;
   AKSI0(2) = subbody(nb_el)%HTA_el (nel+1)

   call beam_kinematics (nbod,nbsub,nel,AKSI0,x,dx,ddx,dxl,ddxl)

!--- nacelle local velocities
   Vnac_wind (1:3) = matmul ( ATnacelle(1:3,1:3), HUB_VEL(1:3) )
   Vnac_rigid(1:3) = matmul ( ATnacelle(1:3,1:3), dx     (1:3) )
   Vnac_tot  (1:3) = Vnac_wind (1:3)-Vnac_rigid (1:3)  

!--- find angle of attack
   alpha = datan2(-Vnac_tot(3),-Vnac_tot(2) ) * R2D
   alpha = -alpha
   if (alpha > 180.d0) alpha = alpha -360.d0
   if (alpha <-180.d0) alpha = alpha +360.d0

!--- find area normal to inflow velocity
   call ROT_MATRIX1 ( 3, -WINDYAW, Ayaw, dAyaw )
   unitv      = 0.d0;
   unitv(1  ) = 1.d0
   tang (1:3) = matmul(ATnacelle,matmul(Ayaw,unitv))

   unitv      = 0.d0;
   unitv(2  ) =-1.d0
   norm (1:3) = matmul(ATnacelle,matmul(Ayaw,unitv))
   area       = ( NAC_LX_el*dabs(norm(2)) + NAC_LY_el*dabs(norm(3)) ) * NAC_LZ_el

!--- find Cl, Cd
   call LIN_INT ( alpha, cl , Alpha_NAC, CL_NAC, NANG_NAC, 900 )
   call LIN_INT ( alpha, cd , Alpha_NAC, CD_NAC, NANG_NAC, 900 )

!--- Calculate local nacelle force and moment
   Fnac_drag  = AIRDEN/2d0 * cd * area  * (Vnac_tot(2)**2+Vnac_tot(3)**2)
   Fnac_lift  = AIRDEN/2d0 * cl * area  * (Vnac_tot(2)**2+Vnac_tot(3)**2)

   F(1:3)     = tang(1:3)*Fnac_drag + norm(1:3)*Fnac_lift

   call  EXTEPR ( RcalcL, F, M )

!--- Output data
# if   ASCII == 1
   open (67,file='nacel.dat',access='append')
    write (67,'(90f15.5)') &
# elif ASCII == 0
   open (67,file='nacel.bin',access='append',form='UNFORMATTED')
    write (67            ) &
# endif
    sngl(TIME          )       ,&
    sngl(alpha         )       ,&
    sngl(cl            )       ,&
    sngl(cd            )       ,&
   (sngl(F          (j)),j=1,3),&
   (sngl(Vnac_rigid (j)),j=1,3),&
   (sngl(Vnac_tot   (j)),j=1,3),&
   (sngl(M          (j)),j=1,3),&
    sngl(Fnac_lift     )       ,& 
    sngl(Fnac_drag     )       ,&
    sngl(area          )
   close (67)

!--- Assembly and Store
   nbod  = NTOWER_el
   nbsub = body(nbod)%NBODSUB_el
   nb_el = body(nbod)%NBODTGLB_el (nbsub)
   nn_el = body(nbod)%NNODTGLB_el (nbsub)
   nel   = subbody(nb_el)%NTEB_el
   nn    =                              subbody(nb_el)%NODMB       (nel,1) !nnod)
   i     = ACCU%NDFTACC_el(nb_el-1)  +  subbody(nb_el)%NDFPNBACC_el(nn-1 )
   HTAL  = 1.d0
   ALLOC = subbody(nb_el)%ALENG_el(nel)
   NDFPE = subbody(nb_el)%NDFPE_el
   NEQPE = subbody(nb_el)%NEQPE_el
   PHIX  = subbody(nb_el)%PHIX_el (nel)
   PHIZ  = subbody(nb_el)%PHIZ_el (nel)

   call SSHAPEFUNC15 ( HTAL, ALLOC, PHIX, PHIZ,& 
                       SHAPE, NDFPE, NEQPE       )

!--- transformation matrix: nacelle to tower
   A(1:3,1:3) = matmul( transf_mat(nn_el)%AT_el(1:3,1:3), Anacelle_el (1:3,1:3) )

!--- Nacelle loads applied on tower top wrt to tower c.s.
   F  = matmul(A,F)
   M  = matmul(A,M)

   AQ (IMDOF_el(1:3)) = F          ( 1:3         )
   AQ (IMDOF_el(4:6)) = M          ( 1:3         )
   SHAPET             = transpose  ( SHAPE       )
   AFLOC0             = matmul     ( SHAPET, AQ  )

!--- STORE_LOC_MATRIX_el
   subbody(nb_el)%AFLOC_el ( nel, 1:NDFPEM_el ) = subbody(nb_el)%AFLOC_el ( nel, 1:NDFPEM_el ) + AFLOC0 ( 1:NDFPEM_el )

!--- ASSEMBLY_MATRIX_el
                  AQ_el    ( i+1 : i+NDFPEM_el) =                AQ_el    ( i+1: i+NDFPEM_el ) + AFLOC0 ( 1:NDFPEM_el )


 END Subroutine DRAG_nac
!----------------------------------------------------------------------
 Subroutine DRAG_tow ( AQ, HTA1, DIAM, nn_el )
!----------------------------------------------------------------------

 Use Cbeam
 Use Craft

   implicit none

   real(8) :: XG(3),UG(3),UL(3),FL(3)
   real(8) :: HTA1,DIAM,Z
   real(8) :: AQ (NEQPEM_el)
   integer :: nn_el


   if (IDRAG_TOW == 0 ) return

   Z  = transf_mat(nn_el)%R_el(3) + transf_mat(nn_el)%A_el(3,2) * HTA1

   if ( Z < Z_DRAG_LIM ) return

   XG(1:2)  = 0.d0;
   XG(  3)  = Z

   call UINFLOW (TIME, XG, UG, 0, 0) !turbulance, tshadow not inlcuded

   UL (1:3) = matmul (transf_mat(nn_el)%AT_el(1:3,1:3), UG (1:3) )
   FL (1  ) = AIRDEN/2d0 * CDRAG_TOW_el * DIAM * dabs(UL(1)) * UL(1)
   FL (2  ) = 0.d0
   FL (3  ) = AIRDEN/2d0 * CDRAG_TOW_el * DIAM * dabs(UL(3)) * UL(3)

   AQ(IMDOF_el(1:3)) = AQ(IMDOF_el(1:3)) - FL (1:3)

!  if (nel ==15) then
!      open(24,file='dragtow15.dat',access='append')
!      write(24,11) TIME,FL(1:3),UL(1:3)
!      close(24)
!11   format ( 150f18.6 )
!  endif


 END Subroutine DRAG_tow
!----------------------------------------------------------------------
 Subroutine VIV_tow ( AQ, HTA1, DIAM, nn_el )
!----------------------------------------------------------------------

 Use Cbeam
 Use Craft

   implicit none

   real(8), intent(inout) :: AQ (NEQPEM_el)
   real(8), intent(in   ) :: HTA1, DIAM
   integer, intent(in   ) :: nn_el

   real(8), save :: Zlow(2), Zup(2), clat, freq, Lzero(2)
   integer, save :: init=0 , imode=1

   integer :: i,j
   real(8) :: Z ,D, rho
   real(8) :: XG(3),UG(3),UL(3),FL(3)


30 if ( imode == 0 ) return

   if (init==0) then
      imode = 0
      open (100,file='vivtow.inp')
      read(100,*,end=30,err=30)  imode
      read(100,*) clat, freq
     do i = 1, imode
      read(100,*) Zlow(i),Zup(i),Lzero(i)
     enddo !i
      init=1
      write(10,'( a,   i3)') 'Eurocode VIV for tower mode', imode
      write(10,'( a,f15.3)') 'Clat     : ', clat
      write(10,'( a,f15.3)') 'Frequency: ', freq
     do i = 1, imode
      write(10,'(a,3f15.3)') 'Zlow/up  : ', Zlow(i),Zup(i),Lzero(i)
     enddo !i
   endif !init


      Z = HTA1 !transf_mat(nn_el)%R_el(3) + transf_mat(nn_el)%A_el(3,2) * HTA1
!  do i = 1, imode
!     if (Z>Zlow(i) .and. Z<Zup(i)) goto 1
!  enddo
      i=1; if (Z>Zlow(1) .and. Z<(Zlow(1)+Zup(1))/2.d0-Lzero(1)/2.d0) goto 1
      i=1; if (Z<Zup (1) .and. Z>(Zlow(1)+Zup(1))/2.d0+Lzero(1)/2.d0) goto 1
      i=2; if (Z>Zlow(2) .and. Z<         Zup(2)      -Lzero(2)     ) goto 1
      return


1  XG(1:2)  = 0.d0;
   XG(  3)  = Z

   call UINFLOW (TIME, XG, UG, 0, 0) !turbulance, tshadow not inlcuded
!qqUG(:)  = 0.0d0;
!qqUG(1)  =29.5d0 
   D      = DIAM    !4.472d0
   rho    = AIRDEN  !1.25d0

   UL (1:3) = matmul (transf_mat(nn_el)%AT_el(1:3,1:3), UG (1:3) ) *min(TIME/30.d0,1.d0)
   FL (:  ) = 0.d0;
   FL (3  ) = rho/2d0 * clat * D * UL(1)**2 * dcos(PI2*freq*TIME + PI*dble(i-1))
   FL (1  ) = rho/2d0 * clat * D * UL(3)**2 * dcos(PI2*freq*TIME + PI*dble(i-1))

   AQ(IMDOF_el(1:3)) = AQ(IMDOF_el(1:3)) - FL (1:3)

   if     (Z> 75.d0.and.Z< 76.d0) then
       open(24,file='vivtow1.dat',access='append')
       write(24,'(8f18.6,i3)') TIME,(FL(j),j=1,3),(UL(j),j=1,3),Z,i
       close(24)
   elseif (Z>133.d0.and.Z<134.d0) then
       open(24,file='vivtow2.dat',access='append')
       write(24,'(8f18.6,i3)') TIME,(FL(j),j=1,3),(UL(j),j=1,3),Z,i
       close(24)
   endif


 END Subroutine VIV_tow
!----------------------------------------------------------------------
 Subroutine Calc_Tow_Shaft_Angles
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   real(8) :: UT    (NDFPEM_el )         
   real(8) :: UT1   (NDFPEM_el )         
   real(8) :: UT2   (NDFPEM_el )         
   integer :: nbod, nbsub,nb_el,nel,idof


!--- Calculate Shaft's end torsion angle, velocity and acceleration
   if (NBODBTWT_el >= NSHAFT_el) then
      nb_el     = body(NBLADE_el+1)%NBODTGLB_el(1)
      nel       = subbody(nb_el)%NTEB_el
!     nn        = subbody(nb_el)%NODMB (nel, NNPE_el)
!     idof      = ACCU%NDFTACC_el(nb_el-1) + subbody(nb_el)%NDFPNBACC_el(nn-1) + 6
      idof      = subbody(nb_el)%NDFPE_el
      call LOCAL_UT_1 ( nb_el, nel, UT, UT1, UT2 )
      UThub_el  = UT_el (idof)
      UThub1_el = UT1_el(idof)
      UThub2_el = UT2_el(idof)
   else
      UThub_el  = 0.d0
      UThub1_el = 0.d0
      UThub2_el = 0.d0
   endif


!--- Calculate Tower's top tilt and yaw angles
!--- tower top twist    angle --> yaw  correction ( y axis)
!--- tower top fore-aft angle --> tilt correction (-z axis)
   if (NBODBTWT_el == NTOWER_el ) then
      nbod      = NTOWER_el
      nbsub     = body(nbod)%NBODSUB_el
      nb_el     = body(nbod)%NBODTGLB_el(nbsub)
      nel       = subbody(nb_el)%NTEB_el
!     nn        = subbody(nb_el)%NODMB (nel, NNPE_el)
!     idof      = ACCU%NDFTACC_el(nb_el-1) + subbody(nb_el)%NDFPNBACC_el(nn-1)
      idof      = subbody(nb_el)%NDFPE_el-6
      call LOCAL_UT_1 ( nb_el, nel, UT, UT1, UT2 )

      yaw_ttop  =   UT_el(idof+IMDOF_el(5))
      tilt_ttop = - UT_el(idof+IMDOF_el(6))
   else
      yaw_ttop  = 0.d0
      tilt_ttop = 0.d0
   endif
 
   if (ICASE_el >= 3) then
      yaw_ttop  = yaw_ttop  + UT_el(NDFBT_el+6)
      tilt_ttop = tilt_ttop + UT_el(NDFBT_el+5)
   endif


 END Subroutine Calc_Tow_Shaft_Angles
