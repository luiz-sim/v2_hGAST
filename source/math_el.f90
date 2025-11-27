!-------------------------------------------------------------------------
!   Subroutine :MATRLOC
! 
!  All Local Matrices [MASS, DAMPING, STIFFNESS, FORCES]
!  with linearly varying properties calculated 
!  with respect to the body FEM local system.
!  Contributions from:
!
!    1. Stiffness terms based on Timoshenko 1st-order beam theory
!    2. Inertial  terms based on multibody dynamics
!    3. Gravity
!
!-------------------------------------------------------------------------
 Subroutine MATRLOC ( AMLOC, ACLOC, AKLOC, AFLOC, AMLOCNQ, ACLOCNQ, AKLOCNQ,&
                      nb_el, nn_el, nel                                       )
!-------------------------------------------------------------------------

 Use Cbeam

   implicit none

   integer                                , intent( in) :: nb_el  , nn_el  , nel
   real(8), dimension(NDFPEM_el,NDFPEM_el), intent(out) :: AKLOC  , ACLOC  , AMLOC
   real(8), dimension(NDFPEM_el,NQS      ), intent(out) :: AMLOCNQ, ACLOCNQ, AKLOCNQ
   real(8), dimension(NDFPEM_el          ), intent(out) :: AFLOC

   real(8), dimension(NDFPEM_el,NDFPEM_el) :: AKNLLOC, AKNLLOC0
   real(8), dimension(NDFPEM_el,NDFPEM_el) :: AKLOC0,ACLOC0,AMLOC0 !!#ELEMENT,AE12, AE12T
   real(8), dimension(NDFPEM_el          ) :: AFLOC0
   real(8), dimension(NDFPEM_el          ) :: UT  ,UT1 ,UT2 
   real(8), dimension(NEQPEM_el          ) :: UTT1,UTT0,DUT0
   real(8), dimension(NEQPEM_el,NEQPEM_el) :: AM , ACM, AKM
   real(8), dimension(NEQPEM_el,NEQPEM_el) :: AK1, AK2, AK3  , AK4  , AK0
   real(8), dimension(NEQPEM_el,NEQPEM_el) ::           AK3NL, AK4NL
   real(8), dimension(NEQPEM_el          ) :: AQM, AQG, AQNL                !G:gravity, NL: nonlinear, M: inertia
   real(8), dimension(NEQPEM_el,NDFPEM_el) :: SHAPE , DSHAPE
   real(8), dimension(NDFPEM_el,NEQPEM_el) :: SHAPET, DSHAPET, SHAPETAM
   real(8) :: ATDAx2(3,3), ATDDA(3,3), ATDDR(3), AT(3,3)          !!#ELEMENT, AE(3,3), AET(3,3)
   real(8) :: ALLOC , PA01 , PA02 , cgauss, F_Grav
   real(8) :: DENS  , AMOMX, AMOMZ, RIXX, RIZZ, RIXZ, POLI, HTA   , HTAL
   real(8) :: TRACT0, XCM  , ZCM  , PHPX, PHPZ, HTA1, HTA0, F0(3)
   real(8) :: Kf(21)
   integer :: IB    , NEQPE, NDFPE, k   , n1  , n2  , IND(6), ii,jj


!--- Initialize local matrices
   AMLOC   = 0.d0;
   ACLOC   = 0.d0;
   AKLOC   = 0.d0;
   AKNLLOC = 0.d0;
   AFLOC   = 0.d0;
   AMLOCNQ = 0.d0;
   ACLOCNQ = 0.d0;
   AKLOCNQ = 0.d0;


!--- Pass Transformation matrices to local arrays
!!!#ELEMENT
!!   AE  (1:3, 1:3 ) = subbody(nb_el)% AE (nel,1:3,1:3)
!!   AET             = transpose(AE);
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
   IND(:)  = IMDOF_el(:)
!!!#ELEMENT
!!   AE12                          = 0.d0;
!!   AE12  (  IND(1:3),  IND(1:3)) = AE (1:3,1:3);
!!   AE12  (  IND(4:6),  IND(4:6)) = AE (1:3,1:3);
!!   AE12  (6+IND(1:3),6+IND(1:3)) = AE (1:3,1:3);
!!   AE12  (6+IND(4:6),6+IND(4:6)) = AE (1:3,1:3);
!!   AE12T                         = transpose(AE12);


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

      call SHAPEFUNC15 ( HTAL, ALLOC, PHPX, PHPZ, SHAPE, DSHAPE, NDFPE, NEQPE )

      UTT1     = matmul ( SHAPE , UT1 )
      UTT0     = matmul ( SHAPE , UT  )
      DUT0     = matmul ( DSHAPE, UT  ) !for TRACT force

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
!------ gravity
      XCM      = PA01 * beam_timo(n1)%XCMD_el   +  PA02 * beam_timo(n2)%XCMD_el
      ZCM      = PA01 * beam_timo(n1)%ZCMD_el   +  PA02 * beam_timo(n2)%ZCMD_el
      F_Grav   =-DENS*GRAV


!------ Inertial terms
!------ local Int {Nt*II*ρ*S*N}
      AM                 = 0.d0;
      AM (IND(1),IND(1)) = DENS
      AM (IND(1),IND(5)) = AMOMX
      AM (IND(2),IND(2)) = DENS
      AM (IND(2),IND(4)) =-AMOMX
      AM (IND(2),IND(6)) = AMOMZ
      AM (IND(3),IND(3)) = DENS
      AM (IND(3),IND(5)) =-AMOMZ
      AM (IND(4),IND(2)) =-AMOMX
      AM (IND(4),IND(4)) = RIXX
      AM (IND(4),IND(6)) =-RIXZ
      AM (IND(5),IND(1)) = AMOMX
      AM (IND(5),IND(3)) =-AMOMZ
      AM (IND(5),IND(5)) = POLI
      AM (IND(6),IND(2)) = AMOMZ
      AM (IND(6),IND(4)) =-RIXZ
      AM (IND(6),IND(6)) = RIZZ

      AQM = 0.d0;
      call IIxAxS           ( ACM, ATDAx2      , DENS, AMOMX, AMOMZ, RIXX, RIXZ, RIZZ,       1.d0 )
      call IIxAxS           ( AKM, ATDDA       , DENS, AMOMX, AMOMZ, RIXX, RIXZ, RIZZ,       1.d0 )
      call IIxAtDDR_AtDDAr0 ( AQM, ATDDA, ATDDR, DENS, AMOMX, AMOMZ, RIXX, RIXZ, RIZZ, HTA1, 1.d0 )


!------ Stiffness
!------ local Int {(N't*K1*N')+(N't*K2*N)+(Nt*K3*N')+(Nt*K4*N)}
!------ K1
      AK0        = 0.d0;
      AK0(1,1:6) = Kf( 1: 6)
      AK0(2,2:6) = Kf( 7:11)
      AK0(3,3:6) = Kf(12:15)
      AK0(4,4:6) = Kf(16:18)
      AK0(5,5:6) = Kf(19:20)
      AK0(6,  6) = Kf(   21)

      do ii      = 2,6
      do jj      = 1,ii-1
      AK0(ii,jj) = AK0(jj,ii)
      enddo!jj
      enddo!ii
      AK1(IND(1:6),IND(1:6)) = AK0(1:6,1:6);
!------ K2
      AK2(1:6   ,1:6   ) = 0.d0;
      AK2(1:6   ,IND(4)) =-AK1(1:6   ,IND(3))
      AK2(1:6   ,IND(6)) = AK1(1:6   ,IND(1))
!------ K3 [-]
      AK3(1:6   ,1:6   ) = 0.d0;
      AK3(IND(4),1:6   ) = AK1(IND(3),1:6   )
      AK3(IND(6),1:6   ) =-AK1(IND(1),1:6   )
!------ K4 [-]
      AK4                = 0.d0;
      AK4(1:6   ,IND(4)) =-AK3(1:6   ,IND(3))
      AK4(1:6   ,IND(6)) = AK3(1:6   ,IND(1))


!------ Nonlinear coupling terms (i.e. Bending-Tension) in K3 and K4
!     TRACT0               = dot_product(AK1(IND( 2),1:6),DUT0(1:6)) &
!                           +dot_product(AK2(IND( 2),1:6),UTT0(1:6))
   do ii                   = 1, 3
      F0(ii)               = dot_product(AK1(IND(ii),1:6),DUT0(1:6)) &
                            +dot_product(AK2(IND(ii),1:6),UTT0(1:6))
   enddo
      TRACT0               = F0(2)

      AK3NL                = 0.d0;
      AK4NL                = 0.d0;
      AQNL                 = 0.d0;
    if (nb_el<= NBODTWT_el.and.IEIG<1) then
      AK3NL(IND(4),1:6   ) =-AK1  (IND(2),1:6   ) * DUT0(IND(3)) !-K(2,:)wo'
      AK3NL(IND(6),1:6   ) = AK1  (IND(2),1:6   ) * DUT0(IND(1)) ! K(2,:)uo'
      AK3NL(IND(4),IND(3)) = AK3NL(IND(4),IND(3)) - TRACT0
      AK3NL(IND(6),IND(1)) = AK3NL(IND(6),IND(1)) + TRACT0
!------ Torsion nonlinear terms
!!    AK3NL(IND(5),1:6   ) = AK1  (IND(1),1:6   ) * DUT0(IND(3))&! K(1,:)wo'-K(3,:)uo'
!!                          -AK1  (IND(3),1:6   ) * DUT0(IND(1))
!!    AK3NL(IND(5),IND(3)) = AK3NL(IND(5),IND(3)) + F0(1)
!!    AK3NL(IND(5),IND(1)) = AK3NL(IND(5),IND(1)) - F0(3)

      AK4NL(IND(4),1:6   ) =-AK2  (IND(2),1:6   ) * DUT0(IND(3)) !-K(2,:)wo'
      AK4NL(IND(6),1:6   ) = AK2  (IND(2),1:6   ) * DUT0(IND(1)) ! K(2,:)uo'
!------ Torsion nonlinear terms
!!    AK4NL(IND(5),1:6   ) = AK2  (IND(1),1:6   ) * DUT0(IND(3))&! K(1,:)wo'-K(3,:)uo'
!!                          -AK2  (IND(3),1:6   ) * DUT0(IND(1))
    endif                
      AQNL(IND(4))         = TRACT0 * DUT0(IND(3))
      AQNL(IND(6))         =-TRACT0 * DUT0(IND(1))
!------ Torsion nonlinear terms
!!    AQNL(IND(5))         = F0(1)  * DUT0(IND(3)) &
!!                          -F0(3)  * DUT0(IND(1))


!------- Gravity [AQ= -II.At.F_Grav]
      AQG(1:6)             = 0.d0;

   if (int(Var_Paper(8))==1) &
      call IIxAtxF_Grav ( AQG, AT, F_Grav, XCM, ZCM )


!------ SF and Virtual Works
      AMLOC0               =   matmul( matmul(  SHAPET,  AM           ),  SHAPE )
      ACLOC0               =   matmul( matmul(  SHAPET,  ACM          ),  SHAPE )
      AKLOC0               =   matmul( matmul( DSHAPET,  AK1          )           &
                             +         matmul(  SHAPET, -AK3          ), DSHAPE ) &
                             + matmul( matmul( DSHAPET,  AK2          )           &
                             +         matmul(  SHAPET, -AK4+AKM      ),  SHAPE )
      AKNLLOC0             =   matmul( matmul(  SHAPET, -AK3NL        ), DSHAPE ) &
                             + matmul( matmul(  SHAPET, -AK4NL        ),  SHAPE )
      AFLOC0               =           matmul(  SHAPET, -AQG-AQNL-AQM )


!------ Gauss integration
      cgauss               = 0.5d0*ALLOC*WGAUSS(k)
      AMLOC                = AMLOC   + cgauss*AMLOC0
      ACLOC                = ACLOC   + cgauss*ACLOC0
      AKLOC                = AKLOC   + cgauss*AKLOC0
      AKNLLOC              = AKNLLOC + cgauss*AKNLLOC0
      AFLOC                = AFLOC   + cgauss*AFLOC0


!------ For all q's [Replace inside]
      call MATRLOC_Q ( AMLOCNQ , ACLOCNQ, AKLOCNQ, SHAPET, UTT0 , UTT1 , F_Grav, &
                       DENS    , AMOMX  , AMOMZ  , RIXX  , RIXZ , RIZZ , HTA1  , &
                       cgauss  , XCM    , ZCM    , nn_el , NDFPE, 1.d0             )
   enddo !k

!!!#ELEMENT
!!!----Rotate element matrices wrt element local c.s.
!!   AMLOC = matmul( matmul( AE12, AMLOC ) , AE12T );
!!   ACLOC = matmul( matmul( AE12, ACLOC ) , AE12T );
!!   AKLOC = matmul( matmul( AE12, AKLOC ) , AE12T );
!!   AFLOC =         matmul( AE12, AFLOC )          ;

!------ Linearization
   AFLOC = AFLOC - matmul(  AMLOC,  UT2 ) &
                 - matmul(  ACLOC,  UT1 ) &
                 - matmul(  AKLOC,  UT  )

   AKLOC = AKLOC + AKNLLOC;


 END Subroutine MATRLOC
!-------------------------------------------------------------------------------
 Subroutine MATRLOC_Q ( AMLOCNQ, ACLOCNQ, AKLOCNQ, SHAPET, UTT0 , UTT1 , F_Grav, &
                        DENS   , AMOMX  , AMOMZ  , RIXX  , RIXZ , RIZZ , HTA1  , &
                        cgauss , XCM    , ZCM    , nn_el , NDFPE, RAT              )
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
   real(8) :: AMLOCNQ (NDFPEM_el,NQS      )
   real(8) :: ACLOCNQ (NDFPEM_el,NQS      )
   real(8) :: AKLOCNQ (NDFPEM_el,NQS      )
   real(8) :: AMLOCNQ0(NDFPEM_el          )
   real(8) :: ACLOCNQ0(NDFPEM_el          )
   real(8) :: AKLOCNQ0(NDFPEM_el          )
   real(8) :: UTT0    (NEQPEM_el          )
   real(8) :: UTT1    (NEQPEM_el          )
   real(8) :: AM      (NEQPEM_el,NEQPEM_el)
   real(8) :: AQ      (NEQPEM_el          )
   real(8) :: SHAPET  (NDFPEM_el,NEQPEM_el)
   integer :: nn_el, n, NDFPE, iq
   real(8) :: HTA1, cgauss, DENS, AMOMX, AMOMZ, RIXX, RIZZ, RIXZ, F_Grav, XCM, ZCM, RAT


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

      AKLOCNQ0 =  matmul( SHAPET  , AQ )

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

      ACLOCNQ0 =  matmul( SHAPET  , AQ )

!---------------
!-- AMLOCNQ --
!---------------

!------ ATDDA2
      call IIxAxS (AM, ATDDA2, DENS, AMOMX, AMOMZ, RIXX, RIXZ, RIZZ, RAT)

      AQ = matmul ( AM, UTT0 )

!------ ATDDR2, ATDDA2
      call IIxAtDDR_AtDDAr0 ( AQ, ATDDA2, ATDDR2, DENS, AMOMX, AMOMZ, RIXX, RIXZ, RIZZ, HTA1, RAT )

      AMLOCNQ0 =  matmul( SHAPET  , AQ )

!------ Replace
      AKLOCNQ(1:NDFPE,n) = AKLOCNQ(1:NDFPE,n) + cgauss * AKLOCNQ0(1:NDFPE)
      ACLOCNQ(1:NDFPE,n) = ACLOCNQ(1:NDFPE,n) + cgauss * ACLOCNQ0(1:NDFPE)
      AMLOCNQ(1:NDFPE,n) = AMLOCNQ(1:NDFPE,n) + cgauss * AMLOCNQ0(1:NDFPE)
   enddo !iq


 END Subroutine MATRLOC_Q
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
 Subroutine FORCLOC (  AFLOC, AKLOCNQ,&
                       nbod , nbsub, nb_el, nn_el, nel )
!-------------------------------------------------------------------------

 Use Cbeam
 Use coupling_cfd_f

   implicit none

   real(8) :: AT      (3,3)
   real(8) :: AT0     (3,3)
   real(8) :: AFLOC   (NDFPEM_el)
   real(8) :: AKLOCNQ (NDFPEM_el,NQS)
   real(8) :: AKLOCNQ0(NDFPEM_el)
   real(8) :: AFLOC0  (NDFPEM_el)
   real(8) :: UT      (NDFPEM_el)
   real(8) :: UT1     (NDFPEM_el)
   real(8) :: UT2     (NDFPEM_el)
   real(8) :: UTT1    (NEQPEM_el)
   real(8) :: UTT0    (NEQPEM_el)
   real(8) :: DUT0    (NEQPEM_el)
   real(8) :: AQ      (NEQPEM_el)
   real(8) :: AIIA    (NEQPEM_el,6)
   real(8) :: FORC    (          6)
   real(8) :: F_ae    (          6)
   real(8) :: SHAPE   (NEQPEM_el,NDFPEM_el)
   real(8) :: DSHAPE  (NEQPEM_el,NDFPEM_el)
   real(8) :: SHAPET  (NDFPEM_el,NEQPEM_el)
   integer :: nbod, nbsub, nb_el, nn_el, nel, IB, NEQPE, NDFPE, k, n1, n2
   real(8) :: ALLOC, PA01 , PA02  , allocwgauss   , XC    , ZC    , F_Grav
   real(8) :: DENS ,                                                HTA   , HTAL
   real(8) ::                                                               XCM
   real(8) ::                               PHPX  , PHPZ  , HTA1  , HTA0  , ZCM, HTA1b
   real(8) :: DIAM
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

      call SHAPEFUNC15 ( HTAL, ALLOC, PHPX, PHPZ, SHAPE, DSHAPE, NDFPE, NEQPE )

      UTT1     = matmul ( SHAPE , UT1 )
      UTT0     = matmul ( SHAPE , UT  )
      DUT0     = matmul ( DSHAPE, UT  ) !for TRACT force

      SHAPET   = transpose (SHAPE )
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
      if (int(Var_Paper(10))==12                    .and.&
                         nel==subbody(nb_el)%NTEB_el.and.&
                       nbsub==   body(nbod )%NBODSUB_el    )then
!--------- flapwise force only in the last subbody and last element
         FX       = FX/ALLOC
         call IIxAtxF_GravX( AQ, AT, FX, XLoa, ZLoa, 1 )          !AQ is added
      endif
   endif


!------ Aerodynamic forces
      if (IB == 1.and.IAERO_el>0) then
        if     (IAERO_el == 1) then

           call LOCAL_FORC0_ae    ( nb_el, nbod, nbsub, HTA1, HTA1b, F_ae, XC, ZC )

        elseif (IAERO_el == 2) then

#         ifdef HAVE_GENUVP
          !call LOCAL_FORC0_gn_ae (        nbod, nbsub,          HTA1, F_ae         )
           call LOCAL_FORC0_gn_ae ( UTT0 , nbod, nbsub, nel , k, HTA1, F_ae, XC, ZC )
#         endif

          !XC = 0.d0 
          !ZC = 0.d0

        elseif (IAERO_el == 3) then

!--------- HAVE_GENUVP or HAVE_cfd

           F_ae(1:3) = matmul ( AT (1:3,1:3), F_cfd(nel  ,1:3,nb_el)*PA01 + &
                                             +F_cfd(nel+1,1:3,nb_el)*PA02     )
           F_ae(4:6) = matmul ( AT (1:3,1:3), F_cfd(nel  ,4:6,nb_el)*PA01 + &
                                             +F_cfd(nel+1,4:6,nb_el)*PA02     )
           XC        = 0.d0
           ZC        = 0.d0

           if (TIME.lt.0.02) &
               F_ae(1:6)= TIME/0.02 * F_ae(1:6)
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
      if (nbod == NTOWER_el) then
         DIAM  = PA01 * substr_mor(n1)%DIAMET_el  +  PA02 * substr_mor(n2)%DIAMET_el
         call VIV_tow  ( AQ, HTA1, DIAM, nn_el )
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
      AFLOC0  =    + matmul( SHAPET, AQ )

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
            if (int(Var_Paper(10))==12                    .and.&
                               nel==subbody(nb_el)%NTEB_el.and.&
                             nbsub==   body(nbod )%NBODSUB_el    )then
!--------------- flapwise force only in the last subbody and last element
               FX       = FX/ALLOC
               call IIxAtxF_GravX( AQ, AT0, FX, XLoa, ZLoa, 1 )          !AQ is added
            endif
         endif
         AKLOCNQ0 =  matmul( SHAPET  , AQ )
!--------- Replace
         AKLOCNQ(1:NDFPE,n) = AKLOCNQ(1:NDFPE,n) + allocwgauss * AKLOCNQ0(1:NDFPE)
      enddo !iq
   enddo !k


 END Subroutine FORCLOC
!-------------------------------------------------------------------------
!   Subroutine :CON_MATRLOC
! 
!  All Local Matrices [MASS, DAMPING, STIFFNESS, FORCES]
!  with concentrated masses calculated 
!  with respect to the body FEM local system.
!
!-------------------------------------------------------------------------
 Subroutine CON_MATRLOC( AMLOC, ACLOC, AKLOC, AFLOC, AMLOCNQ, ACLOCNQ, AKLOCNQ,&
                         ncal , nb_el, nn_el, nel                                )
!-------------------------------------------------------------------------

 Use Cbeam

   implicit none

   real(8) :: ATDAx2  (3,3)
   real(8) :: ATDDA   (3,3)
   real(8) :: ATDDR   (  3)
   real(8) :: AT      (3,3)
   real(8) :: AKLOC   (NDFPEM_el,NDFPEM_el)
   real(8) :: ACLOC   (NDFPEM_el,NDFPEM_el)
   real(8) :: AMLOC   (NDFPEM_el,NDFPEM_el)
   real(8) :: AFLOC   (NDFPEM_el          )
   real(8) :: AMLOCNQ (NDFPEM_el,NQS      )
   real(8) :: ACLOCNQ (NDFPEM_el,NQS      )
   real(8) :: AKLOCNQ (NDFPEM_el,NQS      )
   real(8) :: AKLOC0  (NDFPEM_el,NDFPEM_el)
   real(8) :: ACLOC0  (NDFPEM_el,NDFPEM_el)
   real(8) :: AMLOC0  (NDFPEM_el,NDFPEM_el)
   real(8) :: AFLOC0  (NDFPEM_el )         
   real(8) :: UT      (NDFPEM_el )         
   real(8) :: UT1     (NDFPEM_el )         
   real(8) :: UT2     (NDFPEM_el )         
   real(8) :: UTT1    (NEQPEM_el )         
   real(8) :: UTT0    (NEQPEM_el )         
   real(8) :: AM      (NEQPEM_el,NEQPEM_el)
   real(8) :: AQ      (NEQPEM_el          )
   real(8) :: SHAPE   (NEQPEM_el,NDFPEM_el)
   real(8) :: DSHAPE  (NEQPEM_el,NDFPEM_el)
   real(8) :: SHAPET  (NDFPEM_el,NEQPEM_el)
   real(8) :: DSHAPET (NDFPEM_el,NEQPEM_el)
   real(8) :: SHAPETAM(NDFPEM_el,NEQPEM_el)
   real(8) :: ALLOC, y0  , HTAL, Hta, PHPX, PHPZ, RAT
   real(8) :: Diag(3,3)
   real(8) :: Mass  , Sxoff , Syoff , Szoff
   real(8) :: IxxTot, IyyTot, IzzTot, IxyTot   , IxzTot, Iyztot
   real(8) :: Sy0   , Ixy0  , Iyz0  , Iyy0
   real(8) :: Xoff , Yoff, Zoff
   integer :: nb_el, nn_el, nel, NEQPE, NDFPE, ncal, nbsha_el


!--- Initialize local matrices
!  AMLOC   = 0.d0;
!  ACLOC   = 0.d0;
!  AKLOC   = 0.d0;
!  AFLOC   = 0.d0;
!  AMLOCNQ = 0.d0;
!  ACLOCNQ = 0.d0;
!  AKLOCNQ = 0.d0;

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


   call SHAPEFUNC15 ( HTAL, ALLOC, PHPX, PHPZ, SHAPE, DSHAPE, NDFPE, NEQPE )

   UTT1    = matmul ( SHAPE , UT1 )
   UTT0    = matmul ( SHAPE , UT  )

   SHAPET  = transpose (SHAPE )
   DSHAPET = transpose (DSHAPE)


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

   nbsha_el= 0
   if (NBODBTWT_el>NBLADE_el) &
   nbsha_el= body(NBLADE_el+1)%NBODTGLB_el (1)

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

   SHAPETAM = matmul( SHAPET  , AM    )
   AMLOC0   = matmul( SHAPETAM, SHAPE )
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

   SHAPETAM = matmul( SHAPET  , AM    )
   ACLOC0   = matmul( SHAPETAM, SHAPE )
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

   SHAPETAM = matmul( SHAPET , AM    )
   AKLOC0   = matmul( SHAPETAM, SHAPE )
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

   AFLOC0  = AFLOC0 + matmul( SHAPET, AQ )

!--- Replace
   AKLOC = AKLOC + AKLOC0 
   ACLOC = ACLOC + ACLOC0 
   AMLOC = AMLOC + AMLOC0 
   AFLOC = AFLOC - AFLOC0 

!----------------------------------------------
! QS
!----------------------------------------------

!---- For all q's
   call CON_MATRLOC_Q ( AMLOCNQ, ACLOCNQ, AKLOCNQ       ,&
                        SHAPET , UTT0   , UTT1          ,&
                        Mass   , Sxoff  , Syoff , Szoff ,&
                        IxxTot , IyyTot , IzzTot        ,&
                        IxyTot , IxzTot , Iyztot        ,&
                        Sy0    , Ixy0   , Iyz0  , Iyy0  ,&
                        nn_el  , NDFPE  , RAT              )


 END Subroutine CON_MATRLOC
!-------------------------------------------------------------------------------
 Subroutine CON_MATRLOC_Q ( AMLOCNQ, ACLOCNQ, AKLOCNQ       ,&
                            SHAPET , UTT0   , UTT1          ,&
                            Mass   , Sxoff  , Syoff , Szoff ,&
                            IxxTot , IyyTot , IzzTot        ,&
                            IxyTot , IxzTot , Iyztot        ,&
                            Sy0    , Ixy0   , Iyz0  , Iyy0  ,&
                            nn_el  , NDFPE  , RAT              )
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
   real(8) :: AMLOCNQ (NDFPEM_el,NQS      )
   real(8) :: ACLOCNQ (NDFPEM_el,NQS      )
   real(8) :: AKLOCNQ (NDFPEM_el,NQS      )
   real(8) :: AMLOCNQ0(NDFPEM_el          )
   real(8) :: ACLOCNQ0(NDFPEM_el          )
   real(8) :: AKLOCNQ0(NDFPEM_el          )
   real(8) :: UTT0    (NEQPEM_el          )
   real(8) :: UTT1    (NEQPEM_el          )
   real(8) :: AM      (NEQPEM_el,NEQPEM_el)
   real(8) :: AQ      (NEQPEM_el          )
   real(8) :: SHAPET  (NDFPEM_el,NEQPEM_el)
   integer :: nn_el, n, NDFPE, iq
   real(8) :: RAT
   real(8) :: Mass  , Sxoff , Syoff , Szoff
   real(8) :: IxxTot, IyyTot, IzzTot, IxyTot   , IxzTot, Iyztot
   real(8) :: Sy0   , Ixy0  , Iyz0  , Iyy0


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

      AKLOCNQ0 = matmul( SHAPET  , AQ )
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

      ACLOCNQ0 = matmul( SHAPET  , AQ )
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

      AMLOCNQ0 = matmul( SHAPET  , AQ )

!------ Replace
      AKLOCNQ(1:NDFPE,n) = AKLOCNQ(1:NDFPE,n) + AKLOCNQ0(1:NDFPE)
      ACLOCNQ(1:NDFPE,n) = ACLOCNQ(1:NDFPE,n) + ACLOCNQ0(1:NDFPE)
      AMLOCNQ(1:NDFPE,n) = AMLOCNQ(1:NDFPE,n) + AMLOCNQ0(1:NDFPE)
   enddo !iq


 END Subroutine CON_MATRLOC_Q
!---------------------------------------------------------------------------
!   Subroutine : LOCAL_FORC0_ae
!
!   Local loads communication from the aerodynamic
!-------------------------------------------------------------------------------
 Subroutine LOCAL_FORC0_ae ( nb_el, nbod, nbsub, HLOC, HLOCb, F_ae, XC_ae, ZC_ae )
!-------------------------------------------------------------------------------

 Use Cbeam
 Use Craft

   implicit none

   real(8) :: vec   (NSTRIPM)
   real(8) :: RL    (NSTRIPM)
   real(8) :: F     (  6)
   real(8) :: F_ae  (  6)
!! real(8) :: RB    (  3)
   real(8) :: AB    (3,3)
   real(8) :: ATB   (3,3)
   integer :: nb_el , nbod , nbsub, nboda, nn_el !, nn0_el
   real(8) :: HLOCb
   real(8) :: HLOC  , HLOC0, XC_ae, ZC_ae 
   integer :: i
!! real(8) :: HLOC00


   nn_el            = body(nbod)%NNODTGLB_el (nbsub)
   AB (1:3,1:3)     =         transf_mat(nn_el)%ATPA_el(1:3,1:3)   !from rotor-plane before pitch to blade c.s.
   ATB(1:3,1:3)     = transpose (AB(1:3,1:3))
   HLOC0            = HLOCb + Hhub !for new machi bend-sweep
   nboda            = subbody(nb_el)%IAENUM
   RL    (1:NSTRIP) = strip(nboda,1:NSTRIP)%rBld;


   if (HLOC0 < RL(1)) then
      F_ae(1:6) = 0.d0
      XC_ae     = 0.d0
      ZC_ae     = 0.d0
      return
   endif


                                      vec  (1:NSTRIP) = XCAER_el  (1,nboda,1:NSTRIP);
   call LIN_INT ( HLOC0, XC_ae  , RL, vec   , NSTRIP, NSTRIPM )
                                      vec  (1:NSTRIP) = XCAER_el  (3,nboda,1:NSTRIP);
   call LIN_INT ( HLOC0, ZC_ae  , RL, vec   , NSTRIP, NSTRIPM )

  do i = 1, 6
                                      vec  (1:NSTRIP) = FSTRIP_el (i,nboda,1:NSTRIP);
!!      LIN_INT (  x   ,  y     , xi, yi    , n1    , nnm     )
   call LIN_INT ( HLOC0,  F (i) , RL, vec   , NSTRIP, NSTRIPM )
  enddo !i

   F_ae(1:3) = matmul ( ATB(1:3,1:3), F(1:3)  )
   F_ae(4:6) = matmul ( ATB(1:3,1:3), F(4:6)  )
!old
!  F_ae(4:6) =                        F(4:6)


 END Subroutine LOCAL_FORC0_ae
!----------------------------------------------------------------------
 Subroutine LOADS_TRANS0_NEW ( B2B, A, J1, J2, JMAX )
!----------------------------------------------------------------------
! A*Floc = Fglob
! Floc   = At*Fglob
!----------------------------------------------------------------------
   implicit none

   integer :: J1, J2, JMAX
   real(8) :: B2B(6,1:JMAX), B2Bt(6,J1:J2), A (3,3)


   B2Bt(1  ,J1:J2) = B2B(1,J1:J2)*A(1,1) + B2B(3,J1:J2)*A(1,2) + B2B(4,J1:J2)*A(1,3)
   B2Bt(3  ,J1:J2) = B2B(1,J1:J2)*A(2,1) + B2B(3,J1:J2)*A(2,2) + B2B(4,J1:J2)*A(2,3)
   B2Bt(4  ,J1:J2) = B2B(1,J1:J2)*A(3,1) + B2B(3,J1:J2)*A(3,2) + B2B(4,J1:J2)*A(3,3)

   B2Bt(5  ,J1:J2) = B2B(5,J1:J2)*A(1,1) + B2B(6,J1:J2)*A(1,2) + B2B(2,J1:J2)*A(1,3)
   B2Bt(6  ,J1:J2) = B2B(5,J1:J2)*A(2,1) + B2B(6,J1:J2)*A(2,2) + B2B(2,J1:J2)*A(2,3)
   B2Bt(2  ,J1:J2) = B2B(5,J1:J2)*A(3,1) + B2B(6,J1:J2)*A(3,2) + B2B(2,J1:J2)*A(3,3)

   B2B (1:6,J1:J2) = B2Bt(1:6,J1:J2)

!  do j = J1, J2
!   B2B(IMDOF_el(1:3),j) = matmul ( B2B(IMDOF_el(1:3),j), A(1:3,1:3) ) 
!   B2B(IMDOF_el(4:6),j) = matmul ( B2B(IMDOF_el(4:6),j), A(1:3,1:3) ) 
!  enddo !j

 END Subroutine LOADS_TRANS0_NEW
!----------------------------------------------------------------------
 Subroutine LOADS_CALC_M0 ( B2B, R, J1, J2, JMAX )
!----------------------------------------------------------------------
   implicit none

   integer :: J1, J2, JMAX
   real(8) :: B2B(6,1:JMAX), R(3)


   B2B (5,J1:J2)   = B2B(5,J1:J2) + R(2)*B2B(4,J1:J2) - R(3)*B2B(3,J1:J2)
   B2B (6,J1:J2)   = B2B(6,J1:J2) + R(3)*B2B(1,J1:J2) - R(1)*B2B(4,J1:J2)
   B2B (2,J1:J2)   = B2B(2,J1:J2) + R(1)*B2B(3,J1:J2) - R(2)*B2B(1,J1:J2)
!  do j = J1, J2
!   F (1:3) = B2B(IMDOF_el(1:3),j)
!   call EXTEPR ( R, F, M )
!   B2B(IMDOF_el(4:6),j) = B2B(IMDOF_el(4:6),j) + M(1:3)
!  enddo !j


 END Subroutine LOADS_CALC_M0
!----------------------------------------------------------------------
 Subroutine LOADS_CALC_M00 ( B2B, R, J1, J2, JMAX )
!----------------------------------------------------------------------
   implicit none

   integer :: J1, J2, JMAX
   real(8) :: B2B(6,1:JMAX), R(3)


   B2B (5,J1:J2)   =                R(2)*B2B(4,J1:J2) - R(3)*B2B(3,J1:J2)
   B2B (6,J1:J2)   =                R(3)*B2B(1,J1:J2) - R(1)*B2B(4,J1:J2)
   B2B (2,J1:J2)   =                R(1)*B2B(3,J1:J2) - R(2)*B2B(1,J1:J2)


 END Subroutine LOADS_CALC_M00
!----------------------------------------------------------------------
 Subroutine LOCAL_UT_1 ( nb_el, nel, UT, UT1, UT2 )
!----------------------------------------------------------------------

 Use Cbeam
#ifdef HAVE_MODAL
 Use modal_mod
#endif

   implicit none

   integer, intent(in ) :: nb_el, nel
   real(8), intent(out) :: UT    (NDFPEM_el)
   real(8), intent(out) :: UT1   (NDFPEM_el)
   real(8), intent(out) :: UT2   (NDFPEM_el)

   integer :: nn, i, NDFPE


      NDFPE        = subbody(nb_el)%NDFPE_el
      nn           = subbody(nb_el)%NODMB(nel,1)                                  !nnod)
      i            = ACCU%NDFTACC_el(nb_el-1) + subbody(nb_el)%NDFPNBACC_el(nn-1) !+ndof

   if (IYNmodal_ut==0) then
      UT (1:NDFPE) = UT_el (i+1:i+NDFPE)
      UT1(1:NDFPE) = UT1_el(i+1:i+NDFPE)
      UT2(1:NDFPE) = UT2_el(i+1:i+NDFPE)
#  ifdef HAVE_MODAL
   else
      UT (1:NDFPE) = UT_md (i+1:i+NDFPE)
      UT1(1:NDFPE) = UT1_md(i+1:i+NDFPE)
      UT2(1:NDFPE) = UT2_md(i+1:i+NDFPE)
#  endif
   endif


 END Subroutine LOCAL_UT_1
!---------------------------------------------------------------------
!     Evaluation of the SHAPE FUNCTIONS
!
!
!---------------------------------------------------------------------
 Subroutine SHAPEFUNC15  (HTAL, ALLOC, PHPX, PHPZ, SHAPE, DSHAPE, NDFPE, NEQPE)
!----------------------------------------------------------------------

   implicit none

   integer, intent(in ) :: NDFPE , NEQPE
   real(8), intent(in ) :: HTAL  , ALLOC, PHPX  , PHPZ
   real(8), intent(out) :: SHAPE   (NEQPE,NDFPE)
   real(8), intent(out) :: DSHAPE  (NEQPE,NDFPE)


   if (NDFPE==12) call SF12 (1, HTAL, ALLOC, PHPX, PHPZ, SHAPE, DSHAPE, NDFPE, NEQPE)
   if (NDFPE==15) call SF15 (1, HTAL, ALLOC, PHPX, PHPZ, SHAPE, DSHAPE, NDFPE, NEQPE)


 END Subroutine SHAPEFUNC15
!---------------------------------------------------------------------
 Subroutine SSHAPEFUNC15 (HTAL, ALLOC, PHPX, PHPZ, SHAPE, NDFPE, NEQPE)
!---------------------------------------------------------------------

   implicit none

   integer, intent(in ) :: NDFPE , NEQPE
   real(8), intent(in ) :: HTAL  , ALLOC, PHPX  , PHPZ
   real(8), intent(out) :: SHAPE   (NEQPE,NDFPE)
   real(8)              :: DSHAPE  (NEQPE,NDFPE)


   if (NDFPE==12) call SF12 (0, HTAL, ALLOC, PHPX, PHPZ, SHAPE, DSHAPE, NDFPE, NEQPE)
   if (NDFPE==15) call SF15 (0, HTAL, ALLOC, PHPX, PHPZ, SHAPE, DSHAPE, NDFPE, NEQPE)


 END Subroutine SSHAPEFUNC15
!---------------------------------------------------------------------
 Subroutine SF12L (ider, HTAL, ALLOC, PHPX, PHPZ, SHAPE, DSHAPE, NDFPE, NEQPE)
!----------------------------------------------------------------------

   implicit none

   integer, intent(in ) :: NDFPE , NEQPE, ider
   real(8), intent(in ) :: HTAL  , ALLOC, PHPX , PHPZ
   real(8), intent(out) :: SHAPE  (NEQPE,NDFPE)
   real(8), intent(out) :: DSHAPE (NEQPE,NDFPE)

   real(8) :: S(2)


   S(1) = 1.d0-HTAL
   S(2) = HTAL

   SHAPE(:,:) = 0.d0;
   SHAPE(1,1) = S(1);   SHAPE(1, 7) = S(2);
   SHAPE(2,2) = S(1);   SHAPE(2, 8) = S(2);
   SHAPE(3,3) = S(1);   SHAPE(3, 9) = S(2);  
   SHAPE(4,4) = S(1);   SHAPE(4,10) = S(2);
   SHAPE(5,5) = S(1);   SHAPE(5,11) = S(2);
   SHAPE(6,6) = S(1);   SHAPE(6,12) = S(2);
   
   if (ider==0) return

! -- first order derivatives ---------------------------------

   S(1) = -1.d0/ALLOC
   S(2) =  1.d0/ALLOC

   DSHAPE(:,:) = 0.d0;
   DSHAPE(1,1) = S(1);   DSHAPE(1, 7) = S(2);
   DSHAPE(2,2) = S(1);   DSHAPE(2, 8) = S(2);
   DSHAPE(3,3) = S(1);   DSHAPE(3, 9) = S(2);  
   DSHAPE(4,4) = S(1);   DSHAPE(4,10) = S(2);
   DSHAPE(5,5) = S(1);   DSHAPE(5,11) = S(2);
   DSHAPE(6,6) = S(1);   DSHAPE(6,12) = S(2);


 END Subroutine SF12L
!---------------------------------------------------------------------
 Subroutine SF12 (ider, HTAL, ALLOC, PHPX, PHPZ, SHAPE, DSHAPE, NDFPE, NEQPE)
!----------------------------------------------------------------------

   implicit none

   integer, intent(in ) :: NDFPE , NEQPE, ider
   real(8), intent(in ) :: HTAL  , ALLOC, PHPX , PHPZ
   real(8), intent(out) :: SHAPE  (NEQPE,NDFPE)
   real(8), intent(out) :: DSHAPE (NEQPE,NDFPE)

   real(8) ::                             PHIX1, PHIX2 , PHIZ1 , PHIZ2
   real(8) :: PAX1 , PAX2 , PAX3 , PAX4 , PAX5 , PAX6  , PAX7  , PAX8
   real(8) :: PAZ1 , PAZ2 , PAZ3 , PAZ4 , PAZ5 , PAZ6  , PAZ7  , PAZ8
   real(8) :: PAV1 , PAV2 , PAT1 , PAT2
   real(8) :: DPAX1, DPAX2, DPAX3, DPAX4, DPAX5, DPAX6 , DPAX7 , DPAX8
   real(8) :: DPAZ1, DPAZ2, DPAZ3, DPAZ4, DPAZ5, DPAZ6 , DPAZ7 , DPAZ8
   real(8) :: DPAV1, DPAV2, DPAT1, DPAT2
   real(8) :: HLx2 , HLx3 , HLx4 , HLx45, HLx9 , HLxFX2, HLxFZ2, HLxL , HLx8
   real(8) :: PHIX1x05, PHIZ1x05, PHIX2x2, PHIZ2x2


!  PHIX1 = PHPX/(1.d0+PHPX)
!  PHIX2 = 1.d0/(1.d0+PHPX)
!  PHIZ1 = PHPZ/(1.d0+PHPZ)
!  PHIZ2 = 1.d0/(1.d0+PHPZ)
!  PAX1  =  1.0d0 - HTAL*PHIX1 - 3.0d0*HTAL**2*PHIX2 + 2.0d0*HTAL**3*PHIX2
!  PAX2  =          HTAL*PHIX1 + 3.0d0*HTAL**2*PHIX2 - 2.0d0*HTAL**3*PHIX2
!  PAX3  = (-HTAL +0.5d0*HTAL   *PHIX1  +0.5d0*HTAL**2*PHIX1 +2.0d0*HTAL**2*PHIX2-HTAL**3*PHIX2 )*ALLOC
!  PAX4  = (       0.5d0*HTAL   *PHIX1  -0.5d0*HTAL**2*PHIX1 +      HTAL**2*PHIX2-HTAL**3*PHIX2 )*ALLOC
!  PAX5  = ( 6.0d0*HTAL   *PHIX2 - 6.0d0*HTAL**2*PHIX2 )/ALLOC
!  PAX7  =  1.0d0 - HTAL * PHIX1 - 4.0d0*HTAL * PHIX2 +3.0d0*HTAL**2*PHIX2
!  PAX8  =          HTAL * PHIX1 - 2.0d0*HTAL * PHIX2 +3.0d0*HTAL**2*PHIX2

   PHIX1    = PHPX /(1.d0+PHPX)
   PHIX2    = PHIX1/PHPX
   PHIZ1    = PHPZ /(1.d0+PHPZ)
   PHIZ2    = PHIZ1/PHPZ
   PHIX1x05 = 0.5d0 * PHIX1
   PHIZ1x05 = 0.5d0 * PHIZ1
   PHIX2x2  = 2.0d0 * PHIX2
   PHIZ2x2  = 2.0d0 * PHIZ2
   HLx2     = 2.0d0 * HTAL
   HLx3     = 3.0d0 * HTAL
   HLx4     = 4.0d0 * HTAL
   HLx45    = 4.5d0 * HTAL
   HLx8     = 8.0d0 * HTAL
   HLx9     = 9.0d0 * HTAL
   HLxFX2   = PHIX2 * HTAL
   HLxFZ2   = PHIZ2 * HTAL
   HLxL     = ALLOC * HTAL

   PAX1 =  1d0 + HTAL * ( -PHIX1 + HLxFX2 * (-3d0 + HLx2) )
   PAX2 =  1d0 - PAX1
   PAX3 =  HLxL * (-1d0 + PHIX1x05*(1d0+HTAL) + HLxFX2*(2d0-HTAL))
   PAX4 =  HLxL * (1d0-HTAL)*(PHIX1x05 + HLxFX2)
   PAX5 =  6d0*HTAL*PHIX2/ALLOC * (1d0-HTAL)
   PAX6 = -PAX5
   PAX7 =  1d0 + HTAL * ( -PHIX1 + PHIX2 * (-4d0 + HLx3) )
   PAX8 =        HTAL * (  PHIX1 + PHIX2 * (-2d0 + HLx3) )
   PAZ1 =  1d0 + HTAL * ( -PHIZ1 + HLxFZ2 * (-3d0 + HLx2) )
   PAZ2 =  1d0 - PAZ1
   PAZ3 =  HLxL * (-1d0 + PHIZ1x05*(1d0+HTAL) + HLxFZ2*(2d0-HTAL))
   PAZ4 =  HLxL * (1d0-HTAL)*(PHIZ1x05 + HLxFZ2)
   PAZ5 =  6d0*HTAL*PHIZ2/ALLOC * (1d0-HTAL)
   PAZ6 = -PAZ5
   PAZ7 =  1d0 + HTAL * ( -PHIZ1 + PHIZ2 * (-4d0 + HLx3) )
   PAZ8 =        HTAL * (  PHIZ1 + PHIZ2 * (-2d0 + HLx3) )
   PAV1 = 1.d0-HTAL
   PAV2 = HTAL
   PAT1 = PAV1
   PAT2 = PAV2

   SHAPE (:,:) = 0.d0;
   SHAPE (1,1) = PAX1;   SHAPE (2,1) = PAX5;   SHAPE (3,3) = PAV1;   SHAPE (4, 4) = PAZ1;   SHAPE (5, 4) =-PAZ5;   SHAPE (6, 6) = PAT1;
   SHAPE (1,2) = PAX3;   SHAPE (2,2) = PAX7;   SHAPE (3,9) = PAV2;   SHAPE (4, 5) =-PAZ3;   SHAPE (5, 5) = PAZ7;   SHAPE (6,12) = PAT2;
   SHAPE (1,7) = PAX2;   SHAPE (2,7) = PAX6;                         SHAPE (4,10) = PAZ2;   SHAPE (5,10) =-PAZ6;
   SHAPE (1,8) = PAX4;   SHAPE (2,8) = PAX8;                         SHAPE (4,11) =-PAZ4;   SHAPE (5,11) = PAZ8;

   if (ider==0) return

! -- first order derivatives ---------------------------------

!  DPAX1 =(-              PHIX1 -6.0d0*HTAL   *PHIX2  +6.0d0*HTAL**2*PHIX2 )/ALLOC
!  DPAX2 = -DPAX1
!  DPAX3 = -1.0d0  +0.5d0*PHIX1  + HTAL*PHIX1 + 4.0d0*HTAL*PHIX2 - 3.0d0*HTAL**2*PHIX2 
!  DPAX4 =  0.5d0        *PHIX1  - HTAL*PHIX1 + 2.0d0*HTAL*PHIX2 - 3.0d0*HTAL**2*PHIX2 
!  DPAX5 =(  6.0d0       *PHIX2  -12.0d0*HTAL  *PHIX2 )/ALLOC**2
!  DPAX6 = -DPAX5
!  DPAX7 =(-              PHIX1  -4.0d0        *PHIX2 +6.0d0*HTAL   *PHIX2 )/ALLOC
!  DPAX8 =(               PHIX1  -2.0d0        *PHIX2 +6.0d0*HTAL   *PHIX2 )/ALLOC

   DPAX1 = -PAX5 - PHIX1/ALLOC
   DPAX2 = -DPAX1
   DPAX3 = -PAX7 + PHIX1x05
   DPAX4 = -PAX8 + PHIX1x05
   DPAX5 =  6d0*PHIX2/ALLOC**2 * (1d0 - HLx2)
   DPAX6 = -DPAX5
   DPAX7 =(-PHIX1 + PHIX2x2 *(HLx3-2d0) )/ALLOC
   DPAX8 =( PHIX1 + PHIX2x2 *(HLx3-1d0) )/ALLOC
   DPAZ1 = -PAZ5 - PHIZ1/ALLOC
   DPAZ2 = -DPAZ1
   DPAZ3 = -PAZ7 + PHIZ1x05
   DPAZ4 = -PAZ8 + PHIZ1x05
   DPAZ5 =  6d0*PHIZ2/ALLOC**2 * (1d0 - HLx2)
   DPAZ6 = -DPAZ5
   DPAZ7 =(-PHIZ1 + PHIZ2x2 *(HLx3-2d0) )/ALLOC
   DPAZ8 =( PHIZ1 + PHIZ2x2 *(HLx3-1d0) )/ALLOC
   DPAV1 = -1.d0/ALLOC
   DPAV2 =  1.d0/ALLOC
   DPAT1 = -1.d0/ALLOC
   DPAT2 =  1.d0/ALLOC

   DSHAPE(:,:) = 0.d0;
   DSHAPE(1,1) = DPAX1;  DSHAPE(2,1) = DPAX5;  DSHAPE(3,3) = DPAV1;  DSHAPE(4, 4) = DPAZ1;  DSHAPE(5, 4) =-DPAZ5;  DSHAPE(6, 6) = DPAT1;
   DSHAPE(1,2) = DPAX3;  DSHAPE(2,2) = DPAX7;  DSHAPE(3,9) = DPAV2;  DSHAPE(4, 5) =-DPAZ3;  DSHAPE(5, 5) = DPAZ7;  DSHAPE(6,12) = DPAT2;
   DSHAPE(1,7) = DPAX2;  DSHAPE(2,7) = DPAX6;                        DSHAPE(4,10) = DPAZ2;  DSHAPE(5,10) =-DPAZ6;
   DSHAPE(1,8) = DPAX4;  DSHAPE(2,8) = DPAX8;                        DSHAPE(4,11) =-DPAZ4;  DSHAPE(5,11) = DPAZ8;


 END Subroutine SF12
!---------------------------------------------------------------------
 Subroutine SF15 (ider, HTAL, ALLOC, PHPX, PHPZ, SHAPE, DSHAPE, NDFPE, NEQPE)
!----------------------------------------------------------------------

   implicit none

   integer, intent(in ) :: NDFPE , NEQPE, ider
   real(8), intent(in ) :: HTAL  , ALLOC, PHPX , PHPZ
   real(8), intent(out) :: SHAPE  (NEQPE,NDFPE)
   real(8), intent(out) :: DSHAPE (NEQPE,NDFPE)

   real(8) ::                             PHIX1, PHIX2 , PHIZ1 , PHIZ2
   real(8) :: PAX1 , PAX2 , PAX3 , PAX4 , PAX5 , PAX6  , PAX7  , PAX8
   real(8) :: PAZ1 , PAZ2 , PAZ3 , PAZ4 , PAZ5 , PAZ6  , PAZ7  , PAZ8
   real(8) :: PAV1 , PAV2 , PAV3 , PAV4 , PAT1 , PAT2  , PAT3
   real(8) :: DPAX1, DPAX2, DPAX3, DPAX4, DPAX5, DPAX6 , DPAX7 , DPAX8
   real(8) :: DPAZ1, DPAZ2, DPAZ3, DPAZ4, DPAZ5, DPAZ6 , DPAZ7 , DPAZ8
   real(8) :: DPAV1, DPAV2, DPAV3, DPAV4, DPAT1, DPAT2 , DPAT3
   real(8) :: HLx2 , HLx3 , HLx4 , HLx45, HLx9 , HLxFX2, HLxFZ2, HLxL , HLx8
   real(8) :: HTALm1, HTALm1xHTAL, HTALx2m1, HTALx2m1_2p1x9m10
   real(8) :: PHIX1x05, PHIZ1x05, PHIX2x2, PHIZ2x2


!  PHIX1 = PHPX/(1.d0+PHPX)
!  PHIX2 = 1.d0/(1.d0+PHPX)
!  PHIZ1 = PHPZ/(1.d0+PHPZ)
!  PHIZ2 = 1.d0/(1.d0+PHPZ)
!  PAX1  =  1.0d0 - HTAL*PHIX1 - 3.0d0*HTAL**2*PHIX2 + 2.0d0*HTAL**3*PHIX2
!  PAX2  =          HTAL*PHIX1 + 3.0d0*HTAL**2*PHIX2 - 2.0d0*HTAL**3*PHIX2
!  PAX3  = (-HTAL +0.5d0*HTAL   *PHIX1  +0.5d0*HTAL**2*PHIX1 +2.0d0*HTAL**2*PHIX2-HTAL**3*PHIX2 )*ALLOC
!  PAX4  = (       0.5d0*HTAL   *PHIX1  -0.5d0*HTAL**2*PHIX1 +      HTAL**2*PHIX2-HTAL**3*PHIX2 )*ALLOC
!  PAX5  = ( 6.0d0*HTAL   *PHIX2 - 6.0d0*HTAL**2*PHIX2 )/ALLOC
!  PAX7  =  1.0d0 - HTAL * PHIX1 - 4.0d0*HTAL * PHIX2 +3.0d0*HTAL**2*PHIX2
!  PAX8  =          HTAL * PHIX1 - 2.0d0*HTAL * PHIX2 +3.0d0*HTAL**2*PHIX2
!  PAV1 = 0.125d0*(1.d0-HTAL)*(-10.d0+9.d0*((2.d0*HTAL-1.d0)**2+1.))
!  PAV2 = 4.500d0*(1.d0-HTAL)*HTAL*( 2.d0-3.d0*HTAL)
!  PAV3 = 4.500d0*(1.d0-HTAL)*HTAL*(-1.d0+3.d0*HTAL)
!  PAV4 = 0.125d0*      HTAL *(-10.d0+9.d0*((2.d0*HTAL-1.d0)**2+1.))
!  PAT1 =-     (1.d0-HTAL)     *(2.d0*HTAL-1.d0)
!  PAT2 = 4.d0*(1.d0-HTAL)*HTAL
!  PAT3 =                  HTAL*(2.d0*HTAL-1.d0)

   PHIX1    = PHPX /(1.d0+PHPX)
   PHIX2    = PHIX1/PHPX
   PHIZ1    = PHPZ /(1.d0+PHPZ)
   PHIZ2    = PHIZ1/PHPZ
   PHIX1x05 = 0.5d0 * PHIX1
   PHIZ1x05 = 0.5d0 * PHIZ1
   PHIX2x2  = 2d0   * PHIX2
   PHIZ2x2  = 2d0   * PHIZ2
   HLx2     = 2d0   * HTAL
   HLx3     = 3d0   * HTAL
   HLx4     = 4d0   * HTAL
   HLx45    = 4.5d0 * HTAL
   HLx8     = 8d0   * HTAL
   HLx9     = 9d0   * HTAL
   HLxFX2   = PHIX2 * HTAL
   HLxFZ2   = PHIZ2 * HTAL
   HLxL     = ALLOC * HTAL

   PAX1 =  1d0 + HTAL * ( -PHIX1 + HLxFX2 * (-3d0 + HLx2) )
   PAX2 =  1d0 - PAX1
   PAX3 =  HLxL * (-1d0 + PHIX1x05*(1d0+HTAL) + HLxFX2*(2d0-HTAL))
   PAX4 =  HLxL * (1d0-HTAL)*(PHIX1x05 + HLxFX2)
   PAX5 =  6d0*HTAL*PHIX2/ALLOC * (1d0-HTAL)
   PAX6 = -PAX5
   PAX7 =  1d0 + HTAL * ( -PHIX1 + PHIX2 * (-4d0 + HLx3) )
   PAX8 =        HTAL * (  PHIX1 + PHIX2 * (-2d0 + HLx3) )
   PAZ1 =  1d0 + HTAL * ( -PHIZ1 + HLxFZ2 * (-3d0 + HLx2) )
   PAZ2 =  1d0 - PAZ1
   PAZ3 =  HLxL * (-1d0 + PHIZ1x05*(1d0+HTAL) + HLxFZ2*(2d0-HTAL))
   PAZ4 =  HLxL * (1d0-HTAL)*(PHIZ1x05 + HLxFZ2)
   PAZ5 =  6d0*HTAL*PHIZ2/ALLOC * (1d0-HTAL)
   PAZ6 = -PAZ5
   PAZ7 =  1d0 + HTAL * ( -PHIZ1 + PHIZ2 * (-4d0 + HLx3) )
   PAZ8 =        HTAL * (  PHIZ1 + PHIZ2 * (-2d0 + HLx3) )

   HTALm1       = 1d0 - HTAL
   HTALm1xHTAL  = HTALm1 * HTAL
   HTALx2m1     = HLx2-1d0
   HTALx2m1_2p1x9m10 = 0.125d0*( (HTALx2m1**2 +1.d0)*9d0 - 10d0 )
   
   PAV1 =       HTALm1     * HTALx2m1_2p1x9m10    !0.125d0*(-10.d0+9.d0*HTALx2m1_2p1)
   PAV2 = 4.5d0*HTALm1xHTAL*(  2.d0-HLx3  )
!  PAV3 = 4.5d0*HTALm1xHTAL*(- 1.d0+HLx3  )
   PAV3 = 4.5d0*HTALm1xHTAL - PAV2
   PAV4 =              HTAL* HTALx2m1_2p1x9m10    !0.125d0*(-10.d0+9.d0*HTALx2m1_2p1)
   PAT1 =-     HTALm1     *HTALx2m1
   PAT2 = 4.d0*HTALm1xHTAL
   PAT3 =                  HTAL*HTALx2m1

   SHAPE(:, :)  = 0.d0;
   SHAPE(1, 1)  = PAX1;   SHAPE(2, 1)  = PAX5;   SHAPE(3, 3)  = PAV1;   SHAPE(4, 4)  = PAZ1;   SHAPE(5, 4)  =-PAZ5;   SHAPE(6, 6)  = PAT1;
   SHAPE(1, 2)  = PAX3;   SHAPE(2, 2)  = PAX7;   SHAPE(3, 7)  = PAV2;   SHAPE(4, 5)  =-PAZ3;   SHAPE(5, 5)  = PAZ7;   SHAPE(6, 8)  = PAT2;
   SHAPE(1,10)  = PAX2;   SHAPE(2,10)  = PAX6;   SHAPE(3, 9)  = PAV3;   SHAPE(4,13)  = PAZ2;   SHAPE(5,13)  =-PAZ6;   SHAPE(6,15)  = PAT3;
   SHAPE(1,11)  = PAX4;   SHAPE(2,11)  = PAX8;   SHAPE(3,12)  = PAV4;   SHAPE(4,14)  =-PAZ4;   SHAPE(5,14)  = PAZ8;  

   if (ider==0) return

! -- first order derivatives ---------------------------------

!  DPAX1 =(-              PHIX1 -6.0d0*HTAL   *PHIX2  +6.0d0*HTAL**2*PHIX2 )/ALLOC
!  DPAX2 = -DPAX1
!  DPAX3 = -1.0d0  +0.5d0*PHIX1  + HTAL*PHIX1 + 4.0d0*HTAL*PHIX2 - 3.0d0*HTAL**2*PHIX2 
!  DPAX4 =  0.5d0        *PHIX1  - HTAL*PHIX1 + 2.0d0*HTAL*PHIX2 - 3.0d0*HTAL**2*PHIX2 
!  DPAX5 =(  6.0d0       *PHIX2  -12.0d0*HTAL  *PHIX2 )/ALLOC**2
!  DPAX6 = -DPAX5
!  DPAX7 =(-              PHIX1  -4.0d0        *PHIX2 +6.0d0*HTAL   *PHIX2 )/ALLOC
!  DPAX8 =(               PHIX1  -2.0d0        *PHIX2 +6.0d0*HTAL   *PHIX2 )/ALLOC
!  DPAV1 = ( -5.5d0 + 18d0*HTAL - 13.5d0*HTAL**2 ) / ALLOC
!  DPAV2 = (  9.0d0 - 45d0*HTAL + 40.5d0*HTAL**2 ) / ALLOC
!  DPAV3 = ( -4.5d0 + 36d0*HTAL - 40.5d0*HTAL**2 ) / ALLOC
!  DPAV4 = (  1.0d0 -  9d0*HTAL + 13.5d0*HTAL**2 ) / ALLOC

   DPAX1 = -PAX5 - PHIX1/ALLOC
   DPAX2 = -DPAX1
   DPAX3 = -PAX7 + PHIX1x05
   DPAX4 = -PAX8 + PHIX1x05
   DPAX5 =  6d0*PHIX2/ALLOC**2 * (1d0 - HLx2)
   DPAX6 = -DPAX5
   DPAX7 =(-PHIX1 + PHIX2x2 *(HLx3-2d0) )/ALLOC
   DPAX8 =( PHIX1 + PHIX2x2 *(HLx3-1d0) )/ALLOC
   DPAZ1 = -PAZ5 - PHIZ1/ALLOC
   DPAZ2 = -DPAZ1
   DPAZ3 = -PAZ7 + PHIZ1x05
   DPAZ4 = -PAZ8 + PHIZ1x05
   DPAZ5 =  6d0*PHIZ2/ALLOC**2 * (1d0 - HLx2)
   DPAZ6 = -DPAZ5
   DPAZ7 =(-PHIZ1 + PHIZ2x2 *(HLx3-2d0) )/ALLOC
   DPAZ8 =( PHIZ1 + PHIZ2x2 *(HLx3-1d0) )/ALLOC

   DPAV1 = ( -5.5d0 + HLx45 * ( 4d0 - HLx3) ) / ALLOC
   DPAV2 = (  9.0d0 - HLx45 * (10d0 - HLx9) ) / ALLOC
   DPAV3 = ( -4.5d0 + HLx45 * ( 8d0 - HLx9) ) / ALLOC
   DPAV4 = (  1.0d0 - HLx45 * ( 2d0 - HLx3) ) / ALLOC
   DPAT1 = ( -3d0 + HLx4 ) / ALLOC
   DPAT2 = (  4d0 - HLx8 ) / ALLOC
   DPAT3 = ( -1d0 + HLx4 ) / ALLOC

   DSHAPE(:, :) = 0.d0;
   DSHAPE(1, 1) = DPAX1;  DSHAPE(2, 1) = DPAX5;  DSHAPE(3, 3) = DPAV1;  DSHAPE(4, 4) = DPAZ1;  DSHAPE(5, 4) =-DPAZ5;  DSHAPE(6, 6) = DPAT1;
   DSHAPE(1, 2) = DPAX3;  DSHAPE(2, 2) = DPAX7;  DSHAPE(3, 7) = DPAV2;  DSHAPE(4, 5) =-DPAZ3;  DSHAPE(5, 5) = DPAZ7;  DSHAPE(6, 8) = DPAT2;
   DSHAPE(1,10) = DPAX2;  DSHAPE(2,10) = DPAX6;  DSHAPE(3, 9) = DPAV3;  DSHAPE(4,13) = DPAZ2;  DSHAPE(5,13) =-DPAZ6;  DSHAPE(6,15) = DPAT3;
   DSHAPE(1,11) = DPAX4;  DSHAPE(2,11) = DPAX8;  DSHAPE(3,12) = DPAV4;  DSHAPE(4,14) =-DPAZ4;  DSHAPE(5,14) = DPAZ8;


 END Subroutine SF15
!--------------------------------------------------------------------
 Subroutine  IIxAxS ( AM, A, DENS, AMOMX, AMOMZ, RIXX, RIXZ, RIZZ, RAT )
!--------------------------------------------------------------------

   implicit none

   real(8) :: AM(6,6), A(3,3), DENS, AMOMX, AMOMZ, RIXX, RIXZ, RIZZ, RAT


   AM(1,1) =  DENS *A(1,1)
   AM(3,1) =  DENS *A(2,1)
   AM(4,1) =  DENS *A(3,1)
   AM(5,1) = -AMOMX*A(2,1)
   AM(6,1) =  AMOMX*A(1,1) - AMOMZ*A(3,1)
   AM(2,1) =  AMOMZ*A(2,1)

   AM(1,3) =  DENS *A(1,2)
   AM(3,3) =  DENS *A(2,2)
   AM(4,3) =  DENS *A(3,2)
   AM(5,3) = -AMOMX*A(2,2)
   AM(6,3) =  AMOMX*A(1,2) - AMOMZ*A(3,2)
   AM(2,3) =  AMOMZ*A(2,2)

   AM(1,4) =  DENS *A(1,3)
   AM(3,4) =  DENS *A(2,3)
   AM(4,4) =  DENS *A(3,3)
   AM(5,4) = -AMOMX*A(2,3)
   AM(6,4) =  AMOMX*A(1,3) - AMOMZ*A(3,3)
   AM(2,4) =  AMOMZ*A(2,3)

   AM(1,5) = -AMOMX*A(1,2)
   AM(3,5) = -AMOMX*A(2,2)
   AM(4,5) = -AMOMX*A(3,2)
   AM(5,5) =  RIXX* A(2,2)                           !-massYZ*[A(2,3)+A(3,2)]+massYY*A(3,3)
   AM(6,5) =(-RIXX *A(1,2) + RIXZ *A(3,2))  *  RAT   !+massYZ* A(1,3)        -massXY*A(3,3)
   AM(2,5) = -RIXZ *A(2,2)                           !+massYZ* A(1,2)-massYY*A(1,3)+massXY*A(2,3)

   AM(1,6) = -AMOMZ*A(1,3) + AMOMX*A(1,1)
   AM(3,6) = -AMOMZ*A(2,3) + AMOMX*A(2,1)
   AM(4,6) = -AMOMZ*A(3,3) + AMOMX*A(3,1)
   AM(5,6) =  RIXZ *A(2,3) - RIXX *A(2,1)            !massYZ*A(3,1)-massXY*A(3,3)
   AM(6,6) = (RIZZ *A(3,3) - RIXZ *( A(1,3)+A(3,1) ) + RIXX*A(1,1))  *  RAT
   AM(2,6) = -RIZZ *A(2,3) + RIXZ *A(2,1)            !-massYZ*A(1,1)+massXY*A(1,3)

   AM(1,2) =  AMOMZ*A(1,2)
   AM(3,2) =  AMOMZ*A(2,2)
   AM(4,2) =  AMOMZ*A(3,2)
   AM(5,2) = -RIXZ *A(2,2)                           ! massYZ*A(2,1)-massYY*A(3,1)+massXY*A(3,2)
   AM(6,2) =(-RIZZ *A(3,2) + RIXZ *A(1,2))  *  RAT   !-massYZ*A(1,1)+massXY*A(3,1)
   AM(2,2) =  RIZZ *A(2,2)                           ! massYY*A(1,1)-massXY*[A(1,2)+A(2,1)]


 END Subroutine  IIxAxS
!--------------------------------------------------------------------
 Subroutine IIxAtDDR_AtDDAr0 ( AQ, A, R, DENS, AMOMX, AMOMZ, RIXX, RIXZ, RIZZ, HTA1, RAT )
!--------------------------------------------------------------------

   implicit none

   real(8) :: AQ(6), A(3,3), R(3), DENS, AMOMX, AMOMZ, RIXX, RIXZ, RIZZ, HTA1, RAT


   AQ(1) = AQ(1) + AMOMZ*A(1,1) + DENS*( HTA1*A(1,2)+R(1) ) + AMOMX*  A(1,3)
   AQ(3) = AQ(3) + AMOMZ*A(2,1) + DENS*( HTA1*A(2,2)+R(2) ) + AMOMX*  A(2,3)
   AQ(4) = AQ(4) + AMOMZ*A(3,1) + DENS*( HTA1*A(3,2)+R(3) ) + AMOMX*  A(3,3)
   AQ(5) = AQ(5) - RIXZ *A(2,1) - RIXX       *A(2,3)        - AMOMX*( A(2,2)*HTA1 + R(2) )
   AQ(6) = AQ(6) +(RIXZ * ( A(1,1) - A(3,3) ) + RIXX*A(1,3) - RIZZ*A(3,1) ) * RAT  &
                 + AMOMX*( HTA1*A(1,2)+R(1) ) - AMOMZ*( HTA1*A(3,2)+R(3) )
   AQ(2) = AQ(2) + RIZZ *A(2,1) + RIXZ       *A(2,3)        + AMOMZ*( A(2,2)*HTA1 + R(2) )


 END Subroutine  IIxAtDDR_AtDDAr0 
!--------------------------------------------------------------------
 Subroutine IIxAtxF_Grav ( AQ, A, F_Grav, XCM, ZCM )
!--------------------------------------------------------------------

   implicit none

   real(8) :: AQ(6), A(3,3), F_Grav, XCM, ZCM


   AQ(1) = -F_Grav * A(1,3)
   AQ(3) = -F_Grav * A(2,3)
   AQ(4) = -F_Grav * A(3,3)
   AQ(5) =  F_Grav * A(2,3)*ZCM
   AQ(6) = -F_Grav *(A(1,3)*ZCM - A(3,3)*XCM)
   AQ(2) = -F_Grav * A(2,3)*XCM


 END Subroutine IIxAtxF_Grav
!--------------------------------------------------------------------
 Subroutine IIxAtxF_GravX ( AQ, A, F_Grav, XCM, ZCM, IAX )
!--------------------------------------------------------------------

!jim checks for the paper
   implicit none


   real(8) :: AQ(6), A(3,3), F_Grav, XCM, ZCM
   integer :: IAX !1:x, 2:z


   AQ(1) = AQ(1) + F_Grav * A(1,IAX)
   AQ(3) = AQ(3) + F_Grav * A(2,IAX)
   AQ(4) = AQ(4) + F_Grav * A(3,IAX)
   AQ(5) = AQ(5) - F_Grav * A(2,IAX)*ZCM
   AQ(6) = AQ(6) + F_Grav *(A(1,IAX)*ZCM - A(3,IAX)*XCM)
   AQ(2) = AQ(2) + F_Grav * A(2,IAX)*XCM


 END Subroutine IIxAtxF_GravX
!--------------------------------------------------------------------
 Subroutine  Con_IIxAxS ( AM    , A     , RAT          ,&
                          m     , Sxoff , Syoff , Szoff,&
                          IxxTot, IyyTot, IzzTot       ,&
                          IxyTot, IxzTot, Iyztot       ,&
                          IND                             )
!--------------------------------------------------------------------

   implicit none

   real(8) :: AM(6,6), A(3,3), RAT
   real(8) :: m      , Sxoff , Syoff , Szoff
   real(8) :: IxxTot , IyyTot, IzzTot, IxyTot, IxzTot, Iyztot
   integer :: IND(6)


   AM(IND(1:3),IND(1:3)) =  m     *A(1:3,1:3)
   AM(IND(1:3),IND(4  )) = -Szoff *A(1:3,2  ) + Syoff * A(1:3,3  )
   AM(IND(1:3),IND(5  )) =  Szoff *A(1:3,1  ) - Sxoff * A(1:3,3  )
   AM(IND(1:3),IND(6  )) = -Syoff *A(1:3,1  ) + Sxoff * A(1:3,2  )

   AM(IND(4  ),IND(1:3)) = -Szoff *A(2  ,1:3) + Syoff * A(3  ,1:3)
   AM(IND(4  ),IND(4  )) =  IzzTot*A(2  ,2  ) - IyzTot*(A(2  ,3  ) +        A(3,2))+ IyyTot*A(3,3)
   AM(IND(4  ),IND(5  )) = -IzzTot*A(2  ,1  ) + IxzTot* A(2  ,3  ) + IyzTot*A(3,1) - IxyTot*A(3,3)
   AM(IND(4  ),IND(6  )) =  IyzTot*A(2  ,1  ) - IxzTot* A(2  ,2  ) - IyyTot*A(3,1) + IxyTot*A(3,2)

   AM(IND(5  ),IND(1:3)) =  Szoff *A(1  ,1:3) - Sxoff * A(3  ,1:3)
   AM(IND(5  ),IND(4  )) =(-IzzTot*A(1  ,2  ) + IyzTot* A(1  ,3  ) + IxzTot*A(3,2) - IxyTot*A(3,3)) * RAT
   AM(IND(5  ),IND(5  )) =( IzzTot*A(1  ,1  ) - IxzTot*(A(1  ,3  ) +        A(3,1))+ IxxTot*A(3,3)) * RAT
   AM(IND(5  ),IND(6  )) =(-IyzTot*A(1  ,1  ) + IxzTot* A(1  ,2  ) + IxyTot*A(3,1) - IxxTot*A(3,2)) * RAT

   AM(IND(6  ),IND(1:3)) = -Syoff *A(1  ,1:3) + Sxoff * A(2  ,1:3)
   AM(IND(6  ),IND(4  )) =  IyzTot*A(1  ,2  ) - IyyTot* A(1  ,3  ) - IxzTot*A(2,2) + IxyTot*A(2,3)
   AM(IND(6  ),IND(5  )) = -IyzTot*A(1  ,1  ) + IxyTot* A(1  ,3  ) + IxzTot*A(2,1) - IxxTot*A(2,3)
   AM(IND(6  ),IND(6  )) =  IyyTot*A(1  ,1  ) - IxyTot*(A(1  ,2  ) +        A(2,1))+ IxxTot*A(2,2)


 END Subroutine  Con_IIxAxS
!--------------------------------------------------------------------
 Subroutine Con_IIxAtDDR_AtDDAr0 ( AQ    , A     , R     , RAT   ,&
                                   m     , Sxoff , Syoff , Szoff ,&
                                   IxxTot, IyyTot, IzzTot        ,&
                                   IxyTot, IxzTot, Iyztot        ,&
                                   Sy0   , Ixy0  , Iyz0  , Iyy0  ,&
                                   IND                              )
!--------------------------------------------------------------------

   implicit none

   real(8) :: AQ (6), A(3,3), R(3)  , RAT
   real(8) :: m     , Sxoff , Syoff , Szoff
   real(8) :: IxxTot, IyyTot, IzzTot, IxyTot, IxzTot, Iyztot
   real(8) :: Sy0   , Ixy0  , Iyz0  , Iyy0
   integer :: IND(6)


   AQ(IND(1:3)) = AQ(IND(1:3)) + Sxoff *A(1:3,1) + (Sy0  + Syoff )*A(1:3,2) + Szoff *A(1:3,3) + m*R(1:3)
   AQ(IND(4  )) = AQ(IND(4  )) - IxzTot*A(2  ,1) - (Iyz0 + IyzTot)*A(2  ,2) - IzzTot*A(2  ,3) + IxyTot*A(3,1) + (Iyy0 + IyyTot)*A(3,2) + IyzTot*A(3,3)        - Szoff*R(2) + Syoff*R(3)
   AQ(IND(5  )) = AQ(IND(5  )) +(IxzTot*A(1  ,1) + (Iyz0 + IyzTot)*A(1  ,2) + IzzTot*A(1  ,3) - IxxTot*A(3,1) - (Ixy0 + IxyTot)*A(3,2) - IxzTot*A(3,3)) * RAT + Szoff*R(1) - Sxoff*R(3)
   AQ(IND(6  )) = AQ(IND(6  )) - IxyTot*A(1  ,1) - (Iyy0 + IyyTot)*A(1  ,2) - IyzTot*A(1  ,3) + IxxTot*A(2,1) + (Ixy0 + IxyTot)*A(2,2) + IxzTot*A(2,3)        - Syoff*R(1) + Sxoff*R(2)


 END Subroutine  Con_IIxAtDDR_AtDDAr0 
!--------------------------------------------------------------------
 Subroutine Con_IIxAtxF_Grav (AQ, A, m, Sxoff, Syoff, Szoff, g, IND)
!--------------------------------------------------------------------

   implicit none

   real(8) :: AQ (6), A(3,3), m, Sxoff, Syoff, Szoff, g
   integer :: IND(6)

   AQ(IND(1:3)) =  g*m*A(1:3,3)                         !-- (-) at LHS, at the end it will be subtracted
   AQ(IND(4  )) = -g*( Szoff*A(2,3) - Syoff*A(3,3))
   AQ(IND(5  )) = -g*(-Szoff*A(1,3) + Sxoff*A(3,3))
   AQ(IND(6  )) = -g*( Syoff*A(1,3) - Sxoff*A(2,3))


 END Subroutine Con_IIxAtxF_Grav
!-------------------------------------------------------------------------
 Subroutine Calc_SS ( SS, X0 )
!-------------------------------------------------------------------------

 implicit none

   real(8) :: X0(3), SS(3,6)


   SS(:,:) = 0.d0;
   SS(1,1) = 1.d0
   SS(2,2) = 1.d0
   SS(3,3) = 1.d0
   SS(2,4) = -X0(3)
   SS(1,5) =  X0(3)
   SS(3,5) = -X0(1)
   SS(2,6) =  X0(1)


 END Subroutine Calc_SS
!
!
!---------------------------------------------------------------
!
!-- Output position, velocity and acceleration wrt local and global c.s.
!
!-- example of call for blade tip deflections
!do nbod     = 1, NBLADE_el
!   nbsub    = body   (nbod )%NBODSUB_el
!   nb_el    = body   (nbod )%NBODTGLB_el (nbsub)
!   nel      = subbody(nb_el)%NTEB_el
!   AKSI0(:) = 0.d0;
!   AKSI0(2) = subbody(nb_el)%HTA_el  (nel+1)
!   call beam_kinematics (nbod,nbsub,nel,AKSI0,x,dx,ddx,dxl,ddxl)
!enddo
!--------------------------------------------------------------
 Subroutine beam_kinematics (nbod,nbsub,nel,AKSI0,x,dx,ddx,dxl,ddxl)
!--------------------------------------------------------------

 Use Cbeam

   implicit none

   integer,               intent(in ) :: nbod, nbsub, nel
   real(8), dimension(3), intent(in ) :: AKSI0
   real(8), dimension(3), intent(out) :: x, dx, ddx, dxl,ddxl

   real(8), dimension(NDFPEM_el) :: UT  , UT1 , UT2
   real(8), dimension(NEQPEM_el) :: UTT0, UTT1, UTT2
   real(8), dimension(3,3      ) :: A   , ATDA, ATDAx2, ATDDA
   real(8), dimension(3        ) :: R   , ATDR,         ATDDR
   real(8), dimension(3        ) :: SSU0, SSU1, SSU2
   real(8)                       :: SHAPE(NEQPEM_el, NDFPEM_el)
   real(8)                       :: SS   (3        , NEQPEM_el)

   real(8) :: ALLOC, HTA0 , HTA  , HTAL , PHPX, PHPZ
   integer :: nb_el, nn_el, NEQPE, NDFPE


   nb_el            = body   (nbod )%NBODTGLB_el(nbsub)
   nn_el            = body   (nbod )%NNODTGLB_el(nbsub)
   ALLOC            = subbody(nb_el)%ALENG_el   (nel)
   HTA0             = subbody(nb_el)%HTA_el     (nel)
   PHPX             = subbody(nb_el)%PHIX_el    (nel)
   PHPZ             = subbody(nb_el)%PHIZ_el    (nel)
   NEQPE            = subbody(nb_el)%NEQPE_el
   NDFPE            = subbody(nb_el)%NDFPE_el
   HTA              = AKSI0 (2) - HTA0
   HTAL             = HTA/ALLOC

!--- check
!  if (HTAL>1.d0.or.HTAL<0.d0) then
!     write(*,*) 'error in beam_kinematics. check AKSI0'
!     stop
!  endif
   
   call SSHAPEFUNC15 (HTAL , ALLOC, PHPX, PHPZ, SHAPE, NDFPE, NEQPE)
   call LOCAL_UT_1   (nb_el, nel  , UT  , UT1 , UT2  )

   UTT0             = matmul  ( SHAPE, UT  )
   UTT1             = matmul  ( SHAPE, UT1 )
   UTT2             = matmul  ( SHAPE, UT2 )

   A      (1:3,1:3) = transf_mat(nn_el)%A_el     (1:3,1:3)
   ATDA   (1:3,1:3) = transf_mat(nn_el)%ATDA_el  (1:3,1:3)
   ATDAx2 (1:3,1:3) = transf_mat(nn_el)%ATDA_el  (1:3,1:3) * 2.d0
   ATDDA  (1:3,1:3) = transf_mat(nn_el)%ATDDA_el (1:3,1:3)
   R      (1:3    ) = transf_mat(nn_el)%R_el     (1:3    )
   ATDR   (1:3    ) = transf_mat(nn_el)%ATDR_el  (1:3    )
   ATDDR  (1:3    ) = transf_mat(nn_el)%ATDDR_el (1:3    )
      
   call Calc_SS ( SS, AKSI0 )

   SSU0   (1:3)     = matmul( SS(1:3,1:NEQPE), UTT0(IMDOF_el(1:NEQPE)) ) + AKSI0(1:3)
   SSU1   (1:3)     = matmul( SS(1:3,1:NEQPE), UTT1(IMDOF_el(1:NEQPE)) )
   SSU2   (1:3)     = matmul( SS(1:3,1:NEQPE), UTT2(IMDOF_el(1:NEQPE)) )

!--- local acceleration, velocity
   ddxl   (1:3)     = ATDDR(1:3) + matmul ( ATDDA (1:3,1:3), SSU0 (1:3) ) + &
                                   matmul ( ATDAx2(1:3,1:3), SSU1 (1:3) ) + SSU2(1:3)
   dxl    (1:3)     = ATDR (1:3) + matmul ( ATDA  (1:3,1:3), SSU0 (1:3) ) + SSU1(1:3)
!--- global acceleration, velocity, position
   ddx    (1:3)     =              matmul ( A     (1:3,1:3), ddxl (1:3) )
   dx     (1:3)     =              matmul ( A     (1:3,1:3), dxl  (1:3) )
   x      (1:3)     = R    (1:3) + matmul ( A     (1:3,1:3), SSU0 (1:3) )


 END Subroutine beam_kinematics
!----------------------------------------------------------------------
 Subroutine beam_loads (nbod,nbsub,nel,nod, FGG, FLG, FLL)
!----------------------------------------------------------------------
!
!   Provides the beam loads at node (nbod, nbsub,nel,nod)
!   FGG: Global       c.s
!   FLG: Local Globar c.s. [see below]
!   FLL: Local        c.s. [see below]
!
!           LG             LL        
!           Global         Local      
!   blade : rotor disk   - local-root 
!   shaft : non-rotating - rotating   
!   tower : global c.s.  - local      
!   jacket: global c.s.  - local      
!                                     
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   integer, intent (in   ) :: nbod  , nbsub , nel   , nod
   real(8), intent (  out) :: FGG(6), FLG(6), FLL(6)

   real(8)      :: AG     (3,3), AL     (3,3), AGG     (3,3)
   real(8)      :: AT0G   (3,3), AT0L   (3,3)
   real(8)      :: AFB2BG (6  ), AFB2BL (6  ), AFB2BGG (6  )
   real(8)      :: dir
   integer      :: I1, I2
   integer      :: nn_el , nb_el , nn0_el


   nb_el          = body(nbod)%NBODTGLB_el(nbsub)
   nn_el          = body(nbod)%NNODTGLB_el(nbsub)

   if     (nod==1      ) then; dir = 1.d0           !1. for nod=1, -1 for nod=NNPE_el
   elseif (nod==NNPE_el) then; dir =-1.d0           !1. for nod=1, -1 for nod=NNPE_el
   else                      ; stop "wrong nod input in beam_loads";
   endif

!------ Define global and local c.s. per body
   if     (nbod <= NBLADE_el) then
!------ wrt blade root after cone and pitch, before pre-curved angle
      AT0G(1:3,1:3) = ATP_el  (nbod,1:3,1:3) !Blade c.s. after cone, before pitch 
      AT0L(1:3,1:3) = ATPWr_el(nbod,1:3,1:3) !Blade c.s. after cone, and pitch, before pre-curved angles
!------ wrt global c.s. for the VAWT
      if (IAPPL_el == 1) &
      AT0G(1:3,1:3) = ATsftNoR_el (1:3,1:3)
!!    call DIAGO (3, AT0G)
   elseif (nbod == NSHAFT_el) then
!------ wrt the non-rotating shaft c.s.
      AT0G(1:3,1:3) = ATsftNoR_el (1:3,1:3)
      nn0_el        = body(nbod)%NNODTGLB_el(1)
      AT0L(1:3,1:3) = transf_mat(nn0_el)%AT_el(1:3,1:3)
   else
!------ wrt global c.s. for jacket and tower
      call DIAGO (3, AT0G)
      nn0_el        = body(nbod)%NNODTGLB_el(1)
      AT0L(1:3,1:3) = transf_mat(nn0_el)%AT_el(1:3,1:3)
   endif


   AGG(1:3,1:3)   =                        transf_mat(nn_el)%A_el(1:3,1:3)  !G  -> FGG
   AG (1:3,1:3)   = matmul (AT0G (1:3,1:3),transf_mat(nn_el)%A_el(1:3,1:3)) !LG -> FLG
   AL (1:3,1:3)   = matmul (AT0L (1:3,1:3),transf_mat(nn_el)%A_el(1:3,1:3)) !LL -> FLL

   I1             = subbody(nb_el)%NDFPNACC_el ( nod-1      ) + 1
   I2             = subbody(nb_el)%NDFPNACC_el ( nod        )
   AFB2BGG( 1:6 ) = subbody(nb_el)%AFLOC_el    ( nel, I1:I2 )*dir;
   AFB2BG ( 1:6 ) = subbody(nb_el)%AFLOC_el    ( nel, I1:I2 )*dir;
   AFB2BL ( 1:6 ) = subbody(nb_el)%AFLOC_el    ( nel, I1:I2 )*dir;

!--- Rotate LOADS with respect to body c.s. (0)
   call LOADS_TRANS0_NEW ( AFB2BGG, AGG, 1, 1, 1 )
   call LOADS_TRANS0_NEW ( AFB2BG , AG , 1, 1, 1 )
   call LOADS_TRANS0_NEW ( AFB2BL , AL , 1, 1, 1 )

   FGG    ( 1:6 ) = AFB2BGG(IMDOF_el(1:6))
   FLG    ( 1:6 ) = AFB2BG (IMDOF_el(1:6))
   FLL    ( 1:6 ) = AFB2BL (IMDOF_el(1:6))


 END Subroutine beam_loads
