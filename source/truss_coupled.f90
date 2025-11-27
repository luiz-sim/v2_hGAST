!----------------------------------------------------------------------
 Subroutine MOORINGS_COUPLED_tr ( ISTEP, ii )
!----------------------------------------------------------------------

 Use Cbeam
 Use truss

   implicit none

   integer :: ISTEP, i, I1, I2, ii, ndf, IP


   if ( ITRUSS_el /= 2) return

   goto (1,2,3,4), ISTEP


!--- Initialize the coupled truss code ---
!-   increase the global Matrices [+NDFT_tr] and set (again) initial values
 1 continue

   NDFBT_el = NDFBT_el + NDFT_tr
   NDFT_el  = NDFT_el  + NDFT_tr

   write(* ,*)'***Truss Coupled***'
   write(* ,*)'NDFT_tr ', NDFT_tr
   write(* ,*)'NDFT_el ', NDFT_el
   write(* ,*)'NDFBT_el', NDFBT_el
   write(* ,*)
   write(10,*)'***Truss Coupled***'
   write(10,*)'NDFT_tr ', NDFT_tr
   write(10,*)'NDFT_el ', NDFT_el
   write(10,*)'NDFBT_el', NDFBT_el
   write(10,*)


   Deallocate ( AM_el , AC_el  , AK_el  , AQ_el, INDSYSB,&
                UT_el , UT1_el , UT2_el                 ,&
                UTP_el, UTP1_el, UTP2_el                ,&
                UT0_el, UT01_el, UT02_el                   )

   Allocate   ( AM_el   (NDFT_el,NDFT_el) ,&
                AC_el   (NDFT_el,NDFT_el) ,&
                AK_el   (NDFT_el,NDFT_el) ,&
                AQ_el   (NDFT_el        ) ,&
                INDSYSB (NDFT_el        ) ,&
                UT_el   (NDFT_el), UT1_el (NDFT_el), UT2_el (NDFT_el) ,&
                UTP_el  (NDFT_el), UTP1_el(NDFT_el), UTP2_el(NDFT_el) ,&
                UT0_el  (NDFT_el), UT01_el(NDFT_el), UT02_el(NDFT_el)    )


!- Set the initial values for basic parameters
   call INIT_DEFLECTIONS ( 0.d0 ) !omega_set
 return


!--- Calculate truss Matrices, assembly to global ---
!-      Communicate loads and translations 
 2 continue

   IP = ii

   I1 = ACCU%NDFTACC_el (NBODT_el) + 1
   I2 = NDFBT_el

   UT_tr (1:NDFT_tr) = UT_el (I1:I2);
   UT1_tr(1:NDFT_tr) = UT1_el(I1:I2);
   UT2_tr(1:NDFT_tr) = UT2_el(I1:I2);

   call MATRIX_tr
   call BOUNDC_tr(IP)

   AM_el(I1:I2,I1:I2) = AM_tr (1:NDFT_tr,1:NDFT_tr)
   AC_el(I1:I2,I1:I2) = AC_tr (1:NDFT_tr,1:NDFT_tr)
   AK_el(I1:I2,I1:I2) = AK_tr (1:NDFT_tr,1:NDFT_tr)
   AQ_el(I1:I2      ) = AQ_tr (1:NDFT_tr          )

!- Loads from moorings to the floater
   call truss2float_LOADS
!- Positions or accelerations from floater to connecting points
   call set_connect_tr(IP)         !apply different b.c. at the same lines where the uncoupled code would set the known connection values.
 return


!--- Global matrix Reduction ---
!--  keep lines
 3 continue

   if (IMODAL > 0) return

   do ndf = 1, NDFT_tr
! -For all truss dofs
      ii = ii + 1
      i  = ACCU%NDFTACC_el(NBODT_el) + ndf

      INDSYSB(ii) = i

      AM_el ( ii, 1:NDFT_el ) = AM_el ( i, 1:NDFT_el )
      AC_el ( ii, 1:NDFT_el ) = AC_el ( i, 1:NDFT_el )
      AK_el ( ii, 1:NDFT_el ) = AK_el ( i, 1:NDFT_el )
      AQ_el ( ii            ) = AQ_el ( i            )
   enddo
 return

!--  keep Columns
 4 continue

   if (IMODAL > 0) return

   do ndf = 1, NDFT_tr
! -For all truss dofs
      ii = ii + 1
      i  = ACCU%NDFTACC_el(NBODT_el) + ndf
     
      AM_el ( 1:NDFT0_el, ii ) = AM_el ( 1:NDFT0_el, i )
      AC_el ( 1:NDFT0_el, ii ) = AC_el ( 1:NDFT0_el, i )
      AK_el ( 1:NDFT0_el, ii ) = AK_el ( 1:NDFT0_el, i )
   enddo
 return


 END Subroutine MOORINGS_COUPLED_tr
!----------------------------------------------------------------------
 Subroutine set_connect_tr (IP)
!----------------------------------------------------------------------

 Use truss
 Use Cbeam


   implicit none
  
   real(8) :: RG (3), RG0, DDRG(3), DDRG0,DDRG1,DDRG2
   integer :: IP, nc, nb_tr, ni, inod, i, iq, nq, itr, nrot, irot


   do nc = 1, NCONNECT_tr

      nb_tr   = connect_tr(nc)%IBODCONNECT_tr
      ni      = connect_tr(nc)%INODCONNECT_tr
      inod    = body_tr(nb_tr)%INODE_tr(ni)
      i       = ACCU%NDFTACC_el(NBODT_el) + NDFACC_tr (inod-1)


      AM_el ( i+1:i+3, 1:NDFT_el ) = 0.d0
      AC_el ( i+1:i+3, 1:NDFT_el ) = 0.d0
      AK_el ( i+1:i+3, 1:NDFT_el ) = 0.d0


      if (IP.eq.0) then

         RG(1:3) = matmul (A_float(1:3,1:3), connect_tr(nc)%RLOC_CONNECT_tr(1:3) )

         do nq=1,3

            itr   = i        + nq
            iq    = NDFBT_el + nq

            AK_el(itr, itr ) = 1.d0
            AK_el(itr, iq  ) =-1.d0
            AQ_el(itr      ) = UT_el(iq) - UT_el(itr) + RG(nq) - XG_tr(nq, inod)

            do nrot=1,3
               RG0 = dot_product (A0_float(nrot,nq,1:3), connect_tr(nc)%RLOC_CONNECT_tr(1:3) )
               irot = NDFBT_el+3+nrot
               AK_el( itr, irot ) = -RG0
            enddo

         enddo !nq

      else

         DDRG(1:3) = matmul (DDA_float(1:3,1:3), connect_tr(nc)%RLOC_CONNECT_tr(1:3) )

         do nq=1,3

            itr   = i        + nq
            iq    = NDFBT_el + nq

            AM_el(itr, itr ) = 1.d0
            AM_el(itr, iq  ) =-1.d0
            AQ_el(itr      ) = UT2_el(iq)-UT2_el(itr) + DDRG(nq)

            do nrot=1,3

               irot = NDFBT_el+3+nrot

               DDRG0 = dot_product (DDA0_float(nrot,nq,1:3), connect_tr(nc)%RLOC_CONNECT_tr(1:3) )      !DDA0
               DDRG1 = dot_product ( DA0_float(nrot,nq,1:3), connect_tr(nc)%RLOC_CONNECT_tr(1:3) )*2.d0 !DDA1=2DA0
               DDRG2 = dot_product (  A0_float(nrot,nq,1:3), connect_tr(nc)%RLOC_CONNECT_tr(1:3) )      !DDA2=A0

               AK_el( itr, irot ) = -DDRG0
               AC_el( itr, irot ) = -DDRG1
               AM_el( itr, irot ) = -DDRG2

            enddo

         enddo !nq

      endif

   enddo  !nc


 END Subroutine set_connect_tr
!----------------------------------------------------------------------
 Subroutine truss2float_LOADS
!----------------------------------------------------------------------

 Use truss
 Use Cbeam


   implicit none
  
   real(8) :: AMB2B    (6,NDFPE_tr) ,&
              ACB2B    (6,NDFPE_tr) ,&
              AKB2B    (6,NDFPE_tr) ,&
              AFB2B    (6         ) ,&
              AKB2BNQ  (6,  3     ) ,&
              RG       (3         )
   integer :: nb_tr, ni, nc, nq, i,j, inod, ii, j1


   do nc = 1, NCONNECT_tr


      call LOC_LOADS_tr ( AMB2B, ACB2B, AKB2B, AFB2B, nc ) !loads
      AKB2BNQ(:,:) = 0.d0;


      RG(1:3) = matmul (A_float(1:3,1:3), connect_tr(nc)%RLOC_CONNECT_tr(1:3) )

      call LOADS_CALC_M0_tr ( AMB2B  , RG, 1 , NDFPE_tr, NDFPE_tr )! RG x M
      call LOADS_CALC_M0_tr ( ACB2B  , RG, 1 , NDFPE_tr, NDFPE_tr )! RG x C
      call LOADS_CALC_M0_tr ( AKB2B  , RG, 1 , NDFPE_tr, NDFPE_tr )! RG x K
      call LOADS_CALC_M0_tr ( AFB2B  , RG, 1 , 1       , 1        )! RG x Q

      do nq=1,3

         AKB2BNQ(1:3,nq) = AFB2B(1:3)    !only to pass AQ to sub-routine
         RG     (1:3   ) = matmul ( A0_float(nq,1:3,1:3), connect_tr(nc)%RLOC_CONNECT_tr(1:3) ) !RG0

         call LOADS_CALC_M0_tr ( AKB2BNQ, RG, nq, nq      , 3        )! R0G_nq x Q

         AKB2BNQ(1:3,nq) = 0.d0;         !zero the force terms and only keep the moment

!     write(*,*)'wr0',nq,AKB2BNQ ( 1:6, nq )
      enddo


!-assembly to global matrices
      nb_tr = connect_tr(nc)%IBODCONNECT_tr
      ni    = connect_tr(nc)%INODCONNECT_tr
      inod  = body_tr(nb_tr)%INODE_tr (ni)

      j     = ACCU%NDFTACC_el(NBODT_el) + NDFACC_tr (inod-1)


      do ii = 1,6     !6 floating equations
         i  = NDFBT_el + ii

         AM_el ( i, j+1:j+NDFPE_tr ) = AM_el ( i, j+1:j+NDFPE_tr ) + AMB2B ( ii, 1:NDFPE_tr )
         AC_el ( i, j+1:j+NDFPE_tr ) = AC_el ( i, j+1:j+NDFPE_tr ) + ACB2B ( ii, 1:NDFPE_tr )
         AK_el ( i, j+1:j+NDFPE_tr ) = AK_el ( i, j+1:j+NDFPE_tr ) + AKB2B ( ii, 1:NDFPE_tr )
         AQ_el ( i                 ) = AQ_el ( i                 ) + AFB2B ( ii             )

         do nq = 1, 3    !3 floating rotations

            j1 = NDFBT_el + 3 + nq

            AK_el ( i, j1 ) = AK_el ( i, j1 ) - AKB2BNQ ( ii, nq )

         enddo   !nq
      enddo   !ii

   enddo   !nc


 END Subroutine truss2float_LOADS
!----------------------------------------------------------------------
 Subroutine LOC_LOADS_tr ( AMB2B, ACB2B, AKB2B, AFB2B, nc )
!----------------------------------------------------------------------

 Use truss

   implicit none
  
   real (8) :: AMB2B    (6,NDFPE_tr) ,&
               ACB2B    (6,NDFPE_tr) ,&
               AKB2B    (6,NDFPE_tr) ,&
               AFB2B    (6         )
   integer  :: nb_tr, ni, I1, I2, nc


  AMB2B = 0.d0;
  ACB2B = 0.d0;
  AKB2B = 0.d0;
  AFB2B = 0.d0;

! Set the loads of connection
  nb_tr = connect_tr(nc)%IBODCONNECT_tr 
  ni    = connect_tr(nc)%INODCONNECT_tr
   
  I1    = (ni-1)*NDFPN_tr + 1
  I2    = (ni-1)*NDFPN_tr + 3

  AMB2B  ( 1:3, 1:NDFPE_tr ) = body_tr(nb_tr)%AMLOC_tr ( I1:I2, 1:NDFPE_tr )
  ACB2B  ( 1:3, 1:NDFPE_tr ) = body_tr(nb_tr)%ACLOC_tr ( I1:I2, 1:NDFPE_tr )
  AKB2B  ( 1:3, 1:NDFPE_tr ) = body_tr(nb_tr)%AKLOC_tr ( I1:I2, 1:NDFPE_tr )
  AFB2B  ( 1:3             ) = body_tr(nb_tr)%AQLOC_tr ( I1:I2             )


 END Subroutine LOC_LOADS_tr
!----------------------------------------------------------------------
 Subroutine LOADS_CALC_M0_tr ( B2B, R, J1, J2, JMAX )
!----------------------------------------------------------------------
!
! A*Floc = Fglob
! Floc   = At*Fglob
!
!----------------------------------------------------------------------
   implicit none

   integer :: J1, J2, JMAX
   real(8) :: B2B(6,JMAX), R(3)


   B2B(4,J1:J2) = R(2)*B2B(3,J1:J2) - R(3)*B2B(2,J1:J2)
   B2B(5,J1:J2) = R(3)*B2B(1,J1:J2) - R(1)*B2B(3,J1:J2)
   B2B(6,J1:J2) = R(1)*B2B(2,J1:J2) - R(2)*B2B(1,J1:J2)


 END Subroutine LOADS_CALC_M0_tr
!
!
!
!-------------------------------------------------------------------------------------------------
!
! -- Subroutine :EIGEN_FLOATER
!
!------------------------------------------------------------------------------------------------- 
 Subroutine EIGEN_FLOATER
!------------------------------------------------------------------------------------------------- 

 Use Cbeam

   implicit none

   real(8), allocatable  :: AM_INV (:,:)               !(NDFT0_el,NDFT0_el)
   real(8), allocatable  :: FTI_el (:  ),&
                            FTR1_el(:  ),&
                            FTI1_el(:  )               !(NDFT0_el)
   integer, allocatable  :: INDXE  (:  )               !(NDFT0_el)
   real(8), allocatable  :: ALFRS  (:  ),&
                            ALFIS  (:  ),&
                            RFR    (:  ),&
                            RFI    (:  )               !(2*NDFT0_el)
   integer, allocatable  :: IND    (:  )               !(2*NDFT0_el)
   real(8), allocatable  :: RZR    (:,:),&
                            RZI    (:,:),&
                            RZR1   (:,:),&
                            RZI1   (:,:),&
                            VS     (:,:),&
                            AAS    (:,:)               !(2*NDFT0_el,2*NDFT0_el)
   real(8), allocatable  :: WORK   (:  )               !(LWORK)
   real(8)               :: VL(1)
   integer               :: NDFT2_el, LWORK
   real(8)               :: A1,C1, RATIO, ZCRIT
   integer               :: i, j, m, NUMF, NUMFT, ierr, NOUT1, NDFTR_el
   character*80          :: outfil


   NDFT0_el = 6
   NDFT2_el = 2*NDFT0_el
   LWORK    = 8*NDFT2_el+16

   write (*,*) 'EIGEN NEW is called',NDFT0_el,NDFT2_el

   Allocate  (  AM_INV    (NDFT0_el, NDFT0_el),&
                FTI_el    (NDFT0_el          ),&
                FTR1_el   (NDFT0_el          ),&
                FTI1_el   (NDFT0_el          ),&
                INDXE     (NDFT0_el          )   )

   Allocate  (  ALFRS     (NDFT2_el          ),&
                ALFIS     (NDFT2_el          ),&
                RFR       (NDFT2_el          ),&
                RFI       (NDFT2_el          ),&
                IND       (NDFT2_el          )   )

   Allocate  (  RZR       (NDFT2_el, NDFT2_el),&
                RZI       (NDFT2_el, NDFT2_el),&
                RZR1      (NDFT2_el, NDFT2_el),&
                RZI1      (NDFT2_el, NDFT2_el),&
                VS        (NDFT2_el, NDFT2_el),&
                AAS       (NDFT2_el, NDFT2_el)   )

   Allocate  (  WORK      (LWORK             )   )

   NOUT1 = 98


   open (NOUT1, file='eigen.dat'      )

   AM_INV(1:NDFT0_el,1:NDFT0_el) = AM_el(1+NDFBT0_el:NDFT0_el+NDFBT0_el,1+NDFBT0_el:NDFT0_el+NDFBT0_el);

   call DGETRF ( NDFT0_el, NDFT0_el, AM_INV  , NDFT0_el, INDXE,        ierr )
   call DGETRI ( NDFT0_el, AM_INV  , NDFT0_el, INDXE   , WORK , LWORK, ierr )

   AAS(1:NDFT2_el,1:NDFT2_el) = 0.d0;

   
   do i=1,NDFT0_el

      do j=1,NDFT0_el

         if (i.eq.j) AAS(i,j+NDFT0_el) = 1.d0

         do m=1,NDFT0_el
            AAS(i+NDFT0_el,j         )= AAS(i+NDFT0_el,j         ) - AM_INV(i,m)*AK_el(m+NDFBT0_el,j+NDFBT0_el)
            AAS(i+NDFT0_el,j+NDFT0_el)= AAS(i+NDFT0_el,j+NDFT0_el) - AM_INV(i,m)*AC_el(m+NDFBT0_el,j+NDFBT0_el)
         enddo

      enddo
   enddo


!-- call subroutine which solves the eigen-value problem

   call dgeev ( 'N', 'V', NDFT2_el, AAS, NDFT2_el, ALFRS, ALFIS, VL, 1, VS, &
                          NDFT2_el, WORK, LWORK, ierr )


   NUMF = 0
   do I = 1, NDFT2_el

      A1    = ALFRS(I)
      C1    = ALFIS(I)
      
      if ( C1 < 1d-10 ) cycle   !to avoid the symmetric negative eigen-values

      NUMF = NUMF + 1
      RFR (NUMF) = A1/PI2
      RFI (NUMF) = C1/PI2
    
      do J=1,NDFT0_el 
         RZR  (J,NUMF) = VS(J,I+1)
         RZI  (J,NUMF) = VS(J,I  )
      enddo

   enddo

   NUMFT    = NUMF
   NDFTR_el = 6 !NUMFT  
   

!-- sort
   call SORTRX (RFI, NDFT2_el, NUMFT, IND)


!- Zero out elastic deflections and set initial pitch
   UT_el   = 0.d0;
   UT1_el  = 0.d0;
   UT2_el  = 0.d0;

!  do nbod=1,NBLADE_el
!     NQACC = ACCU%NQSACC_el(nbod-1) 
!     UT_el (NDFBT_el+NQSP+NQACC+5) = PITCH0_el(nbod)
!  enddo


   do I=1,NDFTR_el

!-- Substitute the modal - solution to all active dofs
!     do mm=1,NDFT0_el
!        m = INDSYSB (mm)
!        UT_el(m) = RZR(mm,IND(I))
!     enddo
!     
!     call GET_REDUSED_DOFS


!--  Scale Modes to make deformations visible!
!     UT_el = UT_el * 10.d0; !1500.d0;

!     call ROTMAT_el
      
!-- global deformations
      outfil =  'modeR'
!     call WRITEOUT_WT  (outfil,I,3)
!-- local deformations
      outfil =  'mode'
!     call WRITEOUTEIG  (outfil,I,3)
!!    call WRITEOUT     (outfil, outfil1)

      if (dabs(RFR(IND(I)))>=1d-17) then
         RATIO =  RFI(IND(I))/RFR(IND(I))
         ZCRIT = 1.d0/(dsqrt(1.d0+RATIO**2))
      else
!        RATIO = 1.d10
         ZCRIT = 0.d0
      endif

      write (NOUT1,1000) I,IND(I),RFR(IND(I)),RFI(IND(I)), ZCRIT !, RATIO

   enddo !I

   close (NOUT1)


 100  format (12F15.8)
1000  format (I5,I7,4(F15.5,5x))

   Deallocate  ( AM_INV, FTI_el, FTR1_el, FTI1_el, INDXE      )
   Deallocate  ( ALFRS , ALFIS , RFR    , RFI    , IND        )
   Deallocate  ( RZR   , RZI   , RZR1   , RZI1   , VS   , AAS )
   Deallocate  ( WORK                                         )


 END Subroutine EIGEN_FLOATER
