!-------------------------------------------------------------------------------------------------
!
! -- Subroutine :EIGEN
!
!------------------------------------------------------------------------------------------------- 
 Subroutine EIGEN
!------------------------------------------------------------------------------------------------- 

 Use Cbeam
#ifdef HAVE_MODAL
 Use modal_mod
#endif

   implicit None

   real(8), Allocatable :: AM_INV (:,:)               !(NDFT0_el,NDFT0_el)
   real(8), Allocatable :: FTI_el (:  )
   real(8), Allocatable :: FTR1_el(:  )
   real(8), Allocatable :: FTI1_el(:  )               !(NDFT0_el)
   integer, Allocatable :: INDXE  (:  )               !(NDFT0_el)
   real(8), Allocatable :: ALFRS  (:  )
   real(8), Allocatable :: ALFIS  (:  )
   real(8), Allocatable :: RFR    (:  )
   real(8), Allocatable :: RFI    (:  )               !(2*NDFT0_el)
   integer, Allocatable :: IND    (:  )               !(2*NDFT0_el)
   integer, Allocatable :: IND2   (:  )               !(2*NDFT0_el)
   real(8), Allocatable :: RZR    (:,:)
   real(8), Allocatable :: RZI    (:,:)
   real(8), Allocatable :: RZR1   (:,:)
   real(8), Allocatable :: RZI1   (:,:)
   real(8), Allocatable :: VS     (:,:)
   real(8), Allocatable :: AAS    (:,:)               !(2*NDFT0_el,2*NDFT0_el)
   real(8), Allocatable :: WORK   (:  )               !(LWORK)
   real(8)              :: VL(1)
   integer              :: NDFT2_el, LWORK
   real(8)              :: A1,C1, RATIO, ZCRIT
   integer              :: i, j, m, mm, NUMF, NUMFT, ierr(2), NOUT1, NDFTR_el
   character*80         :: outfil

   integer              :: NFREQ_out
   real(8), Allocatable :: ScaleFact(:)
   integer              :: nbod, NQACC  !for pitch imbalance


   open (105,file='eigen.inp')
     read (105,*) NFREQ_out
     Allocate ( ScaleFact(0:NFREQ_out) )

     ScaleFact(0          ) = 0.d0
     ScaleFact(1:NFREQ_out) = 1.d0

     do i=1,NFREQ_out
        read (105,*)j,ScaleFact(i)
     enddo
   close(105)


   NDFT2_el = 2*NDFT0_el
   LWORK    = 8*NDFT2_el+16

   write (*,*) 'EIGEN NEW is called',NDFT0_el,NDFT2_el

   Allocate  (  AM_INV  (  NDFT0_el,   NDFT0_el) )
   Allocate  (  FTI_el  (  NDFT0_el            ) )
   Allocate  (  FTR1_el (  NDFT0_el            ) )
   Allocate  (  FTI1_el (  NDFT0_el            ) )
   Allocate  (  INDXE   (  NDFT0_el            ) )
   Allocate  (  ALFRS   (  NDFT2_el            ) )
   Allocate  (  ALFIS   (  NDFT2_el            ) )
   Allocate  (  RFR     (  NDFT2_el            ) )
   Allocate  (  RFI     (  NDFT2_el            ) )
   Allocate  (  IND     (  NDFT2_el            ) )
   Allocate  (  IND2    (0:NDFT2_el            ) )
   Allocate  (  RZR     (  NDFT2_el, 0:NDFT2_el) )
   Allocate  (  RZI     (  NDFT2_el,   NDFT2_el) )
   Allocate  (  RZR1    (  NDFT2_el,   NDFT2_el) )
   Allocate  (  RZI1    (  NDFT2_el,   NDFT2_el) )
   Allocate  (  VS      (  NDFT2_el,   NDFT2_el) )
   Allocate  (  AAS     (  NDFT2_el,   NDFT2_el) )
   Allocate  (  WORK    (  LWORK               ) )

         NOUT1 = 98
   open (NOUT1, file='eigen.dat')
   write(NOUT1,'(a)')"#Freq_d[Hz]   Zita[%]  Freq_n[Hz]"


!--- Form the matrix
   AM_INV(1:NDFT0_el,1:NDFT0_el) = AM_el(1:NDFT0_el,1:NDFT0_el);

   call DGETRF ( NDFT0_el, NDFT0_el, AM_INV, NDFT0_el, INDXE,              ierr(1) )
   call DGETRI ( NDFT0_el,           AM_INV, NDFT0_el, INDXE, WORK, LWORK, ierr(2) )

   if     (ierr(1)/=0) then
      write(*,*) 'error in DGETRF in EIGEN',ierr(1)
      stop
   elseif (ierr(2)/=0) then
      write(*,*) 'error in DGETRI in EIGEN',ierr(2)
      stop
   endif

   AAS(1:NDFT2_el,1:NDFT2_el) = 0.d0;
   
   do i=1,NDFT0_el
      do j=1,NDFT0_el
         if (i == j) AAS(i,j+NDFT0_el) = 1.d0

         do m=1,NDFT0_el
            AAS(i+NDFT0_el,j         )= AAS(i+NDFT0_el,j         ) - AM_INV(i,m)*AK_el(m,j)
            AAS(i+NDFT0_el,j+NDFT0_el)= AAS(i+NDFT0_el,j+NDFT0_el) - AM_INV(i,m)*AC_el(m,j)
         enddo
      enddo
   enddo


!--- Solve the eigen-value problem
   call dgeev ( 'N', 'V', NDFT2_el, AAS, NDFT2_el, ALFRS, ALFIS, VL, 1, VS, &
                          NDFT2_el, WORK, LWORK, ierr(1) )

   if     (ierr(1)/=0) then
      write(*,*) 'error in DGEEV in EIGEN',ierr(1)
      stop
   endif

   RZR  (:,0) = 0.d0;
   NUMF = 0
   do I = 1, NDFT2_el
      A1 = ALFRS(I)
      C1 = ALFIS(I)
      
      if ( C1 < 1d-10 ) cycle   !to avoid the symmetric negative eigen-values

      NUMF       = NUMF + 1
      RFR (NUMF) = A1/PI2
      RFI (NUMF) = C1/PI2
    
      do J=1,NDFT0_el 
         RZR  (J,NUMF) = VS(J,I+1)
         RZI  (J,NUMF) = VS(J,I  )
      enddo
   enddo

   write(*,*) 'number of positive eigenvalues',NUMF
   NUMFT    = NUMF
   NDFTR_el = min(NUMFT, NFREQ_out)
   

!--- Sort modes with respect to the frequency
   call SORTRX (RFI, NDFT2_el, NUMFT, IND)
   IND2(0         ) = 0
   IND2(1:NDFT2_el) = IND(1:NDFT2_el)


   do I     = 0, NDFTR_el

      write(*,*) I
!------ set the initial values for basic parameters
      call INIT_DEFLECTIONS ( -OMEGA )


!------ Substitute the modal - solution to all active dofs
      do mm = 1, NDFT0_el
         m  = INDSYSB (mm)
         UT_el(m) = RZR(mm,IND2(I))
      enddo


      if (IYNmodal==0) then
         call GET_REDUSED_DOFS
!--------- Scale Modes to make deformations visible!
         UT_el = UT_el * ScaleFact(I)

!--------- initial blade pitch angles
         do nbod  = 1, NBLADE_el
            NQACC = ACCU%NQSACC_el(nbod-1) 
            UT_el (NDFBT_el+NQSP+NQACC+5) = PITCH0_el(nbod) + PITCH_Imb_el(nbod)
         enddo
# ifdef HAVE_MODAL
      else

         call GET_REDUSED_DOFS_md
!--------- Scale Modes to make deformations visible!
         UT_el = UT_el * ScaleFact(I)

         call modal2fem 
# endif
      endif


      call ROTMAT_el
      
!------ Write global deformations
      outfil = 'modeR'
      if (IYNmodal==0) call WRITEOUT_GLOB (UT_el, outfil,I,3)
#    ifdef HAVE_MODAL
      if (IYNmodal==1) call WRITEOUT_GLOB (UT_md, outfil,I,3)
#    endif
!------ Write local deformations
      outfil = 'mode'
      if (IYNmodal==0) call WRITEOUT_LOC  (UT_el, outfil,I,3)
#    ifdef HAVE_MODAL
      if (IYNmodal==1) call WRITEOUT_LOC  (UT_md, outfil,I,3)
#    endif

      if (I==0) cycle

!------ lamba = -z wn +/- i wn *sqrt(1-z^2) --> ratio = -Re(lamda)/Im(lamda) = z/sqrt(1-z^2)
      RATIO =-RFR(IND2(I))/RFI(IND2(I))
      ZCRIT = RATIO / dsqrt(1.d0+RATIO**2)
!               Omega_d (reduced due to damping), zita [%]    , Omega_n
      write (NOUT1,'(4(F18.8,5x))') RFI(IND2(I)), ZCRIT*100.d0, RFI(IND2(I)) / dsqrt(1.d0-ZCRIT**2)
   enddo !I

   close (NOUT1)

   Deallocate  ( AM_INV, FTI_el, FTR1_el, FTI1_el, INDXE      )
   Deallocate  ( ALFRS , ALFIS , RFR    , RFI    , IND  , IND2)
   Deallocate  ( RZR   , RZI   , RZR1   , RZI1   , VS   , AAS )
   Deallocate  ( WORK                                         )


 END Subroutine EIGEN
!-------------------------------------------------------------------------------------------------
!
! -- Subroutine :MODAL_ANAL   
!
!------------------------------------------------------------------------------------------------- 
 Subroutine MODAL_ANAL
!------------------------------------------------------------------------------------------------- 

 Use Cbeam

   implicit None

   real(8), allocatable :: AK        (:,:)
   real(8), allocatable :: AC        (:,:)
   real(8), allocatable :: AM        (:,:)
   real(8), allocatable :: ALFI      (:  )
   real(8), allocatable :: ALFR      (:  )
   real(8), allocatable :: BETA      (:  )
   real(8), allocatable :: V         (:,:)
   real(8), allocatable :: RFR       (:  )
   real(8), allocatable :: RFI       (:  )
   real(8), allocatable :: RZR       (:,:)
   real(8), allocatable :: Md        (:,:)
   real(8), allocatable :: Cd        (:  )
   real(8), allocatable :: CCRIT_act (:  )
   real(8), allocatable :: Cd_FTTI_el(:,:)
   real(8), allocatable :: FTT_el    (:,:)
   real(8), allocatable :: FTTI_el   (:,:)
   real(8), allocatable :: FREQ_el   (:  )
   integer, allocatable :: IND       (:  )
   integer              :: NDF1,NDF2,nbod,NDFTB,i,j
   integer              :: m, ierr, nb1_el, nb2_el , mmax, I1, Inorm_mb
   real(8)              :: A1, C1 !!, RR, RI, RFI21, RFI22, RFI2, RI1, RI2, RR1, RR2
   real(8)              :: CCRIT
!-- DGGEV
   real(8)              :: VL(1)
   real(8), allocatable :: WORK  (:)
   integer              :: LWORK


   write (* ,*)
   write (* ,*) 'Modal Damping calculation'
   write (* ,*)
   write (10,*)
   write (10,*) 'Modal Damping calculation'
   write (10,*)

   CDAMPSTR = 0.d0;

   open (100, file='modal_anal_eigen.dat' )

!--- for all physical bodies
   do nbod = 1, min(NBODBT_el,NBODBTWT_el+1)

!------ consider the reduction of the zero dofs
!------ foundation is automatically taken into account by selection 
!------ the total number of elastic dofs, because foundation always
!------ is applied at the last body [tower/monopile/jacket]
      if (nbod<=NBODBTWT_el) then

         if (IREDUCT_el == 1) then
            nb1_el = body(nbod)%NBODTGLB_el(1) - 1
            nb2_el = body(nbod)%NBODTGLB_el(body(nbod)%NBODSUB_el)
         else
            nb1_el = 0
            nb2_el = 0
         endif

         NDF1  = ACCU%NDFTBACC_el(nbod-1) - nb1_el * 6
         NDF2  = ACCU%NDFTBACC_el(nbod  ) - nb2_el * 6

!--------- foundation consideration
         if ( ITYPE_FOUND_el >  0 .and.     & !foundation consideration
              ICASE_el       <= 1 .and.     & !foundation for onshore or monopile concept
              nbod           == NTOWER_el ) & !tower is the last body
              NDF2  = NDFBT0_el
      else
         write(*,*)'JACKET modal',nbod

         if (IREDUCT_el == 1) then
            nb1_el = body(nbod)%NBODTGLB_el(1) - 1
         else
            nb1_el = 0
         endif

         NDF1  = ACCU%NDFTBACC_el(nbod-1) - nb1_el * 6
         NDF2  = NDFBT0_el
      endif

      NDFTB = NDF2 - NDF1
      LWORK = 8*NDFTB +16

      Allocate ( AK        (NDFTB, NDFTB) )
      Allocate ( AC        (NDFTB, NDFTB) )
      Allocate ( AM        (NDFTB, NDFTB) )
      Allocate ( ALFI      (NDFTB       ) )
      Allocate ( ALFR      (NDFTB       ) )
      Allocate ( BETA      (NDFTB       ) )
      Allocate ( V         (NDFTB, NDFTB) )

      Allocate ( RFR       (NDFTB       ) )
      Allocate ( RFI       (NDFTB       ) )
      Allocate ( RZR       (NDFTB, NDFTB) )
      Allocate ( IND       (NDFTB       ) )
      Allocate ( FREQ_el   (NDFTB       ) )
      Allocate ( Cd        (NDFTB       ) )
      Allocate ( CCRIT_act (NDFTB       ) )
      Allocate ( Cd_FTTI_el(NDFTB, NDFTB) )

      Allocate ( Md        (NDFTB, NDFTB) )
      Allocate ( FTT_el    (NDFTB, NDFTB) )
      Allocate ( FTTI_el   (NDFTB, NDFTB) )

      Allocate ( WORK      (LWORK       ) )


      write(* ,*)'MODAL_ANAL',nbod,NDF1,NDF2,NDFTB
      write(10,*)'MODAL_ANAL',nbod,NDF1,NDF2,NDFTB

      AK(1:NDFTB,1:NDFTB) = AK_el (NDF1+1:NDF2,NDF1+1:NDF2);
      AM(1:NDFTB,1:NDFTB) = AM_el (NDF1+1:NDF2,NDF1+1:NDF2);

!------ Solve the eigenvalue problem
!     call REIGG(NDFTB, NDFTB   ,AK,AM,ALFR,ALFI,BETA,1,V,ierr)
      call DGGEV('N', 'V', NDFTB, AK , NDFTB, AM, NDFTB, ALFR, ALFI, BETA, VL, 1, V, NDFTB, WORK, LWORK, ierr )

      if (ierr>0) write(*,*) 'Modal Anal - DGGEV', ierr, WORK(1)/LWORK

      AK(1:NDFTB,1:NDFTB) = AK_el (NDF1+1:NDF2,NDF1+1:NDF2);
      AM(1:NDFTB,1:NDFTB) = AM_el (NDF1+1:NDF2,NDF1+1:NDF2);


!------ Check the eigen-values. Set the frequency and the mode matrix.
!------ lamda = (AKFR/BETA) + (ALFI/BETA) yiot = A1 + C1 yiot
      do I = 1, NDFTB

         write(333,'(2i5,5e15.5)') nbod, I, ALFR(I), ALFI(I), BETA(I) !, dsqrt(ALFR(I)/BETA(I))/PI2
!--------- check BETA
         if (dabs(BETA(I)) < 1.d-10) then
            write(10,*) 'beta = 0',I
            write( *,*) 'beta = 0',I

            RFR (  I) = 1.d+10 !0.d0
            RFI (  I) = 0.d0
            RZR (:,I) = 0.d0;
            RZR (I,I) = 1.d0
            cycle
         endif

!--------- check Imaginary Part and Real Part sign
         if (dabs(ALFI(I)) > 1.d-10) then
            write(10,*) 'Complex in modal anal transformed to real'
            write(* ,*) 'Complex in modal anal transformed to real'
            ALFI(I)=0.d0
         endif
         A1    = ALFR(I)/BETA(I)
         C1    = ALFI(I)/BETA(I)

         if (dabs(C1) < 1.d-13) then

!------------ Real Eigen-Value -----------
            if     (A1 <  0.d0) then

!--------------- Real Part is Negative --
               write(10,*) 'lamda REAL, A1 neg',I, A1 !, C1
               write( *,*) 'lamda REAL, A1 neg',I, A1 !, C1

               RFR (  I) = 1.d+10 !0.d0
               RFI (  I) = 0.d0
               RZR (:,I) = 0.d0;
               RZR (I,I) = 1.d0

            elseif (A1 >= 0.d0) then

!--------------- Real Part is Positive [Normal case] --
!w             write( *,*) 'lamda REAL, A1 pos',I, A1, dsqrt(A1)/PI2 !, C1

               RFR (  I) = dsqrt(A1) / PI2
               RFI (  I) = 0.d0
               RZR (:,I) = V(:,I); ! Fi [ndfbxndfb] or [ndfbxndftr] when reduced
            endif !Real pos or Real neg

         else

!------------ Complex Eigen-Value -----------
            write(10,*) 'lamda COMPLEX',I,C1,A1
            write( *,*) 'lamda COMPLEX',I,C1,A1

!           STOP
            RFR (  I) = 1.d+10 !0.d0
            RFI (  I) = 0.d0
            RZR (:,I) = 0.d0;
            RZR (I,I) = 1.d0
            cycle

      !!    RFI21 = 0.5d0*(-A1 + dsqrt(A1**2+C1**2))
      !!    RFI22 = 0.5d0*(-A1 - dsqrt(A1**2+C1**2))
      !!    RFI2  = dmax1(RFI21,RFI22)
      !!    RI1   = dsqrt(RFI2)
      !!    RI2   =-dsqrt(RFI2)
      !!    RR1   = C1/(2.d0*RI1)
      !!    RR2   = C1/(2.d0*RI2)

!     !!   write(*,*)  RFI21, RFI22, RFI2,RI1,RI2
!     !!   write(*,*)  RR1  , RR2

      !!    if (RR1.ge.0.d0) then
      !!        RR = RR1
      !!        RI = RI1
      !!    else
      !!        RR = RR2
      !!        RI = RI2
      !!    endif
      !!    RFR (  I) = RR / PI2
      !!    RFI (  I) = RI / PI2
      !!    RZR (:,I) = V(:,I); ! Fi [ndfbxndfb] or [ndfbxndftr] when reduced
         endif !Real or Complex
      enddo !I


!------ Sort modes with respect to the frequency
      call SORTRX (RFR ,NDFTB,NDFTB,IND)

      FREQ_el (:  ) = 0.d0;
      Md      (:,:) = 0.d0;
      Cd      (:  ) = 0.d0;
      FTT_el  (:,:) = 0.d0;


!------ Normalize Modes
      do I    = 1, NDFTB
         I1   = IND(I)
         mmax = 1
         do m = 2, NDFTB
            if ( dabs(RZR(m,I1)) > dabs(RZR(mmax,I1)) ) mmax=m
         enddo

         FREQ_el    (I) = RFR(I1)

         if (RZR(mmax,I1) == 0.d0) cycle

         FTT_el (1:NDFTB,I) = RZR(1:NDFTB,I1) / RZR(mmax,I1);
      enddo


!------ Inverse FTT_el
      FTTI_el = FTT_el;

      call DGETRF (NDFTB, NDFTB, FTTI_el, NDFTB, IND, ierr)
      if (ierr>0) write(*,*) 'Modal Anal - DGETRF', ierr
      call DGETRI( NDFTB, FTTI_el, NDFTB, IND, WORK, LWORK, ierr )
      if (ierr>0) write(*,*) 'Modal Anal - DGETRI', ierr


!------ Determine modal masses
      Md = matmul ( Transpose(FTT_el), matmul(AM, FTT_el) )


!------ Determine modal damping
      Inorm_mb =  0 !0: normalize, -1: no

 1    do i = 1, NDFTB
                                      CCRIT=CCRIT_el(nbod)
         if ( FREQ_el(i) > FREQ0_el ) CCRIT=CCRIT_el(0)

!!!         if (Inorm_mb == 1) &
!!!         CCRIT    = CCRIT * (CCRIT/ CCRIT_act(i))
         if     (Inorm_mb == 1.and.dabs(CCRIT_act(i))>1.d-15) then
            CCRIT    = CCRIT * (CCRIT/ CCRIT_act(i))
         elseif (Inorm_mb == 1                              ) then
            CCRIT    = 0.d0
         endif

         Cd   (i) = 2.d0*CCRIT*PI2*FREQ_el(i)*Md(i,i)
      enddo

!------ Determine Structural damping
      do i  = 1, NDFTB
         Cd_FTTI_el(i, 1:NDFTB) = Cd(i) * FTTI_el(i, 1:NDFTB)
      enddo

      AC (1:NDFTB, 1:NDFTB) = matmul ( transpose(FTTI_el(1:NDFTB, 1:NDFTB)), Cd_FTTI_el(1:NDFTB, 1:NDFTB) )

      if (Inorm_mb == 0) then
         call EIGEN_MODAL ( AM, AC, AK, NDFTB, CCRIT_act, Inorm_mb )
        if (Inorm_mb == 1)  goto 1
      endif


!------ Assembly to global Structural Damping Matrix
      CDAMPSTR (NDF1+1:NDF1+NDFTB, NDF1+1:NDF1+NDFTB) = AC(1:NDFTB, 1:NDFTB)


!------ Write frequency, Modal mass and Modal damping
      do i=1,NDFTB
         write (100,101) FREQ_el(i), Md(i,i), Cd(i)
      enddo

      write (100,*)
      write (100,*)

      Deallocate ( AK , AM    , ALFI   , ALFR, BETA, V )
      Deallocate ( RFR, RFI   , RZR    , IND           )
      Deallocate ( Md , FTT_el, FTTI_el                )
      Deallocate ( Cd , CCRIT_act, Cd_FTTI_el, AC      )
      Deallocate ( FREQ_el )
      Deallocate ( WORK    )
   enddo!nbod

   close (100)

 101  format ('freq=',F17.3,5x,'modal mass=',F15.3,'modal damp=',F17.3 )


!--- write the total damping matrix
   open(200,file='CSTR.inp') !,form='UNFORMATTED')
    write(200,*) NDFT0_el, NDFT_el
   !write(200  ) NDFT0_el

    do j=1,NDFT0_el
       do i=1,NDFT0_el
          write(200,'(e28.17)') CDAMPSTR(i,j)
         !write(200           ) CDAMPSTR(i,j)
       enddo
    enddo

   close (200)


 END Subroutine MODAL_ANAL
!-------------------------------------------------------------------------------------------------
!
! -- Subroutine :EIGEN
!
!------------------------------------------------------------------------------------------------- 
 Subroutine EIGEN_MODAL ( AM, AC, AK, N, ZCRIT, I_OK )
!------------------------------------------------------------------------------------------------- 

 Use Cbeam

   implicit None

   integer, intent(in)  :: N
   real(8), intent(in)  :: AM     (N,N)
   real(8), intent(in)  :: AC     (N,N)
   real(8), intent(in)  :: AK     (N,N)
   real(8), intent(out) :: ZCRIT  (N)
   integer, intent(out) :: I_OK

   real(8), Allocatable :: AM_INV (:,:)             !(NDFT0_el,NDFT0_el)
   integer, Allocatable :: INDXE  (:  )             !(NDFT0_el)
   real(8), Allocatable :: ALFRS  (:  )
   real(8), Allocatable :: ALFIS  (:  )
   real(8), Allocatable :: RFR    (:  )
   real(8), Allocatable :: RFI    (:  )             !(2*NDFT0_el)
   integer, Allocatable :: IND    (:  )             !(2*NDFT0_el)
   real(8), Allocatable :: VS     (:,:)
   real(8), Allocatable :: AAS    (:,:)             !(2*NDFT0_el,2*NDFT0_el)
   real(8), Allocatable :: WORK   (:  )             !(LWORK)
   real(8)              :: VL(1)
   integer              :: N2, LWORK
   real(8)              :: RATIO
   integer              :: i, j, m, NUMF, ierr


   I_OK  = 1   !by default 
   N2    = 2*N
   LWORK = 8*N2+16

   write (*,*) 'EIGEN_MODAL is called',N, N2

   Allocate ( AM_INV (N , N ) )
   Allocate ( INDXE  (N     ) )
   Allocate ( ALFRS  (N2    ) )
   Allocate ( ALFIS  (N2    ) )
   Allocate ( RFR    (N2    ) )
   Allocate ( RFI    (N2    ) )
   Allocate ( IND    (N2    ) )
   Allocate ( VS     (N2, N2) )
   Allocate ( AAS    (N2, N2) )
   Allocate ( WORK   (LWORK ) )


!--- Form the matrix
   AM_INV(1:N,1:N) = AM(1:N,1:N);

   call DGETRF ( N, N, AM_INV, N, INDXE,              ierr )
   call DGETRI ( N,    AM_INV, N, INDXE, WORK, LWORK, ierr )

   AAS(1:N2,1:N2) = 0.d0;
   
   do i = 1, N
    do j = 1, N
!!!    if (i == j) AAS(i,j+N) = 1.d0

       do m=1,N
          AAS(i+N,j  )= AAS(i+N,j  ) - AM_INV(i,m)*AK(m,j)
          AAS(i+N,j+N)= AAS(i+N,j+N) - AM_INV(i,m)*AC(m,j)
       enddo
    enddo
          AAS(i,i+N) = 1.d0
   enddo


!--- Solve the eigen-value problem
   call dgeev ( 'N', 'V', N2, AAS, N2, ALFRS, ALFIS, VL, 1, VS, &
                          N2, WORK, LWORK, ierr )

   NUMF = 0
   do i = 1, N2
      write(3500,*) i,ALFIS(i)
      if ( ALFIS(i) < 1.d-13 ) cycle   !to avoid the symmetric negative eigen-values
      write(3501,*) i,ALFIS(i)

      NUMF       = NUMF + 1
      RFR (NUMF) = ALFRS(i) / PI2
      RFI (NUMF) = ALFIS(i) / PI2
   enddo

!--- Sort modes with respect to the frequency
   call SORTRX (RFI, N2, NUMF, IND)

  open (105, file='eigen_modal_anal_eigen.dat', access='append')
      write (105,*)'next bondy'

   do i = 1, NUMF
!------ lamba = -z wn +/- i wn *sqrt(1-z^2) --> ratio = -Re(lamda)/Im(lamda) = z/sqrt(1-z^2)
      RATIO    =-RFR(IND(I))/RFI(IND(I))
      ZCRIT(i) = RATIO / dsqrt(1.d0+RATIO**2)
      write (105,'(3(F18.8,5x))') RFI(IND(i)), ZCRIT(i)*100.d0, RFI(IND(I)) / dsqrt(1.d0-ZCRIT(i)**2)
   enddo !i
      write (105,*)
      write (105,*)
  close (105)


!---- check if all eigenvalues are well collected
   if (NUMF /= N) then
      write(*,*) 'error in EIGEN_MODAL',N,NUMF
      I_OK = 0
   endif
   

   Deallocate  ( AM_INV, INDXE                 )
   Deallocate  ( ALFRS , ALFIS, RFR , RFI, IND )
   Deallocate  ( VS    , AAS  , WORK           )


 END Subroutine EIGEN_MODAL
!------------------------------------------------------------------------------------------------------- 
 Subroutine STRUCT_DAMP
!------------------------------------------------------------------------------------------------------- 
  
 Use Cbeam
 Use truss
  

   implicit None


   integer :: i, j, js, NDFB


   if (IMODAL == 0) return

                     NDFB = NDFBT0_el
   if (ITRUSS_el==2) NDFB = NDFBT0_el - NDFT_tr

   do j  = 1, NDFB
      js = INDSYSB (j)

      do i = 1, NDFB
         AC_el(i,j) = AC_el(i,j) + CDAMPSTR(i,j)
         AQ_el(i  ) = AQ_el(i  ) - CDAMPSTR(i,j)*UT1_el(js)
      enddo!i
   enddo!j


 END Subroutine STRUCT_DAMP
!
!
!
!=======================================================================
 SUBROUTINE SORTRX (DATA,NP,N,INDEX)
!=======================================================================
!
!     SORTRX -- SORT, Real input, indeX output
!
!
!     Input:  N     INTEGER
!             DATA  REAL
!
!     Output: INDEX INTEGER (DIMENSION N)
!
! This routine performs an in-memory sort of the first N elements of
! array DATA, returning into array INDEX the indices of elements of
! DATA arranged in ascending order.  Thus,
!
!    DATA(INDEX(1)) will be the smallest number in array DATA;
!    DATA(INDEX(N)) will be the largest number in DATA.
!
! The original data is not physically rearranged.  The original order
! of equal input values is not necessarily preserved.
!
!===================================================================
!
! SORTRX uses a hybrid QuickSort algorithm, based on several
! suggestions in Knuth, Volume 3, Section 5.2.2.  In particular, the
! "pivot key" [my term] for dividing each subsequence is chosen to be
! the median of the first, last, and middle values of the subsequence;
! and the QuickSort is cut off when a subsequence has 9 or fewer
! elements, and a straight insertion sort of the entire array is done
! at the end.  The result is comparable to a pure insertion sort for
! very short arrays, and very fast for very large arrays (of order 12
! micro-sec/element on the 3081K for arrays of 10K elements).  It is
! also not subject to the poor performance of the pure QuickSort on
! partially ordered data.
!
! Created:  15 Jul 1986  Len Moss
!
!===================================================================

      IMPLICIT NONE
      INTEGER      :: NP,N
      INTEGER      :: INDEX(NP)
      REAL(KIND=8) :: DATA(NP)

      INTEGER      :: LSTK(31),RSTK(31),ISTK
      INTEGER      :: L,R,I,J,P,INDEXP,INDEXT
      REAL(KIND=8) :: DATAP

!     QuickSort Cutoff
!
!     Quit QuickSort-ing when a subsequence contains M or fewer
!     elements and finish off at end with straight insertion sort.
!     According to Knuth, V.3, the optimum value of M is around 9.

      INTEGER , PARAMETER ::M=9

!===================================================================
!
!     Make initial guess for INDEX

      DO 50 I=1,N
         INDEX(I)=I
   50    CONTINUE

!     If array is short, skip QuickSort and go directly to
!     the straight insertion sort.

      IF (N.LE.M) GOTO 900

!===================================================================
!
!     QuickSort
!
!     The "Qn:"s correspond roughly to steps in Algorithm Q,
!     Knuth, V.3, PP.116-117, modified to select the median
!     of the first, last, and middle elements as the "pivot
!     key" (in Knuth's notation, "K").  Also modified to leave
!     data in place and produce an INDEX array.  To simplify
!     comments, let DATA[I]=DATA(INDEX(I)).

! Q1: Initialize
      ISTK=0
      L=1
      R=N

  200 CONTINUE

! Q2: Sort the subsequence DATA[L]..DATA[R].
!
!     At this point, DATA[l] <= DATA[m] <= DATA[r] for all l < L,
!     r > R, and L <= m <= R.  (First time through, there is no
!     DATA for l < L or r > R.)

      I=L
      J=R

! Q2.5: Select pivot key
!
!     Let the pivot, P, be the midpoint of this subsequence,
!     P=(L+R)/2; then rearrange INDEX(L), INDEX(P), and INDEX(R)
!     so the corresponding DATA values are in increasing order.
!     The pivot key, DATAP, is then DATA[P].

      P=(L+R)/2
      INDEXP=INDEX(P)
      DATAP=DATA(INDEXP)

      IF (DATA(INDEX(L)) .GT. DATAP) THEN
         INDEX(P)=INDEX(L)
         INDEX(L)=INDEXP
         INDEXP=INDEX(P)
         DATAP=DATA(INDEXP)
      ENDIF

      IF (DATAP .GT. DATA(INDEX(R))) THEN
         IF (DATA(INDEX(L)) .GT. DATA(INDEX(R))) THEN
            INDEX(P)=INDEX(L)
            INDEX(L)=INDEX(R)
         ELSE
            INDEX(P)=INDEX(R)
         ENDIF
         INDEX(R)=INDEXP
         INDEXP=INDEX(P)
         DATAP=DATA(INDEXP)
      ENDIF

!     Now we swap values between the right and left sides and/or
!     move DATAP until all smaller values are on the left and all
!     larger values are on the right.  Neither the left or right
!     side will be internally ordered yet; however, DATAP will be
!     in its final position.

  300 CONTINUE

! Q3: Search for datum on left >= DATAP
!
!     At this point, DATA[L] <= DATAP.  We can therefore start scanning
!     up from L, looking for a value >= DATAP (this scan is guaranteed
!     to terminate since we initially placed DATAP near the middle of
!     the subsequence).

         I=I+1
         IF (DATA(INDEX(I)).LT.DATAP) GOTO 300

  400 CONTINUE

! Q4: Search for datum on right <= DATAP
!
!     At this point, DATA[R] >= DATAP.  We can therefore start scanning
!     down from R, looking for a value <= DATAP (this scan is guaranteed
!     to terminate since we initially placed DATAP near the middle of
!     the subsequence).

         J=J-1
         IF (DATA(INDEX(J)).GT.DATAP) GOTO 400

! Q5: Have the two scans collided?

      IF (I.LT.J) THEN

! Q6: No, interchange DATA[I] <--> DATA[J] and continue

         INDEXT=INDEX(I)
         INDEX(I)=INDEX(J)
         INDEX(J)=INDEXT
         GOTO 300
      ELSE

! Q7: Yes, select next subsequence to sort
!
!     At this point, I >= J and DATA[l] <= DATA[I] == DATAP <= DATA[r],
!     for all L <= l < I and J < r <= R.  If both subsequences are
!     more than M elements long, push the longer one on the stack and
!     go back to QuickSort the shorter; if only one is more than M
!     elements long, go back and QuickSort it; otherwise, pop a
!     subsequence off the stack and QuickSort it.

         IF (R-J .GE. I-L .AND. I-L .GT. M) THEN
            ISTK=ISTK+1
            LSTK(ISTK)=J+1
            RSTK(ISTK)=R
            R=I-1
         ELSE IF (I-L .GT. R-J .AND. R-J .GT. M) THEN
            ISTK=ISTK+1
            LSTK(ISTK)=L
            RSTK(ISTK)=I-1
            L=J+1
         ELSE IF (R-J .GT. M) THEN
            L=J+1
         ELSE IF (I-L .GT. M) THEN
            R=I-1
         ELSE
! Q8: Pop the stack, or terminate QuickSort if empty
            IF (ISTK.LT.1) GOTO 900
            L=LSTK(ISTK)
            R=RSTK(ISTK)
            ISTK=ISTK-1
         ENDIF
         GOTO 200
      ENDIF

  900 CONTINUE

!===================================================================
!
! Q9: Straight Insertion sort

      DO 950 I=2,N
         IF (DATA(INDEX(I-1)) .GT. DATA(INDEX(I))) THEN
            INDEXP=INDEX(I)
            DATAP=DATA(INDEXP)
            P=I-1
  920       CONTINUE
               INDEX(P+1) = INDEX(P)
               P=P-1
               IF (P.GT.0) THEN
                  IF (DATA(INDEX(P)).GT.DATAP) GOTO 920
               ENDIF
            INDEX(P+1) = INDEXP
         ENDIF
  950    CONTINUE

!===================================================================
!
!     All done

 END SUBROUTINE SORTRX
