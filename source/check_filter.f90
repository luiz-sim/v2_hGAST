 Program Check_filter

   implicit none

   integer, parameter :: NQSC0 = 6
   real(8), save      :: QTC_el (NQSC0), QTC1_el (NQSC0), QTC2_el (NQSC0)
   real(8), save      :: QTCP_el(NQSC0), QTCP1_el(NQSC0), QTCP2_el(NQSC0)

   real(8), save      :: TIME, DT_el, PI2, TPERIOD_el
   real(8)            :: input(13001)
   integer            :: ntime

   integer            :: ieq, ieq1, ieq2
   real(8)            :: W,W1,W2, D,D1,D2, VAR, DVAR, DDVAR

   real(8)            :: A(4,4),B(4),C(4),DD, YY,YY1


!--- initialization
   QTC_el (:)=0.d0;   QTCP_el (:)=0.d0;
   QTC1_el(:)=0.d0;   QTCP1_el(:)=0.d0;
   QTC2_el(:)=0.d0;   QTCP2_el(:)=0.d0;

   TPERIOD_el = 6.25d0
   DT_el      = 0.05d0
   TIME       = 0.d0
   PI2        = 2.d0*dacos(-1.d0)

!--- read input signal
   open(1,file="input.dat")
    do ntime = 1, 13001
    read(1,*) var, input(ntime)
    enddo
   close(1)

!--- time domain realization of the filter
   do ntime = 1, 13001
      TIME  = TIME + DT_el 
      QTCP_el (1:NQSC0) = QTC_el (1:NQSC0);
      QTCP1_el(1:NQSC0) = QTC1_el(1:NQSC0);
      QTCP2_el(1:NQSC0) = QTC2_el(1:NQSC0);

!--- Tower Acceleration: 2nd-order low-pass filter at 2.5P
      ieq1             = 1
!     W                = 2.5d0 * (PI2/TPERIOD_el)    ![rad/s]
      W                = 1.875d0 * (PI2/TPERIOD_el)    ![rad/s]
      D                = 0.8d0
      VAR              = input(ntime)
      
      call filter_lowpass2 ( W             , D              , DT_el          , &
                             VAR                                             , &
                             QTCP_el (ieq1), QTCP1_el (ieq1), QTCP2_el (ieq1), &
                             QTC_el  (ieq1), QTC1_el  (ieq1),  QTC2_el (ieq1)    )
      
!--- Tower Acceleration: 2nd-order Notch filter at 3P
      ieq2             = 2
      W1               = 3.0d0 * (PI2/TPERIOD_el)    ![rad/s]
      W2               = W1
      D1               = 0.001d0
      D2               = 0.200d0
      VAR              = QTC_el (ieq1)
      DVAR             = QTC1_el(ieq1)
      DDVAR            = QTC2_el(ieq1)
      
      call filter_notch2 ( W1   , W2   , D1   , D2   , DT_el                , &
                           VAR  , DVAR , DDVAR                              , &
                           QTCP_el (ieq2) , QTCP1_el (ieq2), QTCP2_el (ieq2), &
                           QTC_el  (ieq2) , QTC1_el  (ieq2),  QTC2_el (ieq2)    )
      
!--- low pass filter 0.3Hz 10db   A(4x4), B(4x1), C(1x4), D(1x1)
      A(1,1:4) = (/ -0.7806d0, -1.6650d0,  0.0000d0,  0.0000d0 /)
      A(2,1:4) = (/  1.6650d0,  0.0000d0,  0.0000d0,  0.0000d0 /)
      A(3,1:4) = (/ -0.7806d0,  2.0087d0, -0.0158d0, -1.8830d0 /)
      A(4,1:4) = (/  0.0000d0,  0.0000d0,  1.8830d0,  0.0000d0 /)
      B(1:4  ) = (/  1.8850d0,  0.0000d0,  1.8850d0,  0.0000d0 /)
      C(  1:4) = (/ -0.1310d0,  0.3370d0, -0.0027d0,  0.0046d0 /)
      DD       =     0.3162d0
      VAR      = input(ntime)

!    call filterN ( A,B,C,D , YP,YP1,YP2 ,Y,Y1,Y2, DT, N, XX,XX1, YY,YY1 )
     call filterN ( A,B,C,DD, QTCP_el(3:6),QTCP1_el(3:6),QTCP2_el(3:6), &
                              QTC_el (3:6),QTC1_el (3:6),QTC2_el (3:6), &
                    DT_el, 4, VAR,0.d0, YY,YY1 )

!--- Output
 2    open(27,file='filter.dat',access='append')
       write(27,'(1000e16.7)')              &
         sngl(TIME        )                ,& ! 1
         sngl(input(ntime))                ,& ! 2
        (sngl(QTC_el (ieq)),ieq=1,NQSC0),YY,& ! 3- 9
        (sngl(QTC1_el(ieq)),ieq=1,NQSC0)   ,& !10-15
        (sngl(QTC2_el(ieq)),ieq=1,NQSC0)      !16-21
      close(27)
   enddo !ntime


 END
!--------------------------------------------------------------------------------
!
!  Equation for 2nd-order BandPass Filter : Gain (2zw s) / (s^2 + 2zw s + w^2)
!
!--------------------------------------------------------------------------------
 Subroutine filter_bandpass2 ( Gain  , w     , z     ,&
                               DT    , VAR   , Q0M   ,&
                               QTCP  , QTCP1 , QTCP2 ,&
                               QTC   , QTC1  , QTC2     )
!--------------------------------------------------------------------------------
   implicit none

   real(8), intent(in   ) :: Gain, w, z, DT, VAR
   real(8), intent(in   ) :: QTCP, QTCP1, QTCP2
   real(8), intent(  out) :: QTC , QTC1 , QTC2
   real(8), intent(  out) :: Q0M

   real(8) :: AM, AC, AK, AQ
   real(8) :: BITA, GAMA


   AM   = 1.d0
   AC   = 2.d0*z*w
   AK   =        w**2
   AQ   = 2.d0*Gain*z*w *VAR

   call Newmark_solve1 (AM  , AC   , AK   , AQ, DT,&
                        QTCP, QTCP1, QTCP2,        &
                        QTC , QTC1 , QTC2            )

!--- derivative dQTC/dVAR, based on Newmark
   BITA = 0.25d0
   GAMA = 0.50d0
   Q0M  = 2.d0*Gain*z*w / (AM/BITA/DT**2 + AC*GAMA/BITA/DT + AK) !dQ0/d(VAR)
!  Q0M  = 2.d0*Gain*z*w / (AM/(GAMA*DT) + AC + AK*BITA*DT/GAMA)  !dQ1/d(VAR)


 END Subroutine filter_bandpass2
!--------------------------------------------------------------------------------
!
!  Equation for 2nd-order Notch Filter: (s^2 + 2z1w1 s + w1^2) / (s^2 + 2z2w2 s + w2^2)
!
!--------------------------------------------------------------------------------
 Subroutine filter_notch2 ( W1, W2, Z1, Z2, DT, &
                            VAR , DVAR , DDVAR, &
                            QTCP, QTCP1, QTCP2, &
                            QTC , QTC1 , QTC2     )
!--------------------------------------------------------------------------------
   implicit none

   real(8), intent(in   ) :: W1, W2, Z1, Z2, DT
   real(8), intent(in   ) :: VAR , DVAR , DDVAR
   real(8), intent(in   ) :: QTCP, QTCP1, QTCP2
   real(8), intent(  out) :: QTC , QTC1 , QTC2

   real(8) :: AM, AC, AK, AQ


   AM = 1.d0
   AC = 2.d0*Z2*W2
   AK =         W2**2
   AQ =               DDVAR  &
       +2.d0*Z1*W1   * DVAR  &
       +        W1**2*  VAR

   call Newmark_solve1 (AM, AC, AK, AQ, DT, &
                        QTCP, QTCP1, QTCP2, &
                        QTC , QTC1 , QTC2     )


 END Subroutine filter_notch2
!--------------------------------------------------------------------------------
!
!  Equation for 2nd-order low-pass filter: w^2 / (s^2      + 2z w s + w^2)
!                                         =1   / (s^2 /w^2 + 2z/w s + 1  )
!
!--------------------------------------------------------------------------------
 Subroutine filter_lowpass2 ( W   , Z    , DT   ,&
                              VAR               ,&
                              QTCP, QTCP1, QTCP2,&
                              QTC , QTC1 , QTC2    )
!--------------------------------------------------------------------------------
   implicit none

   real(8), intent(in   ) :: W     , Z     , DT
   real(8), intent(in   ) :: VAR
   real(8), intent(in   ) :: QTCP  , QTCP1 , QTCP2
   real(8), intent(  out) :: QTC   , QTC1  , QTC2

   real(8) :: AM, AC, AK, AQ


   AM = 1.d0
   AC = 2.d0*Z*W
   AK =        W**2
   AQ =        W**2 * VAR

   call Newmark_solve1 (AM, AC, AK, AQ, DT, &
                        QTCP, QTCP1, QTCP2, &
                        QTC , QTC1 , QTC2     )


 END Subroutine filter_lowpass2
!--------------------------------------------------------------------------------
!
!  Equation for 1st-order lowpass filter: 1 / (s (1+s τ))
!
!--------------------------------------------------------------------------------
 Subroutine filter_lowpass1 ( TVLAG, DT   , VAR  , &
                              QTCP , QTCP1, QTCP2, &
                              QTC  , QTC1 , QTC2 , &
                              Q0M                    )
!--------------------------------------------------------------------------------
   implicit none

   real(8), intent(in   ) :: TVLAG , DT
   real(8), intent(in   ) :: VAR
   real(8), intent(in   ) :: QTCP  , QTCP1 , QTCP2
   real(8), intent(  out) :: QTC   , QTC1  , QTC2
   real(8), intent(  out) :: Q0M

   real(8) :: AM, AC, AK, AQ
   real(8) :: BITA, GAMA


   AM   = TVLAG
   AC   = 1.d0
   AK   = 0.d0
   AQ   = VAR

   call Newmark_solve1 ( AM, AC, AK, AQ, DT, &
                         QTCP, QTCP1, QTCP2, &
                         QTC , QTC1 , QTC2     )

!--- derivative dQTC1/dVAR, based on Newmark
   BITA = 0.25d0
   GAMA = 0.50d0
   Q0M  = 1.d0 / (AM/GAMA/DT + AC + AK*BITA*DT/GAMA) !dQTC1/dVAR


 END Subroutine filter_lowpass1
!----------------------------------------------------------------------
 Subroutine ACTUATOR_1 ( TVLAG   , DT             , &
                         MinVAR  , MaxVAR         , &
                         MinDVAR , MaxDVAR        , &
                         MinDDVAR, MaxDDVAR       , &
                         VAR                      , &
                         QTCP    , QTCP1   , QTCP2, &
                         QTC     , QTC1    , QTC2     )
!-----------------------------------------------------------------------

   implicit none

   real(8), intent(in   ) :: TVLAG , DT
   real(8), intent(in   ) :: MinVAR, MaxVAR, MinDVAR, MaxDVAR, MinDDVAR, MaxDDVAR
   real(8), intent(in   ) :: VAR
   real(8), intent(in   ) :: QTCP  , QTCP1 , QTCP2
   real(8), intent(  out) :: QTC   , QTC1  , QTC2

   real(8) :: AM, AC, AK, AQ


   AM   = TVLAG
   AC   = 1.d0
   AK   = 0.d0
   AQ   = VAR

   call Newmark_solve1 ( AM, AC, AK, AQ, DT, &
                         QTCP, QTCP1, QTCP2, &
                         QTC , QTC1 , QTC2     )

   call Newmark_saturate1 ( DT                        , &
                            MinVAR  , MaxVAR          , &
                            MinDVAR , MaxDVAR         , &
                            MinDDVAR, MaxDDVAR        , &
                            QTCP    , QTCP1   , QTCP2 , &
                            QTC     , QTC1    , QTC2      )

 END Subroutine ACTUATOR_1
!--------------------------------------------------------------------------------
!
!  Equation for 2nd-order pitch actuator: (2z w s + w^2) / (s^2       + 2z w s + w^2) =
!                                         (2z/w s + 1  ) / (s^2 / w^2 + 2z/w s + 1  )
!
!  numinator order < denuminator order
!--------------------------------------------------------------------------------
 Subroutine ACTUATOR_2 ( W, D    , DT             , &
                         MinVAR  , MaxVAR         , &
                         MinDVAR , MaxDVAR        , &
                         MinDDVAR, MaxDDVAR       , &
                         VAR     , DVAR           , &
                         QTCP    , QTCP1   , QTCP2, &
                         QTC     , QTC1    , QTC2     )
!--------------------------------------------------------------------------------
   implicit none

   real(8), intent(in   ) :: W, D  , DT
   real(8), intent(in   ) :: MinVAR, MaxVAR, MinDVAR, MaxDVAR, MinDDVAR, MaxDDVAR
   real(8), intent(in   ) :: VAR   , DVAR
   real(8), intent(in   ) :: QTCP  , QTCP1 , QTCP2
   real(8), intent(  out) :: QTC   , QTC1  , QTC2

   real(8) :: AM, AC, AK, AQ


   AM = 1.d0
   AC = 2.d0*D*W
   AK = W**2
   AQ = W**2 * VAR + 2.d0*D*W * DVAR

   call Newmark_solve1 ( AM, AC, AK, AQ, DT, &
                         QTCP, QTCP1, QTCP2, &
                         QTC , QTC1 , QTC2     )

   call Newmark_saturate1 ( DT                        , &
                            MinVAR  , MaxVAR          , &
                            MinDVAR , MaxDVAR         , &
                            MinDDVAR, MaxDDVAR        , &
                            QTCP    , QTCP1   , QTCP2 , &
                            QTC     , QTC1    , QTC2      )


 END Subroutine ACTUATOR_2
!-----------------------------------------------------------------------
 Subroutine Newmark_solve1 (AM  , AC   , AK   , AQ , DT, &
                            QTCP, QTCP1, QTCP2         , &
                            QTC , QTC1 , QTC2             )
!-----------------------------------------------------------------------

   implicit none

   real(8), intent(in ) :: DT
   real(8), intent(in ) :: AM  , AC   , AK   , AQ
   real(8), intent(in ) :: QTCP, QTCP1, QTCP2
   real(8), intent(out) :: QTC , QTC1 , QTC2

   real(8) :: BITA, GAMA, c1,c2,c3,c4
   real(8) :: AKEQ, UPRE, UTPRE, RQ, U


   BITA  = 0.25d0
   GAMA  = 0.50d0
   c1    =(0.50d0 - BITA)*DT**2
   c2    =(1.00d0 - GAMA)*DT
   c3    = 1.00d0 /(BITA*DT**2)
   c4    = GAMA   /(BITA*DT   )

   AKEQ  = AM*c3     + AC*c4     + AK
   UPRE  = QTCP      + QTCP1*DT  + QTCP2*c1
   UTPRE = QTCP1     +             QTCP2*c2
   RQ    = AQ + ( AM*c3 + AC*c4 )*UPRE - AC * UTPRE
   U     = RQ/AKEQ
   QTC   = U
   QTC1  = UTPRE + c4*(U-UPRE)
   QTC2  =         c3*(U-UPRE)


 END Subroutine Newmark_solve1
!-----------------------------------------------------------------------
 Subroutine Newmark_saturate1 ( DT                        , &
                                MinVAR  , MaxVAR          , &
                                MinDVAR , MaxDVAR         , &
                                MinDDVAR, MaxDDVAR        , &
                                QTCP    , QTCP1   , QTCP2 , &
                                QTC     , QTC1    , QTC2      )
!-----------------------------------------------------------------------

   implicit none

   real(8), intent(in   ) :: DT
   real(8), intent(in   ) :: MinVAR, MaxVAR, MinDVAR, MaxDVAR, MinDDVAR, MaxDDVAR
   real(8), intent(in   ) :: QTCP  , QTCP1 , QTCP2
   real(8), intent(inout) :: QTC   , QTC1  , QTC2

   real(8) :: UPRE, U1PRE, BITA, GAMA, c1,c2


      BITA  = 0.25d0
      GAMA  = 0.50d0
      c1    =(0.50d0 - BITA)*DT**2
      c2    =(1.00d0 - GAMA)*DT
!--- Satate acceleration (DDVAR)
   if (QTC2 < MinDDVAR .or. QTC2 > MaxDDVAR ) then
      UPRE  =  QTCP  + QTCP1*DT + QTCP2*c1
      U1PRE =          QTCP1    + QTCP2*c2
      QTC2  = min ( max( QTC2, MinDDVAR ), MaxDDVAR )
      QTC   = UPRE  + BITA * DT**2 * QTC2
      QTC1  = U1PRE + GAMA * DT    * QTC2
   endif

!--- Saturate velocity (DVAR)
   if (QTC1 < MinDVAR .or. QTC1 > MaxDVAR ) then
      UPRE  =  QTCP  + QTCP1*DT + QTCP2*c1
      U1PRE =          QTCP1    + QTCP2*c2
      QTC1  = min ( max( QTC1, MinDVAR  ), MaxDVAR )
      QTC2  = (QTC1-U1PRE) / (GAMA*DT)
      QTC   = UPRE  + BITA  * DT**2 * QTC2
   endif

!--- Saturate position (VAR)
   if (QTC < MinVAR .or. QTC > MaxVAR ) then
      UPRE  =  QTCP  + QTCP1*DT + QTCP2*c1
      U1PRE =          QTCP1    + QTCP2*c2
      QTC   = min ( max( QTC , MinVAR   ), MaxVAR )
      QTC2  = (QTC-UPRE) / (BITA * DT**2)
      QTC1  = U1PRE + GAMA * DT    * QTC2
   endif


 END Subroutine Newmark_saturate1
!----------------------------------------------------------------------
 Subroutine filterN ( A,B,C,D, QTCP,QTCP1,QTCP2 ,QTC,QTC1,QTC2, DT, N, XX,XX1, YY,YY1 )
!----------------------------------------------------------------------
! [A,B,C,D]=func(n,db,[w1,w2],type,'s')
!  func: cheby1, cheby2, butter, ellip
!  n:order (i.e.2)
!  w1,w2 frequency range [rad/s]
!  type='low' for lowpass, 'high' for highpass,'bandpass' for bandpass and 'stop' for bandstop
!  's' for analog filter
!----------------------------------------------------------------------
!---- from lowpass wind
!   N           = 2
!   damp        = 0.20d0 !0.20 default
!   A (1:2,1:2) = 0.d0;   A (1,1) =-damp;    A (2,1) = 1.00d0
!   B (1:2    ) = 0.d0;   B (1  ) =-damp**2; B (2  ) = damp
!   C (1:2    ) = 0.d0;   C (1  ) = 0.00d0;  C (2  ) = 1.00d0
!   D           = 0.d0
!----------------------------------------------------------------------

   implicit none

   integer, intent(in ) :: N
   real(8), intent(in ) :: A (N,N), B  (N), C  (N), D
   real(8), intent(in ) :: QTCP(N), QTCP1(N), QTCP2(N)
   real(8), intent(in ) :: DT
   real(8), intent(in ) :: XX     , XX1    !, XX2
   real(8), intent(out) :: QTC (N), QTC1 (N), QTC2 (N)
   real(8), intent(out) :: YY     , YY1    !, YY2

   real(8)  ::  AM(N,N), AC(N,N), AK(N,N), AQ(N)
   integer  ::  i


   AK(1:N,1:N) = 0.d0;
   AC(1:N,1:N) =-A(1:N,1:N)
   AM(1:N,1:N) = 0.d0; forall(i = 1:N) AM(i,i) = 1.d0
   AQ(1:N    ) = B(1:N)*XX

   call NEWMARKN (QTCP,QTCP1,QTCP2, QTC,QTC1,QTC2 ,DT, AM,AC,AK,AQ, N)

   YY  = dot_product (C(1:N), QTC1(1:N)) + D*XX   !value
   YY1 = dot_product (C(1:N), QTC2(1:N)) + D*XX1  !1st time derivative


 END Subroutine filterN
!----------------------------------------------------------------------
 Subroutine NEWMARKN (QTCP,QTCP1,QTCP2, QTC,QTC1,QTC2, DT, AM,AC,AK,AQ, N) 
!----------------------------------------------------------------------

   implicit none

   real(8), intent(in ) :: DT
   integer, intent(in ) :: N
   real(8), intent(in ) :: AM(N,N), AC(N,N) , AK(N,N), AQ(N)
   real(8), intent(in ) :: QTCP(N), QTCP1(N), QTCP2(N)
   real(8), intent(out) :: QTC (N), QTC1 (N), QTC2 (N)

   integer :: ierr
   real(8) :: RQ(N), AKEQ(N,N)
   integer :: INDXO(N)
   real(8) :: BITA, GAMA, c1,c2,c3,c4
   real(8) :: UPRE(N), UTPRE(N), d


   BITA           = 0.25d0
   GAMA           = 0.50d0
   c1             =(0.50d0 - BITA)*DT**2
   c2             =(1.00d0 - GAMA)*DT
   c3             = 1.00d0 /(BITA*DT**2)
   c4             = GAMA   /(BITA*DT   )

   UPRE (1:N    ) = QTCP  (1:N    )    + QTCP1(1:N)*DT + QTCP2(1:N)*c1
   UTPRE(1:N    ) = QTCP1 (1:N    )    +                 QTCP2(1:N)*c2
   AKEQ (1:N,1:N) = AM    (1:N,1:N)*c3 +         AC (1:N,1:N)*c4 + AK (1:N,1:N)
   RQ   (1:N    ) = AQ    (1:N    )    + matmul( AM (1:N,1:N)*c3 + AC (1:N,1:N)*c4, UPRE (1:N) ) &
                                       - matmul( AC (1:N,1:N)                     , UTPRE(1:N) )

   call LUDCMP (AKEQ, N, N, INDXO, d , ierr); if (ierr==1) then;write(*,*)'LUDCMP in NEWMARKN error',N;stop;endif
   call LUBKSB (AKEQ, N, N, INDXO, RQ      )

   QTC  (1:N    ) = RQ   (1:N)
   QTC1 (1:N    ) = UTPRE(1:N) + c4*(RQ(1:N)-UPRE(1:N))
   QTC2 (1:N    ) =              c3*(RQ(1:N)-UPRE(1:N))


 END Subroutine NEWMARKN
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
