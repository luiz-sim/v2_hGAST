!----------------------------------------------------------------------
 Module Cinwind
!----------------------------------------------------------------------
   integer,              save :: NREC0     , NTOT
   integer,              save :: NGRIDP    , NGRIDP1      , NGRIDC
   integer,              save :: ITYPU     , ITYPV        , ITYPW
   real(8),              save :: TREC0     , DTINF        , VELHUBI
   real(8),              save :: TIME_turb ,  hub_velocity
   real(8), allocatable, save :: UTSR (:,:), VTSR (:,:)   , WTSR (:,:) !(NREC0,NTOT)
   real(8), allocatable, save :: TREC (:  )                            !(NREC0  )
   real(8), allocatable, save :: RINW (:  )                            !(NGRIDC )
   real(8), allocatable, save :: TINW (:  )                            !(NGRIDP1)

 END Module Cinwind
!-----------------------------------------------------------------------
!
!   Wind module: Turbulent wind is used when ITURB=1
!
!                TREC0 = tip_radius/hub_velocity
!                NREC0 = int(2.d0*TREC0/DTINF) + 2
!                TIME1 = TIME + TREC0 - X0(1)/VELHUBT,  X0=0 at Hub
!
!   The Turbulent wind files are given wrt the global c.s, at hub-hight.
!   The disk NYxNZ loop first run the y direction from negative to positive
!   and after the z direction fron negative to positive. So the n = j+(k-1)*NY
!
!   Subroutine: INWINDINT 
!   Initializes code to read data from INWIND
! 
!-----------------------------------------------------------------------
 Subroutine INWINDINT ( windpath, hub_velocity_in, TIME_turb_in, ITYPE, SIG )
!-----------------------------------------------------------------------

 use Cinwind
#ifdef HAVE_MPI
 use MPI
#endif

   Implicit none


   character*256, intent(in) :: windpath
   real(8)      , intent(in) :: hub_velocity_in, TIME_turb_in, SIG(3)
   integer      , intent(in) :: ITYPE

   real(4), allocatable :: RUTS (:), RVTS (:), RWTS (:)
   real(4), allocatable :: RXINW(:), RYINW(:), RZINW(:)
   real(8), allocatable :: XINW (:), YINW (:), ZINW (:)
   real(4)              :: RVELHUB , RFMAX
   real(4)              :: RXLAMBDA, RXMIN, RXMAX, RZROUGH, RRADTIP
   real(4)              :: RZHUB
   real(8)              :: FMAX, DTHET, ZHUBINW
   real(8)              :: PI2
   integer              :: NPTX, i,n,j,k,l
!--- for statistics
   real(8), allocatable :: UU(:), VV(:), WW(:)
   real(8)              :: UMEAN(3),USTD(3),UMIN(3),UMAX(3), SUMD(3),SUMM(3)
   integer              :: ipas
!--- mpi vars
   integer              :: my_rank
#ifdef HAVE_MPI
   integer              :: ierr, nf, mat1, mat2
#endif


#ifdef HAVE_MPI
   call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)
#else
   my_rank = 0
#endif

if (my_rank.eq.0) then
   hub_velocity = hub_velocity_in
   TIME_turb    = TIME_turb_in
   PI2          = dacos(-1.d0)*2.d0

!--- Open Turbulent inflow file
   open (99,file=trim(windpath),form='UNFORMATTED')

! usually Nperiferial > Nradial
!         NZ          <= NY

   read(99)    &
     NPTX,     & !num of components (usually 1)
     NGRIDP,   & !NZ  -  periferial points
     NGRIDC,   & !NY  -  radial     points
     NREC0 ,   & !NX or NT
     RXLAMBDA, &    !Dummy, set zero
     RXMIN,    &    !Dummy, set zero
     RXMAX,    &    !Dummy, set zero
     RZROUGH,  &    !Dummy, set zero
     RRADTIP,  &    !Dummy, set zero
     RZHUB,    &    !Dummy, set zero
     RVELHUB,  & !VELHUB
     RFMAX,    & !1/(2*DT_wind)
     NTOT,     & !NY*NZ
     ITYPU,    & !0 or 1
     ITYPV,    & !0 or 1
     ITYPW       !0 or 1

   Allocate ( RUTS (NTOT), RVTS (NTOT), RWTS (NTOT) )
   Allocate ( RXINW(NTOT), RYINW(NTOT), RZINW(NTOT) )
   Allocate ( XINW (NTOT), YINW (NTOT), ZINW (NTOT) )

   NGRIDP1  = NGRIDP + 1
   ZHUBINW  = dble(RZHUB  )
   VELHUBI  = dble(RVELHUB)
   FMAX     = dble(RFMAX  )
   DTINF    = 1.d0/(2.d0*FMAX)
   TREC0    = dble(NREC0-1)*DTINF

   Allocate ( UTSR (NREC0,NTOT) )
   Allocate ( VTSR (NREC0,NTOT) )
   Allocate ( WTSR (NREC0,NTOT) )
   Allocate ( TREC (NREC0     ) )
   Allocate ( RINW (NGRIDC    ) )
   Allocate ( TINW (NGRIDP1   ) )

   write (10,*)
   write (10,*)
   write (10,*) 'NGRIDP   ' , NGRIDP
   write (10,*) 'NGRIDC   ' , NGRIDC
   write (10,*) 'NREC0    ' , NREC0
   write (10,*) 'RZHUB    ' , RZHUB
   write (10,*) 'RVELHUB  ' , RVELHUB
   write (10,*) 'RFMAX    ' , RFMAX
   write (10,*) 'NTOT     ' , NTOT
   write (10,*) 'ITYPU    ' , ITYPU
   write (10,*) 'ITYPV    ' , ITYPV
   write (10,*) 'ITYPW    ' , ITYPW
   write (10,*)
   write (10,*) 'TIME_turb' , TIME_turb
   write (10,*) 'VELHUBI  ' , VELHUBI
   write (10,*)
   write (10,*) 'TREC0    ' , TREC0
   write (10,*) 'NREC0    ' , NREC0
   write (10,*) 'DTINF    ' , DTINF


   read(99)   (RXINW(i),RYINW(i),RZINW(i),i=1,NTOT)

   open (77,file='velgrd.dat')
    do i       = 1, NTOT
       XINW(i) = dble(RXINW(i))
       YINW(i) = dble(RYINW(i))
       ZINW(i) = dble(RZINW(i))
       write (77,'(3f15.5)') XINW(i),YINW(i),ZINW(i)
    enddo
   close  (77)


!--- Read turbulent box velocities
         UTSR      = 0.d0;
         VTSR      = 0.d0;
         WTSR      = 0.d0;
   do    n         = 1, NREC0
      if (ITYPU == 1) read(99) (RUTS(i),i=1,NTOT)
      if (ITYPV == 1) read(99) (RVTS(i),i=1,NTOT)
      if (ITYPW == 1) read(99) (RWTS(i),i=1,NTOT)

         TREC(n  ) = dble(n-1)*DTINF
      do i         = 1, NTOT
         UTSR(n,i) = (dble(RUTS(i)) - VELHUBI)
         VTSR(n,i) = (dble(RVTS(i))          )
         WTSR(n,i) = (dble(RWTS(i))          )
      enddo
   enddo
   close (99)


!--- Calculate turbulent box statictics (std, mean)
   Allocate (UU(NREC0),VV(NREC0),WW(NREC0))
   open(2,file='turbbox_stats.dat')
   do      ipas    = 1, 2
           SUMD(:) = 0.0
           SUMM(:) = 0.0
      do   i       = 1, NTOT
        do n       = 1, NREC0
           UU(n)   = UTSR(n,i)
           VV(n)   = VTSR(n,i)
           WW(n)   = WTSR(n,i)
        enddo !n
        call STATCO (UU,NREC0, NREC0 ,UMEAN(1),USTD(1),UMIN(1),UMAX(1))
           SUMD(1) = SUMD(1) + USTD (1)/ real(NTOT)                   
           SUMM(1) = SUMM(1) + UMEAN(1)/ real(NTOT)                   
        call STATCO (VV,NREC0, NREC0 ,UMEAN(2),USTD(2),UMIN(2),UMAX(2))
           SUMD(2) = SUMD(2) + USTD (2)/ real(NTOT)                   
           SUMM(2) = SUMM(2) + UMEAN(2)/ real(NTOT)                   
        call STATCO (WW,NREC0, NREC0 ,UMEAN(3),USTD(3),UMIN(3),UMAX(3))
           SUMD(3) = SUMD(3) + USTD (3)/ real(NTOT)
           SUMM(3) = SUMM(3) + UMEAN(3)/ real(NTOT)
         write (2,'(20f12.5)') &
           (USTD (j),j=1,3)   ,&
           (UMEAN(j),j=1,3)   ,&
           (UMIN (j),j=1,3)   ,&
           (UMAX (j),j=1,3)
      enddo !i
         write (2,*)
         write (2,*)

                   write(10,*)
      if (ipas==1) write(10,*) 'original box statistics'
      if (ipas==2) write(10,*) 'modified box statistics'
                   write(10,*) 'umean, ustd: ', SUMM(1), SUMD(1)
                   write(10,*) 'vmean, vstd: ', SUMM(2), SUMD(2)
                   write(10,*) 'wmean, wstd: ', SUMM(3), SUMD(3)
      
      
!------ Scale turbulent velocities if IWINDC=6,7,8
      if (SIG(1)>0.d0.and.SIG(2)>0.d0.and.SIG(3)>0.d0.and.ipas==1) then
            write(10,*)
            write(10,'(a,3f12.5)') ' scale turbulent box velocities', &
                  SIG(1)/SUMD(1),SIG(2)/SUMD(2),SIG(3)/SUMD(3)
         if (SUMD(1)==0.d0.or.SUMD(2)==0.d0.or.SUMD(3)==0.d0) then
            write(*,*)'error in INWINDINT sdv is zero!';stop
         endif
         UTSR(1:NREC0,1:NTOT) = UTSR(1:NREC0,1:NTOT)*SIG(1)/SUMD(1)
         VTSR(1:NREC0,1:NTOT) = VTSR(1:NREC0,1:NTOT)*SIG(2)/SUMD(2)
         WTSR(1:NREC0,1:NTOT) = WTSR(1:NREC0,1:NTOT)*SIG(3)/SUMD(3)
      else
         exit
      endif
   enddo !ipas
   close(2)
   Deallocate (UU,VV,WW)


!--- Define Grid based on ITYPE (Disk, Rectangular)
   if     (ITYPE.eq.1) then
!--- Disk
      DTHET        = PI2/ dble (NGRIDP)
      do i         = 1, NGRIDC
         k         = i*NGRIDP
         RINW(i)   = dsqrt(YINW(k)**2+(ZINW(k)-ZHUBINW)**2)
      enddo
      do j         = 1, NGRIDP
         TINW(j)   = (j-1)*DTHET
      enddo
         j         =    NGRIDP1          !Periodic condition
         TINW(j)   = TINW(NGRIDP)+DTHET
   elseif (ITYPE.eq.2) then
!--- Rectangular
      do i         = 1, NGRIDC !NY
         RINW(i)   = YINW(i)
      enddo
      do j         = 1, NGRIDP !NZ
         l         = (j-1)*NGRIDC+1
         TINW(j)   = ZINW(l)
      enddo
   endif
   Deallocate ( RUTS , RVTS , RWTS  )
   Deallocate ( RXINW, RYINW, RZINW )
   Deallocate ( XINW , YINW , ZINW  )
endif !my_rank=0

#ifdef HAVE_MPI
 call MPI_BCAST(NGRIDP       ,1,MPI_INTEGER          ,0,MPI_COMM_WORLD,ierr)
 call MPI_BCAST(NGRIDP1      ,1,MPI_INTEGER          ,0,MPI_COMM_WORLD,ierr)
 call MPI_BCAST(NGRIDC       ,1,MPI_INTEGER          ,0,MPI_COMM_WORLD,ierr)
 call MPI_BCAST(NREC0        ,1,MPI_INTEGER          ,0,MPI_COMM_WORLD,ierr)
 call MPI_BCAST(ITYPU        ,1,MPI_INTEGER          ,0,MPI_COMM_WORLD,ierr)
 call MPI_BCAST(ITYPV        ,1,MPI_INTEGER          ,0,MPI_COMM_WORLD,ierr)
 call MPI_BCAST(ITYPW        ,1,MPI_INTEGER          ,0,MPI_COMM_WORLD,ierr)
 call MPI_BCAST(NTOT         ,1,MPI_INTEGER          ,0,MPI_COMM_WORLD,ierr)
 call MPI_BCAST(hub_velocity ,1,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,ierr)
 call MPI_BCAST(TIME_turb    ,1,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,ierr)
 call MPI_BCAST(TREC0        ,1,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,ierr)
 call MPI_BCAST(DTINF        ,1,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,ierr)
 
 if (my_rank /= 0) then
   Allocate ( UTSR (NREC0,NTOT) )
   Allocate ( VTSR (NREC0,NTOT) )
   Allocate ( WTSR (NREC0,NTOT) )
   Allocate ( TREC (NREC0     ) )
   Allocate ( RINW (NGRIDC    ) )
   Allocate ( TINW (NGRIDP1   ) )
 endif

 call mpimat2d(mat2,NREC0,NTOT,NREC0,NTOT,0,0)
 call MPI_BCAST(UTSR,1,mat2,0,MPI_COMM_WORLD,ierr)
 call MPI_BCAST(VTSR,1,mat2,0,MPI_COMM_WORLD,ierr)
 call MPI_BCAST(WTSR,1,mat2,0,MPI_COMM_WORLD,ierr)
 call MPI_TYPE_FREE(mat2,ierr)
 
 call mpimat1d(mat1,NREC0,NREC0,0)
 call MPI_BCAST(TREC,1,mat1,0,MPI_COMM_WORLD,ierr)
 call MPI_TYPE_FREE(mat1,ierr)
 
 call mpimat1d(mat1,NGRIDC,NGRIDC,0)
 call MPI_BCAST(RINW,1,mat1,0,MPI_COMM_WORLD,ierr)
 call MPI_TYPE_FREE(mat1,ierr)
 
 call mpimat1d(mat1,NGRIDP1,NGRIDP1,0)
 call MPI_BCAST(TINW,1,mat1,0,MPI_COMM_WORLD,ierr)
 call MPI_TYPE_FREE(mat1,ierr)

 write(*,*) 'my_rank:',my_rank,'Concluding INWINDINT'
#endif



 END Subroutine INWINDINT
!--------------------------------------------------------------------------
! INTERPOLATION OF THE VELOCITIES WITH THE STOCHASTIC PART
!
!     INTERPOLATION OF THE VELOCITIES AT X
!            ---- > VELINWIND   
!   
!     CALCUL OF THE DEFORAMTION 
!            ---- > DEFINWIND     
!
!  SUBROUTINE of the interpolation of the stochastic wind with
!   Polar grid ....
!
!--------------------------------------------------------------------------
 Subroutine VELINWIND (X,U,TIME,ITYPE)
!--------------------------------------------------------------------------

 use Cinwind

   implicit none

   integer, intent(in ) :: ITYPE
   real(8), intent(in ) :: X(3), TIME
   real(8), intent(out) :: U(3)

   integer :: NR1,NR2,I,J,K,NI,NJ
   real(8) :: Time_ramp,TIME1,RAT,R,T
   real(8) :: U2DTS (NGRIDC, NGRIDP1)
   real(8) :: V2DTS (NGRIDC, NGRIDP1)
   real(8) :: W2DTS (NGRIDC, NGRIDP1)


!  if ( TIME < TIME_turb ) return

   TIME1 = TIME  - X(1)/hub_velocity - TIME_turb
!1 if (TIME1< 0.d0) then; TIME1=TIME1+TREC0; goto 1; endif;
!2 if (TIME1>TREC0) then; TIME1=TIME1-TREC0; goto 2; endif;

   if (TIME1< 0.d0) then; U(1:3)=0.d0;       return; endif;
 2 if (TIME1>TREC0) then; TIME1=TIME1-TREC0; goto 2; endif;


      ! LOCATE(NM   ,N    ,A(NM),X    ,I  )
   call LOCATE(NREC0,NREC0,TREC ,TIME1,NR1)

   NR2 = NR1 + 1
   RAT = (TIME1-TREC(NR1))/(TREC(NR2)-TREC(NR1))


   if (ITYPE.eq.1) then
!--- Disk: Polar coordinates (r,T)
      do I          = 1, NGRIDC       ! radial
      do J          = 1, NGRIDP       ! cyrcumferential
         K          = J + (I-1)*NGRIDP
         U2DTS(I,J) = (1.d0-RAT)*UTSR(NR1,K) + RAT*UTSR(NR2,K)
         V2DTS(I,J) = (1.d0-RAT)*VTSR(NR1,K) + RAT*VTSR(NR2,K)
         W2DTS(I,J) = (1.d0-RAT)*WTSR(NR1,K) + RAT*WTSR(NR2,K)
      enddo
      enddo
         J          =    NGRIDP1      ! Periodic condition
      do I          = 1, NGRIDC
         U2DTS(I,J) = U2DTS(I,1)
         V2DTS(I,J) = V2DTS(I,1)
         W2DTS(I,J) = W2DTS(I,1)
      enddo
         NI         = NGRIDC
         NJ         = NGRIDP1
         R          = dsqrt(X(2)**2+X(3)**2)
         T          = datan2(X(3),X(2))
      if (R < RINW(1)) R = RINW(1)
      if (T < 0.d0   ) T = T+2d0*dacos(-1.d0)

      if (R > RINW(NGRIDC)) then
         write (*,'(a,5f15.5)') "POINT OUT OF GRID2",  TIME1,X(2),X(3),R,RINW(NGRIDC)
         U(1:3) = 0.d0;
         return
      endif

   elseif (ITYPE.eq.2) then
!--- Rectangular grid - plane (y,z)
      do K          = 1, NGRIDP         !NZ
      do J          = 1, NGRIDC         !NY
         I          = J + (K-1)*NGRIDC  !NYZ
         U2DTS(J,K) = RAT*UTSR(NR2,I) + (1.d0-RAT)*UTSR(NR1,I)
         V2DTS(J,K) = RAT*VTSR(NR2,I) + (1.d0-RAT)*VTSR(NR1,I)
         W2DTS(J,K) = RAT*WTSR(NR2,I) + (1.d0-RAT)*WTSR(NR1,I)
      enddo
      enddo
         NI         = NGRIDC            !NY
         NJ         = NGRIDP            !NZ
         R          = X(2)              !Y
         T          = X(3)              !Z

      if ((R < RINW(1)).or.(R > RINW(NGRIDC)).or.&
          (T < TINW(1)).or.(T > TINW(NGRIDP))      ) then
          write (*,'(a,5f15.5)') "POINT OUT OF GRID3",  TIME1,X(2),X(3)
          write (*,'(a,5f15.5)') "R,R1,Rn",R,RINW(1), RINW(NGRIDC)
          write (*,'(a,5f15.5)') "T,T1,Tn",T,TINW(1), TINW(NGRIDC)
          U(1:3) = 0.d0;
          return
      endif
   endif !ITYPE

!--- linear interpolation
   if (ITYPU.EQ.1) call BLIN  (RINW,TINW,U2DTS,R,T,U(1),NI,NJ,NGRIDC,NGRIDP1)
   if (ITYPV.EQ.1) call BLIN  (RINW,TINW,V2DTS,R,T,U(2),NI,NJ,NGRIDC,NGRIDP1)
   if (ITYPW.EQ.1) call BLIN  (RINW,TINW,W2DTS,R,T,U(3),NI,NJ,NGRIDC,NGRIDP1)
!--- spline interpolation
!!!if (ITYPU.EQ.1) call SPLIN (RINW,TINW,U2DTS,R,T,U(1),NI,NJ,NGRIDC,NGRIDP1)
!!!if (ITYPV.EQ.1) call SPLIN (RINW,TINW,V2DTS,R,T,U(2),NI,NJ,NGRIDC,NGRIDP1)
!!!if (ITYPW.EQ.1) call SPLIN (RINW,TINW,W2DTS,R,T,U(3),NI,NJ,NGRIDC,NGRIDP1)


!--- Initially linear ramp of turbulent velocities for 10sec
       Time_ramp = 10.d0

   if (TIME <= Time_ramp) then
       U(1:3) = U(1:3)*TIME/Time_ramp
   endif


 END Subroutine VELINWIND
!
!
!
!----------------------------------------------------------------------
! Performs linear interpolation of a 2variable tabulated function
! YA(MM,NM) given at discrite points X1A(MM),X2A(NM).
! The interpolation is performed @ points X1,X2 and the 
! output value is Y.
! M, N are the number of the given values (M<=MM, N<=NM)
!----------------------------------------------------------------------
 Subroutine BLIN (X1A,X2A,YA,X1,X2,Y,M,N,MM,NM)
!----------------------------------------------------------------------

   implicit none

   integer, intent( in) :: M,N, MM,NM
   real(8), intent( in) :: X1A(MM),X2A(NM),YA(MM,NM)
   real(8), intent( in) :: X1, X2
   real(8), intent(out) :: Y

   real(8) :: GKSI, GHTA !contribution of each nearest point
   integer :: M1  , N1   !result from LOCATE


   call LOCATE(MM,M,X1A,X1,M1)
   call LOCATE(NM,N,X2A,X2,N1)

   GKSI = (X1-X1A(M1))/(X1A(M1+1)-X1A(M1))
   GHTA = (X2-X2A(N1))/(X2A(N1+1)-X2A(N1))

   Y = (1.d0-GKSI)*(1.d0-GHTA)*YA(M1  ,N1  ) &
      +      GKSI *(1.d0-GHTA)*YA(M1+1,N1  ) &
      +      GKSI *      GHTA *YA(M1+1,N1+1) &
      +(1.d0-GKSI)*      GHTA *YA(M1  ,N1+1)


 END Subroutine BLIN
!----------------------------------------------------------------------
 Subroutine SPLIN (X1A,X2A,YA,X1,X2,Y,M,N,MM,NM)
!----------------------------------------------------------------------

   implicit none

   real(8) :: X1A(MM),X2A(NM),YA(MM,NM), &
              YTMP(NM),Y2TMP(NM),YYTMP(NM),X1,X2,Y
   integer :: J,M,N,MM,NM


   do J         = 1, M
      YTMP(1:N) = YA(J,1:N)
      CALL SPLINE(X2A,YTMP,Y2TMP,            N,NM)
      CALL SPLINT(X2A,YTMP,Y2TMP,X2,YYTMP(J),N,NM)
   enddo
      CALL SPLINE(X1A,YYTMP,Y2TMP,           M,MM)
      CALL SPLINT(X1A,YYTMP,Y2TMP,X1,Y      ,M,MM)


 END Subroutine SPLIN
!----------------------------------------------------------------------
 Subroutine SPLINT (XA,YA,Y2A,X,Y,N,NM)
!----------------------------------------------------------------------

   implicit none

   real(8) :: XA(NM),YA(NM),Y2A(NM),X,Y,H,A,B
   integer :: N,NM,I,KLO,KHI


   call LOCATE(NM,N,XA,X,I)
   KLO = I
   KHI = I+1
   H   = (XA(KHI)-XA(KLO))
   A   = (XA(KHI)-X)/H
   B   = (X-XA(KLO))/H
   Y   = A*YA(KLO)+B*YA(KHI)+((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.d0


 END Subroutine SPLINT
!----------------------------------------------------------------------
 Subroutine SPLINE (X,Y,Y2,N,NM)
!----------------------------------------------------------------------

   implicit none

   real(8) :: X(NM),Y(NM),Y2(NM),U(NM),SIG,P,QN,UN
   integer :: N,NM,I,K


!  if (NM.LT.N) STOP
      Y2(1) = 0.d0
      U (1) = 0.d0
   do I     = 2, N-1
      SIG   = (X(I)-X(I-1))/(X(I+1)-X(I))
      P     = SIG*Y2(I-1)+2.d0
      Y2(I) = (SIG-1.d0)/P
      U (I) = (6.d0*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))  &
                 /(X(I)-X(I-1)))/(X(I+1)-X(i-1))-SIG*U(I-1))/P
   enddo !I
      QN    = 0.d0
      UN    = 0.d0
      Y2(N) = (UN-QN*U(N-1))/(QN*Y2(N-1)+1.d0)
   do K     = N-1, 1, -1
      Y2(K) = Y2(K)*Y2(K+1)+U(K)
   enddo !K


 END Subroutine SPLINE
!----------------------------------------------------------------------
 Subroutine STATCO (X,NP,N,XMEAN,STD,XMIN,XMAX)
!----------------------------------------------------------------------
   implicit none

   integer, intent(in ) :: NP, N
   real(8), intent(in ) :: X(NP)
   real(8), intent(out) :: XMEAN,STD,XMIN,XMAX

   integer :: I
   real(8) :: SUM


      XMEAN = X(1)
      STD   = 0.d0;
      XMIN  = X(1)
      XMAX  = X(1)
   if (N < 2) return
      SUM   = 0.d0
   do I     = 1, N
      IF (XMAX.LT.X(I)) XMAX=X(I)
      IF (XMIN.GT.X(I)) XMIN=X(I)
      SUM   = SUM + X(I)
   enddo
      XMEAN = SUM/dble(N)
      SUM   = 0.d0
   do I     = 1, N
      SUM   = SUM+(X(I)-XMEAN)**2
   enddo
      STD   =(SUM/dble(N))**0.5d0


 END Subroutine STATCO

#ifdef HAVE_MPI
!--------------------------------------------------------------------------
 Subroutine mpimat2d (mat2d,orig1,orig2,nsize1,nsize2,istart1,istart2)
!--------------------------------------------------------------------------

 use MPI

   implicit none

   integer :: ierr
   integer :: typelist(2)
   integer :: imat(2),mat(2),start(2)
   integer :: istart1,istart2
   integer :: orig1,orig2,nsize1,nsize2,mat2d


   imat    (1) = orig1
   imat    (2) = orig2
   mat     (1) = nsize1
   mat     (2) = nsize2
   start   (1) = istart1
   start   (2) = istart2
   typelist(1) = MPI_DOUBLE_PRECISION
   typelist(2) = MPI_DOUBLE_PRECISION

   call MPI_TYPE_CREATE_SUBARRAY(2,imat,mat,start,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,mat2d,ierr)
   call MPI_TYPE_COMMIT(mat2d,ierr)


 END subroutine mpimat2d
!--------------------------------------------------------------------------
 Subroutine mpimat1d (mat1d,orig,nsize,istart)
!--------------------------------------------------------------------------

 use MPI

   implicit none

   integer :: ierr
   integer :: typelist(1)
   integer :: imat(1),mat(1),start(1)
   integer :: istart
   integer :: orig,nsize,mat1d


   imat    (1) = orig
   mat     (1) = nsize
   start   (1) = istart
   typelist(1) = MPI_DOUBLE_PRECISION

   call MPI_TYPE_CREATE_SUBARRAY(1,imat,mat,start,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,mat1d,ierr)
   call MPI_TYPE_COMMIT(mat1d,ierr)
 

 END subroutine mpimat1d
#endif
