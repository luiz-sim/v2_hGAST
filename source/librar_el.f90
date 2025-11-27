!----------------------------------------------------------------------
 Subroutine DIAGO (NM,A)
!----------------------------------------------------------------------

   implicit none

   integer, intent(in)  :: NM
   real(8), intent(out) :: A(NM,NM)

   integer :: i


      A      = 0.d0;
   do i      = 1, NM
      A(i,i) = 1.d0
   enddo


 END Subroutine DIAGO
!-----------------------------------------------------------------------
!
!     Subroutine :LOCATE(NM,N,A,X,I)
!
!
!     LOCATE locates the interval of a discrete sequence of points
!            A[1:N] within which X falls
!
!     N = maximum legth of A[]
!     I : A(I) < X < A(I+1)
!
!-----------------------------------------------------------------------
 Subroutine LOCATE (nm,n,xa,x,k)
!-----------------------------------------------------------------------

   implicit none

   integer, intent(in ) :: nm,n
   real(8), intent(in ) :: xa(nm), x
   integer, intent(out) :: k

   integer :: khi, klo


   if ((x.lt.xa(1)).or.(x.gt.xa(n))) then
      write (* ,*) 'Problem in sub.LOCATE !!'
      write (* ,*) 'x =',x, xa(1), xa(n), n
      write (10,*) 'Problem in sub.LOCATE !!'
      write (10,*) 'x =',x, xa(1), xa(n), n
      stop
   endif

!   khi = 1
!11 khi = khi+1
!   if (x > xa(khi)) goto 11
!   k   = khi-1

   klo = 1
   khi = n
   if (xa(1) < xa(n)) then
 1    if (khi-klo > 1) then
         k=(khi+klo)/2
         if (xa(k) > x) then
            khi=k
         else
            klo=k
         endif
         goto 1
      endif
   else
 2    if (khi-klo > 1) then
         k=(khi+klo)/2
         if (xa(k) > x) then
            klo=k
         else
            khi=k
         endif
         goto 2
      endif
   endif

   k   = klo


 END Subroutine LOCATE
!--------------------------------------------------------------------
 Subroutine ROTATEP ( IAXIS, TH, XINIT, XFINAL )
!--------------------------------------------------------------------

   implicit none

   real(8), intent(in ) :: TH, XINIT(3)
   integer, intent(in ) :: IAXIS
   real(8), intent(out) :: XFINAL(3)

   real(8) :: COSTH, SINTH, Xtmp(3)


   COSTH = dcos(TH)
   SINTH = dsin(TH)

   goto (1,2,3), IAXIS

!--- Rotation about X axis
 1 Xtmp  (1) =  XINIT(1)
   Xtmp  (2) =  XINIT(2)*COSTH - XINIT(3)*SINTH
   Xtmp  (3) =  XINIT(2)*SINTH + XINIT(3)*COSTH
   XFINAL(:) =  Xtmp (:)
   return

!--- Rotation about Y axis
 2 Xtmp  (1) =  XINIT(1)*COSTH + XINIT(3)*SINTH
   Xtmp  (2) =  XINIT(2)
   Xtmp  (3) = -XINIT(1)*SINTH + XINIT(3)*COSTH
   XFINAL(:) =  Xtmp (:)
   return

!--- Rotation about Z axis
 3 Xtmp  (1) =  XINIT(1)*COSTH - XINIT(2)*SINTH
   Xtmp  (2) =  XINIT(1)*SINTH + XINIT(2)*COSTH
   Xtmp  (3) =  XINIT(3)
   XFINAL(:) =  Xtmp (:)


 END Subroutine ROTATEP
!-----------------------------------------------------------
 Subroutine LIN_INT(x,y,xi,yi, n1, nnm)
!-----------------------------------------------------------

   implicit none

   integer, intent(in ) :: n1, nnm
   real(8), intent(in ) :: x, xi(nnm), yi(nnm)
   real(8), intent(out) :: y

   real(8) :: sl
   integer :: i

  
   if (x.lt.xi(1)) then
      sl =(yi(2)-yi(1))/(xi(2)-xi(1)+1.d-17)
      y  = yi(1)-sl*(x-xi(1))
      return
   endif

   if (x.gt.xi(n1)) then
      sl =(yi(n1)-yi(n1-1))/(xi(n1)-xi(n1-1)+1.d-17)
      y  = yi(n1-1)+sl*(x-xi(n1-1))
      return
   endif

   do i=1,n1-1
      if ( (x.ge.xi(i)).and.(x.le.xi(i+1)) ) then
         y = yi(i) + (yi(i+1)-yi(i))/(xi(i+1)-xi(i))*(x-xi(i))
         return
      endif
   enddo


 END Subroutine LIN_INT
!----------------------------------------------------------------------
 Subroutine INT_2_CHAR ( NP, CNUM, N )
!----------------------------------------------------------------------

! Transform the integer N of NP digits max to character CNUM(NP),
! for opening files

   implicit none

   integer  , intent(in ) :: NP, N
   character, intent(out) :: CNUM(NP)

   integer :: i, i0, i1


      i0      = N
   do i       = 1, NP
      i1      = i0/10**(NP-i)
      i0      = i0-i1*10**(NP-i)
      CNUM(i) = char(48+i1)
   enddo


 END Subroutine INT_2_CHAR
!----------------------------------------------------------------------
 SUBROUTINE spline1 ( x, y, y2, n )
!----------------------------------------------------------------------
!
! Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e.
! y(i)=f(x(i)), with x(1)<x(2)<...<x(n), and given values yp1 and ypn for
! the first derivative of the interpolating function at points 1 and n,
! respectively, this routine returns an array y2(1:n) of length n which
! contains the second derivatives of the interpolating function at the
! tabulated points x(i).  If yp1 and/or ypn are equal to 1.e30 or larger,
! the routine is signaled to set the corresponding boundary condition for
! a natural spline with zero second derivative on that boundary.
! (adopted from Numerical Recipes in FORTRAN 77)
!
   implicit none

   INTEGER, PARAMETER :: DP=KIND(1.0D0)
   INTEGER:: n
   REAL(DP):: x(n), y(n), y2(n) !,yp1, ypn, 
   INTEGER:: i, k
   REAL(DP):: p, qn, sig, un, u(n)

!    if (yp1.gt..99e30) then
        y2(1)=0.d0
        u (1)=0.d0
!    else
!       y2(1)=-0.5
!       u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
!    endif

     do i=2, n-1
        sig   = (x(i)-x(i-1))/(x(i+1)-x(i-1))
        p     = sig*y2(i-1)+2.d0
        y2(i) = (sig-1.d0)/p
        u (i) = (6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/&
                (x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
     enddo

!    if (ypn.gt..99e30) then
        qn=0.d0
        un=0.d0
!    else
!       qn=0.5
!       un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
!    endif

     y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)

     do k=n-1, 1, -1
        y2(k)=y2(k)*y2(k+1)+u(k)
     enddo


 END Subroutine spline1
!----------------------------------------------------------------------
 SUBROUTINE splint1 ( xa, ya, y2a, n, x, y )
!----------------------------------------------------------------------
!
! Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function
! (with the xa(i) in order), and given the array y2a(1:n), which is the output
! from the subroutine spline, and given a value of x, this routine returns a
! cubic spline interpolated value y.
! (adopted from Numerical Recipes in FORTRAN 77)
!
   implicit none

   INTEGER, PARAMETER :: DP = KIND(1.0D0)
   INTEGER:: n
   REAL(DP):: x, y, xa(n), y2a(n), ya(n)
   INTEGER:: k, khi, klo
   REAL(DP):: a, b, h

     klo=1
     khi=n
1   if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if (xa(k).gt.x) then
           khi=k
        else
           klo=k
        endif
        goto 1
     endif

     h=xa(khi)-xa(klo)
     if (h == 0.d0) then; write(*,*)'bad xa input in splint'; stop; endif

     a=(xa(khi)-x)/h
     b=(x-xa(klo))/h
     y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0


 END Subroutine splint1
!----------------------------------------------------------------------
 SUBROUTINE lint1 ( xa, ya, n1, n2, np, x, y )
!----------------------------------------------------------------------
!
! Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function
! (with the xa(i) in order), and given the array y2a(1:n), which is the output
! from the subroutine spline, and given a value of x, this routine returns a
! linear interpolated value y.
!
   implicit none

   integer :: n1, n2, np
   real(8) :: x, y, xa(np), ya(np)
   integer :: k, khi, klo
   real(8) :: a, b, h

     klo=n1
     khi=n2
1   if (khi-klo > 1) then
        k=(khi+klo)/2
        if (xa(k) > x) then
           khi=k
        else
           klo=k
        endif
        goto 1
     endif

     h=xa(khi)-xa(klo)
     if (h == 0.d0) then; write(*,*)'bad xa input in lint1'; stop; endif

     b=(x-xa(klo))/h
     a=1.d0-b              !(xa(khi)-x)/h
     y=a*ya(klo)+b*ya(khi) !+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0


 END Subroutine lint1
!----------------------------------------------------------------------
 Subroutine EXTEPR ( A, B, C )
!----------------------------------------------------------------------

   implicit none

   real(8), intent(in ) :: A(3), B(3)
   real(8), intent(out) :: C(3)

!- External product of two 3D vectors c = a x b

   C(1) =  A(2)*B(3) - A(3)*B(2)
   C(2) = -A(1)*B(3) + A(3)*B(1)
   C(3) =  A(1)*B(2) - A(2)*B(1)


 END Subroutine EXTEPR
!-----------------------------------------------------------
 Subroutine Runinfo (IGO)
!-----------------------------------------------------------
! This subroutine outputs the:
! 1. given title
! 2. git version
! 3. date in the form "mmm dd ccyy" 
! 4. time in the form "hh:mm:ss"
! 5. date-time, paths and omp threads
! 6. exe        path
! 7. simulation path
! 8. OMP threads
!
! for git version the following commads should be added 
! to Makefile:
!
!GITVER := $(shell git describe --dirty --always --tags)
!DIRECT  = -DVER=\"$(GITVER)\"
!-----------------------------------------------------------

 use omp_lib

   implicit none

   integer, intent(in) :: IGO

   integer       :: threads
   character( 8) :: CDate, CurTime
   character(11) :: CTime, CurDate
   character*500 :: title, path_exe, path_sim


   title = 'hGAST hydro-servo-aero-elastic tool'

!--- Call the system date function.
   call DATE_AND_TIME ( CDate      )
   call DATE_AND_TIME ( TIME=CTime )

!--- Parse out the day.
   CurDate(1:3) = CDate(7:8)//'-'

!--- Parse out the month.
   select case     ( CDate  (5:6) )
     case ( '01' );  CurDate(4:6) = 'Jan'
     case ( '02' );  CurDate(4:6) = 'Feb'
     case ( '03' );  CurDate(4:6) = 'Mar'
     case ( '04' );  CurDate(4:6) = 'Apr'
     case ( '05' );  CurDate(4:6) = 'May'
     case ( '06' );  CurDate(4:6) = 'Jun'
     case ( '07' );  CurDate(4:6) = 'Jul'
     case ( '08' );  CurDate(4:6) = 'Aug'
     case ( '09' );  CurDate(4:6) = 'Sep'
     case ( '10' );  CurDate(4:6) = 'Oct'
     case ( '11' );  CurDate(4:6) = 'Nov'
     case ( '12' );  CurDate(4:6) = 'Dec'
   end select

!--- Parse out the year.
   CurDate(7:11) = '-'//CDate(1:4)
   CurTime       = CTime(1:2)//':'//CTime(3:4)//':'//CTime(5:6)


   if (IGO==1) then
!---- Get the Git version
#  ifndef VER
    VER = 'unknown'
#  endif

!---- Get the paths
    call getarg (0,path_exe)     !get the exe        path
    call getcwd (  path_sim)     !get the simulation path
   
!---- Output Version, Date, Time, exe path, input path and OMP threads
    write (10,'( a)'  )                       trim (title   )
    write (10,  *     )
    write (10,'(2a)'  ) 'Version         : ', VER
    write (10,'(2a)'  ) 'Date            : ', trim (CurDate )
    write (10,'(2a)'  ) 'Time            : ', trim (CurTime )
    write (10,'(2a)'  ) 'Exe path        : ', trim (path_exe)
    write (10,'(2a)'  ) 'Sim path        : ', trim (path_sim)
#  ifdef HAVE_OMP
!$omp parallel
    threads = omp_get_num_threads()
!$omp end parallel
#  else
    write(10,'(a)') 'OMP not present'
    threads = 1
#  endif
    write (10,'(a,i4)') 'OMP threads     : ', threads

   elseif (IGO==2) then

    write (10,  *     )
    write (10,'( a)'  ) 'End of simulation'
    write (10,  *     )
    write (10,'(2a)'  ) 'Date            : ', trim (CurDate )
    write (10,'(2a)'  ) 'Time            : ', trim (CurTime )
   endif !IGO


 END Subroutine Runinfo
!!----------------------------------------------------------------------
! Subroutine Polynomial ( X, Y, Np, N, Coeff )
!!----------------------------------------------------------------------
!
!   implicit none
!
!   integer, intent(in ) ::   Np , N
!   real(8), intent(in ) :: X(Np), Y(Np)
!   real(8), intent(out) :: Coeff(N)
!
!   real(8) :: A(N)
!   integer :: Ipivot(N),i,j,ii, info(2)
!
!   do i      = 1, N
!      ii     = 
!   do j      = 1, N
!      A(i,j) = X(i)**(j-1)
!   enddo !j
!      B(i  ) = Y(i)
!   enddo !i
!
!   call DGETRF (     N, N, A, N, Ipivot,       info(1)); if (info(1)/=0) then; write(*,*)'error in DGETRF',info(1);stop
!   call DGETRS ('N', N, 1, A, N, Ipivot, B, N, info(2)); if (info(2)/=0) then; write(*,*)'error in DGETRS',info(2);stop
!
! END Subroutine Polynomial
