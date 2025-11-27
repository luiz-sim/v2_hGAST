include 'module_truss.f90'
!----------------------------------------------------------------------
 Module Cbeam
!----------------------------------------------------------------------
   real(8), save :: UflInit     (6    ) !--->floater%
   real(8), save :: Q_truss     (6    ) !--->floater%
!  real(8), save ::   AT_float  (3,3  )
   real(8), save ::    A_float  (3,3  )
   real(8), save ::   DA_float  (3,3  )
   real(8), save ::  DDA_float  (3,3  )
   real(8), save ::    R_float  (3    )
   real(8), save ::   DR_float  (3    )
   real(8), save ::  DDR_float  (3    )
   real(8), save ::   A0_float  (3,3,3)
   real(8), save ::  DA0_float  (3,3,3)
   real(8), save :: DDA0_float  (3,3,3)

 END Module Cbeam
!----------------------------------------------------------------------
 Program moorings
!----------------------------------------------------------------------

 Use Cbeam
 Use Truss

   implicit None


   Real     (8)  :: AQ0(3,3), AQ01(3,3), AQ02(3,3), AQ03(3,3)
   Real     (8)  :: TIME_MAX,TIME,DT_el,rho,depth,GRAV, RLX
   integer       :: NTIME, NTIMEM, IP, i, int_type
   character*256 :: file_truss


!call mkl_set_dynamic (.false.) !to use all threads (slow)
!call mkl_set_dynamic (.true.)
 call omp_set_num_threads(1)
 call mkl_set_num_threads(1)
 call mkl_domain_set_num_threads(1)


!-- open dfile.inp
   open (1,file='dfile.inp')

    read (1,*) TIME_MAX
    read (1,*) DT_el
    read (1,*) GRAV
    read (1,*) depth
    read (1,*) rho
    read (1,*) file_truss
    read (1,*)

    do i=1,6
       read (1,*) UflInit(i)
    enddo

   close(1)

   write (*,*)'TIME_MAX',TIME_MAX
   write (*,*)'DT_el   ',DT_el
   write (*,*)'GRAV    ',GRAV
   write (*,*)'depth   ',depth
   write (*,*)'rho     ',rho

   TIME           =   0.d0
   NTIMEM         = int(TIME_MAX/DT_el) + 1
   UflInit (4:6)  = UflInit (4:6) * dacos(-1.d0)/180.d0


!-- Set init values
      Q_truss = 0.d0;
      A_float = 0.d0;
     DA_float = 0.d0;
    DDA_float = 0.d0;
      R_float = 0.d0;
     DR_float = 0.d0;
    DDR_float = 0.d0;
     A0_float = 0.d0;
    DA0_float = 0.d0;
   DDA0_float = 0.d0;


!-- calc floater A,R
   R_float(1:3) = UflInit (1:3)
   call DIAGO (3,A_float)
   do i=1,3
      call ROT_MATRIX ( i , UflInit(i+3), AQ0, AQ01, AQ02, AQ03 )
      A_float = matmul(A_float,AQ0)
   enddo


!-- Init truss code as usual
   call INIT_tr ( GRAV, rho , depth, DT_el, file_truss )

!---          dynamic solution    static solution
   IP       =       0           !       0
   int_type =       2           !       1
   RLX      =       1.d0        !       0.7d0

   do NTIME = 1, NTIMEM

      TIME = TIME + DT_el 

      write (* ,*)
      write (* ,*) 'NTIME = ', NTIME, TIME

      call MOORINGS_UPDATE_tr

!--- Moorings contribution for the floater
      call calc_connection_UT_tr                     !- calculate UT,UT1,UT2 at connection points
      call MOORINGS_tr (TIME, int_type, IP, RLX)     !- for static solution RLX=0.8, No substeps, and 51 max iterations!
      call truss2float_loads_uncoupl                 !-

      call WRITELEM_tr(NTIME)

      open (1,file='qstuss.dat',access='append')

         write (1,100)   TIME, Q_truss(1:6)

      close(1)

   enddo !NTIME

 100  format (150f21.8)


 END Program moorings
!----------------------------------------------------------------------
!
! ---  Subroutines from hydroGAST ----
!
!----------------------------------------------------------------------
 Subroutine EXTEPR ( A,B,C )
!----------------------------------------------------------------------

   Implicit None

   Real(Kind=8) :: A(3), B(3), C(3)

!...External product of two 3D vectors a x b = c 

   C(1) =  A(2)*B(3) -A(3)*B(2)
   C(2) = -A(1)*B(3) +A(3)*B(1)
   C(3) =  A(1)*B(2) -A(2)*B(1)

 END Subroutine EXTEPR
!----------------------------------------------------------------------
!   Subroutine :DIAGO
!
!      A(1:NM,1:NM)
!----------------------------------------------------------------------
 Subroutine DIAGO (NM,A)
!----------------------------------------------------------------------

   implicit none

   integer, intent(in)  :: NM
   real(8), intent(out) :: A(NM,NM)
!-- local vars
   integer :: i

   A = 0.d0;

   do i=1,NM
      A(i,i) = 1.d0
   enddo


 END Subroutine DIAGO
!----------------------------------------------------------------------
 Subroutine ROT_MATRIX ( IAX, Q0, AQ0, AQ01, AQ02, AQ03 )
!----------------------------------------------------------------------

   Implicit None

   Integer       :: IAX

   Real (kind=8) :: AQ0(3,3), AQ01(3,3), AQ02(3,3), AQ03(3,3), &
                    CQ0     , SQ0      , Q0


   CQ0 = dcos(Q0)
   SQ0 = dsin(Q0)


   if (IAX == 1) then

      AQ0  (1,1) =  1.d0;   AQ01 (1,1) =  0.d0;   AQ02 (1,1) =  0.d0;   AQ03 (1,1) =  0.d0;
      AQ0  (1,2) =  0.d0;   AQ01 (1,2) =  0.d0;   AQ02 (1,2) =  0.d0;   AQ03 (1,2) =  0.d0;
      AQ0  (1,3) =  0.d0;   AQ01 (1,3) =  0.d0;   AQ02 (1,3) =  0.d0;   AQ03 (1,3) =  0.d0;
      AQ0  (2,1) =  0.d0;   AQ01 (2,1) =  0.d0;   AQ02 (2,1) =  0.d0;   AQ03 (2,1) =  0.d0;
      AQ0  (2,2) =   CQ0;   AQ01 (2,2) = - SQ0;   AQ02 (2,2) = - CQ0;   AQ03 (2,2) =   SQ0;
      AQ0  (2,3) = - SQ0;   AQ01 (2,3) = - CQ0;   AQ02 (2,3) =   SQ0;   AQ03 (2,3) =   CQ0;
      AQ0  (3,1) =  0.d0;   AQ01 (3,1) =  0.d0;   AQ02 (3,1) =  0.d0;   AQ03 (3,1) =  0.d0;
      AQ0  (3,2) =   SQ0;   AQ01 (3,2) =   CQ0;   AQ02 (3,2) = - SQ0;   AQ03 (3,2) = - CQ0;
      AQ0  (3,3) =   CQ0;   AQ01 (3,3) = - SQ0;   AQ02 (3,3) = - CQ0;   AQ03 (3,3) =   SQ0;

 
   elseif (IAX == 2) then
 
      AQ0  (1,1) =   CQ0;   AQ01 (1,1) = - SQ0;   AQ02 (1,1) = - CQ0;   AQ03 (1,1) =   SQ0;
      AQ0  (1,2) =  0.d0;   AQ01 (1,2) =  0.d0;   AQ02 (1,2) =  0.d0;   AQ03 (1,2) =  0.d0;
      AQ0  (1,3) =   SQ0;   AQ01 (1,3) =   CQ0;   AQ02 (1,3) = - SQ0;   AQ03 (1,3) = - CQ0;
      AQ0  (2,1) =  0.d0;   AQ01 (2,1) =  0.d0;   AQ02 (2,1) =  0.d0;   AQ03 (2,1) =  0.d0;
      AQ0  (2,2) =  1.d0;   AQ01 (2,2) =  0.d0;   AQ02 (2,2) =  0.d0;   AQ03 (2,2) =  0.d0;
      AQ0  (2,3) =  0.d0;   AQ01 (2,3) =  0.d0;   AQ02 (2,3) =  0.d0;   AQ03 (2,3) =  0.d0;
      AQ0  (3,1) = - SQ0;   AQ01 (3,1) = - CQ0;   AQ02 (3,1) =   SQ0;   AQ03 (3,1) =   CQ0;
      AQ0  (3,2) =  0.d0;   AQ01 (3,2) =  0.d0;   AQ02 (3,2) =  0.d0;   AQ03 (3,2) =  0.d0;
      AQ0  (3,3) =   CQ0;   AQ01 (3,3) = - SQ0;   AQ02 (3,3) = - CQ0;   AQ03 (3,3) =   SQ0;

   elseif (IAX == 3) then
 
      AQ0  (1,1) =   CQ0;   AQ01 (1,1) = - SQ0;   AQ02 (1,1) = - CQ0;   AQ03 (1,1) =   SQ0;
      AQ0  (1,2) = - SQ0;   AQ01 (1,2) = - CQ0;   AQ02 (1,2) =   SQ0;   AQ03 (1,2) =   CQ0;
      AQ0  (1,3) =  0.d0;   AQ01 (1,3) =  0.d0;   AQ02 (1,3) =  0.d0;   AQ03 (1,3) =  0.d0;
      AQ0  (2,1) =   SQ0;   AQ01 (2,1) =   CQ0;   AQ02 (2,1) = - SQ0;   AQ03 (2,1) = - CQ0;
      AQ0  (2,2) =   CQ0;   AQ01 (2,2) = - SQ0;   AQ02 (2,2) = - CQ0;   AQ03 (2,2) =   SQ0;
      AQ0  (2,3) =  0.d0;   AQ01 (2,3) =  0.d0;   AQ02 (2,3) =  0.d0;   AQ03 (2,3) =  0.d0;
      AQ0  (3,1) =  0.d0;   AQ01 (3,1) =  0.d0;   AQ02 (3,1) =  0.d0;   AQ03 (3,1) =  0.d0;
      AQ0  (3,2) =  0.d0;   AQ01 (3,2) =  0.d0;   AQ02 (3,2) =  0.d0;   AQ03 (3,2) =  0.d0;
      AQ0  (3,3) =  1.d0;   AQ01 (3,3) =  0.d0;   AQ02 (3,3) =  0.d0;   AQ03 (3,3) =  0.d0;

   endif !IAX
 

 END Subroutine ROT_MATRIX
!----------------------------------------------------------------------
 Subroutine INT_2_CHAR ( NP, CNUM, N )
!----------------------------------------------------------------------

! Transform the integer N of NP digits max to character CNUM(NP),
! for opening files

   implicit None

   integer      :: NP, N, i, i0, i1

   Character    :: CNUM(NP)

   i0 = N

   do i=1,NP
      i1      = i0/10**(NP-i)
      i0      = i0-i1*10**(NP-i)
      CNUM(i) = char(48+i1)
   enddo


 END Subroutine INT_2_CHAR
!----------------------------------------------------------------------
 Subroutine UINFLOW_sea ( XG, UG, AG, Icalc )
!----------------------------------------------------------------------

!------------- DUMMY Subroutine-------

   Implicit None

   Real(kind=8) :: XG (3), UG(3), AG (3)
   Integer      :: Icalc


   UG = 0.d0;
   AG = 0.d0;

   return

   write(*,*) XG,Icalc !- to avoid warnings


 END Subroutine UINFLOW_sea

include 'truss.f90'
