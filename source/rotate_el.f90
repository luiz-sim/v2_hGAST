!----------------------------------------------------------------------
 Module Rotate 
!----------------------------------------------------------------------
  type :: type_rot_tmp
!------------------------
    real(8) :: R0     (3  )
    real(8) :: ATDR0  (3  )
    real(8) :: ATDR1  (3  )
    real(8) :: ATDDR0 (3  )
    real(8) :: ATDDR1 (3  )
    real(8) :: ATDDR2 (3  )
    real(8) :: A0     (3,3)
    real(8) :: AT0    (3,3)
    real(8) :: ATDA0  (3,3)
    real(8) :: ATDA1  (3,3)
    real(8) :: ATDDA0 (3,3)
    real(8) :: ATDDA1 (3,3)
    real(8) :: ATDDA2 (3,3)
  END type type_rot_tmp
!------------------------
  type :: type_rotder
!------------------------
    real(8)            , allocatable :: der1 (:,:,:)
  END type type_rotder
!------------------------
  type :: type_rotdder
!------------------------
    type(type_rotder)  , allocatable :: der2 (:)
  END type type_rotdder
!------------------------
  type :: type_rotddder
!------------------------
    type(type_rotdder) , allocatable :: der3 (:)
  END type type_rotddder
!
!------------------------
  type :: type_rotation
!------------------------
    real(8)             :: A  (3,3)
    type(type_rotder)   :: A0
    type(type_rotdder)  :: B0
    type(type_rotddder) :: C0
  END type type_rotation
!------------------------

 type(type_rotation), allocatable :: rotation (:) ! IQrot_T
 type(type_rot_tmp ), allocatable :: rot_tmp  (:) ! NQS


 real(8)             , save :: A         (3,3)
 real(8)             , save :: AT        (3,3)
 real(8)             , save :: DA        (3,3)
 real(8)             , save :: ATDA      (3,3)
 real(8)             , save :: DDA       (3,3)
 real(8)             , save :: ATDDA     (3,3)

 real(8)             , save :: ATP       (3,3)
!!real(8)             , save :: DATP      (3,3)
!!real(8)             , save :: DDATP     (3,3)
 real(8)             , save :: ATPA      (3,3)
!!real(8)             , save :: DATPA     (3,3)
!!real(8)             , save :: DDATPA    (3,3)

 real(8)             , save :: ATDR      (  3)
 real(8)             , save :: ATDDR     (  3)
 real(8)             , save :: R         (  3)
 real(8)             , save :: DR        (  3)
 real(8)             , save :: DDR       (  3)
 real(8)             , save :: St_R      (  3)
 real(8)             , save :: St_DR     (  3)
 real(8)             , save :: St_DDR    (  3)

 real(8), allocatable, save :: A0      (:,:,:) !(IQrot_T,3,3)
 real(8), allocatable, save :: AT0     (:,:,:)
 real(8), allocatable, save :: DA0     (:,:,:)
 real(8), allocatable, save :: ATDA0   (:,:,:)
 real(8), allocatable, save :: ATDA1   (:,:,:)
 real(8), allocatable, save :: DDA0    (:,:,:)
 real(8), allocatable, save :: ATDDA0  (:,:,:)
 real(8), allocatable, save :: ATDDA1  (:,:,:)
 real(8), allocatable, save :: ATDDA2  (:,:,:)
 real(8), allocatable, save :: DDAT0A0 (:,:,:)
 real(8), allocatable, save :: DDAT0A1 (:,:,:)
 real(8), allocatable, save :: DDAT0A2 (:,:,:)

 real(8), allocatable, save :: ATDDR0    (:,:) !(IQ_T,3)
 real(8), allocatable, save :: ATDDR1    (:,:)
 real(8), allocatable, save :: ATDDR2    (:,:)
 real(8), allocatable, save :: R0        (:,:)
 real(8), allocatable, save :: DDR0      (:,:)
 real(8), allocatable, save :: DDR1      (:,:)
 real(8), allocatable, save :: DDR2      (:,:)
 real(8), allocatable, save :: St_R0     (:,:)
 real(8), allocatable, save :: St_DDR0   (:,:)
 real(8), allocatable, save :: St_DDR1   (:,:)
 real(8), allocatable, save :: St_DDR2   (:,:)
 real(8), allocatable, save :: ATDR0     (:,:)
 real(8), allocatable, save :: ATDR1     (:,:)
 real(8), allocatable, save :: DR0       (:,:)
 real(8), allocatable, save :: DR1       (:,:)
 real(8), allocatable, save :: St_DR0    (:,:)
 real(8), allocatable, save :: St_DR1    (:,:)

 integer             , save :: NDFROT    , NDFROT0
 integer             , save :: NDFTRA    , NDFTRA0
 integer             , save :: NDFROT0tow, NDFTRA0tow
 integer             , save :: NDFROT0sha, NDFTRA0sha

 real(8), allocatable, save :: QBr         (:) !(IQrot_T)
 real(8), allocatable, save :: derQBr      (:)  
 real(8), allocatable, save :: dderQBr     (:)   
 real(8), allocatable, save :: IAX_QBr     (:) 
 real(8), allocatable, save :: QBt         (:) !(IQtr_T)
 real(8), allocatable, save :: derQBt      (:)  
 real(8), allocatable, save :: dderQBt     (:)   
 real(8), allocatable, save :: IAX_QBt     (:) 
 real(8), allocatable, save :: AQ      (:,:,:) !(IQrot_T,3,3) 
 real(8), allocatable, save :: AQ1     (:,:,:) 
 real(8), allocatable, save :: AQ2     (:,:,:) 
 real(8), allocatable, save :: AQ3     (:,:,:) 
 real(8), allocatable, save :: AP        (:,:) !(IQtr_T,3) 
 real(8), allocatable, save :: AP1       (:,:) 

 integer             , save :: nsub_tow   , nsub_sh   , nsub_bl
 integer             , save :: IQfl_rot_T , IQfl_tr_T 
 integer             , save :: IQtow_rot_T, IQtow_tr_T
 integer             , save :: IQsh_rot_T , IQsh_tr_T 
 integer             , save :: IQbl_rot_T , IQbl_tr_T 
 integer             , save :: IQrot_T    , IQtr_T    , IQ_T

 integer             , save :: NBODBT, NBLADE


 END Module Rotate
!-----------------------------------------------------------
! --- Subroutine : Init_Rot
!-----------------------------------------------------------
 Subroutine Init_Rot ( NBODSUB_in, NBODBT_in, NBODBTWT_in, NBLADE_in, IQfl_rot_in, IQfl_tr_in )

 use Rotate

   implicit none

   integer :: NBODSUB_in (NBODBT_in), NBODBT_in, NBODBTWT_in, NBLADE_in, n
   integer :: IQfl_rot_in, IQfl_tr_in, IQtow_rot, IQtow_tr
   integer :: IQsh_rot   , IQsh_tr   , IQbl_rot , IQbl_tr


   IQfl_rot_T = IQfl_rot_in
   IQfl_tr_T  = IQfl_tr_in 
   NBLADE     = NBLADE_in
   NBODBT     = NBODBTWT_in

   if (NBODBTWT_in == NBLADE_in+2) then !tow,sh,blades
      n        = NBLADE+2
      nsub_tow = NBODSUB_in(n)
      nsub_sh  = NBODSUB_in(n-1)
   elseif (NBODBTWT_in == NBLADE_in+1) then !no tow, only sh and blades
      nsub_tow = -1
      nsub_sh  = NBODSUB_in(NBODBTWT_in)
   elseif (NBODBTWT_in == NBLADE_in) then !only blades
      nsub_tow = -1
      nsub_sh  = -1
   endif
      nsub_bl  = NBODSUB_in(1)

   IQtow_rot   = (nsub_tow+1)*3
   IQtow_tr    = (nsub_tow+1)*3
   IQsh_rot    = (nsub_sh +1)*3
   IQsh_tr     = (nsub_sh +1)*3
   IQbl_rot    = (nsub_bl +1)*4
   IQbl_tr     = (nsub_bl +1)*3

   IQtow_rot_T = IQfl_rot_T  + IQtow_rot + 1
   IQtow_tr_T  = IQfl_tr_T   + IQtow_tr  + 1 !old
   IQsh_rot_T  = IQtow_rot_T + IQsh_rot  + 5
   IQsh_tr_T   = IQtow_tr_T  + IQsh_tr   + 2
   IQrot_T     = IQsh_rot_T  + IQbl_rot  + 4
   IQtr_T      = IQsh_tr_T   + IQbl_tr       !!+1 Hhub
   IQ_T        = IQrot_T     + IQtr_T


   write(10,*)
   write(10,*)'INIT tr-rot'
   write(10,*)'Tower Tot  ', IQtow_rot_T, IQtow_tr_T
   write(10,*)'Shaft Tot  ', IQsh_rot_T , IQsh_tr_T
   write(10,*)'Blade Tot  ', IQrot_T    , IQtr_T    , IQ_T
   write(10,*)'Subbodies  ', nsub_bl    , nsub_sh   , nsub_tow
   write(10,*)


 END Subroutine Init_Rot
!-----------------------------------------------------------
! --- Subroutine : Alloc_Rot
!----------------------------------------------------------
 Subroutine Alloc_Rot ( NQS_in )

 use Rotate

   implicit none

   integer :: m,n,l, NQS_in


   Allocate ( A0      (IQrot_T,3,3),&
              AT0     (IQrot_T,3,3),&
              DA0     (IQrot_T,3,3),&
              ATDA0   (IQrot_T,3,3),&
              ATDA1   (IQrot_T,3,3),&
              DDA0    (IQrot_T,3,3),&
              ATDDA0  (IQrot_T,3,3),&
              ATDDA1  (IQrot_T,3,3),&
              ATDDA2  (IQrot_T,3,3),&
              DDAT0A0 (IQrot_T,3,3),&
              DDAT0A1 (IQrot_T,3,3),&
              DDAT0A2 (IQrot_T,3,3)  )

   Allocate ( ATDDR0  (IQ_T,3),&
              ATDDR1  (IQ_T,3),&
              ATDDR2  (IQ_T,3),&
              R0      (IQ_T,3),&
              DDR0    (IQ_T,3),&
              DDR1    (IQ_T,3),&
              DDR2    (IQ_T,3),&
              St_R0   (IQ_T,3),&
              St_DDR0 (IQ_T,3),&
              St_DDR1 (IQ_T,3),&
              St_DDR2 (IQ_T,3),&
              ATDR0   (IQ_T,3),&
              ATDR1   (IQ_T,3),&
              DR0     (IQ_T,3),&
              DR1     (IQ_T,3),&
              St_DR0  (IQ_T,3),&
              St_DR1  (IQ_T,3)  )

   Allocate ( QBr    (IQrot_T),&
              derQBr (IQrot_T),&
              dderQBr(IQrot_T),&
              IAX_QBr(IQrot_T),&
              QBt    (IQtr_T ),&
              derQBt (IQtr_T ),&
              dderQBt(IQtr_T ),&
              IAX_QBt(IQtr_T )  )

   Allocate ( AQ     (IQrot_T,3,3),&
              AQ1    (IQrot_T,3,3),&
              AQ2    (IQrot_T,3,3),&
              AQ3    (IQrot_T,3,3),&
              AP     (IQtr_T ,3  ),&
              AP1    (IQtr_T ,3  )  )


 Allocate( rot_tmp  (NQS_in ) )
 Allocate( rotation (IQrot_T) )

 do      m = 1, IQrot_T               !n=1:m
         Allocate( rotation(m)%A0%der1(1:m,3,3) )
         Allocate( rotation(m)%B0%der2(1:m)     )
         Allocate( rotation(m)%C0%der3(1:m)     )

   do    n = 1, m                             !l=n:m
         Allocate( rotation(m)%B0%der2(n)%der1(n:m,3,3) )
         Allocate( rotation(m)%C0%der3(n)%der2(n:m)     )

      do l = n, m                                     !k=l:m
         Allocate( rotation(m)%C0%der3(n)%der2(l)%der1(l:m,3,3) )
      enddo !l
   enddo !n
 enddo !m


 END Subroutine Alloc_Rot
!----------------------------------------------------------
! --- Subroutine : Dealloc_Rot
!----------------------------------------------------------
 Subroutine Dealloc_Rot 

 use Rotate

   implicit none


   Deallocate ( A0     , DA0    , DDA0   ,&
                ATDA0  , ATDA1  , AT0    ,&
                ATDDA0 , ATDDA1 , ATDDA2 ,&
                DDAT0A0, DDAT0A1, DDAT0A2   )

   Deallocate ( ATDDR0 , DDR0   , St_DDR0,&
                ATDDR1 , DDR1   , St_DDR1,&
                ATDDR2 , DDR2   , St_DDR2,&
                         R0     , St_R0     )

   Deallocate ( ATDR0, DR0, St_DR0,&
                ATDR1, DR1, St_DR1   )

   Deallocate ( QBr     , QBt    ,&
                derQBr  , derQBt ,&
                dderQBr , dderQBt,&
                IAX_QBr , IAX_QBt,&
                AQ      , AQ1    ,&
                AQ2     , AQ3    ,&
                AP      , AP1       )

   Deallocate ( rotation )
   Deallocate ( rot_tmp  )


 END Subroutine Dealloc_Rot
!-----------------------------------------------------------
! --- Subroutine : ROT_TRAN_BODYB
!     Elementary matrices of rotation and translation
!----------------------------------------------------------
 Subroutine ROT_TRAN_BODYB ( IROTONLY, mrot1, mrot2, mtr1, mtr2 )

 use Rotate

   implicit none

   integer :: IROTONLY, mrot1, mrot2, mtr1, mtr2
   real(8) :: AQ0(3,3), AQ01(3,3), AQ02(3,3), AQ03(3,3)
   real(8) :: AP0(  3), AP01(  3), Q0
   integer :: m, IAX


!---- Elementary Rotation Matricies
   do m   = mrot1, mrot2 !1, NDFROT
      Q0  =     QBr(m)
      IAX = IAX_QBr(m)
 
      if (IAX==0) cycle !the rotation matrix has been already set

      call ROT_MATRIX (IAX, Q0, AQ0, AQ01, AQ02, AQ03)
 
      AQ ( m,1:3,1:3) =  AQ0 (1:3,1:3)
      AQ1( m,1:3,1:3) =  AQ01(1:3,1:3)
      AQ2( m,1:3,1:3) =  AQ02(1:3,1:3)
      AQ3( m,1:3,1:3) =  AQ03(1:3,1:3)
   enddo

   if (IROTONLY == 1) return

!---- Elementary Translation Matricies
   do m   = mtr1, mtr2 !1, NDFTRA
      Q0  =     QBt(m)
      IAX = IAX_QBt(m)
   
      call TRANS_MATRIX (IAX, Q0, AP0, AP01)
 
      AP (m,1:3) = AP0 (1:3)
      AP1(m,1:3) = AP01(1:3)
   enddo


 END Subroutine ROT_TRAN_BODYB
!-----------------------------------------------------------------
!--- Calculate A  [ rotation(NQ)%A (3,3)                        =    A(m)                ]
!--- Calculate A0 [ rotation(NQ)%A0%der1(n,3,3)                 =   dA(m)/dq_n           ]
!--- Calculate B0 [ rotation(NQ)%B0%der2(n)%der1(l,3,3)         =  ddA(m)/dq_n,dq_l      ]
!--- Calculate C0 [ rotation(NQ)%C0%der3(n)%der2(l)%der1(k,3,3) = dddA(m)/dq_n,dq_l,dq_k ]
!-----------------------------------------------------------------
 Subroutine AA0B0C0 ( NQ0, NQ )
!-----------------------------------------------------------------

 use Rotate
!use omp_lib

   implicit none

   integer :: NQ0, NQ, m, n, l, k 
!  real(8) :: t(0:1)


!t(0) = omp_get_wtime()
   if (NQ0 == 1) then
!    m=n=l=k=1
      rotation(1)%A                      (  1:3,1:3) = AQ (1,1:3,1:3)
      rotation(1)%A0                %der1(1,1:3,1:3) = AQ1(1,1:3,1:3)
      rotation(1)%B0        %der2(1)%der1(1,1:3,1:3) = AQ2(1,1:3,1:3)
      rotation(1)%C0%der3(1)%der2(1)%der1(1,1:3,1:3) = AQ3(1,1:3,1:3)
   endif

   do          m = max(2,NQ0), NQ
!$omp parallel shared  (m,rotation,AQ,AQ1,AQ2,AQ3) &
!$omp          private (n,l,k)
!$omp do 
      do       n = 1, m-1
         do    l = n, m-1
!--A
            do k = l, m-1
                  rotation(m)%C0%der3(n)%der2(l)%der1(k,1:3,1:3) = matmul( rotation(m-1)%C0%der3(n)%der2(l)%der1(k,1:3,1:3), AQ (m,1:3,1:3) )
            enddo !k
!--B
               k = m
                  rotation(m)%B0        %der2(n)%der1(l,1:3,1:3) = matmul( rotation(m-1)%B0        %der2(n)%der1(l,1:3,1:3), AQ (m,1:3,1:3) )
                  rotation(m)%C0%der3(n)%der2(l)%der1(k,1:3,1:3) = matmul( rotation(m-1)%B0        %der2(n)%der1(l,1:3,1:3), AQ1(m,1:3,1:3) )
         enddo !l
!--C
               l = m
               k = m
                  rotation(m)%A0                %der1(n,1:3,1:3) = matmul( rotation(m-1)%A0                %der1(n,1:3,1:3), AQ (m,1:3,1:3) )
                  rotation(m)%B0        %der2(n)%der1(l,1:3,1:3) = matmul( rotation(m-1)%A0                %der1(n,1:3,1:3), AQ1(m,1:3,1:3) )
                  rotation(m)%C0%der3(n)%der2(l)%der1(k,1:3,1:3) = matmul( rotation(m-1)%A0                %der1(n,1:3,1:3), AQ2(m,1:3,1:3) )
      enddo !n
!$omp end do 
!$omp end parallel
!--D
               n = m
               l = m
               k = m
                  rotation(m)%A                      (  1:3,1:3) = matmul( rotation(m-1)%A                        (1:3,1:3), AQ (m,1:3,1:3) )
                  rotation(m)%A0                %der1(n,1:3,1:3) = matmul( rotation(m-1)%A                        (1:3,1:3), AQ1(m,1:3,1:3) )
                  rotation(m)%B0        %der2(n)%der1(l,1:3,1:3) = matmul( rotation(m-1)%A                        (1:3,1:3), AQ2(m,1:3,1:3) )
                  rotation(m)%C0%der3(n)%der2(l)%der1(k,1:3,1:3) = matmul( rotation(m-1)%A                        (1:3,1:3), AQ3(m,1:3,1:3) )
   enddo !m
!t(1) = omp_get_wtime()
!write(*,7)'AA0BOCO',t(1)-t(0)
! 7 format (a7,f12.4)


 END Subroutine AA0B0C0
!-----------------------------------------------------------------------
!
!  For a given rotation NQ..
!
!  Set A, AT, A0 and AT0

!  Calculate linear expansion of dA/dt and d2A/dt2
!    DA , DA0 , (DA1 = A0           ) 
!    DDA, DDA0, (DDA1=2*DA0, DDA2=A0)
!
!-----------------------------------------------------------------------
 Subroutine A_DA_DDA ( NQ )
!-----------------------------------------------------------------------

 use Rotate

   implicit none

   integer :: NQ, m, n, l, k
   real(8) :: sumtr(3,3), sumdiag(3,3), DA0tmp(3,3), DDA0tmp(3,3)


   A0  (:,1:3,1:3) = 0.d0;
   AT0 (:,1:3,1:3) = 0.d0;
   DA0 (:,1:3,1:3) = 0.d0;
   DDA0(:,1:3,1:3) = 0.d0;


   A  (     1:3,1:3) = rotation(NQ)%A      (     1:3,1:3)
   A0 (1:NQ,1:3,1:3) = rotation(NQ)%A0%der1(1:NQ,1:3,1:3)
   AT                = transpose(A)
   do l              = 1, NQ
      AT0(l,1:3,1:3) = transpose(A0(l,1:3,1:3))
   enddo


!--- Calculate DA0(m,n,3,3), DDA0(m,n,3,3)

!  DA1(m,n) = A0(m,n)
! DDA1(m,n) = 2*DA0(m,n)
! DDA2(m,n) = A0(m,n)


!--- DDA0 (C0 terms, diagonal and upper triagonal x 2)
            m                = NQ
!!$omp parallel shared  (m,rotation,derQBr,DDA0) &
!!$omp          private (n,l,k,sumtr,sumdiag)
!!$omp do 
   do       n                = 1, m
            sumtr            = 0.d0
            sumdiag          = 0.d0
      do    l                = 1, n-1 !here use previous C0(m,n,l,k)=C0(m,l,k,n),l<n, k<=n
            k                = l                        !C0(m,n,l,k)=C0(m,l,n,k),l<n, k> n !symmetric matrix of previous n
            sumdiag(1:3,1:3) = sumdiag(1:3,1:3)+rotation(m)%C0%der3(l)%der2(k)%der1(n,1:3,1:3)*derQBr(l)*derQBr(k) !(m,l,k,n)
   
         do k                = l+1, n
            sumtr  (1:3,1:3) = sumtr(1:3,1:3)+rotation(m)%C0%der3(l)%der2(k)%der1(n,1:3,1:3)*derQBr(l)*derQBr(k) !(m,l,k,n)
         enddo !k
   
         do k                = n+1, m
            sumtr  (1:3,1:3) = sumtr(1:3,1:3)+rotation(m)%C0%der3(l)%der2(n)%der1(k,1:3,1:3)*derQBr(l)*derQBr(k) !(m,l,n,k)
         enddo !k
      enddo !l
   
      do    l                = n, m
            k                = l
            sumdiag(1:3,1:3) = sumdiag(1:3,1:3)+rotation(m)%C0%der3(n)%der2(l)%der1(k,1:3,1:3)*derQBr(l)*derQBr(k) !(m,n,l,k)
   
         do k                = l+1, m
            sumtr  (1:3,1:3) = sumtr(1:3,1:3)+rotation(m)%C0%der3(n)%der2(l)%der1(k,1:3,1:3)*derQBr(l)*derQBr(k) !(m,n,l,k)
         enddo !k
      enddo !l
            DDA0 (n,1:3,1:3) = sumdiag(1:3,1:3) + 2.d0*sumtr(1:3,1:3)
   enddo !n
!!$omp end do 
!!$omp end parallel


!--- DA0, DDA0(m,n,3,3) (B0 terms)
!!$omp parallel shared  (m,rotation,derQBr,dderQBr,DA0,DDA0) &
!!$omp          private (n,l,DA0tmp,DDA0tmp)
!!$omp do 
   do    n       = 1, m
         DA0tmp  = 0.d0
         DDA0tmp = 0.d0
      do l       = 1, n-1
         DA0tmp (1:3,1:3) = DA0tmp (1:3,1:3) + rotation(m)%B0%der2(l)%der1(n,1:3,1:3)* derQBr(l) !switch l,n symmetric matrix
         DDA0tmp(1:3,1:3) = DDA0tmp(1:3,1:3) + rotation(m)%B0%der2(l)%der1(n,1:3,1:3)*dderQBr(l) !switch l,n symmetric matrix
      enddo !l
      do l       = n ,m
         DA0tmp (1:3,1:3) = DA0tmp (1:3,1:3) + rotation(m)%B0%der2(n)%der1(l,1:3,1:3)* derQBr(l)
         DDA0tmp(1:3,1:3) = DDA0tmp(1:3,1:3) + rotation(m)%B0%der2(n)%der1(l,1:3,1:3)*dderQBr(l)
      enddo !l
         DA0  (n,1:3,1:3) =                     DA0tmp(1:3,1:3)
         DDA0 (n,1:3,1:3) = DDA0 (n,1:3,1:3) + DDA0tmp(1:3,1:3)
   enddo !n
!!$omp end do 
!!$omp end parallel


!--- DA, DDA (m,3,3) both terms
      DA  = 0.d0;
      DDA = 0.d0;

   do n   = 1, NQ
      DA (1:3,1:3) = DA (1:3,1:3) + A0 (n,1:3,1:3)* derQBr(n)
      DDA(1:3,1:3) = DDA(1:3,1:3) + A0 (n,1:3,1:3)*dderQBr(n) + DA0(n,1:3,1:3)* derQBr(n)
   enddo !n


 END Subroutine A_DA_DDA
!-----------------------------------------------------------------------
 Subroutine ATDA_ATDDA ( NQ,nbod ) 
!-----------------------------------------------------------------------

 use Rotate

   implicit none

   integer :: NQ, l, nbod 


!-- ATDA,ATDDA
   ATDA  = matmul (AT,DA )
   ATDDA = matmul (AT,DDA)

!$omp parallel shared  (NQ,ATDDA0,ATDDA1,ATDDA2,ATDA0,ATDA1,AT,AT0,A0,DA,DA0,DDA,DDA0) &
!$omp          private (l)
!$omp do 
   do l = 1, NQ

!-- ATDDA0,1,2
      ATDDA0 (l,1:3,1:3) = matmul( AT (  1:3,1:3), DDA0(l,1:3,1:3) ) + &
                           matmul( AT0(l,1:3,1:3), DDA (  1:3,1:3) )
      ATDDA1 (l,1:3,1:3) = matmul( AT (  1:3,1:3), DA0 (l,1:3,1:3) )                  ! DDA1(m,n) = 2*DA0(m,n), but multiply after
      ATDDA2 (l,1:3,1:3) = matmul( AT (  1:3,1:3), A0  (l,1:3,1:3) )                  ! DDA2(m,n) = A0(m,n)

!-- ATDA0,1
      ATDA0  (l,1:3,1:3) = ATDDA1 (l,1:3,1:3) + matmul( AT0(l,1:3,1:3), DA(1:3,1:3) )
      ATDA1  (l,1:3,1:3) = ATDDA2 (l,1:3,1:3)

      ATDDA1 (l,1:3,1:3) = 2.d0 * ATDDA1 (l,1:3,1:3)                                  !  DA1(m,n) = A0(m,n)
   enddo
!$omp end do 
!$omp end parallel

!----
   if (nbod <= NBLADE) then
      ATPA   = matmul (ATP  , A)
!!    DATPA  = matmul (DATP , A) + matmul (ATP, DA )
!!    DDATPA = matmul (DDATP, A) + matmul (ATP, DDA) + 2.d0* matmul (DATP, DA)
   endif


 END Subroutine ATDA_ATDDA
!-----------------------------------------------------------------------
 Subroutine R_DR_DDR_INIT (Ibod)
!-----------------------------------------------------------------------

 use Rotate

   implicit none

   integer :: NQ, m1, m2, Ibod


   if (Ibod == NBLADE+2) then !tower

      R   = 0d0;  R0   = 0d0
      DR  = 0d0;  DR0  = 0d0; DR1  = 0d0
      DDR = 0d0;  DDR0 = 0d0; DDR1 = 0d0; DDR2 = 0d0


!------ Give sequence
!------ 3 trans, NO rotation
      NQ  = 0
      m1  = 1
      m2  = IQfl_tr_T

      call DIAGO (3,A)
      DA  = 0d0;
      DDA = 0d0;

!------ Calculate R_float
      call R_DR_DDR ( NQ, m1, m2 )

!------ Calculate A_float
      NQ = IQfl_rot_T

      if (NQ>0) call A_DA_DDA ( NQ )

      call ROTMAT_TRANS_floater !transfer R and A float


!------ if the bottom of the tower isn't at global z=0 (ex. monopile, jacket or floating)
!------ pi/2 and Htow0
      NQ = IQfl_rot_T + 1
      m1 = IQfl_tr_T  + 1
      m2 = IQfl_tr_T  + 1

      call A_DA_DDA ( NQ ) !for 1st R
      call R_DR_DDR ( NQ, m1, m2 )

   elseif (Ibod == NBLADE+1) then !Shaft

!------ Hshof and 3 rot (pi/2, YAW fixed, YAW free movement)
      NQ = IQtow_rot_T + 3
      m1 = IQtow_tr_T  + 1
      m2 = IQtow_tr_T  + 1

      call A_DA_DDA ( NQ )
      call R_DR_DDR ( NQ, m1, m2 )

!------ Hsh (shaft) and 5 rot
!------ if nbod=5 or 4 => Hsh=0. Only rotations take place
      NQ = IQtow_rot_T + 5
      m1 = IQtow_tr_T  + 2
      m2 = IQtow_tr_T  + 2

      call A_DA_DDA ( NQ )
      call R_DR_DDR ( NQ, m1, m2 )

   elseif (Ibod == 1) then !Blade1

!------ Store R,DR,DDR,R0,DDR0,DDR1,DDR2
      St_R    = R
      St_DR   = DR
      St_DDR  = DDR

      St_R0   = R0
      St_DR0  = DR0
      St_DR1  = DR1
      St_DDR0 = DDR0
      St_DDR1 = DDR1
      St_DDR2 = DDR2

!-- -pi/2, pi/2, fi and CONE rotations
! done in rotmat  

   else !1<Ibod<=NBLADE
   
!------ Load R,DR,DDR,R0,DDR0,DDR1,DDR2
      R    = St_R
      DR   = St_DR
      DDR  = St_DDR

      R0   = St_R0
      DR0  = St_DR0
      DR1  = St_DR1
      DDR0 = St_DDR0
      DDR1 = St_DDR1
      DDR2 = St_DDR2

!-- -pi/2, pi/2, fi and CONE rotations
! done in rotmat
   endif


 END Subroutine R_DR_DDR_INIT
!------------------------------------------------------------------
!
!  Calculate:   R, R0, DR, DR0, DR1, DDR, DDR1, DDR2, DDR3
!
!                        (1st order expansions) 
!
!           for m1:m2 translations defined in the same c. s.
!       and for  1:NQ sequencial rotations from global c.s. to local 
!                 (in which m1:m2 translations are defined).
!  
!------------------------------------------------------------------
 Subroutine R_DR_DDR ( NQ, m1, m2 )
!------------------------------------------------------------------

 use Rotate

   implicit none

   integer :: NQ, m1, m2, m, mm, n
   real(8) :: A_nq  (3,3), DA_nq(3,3), DDA_nq(3,3)
   real(8) :: R_m1m2(3)  , Rq1  (3)  , Rq2   (3)  
   real(8) :: A0n   (3,3), DA0n (3,3), DDA0n (3,3)
   real(8) :: AR    (3)  , R0m  (3)           
   real(8) :: ADR1  (3)  , ADR2 (3)  , ADR3  (3)
   real(8) :: ADDR1 (3)  , ADDR2(3)  , ADDR3 (3)  
   real(8) :: ADDR4 (3)  , ADDR5(3)  , ADDR6 (3)


!--- Set ct arrays A_nq, DA_nq, DDA_nq, R_m1m2, Rq1, Rq2
     A_nq(1:3,1:3) =   A(1:3,1:3);
    DA_nq(1:3,1:3) =  DA(1:3,1:3);
   DDA_nq(1:3,1:3) = DDA(1:3,1:3);

   R_m1m2 = 0.d0;
   Rq1    = 0.d0;
   Rq2    = 0.d0;

   do m           = m1, m2
      R_m1m2(1:3) = R_m1m2(1:3) + AP (m,1:3)
      Rq1   (1:3) = Rq1   (1:3) + AP1(m,1:3) *  derQBt(m)
      Rq2   (1:3) = Rq2   (1:3) + AP1(m,1:3) * dderQBt(m)
   enddo !m


!--- Calculate R, DR & DDR and add to previous value
     AR     = matmul(  A_nq, R_m1m2)
    ADR1    = matmul( DA_nq, R_m1m2)
    ADR2    = matmul(  A_nq, Rq1   )
   ADDR1    = matmul(DDA_nq, R_m1m2)
   ADDR2    = matmul( DA_nq, Rq1   )
   ADDR3    = matmul(  A_nq, Rq2   )

     R(1:3) =   R(1:3) +   AR (1:3)
    DR(1:3) =  DR(1:3) +  ADR1(1:3) +  ADR2(1:3)
   DDR(1:3) = DDR(1:3) + ADDR1(1:3) + ADDR2(1:3)*2.d0 + ADDR3(1:3)


!--- Calculate R0, DR0, DR1, DDR0, DDR1, DDR2 of translations 
!--- (no need to add to previous)

!$omp parallel shared  (m1,m2,NQ, IQrot_T, AP1, A_nq,DA_nq,DDA_nq, R_m1m2,Rq1,Rq2, A0,DA0,DDA0, R0, DR0,DR1, DDR0,DDR1,DDR2) &
!$omp          private (m,mm,n, R0m, A0n,DA0n,DDA0n, AR, ADR1,ADR2,ADR3, ADDR1,ADDR2,ADDR3,ADDR4,ADDR5,ADDR6)
!$omp do 
   do m            = m1, m2
      mm           = IQrot_T + m

      R0m(1:3)     = AP1(m,1:3)
      AR           = matmul(  A_nq, R0m)
      ADR1         = matmul( DA_nq, R0m)
      ADR2         = matmul(  A_nq, R0m)
      ADDR1        = matmul(DDA_nq, R0m)
      ADDR2        = matmul( DA_nq, R0m)
      ADDR3        = matmul(  A_nq, R0m)

        R0(mm,1:3) =   AR (1:3)
       DR0(mm,1:3) =  ADR1(1:3)
       DR1(mm,1:3) =  ADR2(1:3)
      DDR0(mm,1:3) = ADDR1(1:3)
      DDR1(mm,1:3) = ADDR2(1:3) * 2.d0
      DDR2(mm,1:3) = ADDR3(1:3)
   enddo !m
!$omp end do 


!--- Calculate R0, DR0, DR1, DDR0, DDR1, DDR2 of rotations
!--- and add to previous value 

!$omp do 
   do n              = 1, NQ
      A0n  (1:3,1:3) = A0  (n,1:3,1:3) !  A0_el1 (NQ,n,i,j)
      DA0n (1:3,1:3) = DA0 (n,1:3,1:3) ! DA0_el1 (NQ,n,i,j)
      DDA0n(1:3,1:3) = DDA0(n,1:3,1:3) !DDA0_el1 (NQ,n,i,j)

      AR             = matmul(A0n  , R_m1m2)
      ADR1           = matmul(DA0n , Rq1   )
      ADR2           = matmul(A0n  , Rq2   )
      ADR3           = matmul(A0n  , R_m1m2)
      ADDR1          = matmul(DDA0n, R_m1m2)
      ADDR2          = matmul(DA0n , Rq1   )
      ADDR3          = matmul(A0n  , Rq2   )
      ADDR4          = matmul(DA0n , R_m1m2)
      ADDR5          = matmul(A0n  , Rq1   )
      ADDR6          = matmul(A0n  , R_m1m2)

      R0   (  n,1:3) = R0  (n,1:3) + AR   (1:3)
      DR0  (  n,1:3) = DR0 (n,1:3) + ADR1 (1:3) + ADR2(1:3)
      DR1  (  n,1:3) = DR1 (n,1:3) + ADR3 (1:3)
      DDR0 (  n,1:3) = DDR0(n,1:3) + ADDR1(1:3) + 2.d0 * ADDR2(1:3) + ADDR3(1:3)
      DDR1 (  n,1:3) = DDR1(n,1:3) + 2.d0 * ( ADDR4(1:3) + ADDR5(1:3) )
      DDR2 (  n,1:3) = DDR2(n,1:3) + ADDR6(1:3)
   enddo!n
!$omp end do 
!$omp end parallel


 END Subroutine R_DR_DDR
!-----------------------------------------------------------------------
 Subroutine ATDR_ATDDR
!-----------------------------------------------------------------------

 use Rotate

   implicit none

   integer :: m


!--- ATDR - ATDDR
   ATDR  = matmul (AT,DR )
   ATDDR = matmul (AT,DDR)

!--- ATDR0,1 - ATDDR0,1,2
   do m              = 1, IQrot_T+NDFTRA
      ATDR0  (m,1:3) = matmul ( AT(1:3,1:3), DR0  (m,1:3) )
      ATDR1  (m,1:3) = matmul ( AT(1:3,1:3), DR1  (m,1:3) )
      ATDDR0 (m,1:3) = matmul ( AT(1:3,1:3), DDR0 (m,1:3) )
      ATDDR1 (m,1:3) = matmul ( AT(1:3,1:3), DDR1 (m,1:3) )
      ATDDR2 (m,1:3) = matmul ( AT(1:3,1:3), DDR2 (m,1:3) )
   enddo

   do m              = 1, NDFROT
      ATDR0  (m,1:3) = ATDR0 (m,1:3) + matmul ( AT0(m,1:3,1:3), DR (1:3) )
      ATDDR0 (m,1:3) = ATDDR0(m,1:3) + matmul ( AT0(m,1:3,1:3), DDR(1:3) )
   enddo


 END Subroutine ATDR_ATDDR
!
!
!----------------------------------------------------------------------
 Subroutine ROT_MATRIX ( IAX, Q0, AQ0, AQ01, AQ02, AQ03 )
!----------------------------------------------------------------------

   implicit none

   integer :: IAX, I(0:4), IAXm, IAXp
   real(8) :: AQ0(3,3), AQ01(3,3), AQ02(3,3), AQ03(3,3)
   real(8) :: CQ0     , SQ0      , Q0


   I(0)=3; I(1)=1; I(2)=2; I(3)=3; I(4)=1;

   IAXm = I(IAX-1)
   IAXp = I(IAX+1)

   CQ0  = dcos(Q0)
   SQ0  = dsin(Q0)

   AQ0  =  0.d0;
   AQ01 =  0.d0;
   AQ02 =  0.d0;
   AQ03 =  0.d0;

   AQ0  (IAX , IAX ) =  1.d0;
   AQ0  (IAXp, IAXp) =   CQ0;   AQ01 (IAXp, IAXp) = - SQ0;   AQ02 (IAXp, IAXp) = - CQ0;   AQ03 (IAXp, IAXp) =   SQ0;
   AQ0  (IAXp, IAXm) = - SQ0;   AQ01 (IAXp, IAXm) = - CQ0;   AQ02 (IAXp, IAXm) =   SQ0;   AQ03 (IAXp, IAXm) =   CQ0;
   AQ0  (IAXm, IAXp) =   SQ0;   AQ01 (IAXm, IAXp) =   CQ0;   AQ02 (IAXm, IAXp) = - SQ0;   AQ03 (IAXm, IAXp) = - CQ0;
   AQ0  (IAXm, IAXm) =   CQ0;   AQ01 (IAXm, IAXm) = - SQ0;   AQ02 (IAXm, IAXm) = - CQ0;   AQ03 (IAXm, IAXm) =   SQ0;

!  goto (1,2,3), IAX

!1 AQ0  (1,1) =  1.d0;
!  AQ0  (2,2) =   CQ0;   AQ01 (2,2) = - SQ0;   AQ02 (2,2) = - CQ0;   AQ03 (2,2) =   SQ0;
!  AQ0  (2,3) = - SQ0;   AQ01 (2,3) = - CQ0;   AQ02 (2,3) =   SQ0;   AQ03 (2,3) =   CQ0;
!  AQ0  (3,2) =   SQ0;   AQ01 (3,2) =   CQ0;   AQ02 (3,2) = - SQ0;   AQ03 (3,2) = - CQ0;
!  AQ0  (3,3) =   CQ0;   AQ01 (3,3) = - SQ0;   AQ02 (3,3) = - CQ0;   AQ03 (3,3) =   SQ0;
! return
!
!
!2 AQ0  (1,1) =   CQ0;   AQ01 (1,1) = - SQ0;   AQ02 (1,1) = - CQ0;   AQ03 (1,1) =   SQ0;
!  AQ0  (1,3) =   SQ0;   AQ01 (1,3) =   CQ0;   AQ02 (1,3) = - SQ0;   AQ03 (1,3) = - CQ0;
!  AQ0  (2,2) =  1.d0;
!  AQ0  (3,1) = - SQ0;   AQ01 (3,1) = - CQ0;   AQ02 (3,1) =   SQ0;   AQ03 (3,1) =   CQ0;
!  AQ0  (3,3) =   CQ0;   AQ01 (3,3) = - SQ0;   AQ02 (3,3) = - CQ0;   AQ03 (3,3) =   SQ0;
! return

!
!3 AQ0  (1,1) =   CQ0;   AQ01 (1,1) = - SQ0;   AQ02 (1,1) = - CQ0;   AQ03 (1,1) =   SQ0;
!  AQ0  (1,2) = - SQ0;   AQ01 (1,2) = - CQ0;   AQ02 (1,2) =   SQ0;   AQ03 (1,2) =   CQ0;
!  AQ0  (2,1) =   SQ0;   AQ01 (2,1) =   CQ0;   AQ02 (2,1) = - SQ0;   AQ03 (2,1) = - CQ0;
!  AQ0  (2,2) =   CQ0;   AQ01 (2,2) = - SQ0;   AQ02 (2,2) = - CQ0;   AQ03 (2,2) =   SQ0;
!  AQ0  (3,3) =  1.d0;
! return
 

 END Subroutine ROT_MATRIX
!----------------------------------------------------------------------
 Subroutine ROT_MATRIX0 ( IAX, Q0, AQ0 )
!----------------------------------------------------------------------

   implicit none

   integer :: IAX, I(0:4), IAXm, IAXp
   real(8) :: AQ0(3,3)
   real(8) :: CQ0     , SQ0      , Q0


   I(0)=3; I(1)=1; I(2)=2; I(3)=3; I(4)=1

   IAXm = I(IAX-1)
   IAXp = I(IAX+1)
   CQ0  = dcos(Q0)
   SQ0  = dsin(Q0)
   AQ0  =  0.d0;

   AQ0  (IAX , IAX ) =  1.d0
   AQ0  (IAXp, IAXp) =   CQ0
   AQ0  (IAXp, IAXm) = - SQ0
   AQ0  (IAXm, IAXp) =   SQ0
   AQ0  (IAXm, IAXm) =   CQ0


 END Subroutine ROT_MATRIX0
!-----------------------------------------------------------------
 Subroutine TRANS_MATRIX (IAX, Q0, AP0, AP01)     
!----------------------------------------------------------------------

   implicit none

   integer :: IAX
   real(8) :: AP0(3), AP01(3), Q0


   AP0  = 0.d0;
   AP01 = 0.d0;

   AP0  (IAX) = Q0
   AP01 (IAX) = 1.d0


 END Subroutine TRANS_MATRIX
!
!
!
!-----------------------------------------------------------
! --- Subroutine : ROT_TRAN_BODYB_ja
!     Elementary matrices of rotation and translation 
!     for jacket
!----------------------------------------------------------
 Subroutine ROT_TRAN_BODYB_ja ( AQ0, AP0, mrot, mtr )

 use Rotate

   implicit none

   integer :: mrot, mtr
   real(8) :: AQ0(3,3), AP0(3)


!--- Elementary Rotation Matrices
   AQ ( mrot,1:3,1:3) = AQ0 (1:3,1:3)
   AQ1( mrot,1:3,1:3) = 0.d0;
   AQ2( mrot,1:3,1:3) = 0.d0;
   AQ3( mrot,1:3,1:3) = 0.d0;

!--- Elementary Translation Matrices
   AP ( mtr ,    1:3) = AP0 (1:3)
   AP1( mtr ,    1:3) = 0.d0;


 END Subroutine ROT_TRAN_BODYB_ja
!-----------------------------------------------------------------------
 Subroutine R_DR_DDR_INIT_ja (ISTEP)
!-----------------------------------------------------------------------

 use Rotate

   implicit none

   integer :: NQ, m1, m2, ISTEP


   goto (1,2), ISTEP

!------ 1st time, calculate again the floater q's
 1    R   = 0d0;  R0   = 0d0
      DR  = 0d0;  DR0  = 0d0; DR1  = 0d0
      DDR = 0d0;  DDR0 = 0d0; DDR1 = 0d0; DDR2 = 0d0

!------ Give sequence
!------ 3 trans, NO rotation
      NQ  = 0
      m1  = 1
      m2  = IQfl_tr_T

      call DIAGO (3,A)
      DA  = 0d0;
      DDA = 0d0;

      call R_DR_DDR ( NQ, m1, m2 )

!------ Store R,DR,DDR,R0,DDR0,DDR1,DDR2
      St_R    = R
      St_DR   = DR
      St_DDR  = DDR

      St_R0   = R0
      St_DR0  = DR0
      St_DR1  = DR1
      St_DDR0 = DDR0
      St_DDR1 = DDR1
      St_DDR2 = DDR2

   return

!------ Load R,DR,DDR,R0,DDR0,DDR1,DDR2
 2    R    = St_R
      DR   = St_DR
      DDR  = St_DDR

      R0   = St_R0
      DR0  = St_DR0
      DR1  = St_DR1
      DDR0 = St_DDR0
      DDR1 = St_DDR1
      DDR2 = St_DDR2

      NQ   = IQfl_rot_T
      m1   = IQfl_tr_T  + 1
      m2   = IQfl_tr_T  + 1

      call A_DA_DDA ( NQ         )
      call R_DR_DDR ( NQ, m1, m2 )


 END Subroutine R_DR_DDR_INIT_ja
