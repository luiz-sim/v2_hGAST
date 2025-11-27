!----------------------------------------------------------------------
 Subroutine GELAST_md (NTIME)
!----------------------------------------------------------------------

 Use Cbeam
 Use Craft
 Use Truss
#ifdef HAVE_MPI
 Use MPI
#endif
 Use modal_mod
!Use omp_lib

   implicit none

   integer       :: NTIME, i, it, ICONV, I1, I2, IP
!! character*80  :: outfil
   integer ::  my_rank,ierr
!  real(8) :: t(0:10)
   integer ::  nbod,nbsub,nb_el,nel


#ifdef HAVE_MPI
   call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)
#else
   ierr    = 0
   my_rank = 0
#endif

 if (my_rank==0) then

!  write (6 ,'(a)',advance='no') char(13)
!  write (* ,'("NTIME =",I8,f10.3)',advance='no')  NTIME, TIME

   write (* ,*)
   write (* ,*) 'NTIME = ', NTIME, TIME
   write (10,*)
   write (10,*) 'NTIME = ', NTIME, TIME


!--- Previous step solution
   do i=1,NDFT_el
      UTP_el (i) = UT_el (i)
      UTP1_el(i) = UT1_el(i)
      UTP2_el(i) = UT2_el(i)
   enddo

   if (ITRUSS_el == 1) &
   call MOORINGS_UPDATE_tr


!--- define external local forces/moments
   call define_external_force


   IP    = 1
   ICONV = 0

   do it = 1, ITERMAX

!qq
!     Debug(:,:,:) = 0.d0;
!t(0) = omp_get_wtime()
      call  modal2fem 
      call  Calc_Tow_Shaft_Angles

!------ Solve the dynamic system of control equations
      call control (NTIME,it,1)


!------ Calculate rotation and translation matrices
!t(1) = omp_get_wtime()

      if ((it/=1).or.(ICMOD>0)) &
       call ROTMAT_el            !1st iter: called before writing deformations
                                 !          or in gelast0 or in recall
!t(2) = omp_get_wtime()

      if (subbody(1)%IBTYP_el /= 0) then
!---------- Calculate aerodynamic loads
         if (IAERO_el == 1) then
!qq jim
            if ( dabs(UT1_el(NDFBT_el + NQSW)) < OMEGAR/3.50d0 ) IPARKED = 1
            if ( dabs(UT1_el(NDFBT_el + NQSW)) > OMEGAR/3.00d0 ) IPARKED = 0

            call RAFT     ( NTIME )
         elseif (IAERO_el == 2) then
!qq         call GENUVP_header( NTIME, ICONV, 1 )
            write(*,*) 'Genuvp for modal is disabled atm'
            stop
         endif
      endif
!t(3) = omp_get_wtime()

!------ Hydrodynamic
      call Hydro_main ( it )
!t(4) = omp_get_wtime()

!------ Integrate loads and store
      call LOC_LOADS_int_store_el

!------ Assembly of matrices
      call MATRIX_md
  !!  call CON_MATRIX_el
!t(5) = omp_get_wtime()

!------ Loads communication between subbodies of same body
  !!  call SBOD2SBOD_LOADS                                    !WT's subbodies(same body)
  !!  call BOD2BOD_LOADSJack                                  !jacket's bodies

!------ Loads communication between bodies
  !!  call BOD2BOD_LOADS (1,NBLADE_el)                        !blades to shaft
  !!  call BOD2BOD_LOADS (NBLADE_el+1,NBLADE_el+1)            !shaft  to tower
  !!  call BOD2BOD_LOADS_WT_JA (NBODBT_el)                    !tower  to jacket
!t(6) = omp_get_wtime()

!------ Springs contribution for the foundation modeling
  !!  call Foundation_contrib

!------ Moorings contribution for the floater
      if     (ITRUSS_el == 1) then
         call calc_connection_UT_tr          !- calculate UT,UT1,UT2 at connection points
         call MOORINGS_tr (TIME, 2, IP, 1.d0)
         call truss2float_loads_uncoupl      !-
      elseif (ITRUSS_el == 2) then
         call calc_connection_UT_tr          !- calculate UT,UT1,UT2 at connection points [just for calculation of the position (writing)]
         call MOORINGS_COUPLED_tr ( 2, IP )  !- coupled with truss code
      endif
!t(7) = omp_get_wtime()

!------ Equations for Qs
      call QS_EQUAT_md (IP)
!t(8) = omp_get_wtime()

!------ Reduction of Matrix
      call MATRIX_REDUCT_md !(IP)
!t(9) = omp_get_wtime()

      if ( (NTIME==1).and.(it==1) ) then
         write(*,*) 'NDFBT_el , NDFT_el ',NDFBT_el , NDFT_el , NQS
         write(*,*) 'NDFBT0_el, NDFT0_el',NDFBT0_el, NDFT0_el, NDFT0_el-NDFBT0_el
      endif

!------ Include Modal damping
      call STRUCT_DAMP

!------ Solve the dynamic system of equations
      RLXs = 1.d0
      call Time_Integrate_md (it, 2, ICONV, MAXERR)
!t(10) = omp_get_wtime()
!write(*,7)'rotmat',t( 2)-t(1)
!write(*,7)'aero  ',t( 3)-t(2)
!write(*,7)'hydro ',t( 4)-t(3)
!write(*,7)'matrix',t( 5)-t(4)
!write(*,7)'bloads',t( 6)-t(5)
!write(*,7)'truss ',t( 7)-t(6)
!write(*,7)'qs equ',t( 8)-t(7)
!write(*,7)'reduct',t( 9)-t(8)
!write(*,7)'newmar',t(10)-t(9)
!write(*,7)'total ',t(10)-t(0)
! 7 format (a6,f12.4)
      if (ICONV == 1) goto 1

   enddo !it


   write(* ,*) 'gelast not converged'
   write(10,*) 'gelast not converged'

 1 continue

!! if (TIME>=2.d0) then 
!!    IEIG = 2
!!    call EIGEN_NEW
!!    stop
!! endif

 endif !my_rank

#ifdef HAVE_MPI
   call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif

   ICONV = 1
!qqif (IAERO_el == 2) call GENUVP_header( NTIME, ICONV, 1 )
   if (IAERO_el == 2) write(*,*) 'Genuvp for modal is disabled atm'
   if (IAERO_el == 2) stop

 if (my_rank==0) then

   if (IAERO_el == 1) call WRITEOUTA

!----- change local element force vectors for output
   if (Integr_md==1) then
      do nbod   = 1, NBODBT_el
      do nbsub  = 1, body(nbod)%NBODSUB_el
         nb_el  = body(nbod)%NBODTGLB_el(nbsub)
      do nel    = 1, subbody(nb_el)%NTEB_el
         subbody(nb_el)%AFLOC_el ( nel, 1:6  ) = subbody_md(nb_el)%LoadIntQ_L ( 1:6, nel   )
        if(subbody(nb_el)%NDFPE_el == 12) &
         subbody(nb_el)%AFLOC_el ( nel, 7:12 ) = subbody_md(nb_el)%LoadIntQ_L ( 1:6, nel+1 )
        if(subbody(nb_el)%NDFPE_el == 15) &
         subbody(nb_el)%AFLOC_el ( nel,10:15 ) = subbody_md(nb_el)%LoadIntQ_L ( 1:6, nel+1 )
      enddo !nel
      enddo !nbsub
      enddo !nbod
   endif

   call write_loads
   call modal2fem 
   call ROTMAT_el
   call write_deform ( UT_md)
   call write_qsdof
   call WRITEOUT_zita

   if     (ITRUSS_el == 1) then

      call MOORINGS_UPDATE_tr
      call WRITELEM_tr(NTIME)

   elseif (ITRUSS_el == 2) then

      TIME_tr = TIME
      I1 = ACCU%NDFTACC_el (NBODT_el) + 1
      I2 = NDFBT_el

      UT_tr (1:NDFT_tr) = UT_el (I1:I2);
      UT1_tr(1:NDFT_tr) = UT1_el(I1:I2);
      UT2_tr(1:NDFT_tr) = UT2_el(I1:I2);

      call WRITELEM_tr(NTIME)

   endif
!!!call Writeout_MORISON

!!                            outfil =  'geometry_'
!! call WRITEOUT_GLOB (UT_md, outfil, NTIME, 5)

   if (ICASE_el == 3) &
      call floater_simulator( 3, it )

!--- Update aerodynamic solution parameters
   if (IAERO_el == 1) &
   call UPDATE_AERO (NTIME)

!--- Update oc3 controller variables
   call control (NTIME,it,2)

 endif !my_rank


 END Subroutine GELAST_md
!----------------------------------------------------------------------
 Subroutine GELAST0_md
!----------------------------------------------------------------------

 Use Cbeam

   implicit none


!--- Modal analysis: Eigen value for each physical body with zero
!--- gravity, in order to determine the modal damping matrix
!  call GELAST0_MODAL     !- After the call IMODAL becomes -1

!--- Eigen value for the whole structure
   call GELAST0_EIG_md

!--- set the initial values for basic parameters
   call INIT_DEFLECTIONS ( -OMEGA )

!--- Compute initial deflections from static solution
!  call GELAST_init_md (1)


   call ROTMAT_el


 END Subroutine GELAST0_md
!----------------------------------------------------------------------
 Subroutine GELAST0_EIG_md
!----------------------------------------------------------------------
!
!  Perform eigen-value calculations for the whole structure, then stop.
!  The damping matrix is also taken into account, and the damping ratio
!  is calculated.
!
!  IEIG  : 0  No eigen value
!  IEIG  : 1  Eigen-value with zero gravity
!  IEIG  : 2  Eigen-value including gravity
!
!  Remarks 1. Atm NO aero or wave forces are taken into account.
!          2. The foundation [springs] is taken into account.
!          3. For the floater gravity and buoyancy are always
!             taken into account in order to reach the equilibrium
!             state.
!          4. Floater's 6 dofs appear in the eigen-frequencies.
!             mass     : i.  floater's structural mass
!                        ii. floater's hydrodynamic added mass at
!                            infinity freq.
!             stiffness: i.  moorings linear stiffness matrix
!                        ii. structural  stiffness (positive)
!                        iii.hydrostatic stiffness (negative)
!             No wave forces.
!          5. If ICMOD = 2 the free-free eigen value is calculated,
!             but the static problem isn't possible to be solved.
!             That's why the SOLVE_el isn't called. [The problem is
!             that the free shaft dof has zero stiffness.]
!
! IP=0 : call from GELAST0_EIG  [bc:AK]
! IP=1 : call from GELAST       [bc:AM]
!
!----------------------------------------------------------------------

 Use Cbeam
 Use Paths
 Use modal_mod

   implicit none

   integer      :: it, ICONV, i, IP, j
   character*80 :: outfil


   if ( IEIG == 0) return

   if ((IEIG == 2).and.(ICMOD == 2)) then
      write(*,*) "atm it's not possible to run static simulation with free-free"
      stop
   endif

      write(*,*)
      write(*,*) "EIGEN modal is called"
      write(*,*)

!--- set the initial values for basic parameters
   call INIT_DEFLECTIONS ( -OMEGA )

!--- Set gravity
   if (IEIG == 1) GRAV  = 0.d0
   if (IEIG == 2) GRAV  = GRAVITY_el
       IEIG  = 1

!--- Disable the aero/hydro dynamic loading
!qqsubbody  (1:NBODT_el)%IBTYP_el = 0;
   IAERO_el                       = 0

   IP    = 0
   ICONV = 0

!--- Previous step solution
   do i=1,NDFT_el
      UTP_el (i) = UT_el (i)
      UTP1_el(i) = UT1_el(i)
      UTP2_el(i) = UT2_el(i)
   enddo

!  if (ITRUSS_el == 1) &
!  call MOORINGS_UPDATE_tr


   do it = 1, 25

      if(NBODBTWT_el > 0) UT1_el(NDFBT_el+NQSW) = UTP1_el(NDFBT_el+NQSW)
      call modal2fem 

!------ Calculate rotation and translation matrices
      call ROTMAT_el

!------ Hydrodynamic
!     call Hydro_main ( it )

!------ Integrate loads and store
      call LOC_LOADS_int_store_el

!------ Assembly of matrices
      call MATRIX_md
!     call CON_MATRIX_el

!------ Loads communication between subbodies of same body
!     call SBOD2SBOD_LOADS                                       !WT's subbodies (same body)
!     call BOD2BOD_LOADSJack                                     !jacket's bodies

!------ Loads communication between bodies
!     call BOD2BOD_LOADS(1,NBLADE_el)                            !blades to shaft
!     call BOD2BOD_LOADS(NBLADE_el+1,NBLADE_el+1)                !shaft to tower
!     call BOD2BOD_LOADS_WT_JA (NBODBT_el)                       !tower to jacket

!------ Moorings contribution for the floater
!     call MOORINGS_COUPLED_tr ( 2, IP )                         !coupled with truss code

!------ Springs contribution for the foundation modelling
!     call Foundation_contrib

!------ Equations for Qs
      call QS_EQUAT_md (IP)

!------ Reduction of Matrix
      call MATRIX_REDUCT_md!(IP)
      if ( it == 1 ) then
         write(*,*) 'NDFBT_el , NDFT_el ',NDFBT_el  , NDFT_el , NQS
         write(*,*) 'NDFBT0_el, NDFT0_el',NDFBT0_el, NDFT0_el, NDFT0_el-NDFBT0_el
      endif
  do i=1,NDFBT0_el
  write(101,'(1000e18.8)') (AK_el(i,j),j=1,NDFBT0_el)
  write(102,'(1000e18.8)') (AC_el(i,j),j=1,NDFBT0_el)
  write(103,'(1000e18.8)') (AM_el(i,j),j=1,NDFBT0_el)
  write(104,'(1000e18.8)')  AQ_el(i  )
  enddo
! stop

!------ Include Modal damping contribution
!     call STRUCT_DAMP
!jim
!goto 1
!------ Only for the free-free eigen analysis, unfortunatelly atm no gravity could be included
      if ( ((IEIG == 1).and.(ICMOD == 2)).or.(ICONV == 2).or.(it == 24) ) goto 1

!------ Solve the static system of equations
      RLXs = 1.d0
      call Time_Integrate_md (it, 1, ICONV, 1.d-11)

      if (ICONV == 1) ICONV = 2
   enddo !it

   write(* ,*) 'gelast0_eig not converged'
   write(10,*) 'gelast0_eig not converged'

 1 continue

!--- write the results from the static equilibrium
   call modal2fem 
   call ROTMAT_el
                              outfil = 'geometry_static_'
   call WRITEOUT_GLOB ( UT_md, outfil, 0, 5)
!  call WRITELOADS
   call write_deform  ( UT_md )
   call write_qsdof

!--- Solve the Eigenvalue problem and Write eigenvalues and modeshapes
   call EIGEN

   stop


 END Subroutine GELAST0_EIG_md
!-------------------------------------------------------------------------
!   Subroutine :STIFLOC_tract
! 
!
!-------------------------------------------------------------------------
 Subroutine STIFLOC_tract_md ( AKLOC, AFLOC, nb_el, nel )
!-------------------------------------------------------------------------

 Use Cbeam
 Use modal_mod

   implicit none

   real(8) :: AKLOC   (NDFPEM_el,NDFPEM_el)
   real(8) :: AFLOC   (NDFPEM_el          )
   real(8) :: AKLOC0  (NDFPEM_el,NDFPEM_el)
   real(8) :: AFLOC0  (NDFPEM_el )
   real(8) :: UT      (NDFPEM_el )
   real(8) :: UT1     (NDFPEM_el )
   real(8) :: UT2     (NDFPEM_el )
   real(8) :: UTT0    (NEQPEM_el )
   real(8) :: DUT0    (NEQPEM_el )
   real(8), dimension(NEQPEM_el,NEQPEM_el) :: AK1, AK2, AK3, AK4
   real(8) :: AQ      (NEQPEM_el          )
   real(8) :: SHAPE   (NEQPEM_el,NDFPEM_el)
   real(8) :: DSHAPE  (NEQPEM_el,NDFPEM_el)
   real(8) :: SHAPET  (NDFPEM_el,NEQPEM_el)
   real(8) :: DSHAPET (NDFPEM_el,NEQPEM_el)
   real(8) :: SHAPETAM(NDFPEM_el,NEQPEM_el)
   integer :: nb_el, nel, NEQPE, NDFPE, k, n1, n2
   real(8) :: ALLOC, PA01 , PA02  , allocwgauss
   real(8) :: HTA   , HTAL, HTA1  , HTA0, PHPX  , PHPZ , TRACT0
   real(8) :: Kf(21)
   integer :: IND(6), ii,jj

!qq
!return

   NEQPE   = subbody(nb_el)%NEQPE_el
   NDFPE   = subbody(nb_el)%NDFPE_el
   n1      = subbody(nb_el)%INODNEL_el(nel,1)
   n2      = subbody(nb_el)%INODNEL_el(nel,2)
   HTA0    = subbody(nb_el)%HTA_el    (nel)
   ALLOC   = subbody(nb_el)%ALENG_el  (nel)
   PHPX    = subbody(nb_el)%PHIX_el   (nel)
   PHPZ    = subbody(nb_el)%PHIZ_el   (nel)
   IND(:)  = IMDOF_el(:)


   call LOCAL_UT_1 (nb_el, nel, UT, UT1, UT2)


   do k        = 1, IORDER
      HTAL     = ( XGAUSS(k) + 1.d0 )/2.d0   ! 0 < HTAL   < 1 HTAL=HTA/L
      HTA      = HTAL*ALLOC                  ! 0 < HTA    < L
      HTA1     = HTA0 + HTA
      PA01     = 1.d0 - HTAL
      PA02     = HTAL

      call SHAPEFUNC15 ( HTAL, ALLOC, PHPX, PHPZ, SHAPE, DSHAPE, NDFPE, NEQPE )

      UTT0     = matmul ( SHAPE , UT  )
      DUT0     = matmul ( DSHAPE, UT  ) !for TRACT force

      SHAPET   = transpose (SHAPE )
      DSHAPET  = transpose (DSHAPE)

      if (Integr_md==1) then
         TRACT0   = PA01 * subbody_md(nb_el)%LoadIntQ_L (IMDOF_el(2), nel  ) +&
                    PA02 * subbody_md(nb_el)%LoadIntQ_L (IMDOF_el(2), nel+1)
         AK3                = 0.d0;
         AK3(IND(4),IND(3)) = - TRACT0
         AK3(IND(6),IND(1)) = + TRACT0
         AK4      = 0.d0;
      else
!--------- stiffness
         Kf(1:21) = PA01 * beam_timo(n1)%K(1:21)   +  PA02 * beam_timo(n2)%K(1:21)
!--------- local Int {(N't*K1*N')+(N't*K2*N)+(Nt*K3*N')+(Nt*K4*N)}   [pure structural terms]
!K1   
         AK1(1:6   ,1:6     ) = 0.d0;
         AK1(IND(1),IND(1:6)) = Kf( 1: 6)
         AK1(IND(2),IND(2:6)) = Kf( 7:11)
         AK1(IND(3),IND(3:6)) = Kf(12:15)
         AK1(IND(4),IND(4:6)) = Kf(16:18)
         AK1(IND(5),IND(5:6)) = Kf(19:20)
         AK1(IND(6),IND(  6)) = Kf(   21)
         do ii=2,6
         do jj=1,ii-1
         AK1(IND(ii),IND(jj)) = AK1(IND(jj),IND(ii))
         enddo
         enddo
!K2   
         AK2(1:6   ,1:6   ) = 0.d0;
         AK2(1:6   ,IND(4)) =-AK1(1:6   ,IND(3))
         AK2(1:6   ,IND(6)) = AK1(1:6   ,IND(1))

!----- Evaluate Axial Force in K3 and K4 [Nonlinear term]
         TRACT0             = dot_product(AK1(IND(2),1:6),DUT0(1:6)) &
                             +dot_product(AK2(IND(2),1:6),UTT0(1:6))
!      if (nb_el<= NBODTWT_el.and.IEIG<1) then
         AK3                = 0.d0;
         AK3(IND(4),1:6   ) =-AK1(IND(2),1:6   ) * DUT0(IND(3)) !-K(2,:)wo'
         AK3(IND(6),1:6   ) = AK1(IND(2),1:6   ) * DUT0(IND(1)) ! K(2,:)uo'
         AK3(IND(4),IND(3)) = AK3(IND(4),IND(3)) - TRACT0
         AK3(IND(6),IND(1)) = AK3(IND(6),IND(1)) + TRACT0
         AK4                = 0.d0;
         AK4(IND(4),1:6   ) =-AK2(IND(2),1:6   ) * DUT0(IND(3)) !-K(2,:)wo'
         AK4(IND(6),1:6   ) = AK2(IND(2),1:6   ) * DUT0(IND(1)) ! K(2,:)uo'
      endif                
         SHAPETAM           =          matmul( SHAPET  , AK3   )
         AKLOC0             =        - matmul( SHAPETAM,DSHAPE )
         SHAPETAM           =          matmul( SHAPET  , AK4   )
         AKLOC0             = AKLOC0 - matmul( SHAPETAM, SHAPE )
     
         AQ                 = 0.d0;
         AQ(IND(4))         = TRACT0 * DUT0(IND(3))
         AQ(IND(6))         =-TRACT0 * DUT0(IND(1))
         AFLOC0             =        + matmul( SHAPET, AQ )

!------ Replace
      allocwgauss = 0.5d0*ALLOC*WGAUSS(k)
      AKLOC = AKLOC + allocwgauss*AKLOC0
      AFLOC = AFLOC - allocwgauss*AFLOC0
   enddo !k


 END Subroutine STIFLOC_tract_md
!----------------------------------------------------------------------
!
!  Subrouine : BOD2BOD_LOADS  ----------------
!
!  Loads Communication between WT's bodies (1st Subbody with last one of the previous body)
!
!----------------------------------------------------------------------
 Subroutine  B_LOADS_md ( nbod, nbsub )
!----------------------------------------------------------------------

 Use Cbeam
 Use modal_mod

   implicit none

   integer              :: nbod, nbsub
   real(8)              :: Aba      (3,3        )
   real(8)              :: Aba0     (3,3,NQS    )
   real(8)              :: Rba      (3          )
   real(8)              :: Rba0     (3  ,NQS    )
   real(8)              :: RGba     (3          )
   real(8), allocatable :: AMB2B    (:  ,:  )  !(6,NDFPEM_el)
   real(8), allocatable :: ACB2B    (:  ,:  )  !(6,NDFPEM_el)
   real(8), allocatable :: AKB2B    (:  ,:  )  !(6,NDFPEM_el)
   real(8)              :: AMB2BNQ  (6,NQS      )
   real(8)              :: ACB2BNQ  (6,NQS      )
   real(8)              :: AKB2BNQ  (6,NQS      )
   real(8)              :: AFB2B    (6          )
   integer              :: nn, i, j, iq, nq
   integer              :: nbodB, nbsubB, nnB_el, nbB_el, nelB, nodB, NmodeB, NDFTBB, NDFPEB, nbodB1, nbodB2
   integer              :: nbodA, nbsubA, nnA_el, nbA_el, nelA, nodA, NmodeA, NDFTBA, nnA2_el
   integer              :: IMAX, JMAX, IQRayl, IcalcR
   real(8)              :: COEFM, COEFK
   real(8)              :: AMLOC    (NMODEM_md,NMODEM_md)
   real(8)              :: ACLOC    (NMODEM_md,NMODEM_md)
   real(8)              :: AKLOC    (NMODEM_md,NMODEM_md)
   real(8)              :: AFLOC    (NMODEM_md          )
   real(8)              :: AMLOCNQ  (NMODEM_md,NQS      )
   real(8)              :: ACLOCNQ  (NMODEM_md,NQS      )
   real(8)              :: AKLOCNQ  (NMODEM_md,NQS      )
   real(8), allocatable :: ftt_t   (:,:)
   real(8), allocatable :: ftt     (:,:)
 !!integer :: nbodB, nnsubB, nbsubB, nnB_el,          nbB_el, nelB, nodB
 !!integer :: nbodA, nn    ,         nnA_el, nnA2_el, nbA_el, nelA, nodA

   if (nbsub==body(nbod)%NBODSUB_el) then
!--------- body loads
      if     (nbod<=NBLADE_el  ) then
         return
      elseif (nbod==NBLADE_el+1) then
         nbodB1 = 1
         nbodB2 = NBLADE_el
         nbsubB = 1
      elseif (nbod==NBLADE_el+2) then
         nbodB1 = NBLADE_el+1
         nbodB2 = NBLADE_el+1
         if (Integr_md==1) &            !Loads calculation flag 0:reaction, 1:integration
         nbodB1 = 1
         nbsubB = 1
      endif
         IQRayl = 2   !0: all qs Rayleigh  , 1: body qs with Rayleigh, 2: no Rayleigh to qs
         IcalcR = 1   !0: no R calculations, 1: R calculations
         COEFM  = 0.d0
         COEFK  = 0.d0
   else
!--------- subbody loads
!qq ------ give expresions for subbodies without RGab
         write(*,*)'stopped in B_LOADS_md, not ready for subbodies'
         stop
         return
         nbodB1 = nbod
         nbodB2 = nbod
         nbsubB = nbsub+1

         IQRayl = 1   !0: all qs Rayleigh  , 1: body qs with Rayleigh, 2: no Rayleigh to qs
         IcalcR = 0   !0: no R calculations, 1: R calculations
         COEFM  = COEFM_el(nbod)
         COEFK  = COEFK_el(nbod)
   endif


   do nbodB   = nbodB1, nbodB2
!------ SubBody (B)
      nelB    = 1
      nodB    = 1
      nnB_el  = body      (nbodB )%NNODTGLB_el(nbsubB)
      nbB_el  = body      (nbodB )%NBODTGLB_el(nbsubB)
      NDFPEB  = subbody   (nbB_el)%NDFPE_el
      NmodeB  = subbody_md(nbB_el)%Nmode_TR_md
      NDFTBB  = subbody   (nbB_el)%NDFTB_el
      j       = 0 !ACCU%NDFTACC_el(nbB_el-1)

!------ SubBody (A)
      nbodA   = nbod
      nbsubA  = nbsub
      nnA2_el = body      (nbodA )%NNODTGLB_el(nbsubA+1 ) ! for R during body loads with offset
      nnA_el  = body      (nbodA )%NNODTGLB_el(nbsubA   )
      nbA_el  = body      (nbodA )%NBODTGLB_el(nbsubA   )
      nelA    = subbody   (nbA_el)%NTEB_el
      nodA    = NNPE_el
      NmodeA  = subbody_md(nbA_el)%Nmode_TR_md
      NDFTBA  = subbody   (nbA_el)%NDFTB_el
      nn      = subbody   (nbA_el)%NODMB(nelA,nodA)
      i       = subbody   (nbA_el)%NDFPNBACC_el(nn-1) !ACCU%NDFTACC_el(nbA_el-1)+subbody(nbA_el)%NDFPNBACC_el(nn-1)

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


      if (Integr_md==0) then         !Loads calculation flag 0:reaction, 1:integration
!------ set local loads of SubBody (B)
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

!!       call LOC_LOADS_int_el   ( AMB2B  , ACB2B  , AKB2B  , AFB2B,&
!!                                 AMB2BNQ, ACB2BNQ, AKB2BNQ       ,&
!!                                 nbodB  , nbsubB , JMAX             )

!--------- load end load [loads transfer]
         AMB2B   = 0.d0;   ACB2B   = 0.d0;   AKB2B   = 0.d0;   AFB2B   = 0.d0;
         AMB2BNQ = 0.d0;   ACB2BNQ = 0.d0;   AKB2BNQ = 0.d0;

         AMB2B  ( 1:6, 1:JMAX ) = subbody_md(nbB_el)%LoadIntM   ( 1:6, 1:JMAX )
         ACB2B  ( 1:6, 1:JMAX ) = subbody_md(nbB_el)%LoadIntC   ( 1:6, 1:JMAX )
         AKB2B  ( 1:6, 1:JMAX ) = subbody_md(nbB_el)%LoadIntK   ( 1:6, 1:JMAX )
         AFB2B  ( 1:6         ) = subbody_md(nbB_el)%LoadIntQ   ( 1:6         )
         do iq = 1, QSnode(nnB_el)%Tot_el (2)
            nq =    QSnode(nnB_el)%IORD_el(iq)
         AMB2BNQ( 1:6, nq     ) = subbody_md(nbB_el)%LoadIntMq  ( 1:6, nq     )
         ACB2BNQ( 1:6, nq     ) = subbody_md(nbB_el)%LoadIntCq  ( 1:6, nq     )
         AKB2BNQ( 1:6, nq     ) = subbody_md(nbB_el)%LoadIntKq  ( 1:6, nq     )
         enddo !iq
      endif

!------ transfer from SubBody (B) to SubBody (A) c.s.
      call LOADS_TRANS_el     ( AMB2B  , ACB2B  , AKB2B  , AFB2B,&
                                AMB2BNQ, ACB2BNQ, AKB2BNQ       ,&
                                Aba0   , Aba    , Rba0   , Rba  ,&
                                         nnB_el , IcalcR        ,&
                                JMAX   , JMAX                      )

!------ modal trunctation
      Allocate ( ftt      (NDFTBB, NmodeB) )
      Allocate ( ftt_t    (NmodeA, NDFTBA) )

      ftt  (1:NDFTBB,1:NmodeB) = subbody_md(nbB_el)%FTT_TR_md   (1:NDFTBB,1:NmodeB)
      ftt_t(1:NmodeA,1:NDFTBA) = subbody_md(nbA_el)%FTT_TR_T_md (1:NmodeA,1:NDFTBA)

      AMLOC   = 0.d0; ACLOC   = 0.d0; AKLOC   = 0.d0; AFLOC = 0.d0;
      AMLOCNQ = 0.d0; ACLOCNQ = 0.d0; AKLOCNQ = 0.d0;

      AMLOC  (1:NmodeA,1:NmodeB) = matmul( ftt_t (1:NmodeA,i+1:i+6), matmul( AMB2B  (1:6,1:JMAX), ftt (j+1:j+JMAX,1:NmodeB) ) )
      ACLOC  (1:NmodeA,1:NmodeB) = matmul( ftt_t (1:NmodeA,i+1:i+6), matmul( ACB2B  (1:6,1:JMAX), ftt (j+1:j+JMAX,1:NmodeB) ) )
      AKLOC  (1:NmodeA,1:NmodeB) = matmul( ftt_t (1:NmodeA,i+1:i+6), matmul( AKB2B  (1:6,1:JMAX), ftt (j+1:j+JMAX,1:NmodeB) ) )
      AFLOC  (1:NmodeA         ) = matmul( ftt_t (1:NmodeA,i+1:i+6),         AFB2B  (1:6       ) )
     do iq = 1, QSnode(nnB_el)%Tot_el(2)
        nq = QSnode(nnB_el)%IORD_el(iq)
      AMLOCNQ(1:NmodeA,      nq) = matmul( ftt_t (1:NmodeA,i+1:i+6),         AMB2BNQ(1:6,    nq) )
      ACLOCNQ(1:NmodeA,      nq) = matmul( ftt_t (1:NmodeA,i+1:i+6),         ACB2BNQ(1:6,    nq) )
      AKLOCNQ(1:NmodeA,      nq) = matmul( ftt_t (1:NmodeA,i+1:i+6),         AKB2BNQ(1:6,    nq) )
     enddo!iq

      Deallocate ( ftt  , ftt_t        )
      DeAllocate ( AMB2B, ACB2B, AKB2B )

!------ assembly loads transfer to global matrices
      IMAX  = NmodeA
      JMAX  = NmodeB
      i     = NDFTACC_md(nbA_el-1)
      j     = NDFTACC_md(nbB_el-1)

      call ASSEMBLY_MATRIX_el  ( AMLOC  , ACLOC  , AKLOC    , AFLOC         ,&
                                 AMLOCNQ, ACLOCNQ, AKLOCNQ  , COEFM , COEFK ,&
                                          nnB_el , i        , j     , IQRayl,&  !IQRayl = 0: all qs Rayleigh, 1: body qs with Rayleigh, 2: no Rayleigh to qs
                                 IMAX   , JMAX   , NMODEM_md, NMODEM_md        )
   enddo !nbod


 END Subroutine B_LOADS_md
!--------------------------------------------------------------------------------
!
!
!
!--------------------------------------------------------------------------------
!
!  Subroutine : QS_EQUAT_md   ------------------
!
!  Equations for the q d.o.f
!
!-- Jimmys check the zero of the qs Q=Qo+dq and when ZER0_DOF is necesary.
!----------------------------------------------------------------------
 Subroutine QS_EQUAT_md (IP)
!----------------------------------------------------------------------

 Use Cbeam
 Use Hydro
 Use modal_mod

   implicit none

   integer :: jj, j, ii, i, NQ_el, nbod, nnsub
   integer :: nel, mnod, mn, IQ, ju, ju1, IP  , nbpre_el, nq, nf
   integer :: imod, jm1


   if (NQS  == 0) return !case of jacket alone


!--- NQ = 1:IQfl
!--- Floater's 6 dof equations
   if (ICASE_el < 3) goto 1

!----------------- FLOATER's EQUATIONS ------------------------------

!--- Communicate Forces from Wind Turbine     to 6 Qs
   call FLOATLOADS_WT26Qs_md
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
      if ( (IEIG == 1).and.(IMODAL<=0).and.(ICMOD == 2) )  then
!--------- drive train: free-free
         INDSYSQ(nq) = 0
         call GENLOAD_md
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
   else !ICMOD 1,2,3,[10-->1]
!------ drive train: free-free
      INDSYSQ(nq) = 0

      call GENLOAD_md
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
!--------- Communicate Jacket's qs to WT
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
      do nnsub    = 2, body(nbod)%NBODSUB_el+1

         nbpre_el = body(nbod)%NBODTGLB_el(nnsub-1)
         nel      = subbody(nbpre_el)%NTEB_el
         mnod     = NNPE_el
         mn       =                               subbody(nbpre_el)%NODMB(nel,mnod)
!new        ju    = ACCU%NDFTACC_el(nbpre_el-1) + subbody(nbpre_el)%NDFPNBACC_el(mn-1)
            ju    =                               subbody(nbpre_el)%NDFPNBACC_el(mn-1)

         NQ_el    = NQSP + ACCU%NQSACC_el (nbod-1) + ACCU%NQSBACC_el(nbod, nnsub-1)
 
         do IQ    = 1, 6
            NQ_el = NQ_el + 1
            i     = NDFBT_el + NQ_el
            ju1   = ju + IMDOF_el(IQ)

!------------ new part
            INDSYSQ_bod(NQ_el) = nbpre_el
!------------ end new part
            INDSYSQ    (NQ_el) = ju1

            if (IP == 0) then
               AK_el(i,i  ) = 1.d0
               AQ_el(i    ) = -UT_el (i)
            else
               AM_el(i,i  ) = 1.d0
               AQ_el(i    ) = -UT2_el(i)
            endif


!------------ new part
            do imod = 1, subbody_md(nbpre_el)%Nmode_TR_md
               jm1  =    NDFTACC_md(nbpre_el-1)+imod

               if (IP == 0) then
                  AK_el(i,jm1) =              -             subbody_md(nbpre_el)%FTT_TR_md (ju1,imod)      ! (NBODT_el,max[NDFTB],max[Nmode])
                  AQ_el(i    ) = AQ_el(i    ) + UT_el (jm1)*subbody_md(nbpre_el)%FTT_TR_md (ju1,imod)
               else
                  AM_el(i,jm1) =              -             subbody_md(nbpre_el)%FTT_TR_md (ju1,imod)      ! (NBODT_el,max[NDFTB],max[Nmode])
                  AQ_el(i    ) = AQ_el(i    ) + UT2_el(jm1)*subbody_md(nbpre_el)%FTT_TR_md (ju1,imod)
               endif
            enddo !imod
!------------ end new part


         enddo !IQ
      enddo !nnsub

   enddo !nbod


 END Subroutine QS_EQUAT_md
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
 Subroutine MATRIX_REDUCT_md!(IP)
!----------------------------------------------------------------------

 Use Cbeam
 Use modal_mod

   implicit none

!  integer :: IP
   integer :: INDQ(NQS), nq, i, j, ju, jq, ii, jj, NQS0
   integer :: nb_el, imod
 

   NDFT0_el  = NDFT_el
   NDFBT0_el = NDFBT_el

   do i = 1, NDFBT0_el
      INDSYSB(i)   = i
   enddo

!----------------------------------------------------------
!-- WT sub-bodies: Add dependent Qs to their contribution
!----------------------------------------------------------
   NQS0     = 0
   do nq    = 1, NQS
      ju    = INDSYSQ    (nq)

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
         jq         = NDFBT_el + nq
         nb_el      = INDSYSQ_bod(nq)
      do i          = 1, NDFBT_el + NQSP !!!NDFT_el
      do imod       = 1, subbody_md(nb_el)%Nmode_TR_md
         j          =    NDFTACC_md (nb_el-1) + imod
         AM_el(i,j) = AM_el(i,j) + AM_el(i,jq)*subbody_md(nb_el)%FTT_TR_md(ju,imod)
         AC_el(i,j) = AC_el(i,j) + AC_el(i,jq)*subbody_md(nb_el)%FTT_TR_md(ju,imod)
         AK_el(i,j) = AK_el(i,j) + AK_el(i,jq)*subbody_md(nb_el)%FTT_TR_md(ju,imod)
      enddo !j
      enddo !i
   enddo !nq

!----------------
!- Keep LINES
!----------------

!--- Dependent Q's ---------
   ii    = NDFBT0_el
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

   jj    = NDFBT0_el
!--- Dependent Q's ---------
   do nq = 1, NQS0
      jj = jj + 1
      j  = NDFBT_el   + INDQ(nq)

      AM_el ( 1:NDFT0_el, jj ) = AM_el ( 1:NDFT0_el, j )
      AC_el ( 1:NDFT0_el, jj ) = AC_el ( 1:NDFT0_el, j )
      AK_el ( 1:NDFT0_el, jj ) = AK_el ( 1:NDFT0_el, j )
   enddo !nq


 END Subroutine MATRIX_REDUCT_md
!----------------------------------------------------------------------
!
!     Performes time integration, solving either the
!     static or the dynamic system of equations.
!       itype=1: the STATIC problem is solved K.u = Q
!       itype=2: the NEWMARK method is used for the 
!                dynamic system               M.ddu + C.du + K.u = Q
!
!----------------------------------------------------------------------
 Subroutine Time_Integrate_md (it, itype, ICONV, ERR)
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   integer, allocatable :: Ipivot  (:)
   real(8)              :: c1, c2, c3, c4, c5, UPRE(NDFT0_el), UTPRE(NDFT0_el), RESDT,RESDTM, UTR1,UTR2, ERR
   integer              :: it, itype, ICONV, i, is, info(2), N, imax


!--- Form the system's matrices, depending on the solution type.
   if     (itype == 2) then
!------------------- Newmark-b method
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
   elseif (itype == 1 ) then
!------------------- Static solution
                   c3 = 0.d0
                   c5 = 0.d0
                   N  = NDFT0_el
      UPRE  (1:N    ) = 0.d0;
      UTPRE (1:N    ) = 0.d0;
   endif


!--- call routines to solve the linear system using LU decomposition
   Allocate ( Ipivot (NDFT0_el)   )

   call DGETRF (     NDFT0_el, NDFT0_el, AK_el, NDFT_el, Ipivot,                 info(1))
   call DGETRS ('N', NDFT0_el,   1     , AK_el, NDFT_el, Ipivot, AQ_el, NDFT_el, info(2))

   if (info(1)>0) write(*,*) 'LU Factorization    ',info(1)
   if (info(2)>0) write(*,*) 'LU Back Substitution',info(2)

   Deallocate ( Ipivot )


!--- calculate the maximum and the mean errors, perform convergence check
   RESDT  = dabs(AQ_el(1))
   RESDTM = dabs(AQ_el(1))
   imax   = 1

   do i = 2, NDFT0_el
      RESDT  = RESDT + dabs(AQ_el(i))
      if (RESDTM<dabs(AQ_el(i))) then
         imax   = i
         RESDTM = dabs(AQ_el(i))
      endif
   enddo

   if (((RESDTM >= 0.d0).or.(RESDTM <= 0.d0)).eqv.(.FALSE.)) then
      write (10,*) '** NaN Error **'
      stop
   endif

   RESDT = RESDT/dble(NDFT0_el)

   if (RESDTM < ERR) ICONV = 1

   write (* ,1) it, RESDT, RESDTM, imax, UT_el(INDSYSB(imax)) + AQ_el(imax)
   write (10,1) it, RESDT, RESDTM, imax, UT_el(INDSYSB(imax)) + AQ_el(imax)

1  format (2x,'it = ', i2,4x,'  RESDT =',2(e23.16,2x),i5,  4x,e23.16)


!--- Substitute the solution for all active dofs
   do i  = 1, NDFT0_el
      is = INDSYSB (i)
   
      UT0_el (is) = AQ_el(i)
      UT_el  (is) = UT_el (is) + AQ_el(i)*RLXs

      UTR2        = c3 * ( UT_el(is)-UPRE(i) )
      UT02_el(is) = UTR2     - UT2_el(is)
      UT2_el (is) = UTR2

      UTR1        = UTPRE(i) + c5 * UTR2
      UT01_el(is) = UTR1     - UT1_el(is)
      UT1_el (is) = UTR1    
   enddo

   call GET_REDUSED_DOFS_md


 END Subroutine Time_Integrate_md
!----------------------------------------------------------------------
 Subroutine GET_REDUSED_DOFS_md
!----------------------------------------------------------------------

 Use Cbeam
 Use modal_mod

   implicit none

   integer :: iq, i, iu, nb_el, imod, icoef


   if (IREDUCT_el == 0) return ! No reduction

!--- WT: Dependent Qs
   do iq    = 1, NQS
      if (INDSYSQ (iq) < 1) cycle

      i     = NDFBT_el + iq
      iu    = INDSYSQ    (iq)
      nb_el = INDSYSQ_bod(iq)

      UT0_el(i)=0.d0; UT01_el(i)=0.d0; UT02_el(i)=0.d0;
      UT_el (i)=0.d0; UT1_el (i)=0.d0; UT2_el (i)=0.d0;

      do imod       = 1, subbody_md(nb_el)%Nmode_TR_md
         icoef      =    NDFTACC_md (nb_el-1) + imod

         UT0_el (i) = UT0_el (i) + UT0_el (icoef)*subbody_md(nb_el)%FTT_TR_md(iu,imod)
         UT01_el(i) = UT01_el(i) + UT01_el(icoef)*subbody_md(nb_el)%FTT_TR_md(iu,imod)
         UT02_el(i) = UT02_el(i) + UT02_el(icoef)*subbody_md(nb_el)%FTT_TR_md(iu,imod)
                      
         UT_el  (i) = UT_el  (i) + UT_el  (icoef)*subbody_md(nb_el)%FTT_TR_md(iu,imod)
         UT1_el (i) = UT1_el (i) + UT1_el (icoef)*subbody_md(nb_el)%FTT_TR_md(iu,imod)
         UT2_el (i) = UT2_el (i) + UT2_el (icoef)*subbody_md(nb_el)%FTT_TR_md(iu,imod)
      enddo !imod
   enddo


 END Subroutine GET_REDUSED_DOFS_md
!--------------------------------------------------------------------------------
!
!  Subrouine : GENLOAD_md  ----------------
!
!  Equation for rotational speed
!
!----------------------------------------------------------------------
 Subroutine GENLOAD_md
!----------------------------------------------------------------------

 Use Cbeam
 Use modal_mod

   implicit none

   real(8), allocatable :: AMB2B    (:,:) !(1,NDFPEM_el)
   real(8), allocatable :: ACB2B    (:,:) !(1,NDFPEM_el)
   real(8), allocatable :: AKB2B    (:,:) !(1,NDFPEM_el)
   real(8)              :: AMB2BNQ  (1,NQS      )
   real(8)              :: ACB2BNQ  (1,NQS      )
   real(8)              :: AKB2BNQ  (1,NQS      )
   real(8)              :: AFB2B    (1          )
   real(8)              :: AMB2B_t  (1,NMODEM_md)
   real(8)              :: ACB2B_t  (1,NMODEM_md)
   real(8)              :: AKB2B_t  (1,NMODEM_md)
   integer              :: i, j, NDFPE, iq, nq
   integer              :: nbod, nbsub, nb_el, nn_el, nel
   real(8), allocatable :: ftt     (:,:)
   integer              :: Nmode, NDFTB, JMAX


   if (NBODBTWT_el == NBLADE_el) return
   if (Integr_md==1) then

      call GENLOAD_md_int
      return
   endif


   nbod  = NBLADE_el + 1
   nbsub = 1
   nel   = 1
   nb_el = body(nbod)%NBODTGLB_el(nbsub)
   nn_el = body(nbod)%NNODTGLB_el(nbsub)
   NDFPE = subbody(nb_el)%NDFPE_el
   i     = IMDOF_el(5)
   j     = 0 !ACCU%NDFTACC_el   ( nb_el-1 )

      
!--- Loads of Wind Turbine
!--- set local loads of SubBody (B)
   JMAX  = subbody(nb_el)%NDFPE_el
      Allocate ( AMB2B (1,JMAX) )
      Allocate ( ACB2B (1,JMAX) )
      Allocate ( AKB2B (1,JMAX) )

      AMB2B  ( 1, 1:JMAX ) = subbody(nb_el)%AMLOC_el    ( nel, i, 1:JMAX )
      ACB2B  ( 1, 1:JMAX ) = subbody(nb_el)%ACLOC_el    ( nel, i, 1:JMAX )
      AKB2B  ( 1, 1:JMAX ) = subbody(nb_el)%AKLOC_el    ( nel, i, 1:JMAX )
      AFB2B  ( 1         ) = subbody(nb_el)%AFLOC_el    ( nel, i         ) &
                            + TGenLSS_el + TMLossLSS_el + TBrakeLSS_el ! QTC1_el(6)*RAT_GEAR
      do iq = 1, QSnode(nn_el)%Tot_el(2)
         nq =    QSnode(nn_el)%IORD_el(iq)

      AMB2BNQ( 1, nq     ) = subbody(nb_el)%AMLOCNQ_el  ( nel, i, iq     )
      ACB2BNQ( 1, nq     ) = subbody(nb_el)%ACLOCNQ_el  ( nel, i, iq     )
      AKB2BNQ( 1, nq     ) = subbody(nb_el)%AKLOCNQ_el  ( nel, i, iq     )
      enddo!iq


!------------------- Damping term for stall-regulated WT
                   nq      = NQSW
      ACB2BNQ ( 1, nq    ) = ACB2BNQ ( 1, nq ) + TGenLSS_D_el


!--- modal trunctation
   Nmode  = subbody_md(nb_el)%Nmode_TR_md
   NDFTB  = subbody   (nb_el)%NDFTB_el
   Allocate ( ftt     (NDFTB, Nmode) )

   ftt  (1:NDFTB,1:Nmode) = subbody_md(nb_el)%FTT_TR_md    (1:NDFTB,1:Nmode)
   AMB2B_t (1,1:Nmode) = matmul( AMB2B  (1,1:JMAX), ftt (j+1:j+JMAX,1:Nmode) )
   ACB2B_t (1,1:Nmode) = matmul( ACB2B  (1,1:JMAX), ftt (j+1:j+JMAX,1:Nmode) )
   AKB2B_t (1,1:Nmode) = matmul( AKB2B  (1,1:JMAX), ftt (j+1:j+JMAX,1:Nmode) )

   Deallocate ( ftt )
   DeAllocate ( AMB2B, ACB2B, AKB2B )


!--- LOADS Communicated to Generator (A)
   i    = NDFBT_el + NQSW !-1+1
   j    = NDFTACC_md(nb_el-1)
   JMAX = Nmode

   AM_el( i, j+1:j+JMAX ) = AM_el( i, j+1:j+JMAX ) + AMB2B_t ( 1, 1:JMAX )
   AC_el( i, j+1:j+JMAX ) = AC_el( i, j+1:j+JMAX ) + ACB2B_t ( 1, 1:JMAX )
   AK_el( i, j+1:j+JMAX ) = AK_el( i, j+1:j+JMAX ) + AKB2B_t ( 1, 1:JMAX )
   AQ_el( i             ) = AQ_el( i             ) + AFB2B   ( 1         )

!--- Qs without Rayleigh Damping
   do iq = 1, QSnode(nn_el)%Tot_el(2)
      nq = QSnode(nn_el)%IORD_el(iq)
      j  = NDFBT_el + nq

      AM_el( i, j ) = AM_el( i, j ) + AMB2BNQ ( 1, nq )
      AC_el( i, j ) = AC_el( i, j ) + ACB2BNQ ( 1, nq )
      AK_el( i, j ) = AK_el( i, j ) + AKB2BNQ ( 1, nq )
   enddo !iq


 END Subroutine GENLOAD_md
!----------------------------------------------------------------------
 Subroutine  GENLOAD_md_int
!----------------------------------------------------------------------

 Use Cbeam
 Use modal_mod

   implicit none

   integer              :: nbod, nbsub
   real(8)              :: Aba      (3,3        )
   real(8)              :: Aba0     (3,3,NQS    )
   real(8)              :: Rba      (3          )
   real(8)              :: Rba0     (3  ,NQS    )
   real(8)              :: RGba     (3          )
   real(8), allocatable :: AMB2B    (:  ,:  )  !(6,NDFPEM_el)
   real(8), allocatable :: ACB2B    (:  ,:  )  !(6,NDFPEM_el)
   real(8), allocatable :: AKB2B    (:  ,:  )  !(6,NDFPEM_el)
   real(8)              :: AMB2BNQ  (6,NQS      )
   real(8)              :: ACB2BNQ  (6,NQS      )
   real(8)              :: AKB2BNQ  (6,NQS      )
   real(8)              :: AFB2B    (6          )
   real(8)              :: AMB2B_t  (6,NMODEM_md)
   real(8)              :: ACB2B_t  (6,NMODEM_md)
   real(8)              :: AKB2B_t  (6,NMODEM_md)
   integer              :: i, j, iq, nq
   integer              :: nbodB, nbsubB, nnB_el, nbB_el, nelB, nodB, NmodeB, NDFTBB, NDFPEB, nbodB1, nbodB2
   integer              :: nbodA, nbsubA, nnA_el, nbA_el, nelA, nodA, NmodeA, NDFTBA, nnA2_el
   integer              :: JMAX, IQRayl, IcalcR
   real(8)              :: COEFM, COEFK
   real(8), allocatable :: ftt     (:,:)
 !!integer :: nbodB, nnsubB, nbsubB, nnB_el,          nbB_el, nelB, nodB
 !!integer :: nbodA, nn    ,         nnA_el, nnA2_el, nbA_el, nelA, nodA


       if (Integr_md/=1) return

       nbod   = NBLADE_el+1
       nbsub  = 1

       nbodB1 =           1
       nbodB2 = NBLADE_el+1
       nbsubB = 1
       IQRayl = 2   !0: all qs Rayleigh  , 1: body qs with Rayleigh, 2: no Rayleigh to qs
       IcalcR = 1   !0: no R calculations, 1: R calculations
       COEFM  = 0.d0
       COEFK  = 0.d0

   do nbodB   = nbodB1, nbodB2
!------ SubBody (B)
      nelB    = 1
      nodB    = 1
      nnB_el  = body      (nbodB )%NNODTGLB_el(nbsubB)
      nbB_el  = body      (nbodB )%NBODTGLB_el(nbsubB)
      NDFPEB  = subbody   (nbB_el)%NDFPE_el
      NmodeB  = subbody_md(nbB_el)%Nmode_TR_md
      NDFTBB  = subbody   (nbB_el)%NDFTB_el
      j       = 0 !ACCU%NDFTACC_el(nbB_el-1)

!------ SubBody (A)
      nbodA   = nbod
      nbsubA  = nbsub
!     nnA2_el = body      (nbodA )%NNODTGLB_el(nbsubA+1 )
      nnA_el  = body      (nbodA )%NNODTGLB_el(nbsubA   )
      nbA_el  = body      (nbodA )%NBODTGLB_el(nbsubA   )
      nnA2_el = body      (nbodA )%NNODTGLB_el(nbsubA   )
      nelA    = subbody   (nbA_el)%NTEB_el
      nodA    = 1
      NmodeA  = subbody_md(nbA_el)%Nmode_TR_md
      NDFTBA  = subbody   (nbA_el)%NDFTB_el
!     nn      = subbody   (nbA_el)%NODMB(nelA,nodA)
!     i       = subbody   (nbA_el)%NDFPNBACC_el(nn-1) !ACCU%NDFTACC_el(nbA_el-1)+subbody(nbA_el)%NDFPNBACC_el(nn-1)

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


!--------- Loads from Integration
      JMAX    = subbody   (nbB_el)%NDFTB_el
         Allocate ( AMB2B (6,JMAX) )
         Allocate ( ACB2B (6,JMAX) )
         Allocate ( AKB2B (6,JMAX) )

         AMB2B   = 0.d0;   ACB2B   = 0.d0;   AKB2B   = 0.d0;   AFB2B   = 0.d0;
         AMB2BNQ = 0.d0;   ACB2BNQ = 0.d0;   AKB2BNQ = 0.d0;

         AMB2B  ( 1:6, 1:JMAX ) = subbody_md(nbB_el)%LoadIntM   ( 1:6, 1:JMAX )
         ACB2B  ( 1:6, 1:JMAX ) = subbody_md(nbB_el)%LoadIntC   ( 1:6, 1:JMAX )
         AKB2B  ( 1:6, 1:JMAX ) = subbody_md(nbB_el)%LoadIntK   ( 1:6, 1:JMAX )
         AFB2B  ( 1:6         ) = subbody_md(nbB_el)%LoadIntQ   ( 1:6         )
         do iq = 1, QSnode(nnB_el)%Tot_el (2)
            nq =    QSnode(nnB_el)%IORD_el(iq)
         AMB2BNQ( 1:6, nq     ) = subbody_md(nbB_el)%LoadIntMq  ( 1:6, nq     )
         ACB2BNQ( 1:6, nq     ) = subbody_md(nbB_el)%LoadIntCq  ( 1:6, nq     )
         AKB2BNQ( 1:6, nq     ) = subbody_md(nbB_el)%LoadIntKq  ( 1:6, nq     )
         enddo !iq


!------ transfer from SubBody (B) to SubBody (A) c.s.
      call LOADS_TRANS_el     ( AMB2B  , ACB2B  , AKB2B  , AFB2B,&
                                AMB2BNQ, ACB2BNQ, AKB2BNQ       ,&
                                Aba0   , Aba    , Rba0   , Rba  ,&
                                         nnB_el , IcalcR        ,&
                                JMAX   , JMAX                      )


!------ modal trunctation
      Allocate ( ftt      (NDFTBB, NmodeB) )

      ftt  (1:NDFTBB,1:NmodeB) = subbody_md(nbB_el)%FTT_TR_md   (1:NDFTBB,1:NmodeB)

      AMB2B_t(1:6,1:NmodeB)     = matmul( AMB2B  (1:6,1:JMAX), ftt (j+1:j+JMAX,1:NmodeB) )
      ACB2B_t(1:6,1:NmodeB)     = matmul( ACB2B  (1:6,1:JMAX), ftt (j+1:j+JMAX,1:NmodeB) )
      AKB2B_t(1:6,1:NmodeB)     = matmul( AKB2B  (1:6,1:JMAX), ftt (j+1:j+JMAX,1:NmodeB) )

      Deallocate ( ftt )
      DeAllocate ( AMB2B, ACB2B, AKB2B )


!------ Generator, Mechanical losses and Brake
      if (nbodB==NBLADE_el+1) then
         AFB2B   ( 6          ) = AFB2B   ( 6     ) + TGenLSS_el + TMLossLSS_el + TBrakeLSS_el         !QTC1_el(6)*RAT_GEAR
!------------------- Damping term for stall-regulated WT
                      nq        = NQSW
         ACB2BNQ ( 6, nq      ) = ACB2BNQ ( 6, nq ) + TGenLSS_D_el
      endif


!------ LOADS Communicated to Generator (A) assembly loads
      i    = NDFBT_el + NQSW !-1+1
      j    = NDFTACC_md(nbB_el-1)
      JMAX = NmodeB
      
      AM_el( i, j+1:j+JMAX ) = AM_el( i, j+1:j+JMAX ) + AMB2B_t ( 6, 1:JMAX )
      AC_el( i, j+1:j+JMAX ) = AC_el( i, j+1:j+JMAX ) + ACB2B_t ( 6, 1:JMAX )
      AK_el( i, j+1:j+JMAX ) = AK_el( i, j+1:j+JMAX ) + AKB2B_t ( 6, 1:JMAX )
      AQ_el( i             ) = AQ_el( i             ) + AFB2B   ( 6         )
      
!------ Qs without Rayleigh Damping
      do iq = 1, QSnode(nnB_el)%Tot_el(2)
         nq = QSnode(nnB_el)%IORD_el(iq)
         j  = NDFBT_el + nq
      
         AM_el( i, j ) = AM_el( i, j ) + AMB2BNQ ( 6, nq )
         AC_el( i, j ) = AC_el( i, j ) + ACB2BNQ ( 6, nq )
         AK_el( i, j ) = AK_el( i, j ) + AKB2BNQ ( 6, nq )
      enddo !iq

      enddo !nbod


 END Subroutine GENLOAD_md_int
!--------------------------------------------------------------------------------
!
!  Subrouine : FLOATLOADS_WT26QS_md ----------------
!
!  Loads Communication between WT and 6 qs [ICASE_el == 3 or 4]
!
!----------------------------------------------------------------------
 Subroutine FLOATLOADS_WT26Qs_md
!----------------------------------------------------------------------

 Use Cbeam
 Use modal_mod

   implicit none

   real(8)              :: Aba      (3,3        )
   real(8)              :: Aba0     (3,3,NQS    )
   real(8)              :: Rba      (3          )
   real(8)              :: Rba0     (3  ,NQS    )
   real(8)              :: Ra       (3          )
   real(8), allocatable :: AMB2B    (:  ,:  )  !(6,NDFPEM_el)
   real(8), allocatable :: ACB2B    (:  ,:  )  !(6,NDFPEM_el)
   real(8), allocatable :: AKB2B    (:  ,:  )  !(6,NDFPEM_el)
   real(8)              :: AMB2BNQ  (6,NQS      )
   real(8)              :: ACB2BNQ  (6,NQS      )
   real(8)              :: AKB2BNQ  (6,NQS      )
   real(8)              :: AFB2B    (6          )
   real(8)              :: AMB2B_t  (6,NMODEM_md)
   real(8)              :: ACB2B_t  (6,NMODEM_md)
   real(8)              :: AKB2B_t  (6,NMODEM_md)
   integer              :: nbod1, nbod2, i, j, NDFPE, iq, nq, iqi
   integer              :: nbod, nnsub, nbsub, nn_el, nb_el, nelB, nodB
   real(8), allocatable :: ftt     (:,:)
   integer              :: Nmode, NDFTB, JMAX


   if ((ICASE_el /= 3).and.(ICASE_el /= 4)) return


   if (NBODBTWT_el == NBLADE_el) then
    nbod1 = 1
    nbod2 = NBLADE_el
   else
    nbod1 = NBODBTWT_el!shaft or tower
    nbod2 = NBODBTWT_el
   endif
   if ( Integr_md == 1 ) then    !Loads calculation flag 0:reaction, 1:integration
    nbod1 = 1
    nbod2 = NBODBTWT_el
   endif


   do nbod  = nbod1, nbod2
!------ Body (B) WT or Jacket
      nnsub = 1
      nbsub = 1

      nn_el = body(nbod)%NNODTGLB_el(nnsub)
      nb_el = body(nbod)%NBODTGLB_el(nbsub)
      nelB  = 1
      nodB  = 1

      NDFPE = subbody(nb_el)%NDFPE_el
      j     = 0 !ACCU%NDFTACC_el   ( nb_el-1 )


!------ Body (A) Floater: Global c.o. System
!------ Aba, Aba0, Rba, Rba0
      Ra  (1:3    ) = UT_el(NDFBT_el+1 : NDFBT_el+3)
      Rba (1:3    ) = transf_mat(nn_el)%R_el (1:3    ) - Ra(1:3)  !transf_mat(nn_el)%A_el (1:3,1:3) x (0,0,Htow0)
      Aba (1:3,1:3) = transf_mat(nn_el)%A_el (1:3,1:3)

      do iq = 1, QSnode(nn_el)%Tot_el(2)
         nq = QSnode(nn_el)%IORD_el(iq)

         Aba0 (1:3,1:3,nq) = transf_mat(nn_el)%A0_el(iq ,1:3,1:3) !(nbsub,nq,i,j)
         Rba0 (1:3,    nq) = 0.d0
      enddo!iq 
         Rba0 (1:3,   4:6) = transf_mat(nn_el)%R0_el(4:6,1:3    ) !(nbsub,nq,i,j)   !! A0 . (0,0,Htow)??


      if (Integr_md==0) then         !Loads calculation flag 0:reaction, 1:integration
!------ set local loads of SubBody (B)
      JMAX    = subbody(nb_el)%NDFPE_el
         Allocate ( AMB2B (6,JMAX) )
         Allocate ( ACB2B (6,JMAX) )
         Allocate ( AKB2B (6,JMAX) )

         call LOC_LOADS_el       ( AMB2B  , ACB2B  , AKB2B  , AFB2B,&
                                   AMB2BNQ, ACB2BNQ, AKB2BNQ       ,&
                                   nb_el  , nn_el  ,nelB   , nodB    )
      else
!--------- Loads from Integration

      JMAX    = subbody   (nb_el)%NDFTB_el
         Allocate ( AMB2B (6,JMAX) )
         Allocate ( ACB2B (6,JMAX) )
         Allocate ( AKB2B (6,JMAX) )

!!       call LOC_LOADS_int_el   ( AMB2B  , ACB2B  , AKB2B  , AFB2B,&
!!                                 AMB2BNQ, ACB2BNQ, AKB2BNQ       ,&
!!                                 nbod   , nbsub  , JMAX             )

!--------- load end load [loads transfer]
         AMB2B   = 0.d0;   ACB2B   = 0.d0;   AKB2B   = 0.d0;   AFB2B   = 0.d0;
         AMB2BNQ = 0.d0;   ACB2BNQ = 0.d0;   AKB2BNQ = 0.d0;

         AMB2B  ( 1:6, 1:JMAX ) = subbody_md(nb_el)%LoadIntM   ( 1:6, 1:JMAX )
         ACB2B  ( 1:6, 1:JMAX ) = subbody_md(nb_el)%LoadIntC   ( 1:6, 1:JMAX )
         AKB2B  ( 1:6, 1:JMAX ) = subbody_md(nb_el)%LoadIntK   ( 1:6, 1:JMAX )
         AFB2B  ( 1:6         ) = subbody_md(nb_el)%LoadIntQ   ( 1:6         )
         do iq = 1, QSnode(nn_el)%Tot_el (2)
            nq =    QSnode(nn_el)%IORD_el(iq)
         AMB2BNQ( 1:6, nq     ) = subbody_md(nb_el)%LoadIntMq  ( 1:6, nq     )
         ACB2BNQ( 1:6, nq     ) = subbody_md(nb_el)%LoadIntCq  ( 1:6, nq     )
         AKB2BNQ( 1:6, nq     ) = subbody_md(nb_el)%LoadIntKq  ( 1:6, nq     )
         enddo !iq
      endif

!------ transfer from SubBody (B) to SubBody (A) c.s.
      call LOADS_TRANS_el     ( AMB2B  , ACB2B  , AKB2B  , AFB2B,&
                                AMB2BNQ, ACB2BNQ, AKB2BNQ       ,&
                                Aba0   , Aba    , Rba0   , Rba  ,&
                                         nn_el  , 1             ,& !IcalcR
                                JMAX   , JMAX                      )


!--- modal trunctation
   Nmode  = subbody_md(nb_el)%Nmode_TR_md
   NDFTB  = subbody   (nb_el)%NDFTB_el
   Allocate ( ftt     (NDFTB, Nmode) )

   ftt  (1:NDFTB,1:Nmode) = subbody_md(nb_el)%FTT_TR_md   (1:NDFTB,1:Nmode)
   AMB2B_t (1:6,1:Nmode) = matmul( AMB2B  (1:6,1:JMAX), ftt (j+1:j+JMAX,1:Nmode) )
   ACB2B_t (1:6,1:Nmode) = matmul( ACB2B  (1:6,1:JMAX), ftt (j+1:j+JMAX,1:Nmode) )
   AKB2B_t (1:6,1:Nmode) = matmul( AKB2B  (1:6,1:JMAX), ftt (j+1:j+JMAX,1:Nmode) )

   Deallocate ( ftt )


!------ Assembly to global Matrices
      j     = NDFTACC_md(nb_el-1)  !ACCU%NDFTACC_el(nb_el-1)
      JMAX  = Nmode

      do iq = 1, IQfl
         i  = NDFBT_el + iq 

         AM_el ( i, j+1:j+JMAX ) = AM_el ( i, j+1:j+JMAX ) + AMB2B_t ( IMDOF_el(iq), 1:JMAX )
         AC_el ( i, j+1:j+JMAX ) = AC_el ( i, j+1:j+JMAX ) + ACB2B_t ( IMDOF_el(iq), 1:JMAX )
         AK_el ( i, j+1:j+JMAX ) = AK_el ( i, j+1:j+JMAX ) + AKB2B_t ( IMDOF_el(iq), 1:JMAX )
         AQ_el ( i             ) = AQ_el ( i             ) + AFB2B   ( IMDOF_el(iq)         )
      enddo

!------ For all q's
      do iqi = 1, IQfl
         i   = NDFBT_el + iqi

         do iq = 1, QSnode(nn_el)%Tot_el(2)
            nq = QSnode(nn_el)%IORD_el(iq)
            j  = NDFBT_el + nq

            AM_el ( i, j ) = AM_el ( i, j ) + AMB2BNQ ( IMDOF_el(iqi), nq )
            AC_el ( i, j ) = AC_el ( i, j ) + ACB2BNQ ( IMDOF_el(iqi), nq )
            AK_el ( i, j ) = AK_el ( i, j ) + AKB2BNQ ( IMDOF_el(iqi), nq )
         enddo!iqj 
      enddo!iqi

      DeAllocate ( AMB2B, ACB2B, AKB2B )
   enddo!nbod


 END Subroutine FLOATLOADS_WT26Qs_md
!-------------------------------------------------------------------------------------------------
! -- Subroutine :MODAL_ANAL_md
!------------------------------------------------------------------------------------------------- 
 Subroutine MODAL_ANAL_md
!------------------------------------------------------------------------------------------------- 

 Use Cbeam
 Use modal_mod

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
   integer              :: NDF1,NDF2,nbod,NDFTB,i,j, Nmode
   integer              :: m, ierr, nb1_el, nb2_el , mmax, I1 !, Inorm_mb, n
   real(8)              :: A1, C1, RR, RI !!, RFI21, RFI22, RFI2, RI1, RI2, RR1, RR2
   real(8)              :: CCRIT
   integer              :: ii, nb_el, nel, nn

!-- DGGEV
   real(8)              :: VL(1)
   real(8), allocatable :: WORK  (:)
   integer              :: LWORK
   real(8), allocatable :: ftt_t    (:,:)
   real(8), allocatable :: ftt      (:,:)
   integer, allocatable :: mode     (:  )

 !!real(8), allocatable :: product_f (:,:)
 !!real(8), allocatable :: Md2       (:,:)
   
   write (* ,*)
   write (* ,*) 'MODAL ANALYSIS is called'
   write (* ,*)
   write (10,*)
   write (10,*) 'MODAL ANALYSIS is called'
   write (10,*)

   CDAMPSTR = 0.d0;

   open (100, file='modal_anal_eigen.dat' )

   !open(201,file='Trancated Modes.dat')

!--- for all physical bodies
   do nbod = 1, min(NBODBT_el,NBODBTWT_el+1)

!------ consider the reduction of the zero dofs
!------ still foundation or free bc is not taken into account
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


! pan AK , AM time0 fixxed
      AK(1:NDFTB,1:NDFTB) = AK_el (NDF1+1:NDF2,NDF1+1:NDF2);
      AM(1:NDFTB,1:NDFTB) = AM_el (NDF1+1:NDF2,NDF1+1:NDF2);


!------ Solve the eigenvalue problem
!     call REIGG(NDFTB, NDFTB   ,AK,AM,ALFR,ALFI,BETA,1,V,ierr)
      call DGGEV('N', 'V', NDFTB, AK , NDFTB, AM, NDFTB, ALFR, ALFI, BETA, VL, 1, V, NDFTB, WORK, LWORK, ierr )

      if (ierr>0) write(*,*) 'Modal Anal - DGGEV', ierr, WORK(1)/LWORK

      AK(1:NDFTB,1:NDFTB) = AK_el (NDF1+1:NDF2,NDF1+1:NDF2);
      AM(1:NDFTB,1:NDFTB) = AM_el (NDF1+1:NDF2,NDF1+1:NDF2);
! pan AK , AM time0 fixxed


!------ Check the eigen-values. Set the frequency and the mode matrix.
!------ lamda = (AKFR/BETA) + (ALFI/BETA) yiot = A1 + C1 yiot
      do I = 1, NDFTB

!--------- check BETA
         if (dabs(BETA(I)) < 1.d-10) then
            write(10,*) 'beta = 0',I
            write( *,*) 'beta = 0',I

            RFR (I) = 1.d+10 !0.d0
            RFI (I) = 0.d0

            do J=1,NDFTB
               RZR (J,I) = 0.d0 ! V(J,I)
            enddo
!              RZR (I,I) = 1.d0
            cycle
         endif

!--------- check Imaginary Part and Real Part sign
         A1    = ALFR(I)/BETA(I)
         C1    = ALFI(I)/BETA(I)

         if (dabs(C1) < 1.d-13) then

!------------ Real Eigen-Value -----------
            if     (A1 <  0.d0) then

!--------------- Real Part is Negative --
               write(10,*) 'lamda REAL, A1 neg',I, A1 !, C1
               write( *,*) 'lamda REAL, A1 neg',I, A1 !, C1

               RFR (I) = 1.d+10 !0.d0
               RFI (I) = 0.d0

               do J=1,NDFTB
                  RZR (J,I) = 0.d0 ! V(J,I)
               enddo
!              RZR (I,I) = 1.d0

               cycle

            elseif (A1 >= 0.d0) then

!--------------- Real Part is Positive [Normal case] --
!w             write( *,*) 'lamda REAL, A1 pos',I, A1, dsqrt(A1)/PI2 !, C1

               RR = dsqrt(A1)
               RI = 0.d0
            endif !Real pos or Real neg

         else

!------------ Complex Eigen-Value -----------
            write(10,*) 'lamda COMPLEX',I,C1,A1
            write( *,*) 'lamda COMPLEX',I,C1,A1

            STOP

            RFR (I) = 1.d+10 !0.d0
            RFI (I) = 0.d0

            do J=1,NDFTB
               RZR (J,I) = 0.d0
            enddo
!              RZR (I,I) = 1.d0

            cycle

      !!    RFI21 = 0.5d0*(-A1 + dsqrt(A1**2+C1**2))
      !!    RFI22 = 0.5d0*(-A1 - dsqrt(A1**2+C1**2))
      !!    RFI2  = dmax1(RFI21,RFI22)
      !!    RI1   = dsqrt(RFI2)
      !!    RI2   =-dsqrt(RFI2)

      !!    RR1 = C1/(2.d0*RI1)
      !!    RR2 = C1/(2.d0*RI2)

!     !!   write(*,*)  RFI21, RFI22, RFI2,RI1,RI2
!     !!   write(*,*)  RR1  , RR2

      !!    if (RR1.ge.0.d0) then
      !!        RR = RR1
      !!        RI = RI1
      !!    else
      !!        RR = RR2
      !!        RI = RI2
      !!    endif

         endif !Real or Complex

         RFR (I) = RR / PI2    ! real part of eigenvalue
         RFI (I) = RI / PI2

         do J=1,NDFTB
            RZR (J,I) = V(J,I) ! Fi [ndfbxndfb] or [ndfbxndftr] when reduced
         enddo
      enddo !I


!------ Sort modes with respect to the frequency
      call SORTRX (RFR ,NDFTB,NDFTB,IND)

      FREQ_el (:  ) = 0.d0;
      Md      (:,:) = 0.d0;
      Cd      (:  ) = 0.d0;
      FTT_el  (:,:) = 0.d0;


!------ Normalize Modes
      do I = 1,NDFTB
         I1   = IND(I)
         mmax = 1
         do m=2,NDFTB
            if ( dabs(RZR(m,I1)) > dabs(RZR(mmax,I1)) ) mmax=m
         enddo

         FREQ_el    (I) = RFR(I1)

         if (RZR(mmax,I1) == 0.d0) cycle

         FTT_el (1:NDFTB,I) = RZR(1:NDFTB,I1) / RZR(mmax,I1);
      enddo !I


!------ Determine modal mass
      Nmode                  = subbody_md(nbod)%Nmode_TR_md

      Allocate ( mode (Nmode) )

       mode(1:Nmode        ) = subbody_md(nbod)%mode (1:Nmode)
      Md   (1:Nmode,1:Nmode) = matmul ( transpose(FTT_el(1:NDFTB,mode(1:Nmode))), matmul(AM, FTT_el(1:NDFTB,mode(1:Nmode))) ) !keep the selected modes [mode(1:Nmode)]

!------ Determine modal damping
      do i = 1, Nmode
         ii=     mode(i)
                                       CCRIT=CCRIT_el(nbod)
         if ( FREQ_el(ii) > FREQ0_el ) CCRIT=CCRIT_el(0)
         
         Cd (i) = 2.d0*CCRIT*PI2*FREQ_el(ii)*Md(i,i)
      enddo


!------ Assembly to global Structural Damping Matrix
      if ( body(nbod)%NBODSUB_el>1) then
         write(*,*) 'modal damping in modal approach with sub-bodies is not supported atm'
         stop
      endif
         nb_el = body(nbod)%NBODTGLB_el(1)
         ii    = NDFTACC_md(nb_el-1)
      do i     = 1, Nmode
         CDAMPSTR (ii+i, ii+i) = Cd(i)
      enddo


!------ Write frequency, Modal mass and Modal damping
      do i = 1, Nmode
         ii=     mode(i)
         write (100,101) FREQ_el(ii), Md(i,i), Cd(i)
      enddo

      write (100,*)
      write (100,*)


!------ store reduced modal matrix and it's transpose
      Nmode = subbody_md(nbod)%Nmode_TR_md
      if (Nmode>NDFTB) then
         write(*,*)"error number of modes can't be bigger than body dofs",nbod,Nmode,NDFTB
         stop
      endif

      Allocate ( ftt   (NDFTB+6,Nmode  ) )
      Allocate ( ftt_t (Nmode  ,NDFTB+6) )

      ftt  (1:      6,1:Nmode  ) = 0.d0;                   ! add the 6 reduced dofs due to zero b.c. 
      ftt  (7:NDFTB+6,1:Nmode  ) = FTT_el(1:NDFTB, mode(1:Nmode))


      if (subbody_md(nbod)%itor_coupl == 0 ) then          ! switch for torsion   coupling  [0:disable, 1:enable]
!!    if (Nmode<=4 .and. nbod /= NBLADE_el+1) then
!------------ remove torsion coupling
            nb_el = body(nbod)%NBODTGLB_el(1)
         do nel   = 1, subbody(nb_el)%NTEB_el
            nn    = subbody(nb_el)%NODMB       (nel,1)!nnod)
            i     = subbody(nb_el)%NDFPNBACC_el(nn-1 )
           if     (subbody(nb_el)%NDFPE_el == 12) then
            ii    = i + 6;  ftt  (ii,1:Nmode) = 0.d0;
            ii    = i +12;  ftt  (ii,1:Nmode) = 0.d0;
           elseif (subbody(nb_el)%NDFPE_el == 15) then
            ii    = i + 6;  ftt  (ii,1:Nmode) = 0.d0;
            ii    = i + 8;  ftt  (ii,1:Nmode) = 0.d0;
            ii    = i +15;  ftt  (ii,1:Nmode) = 0.d0;
           endif
         enddo
      endif
      if (subbody_md(nb_el)%iext_coupl == 0) then          ! switch for extension coupling  [0:disable, 1:enable]
!---------  remove extension coupling
         do nel   = 1, subbody(nb_el)%NTEB_el
            nn    = subbody(nb_el)%NODMB       (nel,1)!nnod)
            i     = subbody(nb_el)%NDFPNBACC_el(nn-1 )
           if     (subbody(nb_el)%NDFPE_el == 12) then
            ii    = i + 3;  ftt  (ii,1:Nmode) = 0.d0;
            ii    = i + 9;  ftt  (ii,1:Nmode) = 0.d0;
           elseif (subbody(nb_el)%NDFPE_el == 15) then
            ii    = i + 3;  ftt  (ii,1:Nmode) = 0.d0;
            ii    = i + 7;  ftt  (ii,1:Nmode) = 0.d0;
            ii    = i + 9;  ftt  (ii,1:Nmode) = 0.d0;
            ii    = i +12;  ftt  (ii,1:Nmode) = 0.d0;
           endif
         enddo
      endif


!------
      ftt_t(1:Nmode  ,1:NDFTB+6) = transpose (ftt (1:NDFTB+6,1:Nmode))

      subbody_md(nbod)%FTT_TR_md   (1:NDFTB+6,1:Nmode  ) = ftt  (1:NDFTB+6,1:Nmode  )
      subbody_md(nbod)%FTT_TR_T_md (1:Nmode  ,1:NDFTB+6) = ftt_t(1:Nmode  ,1:NDFTB+6)

      Deallocate ( ftt, ftt_t, mode )


      Deallocate ( AK , AM    , ALFI   , ALFR, BETA, V )
      Deallocate ( RFR, RFI   , RZR    , IND           )
      Deallocate ( Md , FTT_el, FTTI_el                )
      
      Deallocate ( Cd , CCRIT_act, Cd_FTTI_el, AC      )
      Deallocate ( FREQ_el )
      Deallocate ( WORK    )
!! enddo!nbsub
   enddo!nbod

   !close (201)

   close (100)


 101  format ('freq=',F17.3,5x,'modal mass=',F15.3,'modal damp=',F17.3 )


 END Subroutine MODAL_ANAL_md
