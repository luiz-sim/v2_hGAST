#include "module_truss.f90"
#include "truss.f90"
#include "truss_coupled.f90"
!!#include "flap.f90"
!----------------------------------------------------------------------
 Subroutine GELAST (NTIME)
!----------------------------------------------------------------------

 Use Cbeam
 Use Craft
 Use Truss
#ifdef HAVE_MPI
 Use mpi
#endif
#ifdef HAVE_OMP
 Use omp_lib
#endif

   implicit none

   integer      :: NTIME, it, ICONV, I1, I2, IP, itype
   integer      ::  my_rank,ierr
# ifdef HAVE_OMP
   real(8)      :: t(0:10)
# endif


# ifdef HAVE_MPI
   call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)
# else
   ierr    = 0
   my_rank = 0
# endif

 if (my_rank==0) then

!  write (6 ,'(a)',advance='no') char(13)
!  write (* ,'("NTIME =",I8,f10.3)',advance='no')  NTIME, TIME

   write (* ,*)
   write (* ,*) 'NTIME = ', NTIME, TIME
   write (10,*)
   write (10,*) 'NTIME = ', NTIME, TIME


!--- Previous step solution
   UTP_el (1:NDFT_el) = UT_el (1:NDFT_el);
   UTP1_el(1:NDFT_el) = UT1_el(1:NDFT_el);
   UTP2_el(1:NDFT_el) = UT2_el(1:NDFT_el);

   if (ITRUSS_el == 1) &
   call MOORINGS_UPDATE_tr


!--- define external local forces/moments
   call define_external_force
 endif!my_rank

      IP    = 1
      ICONV = 0
   do it    = 1, ITERMAX
     if (my_rank==0) then

#                                                ifdef HAVE_OMP
                                                 t(0) = omp_get_wtime()
#                                                endif
      call Calc_Tow_Shaft_Angles

!------ Solve the dynamic system of control equations
      call control (NTIME,it,1)

#                                                ifdef HAVE_OMP
                                                 t(1) = omp_get_wtime()
#                                                endif
!------ Calculate rotation and translation matrices
      if ((it/=1).or.(ICMOD>0)) &
       call ROTMAT_el            !1st iter: called before writing deformations
                                 !          or in gelast0 or in recall
#                                                ifdef HAVE_OMP
                                                 t(2) = omp_get_wtime()
#                                                endif

      itype = subbody(1)%IBTYP_el
     endif!my rank
#    ifdef HAVE_MPI
      call MPI_BCAST(itype,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#    endif
!!    if (itype /= 0) then 
!---------- Calculate aerodynamic loads
         if (IAERO_el == 1) then
           if(my_rank.eq.0) then
            if ( dabs(UT1_el(NDFBT_el + NQSW)) < OMEGAR/3.50d0 ) IPARKED = 1
            if ( dabs(UT1_el(NDFBT_el + NQSW)) > OMEGAR/3.00d0 ) IPARKED = 0
            call RAFT            (NTIME    )
           endif!my_rank
         elseif (IAERO_el == 2) then
#         ifdef HAVE_GENUVP
           !call GENUVP_main  ( INIT_gn,TIME_gn,DT_gn,Iter_el,IFull_gn,NTIME_gn,NTIMEM)
            call GENUVP_main  ( 0      ,TIME   ,DT_el, 1     ,ICONV   ,NTIME   ,9999,9999 )
#         ifdef HAVE_MPI
           if(my_rank==0)&
#         endif
            call Aero2Elast
#         endif
         else if (IAERO_el == 3) then
#         ifdef HAVE_cfd
           !call CFD_main (WhatToDo, !=2 perform a time step and provide loads 
           !               NTIME   , !may start from 0
           !               TIME,DT , !time and dt are defined in hGAST
           !               NTIMEM)   !max number of steps in case some allocation needs it
            call broadcast_t
            call kinematics(NTIME)
            call CFD_main (2,NTIME,TIME,DT_el,NTIMEM,-1,it) 
            call Aero2Elast_cfd
#         endif
         endif !IAERO_el
!!    endif !itype

     if(my_rank==0) then
#                                                ifdef HAVE_OMP
                                                 t(3) = omp_get_wtime()
#                                                endif
!------ Hydrodynamic
      call Hydro_main ( it )
#                                                ifdef HAVE_OMP
                                                 t(4) = omp_get_wtime()
#                                                endif

!------ Assembly of matrices
      call MATRIX_el
      call CON_MATRIX_el
#                                                ifdef HAVE_OMP
                                                 t(5) = omp_get_wtime()
#                                                endif

!------ Loads communication between subbodies of same body
      call SBOD2SBOD_LOADS                                    !WT's subbodies(same body)
      call BOD2BOD_LOADSJack                                  !jacket's bodies

!------ Loads communication between bodies
      call BOD2BOD_LOADS     (1,NBLADE_el)                    !blades to shaft
      call BOD2BOD_LOADS     (NBLADE_el+1,NBLADE_el+1)        !shaft  to tower
!!!   call BOD2BOD_LOADS_int (          1,NBLADE_el+1)        !blades+shaft  to tower
      call BOD2BOD_LOADS_WT_JA (NBODBT_el)                    !tower  to jacket
      call BOD2BOD_LOADS_VAWT                                 !blade end to tower
#                                                ifdef HAVE_OMP
                                                 t(6) = omp_get_wtime()
#                                                endif

!------ Springs contribution for the foundation modeling
      call Foundation_contrib

!------ Moorings contribution for the floater
      if     (ITRUSS_el == 1) then
         call calc_connection_UT_tr          !- calculate UT,UT1,UT2 at connection points
         call MOORINGS_tr (TIME, 2, IP, 1.d0)
         call truss2float_loads_uncoupl      !-
      elseif (ITRUSS_el == 2) then
         call calc_connection_UT_tr          !- calculate UT,UT1,UT2 at connection points [just for calculation of the position (writing)]
         call MOORINGS_COUPLED_tr ( 2, IP )  !- coupled with truss code
      endif
#                                                ifdef HAVE_OMP
                                                 t(7) = omp_get_wtime()
#                                                endif

!------ Equations for Qs
      call QS_EQUAT (IP)
#                                                ifdef HAVE_OMP
                                                 t(8) = omp_get_wtime()
#                                                endif

!------ Reduction of Matrix
      call MATRIX_REDUCT (IP)
#                                                ifdef HAVE_OMP
                                                 t(9) = omp_get_wtime()
#                                                endif

      if ( (NTIME==1).and.(it==1) ) then
         write(*,'(a,4i6)') 'NDFBT_el , NDFT_el , NQS ',NDFBT_el , NDFT_el , NQS               , NDFBTWT_el
         write(*,'(a,4i6)') 'NDFBT0_el, NDFT0_el, NQS0',NDFBT0_el, NDFT0_el, NDFT0_el-NDFBT0_el, NDFBTWT0_el
      endif

!------ Include Modal damping
      call STRUCT_DAMP

!------ Solve the dynamic system of equations
      RLXs = 1.d0
      call Time_Integrate (it, 2, ICONV, MAXERR)
#                                                ifdef HAVE_OMP
                                                 t(10) = omp_get_wtime()
                                              if (NTIME==1) then
                                                  write(*,7)'rotmat',t( 2)-t(1)
                                                  write(*,7)'aero  ',t( 3)-t(2)
                                                  write(*,7)'hydro ',t( 4)-t(3)
                                                  write(*,7)'matrix',t( 5)-t(4)
                                                  write(*,7)'bloads',t( 6)-t(5)
                                                  write(*,7)'truss ',t( 7)-t(6)
                                                  write(*,7)'qs equ',t( 8)-t(7)
                                                  write(*,7)'reduct',t( 9)-t(8)
                                                  write(*,7)'newmar',t(10)-t(9)
                                                  write(*,7)'total ',t(10)-t(0)
                                                7 format (a6,f12.4)
                                              endif
#                                                endif
     endif !my_rank

#ifdef HAVE_MPI
 call MPI_BCAST(ICONV,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
      if (ICONV == 1) goto 1
   enddo !it

if(my_rank.eq.0) then
   write(* ,*) 'gelast not converged'
   write(10,*) 'gelast not converged'
 endif !my_rank

 1 continue

!! if (TIME>=2.d0) then 
!!    IEIG = 2
!!    call EIGEN_NEW
!!    stop
!! endif

#ifdef HAVE_MPI
   call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif

   ICONV = 1
   if     (IAERO_el == 2) then
#     ifdef HAVE_GENUVP
        !call GENUVP_main  ( INIT_gn,TIME_gn,DT_gn,Iter_el,IFull_gn,NTIME_gn,NTIMEM)
         call GENUVP_main  ( 0      ,TIME   ,DT_el, 1     ,ICONV   ,NTIME   ,9999, 9999 )
         call write_gn_load( TIME   ,NBLADE_el )
#     endif
   elseif (IAERO_el == 3) then
#     ifdef HAVE_cfd
      !conclude CFD step: write output & prepare next step
      !call CFD_main (WhatToDo, !=3 conlude the time step and prepare next step 
      !               NTIME   , !may start from 0
      !               TIME,DT , !time and dt are defined in hGAST
      !               NTIMEM)   !max number of steps in case some allocation needs it
         call broadcast_t
         call kinematics(NTIME)
         call CFD_main (3,NTIME,TIME,DT_el,NTIMEM,0,it) 
#     endif
   endif

 if (my_rank==0) then

!--- Update BEM solution parameters
   if (IAERO_el == 1) &
   call UPDATE_AERO (NTIME)

!--- Update control variables
   call control (NTIME,it,2)

!--- Update dynamic mooring line code variables
   if     (ITRUSS_el == 1) then

      call MOORINGS_UPDATE_tr

   elseif (ITRUSS_el == 2) then

      TIME_tr = TIME
      I1      = ACCU%NDFTACC_el (NBODT_el) + 1
      I2      = NDFBT_el

      UT_tr (1:NDFT_tr) = UT_el (I1:I2);
      UT1_tr(1:NDFT_tr) = UT1_el(I1:I2);
      UT2_tr(1:NDFT_tr) = UT2_el(I1:I2);
   endif


   call OUTPUT (NTIME)
 endif !my_rank


 END Subroutine GELAST
!----------------------------------------------------------------------
 Subroutine GELAST0
!----------------------------------------------------------------------

 Use Cbeam
#ifdef HAVE_MPI
 Use MPI
#endif

   implicit none

   integer :: my_rank, ierr


#ifdef HAVE_MPI
   call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)
#else
   ierr    = 0
   my_rank = 0
#endif

 if (my_rank==0) then 
!--- Modal analysis: Eigen value for each physical body with zero
!--- gravity, in order to determine the modal damping matrix
   call GELAST0_MODAL     !- After the call IMODAL becomes -1

!--- Eigen value for the whole structure
   call GELAST0_EIG

!--- set the initial values for basic parameters
   call INIT_DEFLECTIONS ( -OMEGA )
 endif !my_rank

!--- Compute initial deflections from static solution
   call GELAST_init (1)


 if(my_rank==0) &
   call ROTMAT_el
!qq flap
!  call controldef


 END Subroutine GELAST0
!----------------------------------------------------------------------
 Subroutine GELAST0_MODAL
!----------------------------------------------------------------------
!
!  Perform eigen-value calculations for each physical body, in order 
!  to estimate the structural [modal based] damping. No gravity or 
!  other external force are included, no rotational effects [omega=zero]
!  and the concentrated masses are not taken into account [i.e. nacelle,
!  hub, gen inertia].
!
!  IMODAL: 0  no modal analysis
!  IMODAL: 1  calculate modal damping
!  IMODAL: 2  read from file
!
!  Remarks 1. Complex numbers maybe appear for big systems or
!             huge stiffness or mass matrices. If this is the
!             case the code STOPS, printing an Error msg.
!          2. The modal damping is calculated wrt the local
!             body coordinate system [so the pitch angle doesn't
!             affect the structural damping]
!          3. It's not possible to reduce the number of the modes,
!             because the inverse of the modal matrix then is wrong!
!          5. For Tower foundation the free bc are not treated. Also 
!             for the jacket's free b.c.

!   IP=0 : call from GELAST0_EIG  [bc:AK]
!   IP=1 : call from GELAST       [bc:AM]
!
!----------------------------------------------------------------------

 Use Cbeam
 Use Hydro
 Use Paths

   implicit none

   real(8),save :: rho_tmp
   integer      :: it, ICONV, NDF, NDF_check, i, j, IP, nf
   integer      :: IBTYP_tmp(NBODT_el), IQfl_active_tmp(6,NFLOATER_el), IAERO_tmp


   if (IMODAL == 0) return

   if (IMODAL == 2) then
!------ Read from file
      IMODAL   =-1
      CDAMPSTR = 0.d0;

      open(200,file=trim(file_cstr)) !,form='UNFORMATTED')

       read(200,*) NDF, NDF_check
      !read(200  ) NDF, NDF_check

       if (NDF_check/=NDFT_el) then
          write(* ,*) 'incompatible CSTR.inp so modal damping is calculated again'
          write(10,*) 'incompatible CSTR.inp so modal damping is calculated again'
          IMODAL = 1
          goto 1
          stop
       endif

       write (* ,*) 'Reading CSTR.inp',NDF
       write (10,*) 'Reading CSTR.inp',NDF

       do j = 1, NDF
       do i = 1, NDF
          read(200,*)CDAMPSTR(i,j)
         !read(200  )CDAMPSTR(i,j)
       enddo
       enddo
      close (200)

      return
   endif

!--- Set the initial values for basic parameters
 1 call INIT_DEFLECTIONS ( 0.d0 ) !omega_set

!--- Disable the aero/hydro dynamic loading
   IBTYP_tmp(1:NBODT_el)          = subbody(1:NBODT_el)%IBTYP_el;
   IAERO_tmp                      = IAERO_el
   subbody  (1:NBODT_el)%IBTYP_el = 0;
   IAERO_el                       = 0
   GRAV                           = 0.d0
   rhoGrav                        = 0.d0
   rho_tmp                        = rho
   rho                            = 0.d0

   IP    = 0
   ICONV = 0

!--- Previous step solution
   do i=1,NDFT_el
      UTP_el (i) = UT_el (i)
      UTP1_el(i) = UT1_el(i)
      UTP2_el(i) = UT2_el(i)
   enddo


   it = 1
!------ Calculate rotation and translation matrices
      call ROTMAT_el

!------ Assembly of matrices
      call MATRIX_el
      call CON_MATRIX_el

!------ Loads communication between subbodies of same body
      call SBOD2SBOD_LOADS                                       !WT's subbodies (same body)
      call BOD2BOD_LOADSJack                                     !jacket's bodies

!------ Springs contribution for the foundation modelling
      call Foundation_contrib

!------ Equations for Qs
      if (ICASE_el >=3 ) then
         do nf = 1, NFLOATER_el
                         IQfl_active_tmp(1:6,nf) = floater(nf)% IQfl_active_el(1:6)      !disable floater's qs because modal analysis is
            floater(nf)% IQfl_active_el (1:6   ) = 0                                     !performed only for each body indipendently
         enddo
      endif

      call QS_EQUAT (IP)

      if (ICASE_el >=3 ) then
         do nf = 1, NFLOATER_el
            floater(nf)% IQfl_active_el(1:6) = IQfl_active_tmp(1:6,nf)
         enddo
      endif

!------ Reduction of Matrix
      call MATRIX_REDUCT (IP)
         write(*,*) 'NDFBT_el , NDFT_el , NQS ',NDFBT_el , NDFT_el , NQS
         write(*,*) 'NDFBT0_el, NDFT0_el, NQS0',NDFBT0_el, NDFT0_el, NDFT0_el-NDFBT0_el


!--- Modal Analysis with zero gravity for each physical body
!--- Determine Structural Damping Matrix
   if (IYNmodal==0) call MODAL_ANAL
#ifdef HAVE_MODAL
   if (IYNmodal==1) call MODAL_ANAL_md
#endif

!--- Restore saved vars
   subbody(1:NBODT_el)%IBTYP_el = IBTYP_tmp(1:NBODT_el);
   IAERO_el                     = IAERO_tmp
   GRAV                         = GRAVITY_el
   rho                          = rho_tmp
   rhoGrav                      = rho*gravity
   IMODAL                       =-1


 END Subroutine GELAST0_MODAL
!----------------------------------------------------------------------
 Subroutine GELAST0_EIG
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
 Use Hydro

   implicit none

   integer      :: it, ICONV, i, IP, ICMOD_tmp, iter
   character*80 :: outfil


   if ( IEIG == 0) return

   if (ICMOD >= 1) &
    write(*,*) "For the static calculation the free-fixed rotor will be simulated"

!--- set the initial values for basic parameters
   call INIT_DEFLECTIONS ( -OMEGA )


   if (IEIG == 1) then
!--- Zero the gravity and Buoyancy loads
         GRAV    = 0.d0     ! structural module
         gravity = 0.d0     ! hydro      module
      do i = 1, nfloater_hyd
         floater(i)% Moorings_weight = 0.d0
         floater(i)% Buoyancy_total  = 0.d0
      enddo

   endif
   if (IEIG == 2) GRAV  = GRAVITY_el
       IEIG  = 1

!qqsubbody  (1:NBODT_el)%IBTYP_el = 0;              !disable the hydrodynamic loading
   IAERO_el                       = 0               !disable the aerodynamic  loading
   ICMOD_tmp                      = ICMOD           !disable the controller/free-free
   ICMOD                          = 0
   IP                             = 0               !Static solution
   ICONV                          = 0
   iter                           = 25

!--- Previous step solution
   do i=1,NDFT_el
      UTP_el (i) = UT_el (i)
      UTP1_el(i) = UT1_el(i)
      UTP2_el(i) = UT2_el(i)
   enddo

   if (ITRUSS_el == 1) &
   call MOORINGS_UPDATE_tr


   do it = 1, iter

      if(NBODBTWT_el > 0) UT1_el(NDFBT_el+NQSW) = UTP1_el(NDFBT_el+NQSW)

!------ Calculate rotation and translation matrices
      call ROTMAT_el

!------ Hydrodynamic
      call Hydro_main ( it )

!------ Assembly of matrices
      call MATRIX_el
      call CON_MATRIX_el

!------ Loads communication between subbodies of same body
      call SBOD2SBOD_LOADS                                       !WT's subbodies (same body)
      call BOD2BOD_LOADSJack                                     !jacket's bodies

!------ Loads communication between bodies
      call BOD2BOD_LOADS(1,NBLADE_el)                            !blades to shaft
      call BOD2BOD_LOADS(NBLADE_el+1,NBLADE_el+1)                !shaft to tower
      call BOD2BOD_LOADS_WT_JA (NBODBT_el)                       !tower to jacket
      call BOD2BOD_LOADS_VAWT                                    !blade end to tower

!------ Moorings contribution for the floater
      call MOORINGS_COUPLED_tr ( 2, IP )                         !coupled with truss code

!------ Springs contribution for the foundation modelling
      call Foundation_contrib

!------ Equations for Qs
      call QS_EQUAT (IP)

!------ Reduction of Matrix
      call MATRIX_REDUCT (IP)
      if ( it == 1 ) then
         write(*,*) 'NDFBT_el , NDFT_el , NQS ',NDFBT_el , NDFT_el , NQS
         write(*,*) 'NDFBT0_el, NDFT0_el, NQS0',NDFBT0_el, NDFT0_el, NDFT0_el-NDFBT0_el
      endif

!------ Include Modal damping contribution
      call STRUCT_DAMP

!------ Only for the free-free eigen analysis, unfortunatelly atm no gravity could be included
      if ( ICONV == 2 ) goto 1
      if ( it == iter ) exit

!------ Solve the static system of equations
      RLXs = 0.80d0 !1.d0
      call Time_Integrate (it, 1, ICONV, 1.d-11)

      if (ICONV == 1) ICMOD = ICMOD_tmp
      if (ICONV == 1) ICONV = 2
   enddo !it

   write(* ,*) 'gelast0_eig not converged'
   write(10,*) 'gelast0_eig not converged'

 1 continue

!--- write the results from the static equilibrium
   call ROTMAT_el
                              outfil =  'geoG_init_eigen'
   call WRITEOUT_GLOB (UT_el, outfil, 0, 0)
                              outfil =  'geoL_init_eigen'
   call WRITEOUT_LOC  (UT_el, outfil, 0, 0)
   call write_loads
   call write_deform  (UT_el)
   call write_qsdof

!--- Solve the Eigenvalue problem and Write eigenvalues and modeshapes
   call EIGEN
!! call EIGEN_FLOATER

   stop


 END Subroutine GELAST0_EIG
!----------------------------------------------------------------------
 Subroutine GELAST_init (NTIME)
!----------------------------------------------------------------------

 Use Cbeam
 Use Hydro
 Use Craft
 Use Truss
#ifdef HAVE_MPI
 Use MPI
#endif

   implicit none

   integer      :: NTIME, i, it, ICONV, IP, nbod,nbsub,nb_el,nel,nn,  I1,I2, nf, itype
   character*80 :: outfil
   integer      :: IQfl_active_tmp(6,NFLOATER_el), ICMOD_tmp, IEIG_tmp, NSubTime_tmp, ITERMAX_tmp, ITER, IYN_Morison_tmp      !-- IYN_Morison_tr    [0: no hydro loads, 1: Morison's eq.]
   real(8)      :: ERROR
   integer      :: my_rank, ierr


#ifdef HAVE_MPI
   call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)
   call MPI_BCAST(ICALCINIT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#else
   ierr    = 0
   my_rank = 0
#endif

   if (ICALCINIT/=1.or.IRERUN==1) then
       ICALCINIT=-1
     return
   endif

 if (my_rank==0) then 

   write (* ,*)
   write (* ,*) 'NTIME = ', NTIME, TIME
   write (* ,*) 'Calc initial deflections'
   write (10,*)
   write (10,*) 'NTIME = ', NTIME, TIME
   write (10,*) 'Calc initial deflections'


!--- Previous step solution
   UTP_el (1:NDFT_el) = UT_el (1:NDFT_el);
   UTP1_el(1:NDFT_el) = UT1_el(1:NDFT_el);
   UTP2_el(1:NDFT_el) = UT2_el(1:NDFT_el);

   if (ITRUSS_el == 1) &
   call MOORINGS_UPDATE_tr

   IP              = 0
   ICONV           = 0
   IEIG_tmp        = IEIG
   ICMOD_tmp       = ICMOD
   IYN_Morison_tmp = IYN_Morison_tr
   IYN_Morison_tr  = 0
   ICMOD           = 0
   ITER            = 50
   if (ITRUSS_el == 2)&
   ITER            = 500
   ERROR           = 1.d-10
   if (IAERO_el  >= 2)&
   ERROR           = 1.d-4
   itype           = subbody(1)%IBTYP_el

 endif !my_rank

#ifdef HAVE_MPI
 call MPI_BCAST(itype,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
 call MPI_BCAST(ITER ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
 call MPI_BCAST(ICONV,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

   do it = 1, ITER
     if (my_rank==0) then 

      if(NBODBTWT_el > 0) UT1_el(NDFBT_el+NQSW) = UTP1_el(NDFBT_el+NQSW)

!!     call  Calc_Tow_Shaft_Angles

!------- Calculate rotation and translation matrices
       call ROTMAT_el
     endif !my_rank

!!    if (itype /= 0) then 
!---------- Calculate aerodynamic loads
         if (IAERO_el == 1) then
           if(my_rank.eq.0) then
!!          if ( dabs(UT1_el(NDFBT_el + NQSW)) < OMEGAR/3.50d0 ) IPARKED = 1
!!          if ( dabs(UT1_el(NDFBT_el + NQSW)) > OMEGAR/3.00d0 ) IPARKED = 0
            call RAFT     ( NTIME )
           endif!my_rank
         elseif (IAERO_el == 2) then
#          ifdef HAVE_GENUVP
             !call GENUVP_main  ( INIT_gn,TIME_gn,DT_gn,Iter_el,IFull_gn,NTIME_gn,NTIMEM)
              call GENUVP_main  ( 0      ,TIME   ,DT_el, 1     ,ICONV   ,NTIME   ,9999, 9999 )
#             ifdef HAVE_MPI
              if(my_rank==0)&
#             endif
              call Aero2Elast
#          endif !HAVE_GENUVP
         else if (IAERO_el == 3) then
#          ifdef HAVE_cfd
             !call CFD_main (WhatToDo, !=1 perform the "initial" step to get loads for the initial 
             !                             deformation / or read from a back-up file
             !               NTIME   , !may start from 0
             !               TIME,DT , !time and dt are defined in hGAST
             !               NTIMEM)   !max number of steps in case some allocation needs it
             if (it==1) then 
              call broadcast_t
              call kinematics(NTIME)
              call CFD_main (1,NTIME,TIME,DT_el,NTIMEM,-1,it) 
              call Aero2Elast_cfd
             endif
#          endif
         endif !IAERO_el
!!    endif !itype

 if (my_rank==0) then 
!------ Hydrodynamic
      call Hydro_main ( it )


!------ Assembly of matrices
      call MATRIX_el
      call CON_MATRIX_el

 
!------ Loads communication between subbodies of same body
      call SBOD2SBOD_LOADS                                    !WT's subbodies(same body)
      call BOD2BOD_LOADSJack                                  !jacket's bodies


!------ Loads communication between bodies
      call BOD2BOD_LOADS (1,NBLADE_el)                        !blades to shaft
      call BOD2BOD_LOADS (NBLADE_el+1,NBLADE_el+1)            !shaft  to tower
      call BOD2BOD_LOADS_WT_JA (NBODBT_el)                    !tower  to jacket
      call BOD2BOD_LOADS_VAWT                                 !blade end to tower


!------ Springs contribution for the foundation modeling
      call Foundation_contrib


!------ Moorings contribution for the floater
      if     (ITRUSS_el == 1) then
         NSubTime_tmp = NSubTime_tr
         NSubTime_tr  =  1
         ITERMAX_tmp  = ITERMAX_tr
         ITERMAX_tr   = 1001
         if (it==1) call calc_connection_UT_tr
         if (it==1) call MOORINGS_tr (TIME, 1, IP, 0.97d0)           !- for static solution RLX=0.8, No substeps, and 1001 max iterations!
                    call truss2float_loads_uncoupl
         NSubTime_tr  = NSubTime_tmp
         ITERMAX_tr   = ITERMAX_tmp
      elseif (ITRUSS_el == 2) then
         call calc_connection_UT_tr                                  !- [just for calculation of the position (writing)]
         call MOORINGS_COUPLED_tr ( 2, IP )                          !- coupled with truss code
      endif


!--- Equations for Qs
   if ((ICASE_el>=3).and.(ITRUSS_el>0)) then
!qqif  (ICASE_el>=3)                    then
      do nf = 1, NFLOATER_el
                      IQfl_active_tmp(1:6,nf) = floater(nf)% IQfl_active_el(1:6)  !disable floater's qs
         floater(nf)% IQfl_active_el (1:6   ) = 0                                 !because atm truss can't perform static solution
      enddo
   endif

      call QS_EQUAT (IP)
     if ((ICASE_el>=3).and.(ITRUSS_el>0)) then
!qq  if  (ICASE_el>=3) then
        do nf = 1, NFLOATER_el
           floater(nf)% IQfl_active_el(1:6) = IQfl_active_tmp(1:6,nf)
        enddo
     endif


!--- Reduction of Matrix
!! IEIG  = 1                     !-- not to keep truss elements dofs, for the initial calculations
      call MATRIX_REDUCT (IP)
!! IEIG  = IEIG_tmp

      if ( (NTIME==1).and.(it==1) ) then
         write(*,*) 'NDFBT_el , NDFT_el , NQS ',NDFBT_el , NDFT_el , NQS
         write(*,*) 'NDFBT0_el, NDFT0_el, NQS0',NDFBT0_el, NDFT0_el, NDFT0_el-NDFBT0_el
      endif

 endif !my_rank

#ifdef HAVE_MPI
      call MPI_BCAST(ICONV,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

      if ( (ICONV == 2).or.(it==ITER-1) ) goto 1
 if(my_rank==0) then
!------ Solve the dynamic system of equations
      RLXs = 1.d0 !0.8d0
!     call Time_Integrate (it, 1, ICONV, 1.d-11)
      call Time_Integrate (it, 1, ICONV, ERROR)
      if (ICONV == 1) ICONV = 2
 endif !my_rank
#   ifdef HAVE_MPI
      call MPI_BCAST(ICONV,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#   endif

   enddo !it

   write(* ,*) 'gelast Init not converged'
   write(10,*) 'gelast Init not converged'
   ICONV = 1

 1 continue


   if     (IAERO_el == 2) then
#     ifdef HAVE_GENUVP
        !call GENUVP_main  ( INIT_gn,TIME_gn,DT_gn,Iter_el,IFull_gn,NTIME_gn,NTIMEM)
         call GENUVP_main  ( 0      ,TIME   ,DT_el, 1     ,1       ,NTIME   ,9999, 9999 )
#     endif
   elseif (IAERO_el == 3) then
#     ifdef HAVE_cfd
        !conclude CFD step: write output & prepare next step
        !call CFD_main (WhatToDo, !=3 conlude the time step and prepare next step 
        !               NTIME   , !may start from 0
        !               TIME,DT , !time and dt are defined in hGAST
        !               NTIMEM)   !max number of steps in case some allocation needs it
         call broadcast_t
         call kinematics(NTIME)
         call CFD_main (3,NTIME,TIME,DT_el,NTIMEM,0,it) 
#     endif
   endif

 if (my_rank==0) then 

   if     (ITRUSS_el == 1) then

      call MOORINGS_UPDATE_tr

   elseif (ITRUSS_el == 2) then

      TIME_tr = TIME
      I1 = ACCU%NDFTACC_el (NBODT_el) + 1
      I2 = NDFBT_el

      UT_tr (1:NDFT_tr) = UT_el (I1:I2);
      UT1_tr(1:NDFT_tr) = UT1_el(I1:I2);
      UT2_tr(1:NDFT_tr) = UT2_el(I1:I2);

   endif
                              outfil =  'geoG_init_static'
   call WRITEOUT_GLOB (UT_el, outfil, NTIME, 0)
                              outfil =  'geoL_init_static'
   call WRITEOUT_LOC  (UT_el, outfil, NTIME, 0)
   call OUTPUT        (NTIME)
!! call control       (NTIME,it,2) !for output

   ICMOD          = ICMOD_tmp
   IYN_Morison_tr = IYN_Morison_tmp

   if (ICMOD==2) then

      nbod = NBLADE_el+1

      do nbsub=1,body(nbod)%NBODSUB_el
         nb_el = body(nbod)%NBODTGLB_el(nbsub)

         do nel = 1, subbody(nb_el)%NTEB_el
            nn       =                                subbody(nb_el)%NODMB       (nel,1)!nnod)
            i        = ACCU%NDFTACC_el(nb_el-1)   +   subbody(nb_el)%NDFPNBACC_el(nn-1 )
            UT_el (i+ 6) = 0.d0
!!!         UT_el (i+1:i+15) = 0.d0;

           if (subbody(nb_el)%NDFPE_el == 15) then
            nn       =                                subbody(nb_el)%NODMB       (nel,3)!nnod)
            i        = ACCU%NDFTACC_el(nb_el-1)   +   subbody(nb_el)%NDFPNBACC_el(nn-1 )
            UT_el (i+ 1) = 0.d0
           endif

            nn       =                                subbody(nb_el)%NODMB       (nel,NNPE_el)!nnod)
            i        = ACCU%NDFTACC_el(nb_el-1)   +   subbody(nb_el)%NDFPNBACC_el(nn-1 )
            UT_el (i+ 6) = 0.d0
         enddo
      enddo

   endif
 endif !my_rank

       ICALCINIT=-1


 END Subroutine GELAST_init
!-----------------------------------------------------------------------
 SUBROUTINE RERUN ( NTIME )
!-----------------------------------------------------------------------

 Use Cbeam
 Use Craft
 Use Hydro
 Use truss
 Use Caero
!Use Ctrl_mod

   implicit none

   integer :: NTIME, ISTRIP, i, j, iblad, nf


   if (mod((NTIME),(NTIMEBACK))/=0) return

   ISTRIP = NSTRIP*NBLADE_el

  open (2,file='RERUN.bk')

      write (2,100)  TIME
      write (2, * )  NTIME
      write (2,100) (UT_el   (  i),i=1,NDFT_el)
      write (2,100) (UT1_el  (  i),i=1,NDFT_el)
      write (2,100) (UT2_el  (  i),i=1,NDFT_el)
!     write (2,100) (QTC_el  (  i),i=1,NQSC0  )
!     write (2,100) (QTC1_el (  i),i=1,NQSC0  )
!     write (2,100) (QTC2_el (  i),i=1,NQSC0  )

   if (IDYNSTALL >  0) then
      write (2,100) (VELEFFP (3,i),i=1,ISTRIP )
      write (2,100) (VELEFFP (2,i),i=1,ISTRIP )
      write (2,100) (VELEFFP (1,i),i=1,ISTRIP )
      write (2,100) (ANGEFFP (3,i),i=1,ISTRIP )
      write (2,100) (ANGEFFP (2,i),i=1,ISTRIP )
      write (2,100) (ANGEFFP (1,i),i=1,ISTRIP )
   endif
                     
   if (IDYNSTALL == 1) then
      write (2,100) (CLOLDA  (1,i),i=1,ISTRIP )
      write (2,100) (CLOLDA  (2,i),i=1,ISTRIP )
      write (2,100) (CLOLDA  (3,i),i=1,ISTRIP )
      write (2,100) (CLOLDB  (1,i),i=1,ISTRIP )
      write (2,100) (CLOLDB  (2,i),i=1,ISTRIP )
      write (2,100) (CLOLDB  (3,i),i=1,ISTRIP )
      write (2,100) (CDOLDB  (1,i),i=1,ISTRIP )
      write (2,100) (CDOLDB  (2,i),i=1,ISTRIP )
      write (2,100) (CDOLDB  (3,i),i=1,ISTRIP )
      write (2,100) (CMOLDB  (1,i),i=1,ISTRIP )
      write (2,100) (CMOLDB  (2,i),i=1,ISTRIP )
      write (2,100) (CMOLDB  (3,i),i=1,ISTRIP )
      write (2,100) (CLIFTA  (1,i),i=1,ISTRIP )
      write (2,100) (CLIFTA  (2,i),i=1,ISTRIP )
      write (2,100) (CLIFTA  (3,i),i=1,ISTRIP )
      write (2,100) (CLIFTB  (1,i),i=1,ISTRIP )
      write (2,100) (CLIFTB  (2,i),i=1,ISTRIP )
      write (2,100) (CLIFTB  (3,i),i=1,ISTRIP )
      write (2,100) (CDRAGB  (1,i),i=1,ISTRIP )
      write (2,100) (CDRAGB  (2,i),i=1,ISTRIP )
      write (2,100) (CDRAGB  (3,i),i=1,ISTRIP )
      write (2,100) (CMOMB   (1,i),i=1,ISTRIP )
      write (2,100) (CMOMB   (2,i),i=1,ISTRIP )
      write (2,100) (CMOMB   (3,i),i=1,ISTRIP )
   elseif (IDYNSTALL == 2) then
      write (2,100) (BXOLD1  (1,i),i=1,ISTRIP )
      write (2,100) (BXOLD1  (2,i),i=1,ISTRIP )
      write (2,100) (BXOLD1  (3,i),i=1,ISTRIP )
      write (2,100) (BXOLD2  (1,i),i=1,ISTRIP )
      write (2,100) (BXOLD2  (2,i),i=1,ISTRIP )
      write (2,100) (BXOLD2  (3,i),i=1,ISTRIP )
      write (2,100) (BXOLD3  (1,i),i=1,ISTRIP )
      write (2,100) (BXOLD3  (2,i),i=1,ISTRIP )
      write (2,100) (BXOLD3  (3,i),i=1,ISTRIP )
      write (2,100) (BXOLD4  (1,i),i=1,ISTRIP )
      write (2,100) (BXOLD4  (2,i),i=1,ISTRIP )
      write (2,100) (BXOLD4  (3,i),i=1,ISTRIP )
   endif

   do iblad=1,NBLADE_el
      write (2,100) (strip(iblad,i)%AINDTB_PRE     , i=1,NSTRIP)
      write (2,100) (strip(iblad,i)%AINDPTB_PRE    , i=1,NSTRIP)
   enddo

   if (ICASE_el >= 3) then
      do nf = 1, NFLOATER_el
      do i  = 1, 6
         write (2,100) (floater(nf)%velocity_buffer (i,j),j=1,retardation_length-1)
      enddo
      enddo
      do i = 1, 6
      enddo

      if (ITRUSS_el == 1) then
         write (2,100) (UT_tr   (  i),i=1,NDFT_tr)
         write (2,100) (UT1_tr  (  i),i=1,NDFT_tr)
         write (2,100) (UT2_tr  (  i),i=1,NDFT_tr)
         write (2,100) ((connect_tr(i)% UT_CONNECT_tr (j), i=1,NCONNECT_tr), j=1,3)
         write (2,100) ((connect_tr(i)%UT1_CONNECT_tr (j), i=1,NCONNECT_tr), j=1,3)
         write (2,100) ((connect_tr(i)%UT2_CONNECT_tr (j), i=1,NCONNECT_tr), j=1,3)
      endif
   endif

   call control (NTIME,1 ,3)
!  write (2,100) TIME_BRAKE

  close(2)


 100 format (90000e28.17)
   

 END SUBROUTINE RERUN
!-----------------------------------------------------------------------
 SUBROUTINE RECALL ( NTIME )
!-----------------------------------------------------------------------

 Use Cbeam
 Use Craft
 Use Hydro
 Use truss
 Use Caero
!Use Ctrl_mod

  implicit none

  integer :: NTIME, ISTRIP, i, j, iblad, lstrip, nf


  if (IRERUN /= 1) then
     NTIME = 1
     return
  endif

  ISTRIP = NSTRIP*NBLADE_el

  write(* ,*)'READING BACKUP started'
  write(10,*)'READING BACKUP started'

  open (2,file='RERUN.bk')

      read (2,100)  TIME
      read (2, * )  NTIME              ! NTIME has been already performed
                    NTIME = NTIME+1
      read (2,100) (UT_el   (  i),i=1,NDFT_el)
      read (2,100) (UT1_el  (  i),i=1,NDFT_el)
      read (2,100) (UT2_el  (  i),i=1,NDFT_el)
!     read (2,100) (QTC_el  (  i),i=1,NQSC0  )
!     read (2,100) (QTC1_el (  i),i=1,NQSC0  )
!     read (2,100) (QTC2_el (  i),i=1,NQSC0  )

   if (IDYNSTALL >  0) then
      read (2,100) (VELEFFP (3,i),i=1,ISTRIP )
      read (2,100) (VELEFFP (2,i),i=1,ISTRIP )
      read (2,100) (VELEFFP (1,i),i=1,ISTRIP )
      read (2,100) (ANGEFFP (3,i),i=1,ISTRIP )
      read (2,100) (ANGEFFP (2,i),i=1,ISTRIP )
      read (2,100) (ANGEFFP (1,i),i=1,ISTRIP )
   endif
                     
   if (IDYNSTALL == 1) then
      read (2,100) (CLOLDA  (1,i),i=1,ISTRIP )
      read (2,100) (CLOLDA  (2,i),i=1,ISTRIP )
      read (2,100) (CLOLDA  (3,i),i=1,ISTRIP )
      read (2,100) (CLOLDB  (1,i),i=1,ISTRIP )
      read (2,100) (CLOLDB  (2,i),i=1,ISTRIP )
      read (2,100) (CLOLDB  (3,i),i=1,ISTRIP )
      read (2,100) (CDOLDB  (1,i),i=1,ISTRIP )
      read (2,100) (CDOLDB  (2,i),i=1,ISTRIP )
      read (2,100) (CDOLDB  (3,i),i=1,ISTRIP )
      read (2,100) (CMOLDB  (1,i),i=1,ISTRIP )
      read (2,100) (CMOLDB  (2,i),i=1,ISTRIP )
      read (2,100) (CMOLDB  (3,i),i=1,ISTRIP )
      read (2,100) (CLIFTA  (1,i),i=1,ISTRIP )
      read (2,100) (CLIFTA  (2,i),i=1,ISTRIP )
      read (2,100) (CLIFTA  (3,i),i=1,ISTRIP )
      read (2,100) (CLIFTB  (1,i),i=1,ISTRIP )
      read (2,100) (CLIFTB  (2,i),i=1,ISTRIP )
      read (2,100) (CLIFTB  (3,i),i=1,ISTRIP )
      read (2,100) (CDRAGB  (1,i),i=1,ISTRIP )
      read (2,100) (CDRAGB  (2,i),i=1,ISTRIP )
      read (2,100) (CDRAGB  (3,i),i=1,ISTRIP )
      read (2,100) (CMOMB   (1,i),i=1,ISTRIP )
      read (2,100) (CMOMB   (2,i),i=1,ISTRIP )
      read (2,100) (CMOMB   (3,i),i=1,ISTRIP )
   elseif (IDYNSTALL == 2) then
      read (2,100) (BXOLD1  (1,i),i=1,ISTRIP )
      read (2,100) (BXOLD1  (2,i),i=1,ISTRIP )
      read (2,100) (BXOLD1  (3,i),i=1,ISTRIP )
      read (2,100) (BXOLD2  (1,i),i=1,ISTRIP )
      read (2,100) (BXOLD2  (2,i),i=1,ISTRIP )
      read (2,100) (BXOLD2  (3,i),i=1,ISTRIP )
      read (2,100) (BXOLD3  (1,i),i=1,ISTRIP )
      read (2,100) (BXOLD3  (2,i),i=1,ISTRIP )
      read (2,100) (BXOLD3  (3,i),i=1,ISTRIP )
      read (2,100) (BXOLD4  (1,i),i=1,ISTRIP )
      read (2,100) (BXOLD4  (2,i),i=1,ISTRIP )
      read (2,100) (BXOLD4  (3,i),i=1,ISTRIP )
   endif

   do iblad=1,NBLADE_el
      read (2,100) (strip(iblad,i)%AINDTB_PRE     , i=1,NSTRIP)
      read (2,100) (strip(iblad,i)%AINDPTB_PRE    , i=1,NSTRIP)
   enddo


   if (ICASE_el >= 3) then
      do nf = 1, NFLOATER_el
      do i  = 1, 6
         read (2,100) (floater(nf)%velocity_buffer (i,j),j=1,retardation_length-1)
      enddo
      enddo

      if (ITRUSS_el == 1) then
         read (2,100) (UT_tr   (  i),i=1,NDFT_tr)
         read (2,100) (UT1_tr  (  i),i=1,NDFT_tr)
         read (2,100) (UT2_tr  (  i),i=1,NDFT_tr)
         read (2,100) ((connect_tr(i)% UT_CONNECT_tr (j), i=1,NCONNECT_tr), j=1,3)
         read (2,100) ((connect_tr(i)%UT1_CONNECT_tr (j), i=1,NCONNECT_tr), j=1,3)
         read (2,100) ((connect_tr(i)%UT2_CONNECT_tr (j), i=1,NCONNECT_tr), j=1,3)
      endif
   endif

   call control (NTIME,1 ,4)
!  read (2,100) TIME_BRAKE

   write(* ,*)'READING BACKUP finished'
   write(10,*)'READING BACKUP finished'


   call ROTMAT_el


   do iblad = 1, NBLADE_el
      do lstrip = 1,NSTRIP
         strip(iblad,lstrip)%AINDTB  = strip(iblad,lstrip)%AINDTB_PRE
         strip(iblad,lstrip)%AINDPTB = strip(iblad,lstrip)%AINDPTB_PRE
      enddo
   enddo

  close(2)
                     
 100 format (90000e28.17)
   

 END SUBROUTINE RECALL
!------------------------------------------------------------------------
 Subroutine get_tower_ref ( RG_tow, AT_tow, DRG_tow )
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   real(8),intent (out) :: RG_tow(3), AT_tow(3,3), DRG_tow(3)


   if (ICASE_el >= 3) then
      RG_tow(1:3    ) =  R_float(1:3    );
     DRG_tow(1:3    ) = DR_float(1:3    );
      AT_tow(1:3,1:3) = AT_float(1:3,1:3);
   else
      RG_tow(1:3    ) =  0.d0;
     DRG_tow(1:3    ) =  0.d0;
     call DIAGO (3,AT_tow)
   endif


 END Subroutine get_tower_ref
