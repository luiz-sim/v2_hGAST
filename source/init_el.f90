 include 'jacket_init.f90'
!-----------------------------------------------------------------------------------
! 
!     Subroutine :INITIA_el
!     
!     Open input FILES, perform initial Calculations
!
!     dofs 1:NDFBT_el are equations of virtual works for the elastic dofs
!                     starting from the 1st body. In case of 15dof elements
!                     with 5 nodes [3 internal nodes] the order is:
!                     {ux,θz,uy,uz,θx,θy}, {uy}, {θy}, {uy}, {ux,θz,uy,uz,θx,θy}
!                           nod1           nod2  nod3  nod4        nod5
!     dofs NDFBT_el+1:NDFTB_el+NQS=NDFT_el are additional equations which 
!                    define the kinematic qs.
!
!          1:NQSP are the primary qs {omega,yaw, {qx,qy,qz,θx,θy,θz}_floater }
!     NQSP+1:NQS  are the subbody endnodes local qs in the order {qx,qy,qz,θx,θy,θz}
!
!
!-----------------------------------------------------------------------------------
 Subroutine INITIA_el
!-----------------------------------------------------------------------------------

 Use Cbeam
 Use Craft
 Use Hydro
 Use Paths
#ifdef HAVE_MPI
 Use MPI
#endif
#ifdef HAVE_OMP
 Use omp_lib
#endif

   implicit none

!--- dfile vars
   integer, parameter   :: IBT=10, ISBT=100
   integer              :: NBODSUB_wt(IBT     ), IDOF_wt  (IBT     )
   integer              :: IBTYP_wt  (IBT,ISBT), NTEB_wt  (IBT,ISBT)
   integer              :: IWRITE    (IBT,2   )
   real(8)              :: PHI0_wt   (IBT     ), PITCH_wt (IBT     ), PITCH_Imb_wt(  IBT)
   real(8)              ::                       PITCHC_wt(IBT     ), PITCHS_wt   (  IBT)
   real(8)              :: COEFM_wt  (IBT     ), COEFK_wt (IBT     ), CCRIT_wt    (0:IBT)

   real(8)              :: rho1, depth1
   integer              :: i, IDT
   character*80         :: outfil
   integer              :: p, ierr, nf, my_rank
   integer, allocatable :: NBODSUB_in(:)



#ifdef HAVE_MPI
   call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)
   call MPI_Comm_size(MPI_COMM_WORLD,p      ,ierr)
#else
   ierr    = 0
   p       = 1
   my_rank = 0
#endif

#ifdef HAVE_OMP
!!call mkl_set_dynamic (.false.) !to use all threads (slow)
!!call mkl_set_dynamic (.true.)
! call omp_set_num_threads(1)
! call mkl_set_num_threads(1)
!!call mkl_domain_set_num_threads(1)
#endif

 if (my_rank==0) then 
!--- open the screen output file
   open (10, file='SCR.out',access='append')
   call Runinfo (1)
   write(10,'(a,i4)') 'MPI processors  : ', p
   write(10,*)

#  ifndef    OS
#     define OS 1
#  endif


!--- read paths, write exe path and working path
   call read_paths

   BIG          = 1.d16
   TIME         = 0.d0
   PI           = dacos(-1.d0)
   PI2          = 2.d0*PI
   PIhalf       = PI/2d0
   R2D          = 180.d0/PI
   RPM2RAD      = PI/30.d0
   GRAVITY_el   = 9.80665d0
   GRAV         = GRAVITY_el
   IYNmodal_ut  = 0           ! switch to FEM solution[0] or modal solution[1] for elastic deflection, velocity, acceleration in sb LOCAL_UT_1


!--- Switch FEM equations/variables (FX,MZ,FY,FZ,MX,MY) to "normal" order (FX,FY,FZ,MX,MY,MZ)
   IMDOF_el(1)  = 1  !u
   IMDOF_el(2)  = 3  !v
   IMDOF_el(3)  = 4  !w
   IMDOF_el(4)  = 5  !ϑx
   IMDOF_el(5)  = 6  !ϑy
   IMDOF_el(6)  = 2  !ϑz


   call read_dfile ( COEFM_wt, COEFK_wt  , CCRIT_wt, PHI0_wt , PITCH_wt, PITCH_Imb_wt, PITCHC_wt, PITCHS_wt, &
                               NBODSUB_wt, IDOF_wt , IBTYP_wt, NTEB_wt , IBT         ,ISBT      , IWRITE      )        !-- Allocate floater_structure

   call read_jacket_nod_bod

   call allocate_main_matrices ( COEFM_wt  , COEFK_wt , CCRIT_wt    ,&
                                 PHI0_wt   , PITCH_wt , PITCH_Imb_wt,&
                                             PITCHC_wt, PITCHS_wt   ,&
                                 NBODSUB_wt, IDOF_wt  , IBTYP_wt    ,&
                                 NTEB_wt   , IBT      , ISBT        , IWRITE )

!--- default values for Var_Paper
   Var_Paper( 1:20) = 0.d0
   Var_Paper( 4: 8) = 1.d0
   Var_Paper(11:14) = 1.d0

   open(30,file='paper.inp')
     read(30,*,end=1,err=1)  ! **Fxct      FxAmpl      omega           ID_TT          ID_BT      ID_Xcg      ID_EAx         ID_GRAV      FxLoc      IAX F         EAx_coeff       EAz_coeff        xcm_coeff    zcm_coeff       EIxz_coeff     EIxy_coeff
     read(30,*,end=1,err=1)  Var_Paper(1),Var_Paper(2),Var_Paper(3),Var_Paper(4),Var_Paper(5),Var_Paper(6),Var_Paper(7),Var_Paper(8),Var_Paper(9), Var_Paper(10), Var_Paper(11), Var_Paper(12), Var_Paper(13), Var_Paper(14), Var_Paper(15), Var_Paper(16)
   close(30)                                                                                                                         !1 :x, 2:y, 3, 10:combined, 20:thy,
                                                                                                                                     !11:tmp gamesa, 12:flap at tip, 13:My(torsion) at tip
   goto 2

1  Write(* ,'(a)') 'file paper.inp not found, so the default values were set'
   Write(10,'(a)') 'file paper.inp not found, so the default values were set'
#  if   OS == 0
       call system ('rm -f paper.inp')
#  elif OS == 1
       call system ('del paper.inp')
#  endif

2  call define_external_force
   call read_machi

!--- JACKET
   call Jacket_trans
   call read_jacket_machi_bc
   call INIT_Atrans

!--- Foundation
   call Foundation_Init
   call manage_con_mass
   call calc_htow_hsh

!--- Initialize (deformations, Forces and their derivatives )----- 
!--- Initialize other elastic VARIABLES
   FSTRIP_el  = 0.d0;
   XCAER_el   = 0.d0;
   CDAMPSTR   = 0.d0;
   UT_el      = 0.d0;  UTP_el  = 0.d0;  UT0_el  = 0.d0;
   UT1_el     = 0.d0;  UTP1_el = 0.d0;  UT01_el = 0.d0;
   UT2_el     = 0.d0;  UTP2_el = 0.d0;  UT02_el = 0.d0;

!--- Set the initial values for basic parameters
   call INIT_DEFLECTIONS ( 0.d0 ) !omega_set

!--- Initialize Control Variables
   call control_init_el

!--- Initialize the rotate package for calculating the Tranformation Matrices
   Allocate (NBODSUB_in(NBODBT_el)); NBODSUB_in(1:NBODBT_el)=body(1:NBODBT_el)%NBODSUB_el;
   call Init_Rot  ( NBODSUB_in, NBODBT_el, NBODBTWT_el, NBLADE_el, IQfl_rot, IQfl_tr )
   call Alloc_Rot ( NQS )
   Deallocate (NBODSUB_in)
 endif !my_rank

#ifdef HAVE_MPI
 call MPI_BCAST(IAERO_el  ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
 call MPI_BCAST(DT_el     ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
 call MPI_BCAST(TPERIOD_el,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
 call MPI_BCAST(TIME      ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
 call MPI_BCAST(NTIMEM    ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
#endif

!--- Initialize AERO / calculate Hub position
     if (my_rank==0) then
      INDXFlap      = 0
      FLCONTR (:,:) = 0.d0;
      call Calc_HubPos
  !!r call INIT_DEFLECTIONS ( -OMEGA ) !omega_set
      call ROTMAT_el
      call Set_undeformed_blade
      if (IAERO_el == 1) call RAFT_init                 !ROTMAT_el needed for ELST2AER_init
     endif !my_rank

#ifdef HAVE_MPI
 call MPI_BARRIER(MPI_COMM_WORLD, IERR)
#endif

   if (IAERO_el == 2) then
      IDT = TPERIOD_el/DT_el
#    ifdef HAVE_GENUVP
     !call GENUVP_main  (INIT_gn,TIME_gn,DT_gn,Iter_el,IFull_gn,NTIME_gn,NTIMEM,IDT)
      call GENUVP_main  (1      ,TIME   ,DT_el,1      , 1      ,0       ,NTIMEM,IDT) !   TIME=0.d0; ntime0=0; INIT_gn=1; ICONV=1
#    else
      write(*,*) 'Compiled without GENUVP support'; stop
#    endif
   else if (IAERO_el == 3) then
#    ifdef HAVE_cfd
     !call CFD_main (WhatToDo, !=0 initialize ... call INIT_cfd 
     !               NTIME   , !may start from 0
     !               TIME,DT , !time and dt are defined in hGAST
     !               NTIMEM)   !max number of steps in case some allocation needs it
      call broadcast_0
      call CFD_main (0,0,TIME,DT_el,NTIMEM,0,0) 
      call kinematics_0
     !in this call CFD is initialized
#    else
      write(*,*) 'Compiled without CFD support'; stop
#    endif
   endif
  !!r call INIT_DEFLECTIONS ( 0.d0   ) !omega_set
  !!r call ROTMAT_el

 if(my_rank==0) then
!--- Initialize Aerodynamic drag for the nacelle and the tower
   HUB_VEL (1:3) = 0.d0

   call  DRAG_INIT_el

!--- Initialize Hydro
   call Hydro_Init (rho1, depth1) !tide_el,Zfloat_offset

!--- Initialize Moorings
   do nf = 1, NFLOATER_el
      floater(nf)% Q_truss (1:6) = 0.d0
   enddo

   if (ITRUSS_el /= 0) &
   call INIT_tr (GRAVITY_el,rho1,depth1,DT_el,file_truss)

!--- Prepare matrices in case of coupled truss code
   call MOORINGS_COUPLED_tr ( 1, i )

!--- Calculate total masses for each component
   call Calc_Masses

!--- Output the initial undeformed geometry
   call ROTMAT_el
   outfil = 'geoG_'
   call WRITEOUT_GLOB ( UT_el, outfil, 0, 5)

   call VAWT_bc_init                               ! Initialize VAWT b.c.
 endif !my_rank


 END Subroutine INITIA_el
!-----------------------------------------------------------
 Subroutine VAWT_bc_init
!-----------------------------------------------------------

 Use Cbeam

   implicit none

   integer :: nbod, nbsub, nb_el, nn_el, nbc, nnT2_el


   if (IAPPL_el /= 1) return!0:HAWT, 1:VAWT

!--- Initialize boundary conditions
   do nbod  = 1, NBLADE_el
      nbsub = body(nbod)%NBODSUB_el
      nb_el = body(nbod)%NBODTGLB_el(nbsub  )
      nn_el = body(nbod)%NNODTGLB_el(nbsub+1)

                        nbc = 2
      boundco (nb_el)%nbcpb = nbc
   
      Deallocate ( boundco (nb_el)%nod   ); Allocate ( boundco (nb_el)%nod    (nbc    ) )
      Deallocate ( boundco (nb_el)%nel   ); Allocate ( boundco (nb_el)%nel    (nbc    ) )
      Deallocate ( boundco (nb_el)%nbcon ); Allocate ( boundco (nb_el)%nbcon  (nbc    ) )
                                            Allocate ( boundco (nb_el)%nelcon (nbc    ) )
                                            Allocate ( boundco (nb_el)%nodcon (nbc    ) )
      Deallocate ( boundco (nb_el)%indx  ); Allocate ( boundco (nb_el)%indx   (nbc,6  ) )
                                            Allocate ( boundco (nb_el)%Atrans1(nbc,3,3) ) !b.c.
!!                                          Allocate ( boundco (nb_el)%Atrans2(nbc,3,3) ) !loads
   
      nbc = 1
                 boundco (nb_el)%nod   (nbc    ) = 1  ! b.c. at the 1st node
                 boundco (nb_el)%nel   (nbc    ) = 1  !   of the 1st element
                 boundco (nb_el)%nbcon (nbc    ) = 0  ! zero condition
                 boundco (nb_el)%indx  (nbc,1:6) = 1  ! 1:b.c., 0:free
      nbc = 2
                 boundco (nb_el)%nod   (nbc    ) = NNPE_el                  ! b.c. at the last node
                 boundco (nb_el)%nel   (nbc    ) = subbody(nb_el)%NTEB_el   !   of the last element
                 boundco (nb_el)%nbcon (nbc    ) = body(NTOWER_el)%NBODTGLB_el ( body(NTOWER_el)%NBODSUB_el )
                 boundco (nb_el)%nelcon(nbc    ) = subbody(body(NTOWER_el)%NBODTGLB_el(body(NTOWER_el)%NBODSUB_el))%NTEB_el
                 boundco (nb_el)%nodcon(nbc    ) = NNPE_el
                 boundco (nb_el)%indx  (nbc,1:6) = 1

                 nnT2_el = body(NTOWER_el)%NNODTGLB_el (body(NTOWER_el)%NBODSUB_el+1 ) !last node of tower: ttop

                 boundco (nb_el)%Atrans1(nbc,1:3,1:3) = matmul (transf_mat(nnT2_el)%AT_el(1:3,1:3), transf_mat(nn_el)%A_el(1:3,1:3))

!     write(*,*         ) 'VAWT Atrans1',nbod
!     write(*,'(3f15.7)') boundco (nb_el)%Atrans1(nbc,1  ,1:3)
!     write(*,'(3f15.7)') boundco (nb_el)%Atrans1(nbc,2  ,1:3)
!     write(*,'(3f15.7)') boundco (nb_el)%Atrans1(nbc,3  ,1:3)
   enddo


 END Subroutine VAWT_bc_init
!-----------------------------------------------------------
 Subroutine INIT_DEFLECTIONS (omega_set)
!-----------------------------------------------------------

 Use Cbeam
 Use Hydro

   implicit none

   real(8) :: omega_set, Ishaft
   integer :: nbod, NQACC, i, nf, ii


!--- Zero out elastic deflections
   UT_el   = 0.d0;
   UT1_el  = 0.d0;
   UT2_el  = 0.d0;

   Ishaft  = 1.d0
   if (IAPPL_el == 1 .or. IAPPL_el == 2) Ishaft  = -1.d0 !0:HAWT, 1:VAWT, 2:Heli

!--- initial rotational speed
   if(NBODBTWT_el > 0) UT1_el(NDFBT_el+NQSW) = omega_set * Ishaft !-OMEGA

!--- initial blade pitch angles
   do nbod  = 1, NBLADE_el
      NQACC = ACCU%NQSACC_el(nbod-1) 
      UT_el (NDFBT_el+NQSP+NQACC+5) = PITCH0_el(nbod) + PITCH_Imb_el(nbod)
   enddo

!--- initial floater's dofs
!--- if truss moorings are enabled special care is needed
!--- in order to set the correct connection points.
   do nf  = 1, NFLOATER_el
    do i  = 1, IQfl
       ii = (nf-1)*6 + i
       UT_el (NDFBT_el+ii) = floater(nf)% UflInit (i)
    enddo
   enddo


 END Subroutine INIT_DEFLECTIONS
!-----------------------------------------------------------
 Subroutine read_paths
!-----------------------------------------------------------
 Use Paths

   implicit none


   open (1 , file='paths.inp')

                                    write(10,*)
    read(1,*) file_dfile_el   ;     write(10,*)'file_dfile_el  :', trim(file_dfile_el  )  !dfile_el.inp
    read(1,*) file_rotorin    ;     write(10,*)'file_rotorin   :', trim(file_rotorin   )  !rotorin.inp
    read(1,*) file_hydro      ;     write(10,*)'file_hydro     :', trim(file_hydro     )  !hydro.inp
    read(1,*) file_jacket     ;     write(10,*)'file_jacket    :', trim(file_jacket    )  !jacket.inp
    read(1,*) file_machi      ;     write(10,*)'file_machi     :', trim(file_machi     )  !machi.inp
    read(1,*) file_geomb      ;     write(10,*)'file_geomb     :', trim(file_geomb     )  !geomb.inp
    read(1,*) file_profilb    ;     write(10,*)'file_profilb   :', trim(file_profilb   )  !profilb.inp
    read(1,*) dir_airfoils    ;     write(10,*)'dir_airfoils   :', trim(dir_airfoils   )  !
    read(1,*) file_drag       ;     write(10,*)'file_drag      :', trim(file_drag      )  !drag.inp
    read(1,*) file_floater    ;     write(10,*)'file_floater   :', trim(file_floater   )  !floater.inp
    read(1,*) file_diffract   ;     write(10,*)'file_diffract  :', trim(file_diffract  )  !diffract.inp
    read(1,*) file_retard     ;     write(10,*)'file_retard    :', trim(file_retard    )  !retard.inp
    read(1,*) file_morison    ;     write(10,*)'file_morison   :', trim(file_morison   )  !morison.inp
    read(1,*) file_truss      ;     write(10,*)'file_truss     :', trim(file_truss     )  !truss.inp
    read(1,*) file_wind       ;     write(10,*)'file_wind      :', trim(file_wind      )  !INWIND.DAT
    read(1,*) dir_wave        ;     write(10,*)'dir_wave       :', trim(dir_wave       )  !
    read(1,*) file_foundation ;     write(10,*)'file_foundation:', trim(file_foundation)  !foundation.inp
    read(1,*) file_cstr       ;     write(10,*)'file_cstr      :', trim(file_cstr      )  !
    read(1,*) file_dllname    ;     write(10,*)'file_dllname   :', trim(file_dllname   )  !
    read(1,*) file_dllcfg     ;     write(10,*)'file_dllcfg    :', trim(file_dllcfg    )  !

   close(1)


 END Subroutine read_paths
!-----------------------------------------------------------
 Subroutine Calc_HubPos
!-----------------------------------------------------------

 Use Cbeam
 Use Hydro

   implicit none

   real(8) :: UflInit_tmp (6,NFLOATER_el)
   integer :: nbod, nbsub, nn_el, i, nf


   if (NBODBTWT_el > NBLADE_el) then

      do nf = 1, NFLOATER_el
                      UflInit_tmp (1:6,nf) = floater(nf)% UflInit (1:6);
         floater(nf)% UflInit     (1:6   ) = 0.d0;
      enddo

      call INIT_DEFLECTIONS ( 0.d0 )

      call ROTMAT_el  !-- only to calculate HUBpos_el

      do nf = 1, NFLOATER_el
         floater(nf)% UflInit     (1:6   ) = UflInit_tmp (1:6,nf)
      enddo

      call INIT_DEFLECTIONS ( 0.d0 )

      nbod    = NSHAFT_el
      nbsub   = body(nbod)%NBODSUB_el+1
      nn_el   = body(nbod)%NNODTGLB_el(nbsub)
      HUBpos_el(1:3) = transf_mat(nn_el)%R_el(1:3)

   else               !-- here the YAW angle isn't taken into account. In any case it's useless, because no tower is considered

      HUBpos_el (1) = -Hsh*dcos(TILT)
      HUBpos_el (2) =  0.d0
      HUBpos_el (3) =  Htow0 + Hshof + Hsh*dsin(TILT)

   endif

   write(10,*) 'HUB Position'
   write(10,*) (HUBpos_el(i),i=1,3)


 END Subroutine Calc_HubPos
!----------------------------------------------------------------------
 Subroutine read_dfile ( COEFM_wt  , COEFK_wt    , CCRIT_wt , PHI0_wt  , &   
                         PITCH_wt  , PITCH_Imb_wt, PITCHC_wt, PITCHS_wt, &
                         NBODSUB_wt, IDOF_wt     , IBTYP_wt , NTEB_wt  , &
                         IBT       , ISBT        , IWRITE                 )
!----------------------------------------------------------------------

 Use Cbeam
 Use Hydro
 Use Paths

   implicit none

   integer :: IBT, ISBT
   integer :: NBODSUB_wt(IBT     ), IDOF_wt  (IBT     )
   integer :: IBTYP_wt  (IBT,ISBT), NTEB_wt  (IBT,ISBT)
   integer :: IWRITE    (IBT,2   )
   real(8) :: PHI0_wt   (IBT     ), PITCH_wt (IBT     ), PITCH_Imb_wt(  IBT)
   real(8) ::                       PITCHC_wt(IBT     ), PITCHS_wt   (  IBT)
   real(8) :: COEFM_wt  (IBT     ), COEFK_wt (IBT     ), CCRIT_wt    (0:IBT)
!--- local vars
   real(8) :: TIME_MAX
   integer :: nbod, nbsub, ncal, i, nf


   open ( 11, file=trim(file_dfile_el) )

   read (11, *)!** title 
   read (11, *)!** Application Parameters
     read (11, *) ICASE_el      !; write(*,*) 'icase_el ', ICASE_el       ! [0:onshore, 1:monopile, 2:tripod/jacket, 3:floating, 4:flex floater] 30,40-->2 floaters-->12dofs
     read (11, *) IAPPL_el      !; write(*,*) 'iappl_el ', IAPPL_el       ! [0:HAWT, 1: VAWT, 2:Helicopter]
     read (11, *) IYNmodal      !; write(*,*) 'iynmodal ', IYNmodal       ! [0:FEM       , 1:MODAL]
     read (11, *) ITRUSS_el     !; write(*,*) 'itruss_el', ITRUSS_el      ! [0:no moorings, 1:moorings implimented with truss elements]
     read (11, *) IEIG          !; write(*,*) 'ieig     ', IEIG           ! eigen-value   flag: 0:no (time domain), 1:eigen value without gravity, 2:eigen value including gravity
     read (11, *) IMODAL        !; write(*,*) 'imodal   ', IMODAL         ! modal damping flag: 0: no, 1:calculate modal damping, 2:read modal damping from file]
     read (11, *) IAERO_el      !; write(*,*) 'iaero_el ', IAERO_el       ! [0:no      , 1:RAFT                 , 2:GENUPV]
     read (11, *) IORDbl_el     !; write(*,*) 'iordbl_el', IORDbl_el      ! blade beam modelling [1:1st order Timoskenko beam, 2:2nd order Euler-Bernoulli beam ]
     read (11, *) ITYPE_FOUND_el!; write(*,*) 'ifound_el', ITYPE_FOUND_el ! foundation type [0:AF, 1:CS, 2:DS, 3:DS p-y]
     read (11, *) ICMOD         !; write(*,*) 'icmod    ', ICMOD          ! Control Mode    [ 0:omega and pitch fixed, 1:controller, 2:WT iddling ]

     if (IORDbl_el/=1) then; write(*,*)'In current implementation only the 1st-order Timoshenko beam is supported';stop;endif

     if (IYNmodal==1) then
#ifdef HAVE_MODAL
#else
     write(*,*) 'in current implementation modal part is disabled'; stop
#endif
     endif

!----- Define the number of floaters
     if (ICASE_el==3.or.ICASE_el==4) then
        NFLOATER_el = 1
     else
        NFLOATER_el = 0
     endif

     if (IAPPL_el == 2) &
        NFLOATER_el = 1

   read (11, *)!** Bodies/Sub-Bodies/Elements Parameters
     read (11, *) NBODBTWT_el   !; write(*,*) 'nbodbtwt_el', NBODBTWT_el ! Total Number of WT's     Bodies [jacket not included]
     read (11, *) NBLADE_el     !; write(*,*) 'nblade_el  ', NBLADE_el   ! Total Number of blades

     if (NBLADE_el==NBODBTWT_el.and.(ICMOD/=0.or.ICMOD/=2)) then;write(10,*)'ICMOD is automatically set to 0';endif;


     IWRITE (:,1) = 0;
     IWRITE (:,2) = 1; !by default output only the local loads/deformations
     NSUBperBM_el = 1
     do nbod = 1, NBODBTWT_el
        if ( (nbod == 1).or.(nbod > NBLADE_el) ) then
           read (11,*,err=1)   &
            NBODSUB_wt(nbod  ),&                                        ! Number of Sub-bodies per Body
            IDOF_wt   (nbod  ),&                                        ! Elasticity flag [0: Elastic dofs disabled, 1: enabled]
            IWRITE    (nbod,1),&                                        ! Output loads, deform in global c.s.
            IWRITE    (nbod,2)                                          ! Output loads, deform in local  c.s.

 1          NSUBperBM_el = max( NSUBperBM_el, NBODSUB_wt(nbod) )
        else
            NBODSUB_wt(nbod    ) = NBODSUB_wt(1    )
            IDOF_wt   (nbod    ) = IDOF_wt   (1    )
            IWRITE    (nbod,1:2) = IWRITE    (1,1:2)
        endif
     enddo

     do nbod=1,NBODBTWT_el
       do nbsub=1,NBODSUB_wt(nbod)

         if ( (nbod == 1).or.(nbod > NBLADE_el) ) then
           read (11, *)           &
            IBTYP_wt(nbod,nbsub), &                                     ! 0 no loads, 1 aero loads, 2 hydro loads
            NTEB_wt (nbod,nbsub)                                        ! Elements per subbody
           if ( (IAERO_el==0).and.(IBTYP_wt(nbod,nbsub)==1) ) &
                                   IBTYP_wt(nbod,nbsub) = 0             ! no aero if IAERO=0
           if ( (ICASE_el==0).and.(IBTYP_wt(nbod,nbsub)==2) ) &
                                   IBTYP_wt(nbod,nbsub) = 0             ! no hydro for onshore WT
         else
           IBTYP_wt(nbod,nbsub) = IBTYP_wt(1,nbsub)
           NTEB_wt (nbod,nbsub) = NTEB_wt (1,nbsub)
         endif

       enddo !nbsub
     enddo !nbod

!----- Jacket's parameters
     if ((ICASE_el==2).or.(ICASE_el==4)) &
           read(11,*)                    &
            IBTYP_wt(NBLADE_el+3, 1)    ,&                              ! ITYPE_ja
            IDOF_wt (NBLADE_el+3   )    ,&                              ! IDOF_ja
            IWRITE  (NBLADE_el+3, 1),&                                  ! Output loads, deform in global c.s.
            IWRITE  (NBLADE_el+3, 2)                                    ! Output loads, deform in local  c.s.


           if ( (ICASE_el==0).and.(IBTYP_wt(NBLADE_el+3, 1)==2) ) &
                                   IBTYP_wt(NBLADE_el+3, 1)= 0          ! no hydro for onshore WT


   read (11, *)!**  Initial Conditions
     read (11, *) ICALCINIT     !; write(*,*) 'icalcinit', ICALCINIT     ! Calculate initial deform state from static solution [ 0:no, 1:yes ]
     if (NBLADE_el > 0) then
        read (11, *) PHI0_wt(1) !; write(*,*) 'phi0     ', PHI0_wt(1)    ! Initial Azimuthal position of Blade 1[deg]
                     PHI0_wt(1) = PHI0_wt(1) / R2D
     endif
                  do nbod=2,NBLADE_el
                     PHI0_wt(nbod) = PHI0_wt(nbod-1) + PI2/NBLADE_el
                     if (PHI0_wt(nbod) >  PI2) PHI0_wt(nbod) = PHI0_wt(nbod) - PI2
                     if (PHI0_wt(nbod) < 0.d0) PHI0_wt(nbod) = PHI0_wt(nbod) + PI2
                  enddo

     do nbod = 1, NBLADE_el 
        if (IAPPL_el<2) then
          PITCHC_wt(:)=0.d0; PITCHS_wt(:)=0.d0;
        read (11, *)          &
          PITCH_wt    (nbod) ,&                                         ! Initial Pitch   angle of Blades  [deg]
          PITCH_Imb_wt(nbod)                                            ! Pitch imbalance angle of Blades  [deg]

          PITCH_wt    (nbod) =-PITCH_wt    (nbod)
          PITCH_Imb_wt(nbod) =-PITCH_Imb_wt(nbod)
        else
!---------- Helicopter
          PITCH_Imb_wt(:) = 0.d0
        read (11, *)          &
          PITCH_wt    (nbod) ,&                                         ! Initial Constant Pitch angle of Blades  [deg]
          PITCHC_wt   (nbod) ,&                                         ! Initial Cosine   Pitch angle of Blades  [deg]
          PITCHS_wt   (nbod)                                            ! Initial   Sine   Pitch angle of Blades  [deg]
        endif

          PITCH_wt    (nbod) = PITCH_wt    (nbod) / R2D   
          PITCHC_wt   (nbod) = PITCHC_wt   (nbod) / R2D   
          PITCHS_wt   (nbod) = PITCHS_wt   (nbod) / R2D   
          PITCH_Imb_wt(nbod) = PITCH_Imb_wt(nbod) / R2D   
     enddo

   if (NBODBTWT_el >0) then
     read (11, *) OMEGA !; write(*,*) 'omega=',OMEGA                     ! Initial Rotor's Speed (LSS) [RPM]
     read (11, *) TREF  !; write(*,*) 'TREF =',TREF                      ! Initial HSS Torque [Nm]
   endif

   if (ICASE_el >= 3.or.IAPPL_el == 2) then

     Allocate ( floater(NFLOATER_el) )

     do i = 1, 6
        read(11,*)                                               &
         ( floater(nf)% UflInit        (i),                      &      ! Floater's Initial position (6 s q's) [m, deg]
           floater(nf)% IQfl_active_el (i), nf = 1, NFLOATER_el )       ! Switch floater's dofs [0:disabled, 1:enabled]
     enddo                                                              

     do nf = 1, NFLOATER_el
        floater(nf)% UflInit (4:6) = floater(nf)% UflInit (4:6) / R2D
     enddo

   endif

   OMEGAR = PI2                                                         ! [rpm] Used when no blades are present, just to set the period

   if (NBODBTWT_el >0) then
   read (11, *)!** WT Characteristics
     read (11, *) Htow                                                  ! Tower top    Height [m]
     read (11, *) Htow0                                                 ! Tower bottom Height [m]
     read (11, *) Hshof                                                 ! Tower off-set from top to nacelle [m]
     read (11, *) Hsh                                                   ! Shaft length [m]
     read (11, *) Hhub                                                  ! distance from hub to root [m]
     read (11, *) YAW                                                   ! YAW  of the nacelle [deg]
     read (11, *) TILT                                                  ! TILT of the nacelle [deg]
     read (11, *) CONE                                                  ! Blade Precone angle [deg]
     read (11, *) OMEGAR                                                ! Rotor's Speed (nominal) [RPM]
     read (11, *) RAT_GEAR                                              ! Gearbox Ratio
     read (11, *) IACT_YAW                                              ! yaw actuator
     read (11, *) CSTIF_YAW                                             ! stiffness of yaw connection [Nm/ rad   ]
     read (11, *) CDAMP_YAW                                             ! damping   of yaw connection [Nm/(rad/s)]
   endif

   read (11, *)!**  Damping coefficients
     if (NBLADE_el   >  0          ) read (11, *) COEFM_wt(1)           ! Blade  mass      Rayleigh damping coefficient
     if (NBLADE_el   >  0          ) read (11, *) COEFK_wt(1)           ! Blade  stiffness Rayleigh damping coefficient
     if (NBODBTWT_el >= NBLADE_el+1) read (11, *) COEFM_wt(NBLADE_el+1) ! Shaft  mass      Rayleigh damping coefficient
     if (NBODBTWT_el >= NBLADE_el+1) read (11, *) COEFK_wt(NBLADE_el+1) ! Shaft  stiffness Rayleigh damping coefficient
     if (NBODBTWT_el == NBLADE_el+2) read (11, *) COEFM_wt(NBLADE_el+2) ! Tower  mass      Rayleigh damping coefficient
     if (NBODBTWT_el == NBLADE_el+2) read (11, *) COEFK_wt(NBLADE_el+2) ! Tower  stiffness Rayleigh damping coefficient
  if ((ICASE_el==2).or.(ICASE_el==4))read (11, *) COEFM_wt(NBLADE_el+3) ! Jacket mass      Rayleigh damping coefficient
  if ((ICASE_el==2).or.(ICASE_el==4))read (11, *) COEFK_wt(NBLADE_el+3) ! Jacket stiffness Rayleigh damping coefficient
     if (NBLADE_el   >  0          ) read (11, *) CCRIT_wt(1)           ! Blade  modal damping ratio for defining structural damping (absolute value, NOT %)
     if (NBODBTWT_el >= NBLADE_el+1) read (11, *) CCRIT_wt(NBLADE_el+1) ! Shaft  modal damping ratio for defining structural damping (absolute value, NOT %)
     if (NBODBTWT_el == NBLADE_el+2) read (11, *) CCRIT_wt(NBLADE_el+2) ! Tower  modal damping ratio for defining structural damping (absolute value, NOT %)
  if ((ICASE_el==2).or.(ICASE_el==4))read (11, *) CCRIT_wt(NBLADE_el+3) ! Jacket modal damping ratio for defining structural damping (absolute value, NOT %)
                                     read (11, *) CCRIT_wt(0), FREQ0_el ! modal  damping ration for modes higher than FREQ0_el

   read (11, *)!**  Simulation parameters
     read (11, *) TIME_MAX                                              ! simulation total length [sec]
     read (11, *) DT_el                                                 ! time-step               [sec]
     read (11, *) BITA_el                                               ! bita coefficient for Newmark Integration
     read (11, *) GAMMA_el                                              ! gama coefficient for Newmark Integration
     read (11, *) IRERUN                                                ! IRERUN = 1: Continue the run reading the last backup
     read (11, *) NTIMEBACK                                             ! every NTIMEBACK timesteps the backup is updated
     read (11, *) MAXERR                                                ! maximum error convergance criterion (not mean error)
     read (11, *) ITERMAX                                               ! maximum number of Iterations per time-step


     COEFM_wt(2:NBLADE_el) = COEFM_wt(1);
     COEFK_wt(2:NBLADE_el) = COEFK_wt(1);
     CCRIT_wt(2:NBLADE_el) = CCRIT_wt(1);

                                         nbod = NBODBTWT_el
     if ((ICASE_el==2).or.(ICASE_el==4)) nbod = NBODBTWT_el+1

     write(10,*)
     write(10,*) 'Rayleigh M',(COEFM_wt(i),i=1,nbod)
     write(10,*) 'Rayleigh K',(COEFK_wt(i),i=1,nbod)
     write(10,*) 'CCrit     ',(CCRIT_wt(i),i=1,nbod)
     write(10,*)

     if ((IACT_YAW == 1).and.(NBODBTWT_el /= NBLADE_el+2)) then
        write(*,*) 'The yaw actuator is disabled, because no tower is present!'
        IACT_YAW = 0
     endif

     YAW        = YAW  / R2D
     TILT       = TILT / R2D
     CONE       = CONE / R2D
     OMEGA      = OMEGA  * RPM2RAD
     NTIMEM     = int(TIME_MAX/DT_el) + 1
     TPERIOD_el = 60.d0 / OMEGAR
     OMEGAR     = OMEGAR * RPM2RAD

     write(10,*)
     write(10,*) 'DT_el     ',DT_el
     write(10,*) 'TIME_MAX  ',TIME_MAX
     write(10,*) 'NTIMEM    ',NTIMEM
     write(10,*) 'TPERIOD_el',TPERIOD_el
     write(10,*)


   read (11,*)!**  Concentrated masses
     read (11,*) NCOMA_el
     read (11,*)!** nb    Hta      Xcm    Ycm    Zcm    Mass[kg]   IPx      IPy     IPz   Ixy   Ixz   Iyz  [local body c.s.]

     if ( NCOMA_el > 0 ) then

        Allocate (conc_mass(NCOMA_el))

        do ncal = 1, NCOMA_el
 
           read (11,*)                &         !-- all properties are defined wrt the local element c.s.
            conc_mass(ncal)%NBCOMA_el,&         !-- body at which mass is attached
            conc_mass(ncal)%Hta      ,&         !-- position along the beam axis at which the mass is attached
            conc_mass(ncal)%Xoff     ,&         !-- X local 
            conc_mass(ncal)%Yoff     ,&
            conc_mass(ncal)%Zoff     ,&         !-- Z local
            conc_mass(ncal)%Mass     ,&         !-- mass
            conc_mass(ncal)%IPx      ,&         !-- read IPx, wrt the x,z offsets --> will become Ixx
            conc_mass(ncal)%IPy      ,&         !-- read IPy, wrt the x,z offsets --> will become Iyy
            conc_mass(ncal)%IPz      ,&         !-- read IPz, wrt the x,z offsets --> will become Izz 
            conc_mass(ncal)%Ixy      ,&
            conc_mass(ncal)%Ixz      ,&
            conc_mass(ncal)%Iyz      
        enddo !ncal
     endif !NCOMA_el>0

   close (11)


 END Subroutine read_dfile
!----------------------------------------------------------------------
 Subroutine allocate_main_matrices ( COEFM_wt  , COEFK_wt , CCRIT_wt    ,& 
                                     PHI0_wt   , PITCH_wt , PITCH_Imb_wt,&
                                                 PITCHC_wt, PITCHS_wt   ,&
                                     NBODSUB_wt, IDOF_wt  , IBTYP_wt    ,&
                                     NTEB_wt   , IBT      , ISBT        , IWRITE )
!----------------------------------------------------------------------

 Use Cbeam
 Use Jacket

   implicit none

   integer :: IBT, ISBT
   integer :: NBODSUB_wt(IBT     ), IDOF_wt  (IBT     )
   integer :: IBTYP_wt  (IBT,ISBT), NTEB_wt  (IBT,ISBT)
   integer :: IWRITE    (IBT,2   )
   real(8) :: PHI0_wt   (IBT     ), PITCH_wt (IBT     ), PITCH_Imb_wt(  IBT)
   real(8) ::                       PITCHC_wt(IBT     ), PITCHS_wt   (  IBT)
   real(8) :: COEFM_wt  (IBT     ), COEFK_wt (IBT     ), CCRIT_wt    (0:IBT)
   integer :: nbod,nbsub,nnsub,nb_el,nn_el,nnod,nnodT,nel,NTNOD
   integer :: nq, nq01, nq02, nqtot,i, ndf,nod, nbc


   NSHAFT_el    = NBLADE_el + 1
   NTOWER_el    = NBLADE_el + 2
   NJACKET_el   = NBLADE_el + 3
   IREDUCT_el   = 1                                                                  ! 0:no Matrix reduction for all b.c and all dependent q's,     1:yes *** default 1 ***
   IREDUCTJa_el = 0                                                                  ! 0:no reduction for all jacket dofs from the system equations 1:yes *** default 1 ***[only if ICASE_el=2]
   NBuoyaTap_el = 0                                                                  ! by default if ICASE_el /= 2 or 4 --> no Buoyancy Tapers
   if (ICASE_el < 3) ITRUSS_el = 0                                                   ! only for the floating case the moorings are supported atm

!qq
   open(30,file='ireduct.inp')
     read(30,*,end=1,err=1) IREDUCT_el
     read(30,*,end=1,err=1) IREDUCTJa_el
   close(30)
   goto 2

1  continue
#  if   OS == 0
       call system ('rm -f ireduct.inp')
#  elif OS == 1
       call system ('del ireduct.inp')
#  endif

2  if (ICASE_el/=2) &
   IREDUCTJa_el = 0


   NBODBT_el    = NBODBTWT_el + NBODBTJA_el

!--- Allocate 1: NBODBT_el Vars
   Allocate ( ACCU%NQSACC_el  (0:NBODBTWT_el    ) )
   Allocate ( ACCU%NDFTBACC_el(0:NBODBT_el      ) )
   Allocate ( ACCU%NQSBACC_el (  NBODBTWT_el,0:NSUBperBM_el+1) )
   Allocate ( body            (  NBODBT_el      ) )
   Allocate ( COEFM_el        (  NBODBT_el      ) )
   Allocate ( COEFK_el        (  NBODBT_el      ) )
   Allocate ( CCRIT_el        (0:NBODBTWT_el+1  ) )
   Allocate ( PITCH0_el       (  NBLADE_el      ) )
   Allocate ( PITCH_Imb_el    (  NBLADE_el      ) )
   Allocate ( PITCHC_el       (  NBLADE_el      ) )
   Allocate ( PITCHS_el       (  NBLADE_el      ) )
   Allocate ( PHI0        (max(  NBLADE_el  ,1  ) ) );   PHI0      (:)   = 0.d0;
   Allocate ( IWRITE_el       (  NBODBTWT_el+1,2) )  ;   IWRITE_el (:,:) = 0;

      PHI0        (1:NBLADE_el        ) = PHI0_wt     (1:NBLADE_el        );
      PITCH0_el   (1:NBLADE_el        ) = PITCH_wt    (1:NBLADE_el        );
      PITCH_Imb_el(1:NBLADE_el        ) = PITCH_Imb_wt(1:NBLADE_el        );
      PITCHC_el   (1:NBLADE_el        ) = PITCHC_wt   (1:NBLADE_el        );
      PITCHS_el   (1:NBLADE_el        ) = PITCHS_wt   (1:NBLADE_el        );
      IWRITE_el   (1:NBODBTWT_el+1,1:2) = IWRITE      (1:NBODBTWT_el+1,1:2);

!------ Damping Coefficients
      COEFM_el(1:NBODBTWT_el) = COEFM_wt(1:NBODBTWT_el);
      COEFK_el(1:NBODBTWT_el) = COEFK_wt(1:NBODBTWT_el);
      CCRIT_el(0:NBODBTWT_el) = CCRIT_wt(0:NBODBTWT_el);

!--- JACKET
   if ((ICASE_el == 2).or.(ICASE_el == 4)) then
      COEFM_el(NBODBTWT_el+1:NBODBT_el) = COEFM_wt(NJACKET_el)
      COEFK_el(NBODBTWT_el+1:NBODBT_el) = COEFK_wt(NJACKET_el)
      CCRIT_el(NBODBTWT_el+1          ) = CCRIT_wt(NJACKET_el)
   endif


   do nbod = 1, NBODBTWT_el
      body(nbod)%IDOF_el    = IDOF_wt   (nbod)                                       ! 0: Elastic dofs disabled or 1: enabled
      body(nbod)%NBODSUB_el = NBODSUB_wt(nbod)                                       ! Sub-bodies per Body
   enddo
!--- JACKET
   do nbod = NBODBTWT_el+1, NBODBT_el
      body(nbod)%IDOF_el    = IDOF_wt   (NJACKET_el)
      body(nbod)%NBODSUB_el = 1
   enddo


!--- Allocate 2:
   ACCU%NQSACC_el(0) = 0
   nb_el             = 0
   nn_el             = 0

   do nbod = 1, NBODBTWT_el
      nbsub = body(nbod)%NBODSUB_el
      Allocate ( body(nbod)%NBODTGLB_el( nbsub   ) )
      Allocate ( body(nbod)%NNODTGLB_el( nbsub +1) )
      if (nbod<=NBLADE_el) then
         Allocate ( body(nbod)%PRECURV    (0:nbsub+1        ) )
                    body(nbod)%PRECURV    (0:nbsub+1        ) = 0.d0;
         Allocate ( body(nbod)%PRECURV_mat(0:nbsub+1,  3,  3) )
                    body(nbod)%PRECURV_mat(0:nbsub+1,1:3,1:3) = 0.d0;
         do i = 1, 3
                    body(nbod)%PRECURV_mat(0:nbsub+1,  i,  i) = 1.d0;
         enddo
      endif

      do nbsub = 1, body(nbod)%NBODSUB_el
         nb_el = nb_el + 1
         body(nbod)%NBODTGLB_el(nbsub) = nb_el                                       ! map subbodies
      enddo
 
      ACCU%NQSBACC_el(nbod,0)=0

      do nnsub=1,body(nbod)%NBODSUB_el+1                                             ! for all subbodies' nodes
         nn_el = nn_el + 1
         body(nbod)%NNODTGLB_el(nnsub) = nn_el                                       ! map nodes
         ACCU%NQSBACC_el  (nbod,nnsub) = ACCU%NQSBACC_el(nbod,nnsub-1) + 6           ! num ACC of Q of subbodies of each body
      enddo

      ACCU%NQSACC_el(nbod) = ACCU%NQSACC_el (nbod-1                          ) + &
                             ACCU%NQSBACC_el(nbod  , body(nbod)%NBODSUB_el+1 )       ! num ACC of Q of bodies
   enddo !nbod


   if (NBODBTWT_el > 0 ) then
      NBODTWT_el = body(NBODBTWT_el)%NBODTGLB_el(body(NBODBTWT_el)%NBODSUB_el  )     ! number of Subbodies Total
      NNODTWT_el = body(NBODBTWT_el)%NNODTGLB_el(body(NBODBTWT_el)%NBODSUB_el+1)     ! number of Nodes     Total
   else
      NBODTWT_el = 0
      NNODTWT_el = 0
   endif


!--- JACKET
   do nbod = NBODBTWT_el+1,NBODBT_el

       Allocate ( body(nbod)%NBODTGLB_el( body(nbod)%NBODSUB_el) )
       Allocate ( body(nbod)%NNODTGLB_el( body(nbod)%NBODSUB_el) )

       do nbsub = 1, body(nbod)%NBODSUB_el
          nb_el = nb_el + 1
          nn_el = nn_el + 1

          body(nbod)%NBODTGLB_el(nbsub) = nb_el                                      ! map subbodies
          body(nbod)%NNODTGLB_el(nbsub) = nn_el                                      ! map nodes
       enddo
 
   enddo

   if     (ICASE_el==2) then
      NBODT_el = body(NBODBT_el)%NBODTGLB_el(body(NBODBT_el)%NBODSUB_el)             ! number of Subbodies Total
      NNODT_el = body(NBODBT_el)%NNODTGLB_el(body(NBODBT_el)%NBODSUB_el)             ! number of Nodes     Total
   elseif (ICASE_el==4) then
      NBODT_el = body(NBODBT_el)%NBODTGLB_el(body(NBODBT_el)%NBODSUB_el)             ! number of Subbodies Total
      NNODT_el = body(NBODBT_el)%NNODTGLB_el(body(NBODBT_el)%NBODSUB_el)             ! number of Nodes     Total
   else
      NBODT_el = NBODTWT_el
      NNODT_el = NNODTWT_el
   endif


   write (10,*) 'ICASE_el     =', ICASE_el
   write (10,*) 'IEIG         =', IEIG
   write (10,*) 'IMODAL       =', IMODAL
   write (10,*) 'IREDUCT_el   =', IREDUCT_el
   write (10,*) 'IREDUCTJA_el =', IREDUCTJA_el
   write (10,*)
   write (10,*) 'NBODBTWT_el  =', NBODBTWT_el
   write (10,*) 'NBLADE_el    =', NBLADE_el
   write (10,*) 'NBODTWT_el   =', NBODTWT_el
   write (10,*) 'NNODTWT_el   =', NNODTWT_el
   write (10,*)
   write (10,*) 'NBODBT_el    =', NBODBT_el
   write (10,*) 'NBLADE_el    =', NBLADE_el
   write (10,*) 'NBODT_el     =', NBODT_el
   write (10,*) 'NNODT_el     =', NNODT_el
   write (10,*) 'NBODSUB_el   =', (body(i)%NBODSUB_el,i=1,NBODBTWT_el)
   write (10,*)


!--- Allocate 3:NBODT_el Vars
   Allocate ( ACCU%NDFTACC_el (0:NBODT_el) )
   Allocate ( subbody         (  NBODT_el) )

   if (NBODBTWT_el > 0) then
      nb_el = body(1)%NBODSUB_el * NBLADE_el
      Allocate ( UTels2aer  (nb_el, 3) )
      Allocate ( UT1els2aer (nb_el, 3) )
      Allocate ( UT2els2aer (nb_el, 3) )
   endif

   do nbod = 1, NBODBT_el
      do nbsub = 1,body(nbod)%NBODSUB_el
         nb_el = body(nbod)%NBODTGLB_el(nbsub)
         Allocate ( subbody(nb_el)%NDFPN_el    (  NNPE_el) )
         Allocate ( subbody(nb_el)%NDFPNACC_el (0:NNPE_el) )
      enddo
   enddo



!--- set subbody:: IBTYP_el,NTEB_el,IAENUM, NDFPN_el,NDFPNACC_el,NDFPE_el, NEQPE_el
!--- set    body:: NTEBB_el
   do nbod=1,NBODBT_el
    body(nbod)%NTEBB_el = 0
    do nbsub=1,body(nbod)%NBODSUB_el

      nb_el = body(nbod)%NBODTGLB_el(nbsub)

      if (nbod<=NBODBTWT_el) then
         subbody(nb_el)%IBTYP_el = IBTYP_wt(nbod,nbsub)                                 ! 0 no loads, 1 aero loads, 2 hydro loads
         subbody(nb_el)%NTEB_el  = NTEB_wt (nbod,nbsub)                                 ! Num of elements of subbody

         if (nbod <= NBLADE_el) subbody(nb_el)%IAENUM = nbod                            ! aerodynamic body correspondance
         if (nbod  > NBLADE_el) subbody(nb_el)%IAENUM = 0                               ! aerodynamic body correspondance
      else
         subbody(nb_el)%IBTYP_el = ITYPE_ja(nb_el-NBODTWT_el)
         if (IBTYP_wt(NJACKET_el,1)==0) &
         subbody(nb_el)%IBTYP_el = 0
         subbody(nb_el)%NTEB_el  = NTEB_ja (nb_el-NBODTWT_el)
         subbody(nb_el)%IAENUM   = 0                                                    ! aerodynamic body correspondance
      endif

      body(nbod)%NTEBB_el = body(nbod)%NTEBB_el+subbody(nb_el)%NTEB_el

      if     (NNPE_el==2) then
!--------- Number of dof per local node of each element - Each subb may have more than 1 element (same type)
         subbody(nb_el)%NDFPN_el   (1) = 6
         subbody(nb_el)%NDFPN_el   (2) = 6

!--------- Accumulative dof per local element node (for each element of the sub-body)
         subbody(nb_el)%NDFPNACC_el(0) = 0
         subbody(nb_el)%NDFPNACC_el(1) = 6
         subbody(nb_el)%NDFPNACC_el(2) = 12

!--------- Total num of dofs per element
         subbody(nb_el)%NDFPE_el       = 12
      elseif (NNPE_el==5) then
!--------- Number of dof per local node of each element - Each subb may have more than 1 element (same type)
         subbody(nb_el)%NDFPN_el   (1) = 6
         subbody(nb_el)%NDFPN_el   (2) = 1
         subbody(nb_el)%NDFPN_el   (3) = 1
         subbody(nb_el)%NDFPN_el   (4) = 1
         subbody(nb_el)%NDFPN_el   (5) = 6

!--------- Accumulative dof per local element node (for each element of the sub-body)
         subbody(nb_el)%NDFPNACC_el(0) = 0
         subbody(nb_el)%NDFPNACC_el(1) = 6
         subbody(nb_el)%NDFPNACC_el(2) = 7
         subbody(nb_el)%NDFPNACC_el(3) = 8
         subbody(nb_el)%NDFPNACC_el(4) = 9
         subbody(nb_el)%NDFPNACC_el(5) = 15

!--------- Total num of dofs per element
         subbody(nb_el)%NDFPE_el       = 15
      endif

!------ Number of equations per element
         subbody(nb_el)%NEQPE_el       = 6
    enddo!nbsub
   enddo!nbod


   if ((ICASE_el==2).or.(ICASE_el==4)) Deallocate ( NTEB_ja, ITYPE_ja )


!--- Set the total number of structural data given at nodes
   INODNELPROPWT_el = 0
   INODNELPROPbl_el = 0
   do nbod = 1, NBODBTWT_el
      if (nbod<=NBLADE_el) INODNELPROPbl_el = INODNELPROPbl_el + body(nbod)%NTEBB_el * 2
                           INODNELPROPWT_el = INODNELPROPWT_el + body(nbod)%NTEBB_el * 2
   enddo
                           INODNELPROP_el   = INODNELPROPWT_el + INODNELPROPJA_el



!--- Allocate 4:
   do    nbod  = 1, NBODBT_el
      do nbsub = 1, body(nbod)%NBODSUB_el
         nb_el = body(nbod)%NBODTGLB_el(nbsub)
         nel   = subbody(nb_el)%NTEB_el
         nod   = nel * (NNPE_el-1) + 1

         Allocate ( subbody(nb_el)%NODMB       (  nel, NNPE_el) )
!!       Allocate ( subbody(nb_el)%NODCB       (  nod         ) )
         Allocate ( subbody(nb_el)%NDFPNBACC_el(0:nod         ) )
         Allocate ( subbody(nb_el)%INODNEL_el  (nel, 2) )
         Allocate ( subbody(nb_el)%ALENG_el    (nel   ) )
         Allocate ( subbody(nb_el)%HTA_el      (nel +1) )
         Allocate ( subbody(nb_el)%PHIX_el     (nel   ) )
         Allocate ( subbody(nb_el)%PHIZ_el     (nel   ) )
         Allocate ( subbody(nb_el)%FCPL_el     (nel, 3) )
         Allocate ( subbody(nb_el)%AMPL_el     (nel, 3) )
!--------- Zero local external forces-moments
         subbody(nb_el)%FCPL_el (1:nel,1:3) = 0.d0;
         subbody(nb_el)%AMPL_el (1:nel,1:3) = 0.d0;
      enddo
   enddo
   do    nbod  = 1, NBLADE_el
      do nbsub = 1, body(nbod)%NBODSUB_el
         nb_el = body(nbod)%NBODTGLB_el(nbsub)
         nel   = subbody(nb_el)%NTEB_el
         nod   = nel * (NNPE_el-1) + 1
         Allocate ( subbody(nb_el)%Rinit (  nel+1,  6) )
                    subbody(nb_el)%Rinit (1:nel+1,1:6) = 0.d0;
      enddo
   enddo



!qq SET CORRECTLY the number of the properties
   Allocate ( beam_timo (INODNELPROP_el) )
   Allocate ( substr_mor(INODNELPROP_el) )

!-- initially set all inner radius zero
  substr_mor(:)%INDIAMET_el = 0.d0;


!--- Correspondance of local element node to global body nodes
!--- Accumulative dof per body node

   do nb_el = 1, NBODT_el  !for all subbodies including JACKET / Tripod

      NTNOD = 0 
      subbody(nb_el)%NDFPNBACC_el(0) = 0

      do nel = 1, subbody(nb_el)%NTEB_el !for all elements of each subbody
 
         do nnod  = 1, NNPE_el       !for all nodes of each element
            NTNOD = NTNOD + 1 
            subbody(nb_el)%NODMB       (nel,nnod)= NTNOD                                    ! map node of element to subbody counting
!!!         subbody(nb_el)%NODCB       (NTNOD)   = nnod                                     ! store element local number of node (1-5)
            subbody(nb_el)%NDFPNBACC_el(NTNOD)   = subbody(nb_el)%NDFPNBACC_el(NTNOD-1) + &
                                                   subbody(nb_el)%NDFPN_el    (nnod   )     ! number ACC of dof per node (of Subb)
         enddo

         NTNOD = NTNOD-1 !last node of an element is the same with the 1st node of the next element 
      enddo

      subbody(nb_el)%NNTB_el  = subbody(nb_el)%NODMB       (subbody(nb_el)%NTEB_el,NNPE_el) ! total num of element nodes (of Subb)
      subbody(nb_el)%NDFTB_el = subbody(nb_el)%NDFPNBACC_el(subbody(nb_el)%NNTB_el        ) ! total num of         dof   (of Subb)
   enddo !nb_el


   ACCU%NDFTACC_el(0) = 0
   do nb_el = 1, NBODT_el !for all subbodies
      ACCU%NDFTACC_el(nb_el) = ACCU%NDFTACC_el(nb_el-1) + subbody(nb_el)%NDFTB_el           ! num ACC of dof of all subbodies
      write (10,*)'NDFTACC_el=', nb_el, ACCU%NDFTACC_el(nb_el)
   enddo

   write (10,*)
   ACCU%NDFTBACC_el(0) = 0
   do nbod  = 1, NBODBT_el
      nbsub = body(nbod)%NBODSUB_el
      nb_el = body(nbod)%NBODTGLB_el(nbsub)

      ACCU%NDFTBACC_el(nbod) = ACCU%NDFTACC_el(nb_el)                                       ! num ACC of dof of all Bodies
      write (10,*)'NDFTBACC_el=', nbod, ACCU%NDFTBACC_el(nbod)!,NQSACC_el(nbod) 
   enddo




!--- add 6 Qs for the floating body 3 rotations and 3 translations
   if (ICASE_el >= 3.or.IAPPL_el == 2) then
       IQfl_tr  = 3
       IQfl_rot = 3
   else
       IQfl_tr  = 0
       IQfl_rot = 0
   endif

   IQfl = IQfl_rot + IQfl_tr

   if (NBODBTWT_el > 0 ) then
      NQSP   = IQfl*NFLOATER_el + 2                            ! Primary Qs (floating, Omega and Yaw Qs)
      NQSW   = IQfl*NFLOATER_el + 1                            ! Omega q
      NQSY   = IQfl*NFLOATER_el + 2                            ! Yaw q for free equation
   else
      NQSP   = 0
      NQSW   = 0
      NQSY   = 0
   endif


   NQS        = NQSP + ACCU%NQSACC_el (NBODBTWT_el)  ! total number of Q of all bodies (Primary + Total Q of all Subbodies)
   NDFBT_el   =        ACCU%NDFTACC_el(NBODT_el   )  ! total number of elastic dof
   NDFT_el    = NDFBT_el + NQS                       ! total Num of dof (elastic and Q)
   NDFBTWT_el =        ACCU%NDFTACC_el(NBODTWT_el )  ! total number of elastic dof of WT bodies


   write(10,*)
   write(10,*) 'NDFT_el', NDFT_el
   write(10,*) 'NQS    ', NQS
   write(10,*)



!--- Allocate main matrices
   Allocate ( AM_el    (NDFT_el,NDFT_el) )
   Allocate ( AC_el    (NDFT_el,NDFT_el) )
   Allocate ( AK_el    (NDFT_el,NDFT_el) )
   Allocate ( AQ_el    (NDFT_el        ) )
   Allocate ( INDSYSB  (NDFT_el        ) )
   Allocate ( INDSYSQ  (NQS            ) )
   Allocate ( UT_el    (NDFT_el), UT1_el (NDFT_el), UT2_el (NDFT_el) ) ! Current iteration deformation, velocity, acceleration
   Allocate ( UTP_el   (NDFT_el), UTP1_el(NDFT_el), UTP2_el(NDFT_el) ) ! Previous Time Step
   Allocate ( UT0_el   (NDFT_el), UT01_el(NDFT_el), UT02_el(NDFT_el) ) ! Perturbation
   Allocate ( CDAMPSTR (NDFT_el,NDFT_el) )


!--- For all q's define the exact matrices
!--- Allocate 8:
   Allocate ( QSnode(NNODT_el) )
!-------------------------------------------------------------------
! for each subbody nb_el
!
!   QSnode(nn_el)%Tot_el  (1) = number of q's from previous bodies
!                               which affect this sub-body.
!                               primary q's are included (no Rayleigh)
!   QSnode(nn_el)%Tot_el  (2) = number of q's total (qs of sub-bodies of
!                               previous bodies + qs of all previous 
!                               sub-bodies which belong to the same body
!                               and qs of current sub-body)
!-------------------------------------------------------------------

   do    nbod  = 1, NBODBTWT_el
      do nnsub = 1, body(nbod)%NBODSUB_el+1

         nn_el = body(nbod)%NNODTGLB_el(nnsub)
         nqtot = 0

         if     (nbod <= NBLADE_el    ) then
            nq = NQSP + ACCU%NQSACC_el(NBODBTWT_el) - ACCU%NQSACC_el(NBLADE_el  ) + ACCU%NQSBACC_el(nbod,nnsub)
         elseif (nbod == NBODBTWT_el  ) then
            nq = NQSP +                                                             ACCU%NQSBACC_el(nbod,nnsub)
         elseif (nbod == NBODBTWT_el-1) then
            nq = NQSP + ACCU%NQSACC_el(nbod+1)      - ACCU%NQSACC_el(nbod)        + ACCU%NQSBACC_el(nbod,nnsub)
         endif

         Allocate ( QSnode(nn_el)%IORD_el(nq) )

!NEW-START-------------------------------------------

!--------- Primary Q's (NQSP)
         do nq = 1, IQfl
            nqtot = nqtot + 1
            QSnode(nn_el)%IORD_el(nqtot) = nq
         enddo
         nq = NQSW
            nqtot = nqtot + 1
            QSnode(nn_el)%IORD_el(nqtot) = nq
         nq = NQSY
            nqtot = nqtot + 1
            QSnode(nn_el)%IORD_el(nqtot) = nq

!NEW-END--------------------------------------------- for the old remove !!
!--------- Primary Q's (NQSP)
   !!    do nq = 1, NQSP
   !!       nqtot = nqtot + 1
   !!       QSnode(nn_el)%IORD_el(nqtot) = nq
   !!    enddo

!--------- Previous Bodies' Qs
!        if (nbod <= NBLADE_el) then
!           nq01 = ACCU%NQSACC_el(NBLADE_el  ) + NQSP + 1
!           nq02 = ACCU%NQSACC_el(NBODBTWT_el) + NQSP
!           do nq = nq01, nq02
!              nqtot = nqtot + 1
!              QSnode(nn_el)%IORD_el(nqtot) = nq
!           enddo
!        elseif (nbod == NBODBTWT_el) then

!        elseif (nbod == NBODBTWT_el-1) then
!           nq01 = ACCU%NQSACC_el(nbod)   + NQSP + 1
!           nq02 = ACCU%NQSACC_el(nbod+1) + NQSP
!           do nq = nq01, nq02
!              nqtot = nqtot + 1
!              QSnode(nn_el)%IORD_el(nqtot) = nq
!           enddo
!        else
!           write(*,*)'fix Q'
!           stop
!        endif

!--------- Previous Bodies' Qs
         if ((nbod <= NSHAFT_el).and.(NBODBTWT_el>=NTOWER_el)) then !tower Qs
               nq01  = ACCU%NQSACC_el(NSHAFT_el  ) + NQSP + 1
               nq02  = ACCU%NQSACC_el(NBODBTWT_el) + NQSP
            do nq    = nq01, nq02
               nqtot = nqtot + 1
               QSnode(nn_el)%IORD_el(nqtot) = nq
            enddo
         endif
         if ((nbod <= NBLADE_el).and.(NBODBTWT_el>=NSHAFT_el)) then !shaft Qs
               nq01  = ACCU%NQSACC_el(NBLADE_el) + NQSP + 1
               nq02  = ACCU%NQSACC_el(NSHAFT_el) + NQSP
            do nq    = nq01, nq02
               nqtot = nqtot + 1
               QSnode(nn_el)%IORD_el(nqtot) = nq
            enddo
         endif

         QSnode(nn_el)%Tot_el(1) = nqtot

!--------- Current Body's Q's [with Rayleigh Damping]
            nq01  = ACCU%NQSACC_el(nbod-1) + NQSP + 1
            nq02  = ACCU%NQSACC_el(nbod-1) + NQSP + ACCU%NQSBACC_el(nbod,nnsub)
         do nq    = nq01, nq02
            nqtot = nqtot + 1
            QSnode(nn_el)%IORD_el(nqtot) = nq
         enddo
         QSnode(nn_el)%Tot_el(2) = nqtot

      enddo !nnsub
   enddo !nbod

   write(10,*         )
   write(10,*         ) '** QS correspondance **'
   write(10,'(1000i3)')0   ,0    , (nq,nq=1,NQS)

   do nbod  = 1,NBODBTWT_el
   do nnsub = 1,body(nbod)%NBODSUB_el+1
      nn_el = body(nbod)%NNODTGLB_el(nnsub)
      write(10,'(1000i3)')nbod,nnsub, (QSnode(nn_el)%IORD_el(nq),nq=1,QSnode(nn_el)%Tot_el(2))
   enddo !nnsub
   enddo !nbod

   write(10,*         )

!--- JACKET
!--- Jacket's bodies have no q's in case of bottom based jacket (ICASE_el=2)
!--- Jacket's bodies have  6 q's in case of floating     jacket (ICASE_el=4)

   if     (ICASE_el == 2) then
      do nn_el=NNODTWT_el+1,NNODT_el
         QSnode(nn_el)%Tot_el(1:2) = 0
      enddo
   elseif (ICASE_el == 4) then
      do nn_el=NNODTWT_el+1,NNODT_el

         Allocate ( QSnode(nn_el)%IORD_el(IQfl) )

         QSnode(nn_el)%Tot_el(1:2) = IQfl

         do nq=1,IQfl
            QSnode(nn_el)%IORD_el(nq) = nq
         enddo
      enddo
   endif !ICASE_el


!--- Allocate Transformation Matrices
   Allocate ( transf_mat (NNODT_el) )

   if (ICASE_el/=4) nnodT = NNODTWT_el
   if (ICASE_el==4) nnodT = NNODT_el

   do nn_el = 1, nnodT
      nq    = QSnode(nn_el)%Tot_el(2)
      Allocate ( transf_mat(nn_el)%A0_el    (nq,3,3) )
      Allocate ( transf_mat(nn_el)%AT0_el   (nq,3,3) )
      Allocate ( transf_mat(nn_el)%ATDA0_el (nq,3,3) )
      Allocate ( transf_mat(nn_el)%ATDA1_el (nq,3,3) )
      Allocate ( transf_mat(nn_el)%ATDDA0_el(nq,3,3) )
      Allocate ( transf_mat(nn_el)%ATDDA1_el(nq,3,3) )
      Allocate ( transf_mat(nn_el)%ATDDA2_el(nq,3,3) )
      Allocate ( transf_mat(nn_el)%R0_el    (nq,  3) )
      Allocate ( transf_mat(nn_el)%ATDDR0_el(nq,  3) )
      Allocate ( transf_mat(nn_el)%ATDDR1_el(nq,  3) )
      Allocate ( transf_mat(nn_el)%ATDDR2_el(nq,  3) )
      Allocate ( transf_mat(nn_el)%ATDR0_el (nq,  3) )
      Allocate ( transf_mat(nn_el)%ATDR1_el (nq,  3) )
   enddo
!--- only for floating jacket
   do nnod = NNODTWT_el+1, nnodT
      Allocate ( transf_mat(nnod)%Aja_el (3, 3) )
      Allocate ( transf_mat(nnod)%Rja_el (   3) )
   enddo

!--- only for blades --- (WT here)
   do nnod = 1, NNODTWT_el
      Allocate ( transf_mat(nnod)%ATPA_el  (3, 3) )
!!    Allocate ( transf_mat(nnod)%DATPA_el (3, 3) )
!!    Allocate ( transf_mat(nnod)%DDATPA_el(3, 3) )
   enddo

!--- no derived data type for these...
   Allocate ( ATP_el    ( NBLADE_el,     3, 3 ) )
!! Allocate ( DATP_el   ( NBLADE_el,     3, 3 ) )
!! Allocate ( DDATP_el  ( NBLADE_el,     3, 3 ) )
   Allocate ( ATPWr_el  ( NBLADE_el,     3, 3 ) )
   Allocate ( ATPWrG_el ( NBLADE_el,     3, 3 ) )
   Allocate ( AT0yaw_el (           NQS, 3, 3 ) )
!! Allocate ( A0yaw_el  (           NQS, 3, 3 ) )



!--- Allocate local element matrices
   do    nbod  = 1, NBODBT_el
      do nbsub = 1, body(nbod)%NBODSUB_el
         nb_el = body(nbod)%NBODTGLB_el(nbsub)
         nn_el = body(nbod)%NNODTGLB_el(nbsub)
         nel   = subbody(nb_el)%NTEB_el
         ndf   = subbody(nb_el)%NDFPE_el
         nq    = QSnode(nn_el)%Tot_el(2)

         if (nbod>NBODBTWT_el) nq=6

         Allocate ( subbody(nb_el)%AMLOC_el  ( nel, ndf, ndf ) )
         Allocate ( subbody(nb_el)%ACLOC_el  ( nel, ndf, ndf ) )
         Allocate ( subbody(nb_el)%AKLOC_el  ( nel, ndf, ndf ) )
         Allocate ( subbody(nb_el)%AFLOC_el  ( nel, ndf      ) )

         if ((ICASE_el/=4).and.(nbod>NBODBTWT_el)) cycle

         Allocate ( subbody(nb_el)%AMLOCNQ_el( nel, ndf, nq  ) )
         Allocate ( subbody(nb_el)%ACLOCNQ_el( nel, ndf, nq  ) )
         Allocate ( subbody(nb_el)%AKLOCNQ_el( nel, ndf, nq  ) )
      enddo
   enddo


!--- Set Gauss Integration parameters
   Allocate ( XGAUSS(IORDER),WGAUSS(IORDER) )

   if (IORDER==6) then
      XGAUSS(1) = -0.932469514203152 
      XGAUSS(2) = -0.661209386466265
      XGAUSS(3) = -0.238619186083197 
      XGAUSS(4) =  0.238619186083197  
      XGAUSS(5) =  0.661209386466265  
      XGAUSS(6) =  0.932469514203152  
      
      WGAUSS(1) =  0.171324492379170
      WGAUSS(2) =  0.360761573048139
      WGAUSS(3) =  0.467913934572691
      WGAUSS(4) =  0.467913934572691
      WGAUSS(5) =  0.360761573048139
      WGAUSS(6) =  0.171324492379170
   elseif (IORDER==8) then
      XGAUSS(1) = -0.9602898564975362317
      XGAUSS(2) = -0.7966664774136267396
      XGAUSS(3) = -0.5255324099163289858
      XGAUSS(4) = -0.18343464249564980494
      XGAUSS(5) =  0.18343464249564980494
      XGAUSS(6) =  0.5255324099163289858
      XGAUSS(7) =  0.7966664774136267396
      XGAUSS(8) =  0.9602898564975362317
      
      WGAUSS(1) =  0.10122853629037626
      WGAUSS(2) =  0.2223810344533745
      WGAUSS(3) =  0.3137066458778873
      WGAUSS(4) =  0.362683783378362
      WGAUSS(5) =  0.362683783378362
      WGAUSS(6) =  0.3137066458778873
      WGAUSS(7) =  0.2223810344533745
      WGAUSS(8) =  0.10122853629037626
   endif

!--- Allocate boundary conditions
   Allocate ( boundco (NBODT_el) )

   do nbod  = 1, NBODBTWT_el
   do nbsub = 1, body(nbod)%NBODSUB_el
      nb_el = body(nbod)%NBODTGLB_el(nbsub)

                        nbc = 1
      boundco (nb_el)%nbcpb = nbc

      Allocate ( boundco (nb_el)%nod   (nbc    ) )
      Allocate ( boundco (nb_el)%nel   (nbc    ) )
      Allocate ( boundco (nb_el)%nbcon (nbc    ) )
      Allocate ( boundco (nb_el)%indx  (nbc,6  ) )

                 boundco (nb_el)%nod   (nbc    ) = 1  ! b.c. at the 1st node
                 boundco (nb_el)%nel   (nbc    ) = 1  !   of the 1st element
                 boundco (nb_el)%nbcon (nbc    ) = 0  ! zero condition
                 boundco (nb_el)%indx  (nbc,1:6) = 1  ! Foundation will modify the index, depending on the dof
   enddo
   enddo


 END Subroutine allocate_main_matrices
!----------------------------------------------------------------------
!-- Open beam structural file
!----------------------------------------------------------------------
 Subroutine read_machi
!----------------------------------------------------------------------

 Use Cbeam
 Use Paths

   implicit none

   real(8) :: EIXXM,EIZZM,GAXM ,GAZM, MassINBAL,diamet(2)
   integer :: nbod,nbsub,nb_el,nel, NUM
   integer :: i, n(2), k, iread, m, j

   real(8) :: XL(3,2,300) !xyz,start-end,nel
   real(8) :: Ex(3),Ey(3),Ez(3),met, A(3,3), AT0(3,3), zer


   write(*,*)'Setting IAX_Ben_Swe to zero'
   IAX_Ben_Swe = 0                                       ! Axis for precurved blades [1:Bend, 3:Sweep]

   if (NBODBTWT_el==0) INODNELPROPWT_el = 0
   if (NBODBTWT_el==0) return

   write(10,*)'Open beam structural file'

   open (22, file=trim(file_machi), STATUS='UNKNOWN')
            read (22, *) !title

            i     = 0
   do       nbod  = 1, NBODBTWT_el
            read (22,*) !** body name
            read (22,*) MassINBAL, iread
            read (22,*) !**num

            if (MassINBAL<0.9d0.or.MassINBAL>1.1d0) then
               write(*,*) 'MassINBAL error input >10%',MassINBAL
               stop
            endif
            write(10,*)
            write(10,*)'nbod',nbod
            write(10,*)
      do    nbsub = 1,body(nbod)%NBODSUB_el
            nb_el =   body(nbod)%NBODTGLB_el(nbsub)
         do nel   = 1,subbody(nb_el)%NTEB_el   !all elements of each subbody 
            i     = i+1; n(1) = i;
            i     = i+1; n(2) = i;

!           write(10,*)'nb_el,nel,n1,n2',nb_el,nel,n1,n2
            subbody(nb_el)%INODNEL_el(nel,1:2) = n(1:2)


                  beam_timo(n(1))%K( :) = 0.d0;
                  beam_timo(n(2))%K( :) = 0.d0;
          if     (iread==0) then
            read (22,*)                                                                        &
             NUM, XL       (2,1,nel)      , XL       (1,1,nel)      , XL       (3,1,nel)     , &
                  beam_timo(n(1))%DENS_el , beam_timo(n(1))%XCMD_el , beam_timo(n(1))%ZCMD_el, &
                  beam_timo(n(1))%RIXX_el , beam_timo(n(1))%RIZZ_el , beam_timo(n(1))%RIXZ_el, &
                  beam_timo(n(1))%K( 7)   , beam_timo(n(1))%K( 9)   , beam_timo(n(1))%K(11)   ,& !EA_el    (1), EAX_el*  (1), EAZ_el  (1)
                  beam_timo(n(1))%K(16)   , beam_timo(n(1))%K(21)   , beam_timo(n(1))%K(18)   ,& !EIXX_el  (1), EIZZ_el  (1), EIXZ_el*(1)
                  beam_timo(n(1))%K(19)                                                       ,& !GIT_el   (1)
                  beam_timo(n(1))%K( 1)   , beam_timo(n(1))%K(12)                             ,& !GAX_el   (1), GAZ_el   (1)
                  beam_timo(n(1))%K( 5)   , beam_timo(n(1))%K(14)                             ,& !GAX_X_el (1), GAZ_Z_el*(1)
                  diamet     (1)                                                                 !  *minus is needed when applied to K
            read (22,*)                                                                        &
                  XL       (2,2,nel)      , XL       (1,2,nel)      , XL       (3,2,nel)     , &
                  beam_timo(n(2))%DENS_el , beam_timo(n(2))%XCMD_el , beam_timo(n(2))%ZCMD_el, &
                  beam_timo(n(2))%RIXX_el , beam_timo(n(2))%RIZZ_el , beam_timo(n(2))%RIXZ_el, &
                  beam_timo(n(2))%K( 7)   , beam_timo(n(2))%K( 9)   , beam_timo(n(2))%K(11)   ,& !EA_el    (2), EAX_el*  (2), EAZ_el  (2)
                  beam_timo(n(2))%K(16)   , beam_timo(n(2))%K(21)   , beam_timo(n(2))%K(18)   ,& !EIXX_el  (2), EIZZ_el  (2), EIXZ_el*(2)
                  beam_timo(n(2))%K(19)                                                       ,& !GIT_el   (2)
                  beam_timo(n(2))%K( 1)   , beam_timo(n(2))%K(12)                             ,& !GAX_el   (2), GAZ_el   (2)
                  beam_timo(n(2))%K( 5)   , beam_timo(n(2))%K(14)                             ,& !GAX_X_el (2), GAZ_Z_el*(2)
                  diamet     (2)                                                                 !  *minus is needed when applied to K
 
                  if (nbod<=NBLADE_el.and.dabs(Var_Paper(15))>1.d-5)&
                  beam_timo(n(1:2))%K(18) = beam_timo(n(1:2))%K(16)*Var_Paper(15) !EIxz=EIxx*Var_Paper(15)

                  if (nbod<=NBLADE_el.and.dabs(Var_Paper(16))>1.d-5)&
                  beam_timo(n(1:2))%K(17) = beam_timo(n(1:2))%K(16)*Var_Paper(16) !EIxy=EIxx*Var_Paper(16)
!qq
!                 if (nbod<=NBLADE_el) &
!                    EIXY_el(1:2) =-0.05d0*EIXX_el(1:2)
          elseif (iread==1) then
            read (22,*)                                                                        &
             NUM, XL       (2,1,nel)      , XL       (1,1,nel)      , XL       (3,1,nel)     , &
                  beam_timo(n(1))%DENS_el , beam_timo(n(1))%XCMD_el , beam_timo(n(1))%ZCMD_el, &
                  beam_timo(n(1))%RIXX_el , beam_timo(n(1))%RIZZ_el , beam_timo(n(1))%RIXZ_el, &
                  beam_timo(n(1))%K( 7)                                                       ,& !EA_el    (1)
                  beam_timo(n(1))%K( 1)   , beam_timo(n(1))%K(12)   , beam_timo(n(1))%K( 3)   ,& !GAX_el   (1), GAZ_el   (1), GAXZ_el (1)
                  beam_timo(n(1))%K(10)   , beam_timo(n(1))%K( 9)   , beam_timo(n(1))%K(11)   ,& !EAY_el   (1), EAX_el*  (1), EAZ_el  (1)
                  beam_timo(n(1))%K( 5)   , beam_timo(n(1))%K(14)                             ,& !GAX_X_el (1), GAZ_Z_el*(1)
                  beam_timo(n(1))%K(16)   , beam_timo(n(1))%K(21)   , beam_timo(n(1))%K(19)   ,& !EIXX_el  (1), EIZZ_el  (1), GIT_el  (1)
                  beam_timo(n(1))%K(18)   , beam_timo(n(1))%K(17)   , beam_timo(n(1))%K(20)   ,& !EIXZ_el* (1), EIXY_el* (1), EIYZ_el*(1)
                  diamet     (1)                                                                 !  *minus is needed when applied to K
            read (22,*)                                                                        &
                  XL       (2,2,nel)      , XL       (1,2,nel)      , XL       (3,2,nel)     , &
                  beam_timo(n(2))%DENS_el , beam_timo(n(2))%XCMD_el , beam_timo(n(2))%ZCMD_el, &
                  beam_timo(n(2))%RIXX_el , beam_timo(n(2))%RIZZ_el , beam_timo(n(2))%RIXZ_el, &
                  beam_timo(n(2))%K( 7)                                                       ,& !EA_el    (2)
                  beam_timo(n(2))%K( 1)   , beam_timo(n(2))%K(12)   , beam_timo(n(2))%K( 3)   ,& !GAX_el   (2), GAZ_el   (2), GAXZ_el (2)
                  beam_timo(n(2))%K(10)   , beam_timo(n(2))%K( 9)   , beam_timo(n(2))%K(11)   ,& !EAY_el   (2), EAX_el*  (2), EAZ_el  (2)
                  beam_timo(n(2))%K( 5)   , beam_timo(n(2))%K(14)                             ,& !GAX_X_el (2), GAZ_Z_el*(2)
                  beam_timo(n(2))%K(16)   , beam_timo(n(2))%K(21)   , beam_timo(n(2))%K(19)   ,& !EIXX_el  (2), EIZZ_el  (2), GIT_el  (2)
                  beam_timo(n(2))%K(18)   , beam_timo(n(2))%K(17)   , beam_timo(n(2))%K(20)   ,& !EIXZ_el* (2), EIXY_el* (2), EIYZ_el*(2)
                  diamet     (2)                                                                 !  *minus is needed when applied to K
!qq
!                 if (nbod<=NBLADE_el) &
!                    EIXY_el(1:2) =-0.05d0*EIXX_el(1:2)
          elseif (iread==2) then
            read (22,*)                                                                        &
             NUM, XL       (2,1,nel)      , XL       (1,1,nel)      , XL       (3,1,nel)     , &
                  beam_timo(n(1))%DENS_el , beam_timo(n(1))%XCMD_el , beam_timo(n(1))%ZCMD_el, &
                  beam_timo(n(1))%RIXX_el , beam_timo(n(1))%RIZZ_el , beam_timo(n(1))%RIXZ_el, &
                 (beam_timo(n(1))%K(j),j=1,21)                                               , &
                  diamet     (1)
            read (22,*)                                                                        &
                  XL       (2,2,nel)      , XL       (1,2,nel)      , XL       (3,2,nel)     , &
                  beam_timo(n(2))%DENS_el , beam_timo(n(2))%XCMD_el , beam_timo(n(2))%ZCMD_el, &
                  beam_timo(n(2))%RIXX_el , beam_timo(n(2))%RIZZ_el , beam_timo(n(2))%RIXZ_el, &
                 (beam_timo(n(2))%K(j),j=1,21)                                               , &
                  diamet     (2)
!qq
!!             do m  = 1, 2
!!                beam_timo(n(m))%K( 2)=0.d0 !K12
!!                beam_timo(n(m))%K( 4)=0.d0 !K14
!!                beam_timo(n(m))%K( 6)=0.d0 !K16
!!                beam_timo(n(m))%K( 8)=0.d0 !K23
!!                beam_timo(n(m))%K(13)=0.d0 !K34
!!                beam_timo(n(m))%K(15)=0.d0 !K36
!!                beam_timo(n(m))%K( 3)=0.d0 !K13-GxzA
!!                beam_timo(n(m))%K(10)=0.d0 !K25-EAy
!!                beam_timo(n(m))%K(17)=0.d0 !K45-EIxy
!!                beam_timo(n(m))%K(20)=0.d0 !K56-EIyz
!!             enddo
          else
            write(*,*)'error in machi-read iread';stop
          endif

          if (iread<2) then
             do m  = 1, 2
                  beam_timo(n(m))%K( 9) = -beam_timo(n(m))%K( 9) !EAX
                  beam_timo(n(m))%K(14) = -beam_timo(n(m))%K(14) !GAZ_Z
                  beam_timo(n(m))%K(17) = -beam_timo(n(m))%K(17) !EIXY
                  beam_timo(n(m))%K(18) = -beam_timo(n(m))%K(18) !EIXZ
                  beam_timo(n(m))%K(20) = -beam_timo(n(m))%K(20) !EIYZ
             enddo
          endif
                  beam_timo(n(1))%XL_el(1:3) = XL(1:3,1,nel)
                  beam_timo(n(2))%XL_el(1:3) = XL(1:3,2,nel)

!--------------------------------------- paper enable-disable Xcm, Xel offsets
                                       if (int(Var_Paper(6))==0) then
                                        do m  = 1, 2
                                           beam_timo(n(m))%XCMD_el = 0.d0
                                           beam_timo(n(m))%ZCMD_el = 0.d0
                                        enddo !m
                                       endif
                                       if (int(Var_Paper(7))==0) then
                                        do m  = 1, 2
                                           beam_timo(n(m))%K( 9)   = 0.d0 !-EAX !(2,4)
                                           beam_timo(n(m))%K(11)   = 0.d0 ! EAZ !(2,6)
                                        enddo !m
                                       endif
                                       
!---------------------------------------- paper scale Xcm, Xel offsets
                                        do m  = 1, 2
                                           beam_timo(n(m))%K( 9)   = beam_timo(n(m))%K( 9)   * Var_Paper(11) !-EAX !(2,4)
                                           beam_timo(n(m))%K(11)   = beam_timo(n(m))%K(11)   * Var_Paper(12) ! EAZ !(2,6)
                                           beam_timo(n(m))%XCMD_el = beam_timo(n(m))%XCMD_el * Var_Paper(13)
                                           beam_timo(n(m))%ZCMD_el = beam_timo(n(m))%ZCMD_el * Var_Paper(14)
                                        enddo

!------------ Stiff Dofs Enabled
            if ( (body(nbod)%IDOF_el == 0 ).or. &
                 (body(nbod)%IDOF_el == 2 ).and.(nel<=6) ) then !for stiff monopile
               do m  = 1, 2
                  beam_timo(n(m))%K( :) = 0.d0;
                  beam_timo(n(m))%K( 1) = BIG   !GAX  !(1,1)
                  beam_timo(n(m))%K( 7) = BIG   !EA   !(2,2)
                  beam_timo(n(m))%K(12) = BIG   !GAZ  !(3,3)
                  beam_timo(n(m))%K(16) = BIG   !EIXX !(4,4)
                  beam_timo(n(m))%K(19) = BIG   !GIT  !(5,5)
                  beam_timo(n(m))%K(21) = BIG   !EIZZ !(6,6)
               enddo !m
            endif

!------------ Apply the mass imbalance
            do m  = 1, 2
               beam_timo(n(m))%DENS_el  = beam_timo(n(m))%DENS_el * MassINBAL
               beam_timo(n(m))%RIXX_el  = beam_timo(n(m))%RIXX_el * MassINBAL
               beam_timo(n(m))%RIZZ_el  = beam_timo(n(m))%RIZZ_el * MassINBAL
               beam_timo(n(m))%RIXZ_el  = beam_timo(n(m))%RIXZ_el * MassINBAL
               
               beam_timo(n(m))%POLI_el  = beam_timo(n(m))%RIXX_el + beam_timo(n(m))%RIZZ_el
               beam_timo(n(m))%AMOMX_el = beam_timo(n(m))%DENS_el * beam_timo(n(m))%ZCMD_el
               beam_timo(n(m))%AMOMZ_el = beam_timo(n(m))%DENS_el * beam_timo(n(m))%XCMD_el
            enddo !m

            if (nbod>NSHAFT_el) &
               substr_mor(n(1:2))%DIAMET_el = diamet(1:2)
         enddo  ! nel


!------------ Define subbody length (ALENGB_el)
            nel   = subbody(nb_el)%NTEB_el
            subbody(nb_el)%ALENGB_el = dsqrt( (XL(1,2,nel)-XL(1,1,1))**2 + (XL(2,2,nel)-XL(2,1,1))**2 + (XL(3,2,nel)-XL(3,1,1))**2 )  !subbody's length


!------------ Define the Precurved rotation matrix(PRECURV_mat) for blades only
         if (nbod<=NBLADE_el) then
            if (nbsub==1)  call DIAGO (3,AT0)
            if (nbsub> 1)  AT0(1:3,1:3) = matmul( transpose(body(nbod)%PRECURV_mat(nbsub-1,1:3,1:3)), AT0(1:3,1:3) )

!------------ Set Ey
            Ey(1:3)  = XL(1:3,2,nel)-XL(1:3,1,1)
            met      = dsqrt(dot_product (Ey,Ey))
            if (met<1.d-10) then
               write(*,*) 'error in node definition of blade, sbod',nbod,nbsub; stop
            endif
            Ey(1:3)  = Ey(1:3)/met

!------------ Set Ez=i x Ey
            Ex(1)    = 1.d0;  Ex(2:3) = 0.d0;
            call EXTEPR ( Ex, Ey, Ez )
            met      = dsqrt(dot_product (Ez,Ez))
            if (met<1.d-15) then
               write(*,*) 'magnitude is zero in blade, sbod',nbod,nbsub; stop
            endif
            Ez(1:3)  = Ez(1:3)/met

!------------ Set Ex=Ey x Ez
            call EXTEPR ( Ey, Ez, Ex )
            met      = dsqrt(dot_product (Ex,Ex))
            Ex(1:3)  = Ex(1:3)/met

            A(1:3,1) = Ex(1:3)
            A(1:3,2) = Ey(1:3)
            A(1:3,3) = Ez(1:3)

            body(nbod)%PRECURV_mat(nbsub,1:3,1:3) = matmul(AT0(1:3,1:3),A(1:3,1:3))
         else
            do k = 1, 3, 2
             if (XL(k,1,nel) /= 0.d0 .or. XL(k,2,nel) /= 0.d0) then !X,Z
              write(*,*) 'at the moment no precurved geometry is supported for the other bodies except the blades',nbod
              stop
             endif
            enddo
         endif !nbod<=BLADE_el

!------------ Define elements length (ALENG_el), and hta coordinate along the elastic axis (HTA_el)
         do nel   = 1,subbody(nb_el)%NTEB_el   !all elements of each subbody 
            n(1:2)=   subbody(nb_el)%INODNEL_el(nel,1:2)
          !!subbody(nb_el)%HTA_el   (nel  ) = subbody(nb_el)%ALENGB_el *                                                                                      &
          !!                                  dsqrt(dot_product(XL(1:3,1,               nel    )-XL(1:3,1,1),XL(1:3,1,               nel    )-XL(1:3,1,1))) / &
          !!                                  dsqrt(dot_product(XL(1:3,2,subbody(nb_el)%NTEB_el)-XL(1:3,1,1),XL(1:3,2,subbody(nb_el)%NTEB_el)-XL(1:3,1,1)))
          !!subbody(nb_el)%HTA_el   (nel+1) = subbody(nb_el)%ALENGB_el *                                                                                      &
          !!                                  dsqrt(dot_product(XL(1:3,2,               nel    )-XL(1:3,1,1),XL(1:3,2,               nel    )-XL(1:3,1,1))) / &
          !!                                  dsqrt(dot_product(XL(1:3,2,subbody(nb_el)%NTEB_el)-XL(1:3,1,1),XL(1:3,2,subbody(nb_el)%NTEB_el)-XL(1:3,1,1)))
            subbody(nb_el)%HTA_el   (nel  ) = subbody(nb_el)%ALENGB_el * (XL(2,1,nel)-XL(2,1,1)) / (XL(2,2,subbody(nb_el)%NTEB_el)-XL(2,1,1))
            subbody(nb_el)%HTA_el   (nel+1) = subbody(nb_el)%ALENGB_el * (XL(2,2,nel)-XL(2,1,1)) / (XL(2,2,subbody(nb_el)%NTEB_el)-XL(2,1,1))
            subbody(nb_el)%ALENG_el (nel  ) = subbody(nb_el)%HTA_el(nel+1) - subbody(nb_el)%HTA_el  (nel)  !element's length

!--------------- Perform checks of input data [should be positive]
               zer = 1.d-10
               if (subbody  (nb_el)%ALENG_el(nel) < zer                            ) then; write(*,'(a,3i4)') 'check HTA  in machi',nbod,nbsub,nel;stop;endif;
               if (NTIMEM>1.or.IEIG>0.or.IMODAL>0) then !otherwise the mass could be zero for static calculations in order to avoid the gravity load
               if (beam_timo(n(1))%DENS_el < zer .or. beam_timo(n(2))%DENS_el < zer) then; write(*,'(a,3i4)') 'check DENS in machi',nbod,nbsub,nel;stop;endif;
               if (beam_timo(n(1))%RIXX_el < zer .or. beam_timo(n(2))%RIXX_el < zer) then; write(*,'(a,3i4)') 'check RIXX in machi',nbod,nbsub,nel;stop;endif;
               if (beam_timo(n(1))%RIZZ_el < zer .or. beam_timo(n(2))%RIZZ_el < zer) then; write(*,'(a,3i4)') 'check RIZZ in machi',nbod,nbsub,nel;stop;endif;
               endif
               if (beam_timo(n(1))%K( 1)   < zer .or. beam_timo(n(2))%K( 1)   < zer) then; write(*,'(a,3i4)') 'check GAX  in machi',nbod,nbsub,nel;stop;endif;
               if (beam_timo(n(1))%K( 7)   < zer .or. beam_timo(n(2))%K( 7)   < zer) then; write(*,'(a,3i4)') 'check EA   in machi',nbod,nbsub,nel;stop;endif;
               if (beam_timo(n(1))%K(12)   < zer .or. beam_timo(n(2))%K(12)   < zer) then; write(*,'(a,3i4)') 'check GAZ  in machi',nbod,nbsub,nel;stop;endif;
               if (beam_timo(n(1))%K(16)   < zer .or. beam_timo(n(2))%K(16)   < zer) then; write(*,'(a,3i4)') 'check EIXX in machi',nbod,nbsub,nel;stop;endif;
               if (beam_timo(n(1))%K(19)   < zer .or. beam_timo(n(2))%K(19)   < zer) then; write(*,'(a,3i4)') 'check GIT  in machi',nbod,nbsub,nel;stop;endif;
               if (beam_timo(n(1))%K(21)   < zer .or. beam_timo(n(2))%K(21)   < zer) then; write(*,'(a,3i4)') 'check EIZZ in machi',nbod,nbsub,nel;stop;endif;

            GAXM  = beam_timo(n(1))%K( 1) + beam_timo(n(2))%K( 1) !GAX  !(1,1)
            GAZM  = beam_timo(n(1))%K(12) + beam_timo(n(2))%K(12) !GAZ  !(3,3)
            EIXXM = beam_timo(n(1))%K(16) + beam_timo(n(2))%K(16) !EIXX !(4,4)
            EIZZM = beam_timo(n(1))%K(21) + beam_timo(n(2))%K(21) !EIZZ !(6,6) 

            subbody(nb_el)%PHIX_el(nel) = 12.d0*EIZZM/(GAXM*subbody(nb_el)%ALENG_el(nel)**2)
            subbody(nb_el)%PHIZ_el(nel) = 12.d0*EIXXM/(GAZM*subbody(nb_el)%ALENG_el(nel)**2)

            write(10,'(a,3i4,2f12.4,2i4)')'nb,nsb,nel,hta1,2, n1, n2',nbod,nbsub,nel,subbody(nb_el)%HTA_el(nel), subbody(nb_el)%HTA_el(nel+1), n(1), n(2)
         enddo  !nel
            write(10,'(a,2i4,1f12.4)')'nb,nsb,ALENGB_el',nbod,nbsub,subbody(nb_el)%ALENGB_el
            write(10,*)
      enddo  !nbsub
      if (nbod<=NBLADE_el) then
         body(nbod)%PRECURV_mat(body(nbod)%NBODSUB_el+1,1:3,1:3) = 0.d0 !unused
      endif
   enddo  !nbod 

   INODNELPROPWT_el = i

   close (22)


 END Subroutine read_machi
!----------------------------------------------------------------------
!
!-- Manage Concentrated Masses
!
!----------------------------------------------------------------------
!
!--- 1. Calculate Ixx, Iyy, Izz, based on the given polar inertias IPx, IPy, IPz
!--- 2. Modify Hta (local position along beam axis) if the mass is not attached
!       at the 1st subbody of the body.
!--- 3. Find the corresponing sub-body and element, based on the given distance
!---    from the begining of the body.
!    
!    Inputs are:
!                 1. the polar inertias IPx, IPy, IPz
!                 2. the distance from the begining of the body
!
! IPx = S dens (y**2+z**2) dA       Ixx = S dens x**2 dA
! IPy = S dens (x**2+z**2) dA       Iyy = S dens y**2 dA
! IPz = S dens (x**2+y**2) dA       Izz = S dens z**2 dA
!
! Iyy + Izz = IPx            Ixx = (IPy+IPz-IPx)/2
! Ixx + Izz = IPy   ====>    Iyy = (IPx+IPz-IPy)/2
! Ixx + Iyy = IPz            Izz = (IPx+IPy-IPz)/2
!
! Attention: Inside code all other FEM vars RIXX assume S dens z**2 dA
! only                                conc_mass()%Ixx = S dens x**2 dA
!------------------------------------------------------------------------
 Subroutine manage_con_mass
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   real(8) :: length, length0
   integer :: nbod,nbsub,nb_el,nel,ncal


   write(10,*)
   write(10,*)'-- Concentrated Masses --', NCOMA_el


!--- Modify hta, if the mass is not attached at the 1st sub-body.
!--- Find the corresponing sub-body and element, based on the given distance
!--- from the begining of the body.
   do    ncal    = 1, NCOMA_el
         nbod    = conc_mass(ncal)%NBCOMA_el
         length  = 0.d0

      do nbsub   = 1, body(nbod)%NBODSUB_el
         nb_el   = body(nbod)%NBODTGLB_el(nbsub)
         length0 = length
         length  = length + subbody(nb_el)%ALENGB_el

         if ( (length0 <= conc_mass(ncal)%Hta).and.(length >= conc_mass(ncal)%Hta) )then

            conc_mass(ncal)%NSBCOMA_el = nbsub
            length                     = length0

            do nel     = 1, subbody(nb_el)%NTEB_el
               length0 = length
               length  = length + subbody(nb_el)%ALENG_el(nel)

               if ( (length0 <= conc_mass(ncal)%Hta).and.(length >= conc_mass(ncal)%Hta) )then

                  conc_mass(ncal)%NELCOMA_el = nel
                  conc_mass(ncal)%Hta        = conc_mass(ncal)%Hta-length0

                  goto 11
               endif
            enddo !nel
         endif
      enddo !nbsub

      write(* ,'(a,i3,3f15.3)') 'error concentrated mass Hta is out of the last node',ncal,length0,length,conc_mass(ncal)%Hta
      write(10,'(a,i3,3f15.3)') 'error concentrated mass Hta is out of the last node',ncal,length0,length,conc_mass(ncal)%Hta
      stop

 11   continue
!     write(* ,'(a,3f15.5)') 'nconmass',conc_mass(ncal)%IPx,conc_mass(ncal)%IPy,conc_mass(ncal)%IPz
!------ Set Ixx, Iyy, Izz using input values for IPx, IPy, IPz
      if     ( (conc_mass(ncal)%IPy /= 0.d0).and.&
               (conc_mass(ncal)%IPx == 0.d0).and.&
               (conc_mass(ncal)%IPz == 0.d0)       )then

         conc_mass(ncal)%Ixx = conc_mass(ncal)%IPy/2.d0
         conc_mass(ncal)%Iyy = 0.d0
         conc_mass(ncal)%Izz = conc_mass(ncal)%IPy/2.d0
      elseif ( (conc_mass(ncal)%IPz /= 0.d0).and.&
               (conc_mass(ncal)%IPx == 0.d0).and.&
               (conc_mass(ncal)%IPy == 0.d0)       )then

         conc_mass(ncal)%Ixx = conc_mass(ncal)%IPz/2.d0
         conc_mass(ncal)%Iyy = conc_mass(ncal)%IPz/2.d0
         conc_mass(ncal)%Izz = 0.d0
      else
         conc_mass(ncal)%Ixx = (-conc_mass(ncal)%IPx + conc_mass(ncal)%IPy + conc_mass(ncal)%IPz )/2d0
         conc_mass(ncal)%Iyy = ( conc_mass(ncal)%IPx - conc_mass(ncal)%IPy + conc_mass(ncal)%IPz )/2d0
         conc_mass(ncal)%Izz = ( conc_mass(ncal)%IPx + conc_mass(ncal)%IPy - conc_mass(ncal)%IPz )/2d0
      endif

      if ( (conc_mass(ncal)%Ixx<0.d0).or.&
           (conc_mass(ncal)%Iyy<0.d0).or.&
           (conc_mass(ncal)%Izz<0.d0)      ) then
         write(* ,'(a,3f15.5)') 'error concentrated inertia is negative',conc_mass(ncal)%Ixx,conc_mass(ncal)%Iyy,conc_mass(ncal)%Izz
         write(10,'(a,3f15.5)') 'error concentrated inertia is negative',conc_mass(ncal)%Ixx,conc_mass(ncal)%Iyy,conc_mass(ncal)%Izz
         stop
      endif


!------ Output concentrated masses
      write(10,100)               &
       ncal                      ,&
       conc_mass(ncal)%NBCOMA_el ,&
       conc_mass(ncal)%NSBCOMA_el,&
       conc_mass(ncal)%NELCOMA_el,&
       conc_mass(ncal)%Hta       ,&
       conc_mass(ncal)%Xoff      ,&
       conc_mass(ncal)%Yoff      ,&
       conc_mass(ncal)%Zoff      ,&
       conc_mass(ncal)%Mass      ,&
       conc_mass(ncal)%IPx       ,&
       conc_mass(ncal)%IPy       ,&
       conc_mass(ncal)%IPz       ,&
       conc_mass(ncal)%Ixx       ,&
       conc_mass(ncal)%Iyy       ,&
       conc_mass(ncal)%Izz       ,&
       conc_mass(ncal)%Ixy       ,&
       conc_mass(ncal)%Ixz       ,&
       conc_mass(ncal)%Iyz      
   enddo !ncal

 100 format(4i5,15e16.6)


 END Subroutine manage_con_mass
!----------------------------------------------------------------------
!
!-- Calculate Htow, Hsh
!
!   Htow     : for Tower Shadow **read
!   Hsh      : for ROTMAT       **read
!   Hshof    : for ROTMAT       **read
!   Htow0    : for ROTMAT       **read
!
!----------------------------------------------------------------------
 Subroutine calc_htow_hsh
!----------------------------------------------------------------------

 Use Cbeam

   implicit none


   if     (NBODBTWT_el == NTOWER_el) then !tow,sh,blades
      Hsh      = 0.d0
   elseif (NBODBTWT_el == NSHAFT_el) then !no tow, only sh and blades
      Hsh      = 0.d0
      Hshof    = Htow  + Hshof - Htow0
   elseif ((NBODBTWT_el == NBLADE_el).and.(NBODBTWT_el>0)) then !only blades
      Hshof    = Htow  + Hshof - Htow0
   else
       write(* ,*)  'attension no Wind Turbine'
       write(10,*)  'attension no Wind Turbine'
       return
   endif

   write(10,*)'Htow0,Htow,Hshof' , Htow0, Htow, Hshof
   write(10,*)'Hsh,Hhub'         , Hsh  , Hhub


 END Subroutine calc_htow_hsh
!-----------------------------------------------------------
! --- Subroutine : FINIT_el
!----------------------------------------------------------
 Subroutine FINIT_el

 Use Cbeam
 Use Hydro

   implicit none


   Deallocate ( body       )
   Deallocate ( subbody    )
   Deallocate ( QSnode     )
   Deallocate ( beam_timo  )
   Deallocate ( substr_mor )

   if ( NCOMA_el >  0) Deallocate ( conc_mass  )

   Deallocate ( AM_el , AC_el  , AK_el  , AQ_el   , &
                UT_el , UTP_el , UT0_el , INDSYSB , &
                UT1_el, UTP1_el, UT01_el, INDSYSQ , &
                UT2_el, UTP2_el, UT02_el, CDAMPSTR    )

   Deallocate ( transf_mat  )
   Deallocate ( AT0yaw_el, ATP_el ) !!, DATP_el, DDATP_el, A0yaw_el
   Deallocate ( COEFM_el, COEFK_el, CCRIT_el, PITCH0_el, PITCH_Imb_el, PHI0, IWRITE_el )
   Deallocate ( XGAUSS, WGAUSS )


   call Dealloc_Rot

   if ((ICASE_el == 2).or.(ICASE_el == 4)) &
   Deallocate ( boundco )

   if (NBuoyaTap_el>0) &
   Deallocate ( buoyancy_tap )

   if ((ICASE_el == 2).or.(ICASE_el == 4)) &
   Deallocate ( IFloodbod_el, RES_hyd )

!--- Deallocate Hydro variables
   if (ICASE_el == 3) &
   call floater_simulator( 2, 0 )

!  Deallocate ( UTSR,VTSR,WTSR,TREC,RINW,TINW )

   call Runinfo (2)

   close(10) !SCR.TOT


 END Subroutine FINIT_el
!-----------------------------------------------------------
! --- Subroutine : DRAG_INIT_el
!----------------------------------------------------------
 Subroutine DRAG_INIT_el

 Use Cbeam
 Use Paths

   implicit none

   integer :: iang


   IDRAG_TOW = 0                  ! Index for activation (1) or not (0) the tower drag
   IDRAG_NAC = 0                  ! Index for activation (1) or not (0) the nacelle drag
   if ((NBODBTWT_el /= NTOWER_el).or.(IAERO_el == 0)) return

   open(16,file=trim(file_drag))  !'drag.inp'

   read(16,*) !** Tower Drag
   read(16,*) CDRAG_TOW_el        !  0.6d0
   read(16,*) Z_DRAG_LIM          !m 5.0d0
   if (dabs(CDRAG_TOW_el)>1.d-05) IDRAG_TOW=1


   read(16,*) !** nacelle Drag
   read(16,*) NAC_LX_el           !m 16.0d0 !10.025d0
   read(16,*) NAC_LY_el           !m  6.0d0 ! 3.400d0
   read(16,*) NAC_LZ_el           !m  6.0d0 ! 3.771d0
   read(16,*) NAC_RX_el           !m  0.d0
   read(16,*) NAC_RY_el           !m  0.d0
   read(16,*) NAC_RZ_el           !m  0.d0
   read(16,*) NANG_NAC
  do iang=1,NANG_NAC
   read(16,*) Alpha_NAC(iang), CL_NAC(iang), CD_NAC(iang)
   if (dabs(CL_NAC(iang))>1.d-05.or.dabs(CD_NAC(iang))>1.d-05) IDRAG_NAC=1 !activation of nacelle drag
  enddo !iang

   close(16)

   write(10,*)
   write(10,*)'!** Tower Drag '
   write(10,*)'CDRAG_TOW_el   ',CDRAG_TOW_el
   write(10,*)'Z_DRAG_LIM     ',Z_DRAG_LIM

   write(10,*)'** Nacelle Drag'
   write(10,*)'NAC_LX_el      ',NAC_LX_el
   write(10,*)'NAC_LY_el      ',NAC_LY_el
   write(10,*)'NAC_LZ_el      ',NAC_LZ_el
   write(10,*)'NAC_FROND_X    ',NAC_RX_el
   write(10,*)'NAC_FROND_Y    ',NAC_RX_el
   write(10,*)'NAC_FROND_Z    ',NAC_RX_el
   write(10,*)'NANG_NAC       ',NANG_NAC
  do iang=1,NANG_NAC
   write(10,'(3f15.5)') Alpha_NAC(iang), CL_NAC(iang), CD_NAC(iang)
  enddo !iang
   write(10,*)


 END Subroutine DRAG_INIT_el
!-----------------------------------------------------------
 Subroutine DateTime
!-----------------------------------------------------------
 ! This subroutine writes one character string encoded 
 !  with the date in the form "dd-mmm-ccyy" and one 
 !  character string encoded with the time in the form "hh:mm:ss".
!-----------------------------------------------------------

   implicit none

!--- Local Variables:
   character( 8) :: CDate, CurTime
   character(11) :: CTime, CurDate


!--- Call the system date function.
   call DATE_AND_TIME ( CDate )
   call DATE_AND_TIME ( TIME=CTime )


!--- Parse out the day.
   CurDate(1:3) = CDate(7:8)//'-'


!--- Parse out the month.
   select case ( CDate(5:6) )
     case ( '01' )
        CurDate(4:6) = 'Jan'
     case ( '02' )
        CurDate(4:6) = 'Feb'
     case ( '03' )
        CurDate(4:6) = 'Mar'
     case ( '04' )
        CurDate(4:6) = 'Apr'
     case ( '05' )
        CurDate(4:6) = 'May'
     case ( '06' )
        CurDate(4:6) = 'Jun'
     case ( '07' )
        CurDate(4:6) = 'Jul'
     case ( '08' )
        CurDate(4:6) = 'Aug'
     case ( '09' )
        CurDate(4:6) = 'Sep'
     case ( '10' )
        CurDate(4:6) = 'Oct'
     case ( '11' )
        CurDate(4:6) = 'Nov'
     case ( '12' )
        CurDate(4:6) = 'Dec'
   end select


!--- Parse out the year.
   CurDate(7:11) = '-'//CDate(1:4)
   CurTime       = CTime(1:2)//':'//CTime(3:4)//':'//CTime(5:6)


!--- Output Date and Time of the code execution
   write (10,*) CurDate
   write (10,*) CurTime
   write (10,*)


 END Subroutine DateTime
!-----------------------------------------------------------
 Subroutine Calc_Masses
!-----------------------------------------------------------

 Use Cbeam
 Use Craft
 Use Hydro

   implicit none

   real(8) :: mass_bl, mass_sh, mass_to, mass_ja, mass_RNA, mass_TOW, mass_JAC,mass_fl, mass_T, massflood, Waddmass
   real(8) :: ALLOC, HTA0, R1, R2, LWET, Y0(3), XG(3)
   integer :: nbod,nbsub,nb_el,nn_el,nel,n1,n2,ncal
   real(8) :: y1,y2,Yc,mass_elem,Ycg_elem,mc1,mc2,Ycg_bl,Ycg_sh,Ycg_to,IXX_bl,IXX_sh,IXX_to, LEN0
   real(8) :: IYY_bl,IYY_sh,IYY_to, IZZ_bl,IZZ_sh,IZZ_to
   real(8) :: Ixx_elem, Izz_elem


   mass_bl = 0.d0;  Ycg_bl  = 0.d0;
   mass_sh = 0.d0;  Ycg_sh  = 0.d0;
   mass_to = 0.d0;  Ycg_to  = 0.d0;
   mass_ja = 0.d0
   IXX_bl  = 0.d0;  IYY_bl  = 0.d0;  IZZ_bl  = 0.d0;
   IXX_sh  = 0.d0;  IYY_sh  = 0.d0;  IZZ_sh  = 0.d0;
   IXX_to  = 0.d0;  IYY_to  = 0.d0;  IZZ_to  = 0.d0;
   mass_RNA  = 0.d0 
   mass_TOW  = 0.d0 
   mass_JAC  = 0.d0 
   massflood = 0.d0
   Waddmass  = 0.d0


   do    nbod  = 1,NBODBT_el
         LEN0  = 0.d0
   do    nbsub = 1,body(nbod)%NBODSUB_el
         nb_el = body(nbod)%NBODTGLB_el(nbsub)
      do nel   = 1,subbody(nb_el)%NTEB_el   !all elements of each subbody 

         n1    = subbody(nb_el)%INODNEL_el(nel,1)
         n2    = subbody(nb_el)%INODNEL_el(nel,2)

         y1    = subbody(nb_el)%HTA_el(nel  ) + LEN0
         y2    = subbody(nb_el)%HTA_el(nel+1) + LEN0

!--------- m = mc1*y + mc2
         mass_elem  = ( beam_timo (n1)%DENS_el+beam_timo (n2)%DENS_el ) * subbody(nb_el)%ALENG_el(nel) / 2.d0 !integration
         Ixx_elem   = ( beam_timo (n1)%RIZZ_el+beam_timo (n2)%RIZZ_el ) * subbody(nb_el)%ALENG_el(nel) / 2.d0 !integration 
         Izz_elem   = ( beam_timo (n1)%RIXX_el+beam_timo (n2)%RIXX_el ) * subbody(nb_el)%ALENG_el(nel) / 2.d0 !integration
         mc1        = ( beam_timo (n2)%DENS_el-beam_timo (n1)%DENS_el ) / subbody(nb_el)%ALENG_el(nel)        !linear interpolation coefficients for m=mc1*y+mc2
         mc2        =   beam_timo (n1)%DENS_el-mc1 * y1  
         Ycg_elem   = ( mc1/3.d0*(y2**3-y1**3) + &
                        mc2/2.d0*(y2**2-y1**2)     )/max(mass_elem,1d-15) !Int{m*y*dy}/Int{m*dy}


            Yc      = 0.d0                                    ! Iyy calculated wrt Yc
            if (nbod<=NBLADE_el) Yc = -Hhub                   ! in order to calculate wrt to Hhub

         if     (nbod >  NBODBTWT_el) then
            mass_ja = mass_ja + mass_elem
         elseif (nbod <= NBLADE_el      ) then
            Ycg_bl  = ( mass_elem*Ycg_elem + Ycg_bl*mass_bl ) / max(mass_elem+mass_bl,1.d-15)
            mass_bl = mass_bl + mass_elem
            IXX_bl  = IXX_bl  +  Ixx_elem
            IYY_bl  = IYY_bl  +  mc1/4.d0*(y2**4-y1**4)  +  (mc2-2*Yc*mc1)/3.d0*(y2**3-y1**3) & !Int{m*(y-yc)^2*ds}
                              +  Yc*(mc1*Yc-2*mc2)/2.d0*(y2**2-y1**2)  + Yc**2 *mc2*(y2-y1)
            IZZ_bl  = IZZ_bl  +  Izz_elem
         elseif (nbod == NBLADE_el+1) then
            Ycg_sh  = ( mass_elem*Ycg_elem + Ycg_sh*mass_sh ) / max(mass_elem+mass_sh,1.d-15)
            mass_sh = mass_sh + mass_elem
            IXX_sh  = IXX_sh  +  Ixx_elem
            IYY_sh  = IYY_sh  +  mc1/4.d0*(y2**4-y1**4)  +  (mc2-2*Yc*mc1)/3.d0*(y2**3-y1**3) &
                              +  Yc*(mc1*Yc-2*mc2)/2.d0*(y2**2-y1**2)  + Yc**2 *mc2*(y2-y1)
            IZZ_sh  = IZZ_sh  +  Izz_elem
         elseif (nbod == NBODBTWT_el) then
            Ycg_to  = ( mass_elem*Ycg_elem + Ycg_to*mass_to ) / max(mass_elem+mass_to,1.d-15)
            mass_to = mass_to + mass_elem
            IXX_to  = IXX_to  +  Ixx_elem
            IYY_to  = IYY_to  +  mc1/4.d0*(y2**4-y1**4)  +  (mc2-2*Yc*mc1)/3.d0*(y2**3-y1**3) &
                              +  Yc*(mc1*Yc-2*mc2)/2.d0*(y2**2-y1**2)  + Yc**2 *mc2*(y2-y1)
            IZZ_to  = IZZ_to  +  Izz_elem
         endif

   !!   write(*,*) nbod,nbsub,nel,n1,n2,mass_to
      enddo !nel
      LEN0 = LEN0 + subbody(nb_el)%ALENGB_el
   enddo !nbsub
   enddo !nbod

   write(10,*)
   write(10,*)'*** Calculate Total body-distributed Masses ***'

!--- IXX=$ρ x^2 dS, IYY=$ρ y^2 dS, IZZ=$ρ z^2 dS ---
!--- den einai polikes ropes!, x,y,z local body  ---
!  if (NBODBTWT_el>0) then
     write(10,'(a)')'the blade properties are calculated wrt Hub and not the blade root'
     write(10,'(a)')'**body       mass[kg]         Ycg[m]     IXX[kg*m2]   IYYr - wrt 0        IZZ       IYY -  wrt Ycg   Sy=mass*Ycg'
   if (NBLADE_el>0)              write(10,101) mass_bl/NBLADE_el, Ycg_bl, IXX_bl/NBLADE_el, IYY_bl/NBLADE_el, IZZ_bl/NBLADE_el, (IYY_bl-mass_bl*Ycg_bl**2)/NBLADE_el, mass_bl/NBLADE_el*Ycg_bl
   if (NBLADE_el>0)              write(10,105) mass_bl          , Ycg_bl, IXX_bl          , IYY_bl          , IZZ_bl          , (IYY_bl-mass_bl*Ycg_bl**2)          , 0.d0
   if (NBODBTWT_el>=NBLADE_el+1) write(10,102) mass_sh          , Ycg_sh, IXX_sh          , IYY_sh          , IZZ_sh          ,  IYY_sh-mass_sh*Ycg_sh**2           , mass_sh          *Ycg_sh
   if (NBODBTWT_el>=NBLADE_el+2) write(10,103) mass_to          , Ycg_to, IXX_to          , IYY_to          , IZZ_to          ,  IYY_to-mass_to*Ycg_to**2           , mass_to          *Ycg_to
!  endif
   write(10,104) mass_ja
   write(10,*)

  101 format ('blade ',7f15.3)
  105 format ('rotor ',7f15.3)
  102 format ('shaft ',7f15.3)
  103 format ('tower ',7f15.3)
  104 format ('jacket',7f15.3)


!--- Calculate concentrated masses
   do ncal  = 1, NCOMA_el
      nbod  = conc_mass(ncal)%NBCOMA_el
      nbsub = conc_mass(ncal)%NSBCOMA_el
      nel   = conc_mass(ncal)%NELCOMA_el
      nb_el = body(nbod)%NBODTGLB_el(nbsub)

      if ( nbod <= NBLADE_el   + 1)                                      mass_RNA = mass_RNA + conc_mass(ncal)%Mass
      if ((nbod == NBODBTWT_el    ).and.(nel == subbody(nb_el)%NTEB_el)) mass_RNA = mass_RNA + conc_mass(ncal)%Mass 
      if ((nbod == NBODBTWT_el    ).and.(nel <  subbody(nb_el)%NTEB_el)) mass_TOW = mass_TOW + conc_mass(ncal)%Mass 
      if ( nbod >  NBODBTWT_el    )                                      mass_JAC = mass_JAC + conc_mass(ncal)%Mass
   enddo !ncal

   write(10,*) '**Concentrated**'
   write(10,*)'RNA   ',mass_RNA
   write(10,*)'tower ',mass_TOW
   write(10,*)'jacket',mass_JAC
   write(10,*)
   write(10,*) '**Total**'

   mass_RNA = mass_RNA + mass_bl  + mass_sh
   mass_TOW = mass_TOW + mass_to
   mass_JAC = mass_JAC + mass_ja
   mass_fl  = 0.d0
   if (ICASE_el==3) & 
   mass_fl  = floater(1)%AM_floater_structural(3,3)
   mass_T   = mass_RNA + mass_TOW + mass_JAC + mass_fl


   write(10,*)'RNA    ',mass_RNA
   write(10,*)'tower  ',mass_TOW
   write(10,*)'jacket ',mass_JAC
   write(10,*)'floater',mass_fl
   write(10,*)'TOTAL  ',mass_T


!--- Calculate Flooded mass
   if ((ICASE_el/=2).and.(ICASE_el/=4)) goto 5

   do nbod  = NTOWER_el, NBODBT_el
      if (IFloodbod_el (nbod) == 0) cycle
   do nbsub = 1, body(nbod)%NBODSUB_el
      nb_el = body(nbod)%NBODTGLB_el(nbsub)
      nn_el = body(nbod)%NNODTGLB_el(nbsub)

      do nel       = 1, subbody(nb_el)%NTEB_el
         n1        = subbody(nb_el)%INODNEL_el(nel,1)
         n2        = subbody(nb_el)%INODNEL_el(nel,2)
         HTA0      = subbody(nb_el)%HTA_el    (nel  )
         ALLOC     = subbody(nb_el)%ALENG_el  (nel  )
         R1        = substr_mor(n1)%INDIAMET_el /2.d0
         R2        = substr_mor(n2)%INDIAMET_el /2.d0

         Y0 (1:3)  = 0.d0
         Y0 (2  )  = HTA0
         XG (1:3)  = transf_mat(nn_el)%R_el(1:3) + matmul ( transf_mat(nn_el)%A_el (1:3,1:3), Y0(1:3) )

         if ( XG(3)>= 0.d0) cycle

         Y0 (2  )  = HTA0 + ALLOC
         XG (1:3)  = transf_mat(nn_el)%R_el(1:3) + matmul ( transf_mat(nn_el)%A_el (1:3,1:3), Y0(1:3) )

!--------- calculate wet length
         if (XG(3)<=0.d0) then
            LWET   = ALLOC
         else
            Y0 (2) = HTA0 !Rloc
            XG (3) = 0.d0 !RG

            call Calc_Yloc ( XG, Y0, LWET, nn_el )
            R2     = R1 * (1d0-LWET/ALLOC) + R2 * (LWET/ALLOC)
         endif

         massflood = massflood + rho * PI/3d0 * LWET * (R1**2 + R1*R2 + R2**2)
      enddo
   enddo
   enddo

 5 write(10,*)'Flooded mass',massflood/1000d0


!--- Calculate water added mass to MSL
   do nbod  = NTOWER_el, NBODBT_el
   do nbsub = 1, body(nbod)%NBODSUB_el
      nb_el = body(nbod)%NBODTGLB_el(nbsub)
      nn_el = body(nbod)%NNODTGLB_el(nbsub)

      if ( subbody(nb_el)%IBTYP_el /=2 ) cycle

      do nel       = 1, subbody(nb_el)%NTEB_el
         n1        = subbody(nb_el)%INODNEL_el(nel,1)
         n2        = subbody(nb_el)%INODNEL_el(nel,2)
         HTA0      = subbody(nb_el)%HTA_el    (nel  )
         ALLOC     = subbody(nb_el)%ALENG_el  (nel  )
         R1        = substr_mor(n1)%DIAMET_el / 2.d0
         R2        = substr_mor(n2)%DIAMET_el / 2.d0

         Y0 (1:3)  = 0.d0
         Y0 (2  )  = HTA0
         XG (1:3)  = transf_mat(nn_el)%R_el(1:3) + matmul ( transf_mat(nn_el)%A_el (1:3,1:3), Y0(1:3) )

         if ( XG(3)>= 0.d0) cycle !tide

         Y0 (2  )  = HTA0 + ALLOC
         XG (1:3)  = transf_mat(nn_el)%R_el(1:3) + matmul ( transf_mat(nn_el)%A_el (1:3,1:3), Y0(1:3) )

!--------- calculate wet length
         if (XG(3)<=0.d0) then !tide
            LWET   = ALLOC
         else
            Y0 (2) = HTA0 !Rloc
            XG (3) = 0.d0 !RG

            call Calc_Yloc ( XG, Y0, LWET, nn_el )
            R2     = R1 * (1d0-LWET/ALLOC) + R2 * (LWET/ALLOC)
         endif

         Waddmass  = Waddmass + rho * PI/3d0 * LWET * (R1**2 + R1*R2 + R2**2) !jim lathos an R1/=R2
      enddo !nel
   enddo !nbsub
   enddo !nbod

   write(10,*)'Water Added mass to MSL',Waddmass/1000d0


 END Subroutine Calc_Masses
!
!
!
!----------------------------------------------------------------------
 Subroutine Calc_Yloc ( Rglob, Rloc, Yloc, nn_el )
!----------------------------------------------------------------------
!
!  For a given Global Z (Rglob(3)), calculate the local y distance of 
!  the beam (to match the Rglob(3)). The yloc is defined wrt the begining of the element.
!
!
! {R0} + [A]*[{Rloc}+(0,y,0)T] = {RG}
!
! (0,y,0)T = [AT]*{RG-R0} - {Rloc} , {Rloc} = (u,v+y0,w)T
!
!
!----------------------------------------------------------------------

 Use Hydro
 Use Cbeam

   implicit none


   real(8) :: A(2,2),B(2),Rglob(3),Rloc(3),Yloc,R0(3),det,det1,det2
   integer :: nn_el


!  Rloc (1:3) = UTT_A(IMDOF_el(1:3))
!  Rloc (2  ) = Rloc(2) + HTA0
   R0   (1:3) = transf_mat(nn_el)%R_el (1:3)

   A    (1,1) = transf_mat(nn_el)%AT_el(1,1)
   A    (1,2) = transf_mat(nn_el)%AT_el(1,2)
   A    (2,1) = transf_mat(nn_el)%AT_el(3,1)
   A    (2,2) = transf_mat(nn_el)%AT_el(3,2)

   B    (1)   = Rloc(1) + transf_mat(nn_el)%AT_el(1,1)*R0(1) + transf_mat(nn_el)%AT_el(1,2)*R0(2) - transf_mat(nn_el)%AT_el(1,3) * ( Rglob(3) - R0(3) )
   B    (2)   = Rloc(3) + transf_mat(nn_el)%AT_el(3,1)*R0(1) + transf_mat(nn_el)%AT_el(3,2)*R0(2) - transf_mat(nn_el)%AT_el(3,3) * ( Rglob(3) - R0(3) )

!--- solve 2x2 and find RG1,RG2
   det        = A(1,1)*A(2,2)-A(1,2)*A(2,1)
   det1       = B(1  )*A(2,2)-A(1,2)*B(2  )
   det2       = A(1,1)*B(2  )-B(1  )*A(2,1)

   Rglob (1)  = det1/det
   Rglob (2)  = det2/det

   Yloc       = transf_mat(nn_el)%AT_el(2,1) * ( Rglob(1) - R0(1) ) + &
                transf_mat(nn_el)%AT_el(2,2) * ( Rglob(2) - R0(2) ) + &
                transf_mat(nn_el)%AT_el(2,3) * ( Rglob(3) - R0(3) ) - &
                Rloc(2)


 END Subroutine Calc_Yloc

!----------------------------------------------------------------------
 Subroutine define_external_force
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   real(8) :: F_tmp
   integer :: nbod, nbsub, nb_el, nel


   if (Var_Paper(9)==0.d0) return
!--- The external loads are already zero in allocate_main_matrices
!--- paper: Set the external forces-moments

         F_tmp = Var_Paper(9) + Var_Paper(2) * dcos(Var_Paper(3)*TIME)
   write(10,*)'magnitude of F in define_external_force',F_tmp

   do    nbod  = 1, NBODBT_el
      do nbsub = 1, body(nbod)%NBODSUB_el
         nb_el = body(nbod)%NBODTGLB_el(nbsub)
         nel   = subbody(nb_el)%NTEB_el
!--------- local Fx enabled (edge-wise)
         if (int(Var_Paper(10))== 1) subbody(nb_el)%FCPL_el (1:nel,1) = F_tmp

!--------- local Fy enabled (extension)
         if (int(Var_Paper(10))== 2) subbody(nb_el)%FCPL_el (1:nel,2) = F_tmp

!--------- local Fz enabled (flap-wise)
         if (int(Var_Paper(10))== 3) subbody(nb_el)%FCPL_el (1:nel,3) = F_tmp

!--------- local combined Fx/Fz (edge-wise/flap-wise)
         if (int(Var_Paper(10))==10) subbody(nb_el)%FCPL_el (1:nel,1) = F_tmp
         if (int(Var_Paper(10))==10) subbody(nb_el)%FCPL_el (1:nel,3) = F_tmp

!--------- local My enabled (torsion)
         if (int(Var_Paper(10))==20) subbody(nb_el)%AMPL_el (1:nel,2) = F_tmp
!--------- local My enabled (torsion) at tip
         if (int(Var_Paper(10))==13.and.&
                          nbsub==   body(nbod )%NBODSUB_el    ) subbody(nb_el)%AMPL_el (nel,2) = F_tmp/0.366d0 !subbody(nb_el)%ALENG_el(nel)
      enddo
   enddo


 END Subroutine define_external_force
