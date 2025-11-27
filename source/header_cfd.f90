!----------------------------------------------------------------------------
 Subroutine kinematics_0
!----------------------------------------------------------------------------

   use grid
   use param, only : PI,NNT,NFacesT
   use heli
   use hgast_interface
   use Variables
   use coupling_cfd_0
   
   Implicit None
   
   integer :: ip,is,ifo,inod,IndComp,IndElem,IndSubb,j,iface
   real(8) :: X0L(3),XLL(3),XCL(3),def(3),veldef(3),UB(3)
   real(8) :: XGr(3),VelGr(3)
!
!  Initialization
!  --------------
!
!  we assume that
!G
!  - the blades are straight and that at t=0 they are placed at 0, 90, 180, 270 deg 
!  - the rotor center is at 0,0,0 
!  - the rotor axis is in 0,0,1 wrt the rotor plane and eventually tilted 
!  - the flight velocity is (Ufl,0,0) wrt the global system 
!
!  ####### t=0 #######
!
      allocate (INDX_subbo_nod(NNT),INDX_Eleme_nod(NNT))
      allocate (XGLL_cfd(3,NNT))
      call calculate_force_points(1)
!   a
!  Define IndSubb, IndElemb for the grid and force_el points ### interface with hGAST
   open (1,file='test')
   do ip=1,NNT               
      IndComp=INDX_blade_nod(ip) 
      X0L(:)=XGL_cfd(:,ip)*dimL
      call DefNodeTopology(IndComp,X0L,IndSubb,IndElem,XLL)
           INDX_Subbo_nod(ip)=IndSubb
           INDX_Eleme_nod(ip)=IndElem
           XGLL_cfd(:,ip)=XLL(:) 
  !   write(1,'(10(e28.17,1x),4(i,1x))') X0L(1:3),XLL(1:3), subbody_HTA_el(IndElem,body_NBODTGLB_el(IndSubb,IndComp)),XGr(1:3),&
  !                                                        IndComp,IndSubb,Indelem,body_NBODTGLB_el(IndSubb,IndComp)
!     note: XLL is defined wrt the 1st node of the element it belongs 
   enddo
  !STOP
!
!  force_el Points are defined similarly to the way cfd particles are generated 
!  e.g. for all surface panels/faces divide ...
!
      allocate (INDX_subbo_lop(nforce_el),INDX_Eleme_lop(nforce_el),INDX_blade_lop(nforce_el))
      allocate (XFLL_cfd(3,nforce_el))
   INDX_blade_lop=INDX_blade_force_el
   do is=1,nforce_el ! for all surface panels
      IndComp=INDX_blade_force_el(is)
     !print *, IndComp
! ?? PAPIS add this part
!      - get the nodes
!      - subdivide into NumFoX x NumFoY panels 
!      - define the center of every sub-panel
      XCL(:)= xcforce_el(1:3,is)
! ?? PAPIS
! ### we may also repeat this part whenever necessary but we need
!     to recalculate XFLL_cfd(:,*) and the INDX_*
      X0L(:)=XCL(:)*dimL
      call DefNodeTopology(IndComp,X0L,IndSubb,IndElem,XLL)
         INDX_Subbo_lop(is)=IndSubb
         INDX_Eleme_lop(is)=IndElem
         XFLL_cfd(:,is)=XLL(:)
   enddo


 END subroutine kinematics_0
!----------------------------------------------------------------------------
 Subroutine kinematics(itime)
!----------------------------------------------------------------------------

   use grid
   use param, only : PI,NNT,NFacesT
   use heli
   use hgast_interface
   use Variables
   use mpi
   
   implicit none
   integer,  intent(IN) :: itime
   
   integer              :: ip,is,ifo,inod,IndComp,IndElem,IndSubb,j,iface
   real(8)              :: X0L(3),XLL(3),XCL(3),def(3),veldef(3),UB(3)
   real(8)              :: XGr(3),VelGr(3)
   real(8), allocatable :: velnod(:,:)
   integer              :: my_rank,ierr,p
   character *15        :: filout
!
!   Initialization
!   --------------
!
!   we assume that
!G
!   - the blades are straight and that at t=0 they are placed at 0, 90, 180, 270 deg 
!   - the rotor center is at 0,0,0 
!   - the rotor axis is in 0,0,1 wrt the rotor plane and eventually tilted 
!   - the flight velocity is (Ufl,0,0) wrt the global system 
!
!   ####### t=0 #######
!
! -----------------------------------------------------------------
! 1 Calculate the rigid positions and velocities of all points
   call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)
   call MPI_Comm_size(MPI_COMM_WORLD,p      ,ierr)
  !write(filout,'(a,i2.2)') 'geom',my_rank
  !open(18,file=filout)

    allocate (velnod(3,NNT))
    do ip=1,NNT
        IndComp=INDX_blade_nod(ip) 
        IndSubb=INDX_subbo_nod(ip) 
        IndElem=INDX_eleme_nod(ip) 
        XLL(:) =XGLL_cfd(:,ip)
!   ### for the current problem we include deformations in XGr, VelGr
        call elastic(IndComp,IndSubb,IndElem,XLL,XGr,VelGr)
        grid_po(ip)%x = XGr(1)/dimL
        grid_po(ip)%y = XGr(2)/dimL
        grid_po(ip)%z = XGr(3)/dimL
        velnod(:,ip)  = VelGr(:)/dimV
  !     write(18,'(6(e28.17,1x))') XGr(1:3)/dimL,XGL_cfd(1:3,ip)*dimL
    enddo
  ! close(18)

 ! write(filout,'(a,i2.2)') 'blsrf_init',my_rank
 ! open(18,file=filout)
 !  do ip=1,nforce_el
 !      write(18,'(9(e28.17,1x))') xcforce_el(1:3,ip)*dimL,XFLL_cfd(1:3,ip)
 !  enddo
 !  close(18)
 
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
! -----------------------------------------------------------------
! 2 Define the current grid point positions
!
    call GEOM_CALC_UNSTRUC(ITIME)
    call wall_comm
 !  call calculate_force_points(1)
 ! write(filout,'(a,i2.2)') 'blsrf',my_rank
 ! open(18,file=filout)
 !  do ip=1,nforce_el
 !      write(18,'(9(e28.17,1x))') xcforce_el(1:3,ip)*dimL,XFLL_cfd(1:3,ip)
 !  enddo
 !  close(18)
! -----------------------------------------------------------------
! 3 Define the face averaged velocities 
!
    do is=1,IBcount(12)
       iface=IBF_index(is)
       UB=0
       do j = 1, grid_fa(iface)%NDP
          inod = grid_fa(iface)%Node(j)
          UB(1:3) = UB(1:3) + velnod(1:3,inod)
       enddo
       UB = UB / grid_fa(iface)%NDP
       flux_al(iface)%UBody(1)= UB(1)!+ Vtransx !Uxbody
       flux_al(iface)%UBody(2)= UB(2)!+ Vtransy !Uybody
       flux_al(iface)%UBody(3)= UB(3)!+ Vtransz
    enddo
   deallocate(velnod)


 END subroutine kinematics
!----------------------------------------------------------------------------
 Subroutine Aero2Elast_cfd
!----------------------------------------------------------------------------

   use hgast_interface
   use Cbeam
   use mod_matrix_vars_el
   
   use coupling_cfd_0
   use coupling_cfd_f
   use MPI
!
!   In "Aero2Elast_cfd" we tranfer to hGAST the loads defined by:
!   ip=1,NumForcP: IndComp,IndSubb,IndElem,forc_po(ip)%[x,y,z], 
!                  ForcG(1:3) (Nt)
!   ### this is done in add_load_cfd(...)


    implicit none


    integer          :: ip,IndComp,IndSubb,IndElem
    double precision :: XPG(3),XPL(3),LoaG(3),MomG(3),ALEN,   AT(3,3)
    double precision,allocatable :: F_tmp(:,:,:)
    integer          :: my_rank,np,ierr,status(MPI_STATUS_SIZE),mat3,i
    integer          :: ic,isubb_first,isubb_last,ielem,isubb,IndSubb_gl,  nn_el

    call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)
    call MPI_Comm_size(MPI_COMM_WORLD,np      ,ierr)

    call calculate_force_points(1)
    F_cfd=0
    do ip=1,nforce_el
       IndComp=INDX_blade_force_el(ip) 
       IndSubb=INDX_subbo_lop(ip) 
       IndElem=INDX_eleme_lop(ip) 
       XPG(1) =xcforce_el(1,ip)*dimL
       XPG(2) =xcforce_el(2,ip)*dimL
       XPG(3) =xcforce_el(3,ip)*dimL
       LoaG(1)=force_el(1,ip)*dimF
       LoaG(2)=force_el(2,ip)*dimF
       LoaG(3)=force_el(3,ip)*dimF
       XPL(:) =XFLL_cfd(1:3,ip)
       MomG(:)=0.
       call add_load_cfd(IndComp,IndSubb,IndElem,XPG,XPL,LoaG,MomG) !F_cfd is added inside
    enddo
    
!   Send to zero
    if (my_rank.eq.0) then 
       allocate(F_tmp(NTEB_el_cfd+1,6,NBODT_el))
        call mpimat3el(mat3,NTEB_el_cfd+1,6,NBODT_el)
       do i=2,np
          call MPI_RECV(F_tmp,1,mat3,i-1,1,MPI_COMM_WORLD,status,ierr)
          F_cfd = F_cfd + F_tmp
       enddo
      deallocate(F_tmp)
    else 
       call mpimat3el(mat3,NTEB_el_cfd+1,6,NBODT_el)
       call MPI_SEND(F_cfd,1,mat3,0,1,MPI_COMM_WORLD,ierr)
    endif
    call MPI_TYPE_FREE(mat3,ierr)
   
    if (my_rank.eq.0) then 
      
       do ic=1,NBLADE_el
       isubb_first =body_NBODTGLB_el(1,ic)                   ! body(ic)%NBODTGLB_el(1)
       isubb_last  =body_NBODTGLB_el(body_NBODSUB_el(ic),ic) ! body(ic)%NBODTGLB_el(body(ic)%NBODSUB_el) 

       nn_el           = body(ic)%NNODTGLB_el (1)
       AT    (1:3,1:3) = transf_mat(nn_el)%AT_el    (1:3,1:3) !Local base of 1st subbody
       open (18,file='forceG'//char(48+ic)//'.dat')
       open (19,file='forceL'//char(48+ic)//'.dat')
       ALEN=0
       do isubb=1,body_NBODSUB_el(ic)
          IndSubb_gl=body_NBODTGLB_el(isubb,ic)

          do ielem=1,subbody_NTEB_el(IndSubb_gl)+1                  ! subbody(IndSubb)%NTEB_el

             write (18,'(7(e28.17,1x))') ALEN+subbody_HTA_el(ielem,IndSubb_gl),(F_cfd(ielem,i,IndSubb_gl),i=1,6)
             write (19,'(7(e28.17,1x))') ALEN+subbody_HTA_el      (ielem,    IndSubb_gl)  ,&
                                         matmul( AT(1:3,1:3),F_cfd(ielem,1:3,IndSubb_gl) ),&
                                         matmul( AT(1:3,1:3),F_cfd(ielem,4:6,IndSubb_gl) )
          enddo
          ALEN=ALEN+subbody(IndSubb_gl)%ALENGB_el   
       enddo
       close(18)
       close(19)
      enddo
    endif


 END subroutine aero2elast_cfd
! -----------------------------------------------------------------------------
!
! hGAST - CFD / GenUVP coupling
! -----------------------------
!
! Procedure ...
!
! hGAST		  	                    CFD
! in "init_el.f90"  INITIA_el
! t=0   call INITIA_el (init_el.f90)
!            read input and initialize
!            call broadcast_0 
!            call CFD_main(0,...)      >>>  read input and initialize
!                                               run kinematics_0
!                                           ### run FIRST in genuvp
!                                      <<<      
!       call GELAST0 (gelast.f90)
!       call GELAST_init(1)
!            call broadcast_t
!            call CFD_main(1,...)      >>>  get an initial set of loads 
!                                               run one or more steps or read back-up
!                                           ### run POTENTC/RESCPV in genuvp
!                                               run aero2elast_cfd
!                                      <<<
!            run hGAST 
!        ... ICONV=1                   >>>  conclude the time step, prepare the next
!                 
! t>0,  NTIME/TIME, ICONV=0
!       call GELAST
!   [*] do elastic iterations
!            prepare iteration           
!            call broadcast_t
!            call CFD_main(2,...)      >>>  perform a step
!                                               run kinematics_t
!                                               run one step
!                                           ### run POTENTC/RESCPV in genuvp
!                                               run aero2elast_cfd
!                                      <<<
!            run hGAST 
!            check error
!            if ICONV=0 go to [*]
!               else                   >>>  conclude the time step
!                                           output results 
!                                           prepare the next step
!                                      <<<          
!       enddo iter_el
!       output elastic results
! enddo t     
!
! ---------------------------------------------------------------------------
! new in CFD
! ---------------------------------------------------------------------------
! EVERYTHING CAN BE DEFINED PER PROCESSOR
! 
!     NNT
! defined in kinematics and used in kinematics (we keep them)
!     INDX_blade_nod(NNT), INDX_Subbo_nod(NNT), INDX_Eleme_nod(NNT)
!
!     XGL_cfd(3,NNT)             !used in kinematics(t=0) and nowhere else
!     XGLL_cfd(3,NNT)            !defined in kinematics and used in kinematics 
!                                !        ### we must keep it 
!     velnod(1:3,1:NNT)          !defined in kinematics and used locally
!        grid_po(NNT)%x,y,z      !defined in kinematics and used in GEOM_CALC 
!        flux_al(iface)%UBody(3) !defined in kinematics  
!
!     NumForcP
! defined in kinematics and used in Aero2Elast_cfd (we keep them)
!     INDX_blade_lod(NumForcP), INDX_Subbo_lod(NumForcP), INDX_Eleme_lod(NumForcP)
!     XFLL_cfd(3,NumForcP)       !defined in kinematics and used in Aero2Elast_cfd
!     forc_po(NumForcP)%x,y,z    !it is not really needed and can be recalculated 
!                                !   by re-dividing the surface panels  
!     ForcG(3,NumForcP)          != [stress}^T*Normal*Surface/NumFPperPan
!`
! -----------------------------------------------------------------------------
! from hGAST
! -----------------------------------------------------------------------------
! + NDFPEM_el     , NEQPEM_el  , NNPE_el    , NBODBT_el, NBODT_el, NDFT_el, NNODT_el
!   NBODSUB_el_cfd, NTEB_el_cfd, NNTB_el_cfd
!
! + UT_el(NDFT_el), UT1_el(NDFT_el)
! 
! + body_NBODSUB_el  (    NBODBT_el)--> body(NBODBT_el)%NBODSUB_el 
! + body_NBODTGLB_el (*  ,NBODBT_el)--> body(NBODBT_el)%NBODTGLB_el(*),      *=body(ib)%NBODSUB_el
! + body_NNODTGLB_el (*+1,NBODBT_el)--> body(NBODBT_el)%NNODTGLB_el(*+1)     *max=NBODSUB_el_cfd
!
! + subbody_NTEB_el  (    NBODT_el) --> subbody(NBODT_el)%NTEB_el
! + subbody_NEQPE_el (    NBODT_el) --> subbody(NBODT_el)%NEQPE_el
! + subbody_NDFPE_el (    NBODT_el) --> subbody(NBODT_el)%NDFPE_el
!
! + subbody_PHIX_el  (*  ,NBODT_el) --> subbody(NBODT_el)%PHIX_el(*)         *=subbody(isub)%NTEB_el
! + subbody_PHIZ_el  (*  ,NBODT_el) --> subbody(NBODT_el)%PHIZ_el(*)         *max=NTEB_el_cfd
! + subbody_ALENG_el (*  ,NBODT_el) --> subbody(NBODT_el)%ALENG_el(*)      
! + subbody_HTA_el   (*+1,NBODT_el) --> subbody(NBODT_el)%HTA_el (*+1)
!   
! subbody_NDFPNBACC_el(0:*,NBODT_el) -> subbody(NBODT_el)%NDFPNBACC(0:*)     *=subbody(isub)%NNTB_el
!                                                                            *max=NNTB_el_cfd=NTEB_el_cfd*(NNPE_el-1)+1
! subbody_NODMB(*,NNPE_el,NBODT_el) --> subbody(NBODT_el)%NODMB(*,NNPE_el)   *=subbody(isub)%NTEB_el
!                                                                            *max=NTEB_el_cfd
!         NDFPNBACC_el(0:#), #=n*(NNPE_el-1)+1
!
! + ACCU_NDFTACC_el (0:NBODT_el)  ---> ACCU%NDFTACC_el(0:NBODT_el) 
!
! + transf_mat_R_el (3,NNODT_el)   ---> transf_mat(NNODT_el)%R_el(3)
!              DR_el(3,NNODT_el)   --->                      DR_el(3))
!              A_el (3,3,NNODT_el) --->                      A_el(3,3)
!              DA_el(3,3,NNODT_el) --->                      DA_el(3,3)
!             AT_el(3,3,NNODT_el)  --->                      AT_el(3,3)
! 
! + F_cfd(#+1,6,NBODT_el)            --->                                   #=NTEB_el_cfd
! 
!----------------------------------------------------------------------------
 Subroutine elastic (IndComp,IndSubb,IndElem,XLL,XGr,VGr)
!----------------------------------------------------------------------------

    use Cbeam,              only: NDFPEM_el,NEQPEM_el,R2D
    use mod_matrix_vars_el, only: UT_el,UT1_el

    use coupling_cfd_0

    implicit none

    integer, intent(IN)  :: IndComp,IndSubb,IndElem
    real(8), intent(IN)  :: XLL(3)
    real(8), intent(OUT) :: XGr(3),VGr(3)

    integer :: inode,NEQPE,NDFPE,i,nn,IndSubb_gl
    real(8) :: AKSI0(3),ALLOC,HTAL,DefL(3),DtDefL(3), &
               PHPX,PHPZ,SS(3,6)
    real(8) :: R_el(3),DR_el(3),A_el(3,3),DA_el(3,3)
    real(8) :: SHAPE (NEQPEM_el,NDFPEM_el)
    real(8) :: UT    (NDFPEM_el ), UTT0(NEQPEM_el)
    real(8) :: UT1   (NDFPEM_el ), UTT1(NEQPEM_el)
    real(8) :: UT2   (NDFPEM_el )
         

    IndSubb_gl     = body_NBODTGLB_el(IndSubb,IndComp)             !nb_el = body(nbod_el)%NBODTGLB_el(nbsub)

    call Calc_SS(SS,XLL)

    ALLOC          = subbody_ALENG_el(IndElem,IndSubb_gl)          !subbody(nb_el)%ALENG_el  (nel)
    HTAL           = XLL(2)/ALLOC                                  ! HTA/ALLOC
    PHPX           = subbody_PHIX_el (IndElem,IndSubb_gl)          !subbody(nb_el)%PHIX_el   (nel)
    PHPZ           = subbody_PHIZ_el (IndElem,IndSubb_gl)          !subbody(nb_el)%PHIZ_el   (nel)
    NEQPE          = subbody_NEQPE_el(        IndSubb_gl)          !subbody(nb_el)%NEQPE_el
    NDFPE          = subbody_NDFPE_el(        IndSubb_gl)          !subbody(nb_el)%NDFPE_el

    call SSHAPEFUNC15 (HTAL,ALLOC,PHPX,PHPZ,SHAPE,NDFPE,NEQPE)

    nn             = subbody_NODMB(IndElem,1,IndSubb_gl)           !subbody(nb_el)%NODMB     (nel,1) !nnod)
    i              = ACCU_NDFTACC_el(IndSubb_gl-1) + subbody_NDFPNBACC_el(nn-1,IndSubb_gl) !ACCU%NDFTACC_el(nb_el-1) + subbody(nb_el)%NDFPNBACC_el(nn-1) !+ndof
    UT (1:NDFPE)   = UT_el (i+1:i+NDFPE)
    UT1(1:NDFPE)   = UT1_el(i+1:i+NDFPE)

    UTT0           = matmul(SHAPE, UT                    )
    UTT1           = matmul(SHAPE, UT1                   )
    DefL           = matmul(SS   , UTT0(IMDOF_elcfd(1:6)))
    DtDefL         = matmul(SS   , UTT1(IMDOF_elcfd(1:6)))

    inode          = body_NNODTGLB_el(IndSubb,IndComp)             !nn_el = body(nbod_el)%NNODTGLB_el(nbsub)
    R_el (1:3    ) = transf_mat_R_el (1:3    ,inode)               !transf_mat(nn_el)% R_el(1:3)
    DR_el(1:3    ) = transf_mat_DR_el(1:3    ,inode)               !transf_mat(nn_el)%DR_el(1:3)
    A_el (1:3,1:3) = transf_mat_A_el (1:3,1:3,inode)               !transf_mat(nn_el)% A_el(1:3,1:3)
    DA_el(1:3,1:3) = transf_mat_DA_el(1:3,1:3,inode)               !transf_mat(nn_el)%DA_el(1:3,1:3)
    AKSI0          = XLL;
    AKSI0(2)       = AKSI0(2) + subbody_HTA_el(IndElem,IndSubb_gl) !HTA0 = subbody(nb_el)%HTA_el(nel)

!---- XG =  RG +  A.[(SS.UTT0)+AKSI0]
!---- VG = DRG + DA.[(SS.UTT0)+AKSI0] + A.[(SS.DUTT0)]
    XGr            =  R_el + matmul(A_el,AKSI0+  DefL)
    VGr            = DR_el + matmul(A_el,      DtDefL) + matmul(DA_el,AKSI0+DefL)


 END subroutine elastic
!----------------------------------------------------------------------------
 Subroutine add_load_cfd(IndComp,IndSubb,IndElem,XPG,XPL,LoaG,MomG)
!----------------------------------------------------------------------------
!
!   add_load_cfd will add LoaG,MomG to subbody(IndSubb)%F_cfd(1:NTEB_el+1,1:6)
!
!   XPG = global position
!   XPL = initial local position wrt the 1st node of element IndElem
!
    use Cbeam,              only: NDFPEM_el,NEQPEM_el
    use mod_matrix_vars_el, only: UT_el,UT1_el

    use coupling_cfd_0
    use coupling_cfd_f

    implicit none

    integer, intent(IN) :: IndComp,IndSubb,IndElem
    real(8), intent(IN) :: XPG(3),XPL(3),LoaG(3),MomG(3)

    integer :: inode,NEQPE,NDFPE,i,nn, IndSubb_gl
    real(8) :: AKSI0(3),SS(3,6),ALLOC,HTAL,PHPX,PHPZ,DefL(3)
    real(8) :: R_el(3),A_el(3,3),XG1(3),VG1(3),RR(3),Mom(3),HTAF,Fun1,Fun2,AT_el(3,3),MomL(3),LoaL(3)


    if (IndSubb>1) then; stop; write(*,*)'only 1 sbod supported atm in add_loads_cfd'; endif

    IndSubb_gl     = body_NBODTGLB_el(IndSubb,IndComp)             !nb_el = body(nbod_el)%NBODTGLB_el(nbsub)
!---- 1. Define XG1 point along the deformed elastic axis where the aero loads will be transfered
    AKSI0    = 0.d0;
    AKSI0(2) = XPL(2)
    ALLOC          = subbody_ALENG_el(IndElem,IndSubb_gl)          !subbody(nb_el)%ALENG_el  (nel)

    call elastic (IndComp,IndSubb,IndElem,AKSI0,XG1,VG1)

!---- 2. Transfer LoadG, MomG from point XPG to XG1
    RR (1:3) = XPG (1:3)-XG1(1:3);

    call EXTEPR (RR,LoaG,Mom)

    Mom(1:3) = MomG(1:3)+Mom(1:3);


!---- 3. Project LoaG, Mom
    HTAF=XPL(2); Fun1=1.d0-HTAF/ALLOC; Fun2=HTAF/ALLOC

    F_cfd(IndElem  , 1:3, IndSubb_gl) = F_cfd(IndElem  , 1:3,IndSubb_gl) + Fun1*LoaG(1:3)/ALLOC !qq consider the deformed length
    F_cfd(IndElem+1, 1:3, IndSubb_gl) = F_cfd(IndElem+1, 1:3,IndSubb_gl) + Fun2*LoaG(1:3)/ALLOC
    F_cfd(IndElem  , 4:6, IndSubb_gl) = F_cfd(IndElem  , 4:6,IndSubb_gl) + Fun1*Mom (1:3)/ALLOC
    F_cfd(IndElem+1, 4:6, IndSubb_gl) = F_cfd(IndElem+1, 4:6,IndSubb_gl) + Fun2*Mom (1:3)/ALLOC

!---- Take into account the other subbodies
!  if     (indElem==1.and.IndSubb_gl>1) then
!   F_cfd(IndElem  , 1:3, IndSubb_gl) = F_cfd(IndElem  , 1:3,IndSubb_gl) + Fun1*LoaG(1:3)/ALLOC !qq consider the deformed length
!   F_cfd(IndElem+1, 1:3, IndSubb_gl) = F_cfd(IndElem+1, 1:3,IndSubb_gl) + Fun2*LoaG(1:3)/ALLOC
!   F_cfd(IndElem  , 4:6, IndSubb_gl) = F_cfd(IndElem  , 4:6,IndSubb_gl) + Fun1*Mom (1:3)/ALLOC
!   F_cfd(IndElem+1, 4:6, IndSubb_gl) = F_cfd(IndElem+1, 4:6,IndSubb_gl) + Fun2*Mom (1:3)/ALLOC
!  elseif (indElem==maxElem.and.IndSubb_gl<MaxSb) then
!   F_cfd(IndElem  , 1:3, IndSubb_gl) = F_cfd(IndElem  , 1:3,IndSubb_gl) + Fun1*LoaG(1:3)/ALLOC !qq consider the deformed length
!   F_cfd(IndElem+1, 1:3, IndSubb_gl) = F_cfd(IndElem+1, 1:3,IndSubb_gl) + Fun2*LoaG(1:3)/ALLOC
!   F_cfd(IndElem  , 4:6, IndSubb_gl) = F_cfd(IndElem  , 4:6,IndSubb_gl) + Fun1*Mom (1:3)/ALLOC
!   F_cfd(IndElem+1, 4:6, IndSubb_gl) = F_cfd(IndElem+1, 4:6,IndSubb_gl) + Fun2*Mom (1:3)/ALLOC
!  endif

!   do isubb    = 1, body_NBODSUB_el(ic)
!      inode    = body_NNODTGLB_el(isubb,ic)                     !body(ic)%NNODTGLB_el(isubb) nn_el
!      isubb_gl = body_NBODTGLB_el(isubb,ic)                                                 !nb_el   (nbod,nbsub,nel)
!      nelem    = subbody_NTEB_el (isubb_gl)                     !subbody(isubb)%NTEB_el


 END subroutine add_load_cfd 
!----------------------------------------------------------------------------
 Subroutine DefNodeTopology(IndComp,XL,IndSubb,IndElem,XLL)
!----------------------------------------------------------------------------
!   
!   Called once at t=0 only for the surface nodes
!                               the force points
!
!   input:  IndComp: index of the component XL belongs to
!           XL     : initial (local) position of a point P
!
!   output: IndSubb: index of sub-body point P belongs
!           IndElem: index of element point P belongs
!           XLL    : local coordinates wrt IndElem 
!
    use coupling_cfd_0

    implicit none

    integer, intent(IN)  :: IndComp
    integer, intent(OUT) :: IndSubb,IndElem
    real(8), intent(IN)  :: XL(3)
    real(8), intent(OUT) :: XLL(3)

    integer :: ic,isubb,isubb_first,isubb_last, &
                  inode,inode_first,            nelem,ielem,i, isubb_gl,IndSubb_gl
    real(8) :: HTA_end,XG(3),R_el(3),A_el(3,3),AT_el(3,3)
!Indcomp points to local subbody
!
    ic=IndComp
!    
!   get the global position XG of P 
    inode_first =body_NNODTGLB_el(1,ic)                   ! body(ic)%NNODTGLB_el(isubb_first)
    R_el(1:3  ) = transf_mat_R_el(1:3,inode_first)
    A_el(1:3,1:3) = transf_mat_A_el(1:3,1:3,inode_first)
    XG = R_el+matmul(A_el,XL)
!   1. Determine IndSubb 

       IndSubb  = 0
    do isubb    = 1, body_NBODSUB_el(ic)
       inode    = body_NNODTGLB_el(isubb,ic)                     !body(ic)%NNODTGLB_el(isubb) nn_el
       isubb_gl = body_NBODTGLB_el(isubb,ic)                                                 !nb_el   (nbod,nbsub,nel)
       nelem    = subbody_NTEB_el (isubb_gl)                     !subbody(isubb)%NTEB_el
       HTA_end  = subbody_HTA_el  (nelem+1,isubb_gl)
!   for every subbody:
!   a) Project XG and obtain XLL         
!         Xg=Rsb1+Asb1.XL
!         Xg=Rsb+Asb.XLL -> XLL=ATsb.(Xg-Rsb) XLL here is defined wrt sb
        R_el (1:3    ) = transf_mat_R_el (1:3,    inode)
        AT_el(1:3,1:3) = transf_mat_AT_el(1:3,1:3,inode)
        XLL            = matmul(AT_el,XG-R_el) 
!   b) Check if the point belongs to isubb 
!      ### The beam axis of a subbody is straight: HTA=[0,HTA_end]
       if (XLL(2) <= HTA_end) then
          IndSubb    = isubb
          IndSubb_gl = body_NBODTGLB_el(IndSubb,ic)
          goto 2
       endif
    enddo
!   c) Correct for the last "subbody"
    if (IndSubb.eq.0) then
!      write(*,*)'DefNodeTopology correction for last sbod'
       IndSubb    = body_NBODSUB_el(ic)
       IndSubb_gl = body_NBODTGLB_el(IndSubb,ic)
    endif 

!   2. Determine IndElem 

 2     IndElem=0
    do ielem=1,subbody_NTEB_el(IndSubb_gl)                    ! subbody(IndSubb)%NTEB_el
       if (XLL(2).lt.subbody_HTA_el(ielem+1,IndSubb_gl)) then ! subbody(IndSubb)%HTA_el(ielem)
          IndElem=ielem
          goto 3
       endif
    enddo
!   correct for the last "element"
    if (IndElem.eq.0) then
       IndElem=subbody_NTEB_el(IndSubb_gl)                    ! subbody(IndSubb)%NTEB_el
    endif
 3  XLL(2)=XLL(2)-subbody_HTA_el(IndElem,IndSubb_gl)          ! subbody(IndSubb)%HTA_el(IndElem)
!  finaly XLL is defined wrt nel
    
!   write(*,*) 'End of DefNodeTopology'
    

 END subroutine DefNodeTopology
!----------------------------------------------------------------------------
  Subroutine broadcast_0
!----------------------------------------------------------------------------
!
! allocates, transfers info to all variables and broadcasts
!
  use Cbeam
  use mod_matrix_vars_el

  use coupling_cfd_0
  use coupling_cfd_f

#ifdef HAVE_MPI
  use MPI
#endif

  Implicit None

  integer :: imax_cfd,i_cfd,nn,k,j,mat2,mat3

#ifdef HAVE_MPI
integer :: my_rank,ierr
call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)
#endif

! allocate and transfer 
! body
  call MPI_BCAST( NBODT_el      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
  call MPI_BCAST( NBODBT_el      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
! call MPI_BCAST( NNPE_el       ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
  allocate (body_NBODSUB_el(NBODBT_el)) 
            imax_cfd=0
            if (my_rank.eq.0) then
                do i_cfd=1,NBODBT_el
                   if (body(i_cfd)%NBODSUB_el.gt.imax_cfd) imax_cfd=body(i_cfd)%NBODSUB_el
                enddo
                NBODSUB_el_cfd=imax_cfd
            endif         

  call MPI_BCAST( NBODSUB_el_cfd,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 

  allocate (body_NBODTGLB_el(NBODSUB_el_cfd  ,NBODBT_el), &
            body_NNODTGLB_el(NBODSUB_el_cfd+1,NBODBT_el)) 
! transfer info
            if (my_rank.eq.0) then 
                body_NBODSUB_el(1:NBODBT_el)=body(1:NBODBT_el)%NBODSUB_el
                do i_cfd=1,NBODBT_el
                     body_NBODTGLB_el(1:body(i_cfd)%NBODSUB_el  ,i_cfd)= &
                          body(i_cfd)%NBODTGLB_el(1:body(i_cfd)%NBODSUB_el  )
                     body_NNODTGLB_el(1:body(i_cfd)%NBODSUB_el+1,i_cfd)= &
                          body(i_cfd)%NNODTGLB_el(1:body(i_cfd)%NBODSUB_el+1)
                enddo
            endif


! subbody
  allocate (subbody_NTEB_el(NBODT_el), subbody_NEQPE_el(NBODT_el), &
            subbody_NDFPE_el(NBODT_el))     
            if (my_rank.eq.0) then 
                imax_cfd=0
                do i_cfd=1,NBODT_el
                   if (subbody(i_cfd)%NTEB_el.gt.imax_cfd) imax_cfd=subbody(i_cfd)%NTEB_el
                enddo
                NTEB_el_cfd=imax_cfd
                NNTB_el_cfd=imax_cfd*(NNPE_el-1)+1
            endif

  call MPI_BCAST( NTEB_el_cfd   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
  call MPI_BCAST( NNTB_el_cfd   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
  allocate (subbody_NODMB(NTEB_el_cfd,NNPE_el,NBODT_el), &
            subbody_NDFPNBACC_el(0:NNTB_el_cfd,NBODT_el), &
            subbody_HTA_el  (NTEB_el_cfd+1,NBODT_el), &
            subbody_ALENG_el(NTEB_el_cfd,  NBODT_el),     &
            subbody_PHIX_el (NTEB_el_cfd,  NBODT_el), &
            subbody_PHIZ_el (NTEB_el_cfd,  NBODT_el))

! transfer info
            if (my_rank.eq.0) then 
                subbody_NTEB_el (1:NBODT_el) = subbody(1:NBODT_el)%NTEB_el
                subbody_NEQPE_el(1:NBODT_el) = subbody(1:NBODT_el)%NEQPE_el
                subbody_NDFPE_el(1:NBODT_el) = subbody(1:NBODT_el)%NDFPE_el
                do i_cfd=1,NBODT_el
                   subbody_HTA_el(1:subbody(i_cfd)%NTEB_el+1,i_cfd) = &
                                  subbody(i_cfd)%HTA_el(1:subbody(i_cfd)%NTEB_el)
                   subbody_ALENG_el(1:subbody(i_cfd)%NTEB_el,i_cfd) = &
                                  subbody(i_cfd)%ALENG_el(1:subbody(i_cfd)%NTEB_el)
                   subbody_PHIX_el(1:subbody(i_cfd)%NTEB_el,i_cfd) = &
                                  subbody(i_cfd)%PHIX_el(1:subbody(i_cfd)%NTEB_el)
                   subbody_PHIZ_el(1:subbody(i_cfd)%NTEB_el,i_cfd) = &
                                  subbody(i_cfd)%PHIZ_el(1:subbody(i_cfd)%NTEB_el)
                nn=subbody(i_cfd)%NTEB_el*(NNPE_el-1)+1
                   subbody_NDFPNBACC_el(0:nn,i_cfd) = subbody(i_cfd)%NDFPNBACC_el(0:nn) 
                   do k=1,NNPE_el
                   subbody_NODMB(1:subbody(i_cfd)%NTEB_el,k,i_cfd)= &
                           subbody(i_cfd)%NODMB(1:subbody(i_cfd)%NTEB_el,k)
                   enddo
                enddo
            endif
               
! bradcast: sizes
! call MPI_BCAST( NDFPEM_el     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
! call MPI_BCAST( NEQPEM_el     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
  call MPI_BCAST( NDFT_el       ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
  call MPI_BCAST( NNODT_el      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
! ACCU

  allocate (ACCU_NDFTACC_el(0:NBODT_el))
           if (my_rank.eq.0) then 
               ACCU_NDFTACC_el(0:NBODT_el)= ACCU%NDFTACC_el(0:NBODT_el)
           endif
   
! transf_mat
  allocate (transf_mat_R_el(3,  NNODT_el), transf_mat_DR_el(  3,NNODT_el), &
            transf_mat_A_el(3,3,NNODT_el), transf_mat_AT_el(3,3,NNODT_el), &
                                           transf_mat_DA_el(3,3,NNODT_el))
  if (my_rank.eq.0) then 
      do nn=1,NNODT_el
         transf_mat_R_el (1:3,nn)=transf_mat(nn)%R_el(1:3)
         transf_mat_DR_el(1:3,nn)=transf_mat(nn)%DR_el(1:3)
         do j=1,3
            transf_mat_A_el (1:3,j,nn)=transf_mat(nn)%A_el (1:3,j)
            transf_mat_AT_el(1:3,j,nn)=transf_mat(nn)%AT_el(1:3,j)
            transf_mat_DA_el(1:3,j,nn)=transf_mat(nn)%DA_el(1:3,j)
         enddo
      enddo
   endif

! Forces: 
  allocate (F_cfd(NTEB_el_cfd+1,6,NBODT_el))
  F_cfd=0.
 
 if (my_rank.ne.0) allocate(UT_el(NDFT_el),UT1_el(NDFT_el))

! bradcast: arrays 

! x              body_NBODSUB_el(  NBODBT_el)) 
  call MPI_BCAST(body_NBODSUB_el(1:NBODBT_el),NBODBT_el ,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

! x     body_NBODTGLB_el(NBODSUB_el_cfd  ,NBODBT_el)
  call mpimat2elint(mat2,  NBODSUB_el_cfd  ,NBODBT_el)
  call MPI_BCAST (body_NBODTGLB_el,1,mat2,0,MPI_COMM_WORLD,ierr)
  call MPI_TYPE_FREE(mat2,ierr)

! x     body_NNODTGLB_el(NBODSUB_el_cfd+1,NBODBT_el) 
  call mpimat2elint(mat2,  NBODSUB_el_cfd+1,NBODBT_el)
  call MPI_BCAST (body_NNODTGLB_el,1,mat2,0,MPI_COMM_WORLD,ierr)
  call MPI_TYPE_FREE(mat2,ierr)

! x              subbody_NTEB_el(NBODT_el)
  call MPI_BCAST(subbody_NTEB_el(1:NBODT_el),NBODT_el,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
! x              subbody_NEQPE_el(NBODT_el)
  call MPI_BCAST(subbody_NEQPE_el(1:NBODT_el),NBODT_el,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
! x              subbody_NDFPE_el(NBODT_el))     
  call MPI_BCAST(subbody_NDFPE_el(1:NBODT_el),NBODT_el,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

! x               subbody_NDFPNBACC_el(0:NNTB_el_cfd,NBODT_el), &

  call mpimat2elint(mat2,NNTB_el_cfd+1,NBODT_el) !CHECK IF IT IS DONE CORRECTLY because it starts from zero
  call MPI_BCAST (subbody_NDFPNBACC_el,1,mat2,0,MPI_COMM_WORLD,ierr)
  call MPI_TYPE_FREE(mat2,ierr)

!                 subbody_NODMB(NTEB_el_cfd,NNPE_el,NBODT_el), &
  call mpimat3elint(mat3,NTEB_el_cfd,NNPE_el,NBODT_el) !CHECK IF IT IS DONE CORRECTLY because it starts from zero
  call MPI_BCAST (subbody_NODMB,1,mat3,0,MPI_COMM_WORLD,ierr)
  call MPI_TYPE_FREE(mat3,ierr)

! x               subbody_HTA_el(NTEB_el_cfd+1,NBODT_el), &
  call mpimat2el(mat2,             NTEB_el_cfd+1,NBODT_el) !CHECK IF IT IS DONE CORRECTLY
  call MPI_BCAST (subbody_HTA_el,1,mat2,0,MPI_COMM_WORLD,ierr)
  call MPI_TYPE_FREE(mat2,ierr)

! x               subbody_ALENG_el(NTEB_el_cfd,NBODT_el),     &
  call mpimat2el(mat2,               NTEB_el_cfd,NBODT_el) !CHECK IF IT IS DONE CORRECTLY
  call MPI_BCAST (subbody_ALENG_el,1,mat2,0,MPI_COMM_WORLD,ierr)
  call MPI_TYPE_FREE(mat2,ierr)

! x               subbody_PHIX_el(NTEB_el_cfd,NBODT_el), &
  call mpimat2el(mat2,              NTEB_el_cfd,NBODT_el) !CHECK IF IT IS DONE CORRECTLY
  call MPI_BCAST (subbody_PHIX_el,1,mat2,0,MPI_COMM_WORLD,ierr)
  call MPI_TYPE_FREE(mat2,ierr)

! x               subbody_PHIZ_el(NTEB_el_cfd,NBODT_el))
  call mpimat2el(mat2,              NTEB_el_cfd,NBODT_el) !CHECK IF IT IS DONE CORRECTLY
  call MPI_BCAST (subbody_PHIZ_el,1,mat2,0,MPI_COMM_WORLD,ierr)
  call MPI_TYPE_FREE(mat2,ierr)
  call MPI_BCAST(ACCU_NDFTACC_el(0:NBODT_el), NBODT_el+1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)


! transf_mat ...
  call mpimat2el(mat2,3,NNODT_el)
  call MPI_BCAST (transf_mat_R_el,1,mat2,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST (transf_mat_DR_el,1,mat2,0,MPI_COMM_WORLD,ierr)
  call MPI_TYPE_FREE(mat2,ierr)

  call mpimat3el(mat3,3,3,NNODT_el)
  call MPI_BCAST (transf_mat_A_el,1,mat3,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST (transf_mat_AT_el,1,mat3,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST (transf_mat_DA_el,1,mat3,0,MPI_COMM_WORLD,ierr)
  call MPI_TYPE_FREE(mat3,ierr)
         
!                UT_el(NDFT_el), UT1_el(NDFT_el)
  if (my_rank.eq.0) IMDOF_elcfd=IMDOF_el
  call MPI_BCAST(IMDOF_elcfd (1:6), 6,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(UT_el (1:NDFT_el), NDFT_el,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(UT1_el(1:NDFT_el), NDFT_el,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

! F_cfd(NTEB_el_cfd+1,6,NBODT_el)
! call mpimat3el(mat3,NTEB_el_cfd+1,6,NBODT_el)
! call MPI_BCAST(F_cfd,1,mat3,0,MPI_COMM_WORLD,ierr)
! call MPI_TYPE_FREE(mat3,ierr)


 END subroutine broadcast_0
!----------------------------------------------------------------------------
  Subroutine broadcast_t
!----------------------------------------------------------------------------
!
! updates variables and broadcasts
!
  use Cbeam
  use mod_matrix_vars_el

  use coupling_cfd_0
  use coupling_cfd_f

#ifdef HAVE_MPI
  use MPI
#endif

  Implicit None

  integer :: imax_cfd,i_cfd,nn,j,mat2,mat3

#ifdef HAVE_MPI
integer :: my_rank,ierr
call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)
#endif

!--------------------
! PHIX_el, PHIZ_el
  if (my_rank.eq.0) then 
     do i_cfd=1,NBODT_el
        subbody_PHIX_el(1:subbody(i_cfd)%NTEB_el,i_cfd) = &
                          subbody(i_cfd)%PHIX_el(1:subbody(i_cfd)%NTEB_el)
        subbody_PHIZ_el(1:subbody(i_cfd)%NTEB_el,i_cfd) = &
                          subbody(i_cfd)%PHIZ_el(1:subbody(i_cfd)%NTEB_el)
     enddo
  endif
           call mpimat2el(mat2,NTEB_el_cfd,NBODT_el) !CHECK IF IT IS DONE CORRECTLY
           call MPI_BCAST (subbody_PHIX_el,1,mat2,0,MPI_COMM_WORLD,ierr)
           call MPI_TYPE_FREE(mat2,ierr)

           call mpimat2el(mat2,NTEB_el_cfd,NBODT_el) !CHECK IF IT IS DONE CORRECTLY
           call MPI_BCAST (subbody_PHIZ_el,1,mat2,0,MPI_COMM_WORLD,ierr)
           call MPI_TYPE_FREE(mat2,ierr)
    
! transf_mat
  if (my_rank.eq.0) then 
      do nn=1,NNODT_el
         transf_mat_R_el(1:3,nn)=transf_mat(nn)%R_el(1:3)
         transf_mat_DR_el(1:3,nn)=transf_mat(nn)%DR_el(1:3)
         do j=1,3
         transf_mat_A_el (1:3,j,nn)=transf_mat(nn)%A_el (1:3,j)
         transf_mat_AT_el(1:3,j,nn)=transf_mat(nn)%AT_el(1:3,j)
         transf_mat_DA_el(1:3,j,nn)=transf_mat(nn)%DA_el(1:3,j)
         enddo
      enddo
  endif

! transf_mat ...
  call mpimat2el(mat2,3,NNODT_el)
  call MPI_BCAST (transf_mat_R_el,1,mat2,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST (transf_mat_DR_el,1,mat2,0,MPI_COMM_WORLD,ierr)
  call MPI_TYPE_FREE(mat2,ierr)

  call mpimat3el(mat3,3,3,NNODT_el)
  call MPI_BCAST (transf_mat_A_el,1,mat3,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST (transf_mat_AT_el,1,mat3,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST (transf_mat_DA_el,1,mat3,0,MPI_COMM_WORLD,ierr)
  call MPI_TYPE_FREE(mat3,ierr)
         
! UT_el(NDFT_el), UT1_el(NDFT_el)
  call MPI_BCAST(UT_el (1:NDFT_el), NDFT_el,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(UT1_el(1:NDFT_el), NDFT_el,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

! Forces: F_cfd(NTEB_el_cfd+1,6,NBODT_el)
 !F_cfd=0.
 !call mpimat3el(mat3,NTEB_el_cfd+1,6,NBODT_el)
 !call MPI_BCAST(F_cfd,1,mat3,0,MPI_COMM_WORLD,ierr)
 !call MPI_TYPE_FREE(mat3,ierr)


 END subroutine broadcast_t
!----------------------------------------------------------------------------
 Subroutine mpimat2elint(mat2,nsize1,nsize2)
!----------------------------------------------------------------------------

    use MPI
    Implicit None
    
    integer, intent(in)::nsize1,nsize2
    integer ierr
    integer::typelist(2)
    integer ::imat(2),mat(2),start(2)
    integer ::mat2
!   allocate(struct%AS_ij(nsize,nsize))
    imat(1)=nsize1
    imat(2)=nsize2
    mat(1)=nsize1
    mat(2)=nsize2
    start(1)=0
    start(2)=0
    
    typelist(1)=MPI_INTEGER
    typelist(2)=MPI_INTEGER
!   write (*,*) nsize
    call MPI_TYPE_CREATE_SUBARRAY(2,imat,mat,start,MPI_ORDER_FORTRAN,MPI_INTEGER,mat2,ierr)
    call MPI_TYPE_COMMIT(mat2,ierr)


 End subroutine mpimat2elint
!----------------------------------------------------------------------------
 Subroutine mpimat2el(mat2,nsize1,nsize2)
!----------------------------------------------------------------------------

    use MPI
    implicit none
    
    integer, intent(in)::nsize1,nsize2
    integer ierr
    integer::typelist(2)
    integer ::imat(2),mat(2),start(2)
    integer ::mat2

!   allocate(struct%AS_ij(nsize,nsize))
    imat(1)=nsize1
    imat(2)=nsize2
    mat(1)=nsize1
    mat(2)=nsize2
    start(1)=0
    start(2)=0
    
    typelist(1)=MPI_DOUBLE_PRECISION
    typelist(2)=MPI_DOUBLE_PRECISION
!   write (*,*) nsize
    call MPI_TYPE_CREATE_SUBARRAY(2,imat,mat,start,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,mat2,ierr)
    call MPI_TYPE_COMMIT(mat2,ierr)


 End subroutine mpimat2el
!----------------------------------------------------------------------------
 Subroutine mpimat3elint(mat3,nsize1,nsize2,nsize3)
!----------------------------------------------------------------------------

    use MPI
    implicit none
    
    integer, intent(in)::nsize1,nsize2,nsize3
    integer :: ierr
    integer :: typelist(3)
    integer :: imat(3),mat(3),start(3)
    integer :: mat3

!   allocate(struct%AS_ij(nsize,nsize))
    imat(1)=nsize1
    imat(2)=nsize2
    imat(3)=nsize3
    mat(1)=nsize1
    mat(2)=nsize2
    mat(3)=nsize3
    start(1)=0
    start(2)=0
    start(3)=0
    
    typelist(1)=MPI_INTEGER
    typelist(2)=MPI_INTEGER
    typelist(3)=MPI_INTEGER
!   write (*,*) nsize
    call MPI_TYPE_CREATE_SUBARRAY(3,imat,mat,start,MPI_ORDER_FORTRAN,MPI_INTEGER,mat3,ierr)
    call MPI_TYPE_COMMIT(mat3,ierr)


 End subroutine mpimat3elint
!----------------------------------------------------------------------------
 Subroutine mpimat3el(mat3,nsize1,nsize2,nsize3)
!----------------------------------------------------------------------------

   use MPI
   Implicit None
   
   integer, intent(in)::nsize1,nsize2,nsize3
   integer ierr
   integer::typelist(3)
   integer ::imat(3),mat(3),start(3)
   integer ::mat3

!  allocate(struct%AS_ij(nsize,nsize))
   imat(1)=nsize1
   imat(2)=nsize2
   imat(3)=nsize3
   mat(1)=nsize1
   mat(2)=nsize2
   mat(3)=nsize3
   start(1)=0
   start(2)=0
   start(3)=0
   
   typelist(1)=MPI_DOUBLE_PRECISION
   typelist(2)=MPI_DOUBLE_PRECISION
   typelist(3)=MPI_DOUBLE_PRECISION
!  write (*,*) nsize
   call MPI_TYPE_CREATE_SUBARRAY(3,imat,mat,start,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,mat3,ierr)
   call MPI_TYPE_COMMIT(mat3,ierr)


 End subroutine mpimat3el
