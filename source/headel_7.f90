!--- TO DO ---
! 1. F2SD, M2DS global [rescpv]
!    ### spyros all forces are in global coordinates
!               the pitching moment is at present a scalar OK
! 2. nbod_el ?? vs nbod_ae. Is there any case that these should be different?
! 3. GetGeometry_ae for the rest aerodynamic bodies OK check with voutsi
!    ?? Forc_str(str0+nsp)%PosSection(:)
!       Forc_str(str0+nsp)%VelSection(:)
! OK
!    ?? Forc_all(nbod_ae)%BodyAxis(:) = BodyAx(:)
!       Forc_all(nbod_ae)%CircAxis(:) = CircAx(:)
!       Forc_all(nbod_ae)%OutPlAxis(:)= OutPlAx(:)
! 4. Flap motion [modify XGL] should be a separate subroutine, which could also be called for
!    the aerodynamic simulation in the same way... OK
! 5. Deside if the geometry is given with or without Hhub [if the hub is included then the
!    input of the aerodynamic and the aeroelastic simulations is the same]. Check only the
!    orientation.
!
!  DefineGeo  is moved to header [both header_GENUVP.f90 or header.f90] because is different
!
!----------------------------------------------------------------------------
!- Subroutine :GetGeometry
!
!  Called at every time step. 
!  Entering Get Geometry the following is known:
!           XGS(node_ae)%XGL
!  Defines: 
!           XGS(node_ae)%XG(3): Global position of every gnvp grid node
!           XGS(node_ae)%VelG(3): Global velocity of every gnvp grid node 
!                              (rigid body motions + elastic velocities)
!           Forc_all(nbod_ae)%BodyAxis(3): y-axis
!           Forc_all(nbod_ae)%CircAxis(3):  
!           Forc_all(nbod_ae)%OutPlAxis(3):
!           Forc_str(istr)%PosSection(3): 
!           Forc_str(istr)%VelSection(3):
!           ROTOR_center(3): 
!           ROTOR_axis(3):
!  
!----------------------------------------------------------------------------
   SUBROUTINE GetGeometry (IAXISRF,NTIME,TIME4,UREF4,TimStep)
!----------------------------------------------------------------------------

   Use Cbeam
   Use Conf_all
   Use Geom_all
   Use Grid
   Use Forces
   Use Move
   Use Wake
#ifdef HAVE_MPI
   Use MPI
#endif
   
   implicit none

   integer, intent(IN) :: IAXISRF,NTIME
   real(4), intent(IN) :: UREF4,TIME4,TimStep !=DT

   integer :: nbod_ae,nb,nbod_el,nbsub,nn_el,nsp,nch,nb_el,nel,ni,NDFPE, NEQPE

   real(8) :: SHAPE   (NEQPEM_el,NDFPEM_el)
   real(8) :: UT      (NDFPEM_el )
   real(8) :: UT1     (NDFPEM_el )
   real(8) :: UT2     (NDFPEM_el )
   real(8) :: UTT1    (NEQPEM_el )
   real(8) :: UTT0    (NEQPEM_el )
   real(8) :: AKSI0(3), DAKSI0(3), SS(3,6),Xnod(3),Vnod(3)
   real(8) :: HTA,HTAL,HTA0,HTA1,PHPX,PHPZ,ALLOC
   integer :: str,nodw,nod1,nod2,nodww
   real(4) :: XX(3),VV(3),BodyAx(3),CircAx(3),OutPlAx(3), RotorAx(3) 
   real(4) :: RotorCnt(3)
   real(4) :: vector(3),Pos_1(3),DPos_1(3),RotMat_1(3,3),DRotMat_1(3,3)  
!--- for mpi
   integer :: nimax
   real(4) :: xnod1,xnod2,xnod3  
   real(4) :: circ(3),bodax(3),outax(3),rotcenter(3),rotax(3)
   real(4), allocatable :: XG_loc(:,:), VelG_loc(:,:)
   integer :: my_rank,ierr,mat2                                                          


#ifdef HAVE_MPI
 call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)
#else
 my_rank=0; ierr=0;
#endif
!-----------------------------------------------------------------
!spyros ...Use XGL in order to define the chord ...  
!
!=======================
!------ Set XGRef ------
!=======================
!--- For all the bodies set the global Reference position on the GENUVP grid
         str     = 0
   do    nbod_ae = 1, NBod_all%NBODT
         nb      = nbod_ae
         nbod_el = Body_all(nbod_ae)%ElastBod

      if (nbod_el == 0) then
!=====================================================
!---------- A E R O D Y N A M I C -------------------
!=====================================================
!        caclulate XG, R_el, DR_el, A_el, DA_el, AT_el
         call MOVEBOD1(NTIME,TIME4,nbod_ae,Pos_1,DPos_1,RotMat_1,DRotMat_1)

         vector(1:3)=Body_all(nb)%ROTOR_axis0(:)
             Body_all(nb)%ROTOR_axis(:) = matmul(RotMat_1,vector)
         vector(1:3)=Body_all(nb)%ROTOR_center0(:)
         RotorCnt(:)=Body_all(nb)%ROTOR_center0(:)
             Body_all(nb)%ROTOR_center(:) = Pos_1+matmul(RotMat_1,RotorCnt)
         vector(1:3)=0.; vector(2)=1.
             Forc_all(nb)%BodyAxis(:) =  matmul(RotMat_1,BodyAx)
         vector(1:3)=0.; vector(1)=1.
             Forc_all(nb)%CircAxis(:) = matmul(RotMat_1,CircAx)
         vector(1:3)=0.; vector(3)=1.
             Forc_all(nb)%OutPlAxis(:) = matmul(RotMat_1,OutPlAx)

         do    nsp  = 1, Body_all(nbod_ae)%NCWB
            do nch  = 1, Body_all(nbod_ae)%NNBB
               ni   = Body_acu(nbod_ae-1)%NNTBACC   + (nsp-1)*Body_all(nbod_ae)%NNBB + nch 
               XX(:)=XGS(ni)%XGL0
               XGS(ni)%XGRef(:) = Move_all(nbod_ae)% R_el(:) + matmul(Move_all(nbod_ae)% A_el,XX)
            enddo !nch
         enddo !nsp
      else
!=====================================================
!---------- A E R O E L A S T I C -------------------
!=====================================================
       if(my_rank.eq.0)then
         do nsp     = 1, Body_all(nbod_ae)%NCWB
            nbsub   = Forc_elast(nsp,nbod_el)%NSUBe2a
            nel     = Forc_elast(nsp,nbod_el)%NELe2a
            HTA     = Forc_elast(nsp,nbod_el)%HTAe2a
            nb_el   = body(nbod_el)%NBODTGLB_el(nbsub)
            nn_el   = body(nbod_el)%NNODTGLB_el(nbsub)
      
            ALLOC   = subbody(nb_el)%ALENG_el  (nel)
            HTA0    = subbody(nb_el)%HTA_el    (nel)
            HTA1    = HTA0 + HTA
            HTAL    = HTA/ALLOC
            PHPX    = subbody(nb_el)%PHIX_el   (nel)
            PHPZ    = subbody(nb_el)%PHIZ_el   (nel)
            NEQPE   = subbody(nb_el)%NEQPE_el
            NDFPE   = subbody(nb_el)%NDFPE_el
      
            call SSHAPEFUNC15 (HTAL , ALLOC, PHPX, PHPZ, SHAPE, NDFPE, NEQPE)
            call LOCAL_UT_1   (nb_el, nel  , UT  , UT1 , UT2  )

            UTT0    =       matmul  ( SHAPE, UT  )
            UTT1    =       matmul  ( SHAPE, UT1 )
      
!---------- Set deformed grid and velocities
            do nch  = 1, Body_all(nbod_ae)%NNBB
               ni   = Body_acu(nbod_ae-1)%NNTBACC + (nsp-1)*Body_all(nbod_ae)%NNBB + nch 

               AKSI0     (1:3) = dble(XGS(ni)%XGL(1:3))
!old           AKSI0     (2)   = HTA1                       !qq!dble(XGL(2,ni))
!qq check HTA1, AKSI0(2)
               call Calc_SS ( SS, AKSI0 )
      
               Xnod         (1:3) =           transf_mat(nn_el)% R_el(1:3) &
                                     + matmul(transf_mat(nn_el)% A_el(1:3,1:3),matmul(SS(1:3,1:6),UTT0(IMDOF_el(1:6)) ) + AKSI0(1:3))
               XGS(ni)%XGRef(1:3) = sngl(Xnod(1:3))
            enddo ! nch(NNB)
         enddo ! nsp(NCW)
!--------------------------------------------------------------------------
       endif !my_rank
      endif !(nbod_el==0)
   enddo ! nbod_ae

!  XGL is the local/reference + flap we have XGL_fl and we then move
!qqXGS(:)%XGL_fl(1)=XGS(:)%XGL(1); XGS(:)%XGL_fl(2)=XGS(:)%XGL(2); XGS(:)%XGL_fl(3)=XGS(:)%XGL(3)
#ifdef HAVE_MPI
!------ send to other ranks
      do    nbod_ae = 1, NBod_all%NBODT
         do nsp     = 1, Body_all(nbod_ae)%NCWB
            do nch  = 1, Body_all(nbod_ae)%NNBB
               ni   = Body_acu(nbod_ae-1)%NNTBACC + (nsp-1)*Body_all(nbod_ae)%NNBB + nch 
             if(my_rank.eq.0)then
              xnod1 =  XGS(ni)%XGRef(1)   
              xnod2 =  XGS(ni)%XGRef(2)   
              xnod3 =  XGS(ni)%XGRef(3)   
             endif !my_rank
             call MPI_BCAST(xnod1,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
             call MPI_BCAST(xnod2,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
             call MPI_BCAST(xnod3,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
              XGS(ni)%XGRef(1) = xnod1   
              XGS(ni)%XGRef(2) = xnod2   
              XGS(ni)%XGRef(3) = xnod3   
            enddo
         enddo
      enddo
#endif

!=======================
!---- Set XG,VelG ----
!=======================
   XGS(:)%XGL_fl (1)=XGS(:)%XGL0(1); XGS(:)%XGL_fl (2)=XGS(:)%XGL0(2); XGS(:)%XGL_fl (3)=XGS(:)%XGL0(3)
   XGS(:)%VelL_fl(1)=0.            ; XGS(:)%VelL_fl(2)=0.            ; XGS(:)%VelL_fl(3)=0.
!
!  open(1001,file='xglast.dat')
!  open(1002,file='strippalast.dat')
!
!--- For all the bodies set the global positions and velocities on the GENUVP grid
         str     = 0
   do    nbod_ae = 1, NBod_all%NBODT
         nbod_el = Body_all(nbod_ae)%ElastBod
         if (Move_all(nbod_ae)%IYNFlap.ne.0) call Move_Flap(nbod_ae,TIME4) !modifies XGL_fl, VelL_fl based on flap motion

      if (nbod_el == 0) then
!=====================================================
!---------- A E R O D Y N A M I C -------------------
!=====================================================
!        R_el, DR_el, A_el, DA_el, AT_el have been calculated !!!
         do    nsp  = 1, Body_all(nbod_ae)%NCWB
            do nch  = 1, Body_all(nbod_ae)%NNBB
               ni   = Body_acu(nbod_ae-1)%NNTBACC   + (nsp-1)*Body_all(nbod_ae)%NNBB + nch 
               XX(:)=XGS(ni)%XGL_fl (:)  ! after flap deflection
               VV(:)=XGS(ni)%VelL_fl(:)  ! flap velocity
               XGS(ni)%XG(:)   = Move_all(nbod_ae)% R_el(:) + matmul(Move_all(nbod_ae)% A_el,XX)
               XGS(ni)%VelG(:) = Move_all(nbod_ae)%DR_el(:) + matmul(Move_all(nbod_ae)%DA_el,XX) &
                                                            + matmul(Move_all(nbod_ae)% A_el,VV)
            enddo !nch
!           (str)%PosSection, %VelSection position and velocity along
!                             the beam axis (SectionLoads)
            XX(1)=0.; XX(3)=0.
            str  = str + 1
            Forc_str(str)%PosSection(:) = Move_all(nbod_ae)% R_el(:) + matmul(Move_all(nbod_ae)% A_el,XX)
            Forc_str(str)%VelSection(:) = Move_all(nbod_ae)%DR_el(:) + matmul(Move_all(nbod_ae)%DA_el,XX)
         enddo !nsp
      else
!=====================================================
!---------- A E R O E L A S T I C -------------------
!=====================================================
       if(my_rank.eq.0)then
        if (NBODBT_el > NBLADE_el) then
!--------- last subbody node of the shaft = hub center
         nbsub = body(NSHAFT_el)%NBODSUB_el
         nn_el = body(NSHAFT_el)%NNODTGLB_el(nbsub+1)
         
!--------- Rotor_center coincides with the instantaneous hub position
         Body_all(nbod_ae)%ROTOR_center(1:3) = sngl(transf_mat(nn_el)%R_el(1:3)  )
         
!--------- Rotor_axis should point along ncw [VAWT]
!--------- Rotor_axis should point along the wind direction [HAWT]
         Body_all(nbod_ae)%ROTOR_axis  (1:3) = -sngl(transf_mat(nn_el)%A_el(1:3,2))
        else
         Body_all(nbod_ae)%ROTOR_center(1:3) = sngl(HUBpos_el(1:3) )
         Body_all(nbod_ae)%ROTOR_axis  (1:3) = 0.
!qq check if it the opposite
         if (IAPPL_el == 0)                       &
         Body_all(nbod_ae)%ROTOR_axis  (  1) = 1.
         if (IAPPL_el == 1 .or. IAPPL_el == 2)    &
         Body_all(nbod_ae)%ROTOR_axis  (  3) =-1.
        endif !NBODBT_el>NBLADE_el
         Forc_all(nbod_ae)%CircAxis    (1:3) = sngl(-ATP_el(nbod_el,1,1:3)) !Here the cone is considered, while in purely aerodynamic case usually not [only if the  cone is defined in the input grid]
         Forc_all(nbod_ae)%BodyAxis    (1:3) = sngl( ATP_el(nbod_el,2,1:3))
         Forc_all(nbod_ae)%OutPlAxis   (1:3) = sngl( ATP_el(nbod_el,3,1:3))
!w       write(*,*)'rotor center',Body_all(nbod_ae)%ROTOR_center(1:3)
!w       write(*,*)'rotor axis  ',Body_all(nbod_ae)%ROTOR_axis  (1:3)
!w       write(*,*)'nb',nbod_el
!w       write(*,*)'bodyaxis ', sngl(ATP_el(nbod_el,2,1:3))
!w       write(*,*)'circaxis ',-sngl(ATP_el(nbod_el,1,1:3))
!w       write(*,*)'outPlaxis', sngl(ATP_el(nbod_el,3,1:3))
       endif !my_rank
#    ifdef HAVE_MPI
      if(my_rank.eq.0)then
         circ(1:3)     = Forc_all(nbod_ae)%CircAxis   (1:3)
         bodax(1:3)    = Forc_all(nbod_ae)%BodyAxis   (1:3)
         outax(1:3)    = Forc_all(nbod_ae)%OutPlAxis  (1:3)
         rotcenter(1:3)= Body_all(nbod_ae)%ROTOR_axis (1:3)
         rotax(1:3)    = Body_all(nbod_ae)%ROTOR_axis (1:3)
      endif !my_rank
      
      call MPI_BCAST(circ     ,3,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(bodax    ,3,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(outax    ,3,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rotcenter,3,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rotax    ,3,MPI_REAL,0,MPI_COMM_WORLD,ierr)
 
      Forc_all(nbod_ae)%CircAxis   (1:3) = circ(1:3)
      Forc_all(nbod_ae)%BodyAxis   (1:3) = bodax(1:3)
      Forc_all(nbod_ae)%OutPlAxis  (1:3) = outax(1:3)
      Body_all(nbod_ae)%ROTOR_axis (1:3) = rotcenter(1:3)
      Body_all(nbod_ae)%ROTOR_axis (1:3) = rotax(1:3)
#    endif

      if(my_rank.eq.0)then
         do nsp     = 1, Body_all(nbod_ae)%NCWB
            nbsub   = Forc_elast(nsp,nbod_el)%NSUBe2a
            nel     = Forc_elast(nsp,nbod_el)%NELe2a
            HTA     = Forc_elast(nsp,nbod_el)%HTAe2a
            nb_el   = body(nbod_el)%NBODTGLB_el(nbsub)
            nn_el   = body(nbod_el)%NNODTGLB_el(nbsub)
      
            ALLOC   = subbody(nb_el)%ALENG_el  (nel)
            HTA0    = subbody(nb_el)%HTA_el    (nel)
            HTA1    = HTA0 + HTA
            HTAL    = HTA/ALLOC
            PHPX    = subbody(nb_el)%PHIX_el   (nel)
            PHPZ    = subbody(nb_el)%PHIZ_el   (nel)
            NEQPE   = subbody(nb_el)%NEQPE_el
            NDFPE   = subbody(nb_el)%NDFPE_el
      
            call SSHAPEFUNC15 (HTAL , ALLOC, PHPX, PHPZ, SHAPE, NDFPE, NEQPE)
            call LOCAL_UT_1   (nb_el, nel  , UT  , UT1 , UT2  )

            UTT0    =       matmul  ( SHAPE, UT  )
            UTT1    =       matmul  ( SHAPE, UT1 )
      
!---------- Set deformed grid and velocities
            do nch  = 1, Body_all(nbod_ae)%NNBB
               ni   = Body_acu(nbod_ae-1)%NNTBACC + (nsp-1)*Body_all(nbod_ae)%NNBB + nch 
!qq no XGL_fl here qq Needs special treatment in case of subbodies....jim
!old           AKSI0     (1:3) = dble(XGS(ni)%XGL(1:3))
!old           AKSI0     (2)   = HTA1                       !qq!dble(XGL(2,ni))
!qq check HTA1, AKSI0(2)
               AKSI0     (1:3) = dble(matmul ( ATsbAsb1     (nbsub,1:3,1:3),  XGS(ni)%XGL_fl (1:3) ) &
                                             - ATsbRsbmRsb1 (nbsub,    1:3)                         )    ! local undeformed geometry, including modifications by flap motion wrt each sub-body local c.s.
               DAKSI0    (1:3) = dble(matmul ( ATsbAsb1     (nbsub,1:3,1:3),  XGS(ni)%VelL_fl(1:3) ))    ! flap velocity wrt each sub-body local c.s.
      
               call Calc_SS ( SS, AKSI0 )
      
               Xnod(1:3) = transf_mat(nn_el)% R_el(1:3) + matmul ( transf_mat(nn_el)% A_el(1:3,1:3), matmul( SS(1:3,1:6), UTT0(IMDOF_el(1:6)) ) +  AKSI0(1:3) )
               Vnod(1:3) = transf_mat(nn_el)%DR_el(1:3) + matmul ( transf_mat(nn_el)%DA_el(1:3,1:3), matmul( SS(1:3,1:6), UTT0(IMDOF_el(1:6)) ) +  AKSI0(1:3) ) &
                                                        + matmul ( transf_mat(nn_el)% A_el(1:3,1:3), matmul( SS(1:3,1:6), UTT1(IMDOF_el(1:6)) ) + DAKSI0(1:3) )
               XGS(ni)%XG(1:3)   = sngl(Xnod(1:3))
               XGS(ni)%VelG(1:3) = sngl(Vnod(1:3))

!              if (nch==1.and.nsp==20) write(666,'(a,6e25.12)')'x,u',Xnod(1:3),Vnod(1:3)
!              if (ni<1989)then                                                  
!              write(2000,'(i8,6e25.12)')ni,XGS(ni)%XG(1:3),XGS(ni)%VG(1:3)
!              write(2001,'(i8,3e25.12)')ni,XGS(ni)%VG(1:3)
!              endif
!
!              if (ni==1580) write(*,*)'bb1', Xnod(1:3)
!              if (ni==1580) write(*,*)'bb2', Vnod(1:3)

!              write(1001,'(3f15.5)') Xnod(1:3)
            enddo ! nch(NNB)
!------------- Set elastic axis points at which loads are communicated
!old           AKSI0     (1:3) = 0.d0;
!old           AKSI0     (2)   = HTA1                       !qq!dble(XGL(2,ni))
               AKSI0     (1)   = 0.d0; !AKSI0(2) known from previous loop
               AKSI0     (3)   = 0.d0;
      
               call Calc_SS ( SS, AKSI0 )
      
               Xnod(1:3) = transf_mat(nn_el)% R_el(1:3) + matmul ( transf_mat(nn_el)% A_el(1:3,1:3), matmul( SS(1:3,1:6), UTT0(IMDOF_el(1:6)) ) + AKSI0(1:3) )
               Vnod(1:3) = transf_mat(nn_el)%DR_el(1:3) + matmul ( transf_mat(nn_el)%DA_el(1:3,1:3), matmul( SS(1:3,1:6), UTT0(IMDOF_el(1:6)) ) + AKSI0(1:3) ) &
                                                        + matmul ( transf_mat(nn_el)% A_el(1:3,1:3), matmul( SS(1:3,1:6), UTT1(IMDOF_el(1:6)) )              )
               Forc_elast(nsp,nbod_el)%STRIPPA(1:3) = sngl(Xnod(1:3))

                        str                  = str + 1
               Forc_str(str)%PosSection(1:3) = sngl(Xnod(1:3))
               Forc_str(str)%VelSection(1:3) = sngl(Vnod(1:3))
!              write(1002,'(3f15.5)') Xnod(1:3); write(1001,*); write(1002,*)
         enddo ! nsp(NCW)
!              write(1001,*); write(1001,*); write(1002,*); write(1002,*)

!              !if (ni==1580) write(*,*)'bb1', Xnod(1:3)
!              !if (ni==1580) write(*,*)'bb2', Vnod(1:3)
!              !if (ni<1989)then
!              !write(3000,'(i8,6e25.12)')ni,XGS(ni)%XG(1:3),XGS(ni)%VG(1:3)
!              !write(3001,'(i8,3e25.12)')ni,XGS(ni)%VG(1:3)
!              !endif
!--------------------------------------------------------------------------
      endif !(nbod_el == 0)
      endif!my_rank
   enddo ! nbod_ae
!              close(1001); close(1002)
!              stop

#  ifdef HAVE_MPI
      nimax=NBod_all%NNT
      call MPI_BCAST(nimax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      allocate ( XG_loc(3,nimax),VelG_loc(3,nimax) )
    if(my_rank.eq.0)then
      XG_loc  (1,:) = XGS(:)%XG(1)
      XG_loc  (2,:) = XGS(:)%XG(2)
      XG_loc  (3,:) = XGS(:)%XG(3)
      VelG_loc(1,:) = XGS(:)%VelG(1)
      VelG_loc(2,:) = XGS(:)%VelG(2)
      VelG_loc(3,:) = XGS(:)%VelG(3)
    endif !my_rank

      call mpimat2(mat2,3,nimax,3,nimax,0,0)
      call MPI_BCAST(XG_loc  ,1,mat2,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(VelG_loc,1,mat2,0,MPI_COMM_WORLD,ierr)
      call MPI_TYPE_FREE(mat2,ierr)


      XGS(:)%XG(1)   = XG_loc  (1,:)
      XGS(:)%XG(2)   = XG_loc  (2,:)  
      XGS(:)%XG(3)   = XG_loc  (3,:)  
      XGS(:)%VelG(1) = VelG_loc(1,:)                                             
      XGS(:)%VelG(2) = VelG_loc(2,:)                                             
      XGS(:)%VelG(3) = VelG_loc(3,:)                                             
                       
      deallocate (XG_loc,VelG_loc)
#  endif

!qq jim what is this?
   do nbod_ae=1,NBod_all%NBODT
   if(Body_all(nbod_ae)%YNLiftB.eq.0) cycle
      if (NTIME.eq.0) then
         do nodw = wake_acc1(nbod_ae-1)%NUEMACC+1, wake_acc1(nbod_ae)%NUEMACC
            nod1=wake_oth(nodw)%NODEMIT(1)
            nod2=wake_oth(nodw)%NODEMIT(2)
            XGWS(nodw)%XWAKE(:) = 0.5*(XGS(nod1)%XG(:)+XGS(nod2)%XG(:))
            nodww = nodw + NBod_all%NUEMNO
            XGWS(nodww)%XWAKE(:) = XGWS(nodw)%XWAKE(:)
            XGWS(nodww)%XWAKE(IAXISRF) = XGWS(nodw)%XWAKE(IAXISRF) + UREF4*TimStep
         enddo ! nodw
      else
         do nodw = wake_acc1(nbod_ae-1)%NUEMACC+1, wake_acc1(nbod_ae)%NUEMACC
            nodww = nodw + NBod_all%NUEMNO
            XGWS(nodww)%XWAKE(:) = XGWS(nodw)%XWAKE(:)
            nod1=wake_oth(nodw)%NODEMIT(1)
            nod2=wake_oth(nodw)%NODEMIT(2)
            XGWS(nodw)%XWAKE(:) = 0.5*(XGS(nod1)%XG(:)+XGS(nod2)%XG(:))
!>                         XWAKE (1,NOD) = 0.5*(XG(1,NODEMIT(NOD,1,1))+
!>         .                                    XG(1,NODEMIT(NOD,2,1))  )
!>                         XWAKE (2,NOD) = 0.5*(XG(2,NODEMIT(NOD,1,1))+
!>         .                                    XG(2,NODEMIT(NOD,2,1))  )
!>                         XWAKE (3,NOD) = 0.5*(XG(3,NODEMIT(NOD,1,1))+
!>         .                                    XG(3,NODEMIT(NOD,2,1))  )
         enddo ! nodw
      endif
   enddo


   END Subroutine GetGeometry
!------------------------------------------------------------------------------------------
!
!  Subroutine :Aero2Elast
!
!  Define distributed aerodynamic loads for aeroelastic simulations at the grid nodes [NCW].
!  Only called from hGAST and not needed in purely aerodynamic cases.
!
!------------------------------------------------------------------------------------------
   Subroutine Aero2Elast  !GetLoads
!------------------------------------------------------------------------------------------

   use Conf_all
   use Geom_all
   use Forces
   use Run_def !velhub_gn_dp
   use Craft   !HUB_VEL

   implicit none

   integer :: I, nbod_ae, nbod_el, nsp
   real(4) :: R(3), F(3), M(3)


   HUB_VEL (1:3)  = velhub_gn_dp(1:3)

   do     nbod_ae = 1, NBod_all%NBODT
          nbod_el = Body_all(nbod_ae)%ElastBod
      if (nbod_el == 0) cycle
      do nsp      = 1, Body_all(nbod_ae)%NCWB1
         F(1:3)   = 0.;
         M(1:3)   = 0.;
         I        = Body_acu(nbod_ae-1)%NSTRIPACC + nsp
         R(1:3)   = Forc_str(I)%XGFC4(1:3) - (Forc_elast(nsp+0,nbod_el)%STRIPPA(1:3)+Forc_elast(nsp+1,nbod_el)%STRIPPA(1:3))/2.
if (nsp.gt.Body_all(nbod_ae)%ROOTcut) then
        !F(1:3)   = Forc_str(I)%FSTRGN  (1:3) / Forc_str(I)%DSPAN
        !F(1:3)   = Forc_str(I)%FSTR2D  (1:3) / Forc_str(I)%DSPAN
         F(1:3)   = Forc_str(I)%FSTRDS2D(1:3) / Forc_str(I)%DSPAN
else
         F(1:3)   = Forc_str(I)%FSTR2D  (1:3) / Forc_str(I)%DSPAN
endif
      
         call EXTEPR_gn ( R(1), R(2), R(3),&
                          F(1), F(2), F(3),&
                          M(1), M(2), M(3)  )
if (nsp.gt.Body_all(nbod_ae)%ROOTcut) then
        !M(1:3)   = M(1:3) + Forc_str(I)%MSTRGN  (1:3) / Forc_str(I)%DSPAN
        !M(1:3)   = M(1:3) + Forc_str(I)%MSTR2D  (1:3) / Forc_str(I)%DSPAN
         M(1:3)   = M(1:3) + Forc_str(I)%MSTRDS2D(1:3) / Forc_str(I)%DSPAN
else
         M(1:3)   = M(1:3) + Forc_str(I)%MSTR2D  (1:3) / Forc_str(I)%DSPAN
endif
    
         Forc_elast(nsp,nbod_el)%FSTRIPG(1:3) = dble(F(1:3))
         Forc_elast(nsp,nbod_el)%MSTRIPG(1:3) = dble(M(1:3))
      enddo !nsp
   enddo !nbod_ae


   END Subroutine Aero2Elast 
!----------------------------------------------------------------------------
   Subroutine write_gn_load(TIMED,NBLADE_el)
!----------------------------------------------------------------------------
   use Conf_all
   use Geom_all
   use Forces

   implicit none

   integer, intent(IN)  :: NBLADE_el
   real(8), intent(IN)  :: TIMED  

   integer      :: nbod, nsp, j
   character    :: CNUM1(2), CNUM2(2)
   character*80 :: outfil

!qq
   return
   do nbod = 1, NBLADE_el
      call INT_2_CHAR ( 2, CNUM1, nbod )
   do nsp  = 1, Body_all(nbod)%NCWB - 1
      call INT_2_CHAR ( 2, CNUM2, nsp  )

      outfil = 'aeroloa'//CNUM1(1)//CNUM1(2)//'_'//CNUM2(1)//CNUM2(2)//'.dat'
      open (1,file=outfil, access='append')! ,form='UNFORMATTED')
       write (1,100)      TIMED                                 ,& ! 1
                         (Forc_elast(nsp,nbod)%FSTRIPG(j),j=1,3),& ! 2-4
                         (Forc_elast(nsp,nbod)%MSTRIPG(j),j=1,3),& ! 5-7
                          Forc_elast(nsp,nbod)%RSTRIP              ! 8

      close (1)
      
   enddo !nsp
   enddo !nbod

 100  format (150e24.8)


   END Subroutine write_gn_load
!----------------------------------------------------------------------------
!- Subroutine :LOCAL_FORC0_gn_ae
!
!  Set the local distributed aerodynamic force on the blade beam.
!  Called from math [FEM] subroutines in order to consider the virtual work
!  of the aerodynamic forces.
!----------------------------------------------------------------------------
 Subroutine LOCAL_FORC0_gn_ae ( UTT0, nbod, nbsub, nel , k, HLOC, F_ae, XC, ZC )
!Subroutine LOCAL_FORC0_gn_ae (       nbod, nbsub,          HLOC, F_ae )
!----------------------------------------------------------------------------

   use Cbeam
   use Geom_all !Body_all
   use Forces

   implicit none

   integer, intent(in ) :: nel, k, nbod, nbsub   !DUMMY
   real(8), intent(in ) :: UTT0(NEQPEM_el), HLOC !DUMMY
   real(8), intent(out) :: F_ae(6), XC, ZC

   integer :: nn_el, ib_el, ibsub, ncw
   real(8) :: Sloc, F(3), AM(3), ATB(3,3)

   real(8), allocatable :: FX(:), FY(:), FZ(:), MX(:), MY(:), MZ(:), S_pa(:)


!--- in this approach the loads have been already transfered to the pitch axis
   F_ae(1:6) = 0.d0; XC=0.d0; ZC=0d0;

   if (Body_all(nbod)%ElastBod == 0) then
      write(*,*)'elast body is zero in LOCAL_FORC0_gn_ae for nbod',nbod
      write(*,*)'zero aerodyamic force was set'
      return
   endif
   if (nbod/=Body_all(nbod)%ElastBod) then
      write(*,*)'error in elastic and aerodynamic bodies correspondance'
      stop
   endif

   ncw = Body_all(nbod)%NCWB - 1

   Allocate (FX(ncw),FY(ncw),FZ(ncw),MX(ncw),MY(ncw),MZ(ncw),S_pa(ncw))

   nn_el        = body(nbod)%NNODTGLB_el(nbsub)
   ATB(1:3,1:3) = transf_mat(nn_el)%AT_el(1:3,1:3)         !from global to local

!--- Calculate Sloc of blade
      Sloc  = HLOC
   do ibsub = 1, nbsub-1
      ib_el = body(nbod)%NBODTGLB_el(ibsub)
      Sloc  = Sloc + subbody(ib_el)%ALENGB_el
   enddo

   S_pa(1:ncw) = Forc_elast(1:ncw,nbod)%RSTRIP

   if (Sloc < S_pa (1)) then
      return
   endif

   FX (1:ncw) = Forc_elast(1:ncw,nbod)%FSTRIPG(1);
   FY (1:ncw) = Forc_elast(1:ncw,nbod)%FSTRIPG(2);
   FZ (1:ncw) = Forc_elast(1:ncw,nbod)%FSTRIPG(3);
   MX (1:ncw) = Forc_elast(1:ncw,nbod)%MSTRIPG(1);
   MY (1:ncw) = Forc_elast(1:ncw,nbod)%MSTRIPG(2);
   MZ (1:ncw) = Forc_elast(1:ncw,nbod)%MSTRIPG(3);

!!      LIN_INT (  x  ,  y    , xi  , yi, n1 , nnm )
   call LIN_INT ( Sloc,  F (1), S_pa, FX, ncw, ncw )
   call LIN_INT ( Sloc,  F (2), S_pa, FY, ncw, ncw )
   call LIN_INT ( Sloc,  F (3), S_pa, FZ, ncw, ncw )
   call LIN_INT ( Sloc, AM (1), S_pa, MX, ncw, ncw )
   call LIN_INT ( Sloc, AM (2), S_pa, MY, ncw, ncw )
   call LIN_INT ( Sloc, AM (3), S_pa, MZ, ncw, ncw )

   F         = matmul ( ATB, F  )
   AM        = matmul ( ATB, AM )
   F_ae(1:3) = F  (1:3)
   F_ae(4:6) = AM (1:3)

   deallocate (FX, FY, FZ, MX, MY, MZ, S_pa)

   END Subroutine LOCAL_FORC0_gn_ae



!----------------------------------------------------------------------------
!
!  Subroutine :DefineGeo
!
!  Defines the initial local geometry [grid] of aerodynamic body XGL.
!  Only called once from the initia.
!
!----------------------------------------------------------------------------
   SUBROUTINE DefineGeo(nbod_ae,nbod_el,File_geo,NNBB,NCWB,NNT)
!----------------------------------------------------------------------------

   use Cbeam
   use Grid
   use Forces
   use Conf_all
   use Geom_all

   implicit none

   integer, intent(IN)  :: nbod_ae,nbod_el,NNBB,NCWB,NNT
   character (LEN=12)   :: File_geo

   integer :: nsp,ibsub,nch,nod,j ,nn0_el, nn_el, NCWmax
   real(8) :: A_sb1(3,3),R_sb1(3),AT_sb(3,3),R_sb(3),R_read(3)
   integer :: nbsub, nb_el, ib_el
   real(8) :: L
   character (LEN=12) :: File_strip
   integer :: nbod1_ae,nbod1_el, I
   real(8) :: AKSI0(3),Xnod(3), AT_sb1(3,3), Rref(3) !qq new for prebend and XGL0

!qq
   if (nbod_ae==1) then
      NCWmax=NCWB
      do nbod1_ae   = 1, NBod_all%NBODT
         if (NCWmax.lt.Body_all(nbod1_ae)%NCWB) NCWmax=Body_all(Nbod1_ae)%NCWB
      enddo
!qq   allocate (Forc_elast(NCWmax,NBod_all%NBODT))    !  maxNWCB, NBLADE_el
      allocate (Forc_elast(NCWmax,NBLADE_el     ))    !  maxNWCB, NBLADE_el
!qq   do nbod1_el   = 1, NBod_all%NBODT
      do nbod1_el   = 1, NBLADE_el
         do I       = 1, NCWmax
            Forc_elast(I,nbod1_el)%FSTRIPG(1:3) = 0.d0;
            Forc_elast(I,nbod1_el)%MSTRIPG(1:3) = 0.d0;
         enddo
      enddo

!------- for flap
                              ibsub = body(1)%NBODSUB_el
      allocate (ATsbAsb1     (ibsub,3,3))
      allocate (ATsbRsbmRsb1 (ibsub,  3))

         nbod1_el        = 1
         nn0_el          = body(nbod1_el)%NNODTGLB_el (1)
          A_sb1(1:3,1:3) = transf_mat(nn0_el)%A_el(1:3,1:3)
          R_sb1(1:3    ) = transf_mat(nn0_el)%R_el(1:3    )
         AT_sb1(1:3,1:3) = transpose(A_sb1(1:3,1:3))
      do ibsub           = 1, body(1)%NBODSUB_el
         nn_el           = body(nbod_el)%NNODTGLB_el (ibsub)
         AT_sb(1:3,1:3)  = transf_mat(nn_el)%AT_el(1:3,1:3)
         R_sb (1:3    )  = transf_mat(nn_el)%R_el (1:3    )

         ATsbAsb1     (ibsub,1:3,1:3) = sngl( matmul ( AT_sb(1:3,1:3), A_sb1(1:3,1:3)      ) )
         ATsbRsbmRsb1 (ibsub,    1:3) = sngl( matmul ( AT_sb(1:3,1:3), R_sb(1:3)-R_sb1(1:3)) )
      enddo !ibsub
   endif


   if (nbod_el==0) then
!=====================================================
!---------- A E R O D Y N A M I C -------------------
!=====================================================
     open (70,file=File_geo)
     read (70,*)  !strip_gn.inp
      do nsp = 1, NCWB
      do nch = 1, NNBB
         nod = NNT+(nsp-1)*NNBB + nch 
         read(70,*)  (R_read(j),j=1,3)
         XGS(nod)%XGL0(1:3) = sngl(R_read(1:3))
      enddo !nch
         read(70,*)
      enddo !nsp
     close (70)

   else !if (nbod_el /= 0) then
!=====================================================
!---------- A E R O E L A S T I C -------------------
!=====================================================
     open (70,file=File_geo)
     read (70,*) File_strip !strip_gn.inp
          open (71,file=File_strip)
           do nsp=1,NCWB
              read(71,*)                         &
               Forc_elast(nsp,nbod_el)%NSUBe2a  ,&
               Forc_elast(nsp,nbod_el)%S_elastic
           enddo !nsp
          close (71)

     call ELST2AERinit(nbod_ae) !set for each strip: NELe2a_gn, HTAe2a_gn
                                !qq could also find NSUBe2a, S_elastic!

!--- Read the blade geometry 
!---- LOCAL COORDINATES wrt blade root at the local blade c.s. of hGAST.
!---- The total length of the blade geometry is without the Hhub distance.
!---- R_input = AT_sb1 . (Rg_nod-Rsb1)
!---- NNBxNCW NNB: starts from LE to TE
!----         NCW: starts from root to tip

   !!    nn0_el          = body(nbod_el)%NNODTGLB_el (1)
   !!     A_sb1(1:3,1:3) = transf_mat(nn0_el)%A_el(1:3,1:3)
   !!     R_sb1(1:3    ) = transf_mat(nn0_el)%R_el(1:3    )
   !!    AT_sb1(1:3,1:3) = transpose(A_sb1(1:3,1:3))
      do nsp             = 1, NCWB
         ibsub           = Forc_elast(nsp,nbod_el)%NSUBe2a
   !!    nn_el           = body(nbod_el)%NNODTGLB_el (ibsub)
   !!    AT_sb(1:3,1:3)  = transf_mat(nn_el)%AT_el(1:3,1:3)
   !!    R_sb (1:3    )  = transf_mat(nn_el)%R_el (1:3    )
      do nch             = 1, NNBB
         nod             = NNT+(nsp-1)*NNBB + nch  !qq NNT??
         read(70,*)        (R_read(j),j=1,3)

!qq      R_read(2)       =  R_read(2) - Hhub !geometry given wrt Rhub
!qq      R_read(2)       =  R_read(2)        !geometry given wrt Rroot

!--------- XGL_gn       = AT_sb*A_sb1*(R_read)-AT_sb*(R_sb-R_sb1)

!qq------- valid for straight blade or if the prebend blade geometry has been given
         XGS(nod)%XGL0(1:3) = sngl(R_read(1:3))
         XGS(nod)%XGL (1:3) = matmul ( ATsbAsb1     (ibsub,1:3,1:3),  sngl(R_read(1:3)) ) & !!sngl( matmul ( matmul(AT_sb(1:3,1:3), A_sb1(1:3,1:3)),  R_read(1:3) ) &
                                     - ATsbRsbmRsb1 (ibsub,    1:3)                         !!     -matmul (        AT_sb(1:3,1:3), R_sb(1:3)-R_sb1(1:3)          )  )

!qq------- transform the straight blade into prebent [valid only if all nnb points of each section have the same y coordinate]
!!!      XGS(nod)%XGL (1:3) = sngl( R_read(1:3) )
!!!      XGS(nod)%XGL (  2) = sngl( Forc_elast(nsp,nbod_el)%S_elastic)
!!!      AKSI0        (1:3) = dble(XGS(nod)%XGL(1:3))
!!!      Rref         (1:3) = R_sb1(1:3)  !HUBpos_el(1:3)   !reference for the calculation of XGL0 (needed for RSPAN [ClCdCm] and flap moton
!!!      Xnod         (1:3) = matmul (AT_sb1(1:3,1:3), transf_mat(nn_el)% R_el(1:3) - Rref(1:3) + matmul ( transf_mat(nn_el)% A_el(1:3,1:3),  AKSI0(1:3) ) )
!!!      XGS(nod)%XGL0(1:3) = sngl( Xnod(1:3) )

!        write(*,'(a,9f12.3)') 'c', HTA1, Forc_elast(nsp,nbod_el)%S_elastic

         write(*,'(a,3i5,9f12.3)') 'AerogeoX', nbod_el, nsp,nch, R_read(1), XGS(nod)%XGL(1), XGS(nod)%XGL0(1)
         write(*,'(a,3i5,9f12.3)') 'AerogeoY', nbod_el, nsp,nch, R_read(2), XGS(nod)%XGL(2), XGS(nod)%XGL0(2), Forc_elast(nsp,nbod_el)%S_elastic
         write(*,'(a,3i5,9f12.3)') 'AerogepZ', nbod_el, nsp,nch, R_read(3), XGS(nod)%XGL(3), XGS(nod)%XGL0(3)
      enddo !nch
         read(70,*)
      enddo !nsp
     close (70)

!==========================================================
!--------- Set RSTRIP
     do nsp      = 1, NCWB - 1
        ibsub    = Forc_elast(nsp,nbod_el)%NSUBe2a
        ib_el    = body(nbod_el)%NBODTGLB_el (ibsub)
           L     = (Forc_elast(nsp  ,nbod_el)%S_elastic + &
                    Forc_elast(nsp+1,nbod_el)%S_elastic     )/2.d0
        do nbsub = 1, ibsub-1
           nb_el = body(nbod_el)%NBODTGLB_el (nbsub)
           L     = L + subbody(nb_el)%ALENGB_el
        enddo
        Forc_elast(nsp,nbod_el)%RSTRIP = L
        write(*,'(a,2i5,f15.5)')'RSTRIP_gn',nbod_el, nsp, L
     enddo !nsp
!==========================================================
!--------- Set RSTRIP for Aero2Elast0,1
!    do nsp      = 1, NCWB
!       ibsub    = Forc_elast(nsp,nbod_el)%NSUBe2a
!       ib_el    = body(nbod_el)%NBODTGLB_el (ibsub)
!          L     = Forc_elast(nsp,nbod_el)%S_elastic
!--- set RSTRIP, if S_elastic <0 or S_elastic>ALENGB_el do nothing set RSTRIP==1st or last distance
!          L     = min(L,subbody(ib_el)%ALENGB_el)
!          L     = max(L,0.d0                    )
!       do nbsub = 1, ibsub-1
!          nb_el = body(nbod_el)%NBODTGLB_el (nbsub)
!          L     = L + subbody(nb_el)%ALENGB_el
!       enddo
!       Forc_elast(nsp,nbod_el)%RSTRIP = L
!       write(*,'(a,2i5,f15.5)')'RSTRIP_gn',nbod_el, nsp, L
!    enddo !nsp
!==========================================================
!---- Define the spanwise indices at which the aerodynamic parts started [for VAWT]
!!    do nsp      = 1, NCWB
!!       ibsub    = Forc_elast(nsp,nbod_el)%NSUBe2a
!!       ib_el    = body(nbod_el)%NBODTGLB_el (ibsub)
!!          L     = Forc_elast(nsp,nbod_el)%S_elastic
!!-------- set RSTRIP_gn, if S_elastic <0 or S_elastic>ALENGB_el do nothing set RSTRIP==1st or last distance
!!          if (L<0.d0                    ) nsp_ext_gn(1)=nsp+1
!!          if (L>subbody(ib_el)%ALENGB_el) nsp_ext_gn(2)=nsp-1
!!          if (L>subbody(ib_el)%ALENGB_el) exit
!!    enddo !nsp
!!    write(*,*) 'external1-2',nsp_ext_gn(1),nsp_ext_gn(2)

   endif !nbod_el==0

   END Subroutine DefineGeo 



!----------------------------------------------------------------------------
!
!  Subroutine :ELST2AERinit
!
!  Defines the strip correspondance between the structural and the
!  aerodynamic grid, stored in NELe2a, HTAe2a.
!  Needed for the calculation of the deformed grid and the velocities.
!  Only called once from the initia and only for aeroelastic cases.
!
!  Interpolations are based on s and not r...
!
!---------------------------------------------------------------------------
   SUBROUTINE ELST2AERinit(nbod_ae)
!---------------------------------------------------------------------------

   use Cbeam
   use Forces
   use Conf_all !NBod_all
   use Geom_all !Body_all, Body_acu, Body_res

   implicit none

   integer,intent(in) :: nbod_ae

   integer :: nbod_el, nsp, ibsub,ib_el,iel
   real(8) :: HLOC, HLOC1, HLOC2


   nbod_el  = Body_all(nbod_ae)%ElastBod
        write (*,*) 'nbod_el  @ELST2AERinit:' , nbod_el
   if (nbod_el == 0) return

   do nsp      = 1, Body_all(nbod_ae)%NCWB
      HLOC     = Forc_elast(nsp,nbod_el)%S_elastic                 !defines the local position of the "airfoil"(ncw not ncw1) wrt the elastic axis. local distance from sbod start
      ibsub    = Forc_elast(nsp,nbod_el)%NSUBe2a
      ib_el    = body(nbod_el)%NBODTGLB_el (ibsub)
!------ Search which element HLOC lies in
      do iel   = 1, subbody(ib_el)%NTEB_el
         HLOC1 =    subbody(ib_el)%HTA_el(iel  ) !based on s
         HLOC2 =    subbody(ib_el)%HTA_el(iel+1)

!------ define NELe2a and HTAe2a for outer parts [upper] [for VAWT]
!!       if ( (iel==1.).and.(HLOC < 0.d0) ) then
!!          Forc_elast(nsp,nbod_el)%NELe2a = iel
!!          Forc_elast(nsp,nbod_el)%HTAe2a = 0.d0
!!          write(*,'(a,2i5,f19.5)')'upper ext: ibsub,iel,HTAe2a_gn',ibsub, iel, Forc_elast(nsp,nbod_el)%HTAe2a
!!          goto 1
!!       endif

         if ( (HLOC >= HLOC1).and.(HLOC <= HLOC2) ) then
         !!!HTA          = (HLOC-HLOC1)/(HLOC2-HLOC1)*ALLOC
         !!!HTAe2a_gn (nspt) = HTA
            Forc_elast(nsp,nbod_el)%NELe2a = iel
            Forc_elast(nsp,nbod_el)%HTAe2a = HLOC - HLOC1
            write(*,'(a,2i5,f19.5)')'ibsub,iel,HTAe2a_gn',ibsub, iel, Forc_elast(nsp,nbod_el)%HTAe2a
            goto 1
         endif
      enddo !iel

!------ define NELe2a and HTAe2a for outer parts [lower] [for VAWT]
!!          Forc_elast(nsp,nbod_el)%NELe2a = subbody(ib_el)%NTEB_el
!!          Forc_elast(nsp,nbod_el)%HTAe2a = HLOC2 - HLOC1
!!          write(*,'(a,2i5,f19.5)')'lower ext: ibsub,iel,HTAe2a_gn',ibsub, iel, Forc_elast(nsp,nbod_el)%HTAe2a
        
!     write(10,*) 'HLOC not found', HLOC, HLOC1,HLOC2
!     write(* ,*) 'HLOC not found', HLOC, HLOC1,HLOC2
!     stop
 1 enddo !nsp

 END Subroutine ELST2AERinit


!----------------------------------------------------------------------
 Subroutine Getfile_wind (file_wind_out)
!----------------------------------------------------------------------
!
! communicate the file_wind
!
!----------------------------------------------------------------------

 Use Paths

   implicit none

   character*256 :: file_wind_out

   file_wind_out = file_wind

 END Subroutine Getfile_wind
