!---------------------------------------------------------------------!
! REMARKS:
! ========
! 1. Elements orientation is from sea bed to tower bottom [Morison]
!    so if nodes are read from sea bed to MSL, the 1st node number of
!    the element MUST be less than the last. [jacket.inp file]
! 2. The last body of the jacket must be vertical, so the deflections of the
!    upper (last) node match the 1st q's of the tower (bottom). This is
!    helpful for the reduction of the qs without the need of matrices.
! 3. communication between jacket and tower is assumed to be performed to the 
!    last body of the jacket (last node) with the 1st node of the tower
!----------------------------------------------------------------------
 Module Jacket
!----------------------------------------------------------------------
    integer,              save :: NNT_ja
    real(8), allocatable, save :: XN_ja     (:,:) ! (3,NNT_ja)
    integer, allocatable, save :: Nod_ja    (:,:) ! (2,NBODBTJA_el)
    integer, allocatable, save :: NTEB_ja   (:  ) ! (  NBODBTJA_el)
    integer, allocatable, save :: NBC_ja    (:  ) ! (  NBODBTJA_el)
    integer, allocatable, save :: ITYPE_ja  (:  ) ! (  NBODBTJA_el)

 END Module Jacket
!----------------------------------------------------------------------
 Subroutine Jacket_trans
!----------------------------------------------------------------------

 Use Cbeam
 Use Jacket

   implicit none

   real(8) :: X(3,2), Ex(3),Ey(3),Ez(3),met
   real(8) :: A(3,3), AT(3,3) !t,Y(3),ALEN
   integer :: nbod, nbsub, nn_el, m1, m2, nb_el


   write(10,*)
   write(10,*)'direction cosines'

   do    nbod     = 1+NBODBTWT_el, NBODBT_el
      do nbsub    = 1, body(nbod)%NBODSUB_el
 
         nb_el    = body(nbod)%NBODTGLB_el(nbsub)
         nn_el    = body(nbod)%NNODTGLB_el(nbsub)

         m1       = Nod_ja(1,nbod-NBODBTWT_el)
         m2       = Nod_ja(2,nbod-NBODBTWT_el)

         X(1:3,1) = XN_ja(1:3,m1)
         X(1:3,2) = XN_ja(1:3,m2)

         Ey(1:3)  = X(1:3,2)-X(1:3,1)
         met      = dsqrt(dot_product (Ey,Ey))
!t       ALEN=met
!!!      write(*,*)'Ey metro', met,Ey
         Ey(1:3)  = Ey(1:3)/met
!!!      write(*,*)'Ey fin', Ey

!-- Ex(1:3) are chosen to be same, as local tower c.system,
!-- in order to set easily the dependent q's

         Ex(1)    = Ey(3)  ! -Ey(2)
         Ex(2)    = 0.d0   !  Ey(1)
         Ex(3)    =-Ey(1)  !  0.d0
         met      = dsqrt(dot_product (Ex,Ex))

         if(met == 0) then
!!!         write(*,*)'Ex0 metro', nbod, met
!!!         write(*,*) Ex
            Ex(1) = 0.d0
            Ex(2) =-Ey(3)
            Ex(3) = Ey(2)
            met   = dsqrt(dot_product (Ex,Ex))
         endif

         if (met == 0) then !ckeck
            write(*,*) 'metro = 0 in Jacket_trans'
            stop
         endif
!!!      write(*,*)'Ex metro', nbod, met
!!!      write(*,*) Ex
         Ex(1:3)  = Ex(1:3)/met

         call EXTEPR ( Ex, Ey, Ez )
         met      = dsqrt(dot_product (Ez,Ez))
         Ez       = Ez/met

         AT(1,1:3)= Ex(1:3)
         AT(2,1:3)= Ey(1:3)
         AT(3,1:3)= Ez(1:3)

         A        = transpose(AT)

!t       Y=0.d0
!t       Y(2)=ALEN
!t       write(*,*)'given', X(1:3,2)
!t       write(*,*)'calc ', X(1:3,1)+matmul(A(1:3,1:3),Y(1:3))

         transf_mat(nn_el)%A_el     (1:3,1:3) = A (1:3,1:3)
         transf_mat(nn_el)%AT_el    (1:3,1:3) = AT(1:3,1:3)
         transf_mat(nn_el)%R_el     (1:3    ) = X (1:3,1  )

!cut start
!qq      if (ICUT_ja == 1) then
!qq         Y    = 0.d0
!qq         Y(2) = Lcut_ja(nb_el,1)
!qq         transf_mat(nn_el)%R_el (1:3) = X (1:3,1) + matmul ( A(1:3,1:3), Y(1:3) )
!qq      endif
!cut end

         transf_mat(nn_el)%ATDA_el  (1:3,1:3) = 0.0;
         transf_mat(nn_el)%ATDDA_el (1:3,1:3) = 0.0;
         transf_mat(nn_el)%ATDR_el  (1:3    ) = 0.0;
         transf_mat(nn_el)%ATDDR_el (1:3    ) = 0.0;

         transf_mat(nn_el)%DA_el    (1:3,1:3) = 0.0;
         transf_mat(nn_el)%DR_el    (1:3    ) = 0.0;
         write(10,100) nbod, Ex(1:3), Ey(1:3), Ez(1:3)
      enddo!nbsub
   enddo!nbod


   if (ICASE_el == 4) then
      do nn_el = NNODTWT_el+1 , NNODT_el
         transf_mat(nn_el)%Aja_el (1:3,1:3 ) = transf_mat(nn_el)%A_el (1:3,1:3);
         transf_mat(nn_el)%Rja_el (    1:3 ) = transf_mat(nn_el)%R_el (    1:3);
      enddo
   endif

 100 format(i5,10f16.10)


 END Subroutine Jacket_trans
!-----------------------------------------------------------------------
!
!-- Subroutine :INIT_Atrans
!
!   Rotation Matrices needed for B.C.
!
!-----------------------------------------------------------------------
 Subroutine INIT_Atrans
!-----------------------------------------------------------------------

 Use Cbeam
 Use Jacket

   implicit none

   real(8) :: X0L(3), X0L0(3), AG1(6,6),AG2(6,6)
   integer :: nbod,nbsub,nb_el,nbc,ibcon,iel,inod,i
   integer :: ibsubcon,isbcon,iecon,incon,icon
     

   do    nbod  = 1+NBODBTWT_el, NBODBT_el
   do    nbsub = 1, body(nbod)%NBODSUB_el
         nb_el = body(nbod)%NBODTGLB_el(nbsub)
      do nbc   = 1, boundco (nb_el)%nbcpb

         ibcon    = boundco (nb_el)%nbcon (nbc)
         if (ibcon == 0) cycle

         iel      = boundco (nb_el)%nel (nbc)
         inod     = boundco (nb_el)%nod (nbc)
         i        = inod/NNPE_el

         ibsubcon = 1
         isbcon   = body    (ibcon)%NBODTGLB_el(ibsubcon)
         iecon    = boundco (nb_el)%nelcon(nbc)
         incon    = boundco (nb_el)%nodcon(nbc)
         icon     = incon/NNPE_el

            X0L (:) = 0.d0;  X0L (2) = subbody(nb_el )%HTA_el(iel  +i   )
            X0L0(:) = 0.d0;  X0L0(2) = subbody(isbcon)%HTA_el(iecon+icon)

            call BOUND_TRANF (nbod, ibcon, X0L, X0L0, AG1, AG2)

         boundco (nb_el)%Atrans1  (nbc,1:6,1:6) = AG1(1:6,1:6);
         boundco (nb_el)%Atrans2  (nbc,1:6,1:6) = AG2(1:6,1:6);
      enddo !nbc

   enddo !nbsub
   enddo !nbod


 END Subroutine INIT_Atrans
!-----------------------------------------------------------------------
!
!-- Subroutine :BOUND_TRANF
!   Transformation operations for the implementation of boundary    
!   condition
!-----------------------------------------------------------------------
 Subroutine BOUND_TRANF (nbod1, nbod2, X0L1, X0L2, AG_bc, AG_load)
!-----------------------------------------------------------------------

 Use Cbeam
 Use Jacket

   implicit none

   real(8) :: AG_bc(6,6), AG_load(6,6)
   real(8) :: A1(3,3), A2(3,3), AT1(3,3), AT2(3,3), RG1(3), RG2(3)
   real(8) :: X0L1(3), X0L2(3), X2Loc(3), A12(3,3), A21(3,3)
   integer :: nbod1, nbod2, nn1_el, nn2_el


!--- body 1 give  loads         to body 2
!--- body 2 set the deflections to body 1
   nn1_el = body(nbod1)%NNODTGLB_el(1)
   nn2_el = body(nbod2)%NNODTGLB_el(1)

   A1     = transf_mat(nn1_el)%A_el
   A2     = transf_mat(nn2_el)%A_el
   AT1    = transf_mat(nn1_el)%AT_el
   AT2    = transf_mat(nn2_el)%AT_el
   RG1    = transf_mat(nn1_el)%R_el  + matmul(A1,X0L1)
   RG2    = transf_mat(nn2_el)%R_el  + matmul(A2,X0L2)

   X2Loc  = matmul   (AT2, RG1-RG2)
   A21    = matmul   (AT1, A2     )
   A12    = transpose(A21         )

   AG_bc ( 1:6, 1:6 )                     = 0.d0;
   AG_bc ( IMDOF_el(1:3), IMDOF_el(1:3) ) = A21 ( 1:3, 1:3 )
   AG_bc ( IMDOF_el(4:6), IMDOF_el(4:6) ) = A21 ( 1:3, 1:3 )
   AG_load                                = transpose(AG_bc)

!qq jim REMOVE the extra terms!
!  X2Loc(1:3) = 0.d0;
   write(1001,*) nbod1-5,nbod2-5
   write(1001,*) X2Loc(1:3)
   write(1001,*) X0L2(2)+X2Loc(2)

!qq permit the y offset!
!  X2Loc(2) = 0.d0

!--- Add terms in case of offset connections
   AG_bc  ( IMDOF_el(1:3),IMDOF_el(4) ) = X2Loc(2)*A21(1:3,3)-X2Loc(3)*A21(1:3,2)
   AG_bc  ( IMDOF_el(1:3),IMDOF_el(5) ) = X2Loc(3)*A21(1:3,1)-X2Loc(1)*A21(1:3,3)
   AG_bc  ( IMDOF_el(1:3),IMDOF_el(6) ) = X2Loc(1)*A21(1:3,2)-X2Loc(2)*A21(1:3,1)

   AG_load( IMDOF_el(4),IMDOF_el(1:3) ) = X2Loc(2)*A12(3,1:3)-X2Loc(3)*A12(2,1:3)
   AG_load( IMDOF_el(5),IMDOF_el(1:3) ) = X2Loc(3)*A12(1,1:3)-X2Loc(1)*A12(3,1:3)
   AG_load( IMDOF_el(6),IMDOF_el(1:3) ) = X2Loc(1)*A12(2,1:3)-X2Loc(2)*A12(1,1:3)


 END Subroutine BOUND_TRANF
!
!
!
!----------------------------------------------------------------------
 Subroutine read_jacket_nod_bod
!----------------------------------------------------------------------

 Use Cbeam
 Use Jacket
 Use Paths

   implicit none

   integer :: i,n


   if ((ICASE_el /= 2).and.(ICASE_el /= 4)) then
      NBODBTJA_el = 0
      NBODTJA_el  = 0
      return
   endif

   open(31,file=trim(file_jacket)) !'jacket.inp'


!--- Read jacket NODES. 3D position wrt global c.o. system
!--- BC_FLAG_ja :  0: No BC
!                  1: Zero BC
!                  2: Connection with tower  *************[NOT WORKING ATM]
   read(31,*) !title
   read(31,*) !** nodes
   read(31,*) NNT_ja
   read(31,*) !**  n       X         Y        Z      b.c.


     Allocate ( XN_ja     (3,NNT_ja  ) ) 

   do i=1,NNT_ja
      read(31,*) n, XN_ja(1,i), XN_ja(2,i), XN_ja(3,i)
   enddo


!--- Read jacket's BODIES given in undeformed condition by setting 1st and last node
   read(31,*) !** bodies
   read(31,*) NBODBTJA_el
   read(31,*) !** nb   n1   n2   NTEB  NBC  Iflooded

!df :: IFloodbod -->morison_typ

   NBODBT_el = NBODBTWT_el+NBODBTJA_el

   Allocate (IFloodbod_el(  NBODBT_el  ))
   Allocate (Nod_ja      (2,NBODBTJA_el))
   Allocate (NTEB_ja     (  NBODBTJA_el))
   Allocate (NBC_ja      (  NBODBTJA_el))
   Allocate (ITYPE_ja    (  NBODBTJA_el))

   IFloodbod_el=0;

   do i=1,NBODBTJA_el
      read(31,*)                    &
       n                           ,&       !-- Body's index number[jacket only]
       Nod_ja(1,i)                 ,&       !-- Body's start node
       Nod_ja(2,i)                 ,&       !-- Body's end   node
       NTEB_ja (i)                 ,&       !-- Body's number of total elements
       NBC_ja  (i)                 ,&       !-- Body's number of boundary conditions
       IFloodbod_el (i+NBODBTWT_el),&       !-- Body's free-flooded flag [0:no, 1:yes]
       ITYPE_ja(i)                          !-- Body's type of loads     [0:no, 2:hydro]

!------ Perform checks
      if (i/=n) then
        write(*,*) 'CHECK jacket input file'
        stop
      endif
      if (XN_ja(3,Nod_ja(2,i)) < XN_ja(3,Nod_ja(1,i))) then
          write(*,*) 'Warning, Jacket nod2 is upper than nod1',n,XN_ja(3,Nod_ja(1,i)), XN_ja(3,Nod_ja(2,i))
      endif
   enddo


!--- machi for jacket
   read(31,*) !** Number of different structural-hydrodynamic properties
   read(31,*) INODNELPROPJA_el


 END Subroutine read_jacket_nod_bod
!----------------------------------------------------------------------
 Subroutine read_jacket_machi_bc
!----------------------------------------------------------------------

 Use Cbeam
 Use Jacket
 Use Paths

   implicit none

   integer :: nbod, nbsub, nb_el,nbc, nbod_ja, nbc_read, nbcon, nbcon_ja, isbcon, j
   real(8) :: EIXXM,EIZZM,GAXM ,GAZM, zer
   integer :: NUM, nel, nel_ja, n, n1, n2


   if ((ICASE_el /= 2).and.(ICASE_el /= 4)) return


!--- Read structural-hydrodynamic properties
   read(31,*) !**num...


   do n = 1+INODNELPROPWT_el,INODNELPROP_el
       beam_timo (n)%K( :) = 0.d0;
      read (31,*)    NUM      ,                                             &
       beam_timo (n)%DENS_el  , beam_timo(n)%XCMD_el , beam_timo(n)%ZCMD_el,&
       beam_timo (n)%RIXX_el  , beam_timo(n)%RIZZ_el , beam_timo(n)%RIXZ_el,&
       beam_timo (n)%K( 7)    , beam_timo(n)%K( 9)   , beam_timo(n)%K(11)  ,& !EA_el    , EAX_el*  , EAZ_el
       beam_timo (n)%K(16)    , beam_timo(n)%K(21)   , beam_timo(n)%K(18)  ,& !EIXX_el  , EIZZ_el  , EIXZ_el*
       beam_timo (n)%K(19)                                                 ,& !GIT_el
       beam_timo (n)%K( 1)    , beam_timo(n)%K(12)                         ,& !GAX_el   , GAZ_el   
       beam_timo (n)%K( 5)    , beam_timo(n)%K(14)                         ,& !GAX_X_el , GAZ_Z_el*
       substr_mor(n)%DIAMET_el, substr_mor(n)%INDIAMET_el                  ,&
       substr_mor(n)%Cm_el    , substr_mor(n)%Cd_el

!------- change sign because the stiffness matrix is set and not the elastic properties [vars with *]
       beam_timo (n)%K( 9)    = -beam_timo(n)%K( 9) !EAX
       beam_timo (n)%K(14)    = -beam_timo(n)%K(14) !GAZ_Z
      !beam_timo (n)%K(17)    = -beam_timo(n)%K(17) !EIXY
       beam_timo (n)%K(18)    = -beam_timo(n)%K(18) !EIXZ
      !beam_timo (n)%K(20)    = -beam_timo(n)%K(20) !EIYZ

       beam_timo (n)%POLI_el  = beam_timo(n)%RIXX_el + beam_timo(n)%RIZZ_el
       beam_timo (n)%AMOMX_el = beam_timo(n)%DENS_el * beam_timo(n)%ZCMD_el
       beam_timo (n)%AMOMZ_el = beam_timo(n)%DENS_el * beam_timo(n)%XCMD_el

!------ Stiff Dofs Enabled for all jacket properties
      if (body(NBODBT_el)%IDOF_el == 0 ) then
       beam_timo (n)%K( :) = 0.d0;
       beam_timo (n)%K( 1) = BIG   !GAX  !(1,1)
       beam_timo (n)%K( 7) = BIG   !EA   !(2,2)
       beam_timo (n)%K(12) = BIG   !GAZ  !(3,3)
       beam_timo (n)%K(16) = BIG   !EIXX !(4,4)
       beam_timo (n)%K(19) = BIG   !GIT  !(5,5)
       beam_timo (n)%K(21) = BIG   !EIZZ !(6,6)
      endif

   enddo !i


!--- Read Body elements length, corresponding property
!--- Set subbody:: HTA_el, ALENG_el, ALENGB_el, PHIX_el, PHIZ_el, INODNEL_el
   read(31,*) !**  nb  nel          hta1             prop1        hta2             prop2

   do nbod  = NBODBTWT_el+1, NBODBT_el
   do nbsub = 1, body(nbod)%NBODSUB_el
      nb_el = body(nbod)%NBODTGLB_el(nbsub)

      do nel = 1, subbody(nb_el)%NTEB_el

         if (nel==1)                     &
          read(31,*)                     &
            nbod_ja                     ,&
            nel_ja                      ,&
            subbody(nb_el)%HTA_el(nel  ),&
            n1                          ,&
            subbody(nb_el)%HTA_el(nel+1),&
            n2   

         if (nel> 1)                     &
          read(31,*)                     &
            nel_ja                      ,&
            subbody(nb_el)%HTA_el(nel  ),&
            n1                          ,&
            subbody(nb_el)%HTA_el(nel+1),&
            n2   

         if ((nbod_ja+NBODBTWT_el/=nbod).or.(nel_ja/=nel)) then
            write(*,*) 'error reading jacket body properties'
            stop
         endif

         n1 = n1 + INODNELPROPWT_el
         n2 = n2 + INODNELPROPWT_el
         subbody(nb_el)%INODNEL_el(nel,1) = n1
         subbody(nb_el)%INODNEL_el(nel,2) = n2
         subbody(nb_el)%ALENG_el  (nel  ) = subbody(nb_el)%HTA_el(nel+1) - subbody(nb_el)%HTA_el  (nel)  !element's length
         subbody(nb_el)%ALENGB_el         = subbody(nb_el)%ALENGB_el     + subbody(nb_el)%ALENG_el(nel)  !subbody's length

!------------ Perform checks of input data [should be positive]
            zer = 1.d-10
            if (subbody  (nb_el)%ALENG_el(nel) < zer                        ) then; write(*,'(a,3i4)') 'check HTA  in machi',nbod,nbsub,nel;stop;endif;
            if (beam_timo(n1)%DENS_el < zer .or. beam_timo(n2)%DENS_el < zer) then; write(*,'(a,3i4)') 'check DENS in machi',nbod,nbsub,nel;stop;endif;
            if (beam_timo(n1)%RIXX_el < zer .or. beam_timo(n2)%RIXX_el < zer) then; write(*,'(a,3i4)') 'check RIXX in machi',nbod,nbsub,nel;stop;endif;
            if (beam_timo(n1)%RIZZ_el < zer .or. beam_timo(n2)%RIZZ_el < zer) then; write(*,'(a,3i4)') 'check RIZZ in machi',nbod,nbsub,nel;stop;endif;
            if (beam_timo(n1)%K( 1)   < zer .or. beam_timo(n2)%K( 1)   < zer) then; write(*,'(a,3i4)') 'check GAX  in machi',nbod,nbsub,nel;stop;endif;
            if (beam_timo(n1)%K( 7)   < zer .or. beam_timo(n2)%K( 7)   < zer) then; write(*,'(a,3i4)') 'check EA   in machi',nbod,nbsub,nel;stop;endif;
            if (beam_timo(n1)%K(12)   < zer .or. beam_timo(n2)%K(12)   < zer) then; write(*,'(a,3i4)') 'check GAZ  in machi',nbod,nbsub,nel;stop;endif;
            if (beam_timo(n1)%K(16)   < zer .or. beam_timo(n2)%K(16)   < zer) then; write(*,'(a,3i4)') 'check EIXX in machi',nbod,nbsub,nel;stop;endif;
            if (beam_timo(n1)%K(19)   < zer .or. beam_timo(n2)%K(19)   < zer) then; write(*,'(a,3i4)') 'check GIT  in machi',nbod,nbsub,nel;stop;endif;
            if (beam_timo(n1)%K(21)   < zer .or. beam_timo(n2)%K(21)   < zer) then; write(*,'(a,3i4)') 'check EIZZ in machi',nbod,nbsub,nel;stop;endif;

         GAXM  = beam_timo(n1)%K( 1) + beam_timo(n2)%K( 1) !GAX  !(1,1)
         GAZM  = beam_timo(n1)%K(12) + beam_timo(n2)%K(12) !GAZ  !(3,3)
         EIXXM = beam_timo(n1)%K(16) + beam_timo(n2)%K(16) !EIXX !(4,4)
         EIZZM = beam_timo(n1)%K(21) + beam_timo(n2)%K(21) !EIZZ !(6,6) 

         subbody(nb_el)%PHIX_el  (nel) = 12.d0*EIZZM/(GAXM*subbody(nb_el)%ALENG_el(nel)**2)
         subbody(nb_el)%PHIZ_el  (nel) = 12.d0*EIXXM/(GAZM*subbody(nb_el)%ALENG_el(nel)**2)
      enddo !nel
   enddo !nbsub
   enddo !nbod



!--- Read Boundary Conditions --> Connections
   read(31,*) !** Boundary Conditions
   read(31,*) !** nb   nbc  nel  nod   nbcon nel  nod   u  v  w θx θy θz !local, 0:free, 1:bc applied

   do nbod    = NBODBTWT_el+1, NBODBT_el
      nbod_ja = nbod - NBODBTWT_el 
      nbsub   = 1
      nb_el   = body(nbod)%NBODTGLB_el(nbsub)

      boundco (nb_el)%nbcpb = NBC_ja(nbod_ja)

      if (boundco (nb_el)%nbcpb==0) cycle

      nbc = boundco (nb_el)%nbcpb

      Allocate ( boundco (nb_el)%nel     (nbc    ) )
      Allocate ( boundco (nb_el)%nod     (nbc    ) )
      Allocate ( boundco (nb_el)%nbcon   (nbc    ) )
      Allocate ( boundco (nb_el)%nelcon  (nbc    ) )
      Allocate ( boundco (nb_el)%nodcon  (nbc    ) )
      Allocate ( boundco (nb_el)%indx    (nbc,6  ) )
      Allocate ( boundco (nb_el)%Atrans1 (nbc,6,6) )
      Allocate ( boundco (nb_el)%Atrans2 (nbc,6,6) )

      isbcon = 1

      do nbc = 1, boundco (nb_el)%nbcpb

         if (nbc==1)                                 &
          read(31,*)                                 &
            nbod_ja                                 ,&
            nbc_read                                ,&
            boundco (nb_el)%nel   (nbc             ),&
            boundco (nb_el)%nod   (nbc             ),&      !1 or NNPE_el
            nbcon_ja                                ,&
            boundco (nb_el)%nelcon(nbc             ),&
            boundco (nb_el)%nodcon(nbc             ),&      !1 or NNPE_el
           (boundco (nb_el)%indx  (nbc, IMDOF_el(j)),j=1,6)

         if (nbc>1)                                  &
          read(31,*)                                 &
            nbc_read                                ,&
            boundco (nb_el)%nel   (nbc             ),&
            boundco (nb_el)%nod   (nbc             ),&      !1 or NNPE_el
            nbcon_ja                                ,&
            boundco (nb_el)%nelcon(nbc             ),&
            boundco (nb_el)%nodcon(nbc             ),&      !1 or NNPE_el
           (boundco (nb_el)%indx  (nbc, IMDOF_el(j)),j=1,6)

!---------- manage the different number of element nodes
         if(boundco (nb_el)%nod   (nbc)>1)&
            boundco (nb_el)%nod   (nbc) = NNPE_el
         if(boundco (nb_el)%nodcon(nbc)>1)&
            boundco (nb_el)%nodcon(nbc) = NNPE_el

         if (nbcon_ja==0) then
                    boundco (nb_el)%nbcon (nbc) = 0          !-- remains zero, when zero b.c. are applied 
                                                             !-- nelcon and nodcon are not used
         else
                    nbcon = nbcon_ja+NBODBTWT_el
                    boundco (nb_el)%nbcon (nbc) = nbcon
         endif

         if ((nbod/=nbod_ja+NBODBTWT_el).or.(nbc/=nbc_read)) then
            write(*,*) 'error in reading boundary conditions'
            stop
         endif

         write(*,1) nbod_ja                                 ,&
                    nbc_read                                ,&
                    boundco (nb_el)%nel   (nbc             ),&
                    boundco (nb_el)%nod   (nbc             ),&
                    nbcon_ja                                ,&
                    boundco (nb_el)%nelcon(nbc             ),&
                    boundco (nb_el)%nodcon(nbc             ),&
                   (boundco (nb_el)%indx  (nbc, IMDOF_el(j)),j=1,6)
      enddo !nbc
   enddo !nbod

 1 format (15i4)

   Deallocate ( NBC_ja )



!--- Read where to calculate Buoyancy on tapers (side surfaces)
!--- Careful: Give the correct jacket member, without adding WT bodies
   read(31,*) !** Buoyancy Caps
   read(31,*) NBuoyaTap_el

   if (NBuoyaTap_el>0) then
      read(31,*) !** nb   iend(1,2)
      Allocate (buoyancy_tap(NBuoyaTap_el))
   endif

   do j = 1, NBuoyaTap_el
      read(31,*)                   &
        buoyancy_tap(j)%nbod     , &  !jacket bodies without WT
        buoyancy_tap(j)%nel      , &
        buoyancy_tap(j)%nod      , &  !1 or NNPE_el
        buoyancy_tap(j)%itype    , &
        buoyancy_tap(j)%area     , &
        buoyancy_tap(j)%ca       , &
        buoyancy_tap(j)%cd       , &
        buoyancy_tap(j)%ifroude  , &
        buoyancy_tap(j)%idyn_pres, &
        buoyancy_tap(j)%irel

        buoyancy_tap(j)%nbod = buoyancy_tap(j)%nbod + NBODBTWT_el
!---------- manage the different number of element nodes
     if(buoyancy_tap(j)%nod>1)&
        buoyancy_tap(j)%nod  = NNPE_el
   enddo

   close(31)

   write(10,*)
   write(10,*) 'NBuoyaTap_el=',NBuoyaTap_el
   do j = 1, NBuoyaTap_el
      write(10,2)                  &
        buoyancy_tap(j)%nbod     , &  !jacket bodies without WT
        buoyancy_tap(j)%nel      , &
        buoyancy_tap(j)%nod      , &  !1 or NNPE_el
        buoyancy_tap(j)%itype    , &
        buoyancy_tap(j)%area     , &
        buoyancy_tap(j)%ca       , &
        buoyancy_tap(j)%cd       , &
        buoyancy_tap(j)%ifroude  , &
        buoyancy_tap(j)%idyn_pres, &
        buoyancy_tap(j)%irel
   enddo
   write(10,*)

 2  format (4i5,3f14.4,3i5)


 END Subroutine read_jacket_machi_bc
