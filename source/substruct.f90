!-------------------------------------------------------------------------------
!   Subroutine :FHYDLOC_mor
!
!  Calculate Local Hydrodynamic Matrices from Morison and Buoyancy
!
!-------------------------------------------------------------------------------
 Subroutine Morison_hyd ( AMLOC, ACLOC, AKLOC, AFLOC, AMLOCNQ, ACLOCNQ, AKLOCNQ,&
                          nbod , nbsub, nb_el, nn_el, nel                         )
!-------------------------------------------------------------------------------

 Use Cbeam
 Use Hydro

   implicit none

   integer,                                 intent(in   ) :: nbod, nbsub, nb_el, nn_el, nel
   real(8), dimension(NDFPEM_el,NDFPEM_el), intent(inout) :: AKLOC  , ACLOC  , AMLOC
   real(8), dimension(NDFPEM_el,NQS      ), intent(inout) :: AKLOCNQ, ACLOCNQ, AMLOCNQ
   real(8), dimension(NDFPEM_el          ), intent(inout) :: AFLOC

   real(8), dimension(NDFPEM_el,NDFPEM_el) :: AKLOC0  , ACLOC0  , AMLOC0
   real(8), dimension(NDFPEM_el          ) :: AKLOCNQ0, ACLOCNQ0, AMLOCNQ0
   real(8), dimension(NDFPEM_el          ) :: AFLOC0     
   real(8), dimension(NDFPEM_el)           :: UT   , UT1  , UT2
   real(8), dimension(NEQPEM_el)           :: UTT_A, UTT_B    

   real(8) :: SHAPET    (NDFPEM_el,NEQPEM_el)
   real(8) :: SHAPE     (NEQPEM_el,NDFPEM_el)
   real(8) :: SHAPETAM  (NDFPEM_el,NEQPEM_el)
   real(8) :: AM        (NEQPEM_el,NEQPEM_el)
   real(8) :: AQ        (NEQPEM_el          )

   real(8) :: XG(3), Y0(3), DL

!--- Vars Passed to Morison
   real(8), dimension (NEQPEM_el) :: UTT, UTT1, UTT2  
!! real(8) :: HTA1
   real(8) :: L_Unit_t(3)
   real(8) :: DIAM, CDrag, CMass !!, Fi
   real(8) :: AKL     (3  ,3)
   real(8) :: ACL     (3  ,3)
   real(8) :: AML     (3  ,3)
   real(8) :: AKLnq   (NQS,3)
   real(8) :: ACLnq   (NQS,3)
   real(8) :: AMLnq   (NQS,3)
   real(8) :: AFL     (3    )

   integer :: k, n1, n2, NDFPE, NEQPE, it
   real(8) :: ALLOC, PA01, PA02, HTA, HTAL, HTA1, HTA0, PHIX, PHIZ
   real(8) ::        zita, L1  , L2
   integer :: iq, n

   integer :: N_low, N_high, Ihalfwet, N_low_tmp
   real(8) :: Z_Calc, Length, Length0, lamda, XG_nod(3,2), Ey(3), LengthP
   real(8) :: Z_Calc_high, Z_Calc_low


!--- hydrodynamic forces not included for this subbody
   if (subbody(nb_el)%IBTYP_el /= 2) return


   HTA0    = subbody(nb_el)%HTA_el    (nel  )
   ALLOC   = subbody(nb_el)%ALENG_el  (nel  )
   n1      = subbody(nb_el)%INODNEL_el(nel,1)
   n2      = subbody(nb_el)%INODNEL_el(nel,2)
   NDFPE   = subbody(nb_el)%NDFPE_el
   NEQPE   = subbody(nb_el)%NEQPE_el
   PHIX    = subbody(nb_el)%PHIX_el   (nel  )
   PHIZ    = subbody(nb_el)%PHIZ_el   (nel  )


   call LOCAL_UT_1   (nb_el, nel, UT, UT1, UT2)


!--- 1st node of the elemenet
!! call SSHAPEFUNC15 ( 0.d0, ALLOC, PHIX, PHIZ,&  !-- HTAL=0.d0
!!                     SHAPE, NDFPE, NEQPE      )

   Y0    (1:3  ) = 0.d0;
   Y0    (2    ) = HTA0
!! UTT           = matmul (SHAPE, UT )
!! XG_nod(1:3,1) = transf_mat(nn_el)%R_el(1:3) + matmul ( transf_mat(nn_el)%A_el(1:3,1:3), Y0(1:3) + UTT(IMDOF_el(1:3))  )
   XG_nod(1:3,1) = transf_mat(nn_el)%R_el(1:3) + matmul ( transf_mat(nn_el)%A_el(1:3,1:3), Y0(1:3) )

!--- last node of the elemenet
!! call SSHAPEFUNC15 ( 1.d0, ALLOC, PHIX, PHIZ,&   !-- HTAL=1.d0
!!                     SHAPE, NDFPE, NEQPE       )

   Y0    (2    ) = HTA0 + ALLOC
!! UTT           = matmul (SHAPE, UT )
!! XG_nod(1:3,2) = transf_mat(nn_el)%R_el(1:3) + matmul ( transf_mat(nn_el)%A_el(1:3,1:3), Y0(1:3) + UTT(IMDOF_el(1:3)) )
   XG_nod(1:3,2) = transf_mat(nn_el)%R_el(1:3) + matmul ( transf_mat(nn_el)%A_el(1:3,1:3), Y0(1:3) )


!--- check which node is lower
   if (XG_nod(3,2)>=XG_nod(3,1)) then
      N_low  = 1
      N_high = 2
   else
      N_low  = 2
      N_high = 1
   endif


   if (XG_nod(3,N_low ) >=  Water_upper_lim) return            !-- both nodes out of water
   if (XG_nod(3,N_high) <= -depth          ) return            !-- the element is below mudline (foundation)
                                             Ihalfwet = 0      !-- indicates if the body is half wet, half dry [1]
   if (XG_nod(3,N_high) >   Water_lower_lim) Ihalfwet = 1


   Ey(1:3)  = XG_nod(1:3,N_high) - XG_nod(1:3,N_low)           !-- vector tangential
   Length   = dsqrt(dot_product (Ey,Ey))
   Length0  = Length

   if     (Ihalfwet==0) then
      L1        = 0.d0
      L2        = Length0
!------ hydro Results
      Res_hyd (nbod, nel, 1, 31) = Length
   elseif (Ihalfwet==1) then
!------ hydro Results
      Res_hyd (nbod, nel, 1, 31) = 0.d0                        !-- record the value in case of return
!------ X-X1=λ(X2-X1) ---> λ = (X(3)-X1(3)) / (X2(3)-X1(3)) 
!------               ---> X(1:2) = X1(1:2) + λ (X2(1:2)-X1(1:2))

      if     ( StretchMod == 0 ) then
         Z_Calc             = tide

         lamda              = (Z_calc-XG_nod(3,N_low)) / Ey(3)
         XG_nod(1:2,N_high) = XG_nod(1:2,N_low) + lamda*Ey(1:2)
         XG_nod(  3,N_high) = Z_calc

         Ey(1:3)            = XG_nod(1:3,N_high) - XG_nod(1:3,N_low)
         Length             = dsqrt(dot_product (Ey,Ey))
      elseif ( StretchMod == 1 ) then
!!!      write(*,*)
!!!      write(*,*)'bike',nbod,nel,N_high,Length0

!--------- Lower Node with water elevation
         XG(1:3)    = XG_nod(1:3,N_low)
         call ZELEV_sea ( XG, zita )
         Z_Calc_low = zita + tide

!--------- Higher Node with water elevation
         XG(1:3)     = XG_nod(1:3,N_high)
         call ZELEV_sea ( XG, zita  )
         Z_Calc_high = zita + tide

!--------- both nodes out of water and no hydrodynamic force is calculated
         if ( (XG_nod(3,N_low ) >=  Z_Calc_low ).and.&
              (XG_nod(3,N_high) >=  Z_Calc_high)       ) return

!--------- both nodes are inside water, the whole element is wet
         if ( (XG_nod(3,N_high) <=  Z_calc_high).and.&
              (XG_nod(3,N_low ) <=  Z_calc_low )       ) then
            Length   = Length0
         else
!--------- The element is half wet - half dry.
!------------ switch N_low, N_high if the lower node is dry and the higher is wet
            if ( (XG_nod(3,N_low ) >=  Z_Calc_low ).and.&
                 (XG_nod(3,N_high) <=  Z_Calc_high)       ) then
!               write(* ,*)'SOS switching lower-higher node',nbod,nel
                write(10,*)'SOS switching lower-higher node',nbod,nel
                N_low_tmp = N_low
                N_low     = N_high
                N_high    = N_low_tmp
                Ey(1:3)   =-Ey(1:3)
            endif

            if (dabs(Ey(3))<1.d-5) then
!              write(* ,*)'SOS Ey is zero',Ey(3),nbod,nel
               write(10,*)'SOS Ey is zero',Ey(3),nbod,nel
            endif

!------------ Iteration for determining the actual wet length, when integrating up to current water elevation
            do it = 1, 101
               LengthP  = Length
               XG(1:3)  = XG_nod(1:3,N_high)
               call ZELEV_sea ( XG, zita )
               Z_Calc   = zita + tide

               lamda              = (Z_calc-XG_nod(3,N_low)) / Ey(3)
               XG_nod(1:2,N_high) = XG_nod(1:2,N_low) + lamda*Ey(1:2)
               XG_nod(  3,N_high) = Z_calc

               Ey(1:3)            = XG_nod(1:3,N_high) - XG_nod(1:3,N_low)
               Length             = dsqrt(dot_product (Ey,Ey))

!!!            write(*,*)'it,XG,XG_nod',it,Length
!!             write(*,*)'XG,XG_nod',Length
!!             write(*,*)XG    (1:3)
!!             write(*,*)XG_nod(1:3,N_high)
!!             write(*,*) dabs( LengthP - Length )
               
               if (dabs( LengthP - Length ) < 1.d-5) goto 1
            enddo

            write(* ,*)"SOS wet length wasn't caclulated",dabs( LengthP - Length ), Length
            write(10,*)"SOS wet length wasn't caclulated",dabs( LengthP - Length ), Length
!qq
      Res_hyd (nbod, nel, 1, 31) = Length
            return

   1        continue
         endif
      endif !StretchMod

!------ hydro Results
      Res_hyd (nbod, nel, 1, 31) = Length

      if (N_low==1) then
         L1     = 0.d0
         L2     = Length
      else
         L1     = Length0 - Length
         L2     = Length0
      endif
   endif !Ihalfwet

!qq
   L2 = min(L2,ALLOC)

!! Fi = datan2(  (substr_mor(n2)%DIAMET_el-substr_mor(n1)%DIAMET_el), (2.d0*ALLOC)  )
   DL = ( L2 - L1 ) / NX_hyd

   if (DL<=0.d0) then
      write(* ,*) 'ZERO Lenth in Morison_hyd'
      write(10,*) 'ZERO Lenth in Morison_hyd'
      return
   endif


!--- Integrate normal drag and inertial Forces from Morison's equation
!!    HTAL           =  L1 /ALLOC

!!    call SSHAPEFUNC15 ( HTAL, ALLOC, PHIX, PHIZ,& 
!!                        SHAPE, NDFPE, NEQPE       )

!!    SHAPET         = transpose (SHAPE )
!!    UTT_B          = matmul    (SHAPE, UT )             !UTT_B = UTT_A;


   do k              = 1, NX_hyd
!------ local deformed Unit Vector parallel to beam axis
!!    UTT_A          = UTT_B
!!    HTAL           = (L1 + DL * dble(k))/ALLOC

!!    call SSHAPEFUNC15 ( HTAL, ALLOC, PHIX, PHIZ,& 
!!                        SHAPE, NDFPE, NEQPE       )

!!    SHAPET         = transpose (SHAPE )
!!    UTT_B          = matmul (SHAPE, UT )

!!    L_Unit_t (1:3) = UTT_B(IMDOF_el(1:3))-UTT_A(IMDOF_el(1:3))
!!    L_Unit_t (2  ) = L_Unit_t (2  ) + DL
!!    L_Unit_t (1:3) = L_Unit_t (1:3) / dsqrt ( dot_product ( L_Unit_t(1:3), L_Unit_t(1:3) ) )
!------ local undeformed Unit Vector
      L_Unit_t (1:3) = 0.d0;
      L_Unit_t (  2) = 1.d0;


!------ Mid point
      HTA      = L1 + dble(2*k-1)/2.d0 * DL
      HTAL     = HTA / ALLOC
      HTA1     = HTA0 + HTA
      PA01     = 1.d0 - HTAL
      PA02     = HTAL

      call SSHAPEFUNC15 ( HTAL, ALLOC, PHIX, PHIZ,& 
                          SHAPE, NDFPE, NEQPE       )

      SHAPET   = transpose (SHAPE )

      CMass    =  PA01 * substr_mor(n1)%Cm_el      +  PA02 * substr_mor(n2)%Cm_el
      CDrag    =  PA01 * substr_mor(n1)%Cd_el      +  PA02 * substr_mor(n2)%Cd_el
      DIAM     =  PA01 * substr_mor(n1)%DIAMET_el  +  PA02 * substr_mor(n2)%DIAMET_el

      UTT      = matmul (SHAPE, UT )
      UTT1     = matmul (SHAPE, UT1)
      UTT2     = matmul (SHAPE, UT2)

!qq
!!    UTT=0.d0; UTT1=0.d0; UTT2=0.d0; L_Unit_t=0.d0; L_Unit_t(2)=1

!------ Hydrodynamic forces from Morison's equation
      call MORISON ( UTT, UTT1, UTT2, HTA1 , L_Unit_t, DIAM , CDrag, CMass,       & !!Fi,
                     AKL, ACL , AML , AKLnq, ACLnq   , AMLnq, AFL  , nbod , nbsub, nel, k )


      AM                              = 0.d0;
      AM(IMDOF_el(1:3),IMDOF_el(1:3)) = AKL(1:3,1:3);
      SHAPETAM                        = matmul( SHAPET  , AM    )
      AKLOC0                          = matmul( SHAPETAM, SHAPE )
      AM                              = 0.d0;
      AM(IMDOF_el(1:3),IMDOF_el(1:3)) = ACL(1:3,1:3);
      SHAPETAM                        = matmul( SHAPET  , AM    )
      ACLOC0                          = matmul( SHAPETAM, SHAPE )
      AM                              = 0.d0;
      AM(IMDOF_el(1:3),IMDOF_el(1:3)) = AML(1:3,1:3);
      SHAPETAM                        = matmul( SHAPET  , AM    )
      AMLOC0                          = matmul( SHAPETAM, SHAPE )
      AQ                              = 0.d0;
      AQ (IMDOF_el(1:3))              = AFL(1:3)
      AFLOC0                          = matmul( SHAPET  , AQ    )

      do iq = 1, QSnode(nn_el)%Tot_el(2)
         n  = QSnode(nn_el)%IORD_el(iq)
         AQ                           = 0.d0;
         AQ (IMDOF_el(1:3))           = AKLnq (iq    , 1:3)
         AKLOCNQ0                     = matmul(SHAPET, AQ )
         AQ                           = 0.d0;
         AQ (IMDOF_el(1:3))           = ACLnq (iq    , 1:3)
         ACLOCNQ0                     = matmul(SHAPET, AQ )
         AQ                           = 0.d0;
         AQ (IMDOF_el(1:3))           = AMLnq (iq    , 1:3)
         AMLOCNQ0                     = matmul(SHAPET, AQ )

         AKLOCNQ(1:NDFPE,n)           = AKLOCNQ(1:NDFPE,n) + AKLOCNQ0(1:NDFPE) * DL
         ACLOCNQ(1:NDFPE,n)           = ACLOCNQ(1:NDFPE,n) + ACLOCNQ0(1:NDFPE) * DL
         AMLOCNQ(1:NDFPE,n)           = AMLOCNQ(1:NDFPE,n) + AMLOCNQ0(1:NDFPE) * DL
      enddo !iq

!qq
     if (IREDUCTJa_el==0.or.IEIG>0.or.IMODAL>0.or.ICALCINIT>0) then
!!----- Add Morison contribution to elastic matrices
      AKLOC = AKLOC + AKLOC0 * DL
      ACLOC = ACLOC + ACLOC0 * DL
      AMLOC = AMLOC + AMLOC0 * DL
     endif !LHS
      AFLOC = AFLOC + AFLOC0 * DL
   enddo !k


 END Subroutine Morison_hyd
!----------------------------------------------------------------------
!
!   Subroutine :MORISON
!
! Implementation of Morison equation for sub-structrures 
! (ICASE_el=1, 2 or 4) --> monopile, tripod or jacket, floating flexible body
!
! Local calculation (U, A tranformed from gl.c.s. to local)
!
! Output : Local Force (inertia and drag terms normal to local beam axis)
! -------  Local Added mass terms
!          Local Added damping terms
!
!         Froud Krylov       Added mass and zero freq. diffraction         Drag Force
! dF/dl =    ρ dVo a_norm + Ca ρ dVo (a-d2q/dt2)_norm + Cd ρ/2 dA |U-dq/dt|_norm (U-dq/dt)_norm [truss]
! dF/dl = Cm ρ dVo a_norm - Ca ρ dVo   (d2q/dt2)_norm + Cd ρ/2 dA |U-dq/dt|_norm (U-dq/dt)_norm [beam ]
!
!     U    , A       : flow    velocity and acceleration (from wave / current)
!     dq/dt, d2q/dt2 : elastic velocity and acceleration
!
!     Cm   = 1 + Ca
!     dVo  = π R ^2 dl
!     dA   = 2 R    dl , l: 0:L beam axis
!
! To do : Find if cone angle affects the force direction
!
!----------------------------------------------------------------------
 Subroutine MORISON ( UTT, UTT1, UTT2, HTA1 , L_Unit_t, DIAM , CDrag, CMass,       & !!Fi,
                      AKL, ACL , AML , AKLnq, ACLnq   , AMLnq, AFL  , nbod , nbsub, nel, k_nhyd )
!----------------------------------------------------------------------

 Use Hydro
 Use Cbeam

   implicit none

   real(8), dimension(NEQPEM_el), intent(in ) :: UTT  , UTT1  , UTT2
   real(8),                       intent(in ) :: HTA1 , L_Unit_t(3), DIAM, CDrag, CMass !!, Fi
   integer,                       intent(in ) :: nbod , nbsub , nel, k_nhyd
   real(8), dimension(3  ,3)    , intent(out) :: AKL  , ACL   , AML
   real(8), dimension(NQS,3)    , intent(out) :: AKLnq, ACLnq , AMLnq 
   real(8)                      , intent(out) :: AFL(3)

   real(8), dimension(3,3) :: ATDA , ATDDA , A     , AT
   real(8), dimension(3  ) :: ATDR , ATDDR , R
   real(8), dimension(3,3) :: ATDA0, ATDA1 , ATDDA0, ATDDA1, ATDDA2,AT0
   real(8), dimension(3  ) :: ATDR0, ATDR1 , ATDDR0, ATDDR1, ATDDR2
   real(8), dimension(3  ) :: Fdrag, Finer1, Finer2, FG
   real(8)                 :: XG    (3), Y0    (3)
   real(8)                 :: UG    (3), AG    (3), pdyn
   real(8)                 :: Uinf  (3), Ainf  (3)
   real(8)                 :: Uelst (3), Aelst (3)
   real(8)                 :: Ur    (3), Ur_n  (3), Ur_n_Meas
   real(8)                 :: CD       , CM1      , CM2
   real(8)                 :: LLt (3,3)
   real(8)                 :: DIAGmLLt(3,3)
   integer                 :: nn_el, i, j, iq


   nn_el      = body(nbod)%NNODTGLB_el (nbsub)

!--- Zero Local Matrices
   AKL  (:,:) = 0.d0;
   ACL  (:,:) = 0.d0;
   AML  (:,:) = 0.d0;
   AKLnq(:,:) = 0.d0;
   ACLnq(:,:) = 0.d0;
   AMLnq(:,:) = 0.d0;
   AFL  (:  ) = 0.d0;

!--- Pass Transformation matrices to local arrays
   A     (1:3,1:3) = transf_mat(nn_el)%A_el     (1:3,1:3)
   AT    (1:3,1:3) = transf_mat(nn_el)%AT_el    (1:3,1:3)
   ATDA  (1:3,1:3) = transf_mat(nn_el)%ATDA_el  (1:3,1:3)
   ATDDA (1:3,1:3) = transf_mat(nn_el)%ATDDA_el (1:3,1:3)
   R     (1:3    ) = transf_mat(nn_el)%R_el     (1:3    )
   ATDR  (1:3    ) = transf_mat(nn_el)%ATDR_el  (1:3    )
   ATDDR (1:3    ) = transf_mat(nn_el)%ATDDR_el (1:3    )


!--- 1. L_Unit_t: Local Unit vector parallel to beam axis [almost (0,1,0)]
!--- 2. XG      : global coordinates of the evaluation point
   Y0 (1:3) = 0.d0;
   Y0 (2  ) = HTA1
   XG (1:3) = R(1:3) + matmul ( A(1:3,1:3), Y0(1:3)+UTT(IMDOF_el(1:3)) )

!--- 3. GLOBAL U, A due to wave / current
   call UINFLOW_sea ( XG, UG, AG, pdyn, 1 )

!--- 4. LOCAL Uinf, Ainf due to wave / current
   Uinf (1:3)    = matmul ( AT(1:3,1:3), UG(1:3) )
   Ainf (1:3)    = matmul ( AT(1:3,1:3), AG(1:3) )

!--- 5. LLt[3,3] = L_Unit_t.transpose(L_Unit_t)
   do j=1,3
    do i=1,3
       LLt(i,j)  = L_Unit_t(i)*L_Unit_t(j)
    enddo
   enddo

   call DIAGO (3,DIAGmLLt)

   DIAGmLLt = DIAGmLLt - LLt;

!--- 6. elastic/rigid body velocity and acceleration wrt the local c.s.
   Uelst(1:3) =               ATDR (1:3    )                                + &
                     matmul ( ATDA (1:3,1:3), UTT (IMDOF_el(1:3))+Y0(1:3) ) + &
                                              UTT1(IMDOF_el(1:3))
  if (IADDMASSELAST == 0) then
   Aelst(1:3) =               ATDDR(1:3    )                                + &
                     matmul ( ATDDA(1:3,1:3),                    +Y0(1:3) )
  else
   Aelst(1:3) =               ATDDR(1:3    )                                + &
                     matmul ( ATDDA(1:3,1:3), UTT (IMDOF_el(1:3))+Y0(1:3) ) + &
                2.d0*matmul ( ATDA (1:3,1:3), UTT1(IMDOF_el(1:3))         ) + &
                                              UTT2(IMDOF_el(1:3))
  endif

!--- Morison coefficients
   Ur   (1:3)          = Uinf (1:3) - Uelst(1:3)
   Ur_n (1:3)          = matmul( DIAGmLLt(1:3,1:3), Ur(1:3) )
   Ur_n_Meas           = dsqrt ( dot_product ( Ur_n, Ur_n ) )
   CD                  = CDrag       * rho/2d0  * DIAM * Ur_n_Meas
   CM1                 = CMass       * rho * PI * DIAM**2/4d0
   CM2                 =(CMass-1.d0) * rho * PI * DIAM**2/4d0

!--- I  . DRAG term  ---------
   AKL                 =       CD  * matmul( DIAGmLLt, ATDA  )
   ACL                 =       CD  *         DIAGmLLt
   Fdrag(1:3)          =       CD  * Ur_n(1:3)
!--- II . INERTIA term 1 [ Cm*ρ*Vo*An, Vo=pi*R**2 ] -----
   Finer1(1:3)         =            CM1 * matmul( DIAGmLLt(1:3,1:3), Ainf(1:3) )
!--- III. INERTIA term 2 [ -Ca*ρ*Vo* d2q/dt2n ,Ca=Cm-1 ] ------
   AKL                 = AKL + CM2 * matmul( DIAGmLLt, ATDDA )
   ACL                 = ACL + CM2 * matmul( DIAGmLLt, ATDA  ) * 2.d0
   AML                 = AML + CM2 *         DIAGmLLt
   Finer2(1:3)         =     - CM2 * matmul( DIAGmLLt(1:3,1:3), Aelst(1:3) )

   AFL   (1:3)         = Fdrag(1:3) + Finer1(1:3) + Finer2(1:3)

   do iq = 1, QSnode(nn_el)%Tot_el(2)
      AT0    (1:3,1:3) = transf_mat(nn_el)%AT0_el    (iq,1:3,1:3)
      ATDA0  (1:3,1:3) = transf_mat(nn_el)%ATDA0_el  (iq,1:3,1:3)
      ATDA1  (1:3,1:3) = transf_mat(nn_el)%ATDA1_el  (iq,1:3,1:3)
      ATDR0  (1:3    ) = transf_mat(nn_el)%ATDR0_el  (iq,1:3    )
      ATDR1  (1:3    ) = transf_mat(nn_el)%ATDR1_el  (iq,1:3    )
      ATDDA0 (1:3,1:3) = transf_mat(nn_el)%ATDDA0_el (iq,1:3,1:3)
      ATDDA1 (1:3,1:3) = transf_mat(nn_el)%ATDDA1_el (iq,1:3,1:3)
      ATDDA2 (1:3,1:3) = transf_mat(nn_el)%ATDDA2_el (iq,1:3,1:3)
      ATDDR0 (1:3    ) = transf_mat(nn_el)%ATDDR0_el (iq,1:3    )
      ATDDR1 (1:3    ) = transf_mat(nn_el)%ATDDR1_el (iq,1:3    )
      ATDDR2 (1:3    ) = transf_mat(nn_el)%ATDDR2_el (iq,1:3    )
!------ I  . DRAG term  ---------
      AKLnq  (iq, 1:3) =                  CD  * matmul (  DIAGmLLt, ATDR0 (1:3) + matmul( ATDA0 (1:3,1:3), Y0(1:3)+UTT(IMDOF_el(1:3)) ) - matmul(AT0(1:3,1:3), UG(1:3))  )
      ACLnq  (iq, 1:3) =                  CD  * matmul (  DIAGmLLt, ATDR1 (1:3) + matmul( ATDA1 (1:3,1:3), Y0(1:3)+UTT(IMDOF_el(1:3)) )                                  )
!------ II . INERTIA term 1 [ Cm*ρ*Vo*An, Vo=pi*R**2 ] -----
      AKLnq  (iq, 1:3) = AKLnq(iq, 1:3) - CM1 * matmul(  DIAGmLLt, matmul(AT0(1:3,1:3), AG(1:3))  )
!------ III. INERTIA term 2 [ -Ca*ρ*Vo* d2q/dt2n ,Ca=Cm-1 ] ------
      AKLnq  (iq, 1:3) = AKLnq(iq, 1:3) + CM2 * matmul (  DIAGmLLt, ATDDR0(1:3) + matmul( ATDDA0(1:3,1:3), Y0(1:3)+UTT(IMDOF_el(1:3)) ) + 2.d0 * matmul( ATDA0(1:3,1:3), UTT1(IMDOF_el(1:3)) )  )
      ACLnq  (iq, 1:3) = ACLnq(iq, 1:3) + CM2 * matmul (  DIAGmLLt, ATDDR1(1:3) + matmul( ATDDA1(1:3,1:3), Y0(1:3)+UTT(IMDOF_el(1:3)) ) + 2.d0 * matmul( ATDA1(1:3,1:3), UTT1(IMDOF_el(1:3)) )  )
      AMLnq  (iq, 1:3) = AMLnq(iq, 1:3) + CM2 * matmul (  DIAGmLLt, ATDDR2(1:3) + matmul( ATDDA2(1:3,1:3), Y0(1:3)+UTT(IMDOF_el(1:3)) )                                                         )
   enddo

!--- hydro Results
   FG (1:3)                        = matmul ( AT (1:3,1:3), AFL(1:3) )
   Res_hyd (nbod,nel,k_nhyd, 1: 3) = XG    (1:3)
   Res_hyd (nbod,nel,k_nhyd, 4: 6) = UG    (1:3)
   Res_hyd (nbod,nel,k_nhyd, 7: 9) = AG    (1:3)
   Res_hyd (nbod,nel,k_nhyd,10:12) = FDrag (1:3)
   Res_hyd (nbod,nel,k_nhyd,13:15) = FIner1(1:3)
   Res_hyd (nbod,nel,k_nhyd,16:18) = FIner2(1:3)
   Res_hyd (nbod,nel,k_nhyd,19:21) = FG    (1:3)


 END Subroutine MORISON
!-------------------------------------------------------------------------------
!   Subroutine :Buoyancy
!
!  Calculate Local Hydrodynamic Matrices from Buoyancy
!
!-------------------------------------------------------------------------------
 Subroutine Buoyancy_hyd ( AFLOC, AKLOCNQ, nbod, nbsub, nb_el, nn_el, nel )
!-------------------------------------------------------------------------------

 Use Cbeam
 Use Hydro

   implicit none

   integer, intent(in   ) :: nbod, nbsub, nb_el, nn_el, nel
   real(8), intent(inout) :: AKLOCNQ   (NDFPEM_el,NQS      )
   real(8), intent(inout) :: AFLOC     (NDFPEM_el          )

   real(8), dimension(NDFPEM_el) :: UT, UT1, UT2
   real(8) :: UTT       (NEQPEM_el)
   real(8) :: AFLOC0    (NDFPEM_el)
   real(8) :: AKLOCNQ0  (NDFPEM_el)
   real(8) :: AQ        (NEQPEM_el)
   real(8) :: XG        (3)
   real(8) :: SHAPET    (NDFPEM_el,NEQPEM_el)
   real(8) :: SHAPE     (NEQPEM_el,NDFPEM_el)
   real(8) :: RG(3),Y0(3), DL
   real(8) :: FB(6)
   real(8) :: FBKnq(NQS,6)
   integer :: k, n1, n2, NDFPE, NEQPE 
   real(8) :: ALLOC, PA01, PA02, HTA, HTAL, HTA1, HTA0
   real(8) :: DIAM , Fi, RADIN, L1, L2, PHIX, PHIZ
   integer :: iq, n
   integer :: N_low, N_high, Ihalfwet
   real(8) :: Z_calc, Length, Length0, lamda, XG_nod(3,2), Ey(3)
   real(8) :: FG(3), MG(3)


!--- hydrodynamic forces not included for this subbody
   if (subbody(nb_el)%IBTYP_el /= 2) return

   Z_calc  = tide

   HTA0    = subbody(nb_el)%HTA_el    (nel  )
   ALLOC   = subbody(nb_el)%ALENG_el  (nel  )
   n1      = subbody(nb_el)%INODNEL_el(nel,1)
   n2      = subbody(nb_el)%INODNEL_el(nel,2)
   NDFPE   = subbody(nb_el)%NDFPE_el
   NEQPE   = subbody(nb_el)%NEQPE_el
   PHIX    = subbody(nb_el)%PHIX_el   (nel  )
   PHIZ    = subbody(nb_el)%PHIZ_el   (nel  )


   call LOCAL_UT_1   (nb_el, nel, UT, UT1, UT2)


!--- 1st node of the elemenet
   call SSHAPEFUNC15 ( 0.d0, ALLOC, PHIX, PHIZ,&  !-- HTAL=0.d0
                       SHAPE, NDFPE, NEQPE      )

   UTT           = matmul (SHAPE, UT )
   Y0    (1:3  ) = 0.d0;
   Y0    (2    ) = HTA0
!! XG_nod(1:3,1) = transf_mat(nn_el)%R_el(1:3) + matmul ( transf_mat(nn_el)%A_el(1:3,1:3), Y0(1:3) + UTT(IMDOF_el(1:3))  )
   XG_nod(1:3,1) = transf_mat(nn_el)%R_el(1:3) + matmul ( transf_mat(nn_el)%A_el(1:3,1:3), Y0(1:3) )

!--- last node of the elemenet
   call SSHAPEFUNC15 ( 1.d0, ALLOC, PHIX, PHIZ,&   !-- HTAL=1.d0
                       SHAPE, NDFPE, NEQPE       )

   UTT           = matmul (SHAPE, UT )
   Y0    (2    ) = HTA0 + ALLOC
!! XG_nod(1:3,2) = transf_mat(nn_el)%R_el(1:3) + matmul ( transf_mat(nn_el)%A_el(1:3,1:3), Y0(1:3) + UTT(IMDOF_el(1:3)) )
   XG_nod(1:3,2) = transf_mat(nn_el)%R_el(1:3) + matmul ( transf_mat(nn_el)%A_el(1:3,1:3), Y0(1:3) )

!--- check which node is lower
   if (XG_nod(3,2)>=XG_nod(3,1)) then
      N_low  = 1
      N_high = 2
   else
      N_low  = 2
      N_high = 1
   endif

   if (XG_nod(3,N_low ) >   Z_calc) return            !-- both nodes out of water
   if (XG_nod(3,N_high) <= -depth ) return            !-- the element is below mudline (foundation)
                                    Ihalfwet = 0      !-- indicates if the body is half wet, half dry
   if (XG_nod(3,N_high) >   Z_calc) Ihalfwet = 1


   Ey(1:3)  = XG_nod(1:3,N_high) - XG_nod(1:3,N_low)  !-- vector tangential
   Length   = dsqrt(dot_product (Ey,Ey))
   Length0  = Length

   if     (Ihalfwet==0) then
      L1        = 0.d0
      L2        = Length0
   elseif (Ihalfwet==1) then
!------ X-X1=λ(X2-X1) ---> λ = (X(3)-X1(3)) / (X2(3)-X1(3)) 
!------               ---> X(1:2) = X1(1:2) + λ (X2(1:2)-X1(1:2))
      lamda              = (Z_calc-XG_nod(3,N_low)) / Ey(3)
      XG_nod(1:2,N_high) = XG_nod(1:2,N_low) + lamda*Ey(1:2)
      XG_nod(  3,N_high) = Z_calc

      Ey(1:3)   = XG_nod(1:3,N_high) - XG_nod(1:3,N_low)
      Length    = dsqrt(dot_product (Ey,Ey))

      if (N_low==1) then
         L1     = 0.d0
         L2     = Length
      else
         L1     = Length0 - Length
         L2     = Length0
      endif
   endif

   DL = ( L2 - L1 ) / NX_hyd
   Fi = datan2(  (substr_mor(n2)%DIAMET_el-substr_mor(n1)%DIAMET_el), (2.d0*ALLOC)  )

!--- Integrate Buoyancy loading for all the wet side surfaces
   do k        = 1, NX_hyd
!------ Mid point
      HTA      = L1 + dble(2*k-1)/2.d0 * DL
      HTAL     = HTA / ALLOC
      HTA1     = HTA0 + HTA
      PA01     = 1.d0 - HTAL
      PA02     = HTAL

      call SSHAPEFUNC15 ( HTAL, ALLOC, PHIX, PHIZ,& 
                          SHAPE, NDFPE, NEQPE       )

      SHAPET   = transpose (SHAPE)
      DIAM     =  PA01 * substr_mor(n1)%DIAMET_el   +  PA02 * substr_mor(n2)%DIAMET_el
      RADIN    = (PA01 * substr_mor(n1)%INDIAMET_el +  PA02 * substr_mor(n2)%INDIAMET_el)/2.d0
!!!   UTT      = matmul (SHAPE, UT )


      call BuoyancySide ( HTA1, DIAM, RADIN, Fi, FB, FBKnq, nbod, nbsub, nel, k )

!------ Total buoyancy
      FG(1:3)                = matmul(transf_mat(nn_el)%A_el(1:3,1:3), FB(1:3)) * DL
      XG(1:3)                = transf_mat(nn_el)%R_el(1:3) + transf_mat(nn_el)%A_el(1:3,2) * HTA1
      RG(1:3)                = XG(1:3) - R_float  (1:3)

      call EXTEPR ( RG, FG, MG )

      BuoyancyTside_el(1:3)  = BuoyancyTside_el(1:3) + FG(1:3)
      BuoyancyTside_el(4:6)  = BuoyancyTside_el(4:6) + MG(1:3)


!------ Using shape functions, add to local element matrices
      AQ                 = 0.d0;
      AQ (IMDOF_el(1:6)) = FB (1:6)
      AFLOC0             = matmul( SHAPET, AQ )
      AFLOC              = AFLOC + AFLOC0 * DL

!------ for all q's
      do iq = 1, QSnode(nn_el)%Tot_el(2)
         n  = QSnode(nn_el)%IORD_el(iq)

         AQ                 = 0.d0;
         AQ (IMDOF_el(1:6)) = FBKnq (iq,1:6)
         AKLOCNQ0           = matmul( SHAPET, AQ )
         AKLOCNQ(1:NDFPE,n) = AKLOCNQ(1:NDFPE,n) + AKLOCNQ0(1:NDFPE) * DL
      enddo !iq
   enddo !k


 END Subroutine Buoyancy_hyd
!----------------------------------------------------------------------
 Subroutine BuoyancySide ( HTA1, DIAM, RADIN, Fi, FB, FBKnq, nbod, nbsub, nel, k_nhyd )
!----------------------------------------------------------------------

 Use Hydro
 Use Cbeam

   implicit none

   integer, intent(in ) :: nbod, nbsub, nel, k_nhyd
   real(8), intent(in ) :: HTA1, DIAM, RADIN, Fi
   real(8), intent(out) :: FBKnq(NQS,6)
   real(8), intent(out) :: FB   (    6)

   real(8) :: XG(3), UG(3  ), AG(3)
   real(8) :: R (3), A (3,3)
   real(8) :: R0(3), A0(3,3)        
   real(8) :: Fct1 , Fct2  , Fct3
   real(8) :: area , perim , areaR
   real(8) :: Z0g  , pdyn
   integer :: nn_el, iq    , idyn


!--- enable/disable the dynamic presure term in case of cone members
   idyn = 1
!--- Zero matrices
   FB   (      1:6) = 0.d0;
   FBKnq(1:NQS,1:6) = 0.d0;

   if (IBUOYANCY == 0) return

!--- Rotation Matrix
   nn_el        = body(nbod)%NNODTGLB_el (nbsub)
   A  (1:3,1:3) = transf_mat(nn_el)%A_el(1:3,1:3)
   R  (    1:3) = transf_mat(nn_el)%R_el(    1:3)
!--- Zoglob
   XG (    1:3) = R(1:3) + A(1:3,2) * HTA1
   Z0g          = XG(3) - tide

   if     (IBUOYANCY == 1) then
!------ Presure integration method
      pdyn = 0.d0
      if ((Fi/=0.d0).and.(idyn==1)) then
         call UINFLOW_sea ( XG, UG, AG, pdyn, 1 )

        !pdyn = 0.5d0*rho*dot_product(UG,UG)
      endif

      if (IFloodbod_el(nbod) == 1) then
         area  = (DIAM**2 /4d0 -     RADIN**2) * PI
         perim = (DIAM         - 2d0*RADIN   ) * PI
         areaR = (DIAM**3 /8d0 -     RADIN**3) * PI

         Fct1  =  rhoGrav      * area  * dcos(Fi)
         Fct2  = -rhoGrav      * perim * dsin(Fi)                    !-- axial     force  due to cone
         Fct3  =  rhoGrav      * areaR * dsin(Fi)                    !-- restoring moment due to cone
         perim =  DIAM         * PI                                  !-- for dynamic presure only the outside water contributes
      else
         area  =  DIAM**2 /4d0 * PI
         perim =  DIAM         * PI
         areaR =  area * DIAM/2d0

         Fct1  =  rhoGrav      * area  * dcos(Fi)
         Fct2  = -rhoGrav      * perim * dsin(Fi)                    !-- axial     force  due to cone
         Fct3  =  rhoGrav      * areaR * dsin(Fi)                    !-- restoring moment due to cone
      endif

      FB(1) =  Fct1 * A(3,1)
      FB(2) =  Fct2 * Z0g              &                             !-- axial     force  due to cone
              -pdyn * perim * dsin(Fi)
      FB(3) =  Fct1 * A(3,3)
      FB(4) =  Fct3 * A(3,3)                                         !-- restoring moment due to cone
      FB(6) = -Fct3 * A(3,1)                                         !-- restoring moment due to cone

      do iq = 1, QSnode(nn_el)%Tot_el(2)
         R0  (    1:3) = transf_mat(nn_el)%R0_el (iq,    1:3)
         A0  (1:3,1:3) = transf_mat(nn_el)%A0_el (iq,1:3,1:3)
!--------- (-) because terms are moved to LHS
         FBKnq (iq ,1) = -Fct1 *  A0(3,1)
         FBKnq (iq ,2) = -Fct2 * (A0(3,2)*HTA1+R0(3))                !-- axial     force  due to cone
         FBKnq (iq ,3) = -Fct1 *  A0(3,3)
         FBKnq (iq ,4) = -Fct3 *  A0(3,3)                            !-- restoring moment due to cone
         FBKnq (iq ,6) =  Fct3 *  A0(3,1)                            !-- restoring moment due to cone
      enddo

   elseif (IBUOYANCY == 2) then
!------ Volume method
!------ Atm no conicity is taken into account
      if (IFloodbod_el(nbod) == 1) then
         Fct1  =  ( DIAM**2/4d0 - RADIN**2 ) * PI * rhoGrav
      else
         Fct1  =    DIAM**2/4d0              * PI * rhoGrav
      endif

      FB(1:3)  =  Fct1 * A (3,1:3) !!AT(1:3,3)

      do iq = 1, QSnode(nn_el)%Tot_el(2)
         A0    (1:3,1:3) = transf_mat(nn_el)%A0_el (iq,1:3,1:3)

         FBKnq (iq ,1:3) = -Fct1 * A0(3,1:3)
      enddo
   endif

!--- hydro Results
   Res_hyd (nbod,nel,k_nhyd,22:24) = XG (1:3)
   Res_hyd (nbod,nel,k_nhyd,25:30) = FB (1:6)


 END Subroutine BuoyancySide
!----------------------------------------------------------------------
 Subroutine BuoyancyTap
!----------------------------------------------------------------------

!--- Buoyancy Taps defined by the user


!--- MSL for rho g z, stretch for morison and dyn pres
!--- Atm the search doen't consider the z_elevation which means that
!--- if the tap is below z_calc_low the stretching is applied,
!--- but it the taper is near free surface only the MSL calculation is consistent.
!--- morison term2


 Use Hydro
 Use Cbeam

   implicit none

   real(8) :: ALLOC, FL
   real(8) :: n_unit, HTAL, HTA0, PHIX, PHIZ, area, ca,cd
   integer :: nbod,nbsub,nb_el,nn_el,nel, itap,nod, i,nn, NDFPE, NEQPE, ityp, idyn, irel, ifroude
   real(8) :: AQ      (NEQPEM_el          )
   real(8) :: SHAPET  (NDFPEM_el,NEQPEM_el)
   real(8) :: SHAPE   (NEQPEM_el,NDFPEM_el)
   real(8), dimension (NDFPEM_el,NQS      ) :: AKLOCNQ , ACLOCNQ , AMLOCNQ
   real(8), dimension (NDFPEM_el          ) :: AKLOCNQ0, ACLOCNQ0, AMLOCNQ0
   real(8) :: AFLOC0  (NDFPEM_el          )
   real(8) :: FLKnq(NQS)
   real(8) :: Fhyd, Fdyn, Fdrag, Faddmass, Finer, Frel
   real(8) :: Fct
   integer :: iq, nq, j
   real(8) :: FG(3), MG(3), RG(3), Z0g
   real(8) :: XG(3), UG(3), AG(3), pdyn
!--- for Morison
   real(8) :: Y0    (3)
   real(8) :: Uinf  (3), Ainf  (3)
   real(8) :: Uelst (3), Aelst (3)
   real(8) :: Ur    (3), Ur_n     , Ur_n_Meas
   real(8) :: Ar    (3), Ar_n     , A_n
   real(8) :: CD_m     , CM_m

   real(8), dimension(3,3) :: ATDA , ATDDA , A     , AT
   real(8), dimension(3  ) :: ATDR , ATDDR , R
   real(8), dimension(3,3) :: ATDA0, ATDA1 , ATDDA0, ATDDA1, ATDDA2, A0, AT0
   real(8), dimension(3  ) :: ATDR0, ATDR1 , ATDDR0, ATDDR1, ATDDR2, R0

!! real(8), dimension(3  ) :: AKL, ACL, AML
   real(8), dimension(NQS) :: AKLnq, ACLnq, AMLnq
   real(8) :: AFL


   do itap    = 1, NBuoyaTap_el
      nbod    = buoyancy_tap(itap)%nbod
      nel     = buoyancy_tap(itap)%nel
      nod     = buoyancy_tap(itap)%nod             !-- 1 or NNPE_el
      ityp    = buoyancy_tap(itap)%itype
      area    = buoyancy_tap(itap)%area
      ca      = buoyancy_tap(itap)%ca
      cd      = buoyancy_tap(itap)%cd
      idyn    = buoyancy_tap(itap)%idyn_pres
      irel    = buoyancy_tap(itap)%irel            !-- enable/disable the relative  term from Bernoulli eq.
      ifroude = buoyancy_tap(itap)%ifroude         !-- enable/disable Froude-Krylov term from Morison's eq.
      nbsub   = 1                                  !-- always for jacket bodies 1 sub-body!
      nb_el   = body(nbod)%NBODTGLB_el(nbsub)
      nn_el   = body(nbod)%NNODTGLB_el(nbsub)

      if (subbody(nb_el)%IBTYP_el/= 2) cycle     !-- hydrodynamic forces not included for this subbody

      if (nod == 1) then
         HTA0   = subbody(nb_el)%HTA_el(nel  )
         n_unit =-1.d0                           !-- local unit vector pointing the outside flow (0,-1,0)
         HTAL   = 0.d0
      else !nod==NNPE_el
         HTA0   = subbody(nb_el)%HTA_el(nel+1)
         n_unit = 1.d0                           !-- local unit vector pointing the outside flow (0, 1,0)
         HTAL   = 1.d0
      endif

      if (ityp==2) n_unit = -n_unit              !-- buoyancy is opposite for flooded members [ityp=2]


      ALLOC = subbody(nb_el)%ALENG_el(nel)
      NDFPE = subbody(nb_el)%NDFPE_el
      NEQPE = subbody(nb_el)%NEQPE_el
      PHIX  = subbody(nb_el)%PHIX_el (nel)
      PHIZ  = subbody(nb_el)%PHIZ_el (nel)

!------ Rotation Matrix
      A  (1:3,1:3) = transf_mat(nn_el)%A_el(1:3,1:3)
      R  (    1:3) = transf_mat(nn_el)%R_el(    1:3)
!------ Zoglob
      XG (    1:3) = R(1:3) + A(1:3,2) * HTA0
!!    XG (    1:3) = R(1:3) + matmul ( A(1:3,1:3), Y0(1:3)+UTT(IMDOF_el(1:3)) )
      Z0g          = XG(3) - tide
      FL           = 0.d0

!qq Here we should also check the z elevation, in order to be consistent with the instantaneous position and wheeler streatching...
!qq Also the check could be splited in 2 parts: a. buoyancy [MSL check], b, Morison [MSL or IW, based of Streatch option]
      if (Z0g>=0.d0) then
!--------- taper is out of the water
         Fhyd      = 0.d0
         Fdyn      = 0.d0
         Frel      = 0.d0
         Fdrag     = 0.d0
         Finer     = 0.d0
         Faddmass  = 0.d0
         goto 2
      endif


!------ Morison's equation for Plates [drag, inertia]
!!    call LOCAL_UT_1   (nb_el, nel, UT, UT1, UT2)

!!    call SSHAPEFUNC15 ( HTAL, ALLOC, PHIX, PHIZ,& 
!!                        SHAPE, NDFPE, NEQPE       )

!!    SHAPET   = transpose (SHAPE )

!!    UTT      = matmul (SHAPE, UT )
!!    UTT1     = matmul (SHAPE, UT1)
!!    UTT2     = matmul (SHAPE, UT2)

!------ Pass Transformation matrices to local arrays
      AT    (1:3,1:3) = transf_mat(nn_el)%AT_el    (1:3,1:3)
      ATDA  (1:3,1:3) = transf_mat(nn_el)%ATDA_el  (1:3,1:3)
      ATDDA (1:3,1:3) = transf_mat(nn_el)%ATDDA_el (1:3,1:3)
      ATDR  (1:3    ) = transf_mat(nn_el)%ATDR_el  (1:3    )
      ATDDR (1:3    ) = transf_mat(nn_el)%ATDDR_el (1:3    )


!------ 1. XG      : global coordinates of the evaluation point
      Y0 (1:3) = 0.d0;
      Y0 (2  ) = HTA0

!------ 2. GLOBAL U, A due to wave / current
      call UINFLOW_sea ( XG, UG, AG, pdyn, 1 )

!------ 3. LOCAL Uinf, Ainf due to wave / current
      Uinf (1:3)    = matmul ( AT(1:3,1:3), UG(1:3) )
      Ainf (1:3)    = matmul ( AT(1:3,1:3), AG(1:3) )


!------ 4. rigid body velocity and acceleration wrt the local c.s. [elastic velocity is not considered]
      Uelst(1:3) =               ATDR (1:3    )            + &
                        matmul ( ATDA (1:3,1:3), Y0(1:3) )

      Aelst(1:3) =               ATDDR(1:3    )            + &
                        matmul ( ATDDA(1:3,1:3), Y0(1:3) )

!!    Uelst(1:3) =               ATDR (1:3    )                                + &
!!                      matmul ( ATDA (1:3,1:3), UTT (IMDOF_el(1:3))+Y0(1:3) ) + &
!!                                               UTT1(IMDOF_el(1:3))
!!    Aelst(1:3) =               ATDDR(1:3    )                                + &
!!                      matmul ( ATDDA(1:3,1:3), UTT (IMDOF_el(1:3))+Y0(1:3) ) + &
!!                 2.d0*matmul ( ATDA (1:3,1:3), UTT1(IMDOF_el(1:3))         ) + &
!!                                               UTT2(IMDOF_el(1:3))

!------ Zero matrices
!!    AKL   = 0.d0;
!!    ACL   = 0.d0;
!!    AML   = 0.d0;
      AFL   = 0.d0
      AKLnq = 0.d0;
      ACLnq = 0.d0;
      AMLnq = 0.d0;

!------ I . DRAG term [ ρ/2*Cd*S*|Ur_n| * *Ur_n ] ---------
      Ur   (1:3)          = Uinf (1:3) - Uelst(1:3)
      Ur_n                = Ur(2)                                         !-- only the y component
      Ur_n_Meas           = dabs(Ur_n)
      CD_m                = area * cd * 0.5d0 * rho * Ur_n_Meas
      AFL                 = CD_m * Ur_n
!!    AKL  (1:3)          = CD_m * ATDA(2,1:3)
!!    ACL  (1:3)          = CD_m

      do iq = 1, QSnode(nn_el)%Tot_el(2)
         AT0    (1:3,1:3) = transf_mat(nn_el)%AT0_el    (iq,1:3,1:3)
         ATDA0  (1:3,1:3) = transf_mat(nn_el)%ATDA0_el  (iq,1:3,1:3)
         ATDA1  (1:3,1:3) = transf_mat(nn_el)%ATDA1_el  (iq,1:3,1:3)
         ATDR0  (1:3    ) = transf_mat(nn_el)%ATDR0_el  (iq,1:3    )
         ATDR1  (1:3    ) = transf_mat(nn_el)%ATDR1_el  (iq,1:3    )

!!       AKLnq  (iq) = CD_m * ( ATDR0(2) + dot_product( ATDA0(2,1:3), Y0(1:3)+UTT(IMDOF_el(1:3)) ) - dot_product(AT0(2,1:3), UG(1:3))  )
!!       ACLnq  (iq) = CD_m * ( ATDR1(2) + dot_product( ATDA1(2,1:3), Y0(1:3)+UTT(IMDOF_el(1:3)) )                                     )

         AKLnq  (iq) = CD_m * ( ATDR0(2) + dot_product( ATDA0(2,1:3), Y0(1:3) ) - dot_product( AT0(2,1:3), UG(1:3) ) )
         ACLnq  (iq) = CD_m * ( ATDR1(2) + dot_product( ATDA1(2,1:3), Y0(1:3) )                                      )
      enddo

      Fdrag               = CD_m * Ur_n
      Finer               = 0.d0
      Faddmass            = 0.d0

      if (ifroude /= 0) then
!--------- II . FROUDE-KRYLOV term [ ρ*Vo*An, Vo=pi*R**2 ] -----
         A_n                 = Ainf(2)
         CM_m                = area * rho
         Finer               =       CM_m * A_n
         AFL                 = AFL + CM_m * A_n

         do iq = 1, QSnode(nn_el)%Tot_el(2)
            AT0    (1:3,1:3) = transf_mat(nn_el)%AT0_el (iq,1:3,1:3)

            AKLnq  (iq)      = AKLnq(iq) - CM_m * dot_product( AT0(2,1:3), AG(1:3) )
         enddo
      endif

!------ III. ADDED MASS term [ ρ*Ca*Vo* Ar_n, Ca=A33(0)/(ρ*Vo), Vo=S->area ] ------
      Ar   (1:3)          = Ainf (1:3) - Aelst(1:3)
      Ar_n                = Ar(2)                                         !-- only the y component
      CM_m                = area * ca * rho
      Faddmass            = CM_m * Ar_n
      AFL                 = AFL + CM_m * Ar_n
!!    AKL  (1:3)          = AKL(1:3) + CM_m *        ATDDA(2,1:3)
!!    ACL  (1:3)          = ACL(1:3) + CM_m * 2.d0 * ATDA (2,1:3)
!!    AML  (1:3)          = AML(1:3) + CM_m

      do iq = 1, QSnode(nn_el)%Tot_el(2)
         AT0    (1:3,1:3) = transf_mat(nn_el)%AT0_el    (iq,1:3,1:3)
         ATDA0  (1:3,1:3) = transf_mat(nn_el)%ATDA0_el  (iq,1:3,1:3)
         ATDA1  (1:3,1:3) = transf_mat(nn_el)%ATDA1_el  (iq,1:3,1:3)
         ATDDA0 (1:3,1:3) = transf_mat(nn_el)%ATDDA0_el (iq,1:3,1:3)
         ATDDA1 (1:3,1:3) = transf_mat(nn_el)%ATDDA1_el (iq,1:3,1:3)
         ATDDA2 (1:3,1:3) = transf_mat(nn_el)%ATDDA2_el (iq,1:3,1:3)
         ATDDR0 (1:3    ) = transf_mat(nn_el)%ATDDR0_el (iq,1:3    )
         ATDDR1 (1:3    ) = transf_mat(nn_el)%ATDDR1_el (iq,1:3    )
         ATDDR2 (1:3    ) = transf_mat(nn_el)%ATDDR2_el (iq,1:3    )

         AKLnq  (iq) = AKLnq(iq) + CM_m * ( ATDDR0(2) + dot_product( ATDDA0(2,1:3), Y0(1:3) ) - dot_product( AT0(2,1:3), AG(1:3) ) )
         ACLnq  (iq) = ACLnq(iq) + CM_m * ( ATDDR1(2) + dot_product( ATDDA1(2,1:3), Y0(1:3) )                                      )
         AMLnq  (iq) = AMLnq(iq) + CM_m * ( ATDDR2(2) + dot_product( ATDDA2(2,1:3), Y0(1:3) )                                      )
!!       AKLnq  (iq) = AKLnq(iq) + CM_m * ( ATDDR0(2) + dot_product( ATDDA0(2,1:3), Y0(1:3)+UTT(IMDOF_el(1:3)) ) + 2.d0 * dot_product( ATDA0(2,1:3), UTT1(IMDOF_el(1:3)) )  )
!!       ACLnq  (iq) = ACLnq(iq) + CM_m * ( ATDDR1(2) + dot_product( ATDDA1(2,1:3), Y0(1:3)+UTT(IMDOF_el(1:3)) ) + 2.d0 * dot_product( ATDA1(2,1:3), UTT1(IMDOF_el(1:3)) )  )
!!       AMLnq  (iq) = AMLnq(iq) + CM_m * ( ATDDR2(2) + dot_product( ATDDA2(2,1:3), Y0(1:3)+UTT(IMDOF_el(1:3)) )                                                         )
      enddo

!-------------------------------------------------------------------------------------------------------------
! p = -rho*g*z - rho*dfi/dt =-rho*g*z + pdyn [for nonlinear theories 0.5grad{fi}^2 is incorporated in pdyn]
! F = -Int{p n ds}, n:points outwards
! F = n*rho*g*z*S -n*pdyn*S
!-------------------------------------------------------------------------------------------------------------
!------ Local Buoyancy Force applied on current tap
      Fct   = n_unit * area * rhoGrav
      Fhyd  = Fct * Z0g
      Fdyn  = 0.d0
      Frel  = 0.d0
      if (irel==1) &
      Frel  =-n_unit * area         * rho * dot_product(UG,Uelst)  !-- Bernoulli equ. term for the relative system
      if (idyn==1) &
      Fdyn  =-n_unit * area * pdyn
!     Fdyn  =-n_unit * area * 0.5d0 * rho * dot_product(UG,UG)
      FL    = Fhyd + Fdyn + Frel

      do iq = 1, QSnode(nn_el)%Tot_el(2)
         R0    (    1:3) = transf_mat(nn_el)%R0_el (iq,    1:3)
         A0    (1:3,1:3) = transf_mat(nn_el)%A0_el (iq,1:3,1:3)

         FLKnq (iq) = - Fct * (R0(3)+A0(3,2)*HTA0)                 !-- (-) to LHS
      enddo


!------ Total buoyancy
      FG(1:3)               = A (1:3,2)*FL
      RG(1:3)               = XG(1:3) - R_float(1:3)

      call EXTEPR ( RG, FG, MG )

      BuoyancyTtap_el(1:3)  = BuoyancyTtap_el(1:3) + FG(1:3)
      BuoyancyTtap_el(4:6)  = BuoyancyTtap_el(4:6) + MG(1:3)


!------ Add to global matrices
      call SSHAPEFUNC15 ( HTAL, ALLOC, PHIX, PHIZ,& 
                          SHAPE, NDFPE, NEQPE       )

!! ------ also  -------------------
!     AKLOC0
!     ACLOC0
!     AMLOC0
      SHAPET           = transpose ( SHAPE      )
      AQ (       1:6 ) = 0.d0;
      AQ (IMDOF_el(2)) = FL + AFL
      AFLOC0           = matmul    ( SHAPET, AQ )
      nn               =                                   subbody(nb_el)%NODMB       (nel, 1)!nnod)
      i                = ACCU%NDFTACC_el ( nb_el-1 )   +   subbody(nb_el)%NDFPNBACC_el(nn-1  )

!--------- for all q's
         AKLOCNQ(1:NDFPE,:) = 0.d0;
         ACLOCNQ(1:NDFPE,:) = 0.d0;
         AMLOCNQ(1:NDFPE,:) = 0.d0;

      do iq = 1, QSnode(nn_el)%Tot_el(2)
         nq = QSnode(nn_el)%IORD_el(iq)

         AQ (       1:6 )    = 0.d0;
         AQ (IMDOF_el(2))    = FLKnq (iq) + AKLnq (iq)
         AKLOCNQ0            = matmul(SHAPET, AQ )
         AQ (       1:6 )    = 0.d0;
         AQ (IMDOF_el(2))    = ACLnq (iq)
         ACLOCNQ0            = matmul(SHAPET, AQ )
         AQ (       1:6 )    = 0.d0;
         AQ (IMDOF_el(2))    = AMLnq (iq)
         AMLOCNQ0            = matmul(SHAPET, AQ )
!--------- Replace
         AKLOCNQ(1:NDFPE,nq) = AKLOCNQ0(1:NDFPE)
         ACLOCNQ(1:NDFPE,nq) = ACLOCNQ0(1:NDFPE)
         AMLOCNQ(1:NDFPE,nq) = AMLOCNQ0(1:NDFPE)
      enddo !iq


!--------- STORE_LOC_MATRIX_el
!!       subbody(nb_el)%AKLOC_el   ( nel, 1:NDFPE, 1:NDFPE) = subbody(nb_el)%AKLOC_el   ( nel, 1:NDFPE, 1:NDFPE) + &
!!                                                                           AKLOC0     (      1:NDFPE, 1:NDFPE)
!!       subbody(nb_el)%ACLOC_el   ( nel, 1:NDFPE, 1:NDFPE) = subbody(nb_el)%ACLOC_el   ( nel, 1:NDFPE, 1:NDFPE) + &
!!                                                                           ACLOC0     (      1:NDFPE, 1:NDFPE)
!!       subbody(nb_el)%AMLOC_el   ( nel, 1:NDFPE, 1:NDFPE) = subbody(nb_el)%AMLOC_el   ( nel, 1:NDFPE, 1:NDFPE) + &
!!                                                                           AMLOC0     (      1:NDFPE, 1:NDFPE)
         subbody(nb_el)%AFLOC_el   ( nel, 1:NDFPE         ) = subbody(nb_el)%AFLOC_el   ( nel, 1:NDFPE         ) + &
                                                                             AFLOC0     (      1:NDFPE         )
      do iq = 1, QSnode(nn_el)%Tot_el(2)
         nq = QSnode(nn_el)%IORD_el(iq)

         subbody(nb_el)%AKLOCNQ_el ( nel, 1:NDFPEM_el, iq ) = subbody(nb_el)%AKLOCNQ_el ( nel, 1:NDFPEM_el, iq ) + &
                                                                             AKLOCNQ    (      1:NDFPEM_el, nq )
         subbody(nb_el)%ACLOCNQ_el ( nel, 1:NDFPEM_el, iq ) = subbody(nb_el)%ACLOCNQ_el ( nel, 1:NDFPEM_el, iq ) + &
                                                                             ACLOCNQ    (      1:NDFPEM_el, nq )
         subbody(nb_el)%AMLOCNQ_el ( nel, 1:NDFPEM_el, iq ) = subbody(nb_el)%AMLOCNQ_el ( nel, 1:NDFPEM_el, iq ) + &
                                                                             AMLOCNQ    (      1:NDFPEM_el, nq )
      enddo !iq

!--------- ASSEMBLY_MATRIX_el
!! ------- also -------------------
!        AKLOC0
!        ACLOC0
!        AMLOC0 
         AQ_el( i+1: i+NDFPE   ) = AQ_el( i+1:i+NDFPE    ) + AFLOC0  ( 1:NDFPE     )

      do iq = 1, QSnode(nn_el)%Tot_el(2)
         nq = QSnode(nn_el)%IORD_el(iq)
         j  = NDFBT_el + nq

         AK_el( i+1:i+NDFPE, j ) = AK_el( i+1:i+NDFPE, j ) + AKLOCNQ ( 1:NDFPE, nq )
         AC_el( i+1:i+NDFPE, j ) = AC_el( i+1:i+NDFPE, j ) + ACLOCNQ ( 1:NDFPE, nq )
         AM_el( i+1:i+NDFPE, j ) = AM_el( i+1:i+NDFPE, j ) + AMLOCNQ ( 1:NDFPE, nq )
      enddo !iq

!------ Store results for output
  2   buoyancy_tap(itap)%Output_res(1) = Fhyd
      buoyancy_tap(itap)%Output_res(2) = Fdyn
      buoyancy_tap(itap)%Output_res(3) = Frel
      buoyancy_tap(itap)%Output_res(4) = Fdrag
      buoyancy_tap(itap)%Output_res(5) = Finer
      buoyancy_tap(itap)%Output_res(6) = Faddmass
      buoyancy_tap(itap)%Output_res(7) = Z0g
   enddo !itap


 END Subroutine BuoyancyTap
!----------------------------------------------------------------------
 Subroutine Writeout_MORISON
!----------------------------------------------------------------------

 Use Hydro
 Use Cbeam

   implicit none

   integer      :: nbod, nb_el, nel, k
   character    :: CNUM1(3), CNUM2(3), CNUM3(3)
   character*80 :: outfil


   if ( (ICASE_el/=2).and.(ICASE_el/=4) ) return


   do    nbod = NBODBTWT_el+1, NBODBT_el
                                          call INT_2_CHAR ( 3, CNUM1, nbod )
         nb_el= body(nbod)%NBODTGLB_el(1)
      do nel  = 1, subbody(nb_el)%NTEB_el
                                          call INT_2_CHAR ( 3, CNUM2, nel  )
         do k = 1, NX_hyd
                                          call INT_2_CHAR ( 3, CNUM3, k    )

#    if   ASCII == 1
            outfil = 'hydro_morison'//CNUM1(1)//CNUM1(2)//CNUM1(3)//'_'//CNUM2(1)//CNUM2(2)//CNUM2(3)//'_'//CNUM3(1)//CNUM3(2)//CNUM3(3)//'.dat'
            open (201, file=outfil,access='append')
              write(201,101)                     &
#    elif ASCII == 0
            outfil = 'hydro_morison'//CNUM1(1)//CNUM1(2)//CNUM1(3)//'_'//CNUM2(1)//CNUM2(2)//CNUM2(3)//'_'//CNUM3(1)//CNUM3(2)//CNUM3(3)//'.bin'
            open (201, file=outfil,access='append',form='UNFORMATTED')
              write(201    )                     &
#    endif
               sngl(TIME                      ), & !              1
               sngl(Res_hyd (nbod,nel,k, 1: 3)), & ! XG    (1:3)  2-4
               sngl(Res_hyd (nbod,nel,k, 4: 6)), & ! UG    (1:3)  5-7
               sngl(Res_hyd (nbod,nel,k, 7: 9)), & ! AG    (1:3)  8-10
               sngl(Res_hyd (nbod,nel,k,10:12)), & ! FDrag (1:3) 11-13
               sngl(Res_hyd (nbod,nel,k,13:15)), & ! FIner1(1:3) 14-16
               sngl(Res_hyd (nbod,nel,k,16:18)), & ! FIner2(1:3) 17-19
               sngl(Res_hyd (nbod,nel,k,19:21)), & ! FG    (1:3) 20-22
               sngl(Res_hyd (nbod,nel,k,22:24)), & ! XG    (1:3) 23-25
               sngl(Res_hyd (nbod,nel,k,25:30)), & ! FB    (1:6) 26-31
               sngl(Res_hyd (nbod,nel,1,31   ))    ! Length      32
            close(201)

         enddo !i
      enddo !nel
   enddo !nbod

101 format ( 150f15.4 )


 END Subroutine Writeout_MORISON
!
!
!
!-------------------------------------------------------------------------
!   Subroutine :FLOODMATRLOC
! 
!  Local Inertia Matrix with linearly varying properties is calculated 
!  with respect to the body FEM local system for the free flooded members.
!
!  Rem: For jacket's members ATDA=ATDDA=ATDDR=0 and nqs=0 --> only AMLOC [15x15]
!
!-------------------------------------------------------------------------
 Subroutine FLOODMATRLOC (  AMLOC, ACLOC, AKLOC, AFLOC, AMLOCNQ, ACLOCNQ, AKLOCNQ,&
                            nbod , nbsub, nb_el, nn_el, nel                         )
!-------------------------------------------------------------------------

 Use Cbeam
 Use Hydro

   implicit none

   integer,                                 intent(in   ) :: nbod, nbsub, nb_el, nn_el, nel
   real(8), dimension(NDFPEM_el,NDFPEM_el), intent(inout) :: AKLOC  , ACLOC  , AMLOC
   real(8), dimension(NDFPEM_el,NQS      ), intent(inout) :: AKLOCNQ, ACLOCNQ, AMLOCNQ
   real(8), dimension(NDFPEM_el          ), intent(inout) :: AFLOC

   real(8), dimension(NDFPEM_el,NDFPEM_el) :: AMLOC0 !j, ACLOC0, AKLOC0
   real(8), dimension(NDFPEM_el)           :: UT    , UT1  , UT2
   real(8), dimension(NEQPEM_el,NDFPEM_el) :: SHAPE , DSHAPE 
   real(8), dimension(NDFPEM_el,NEQPEM_el) :: SHAPET, DSHAPET, SHAPETAM
!j real(8), dimension(NEQPEM_el)           :: UTT0  , UTT1
!j real(8) :: AT(3,3), ATDAx2(3,3), ATDDA(3,3), ATDDR(3)
   real(8) :: AFLOC0  (NDFPEM_el )         
   real(8) :: AM      (NEQPEM_el,NEQPEM_el)
!j real(8) :: AQ      (NEQPEM_el          )
   integer :: NEQPE, NDFPE, k, n1, n2
   real(8) :: ALLOC, PA01 , PA02 , allocwgauss,         L1  , L2  , HTA
   real(8) :: DENS , AMOMX, AMOMZ, RIXX       , RIZZ  , RIXZ, POLI, HTAL
   real(8) :: PHPX , PHPZ , HTA1 , HTA0       ,               R   , Z0
!j real(8) :: F_Grav,XCM  , ZCM
   real(8) :: RG(3), Rloc(3)


   if (subbody(nb_el)%IBTYP_el /= 2) return     !hydrodynamic forces not included for this subbody
   if (ICASE_el                /= 2) return     !only jacket case has flooded members atm
   if (IFloodbod_el(nbod)      /= 1) return     !if member not free flooded


!j ATDAx2(1:3,1:3) = 2.d0*transf_mat(nn_el)%ATDA_el (1:3,1:3)
!j ATDDA (1:3,1:3) =      transf_mat(nn_el)%ATDDA_el(1:3,1:3)
!j AT    (1:3,1:3) =      transf_mat(nn_el)%AT_el   (1:3,1:3)
!j ATDDR (1:3    ) =      transf_mat(nn_el)%ATDDR_el(1:3    )


   HTA0    = subbody(nb_el)%HTA_el (nel)

!--- check position wrt MSL
   Z0      = transf_mat(nn_el)%R_el(3) + transf_mat(nn_el)%A_el(3,2) * HTA0

   if (Z0 >= tide) return     !the MSL is below the 1st nod of current element

   HTA     = subbody(nb_el)%HTA_el (nel+1)
   Z0      = transf_mat(nn_el)%R_el(3) + transf_mat(nn_el)%A_el(3,2) * HTA

!--- set L1, L2
   if (Z0 <= tide) then

      L2   = subbody(nb_el)%ALENG_el (nel)

   else

      Rloc (1:3) = 0.d0
      Rloc (2)   = HTA0
      RG(3)      = tide

      call Calc_Yloc ( RG, Rloc, L2, nn_el )

   endif

   L1      = 0.d0

   n1      = subbody(nb_el)%INODNEL_el(nel,1)
   n2      = subbody(nb_el)%INODNEL_el(nel,2)
   NEQPE   = subbody(nb_el)%NEQPE_el
   NDFPE   = subbody(nb_el)%NDFPE_el
   ALLOC   = subbody(nb_el)%ALENG_el  (nel  )
   PHPX    = subbody(nb_el)%PHIX_el   (nel  )
   PHPZ    = subbody(nb_el)%PHIZ_el   (nel  )


   call LOCAL_UT_1 (nb_el, nel, UT, UT1, UT2)

!----------------------------!
! IIxAxS           : set AM  !
! IIxAtxF_Grav     : set AQ  !
! IIxAtDDR_AtDDAr0 : sum AQ  !
!----------------------------!

   do k        = 1, IORDER
      HTA      = ( (L2-L1)*XGAUSS(k) + L1+L2 )/2.d0    ! 0 < HTA    < ALLOC_Force
      HTAL     = HTA / ALLOC                           ! 0 < HTAL   < 1           , HTAL=HTA/ALLOC
      HTA1     = HTA0 + HTA

   !!!HTA      = ( XGAUSS(k) + 1.d0 )/2.d0 * ALLOC_Force   ! 0 < HTA    < ALLOC_Force
   !!!HTAL     = HTA / ALLOC                               ! 0 < HTAL   < 1           , HTAL=HTA/ALLOC
   !!!HTA1     = HTA0 + HTA

      PA01     = 1.d0 - HTAL
      PA02     = HTAL

      call SHAPEFUNC15 ( HTAL, ALLOC, PHPX, PHPZ, SHAPE, DSHAPE, NDFPE, NEQPE )

!j    UTT1     = matmul ( SHAPE , UT1 )
!j    UTT0     = matmul ( SHAPE , UT  )

      SHAPET   = transpose (SHAPE )
      DSHAPET  = transpose (DSHAPE)

      R        = (PA01 * substr_mor(n1)%INDIAMET_el  +  PA02 * substr_mor(n2)%INDIAMET_el) /2.d0
      DENS     = PI * R**2 * rho
      AMOMX    = 0.d0
      AMOMZ    = 0.d0
      RIXX     = PI/4d0 * R**4 * rho
      RIZZ     = RIXX
      RIXZ     = 0.d0
      POLI     = RIXX + RIZZ

!j    F_Grav   = 0.d0 !pass to matrloc_q
!j    XCM      = 0.d0
!j    ZCM      = 0.d0

!-------------------------------------------
!  calculation of the local MASS matrix [M]
!
!
!  Nt*(rdA)*II*S*N*u
!-------------------------------------------
      AM       = 0.d0;
      AM(1,1)  = DENS
!     AM(6,1)  = AMOMX
      AM(2,2)  = RIZZ 
!     AM(3,2)  = AMOMZ
      AM(5,2)  =-RIXZ
      AM(2,3)  = AMOMZ
!     AM(3,3)  = DENS
      AM(5,3)  =-AMOMX
      AM(4,4)  = DENS
!     AM(6,4)  =-AMOMZ
      AM(2,5)  =-RIXZ
!     AM(3,5)  =-AMOMX
      AM(5,5)  = RIXX
      AM(1,6)  = AMOMX
      AM(4,6)  =-AMOMZ
!     AM(6,6)  = POLI

      SHAPETAM = matmul( SHAPET  , AM    )
      AMLOC0   = matmul( SHAPETAM, SHAPE )
      AFLOC0   = matmul( AMLOC0  , UT2   )

!----------------------------------------------
!
!  calculation of the local DAMPING matrix [C]
!
!----------------------------------------------
!j    call IIxAxS (AM, ATDAx2, DENS, AMOMX, AMOMZ, RIXX, RIXZ, RIZZ, 1.d0)

!j    SHAPETAM = matmul( SHAPET  , AM    )
!j    ACLOC0   = matmul( SHAPETAM, SHAPE )
!j    AFLOC0   = AFLOC0 + matmul( ACLOC0  , UT1 )

!----------------------------------------------
!
!  calculation of the local STIFFNESS matrix [K]
!
!----------------------------------------------
!UTT--> UTT0
!A  --> ATDDA
!j    call IIxAxS (AM, ATDDA, DENS, AMOMX, AMOMZ, RIXX, RIXZ, RIZZ, 1.d0)

!j    SHAPETAM = matmul( SHAPET , AM    )
!j    AKLOC0   = AKLOC0 + matmul( SHAPETAM, SHAPE )
!j    AFLOC0   = AFLOC0 + matmul( AKLOC0  , UT )


!- F
! R --> ATDDR
! A --> ATDDA

!j    call IIxAtDDR_AtDDAr0 ( AQ, ATDDA, ATDDR, DENS, AMOMX, AMOMZ, RIXX, RIXZ, RIZZ, HTA1, 1.d0 )

!----
!j    AFLOC0  = AFLOC0 + matmul( SHAPET, AQ )

!---Replace (Add without Zeroing)

      allocwgauss = 0.5d0*(L2-L1)*WGAUSS(k)

!j    AKLOC = AKLOC + allocwgauss*AKLOC0 
!j    ACLOC = ACLOC + allocwgauss*ACLOC0 
      AMLOC = AMLOC + allocwgauss*AMLOC0 
      AFLOC = AFLOC - allocwgauss*AFLOC0 

!----------------------------------------------
! QS
!----------------------------------------------

!---- For all q's

!j    call MATRLOC_Q ( AMLOCNQ    , ACLOCNQ, AKLOCNQ, SHAPET, UTT0 , UTT1 , F_Grav, &
!j                     DENS       , AMOMX  , AMOMZ  , RIXX  , RIXZ , RIZZ , HTA1  , &
!j                     allocwgauss, XCM    , ZCM    , nn_el , nb_el, NDFPE, 1.d0     )


   enddo !k (IORDER)


 END Subroutine FLOODMATRLOC
!--------------------------------------------------------------------------------------------------------
!
!----- hydro-wave module
! 1. set upper-lower limits depend on the wave case if it is possible
! 2. fenton 1 wave
! 3. fenton irregural [Belibasakis idea]
!
!----- Morison
! 1. no elasticity is considered during the search procedure for the wet/dry length
!    this is the only part, everywhere else elasticity is introduced [ Uelst, L_Unit_t]
!    ***** Elasticity now is introduced ******
! SOS --->  integrate at the deformed length
!
!---- Buoyancy Side
! 1. no elasticity is considered during the search procedure for the wet/dry length
! 2. no elasticity: buoyancy is applied at the undeformed geometry [only consider the 6 float q's]
! 3. check cone --> show to Spyros [normal vectors, Pdyn sign, Restoring Moments] 
!    ***** Pdyn correct *****
!    ***** Spyros has added the relative system term *****
!
!---- Buoyancy Tap
! 1. no elasticity: buoyancy is applied at the undeformed geometry [only consider the 6 float q's]
! 2. Morison      : no elasticity [undeformed geometry, Uelst only due to 6 float q's]
! 3. THINKING     : maybe zero the opposite Morison force ????
!    ***** most probably NO *****
! 4. Morison froude-krylov force added, even worst...
! 5. stretching needs actual wet length, while buoyancy MSL, Also think how to treat Pdyn - Prel when stretch is applied
!--------------------------------------------------------------------------------------------------------
