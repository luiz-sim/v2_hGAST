!-------------------------------------------------------------------------------------------------------
 Subroutine foilfs_init (NTIME)
!-------------------------------------------------------------------------------------------------------

 Use Craft  !NBLADE_ae, strip
 Use Foilfs_mod

   implicit none

   integer :: NTIME, iblade, lstrip, i
!  real(8) :: psi

   real(8) :: chord,twist
   real(8) :: WEFF,PHI,ALPHA0
   real(8) :: AP, A
   real(8) :: rRD, rBld
   real(8) :: AA(3,3), AT(3,3)
   real(8) :: WEFF_RD(3), WEFFP_RD(3), WEFFP_Bld(3)


   if (NTIME>1.or.INDXFlap/=1) return


   do iblade       = 1, NBLADE_ae
   do lstrip       = 1, NSTRIP

!------ Wake initialization for flaps
      GWAKE(iblade,lstrip,1:NWAKE)= 0.d0

      AP           = 0.d0
      A            = 0.d0
      AA (1:3,1:3) = strip(iblade,lstrip)%A(1:3,1:3)
      AT (1:3,1:3) = transpose(AA(1:3,1:3))
      rRD          = strip(iblade,lstrip)%rRD
      rBld         = strip(iblade,lstrip)%rBld
      chord        = strip(iblade,lstrip)%chord
      twist        = strip(iblade,lstrip)%twist
      WEFF_RD(1:3) = strip(iblade,lstrip)%UWINDTB(1:3) - strip(iblade,lstrip)%UBODY(1:3)

!------ limit2 WEFF_RD(3)
      if (IPARKED==0.and.INOINDUCT==0)&
      WEFF_RD  (3  ) = max(WEFF_RD(3),0.1d0)
      WEFFP_RD (1  ) = WEFF_RD(1)*(1.d0+AP)
      WEFFP_RD (2  ) = WEFF_RD(2)
      if (A<=1.d0) &
      WEFFP_RD (3  ) = WEFF_RD(3)*(1.d0-A )
      if (A >1.d0) &
      WEFFP_RD (3  ) =-WEFF_RD(3)*(1.d0-A )
      WEFFP_Bld(1:3) = matmul(AT(1:3,1:3), WEFFP_RD(1:3))
!------ limit3 WEFFP_Bld(3)
!     if (IPARKED==0.and.INOINDUCT==0)&
!     WEFFP_Bld(3  ) = dmax1(WEFFP_Bld(3),1.0d0)
!     WEFFP_RD (1:3) = matmul(AA(1:3,1:3), WEFFP_Bld(1:3))

!------ find AoA and Weff wrt chord c.s.
      WEFF   = dsqrt (WEFFP_Bld(1)**2+WEFFP_Bld(3)**2)
      PHI    = datan2(WEFFP_Bld(3),   WEFFP_Bld(1)   )
      ALPHA0 = PHI-twist
      ALPHAFOIL(iblade, lstrip) = ALPHA0           

!     PHI    = datan2(WEFFZ,WEFFX)
!     WEFF   = dsqrt (WEFFX**2+WEFFZ**2)
!     ALPHA0 = PHI-THET-THY0
!     ALPHAFOIL(iblade, lstrip) = ALPHA0           

!-- Calculation of initial circulation for each flap
!-- in body iblade and strip lstrip
      if ( IFLAPLOC(iblade, lstrip)/=1) cycle

         NNODE       = NPOINTS(iblade,1)
      do i           = 1, NNODE
         XNODE (1,i) = XNSTR(iblade,lstrip,i)*chord   
         XNODE (2,i) = YNSTR(iblade,lstrip,i)*chord   
         XNODEL(1,i) = XNODE(1,i)               
         XNODEL(2,i) = XNODE(2,i)               
         YNODEOLDN (iblade, lstrip, i) = XNODE(2,i)*chord
      enddo !i

      call GEOM
      call ZERO ( WEFF, ALPHA0, CHORD, iblade, lstrip )

      do i = 1, NNODE-1
         circ_old(iblade, lstrip, i) = circ_int(iblade, lstrip, i)
      enddo
!qq
!      if ( IFLAPLOC(iblade, lstrip)==1) then
!           psi                    = -UT_el(NDFB_el+NQSW)-UThub_el !+PHI0(nbod)
!           FLCONTR(iblade,lstrip) = 10.d0*dsin(psi)
!            if (lstrip==35)&
!            write(330+iblade,'(100e15.5)') TIME,psi, FLCONTR(iblade,lstrip)
!      endif 
   enddo !lstrip
   enddo !iblade


 END Subroutine foilfs_init
!-------------------------------------------------------------------------------------------------------
 Subroutine UPDATE_FOILFS
!-------------------------------------------------------------------------------------------------------

 Use Craft  !NBLADE_ae, strip
 Use Foilfs_mod

   implicit none

   integer :: iblad, lstrip, i


   if (INDXFlap/=1) return

!--- Foilfs update parameters
   do    iblad  = 1, NBLADE_el
      do lstrip = 1, NSTRIP

          do i  = NWAKE, 2, -1
             GWAKE    (iblad,lstrip,i) = GWAKE(iblad, lstrip,i-1)
          enddo

          do i = 1, NPOINTS(iblad,1)
             YNODEOLDN(iblad,lstrip,i) = YTEMP(iblad,lstrip,i)
             circ_old (iblad,lstrip,i) = circ_int(iblad,lstrip,i)
          enddo
          
          do i = 1, NPOINTS(iblad,1)-1
             XC4N     (iblad,lstrip,i) = (XTEMP(iblad,lstrip,i+1)+XTEMP(iblad,lstrip,i))/2.d0
          enddo
      enddo !lstrip
   enddo !iblad


 END Subroutine UPDATE_FOILFS
!-------------------------------------------------------------------------------------------------------
!  Subroutine : AIRGEOMETRY (iblade, RLD, lstrip, RCB2, XAIRNODE, YAIRNODE, XNSTR, YNSTR, NSPANB2, NPOINTS)
!
!
!  AIRGEOMETRY - Determines the airfoil geometry  at
!                a given radius along the span by means of
!                linear approximation.
!           
!               RLD      = Strip radial position, not dimensionless(=RADSTRIP(lstrip)*RTIP
!               iblade   = body of interpolation
!               lstrip   = strip of interpolation
!               XNSTR    = Coordinate after interpolation
!               YNSTR    =  -//-              -//-
!               XAIRNODE = Cordination before interponation
!               YAIRNODE =  -//-              -//-
!
!  Interpolation is done at x, y nodes of the 1st airfoil given
!
!-------------------------------------------------------------------------------------------------------
 Subroutine AIRGEOMETRY (iblade, RS, lstr, RB, XNODEAG, YNODEAG, XNS, YNS, SPANB2, POINTS, NST)
!------------------------------------------------------------------------------------------------------

 Use Foilfs_mod

   implicit none
   integer :: i, ispan, k, iblade, lstr, POINTS, NST
   real(8) :: wgt1, wgt2, RS
   real(8) :: RB(MBLADE, NSPANM)
   real(8) :: XNODEAG(NBLADE_el, NST, NNODEM), YNODEAG(NBLADE_el, NST, NNODEM)
   real(8) :: XNS(NBLADE_el, NST, NNODEM), YNS(NBLADE_el, NST, NNODEM)
   integer :: SPANB2(MBLADE) 
 

   do i  = 1, SPANB2(iblade)  
      if (RS.le.RB(iblade, i)) then
         ispan = i-1 
         goto 11
      endif 

   enddo !i

   write(*,*) 'PROBLEM!!!!!!'

11 continue

   wgt1=(RS-RB(iblade, ispan))  /(RB(iblade, ispan+1)-RB(iblade, ispan)) ! weight 1
   wgt2=(RB(iblade, ispan+1)-RS)/(RB(iblade, ispan+1)-RB(iblade, ispan)) ! weight 2

   do k=1, POINTS
      YNS(iblade,lstr,k)=YNODEAG(iblade,ispan,k)*(1-wgt1)+YNODEAG(iblade,ispan+1,k)*(1-wgt2)
      XNS(iblade,lstr,k)=XNODEAG(iblade,ispan,k)
   enddo !k
  

 END Subroutine AIRGEOMETRY
!-------------------------------------------------------------------------------------------------------
 Subroutine GEOM
!------------------------------------------------------------------------------------------------------

 Use Foilfs_mod

   implicit none

   integer :: i
   real(8) :: XDERIV(NNODEM), YDERIV(NNODEM), DYX(NNODEM)


   do i         = 1, NNODE
      XDERIV(i) = XNODE(1,i)
      YDERIV(i) = XNODE(2,i)
   enddo

!w write(*,*)'GEOM-NNODE,NNODEM',NNODE,NNODEM
!--- calculates derivative DYX of the airfoil
   call DERIVE1 (NNODEM, NNODE, XDERIV, YDERIV, DYX)

   do  i        = 1, NNODE-1
       XCP  (i) = (XNODE (1,i+1) + XNODE (1,i)) / 2.d0
       YCP  (i) = (XNODE (2,i+1) + XNODE (2,i)) / 2.d0
       XCPL (i) = (XNODEL(1,i+1) + XNODEL(1,i)) / 2.d0
       YCPL (i) = (XNODEL(2,i+1) + XNODEL(2,i)) / 2.d0
       DX   (i) =  XNODE (1,i+1) - XNODE (1,i)
       DYXCP(i) = (DYX(i)+DYX(i+1))/2.d0
   enddo !i


 END Subroutine GEOM       
!-------------------------------------------------------------------------------------------------------
 Subroutine ZERO (Uinf,alpha,CHORD, nbod, istrip)
!------------------------------------------------------------------------------------------------------

 Use Foilfs_mod

   implicit none
   integer :: i, nbod, istrip
   real(8) :: Uinf, alpha, CHORD


      alpha    = 0.d0
   do i        = 1, NNODE-1              
      WWASH(i) = Uinf*dcos(alpha)*DYXCP(i) - Uinf*dsin(alpha)
   enddo

   call FOURIER (Uinf,CHORD)

!---- Calculation of AZER0
!   AZER0pot_ae(nbod,istrip) =  -((AFOUR(1)-alpha) + AFOUR(2)/2 )

   call CIRCULATION (0,Uinf,CHORD, nbod, istrip)

!---- calculate potential Cm/4 value
   call MOMPOT0     (  Uinf,alpha,CHORD, nbod, istrip)


 END Subroutine ZERO
!--------------------------------------------------------------------
!  Subroutine : FOURIER
!
!  calculates [A] elements of Fourier Series
!  careful A(1) corresponds to A0 element of series etc.
!  giving WWASH vertical speed
!--------------------------------------------------------------------
 Subroutine FOURIER (Uinf, CHORD)
!--------------------------------------------------------------------

 Use Foilfs_mod

   implicit none
   integer :: n,i
   real(8) :: X0, Y0, X0L, Y0L, theta, Dtheta, CHORD, Uinf


   do    n        = 1, NFOUR
         AFOUR(n) = 0.d0
      do i        = 1, NNODE-1
         X0       = XCP (i)
         Y0       = YCP (i)
         X0L      = XCPL(i)
         Y0L      = YCPL(i)
         theta    = dacos(1.d0-2.d0*X0L/CHORD) 
         Dtheta   = dacos(1.d0-2.d0*XNODEL(1,i+1)/CHORD) &
                   -dacos(1.d0-2.d0*XNODEL(1,i  )/CHORD)

         if (n.eq.1) then
            AFOUR(n) = AFOUR(n)&
                       -1./pi*WWASH(i)/Uinf*Dtheta
         else
            AFOUR(n) = AFOUR(n)&
                       +2./pi*WWASH(i)/Uinf*dcos((n-1)*theta)*Dtheta
         endif
      enddo ! control points
   enddo ! fourier Series elements


 END Subroutine FOURIER       
!------------------------------------------------------------------------
! Subroutine : CIRCULATION
!
! calculates total circilation at a specific time using fourrier elements.
!------------------------------------------------------------------------
 Subroutine CIRCULATION (NT, Uinf, CHORD, nb, istr)
!------------------------------------------------------------------------

 Use Foilfs_mod

   implicit none

   integer :: i, n, NT, nb, istr
   real(8) :: Uinf, CHORD, theta, Dtheta, circth


      circ_int(nb, istr,  1) = 0.d0
      circ    (nb, istr, NT) = 0.d0
   do i      = 1, NNODE-1
      theta  = dacos(1.d0-2.d0*XCPL  (  i  )/CHORD)
      Dtheta = dacos(1.d0-2.d0*XNODEL(1,i+1)/CHORD)&
              -dacos(1.d0-2.d0*XNODEL(1,i  )/CHORD)
      circth = 0.d0

      do n = 1, NFOUR
        if (n.eq.1) then
           circth = circth + 2.d0*Uinf*AFOUR(n)*1.d0/dtan(theta/2.d0)
        else
           circth = circth + 2.d0*Uinf*AFOUR(n)*dsin((n-1)*theta)
        endif
      enddo !n
           circ    (nb, istr, NT ) =  circ    (nb, istr, NT) + 0.5d0*CHORD*circth*dsin(theta)*Dtheta
           circ_int(nb, istr, i+1) =  circ_int(nb, istr, i ) + 0.5d0*CHORD*circth*dsin(theta)*Dtheta
   enddo !i


 END Subroutine CIRCULATION   
!------------------------------------------------------------------------
!  Subroutine : MOMPOT0
!------------------------------------------------------------------------
 Subroutine MOMPOT0 (Uinf, alpha, CHORD, nb, istr)
!------------------------------------------------------------------------

 Use Foilfs_mod

   implicit none
   integer :: i, n, nb, istr
   real(8) :: Uinf, alpha, CHORD, theta, Dtheta, circth
   real(8) :: circtheta, XC4


         cmpot(nb, istr) = 0.d0
   do    i       = 1, NNODE-1
         theta   = dacos(1.d0-2.d0*XCPL  (  i  )/CHORD)
         Dtheta  = dacos(1.d0-2.d0*XNODEL(1,i+1)/CHORD)&
                  -dacos(1.d0-2.d0*XNODEL(1,i  )/CHORD)
         circth  = 0.d0
      do n       = 1, Nfour
         if (n.eq.1) then
          circth = circth + 2.d0*Uinf* AFOUR(n)*1.d0/dtan(theta/2.d0)
         else
          circth = circth + 2.d0*Uinf* AFOUR(n)     *dsin((n-1)*theta)
         endif
      enddo
     
      circtheta  = 0.5d0*CHORD*circth*dsin(theta)*Dtheta
      XC4        = 0.5d0*(XNODE(1,i+1)+ XNODE(1,i)) - 0.25d0*CHORD

      cmpot(nb, istr) = cmpot(nb, istr) - AIRDEN*Uinf*dcos(alpha)*circtheta*XC4
   enddo
      cmpot(nb, istr) = cmpot(nb, istr)/(0.5d0*AIRDEN*Uinf**2*CHORD**2)


 END Subroutine MOMPOT0        
!---------------------------------------------------------------------
!
!  Subroutine : FOILFS
!
!  Implementation of FOILFS code
!------------------------------------------------------------------------------
 Subroutine FOILFS (iblade, lstrip, NTIME, ANGMOVE, AFLAP, ALPHA0, WEFF, CHORD)
!------------------------------------------------------------------------------

 Use Foilfs_mod

   implicit none

   integer, intent(in) :: iblade, lstrip, NTIME
   real(8), intent(in) :: ANGMOVE, AFLAP, ALPHA0, WEFF, CHORD

   real(8) :: RLX, diff, gwake_new
   integer :: ITER, ic


   do ic           = 1, NNODE
      XNODE (1,ic) = XNSTR(iblade,lstrip,ic)*CHORD
      XNODE (2,ic) = YNSTR(iblade,lstrip,ic)*CHORD
      XNODEL(1,ic) = XNODE(1,ic)               
      XNODEL(2,ic) = XNODE(2,ic)               
      YNODEOLD(ic) = YNODEOLDN (iblade, lstrip, ic)
   enddo

   ALPHAFOIL(Iblade,lstrip) = ALPHA0 !foilfs


   call MOVEBO(ANGMOVE, CHORD, AFLAP, lstrip)

   do ic        = 1, NNODE-1
      WWASH(ic) = WEFF*DYXCP(ic)
   enddo !ic

   call FOURIER (WEFF,CHORD)
   
   AZER_ae0(iblade,lstrip)  = -((AFOUR(1)-ANGMOVE) + AFOUR(2)/2 ) 
!  cmpot(iblade, lstrip) = pi/4*(AFOUR(3)-AFOUR(2))

   call WAKEGEOM    (iblade, lstrip, WEFF)
   call DOWNWASH    (NTIME , iblade, lstrip, CHORD, WEFF)

!--- initialize near wake for iterations
   GWAKE(iblade, lstrip, 1) = GWAKE(iblade, lstrip, 2)

  do  ITER   = 1, 30
   call DOWNWASH_1  (iblade, lstrip, CHORD                )
   call FOURIER     (        WEFF  , CHORD                )
   call CIRCULATION (NTIME , WEFF  , CHORD, iblade, lstrip)

!----- convergence for gwake
   RLX       = 0.75d0
   gwake_new = -( circ(iblade, lstrip, NTIME) - circ(iblade, lstrip, NTIME-1) )&
                / DXW(iblade, lstrip, 1)
   diff      = dabs( gwake_new - GWAKE(iblade, lstrip, 1) )
   GWAKE(iblade, lstrip, 1) = (1-RLX)*gwake_new + RLX * GWAKE(iblade, lstrip, 1) 

   if ((diff).lt.(1.0d-1)) goto 11

!CFan_START
!  if ((iblade==1).and.(lstrip==40)) then
!     write(40,105) TIME,iblade,lstrip, GWAKE(iblade, lstrip, 1),diff,J 
!105   format   (f20.10,2X,i2,2X,i2,2X,2f20.10,2X,i2)
!  endif
!CFan_END
!qq check
!! if (iblade==1.and.lstrip==29)&
!! write(40,'(7e15.5,i5)') TIME, ANGMOVE, AFLAP, CHORD, WEFF, GWAKE(iblade, lstrip, 1),diff,J
  enddo !ITER

   write(20,*) 'FOILFS not converged, blade=',iblade,' strip=',lstrip

11 call FORCE (CHORD, WEFF, iblade, lstrip, NTIME)

   do ic                       = 1, NNODE
     YTEMP(iblade, lstrip, ic) = XNODE(2,ic)
     XTEMP(iblade, lstrip, ic) = XNODE(1,ic)
   enddo


 END Subroutine FOILFS
!---------------------------------------------------------------------
!  Subroutine : DERIVE1(NM,NX,X,Z,DZX)
!
!
!  DERIVE1 calculates the discrete sequence of derivatives DZX[1:NX]
!          of a function z(x) corresponding to the discrete data
!          X[1:NX], Z[1:NX]
!
!          NM = maximum values of the NX
!---------------------------------------------------------------------
 Subroutine DERIVE1 ( NM , NX , X , Z , DZX )
!---------------------------------------------------------------------
   implicit none

   integer :: NX, NM, I1, I2, I3, I
   real(8) :: DX1,DY1,DX2,DY2
   real(8) :: X(NM),Z(NM),DZX(NM)

!--- i = 1
      I1       = 1
      I2       = 2
      I3       = 3
      DX1      = X(I2) - X(I1)
      DY1      = Z(I2) - Z(I1)
      DX2      = X(I3) - X(I2)
      DY2      = Z(I3) - Z(I2)
      DZX (I1) = DY1/DX1+(DY1*DX2-DX1*DY2)/DX2/(DX2+DX1)


!--- i = 2 : NX-1
   do I        = 2, NX-1
      I1       = I - 1
      I2       = I
      I3       = I + 1
      DX1      = X(I2) - X(I1)
      DY1      = Z(I2) - Z(I1)
      DX2      = X(I3) - X(I2)
      DY2      = Z(I3) - Z(I2)
     if (DX1.eq.-DX2) then
      DZX (I2) = 1.d0
      cycle
     endif
      DZX (I2) = (DY2*DX1/DX2+DY1*DX2/DX1)/(DX1+DX2)
   enddo


!--- i = NX
      I1       = NX - 2
      I2       = NX - 1
      I3       = NX
      DX1      = X(I2) - X(I1)
      DY1      = Z(I2) - Z(I1)
      DX2      = X(I3) - X(I2)
      DY2      = Z(I3) - Z(I2)
      DZX (I3) = DY2/DX2+(DY2*DX1-DX2*DY1)/DX1/(DX2+DX1)


 END Subroutine DERIVE1        
!-----------------------------------------------------------------------
!  Subroutine : MOVEBO 
!
!  Sets Airfoil in current step potition and deployes flap 
!-----------------------------------------------------------------------
 Subroutine MOVEBO(ANGMOV, CHORD, BFLAP, lstrip)
!-----------------------------------------------------------------------
 
 Use Foilfs_mod

   implicit none
 
   integer :: i, lstrip
   real(8) :: XDERIV(NNODEM), YDERIV(NNODEM), DYX(NNODEM)
   real(8) :: ANGMOV, CHORD !, DANGMOV, DDANGMOV
   real(8) :: XEXT, RCC, YB, YFLAP, BFLAP, DT !, FLAP
   real(8) :: ksi, hmax


   DT   = DT_el
   XEXT = 1.d0-XFLAP(lstrip)

!--- Flap10%
   if     (FLAPMODE (1) ==1) then;  RCC  = ((XEXT*dsin(pi/180.d0))**2 + XEXT**2)/(2.d0*XEXT*dsin(pi/180.d0))
!--- Flap 30%
   elseif (FLAPMODE (1) ==2) then;  RCC = 3.d0
   endif

   do i          = 1, NNODE
      XNODE(1,i) = XNODEL(1,i)
      XNODE(2,i) = XNODEL(2,i)
      YB         = 0.d0

      if ( (XNODEL(1,i)/CHORD) .lt. XFLAP(lstrip)) goto 1

!--- FLap 10%
   if     (FLAPMODE (1) ==1) then; YB  = dsqrt(RCC**2-(XNODE(1,i)/CHORD-(1.-XEXT))**2)-RCC
!--- Flap 30%
   elseif (FLAPMODE (1) ==2) then; YB = (0.0950*RCC - 0.2087)+(0.778 - 0.3665*RCC)*(XNODE(1,i)/CHORD)     +  &
                                        (0.4654*RCC - 0.9445)*((XNODE(1,i)/CHORD)**2)+(0.3703 - 0.19391*RCC) &
                                       *((XNODE(1,i)/CHORD)**3)
!!--- Flap 10% Avatar
    elseif (FLAPMODE (1) ==3) then
        ksi=((XNODEL(1,i)/CHORD)- XFLAP(lstrip))/(1.d0-XFLAP(lstrip))
        hmax=(1.d0-XFLAP(lstrip))*dtan(1.d0*pi/180.d0)
        YB =-hmax*(ksi**2*(3.d0-ksi)/2.d0)
   endif

      YB    = YB * CHORD

      YFLAP = YB*BFLAP  
      XNODE(1,i) = XNODE(1,i)
      XNODE(2,i) = XNODE(2,i) + YFLAP

 1    continue

      XNODE(2,i) = XNODE(2,i) - ANGMOV*(XNODEL(1,i)-XCLPI*CHORD)
   enddo


   do i         = 1, NNODE
      XDERIV(i) = XNODE(1,i)
      YDERIV(i) = XNODE(2,i)
   enddo

!--- calculates derivative DYX of the airfoil
   call DERIVE1 (NNODEM, NNODE, XDERIV, YDERIV, DYX)

!--- spatial and time derivative of geometry
   do i        = 1, NNODE-1
      XCP  (i) = (XNODE (1,i+1) + XNODE (1,i)) / 2.d0
      YCP  (i) = (XNODE (2,i+1) + XNODE (2,i)) / 2.d0
      XCPL (i) = (XNODEL(1,i+1) + XNODEL(1,i)) / 2.d0
      YCPL (i) = (XNODEL(2,i+1) + XNODEL(2,i)) / 2.d0
      DX   (i) =  XNODE (1,i+1) - XNODE (1,i)
      DYXCP(i) = (DYX(i)+DYX(i+1))/2.d0
   enddo

   do i        = 1, NNODE
      DYT  (i) = (XNODE(2,i)-YNODEOLD(i))/DT
   enddo
 

 END Subroutine MOVEBO
!-----------------------------------------------------------------------
!  Subroutine : WAKEGEOM
!  
!  Generates wake geometry
!-----------------------------------------------------------------------
 Subroutine WAKEGEOM (nb, istr, Uinf)
!-----------------------------------------------------------------------

 Use Foilfs_mod

   implicit none
    
   integer  :: nb, istr, i
   real(8)  :: alpha, DT, Uinf


      alpha             = ALPHAFOIL (nb, istr)      
      DT                = DT_el
      DXW  (nb, istr,1) = Uinf*DT
      XWAKE(nb, istr,1) = XNODE(1,NNODE)
      YWAKE(nb, istr,1) = XNODE(2,NNODE)
   do i                 = 2, NWAKE
      DXW  (nb, istr,i) = DXW  (nb,istr,i-1)
      XWAKE(nb, istr,i) = XWAKE(nb,istr,i-1) + DXW(nb,istr,i)*dcos(alpha)
      YWAKE(nb, istr,i) = YWAKE(nb,istr,i-1) + DXW(nb,istr,i)*dsin(alpha)
   enddo


 END Subroutine WAKEGEOM
!-------------------------------------------------------------------- 
!  Subroutine : DOWNWASH
!
!  calculates downwash
!-------------------------------------------------------------------- 
    Subroutine DOWNWASH (NT, nb, istr, CHORD, WEFF)
!-------------------------------------------------------------------- 

 Use Foilfs_mod

   implicit none
   integer  :: NT, nb, istr, i
   real(8)  :: X0, Y0, X0L, Y0L, theta, Dtheta
   real(8)  :: CHORD, alpha, WEFF, WX, WY, WW, UW


      alpha  = ALPHAFOIL (nb, istr)
      WX     = WEFF * dcos (alpha) 
      WY     = WEFF * alpha         !WEFF * dsin (alpha)
   do i      = 1, NPOINTS(nb,1)-1
      X0     = XCP (i)
      Y0     = YCP (i)
      X0L    = XCPL(i)
      Y0L    = YCPL(i)
      theta  = dacos(1.d0-2.d0*X0L          /CHORD)
      Dtheta = dacos(1.d0-2.d0*XNODEL(1,i+1)/CHORD)&
              -dacos(1.d0-2.d0*XNODEL(1,i  )/CHORD)

      WWASH_1(i) =  WX*DYXCP(i) - WY + (DYT(i)+DYT(i+1))/2.d0
      WW     = 0.d0

      if (NT.ge.2) &
      call WAKE_CONTRIB(X0,Y0, UW, WW, nb, istr, alpha)

      WWASH_1(i) = WWASH_1(i) - WW
   enddo ! control points (i)


 END Subroutine DOWNWASH
!-------------------------------------------------------------------------
!  Subroutine : WAKE_CONTRIB
!
!  Far wake contribution         
!-------------------------------------------------------------------------
 Subroutine WAKE_CONTRIB (X0, Y0, UW, WW, iblade, lstr, alpha)
!-------------------------------------------------------------------------

 Use Foilfs_mod

   implicit none

   integer  :: iblade, lstr, j
   real(8)  :: X0, Y0, UW, WW, XW1, YW1, DXW1, alpha
   real(8)  :: CS, SN,U1S, V1S, U1D, V1D    


      UW   = 0.d0
      WW   = 0.d0
      CS   = dcos(alpha)
      SN   = dsin(alpha)
   do j    = 2, NWAKE
      XW1  = XWAKE  (iblade, lstr, j)
      YW1  = YWAKE  (iblade, lstr, j)
      DXW1 = DXW    (iblade, lstr, j)
      
      call SOURV (X0, Y0, XW1, YW1, DXW1, CS, SN, U1S, V1S, U1D, V1D)

      UW   = UW+U1D*GWAKE(iblade, lstr, j)
      WW   = WW+V1D*GWAKE(iblade, lstr, j)
   enddo


 END Subroutine WAKE_CONTRIB       
!----------------------------------------------------------------------
!  Subroutine : SOURV
!
!  Velocity calculation from a piecewise constant
!  SOURCE and VORTICITY distributions
!----------------------------------------------------------------------
 Subroutine SOURV (X0, Y0, X1, Y1, DS, CB, SB, U1S, V1S, U1D, V1D) 
!----------------------------------------------------------------------

 Use Foilfs_mod
 
   implicit none

   real(8) :: X0L, Y0L, X0, Y0, X1, Y1, DS, CB, SB, U1S, V1S, U1D, V1D
   real(8) :: TA1, TA2, V1L, U1L


      X0L = (X0-X1)*CB+(Y0-Y1)*SB
      Y0L =-(X0-X1)*SB+(Y0-Y1)*CB
      TA1 = X0L
      TA2 = X0L-DS

   if (DS.lt.1.d-15) then
      U1S = 0.d0
      V1S = 0.d0
      U1D = 0.d0
      V1D = 0.d0
      return
   endif

   if     (DABS(Y0L).gt.1.d-15) then

      V1L =-1.d0/(2.d0*PI)*(DATAN(TA2/Y0L)-DATAN(TA1/Y0L))
      U1L =-1.d0/(2.d0*PI)*DLOG(DSQRT((TA2**2+Y0L**2)/(TA1**2+Y0L**2)))

   else if (  (Y0L.lt.1.d-15).and.  &
              (Y0L.gt.0.d0  ).and.  & 
              (X0L.gt.0.d0  ).and.  &
              (X0L.lt.   DS)        ) then

      V1L = 0.5d0
      U1L =-1.d0/(2.d0*PI)*DLOG(DSQRT(TA2**2/TA1**2))

   else if (  (Y0L.gt.-1.d-15).and. &
              (Y0L.lt. 0.d0  ).and. &
              (X0L.gt. 0.d0  ).and. &
              (X0L.lt.    DS)       ) then
      V1L =-0.5d0
      U1L =-1.d0/(2.d0*PI)*DLOG(DSQRT(TA2**2/TA1**2))

   else if (  (DABS(Y0L).lt.1.0d-15).and. &
              ((X0L.gt.DS).or.(X0L.lt.0.d0))   ) then
      V1L = 0.d0
      U1L =-1.d0/(2.d0*PI)*DLOG(DSQRT(TA2**2/TA1**2))

   end if

      U1S = U1L*CB-V1L*SB
      V1S = U1L*SB+V1L*CB
      U1D = V1S
      V1D =-U1S


 END Subroutine SOURV
!--------------------------------------------------------------------
!  Subroutine : DOWNWASH_1
!
!  calculates downwash adding the XW(1) circulation
!--------------------------------------------------------------------
 Subroutine DOWNWASH_1 (nb, istr, CHORD)
!--------------------------------------------------------------------

 Use Foilfs_mod

   implicit none

   integer  :: nb, istr, i
   real(8)  :: X0, Y0, X0L, Y0L, theta, Dtheta, CHORD
   real(8)  :: UW, WW, alpha


      alpha  = ALPHAFOIL( nb, istr)
   do i      = 1, NPOINTS(nb,1)-1
      X0     = XCP (i)
      Y0     = YCP (i)
      X0L    = XCPL(i)
      Y0L    = YCPL(i)
      theta  = dacos(1.d0-2.d0 * X0L          /CHORD) 
      Dtheta = dacos(1.d0-2.d0 * XNODEL(1,i+1)/CHORD)&
              -dacos(1.d0-2.d0 * XNODEL(1,i  )/CHORD)

      call WAKE_CONTRIB_1 (X0,Y0, UW, WW, nb, istr, alpha)

      WWASH(i) = WWASH_1(i) - WW
   enddo ! control points


 END Subroutine DOWNWASH_1
!-------------------------------------------------------------------------
!  Subroutine : WAKE_CONTRIB_1
!
!  Near wake contribution
!-------------------------------------------------------------------------
 Subroutine WAKE_CONTRIB_1 (X0, Y0, UW, WW, iblade, lst, alpha)
!-------------------------------------------------------------------------

 Use Foilfs_mod

   implicit none
   
   integer  :: iblade, lst
   real(8)  :: X0, Y0, UW, WW, XW1, YW1, DXW1, alpha, CS, SN
   real(8)  :: U1S, V1S, U1D, V1D


   UW   = 0.d0
   WW   = 0.d0
   XW1  = XWAKE (iblade, lst, 1)
   YW1  = YWAKE (iblade, lst, 1)
   DXW1 = DXW   (iblade, lst, 1)
   CS   = dcos(alpha)
   SN   = dsin(alpha)

   call SOURV (X0,Y0, XW1,YW1,DXW1, CS, SN, U1S ,V1S , U1D, V1D)

   UW   = UW+U1D*GWAKE(iblade, lst, 1)
   WW   = WW+V1D*GWAKE(iblade, lst, 1)


 END Subroutine WAKE_CONTRIB_1       
!--------------------------------------------------------------------
 Subroutine FORCE (CHORD, WEFF, nb, istr, NTI)
!--------------------------------------------------------------------

 Use Foilfs_mod

   implicit None

   integer :: nb, istr, NTI, i, n
   real(8) :: Dxx, XC4, CHORD, WEFF, WX, alpha
   real(8) :: DT, circtheta,circth,theta, Dtheta
   real(8) :: X0, Y0, UW, WW, w_wake1


   DT    = DT_el
   alpha = ALPHAFOIL(nb, istr)
   WX    = WEFF * dcos (alpha)
   
   clf (nb, istr, NTI) = 0.d0
   cmf (nb, istr, NTI) = 0.d0

 
   do i      = 1, NNODE-1
      theta  = dacos(1.d0-2.d0*XCPL  (  i  )/CHORD)
      Dtheta = dacos(1.d0-2.d0*XNODEL(1,i+1)/CHORD)&
              -dacos(1.d0-2.d0*XNODEL(1,i  )/CHORD)
      circth = 0.d0
      do n   = 1, NFOUR
         if (n.eq.1) then
            circth = circth + 2.d0*WEFF* AFOUR(n)*1.d0/dtan(theta/2.d0)
         else
            circth = circth + 2.d0*WEFF* AFOUR(n)     *dsin((n-1)*theta)
         endif
      enddo

      circtheta = 0.5d0*CHORD*circth*dsin(theta)*Dtheta

      if (NTI.lt.4) cycle

      Dxx  = XNODE(1,i+1)- XNODE(1,i)
      XC4  = XC4N(nb, istr, i) - 0.25d0 * CHORD

      if (i.eq.1) then
         clf (nb, istr, NTI) = clf (nb, istr, NTI) + ( circ_int(nb, istr, i+1) -&
                               circ_old(nb, istr, i+1)   )/DT*Dxx*AIRDEN
         cmf (nb, istr, NTI) = cmf (nb, istr, NTI) - ( circ_int(nb, istr, i+1) -&
                               circ_old(nb, istr, i+1)   )/DT*Dxx*AIRDEN*XC4 &
                               - AIRDEN*WEFF*circtheta*XC4 !- AIRDEN*WX*circtheta*XC4
      else

         clf (nb, istr, NTI) = clf   (nb, istr, NTI)+(0.5d0* (circ_int(nb, istr, i+1)+circ_int(nb, istr, i)) -&
                               0.5d0*(circ_old(nb, istr, i+1)+circ_old(nb, istr, i)))/DT &
                               *Dxx*AIRDEN
         cmf (nb, istr, NTI) = cmf   (nb, istr, NTI)-(0.5d0* (circ_int(nb, istr, i+1)+circ_int(nb, istr, i)) -&
                               0.5d0*(circ_old(nb, istr, i+1)+circ_old(nb, istr, i)))/DT &
                               *Dxx*AIRDEN*XC4- AIRDEN*WEFF*circtheta*XC4  !*Dxx*AIRDEN*XC4- AIRDEN*WX*circtheta*XC4 
     endif
   enddo !i

   clf (nb, istr, NTI) = clf (nb, istr, NTI)+ AIRDEN*WEFF*circ(nb, istr, NTI)  !AIRDEN*WX*circ(nb, istr, NTI)
   clf (nb, istr, NTI) = clf (nb, istr, NTI)/(0.5d0*AIRDEN*WEFF**2*CHORD)
   cmf (nb, istr, NTI) = cmf (nb, istr, NTI)/(0.5d0*AIRDEN*WEFF**2*CHORD**2)

!--- calculation of downwash angle -----------------------------------
   X0      = 0.25d0
   Y0      = 0.00d0
   w_wake1 = 0.00d0

   call WAKE_CONTRIB   (X0, Y0, UW, WW, nb, istr, alpha)
   w_wake1 = w_wake1 - WW

   call WAKE_CONTRIB_1 (X0, Y0, UW, WW, nb, istr, alpha)
   w_wake1 = w_wake1 - WW

!--- calculation of drag coefficient ---------------------------------
   cdf (nb, istr, NTI)  = clf (nb, istr, NTI)*dsin( datan(w_wake1/WX ))

!qq check
!! if (nb==1.and.istr==29) &
!! write(51,'(20e15.5)')   & 
!!   TIME               ,  & !1
!!   WEFF               ,  & !2
!!   theta              ,  & !3
!!  Dtheta              ,  & !4
!!   alpha              ,  & !5
!!   clf (nb, istr, NTI),  & !6
!!   cdf (nb, istr, NTI),  & !7
!!   cmf (nb, istr, NTI),  & !8
!!   circ(nb, istr, NTI),  & !9
!!   circtheta          ,  & !10
!!   w_wake1            ,  & !11
!!   GWAKE(nb, istr,  1)     !12


 END Subroutine FORCE
