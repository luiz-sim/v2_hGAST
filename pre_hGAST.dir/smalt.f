      PROGRAM BM1_FE_DMP
c
c     GAMESA Coordinate System
c     Linear Geometry Version
c     Properties are calculated along the mean-thickness line
c     New shear center calculation
c
      IMPLICIT NONE
      include '/home/seraf/hGAST/pre_hGAST.dir/paramet.for'
 
      REAL*8 wk
c
      INTEGER nl
      INTEGER im(nlayer, mlamin)
      INTEGER nlay(mlamin), lamsg(mseg, msct)
      INTEGER npts
      INTEGER nsct,  nsegm(msct), nskseg(msct), isgpt(2, mseg, msct)
      INTEGER sknseg(mseg, msct), webseg(mwebs, msct), nwebs(msct)
      INTEGER lpseg (mseg, mloops, msct)
      INTEGER nloops(msct),nlpsg(mloops, msct)
      INTEGER icase
c
      REAL*8 hl (nlayer,mlamin), angl(nlayer, mlamin)
      REAL*8 pm      (mater, mater)
      REAL*8 strength(mater, mater)
      REAL*8 cl(5, 5, mater),cc(3,3,nlayer) 
      REAL*8 h(nlayer),zbar(nlayer),rho(nlayer),thk(mlamin)
      REAL*8 r1inv, r2inv
      REAL*8 W1(nlayer), W3(nlayer), W4(nlayer), W5(nlayer)
      REAL*8 A0(3, 3), D0(3, 3), B0(3, 3), thks(mseg), A(3, 3, mseg), 
     &       D(3, 3, mseg), B(3, 3, mseg), as(2, 2, mseg)
      REAL*8 A0T(3, 3), B0T(3, 3), D0T(3, 3)
      REAL*8 AM0(3, 3), BM0(3, 3), DM0(3, 3), am(mseg), bm(mseg), 
     &       dm(mseg)
      REAL*8 xpnt(maxpoint, msct), ypnt(maxpoint, msct)
      REAL*8 xpno(maxpoint, msct), ypno(maxpoint, msct)
      REAL*8 aahstr(mloops), area0(mloops), aa(mseg), area(mseg)
      REAL*8 rshr(2)
 
      REAL*8 Loads(6, msct), EpsCur(6, msct)
      REAL*8 ezz(maxpoint, nlayer), ezs(maxpoint, nlayer)
      REAL*8 ezj(maxpoint, nlayer)
      REAL*8 szz(maxpoint, nlayer), szs(maxpoint, nlayer)
      REAL*8 szj(maxpoint, nlayer)
      
      REAL*8 qshear(mseg,maxpoint)
c      
      REAL*8 corweb
      REAL*8 xmem,ymem
c
c     Local variables
      INTEGER I, ICOSY, ILAM, IP, ISCT, ISEG, J
      INTEGER N, NLAM, NMAT
      REAL*8     XCEN, YCEN
      REAL*8   FCRIT1,FCRIT2
      CHARACTER*300 title
      CHARACTER*300 inpstr, inpgeo, inpload
      CHARACTER*300 n_p
      CHARACTER*350 p1, p2, p3, p4, p5
c
c     ... input files
      read(*,"(a)") n_p
      inpgeo=trim(adjustl(n_p))//"/other.dir/geometry.dat"   ! READ (*, '(a)') inpgeo                 !=========================================> seraf changes
      inpstr=trim(adjustl(n_p))//"/other.dir/structury.dat"  ! READ (*, '(a)') inpstr                 !
      OPEN (1, FILE=inpgeo, STATUS='old')
      OPEN (2, FILE=inpstr, STATUS='old')
c
      p1=trim(adjustl(n_p))//"/other.dir/section.out"
      OPEN (3, FILE=p1, STATUS='unknown')
c
c     ----------------------------------------------------------
c     ....read input data
c     ----------------------------------------------------------
c
      READ (2, 99001) title
99001 FORMAT (a80)
c
      READ (2, *) nmat, nlam
      CALL laminp(nmat, nlam, nlay, pm, cl, thk, im, hl, 
     &            angl, strength)
c     
      CALL sctinp(nsct , nsegm , nskseg, lamsg, isgpt, sknseg, 
     &            nwebs, webseg, nloops, lpseg, nlpsg)
c
      CALL input (nsct , npts  ,xpnt   ,ypnt)
 
c     ------------------------------------------------------------
c
       corweb = 0.0d0             ! No correction
c      corweb = 0.5d0             ! Half thickness correction
c      corweb = 1.0d0             ! Full thickness correction
c      
      DO isct = 1, nsct
c
         DO   n = 1,npts
         xpno(n,isct)=xpnt(n,isct)
         ypno(n,isct)=ypnt(n,isct)
         END DO 
         
         DO iseg = 1,nsegm(isct)
            ilam =   lamsg(iseg, isct)
            thks(iseg) = thk(ilam)
            NL = nlay(ilam)
               DO  ip  = 1, NL
               h  (ip) = hl(      ip,ilam)
               rho(ip) = pm(1, im(ip,ilam))
               CALL offply (im(ip,ilam), ip, angl(ip,ilam), cl, cc)
c               CALL offply1(im(ip,ilam), ip, angl(ip,ilam), cl, cc)
               END DO       
            CALL lamstf(nl, thk(ilam), h, zbar, cc, w4, w5, 
     &           a(1:3,1:3,iseg), b(1:3,1:3,iseg), d(1:3,1:3,iseg))
            CALL lamdns(nl, thk(ilam), h, zbar, rho, w1, w3, am(iseg), 
     &           bm(iseg), dm(iseg))
         END DO
         write(*,*) 'THK(3) = ',thks(3)
         write(*,*) 'THK(7) = ',thks(7)
         corweb   =  corweb*(thks(3)+thks(7))
         write(*,*) 'CorWeb = ',corweb
c
         rshr(1) = 0.
         rshr(2) = 0.
         CALL secstf(isct, thks, h, zbar, a, b, d, as, A0, B0, D0, 
     &               xpnt, ypnt,  
     &               nskseg(isct), lamsg(1:mseg,isct), 
     &               isgpt(1:2,1:mseg,isct),sknseg, nwebs(isct), 
     &               webseg, nloops(isct), 
     &               lpseg(1:mseg,1:mloops,isct), nlpsg(1:mloops,isct), 
     &               aahstr, area0, aa, area, corweb, n_p)
     
         xcen = -B0(1, 1)/A0(1, 1)
         ycen =  B0(1, 2)/A0(1, 1)
         WRITE (*, *) ' Centroid Coordinates', xcen, ycen
c        
c        TRANSLATE C.S. TO CENTROID
c       
         do   n = 1,npts
         xpnt(n,isct)=xpno(n,isct) - xcen
         ypnt(n,isct)=ypno(n,isct) - ycen
         enddo 
         CALL secstf(isct, thks, h, zbar, a, b, d, as, A0, B0, D0, 
     &               xpnt, ypnt,  
     &               nskseg(isct), lamsg(1:mseg,isct), 
     &               isgpt(1:2,1:mseg,isct),sknseg, nwebs(isct), 
     &               webseg, nloops(isct), 
     &               lpseg(1:mseg,1:mloops,isct), nlpsg(1:mloops,isct), 
     &               aahstr, area0, aa, area, corweb, n_p)
     
         CALL SHEARFLO(isct, a , b, A0, B0, D0, xpnt, ypnt,  
     &                 nskseg(isct), lamsg(1:mseg,isct), 
     &                 isgpt(1:2,1:mseg,isct), sknseg,  
     &                 nwebs(isct),webseg, nloops(isct), 
     &               lpseg(1:mseg,1:mloops,isct), nlpsg(1:mloops,isct), 
     &                 corweb,0.d0,1.d0,rshr(1),qshear, n_p)
c   
         CALL SHEARFLO(isct, a , b, A0, B0, D0, xpnt, ypnt,  
     &                 nskseg(isct), lamsg(1:mseg,isct), 
     &                 isgpt(1:2,1:mseg,isct), sknseg,  
     &                 nwebs(isct),webseg, nloops(isct), 
     &               lpseg(1:mseg,1:mloops,isct), nlpsg(1:mloops,isct), 
     &                 corweb,-1.d0,0.d0,rshr(2),qshear,n_p)
 
         WRITE (*, *) ' Shear Center Coordinates w.r.t Centroid',
     &                  rshr(1), rshr(2)
c
         WRITE (*, *) ' Enter Reference Coordinate System '
         WRITE (*, *) ' 0 for Section C.S.'
         WRITE (*, *) ' 1 for Shear Center C.S.'
         WRITE (*, *) ' 2 for Centroid C.S.'
         icosy=0          ! READ (*, *) icosy   !=====================================================> seraf changes
c
         IF ( icosy==2 ) THEN
         END IF
         IF ( icosy==1 ) THEN
            DO i = 1, npts
               xpnt(i, isct) = xpnt(i, isct) - rshr(1)
               ypnt(i, isct) = ypnt(i, isct) - rshr(2)
            END DO
            xmem = rshr(1) + xcen
            ymem = rshr(2) + ycen 
            rshr(1) = 0.
            rshr(2) = 0.         
         END IF
         IF ( icosy==0 ) THEN
            DO i = 1, npts
               xpnt(i, isct) = xpnt(i, isct) + xcen
               ypnt(i, isct) = ypnt(i, isct) + ycen
            END DO
            rshr(1) = rshr(1) + xcen
            rshr(2) = rshr(2) + ycen
            xmem  = - rshr(1)
            ymem  = - rshr(2)
         END IF
        
         CALL secstf(isct, thks,h,zbar, a, b, d, as, A0, B0, D0, 
     &               xpnt, ypnt,  
     &               nskseg(isct), lamsg(1:mseg,isct), 
     &               isgpt(1:2,1:mseg,isct),sknseg, nwebs(isct), 
     &               webseg, nloops(isct), 
     &               lpseg(1:mseg,1:mloops,isct), nlpsg(1:mloops,isct), 
     &               aahstr, area0, aa, area, corweb,n_p)
         
         CALL secdns(isct,thks,h,zbar, am, bm, dm, AM0, BM0, DM0,
     &               xpnt, ypnt,  
     &               nskseg(isct), lamsg(1:mseg,isct), 
     &               isgpt(1:2,1:mseg,isct), sknseg, nwebs(isct),
     &               webseg, nloops(isct), 
     &               lpseg(1:mseg,1:mloops,isct),nlpsg(1:mloops,isct), 
     &               aahstr, area0, aa, corweb)
c
         xcen = -B0(1, 1)/A0(1, 1)
         ycen =  B0(1, 2)/A0(1, 1)
         WRITE (*, *) ' Centroid Coordinates', xcen, ycen
         WRITE (*, *) ' Shear Center Coordinates', rshr(1), rshr(2)
         B0(1,3) =  rshr(1)*A0(3,1)-rshr(2)*A0(2,1) + B0(1,3)
         B0(2,3) =  rshr(1)*A0(3,2)-rshr(2)*A0(2,2)
         B0(3,3) =  rshr(1)*A0(3,3)-rshr(2)*A0(2,3)
         D0(1,3) =  rshr(1)*B0(3,1)-rshr(2)*B0(2,1) + D0(1,3)
         D0(2,3) =  rshr(1)*B0(3,2)-rshr(2)*B0(2,2) + D0(2,3)
         D0(3,3) =  rshr(1)*B0(3,3)-rshr(2)*B0(2,3) + D0(3,3)
         D0(3,1) =  D0(1,3)
         D0(3,2) =  D0(2,3)
         WRITE (*, *) ' B0(2,3),B0(3,3)'
         WRITE (*, *)   B0(2,3),B0(3,3) 
         
         call TRLFSM(A0, B0, D0, A0T, B0T, D0T, xmem, ymem)
c
c        ... output beam section properties here
c
         WRITE (3, *) '**SECTION =', isct, '   **IELM =', n
         WRITE (3, *) 'STIFFNESS MATRICES'
         WRITE (3, *) '(I,J)      A0  ', '      B0  ', '        D0  '
         WRITE(3,99002) 1, 1,     A0(1,1), B0(1,1), D0(1,1) !
         WRITE(3,99002) 1, 2,     A0(1,2), B0(1,2), D0(1,2) !
         WRITE(3,99002) 1, 3,     A0(1,3), B0(1,3), D0(1,3) !
         WRITE(3,99002) 2, 1,     A0(2,1), B0(2,1), D0(2,1) !
         WRITE(3,99002) 2, 2,     A0(2,2), B0(2,2), D0(2,2) !=======================> seraf change
         WRITE(3,99002) 2, 3,     A0(2,3), B0(2,3), D0(2,3) !
         WRITE(3,99002) 3, 1,     A0(3,1), B0(3,1), D0(3,1) !
         WRITE(3,99002) 3, 2,     A0(3,2), B0(3,2), D0(3,2) !
         WRITE(3,99002) 3, 3,     A0(3,3), B0(3,3), D0(3,3) !
c         DO i = 1, 3
c            DO j = 1, 3
c               WRITE (3, 99002) i, j, A0(i, j), B0(i, j), D0(i, j)
c            END DO
c         END DO
         WRITE (3, *)
         WRITE (3, *) '**SECTION =', isct, '   **IELM =', n
         WRITE (3, *) 'STIFFNESS MATRICES TRANSLATED'
         WRITE (3, *) '(I,J)      A0  ', '      B0  ', '        D0  '
         DO i = 1, 3
            DO j = 1, 3
               WRITE (3, 99002) i, j, A0T(i, j), B0T(i, j), D0T(i, j)
            END DO
         END DO
         WRITE (3, *)
         WRITE (3, *) 'MASS MATRICES'
         WRITE (3, *) '(I,J)      Am0  ', '     Bm0  ', '    Dm0  '
         DO i = 1, 3
            DO j = 1, 3
               WRITE (3, 99002) i, j, Am0(i, j), Bm0(i, j), Dm0(i, j)
            END DO
         END DO
         WRITE (3, *)
         WRITE (3, *) 'MATERIALs REDUCED Q MATRIX'
         DO ilam = 1, nlam
            WRITE (3, *) 'Laminate No =', ilam
            nl = nlay(ilam)
            DO ip = 1, nl
               CALL offply (im(ip,ilam), ip, angl(ip,ilam),cl,cc)
c               CALL offply1(im(ip,ilam), ip, angl(ip,ilam),cl,cc)
               WRITE (3, *) 'Ply No =', ip
               DO i = 1, 3
                  WRITE (3, 99002) i, i, cc(i, 1, ip), cc(i, 2, ip), 
     &                             cc(i, 3, ip)
               END DO
            END DO
         END DO

        p2=trim(adjustl(n_p))//"/other.dir/section_matrix_hGAST.dat"
        OPEN(UNIT=43, FILE=p2)                           !
           WRITE(43,"(a)") "stiffness matrix"            !
           WRITE(43,735) A0(3,3), A0(3,1), A0(3,2),      !
     &                   B0(3,1), B0(3,3), B0(3,2)       !
           WRITE(43,735) A0(1,3), A0(1,1), A0(1,2),      !
     &                   B0(1,1), B0(1,3), B0(1,2)       !
           WRITE(43,735) A0(2,3), A0(2,1), A0(2,2),      !
     &                   B0(2,1), B0(2,3), B0(2,2)       ! 
           WRITE(43,735) 0.00, 0.00, 0.00,               !  
     &     D0(1,1), D0(1,3), D0(1,2)                     ! 
           WRITE(43,735) 0.00, 0.00, 0.00,               !   
     &     D0(3,1), D0(3,3), D0(3,2)                     !
           WRITE(43,735) 0.00, 0.00, 0.00,               ! 
     &     D0(2,1), D0(2,3), D0(2,2)                     !
           WRITE(43,"(a)") "mass mtrix"                  !========================> seraf change 
           WRITE(43,735) Am0(3,3), Am0(3,1), Am0(3,2),   !  
     &     Bm0(3,1), Bm0(3,3), Bm0(3,2)                  !
           WRITE(43,735) Am0(1,3), Am0(1,1), Am0(1,2),   !  
     &     Bm0(1,1), Bm0(1,3), Bm0(1,2)                  ! 
           WRITE(43,735) Am0(2,3), Am0(2,1), Am0(2,2),   !    
     &     Bm0(2,1), Bm0(2,3), Bm0(2,2)                  ! 
           WRITE(43,735) Bm0(3,1), Bm0(1,1), Bm0(2,1),   !
     &     Dm0(1,1), Dm0(1,3), Dm0(1,2)                  !
           WRITE(43,735) Bm0(3,3), Bm0(1,3), Bm0(2,3),   !
     &     Dm0(3,1), Dm0(3,3), Dm0(3,2)                  !
           WRITE(43,735) Bm0(3,2), Bm0(1,2), Bm0(2,2),   !
     &     Dm0(2,1), Dm0(2,3), Dm0(2,2)                  ! 
        CLOSE (43)                                       !
 735    FORMAT(6(F30.5,3X))                              !   

        p3="/home/seraf/hGAST/pre_hGAST.dir/elstif.dat"
        open (4,file=p3,status='unknown')
        do iseg     = 1, nsegm(isct)
        write(4,*) ' ISEG = ',iseg
        do i= 1,3
        do j= 1,3
        write(4,2222) i,j,a(i,j,iseg),b(i,j,iseg),d(i,j,iseg)
        enddo
        enddo
        enddo
        close(4)
c
         WRITE (*, *) 
     &          'ENTER O to STOP or 1 to continue with strains/stresses'
         icase=0   ! READ (*, *) icase   !===============================================> seraf changes
         IF ( icase==1 ) THEN
            WRITE (*, *) '  ENTER LOADS FILE NAME'
            p4=trim(adjustl(n_p))//"/other.dir/loads.dat" !READ (*, '(a)') inpload
            OPEN (8, FILE=p4, STATUS='unknown')
c           Ordering of Loads: {Nz,Qx,Qy,My,Mx,Mz}
            DO i = 1, 6
               READ (8, *) Loads(i, isct)
            END DO
            CLOSE (8)
            
            CALL strains(isct,thks,h,zbar,A0, B0, D0, a, b, d, 
     &                   Loads(1:6,isct), EpsCur(1:6,isct), 
     &                   xpnt, ypnt, 
     &                   nskseg(isct), lamsg(1:mseg,isct), 
     &                   isgpt(1:2,1:mseg,isct), sknseg, 
     &                   nwebs(isct), webseg, nloops(isct), 
     &                   lpseg(1:mseg,1:mloops,isct), 
     &                   nlpsg(1:mloops,isct), rshr, 
     &                   nlay, hl, angl, aahstr, cl, cc, ezz, ezs, 
     &                   ezj, szz, szs, szj, im, strength,
     &                   fcrit1,fcrit2, corweb, n_p)
            write(*,*)'FAILURE if CRIT 1 >0 ',fcrit1
            write(*,*)'FAILURE if CRIT 2 <1 ',fcrit2   
            DO i = 1, 6
               WRITE (*, *) EpsCur(i, isct)
            END DO 
         END IF
             
c
      END DO
      CLOSE (1)
      CLOSE (2)
      CLOSE (3)
99002 FORMAT (1x, '(', i1, ',', i1, ')', 1x, 3(',',E12.4))
 2222 FORMAT (2i4,3(2x,E12.4))
      END PROGRAM BM1_FE_DMP
c
 
 
 
c --------------------------------------------------------------------
c        scope: this subroutine reads input required for laminate definition
c                - problem parameters
c                - material properties
c                - laminate configuration
c                - forms property matrices on material coords.
c
c
      SUBROUTINE laminp(nmat, nlam, nlay, pm, cl, thk, 
     &                  im, hl, angl, strength)
      IMPLICIT NONE      
      include '/home/seraf/hGAST/pre_hGAST.dir/paramet.for'
      INTEGER I, IL, ILAM, IP, J, K, KLAM, NLAM, NMAT
c
      INTEGER im(nlayer, mlamin),nlay(mlamin)
      REAL*8  pm(mater, mater)
      REAL*8  cl(5,  5, mater)
      REAL*8  thk(mlamin),hl(nlayer,mlamin),angl(nlayer, mlamin)
      REAL*8  strength(mater, mater)
c
c     .. read material properties (in material system 0x1x2x3))
c
      READ (2, *)
c     pm(1) - rho
c     pm(2) - Yl11
c     pm(3) - yl22
c     pm(4) - nul12
c     pm(5) - Yl12
c     pm(6) - Yl23
c     pm(7) - Yl33
c     pm(8) - nul23
c     pm(9) - Yl13
c     pm(10) - nul13
c
c     strength(1)= Tension strength in direction 1 (along the fibres)   :  XT 
c     strength(2)= Compressive strength in direction 1 (along the fibres): XC 
c     strength(3)= Tension strength in direction 2 (Transverse to the fibres): YT
c     strength(4)= Compressive strength in direction 2 (Transverse to the fibres): YC 
c     strength(5)= Tension strength in direction 3 (through the thickness): ZT can be neglected 
c     strength(6)= Compressive strength in direction 3 (Through the thickness): ZC ' can be neglected 
c     strength(7)= Transverse Shear strength through the thickness: S23 
c     strength(8)= Transverse Shear strength: S13 
c     strength(9)= Shear strengh in plane: S12
c
      DO k = 1, nmat
         READ (2, *) i, pm(1, i), pm(2, i), pm(3, i), pm(4, i), pm(5, i)
         READ (2, *) pm(6, i), pm(7, i), pm(8, i), pm(9, i), pm(10, i)
         READ (2, *) (strength(j,i), j=1, 5)
         READ (2, *) (strength(j,i), j=6, 9)
         READ (2, *)
      END DO
c
      CALL mtcnst(nmat, pm, cl)
c
      READ (2, *)
      DO klam = 1, nlam
         READ (2, *) ilam, nlay(ilam)
         thk(ilam) = 0.
         DO il = 1, nlay(ilam)
            READ (2, *) ip, im(ip, ilam), hl(ip, ilam), angl(ip, ilam)
            thk(ilam) = thk(ilam) + hl(ip, ilam)
         END DO
      END DO
c
      END SUBROUTINE LAMINP
c --------------------------------------------------------------------
c        scope:  reads input required for z-section definition
c        status: develop
c                - subsections, laminations
c                - reads geometry configuration of each subsection
c                - forms parametric representation of section geometry
c
      SUBROUTINE sctinp(nsct  , nsegm, nskseg, lamsg , isgpt, 
     &                  sknseg, nwebs, webseg, nloops, 
     &                  lpseg , nlpsg)
      IMPLICIT NONE
      include '/home/seraf/hGAST/pre_hGAST.dir/paramet.for'
c
      INTEGER nsegm (msct),nskseg(msct),isgpt(2, mseg, msct), 
     &        lamsg (mseg, msct)
      INTEGER sknseg(mseg, msct), webseg(mwebs, msct), nwebs(msct)
      INTEGER lpseg (mseg, mloops, msct), nloops(msct), 
     &        nlpsg (mloops, msct)
      INTEGER NSCT
      INTEGER ILOOP, ISCT, ISEG, IWEB, JSKN, K, KSCT, KSEG 
      CHARACTER*10 dummy10
c
      READ (2, *)
      READ (2, *) nsct
                      ! No of sections
      DO ksct = 1, nsct
         READ (2, *) dummy10, isct, nsegm(isct)
         READ (2, *) dummy10, nskseg(isct)
         DO jskn = 1, nskseg(isct)
            READ (2, *) k, isgpt(1, k, isct), isgpt(2, k, isct), 
     &                  lamsg(k, isct)
            sknseg(jskn, isct) = k
         END DO
         READ (2, *) dummy10, nwebs(isct)
         DO iweb = 1, nwebs(isct)
            READ (2, *) k, isgpt(1, k, isct), isgpt(2, k, isct), 
     &                  lamsg(k, isct)
            webseg(iweb, isct) = k
         END DO
         READ (2, *) dummy10, nloops(isct)
         DO iloop = 1, nloops(isct)
            READ (2, *) k, kseg
            READ (2, *) (lpseg(iseg,k,isct), iseg=1, kseg)
            nlpsg(k, isct) = kseg
         END DO
      END DO
c
      END SUBROUTINE SCTINP
c
c  ---------------------------------------------------------------------
      SUBROUTINE input(nsct,npts,xpnt,ypnt)
      IMPLICIT NONE
      include '/home/seraf/hGAST/pre_hGAST.dir/paramet.for'
           
      INTEGER nsct,npts
      REAL*8  xpnt(maxpoint,msct),ypnt(maxpoint,msct)
      INTEGER I, ISCT, K, KSCT
 
      READ (1, *)
      READ (1, *)              npts
      WRITE(*, *)  ' NPTS = ', npts
      READ (1, *)
 
      WRITE(*, *)  ' NSCT = ', nsct
      DO isct = 1,   nsct
      READ (1, *)    ksct
      DO i = 1,      npts
      READ (1, *) k, xpnt(k, ksct), ypnt(k, ksct)
      END DO
      END DO
 
      END SUBROUTINE INPUT
c  ------------------------------------------------------------------
c               pm(1) - rho
c               pm(2) - Yl11
c               pm(3) - yl22
c               pm(4) - nul12
c               pm(5) - Yl12
c               pm(6) - Yl23
c               pm(7) - Yl33
c               pm(8) - nul23
c               pm(9) - Yl13
c               pm(10) - nul13
      SUBROUTINE mtcnst(nmat, pm, CL)
      IMPLICIT NONE
      include '/home/seraf/hGAST/pre_hGAST.dir/paramet.for'
      INTEGER IMAT,NMAT
      REAL*8 pm(mater, mater)
      REAL*8 CL(5,  5, mater)
c
      DO imat = 1, nmat
c        ..formulate matrices in material system
c        ... CL stands for compliance matrix
         CL(1, 1, imat) = 1./pm(2, imat)
         CL(1, 2, imat) = -pm(4, imat)/pm(2, imat)
         CL(2, 1, imat) =    CL(1, 2, imat)
         CL(2, 2, imat) = 1./pm(3, imat)
         CL(3, 3, imat) = 1./pm(6, imat)
         CL(4, 4, imat) = 1./pm(9, imat)
         CL(5, 5, imat) = 1./pm(5, imat)
      END DO
c
      END SUBROUTINE MTCNST
c  -------------------------------------------------------------
c     Scope: Calculation of off-axis ply moduli
c            this version is limited to Yc11, Yc13 only
c
      SUBROUTINE offply(imat, ip, ang, cl, cc)
      IMPLICIT NONE
      include '/home/seraf/hGAST/pre_hGAST.dir/paramet.for'
      REAL*8  ANG
      INTEGER IMAT,IP
c
      REAL*8 cl(5, 5 , mater), cc(3,3,nlayer)
      REAL*8 R2(5, 5), R2INV(5, 5), R2TR(5, 5), R2INTR(5, 5), w2(5, 5), 
     &       R1(3, 3), R1TR (3, 3), w1  (3, 5), w3    (3, 3), sc(5, 5)
      REAL*8 DD
c
      CALL TMTR2D(ang, R2, R2INV, R2TR, R2INTR, R1, R1TR)
c     ...OFF-AXIS PLY Compliance [S]c (5x5)
      CALL MULT(R2TR, CL(1:5,1:5,IMAT), W2, 5, 5, 5)
      CALL MULT(W2, R2, SC, 5, 5, 5)
c     equivalent off-axis stiffness [C]c (3x3)
      dd =   sc(1, 1)*sc(5, 5)-sc(1, 5)**2
      cc(1, 1, ip) =  sc(5, 5)/dd
      cc(3, 3, ip) =  sc(1, 1)/dd
      cc(1, 3, ip) = -sc(1, 5)/dd
      cc(3, 1, ip) =  cc(1, 3, ip)
      cc(2, 2, ip) =  1./ sc(4, 4)
c
      END SUBROUTINE OFFPLY
      
c  -------------------------------------------------------------
c     Scope: Calculation of off-axis ply moduli
c            this version is limited to Yc11, Yc13 only
c
      SUBROUTINE offply1(imat, ip, ang, cl, cc)
      IMPLICIT NONE
      include '/home/seraf/hGAST/pre_hGAST.dir/paramet.for'
      REAL*8  ANG
      INTEGER IMAT,IP,I,J
      INTEGER kpvt(5), info, inert(3)
c
      REAL*8 cl(5, 5 , mater), cc(3,3,nlayer)
      REAL*8 R2(5, 5), R2INV(5, 5), R2TR(5, 5), R2INTR(5, 5), w2(5, 5), 
     &       R1(3, 3), R1TR (3, 3), w1  (3, 5), w3    (3, 3), sc(5, 5)
      REAL*8 DD
      REAL*8 wk(5),det1(2)
c
      CALL TMTR2D(ang, R2, R2INV, R2TR, R2INTR, R1, R1TR)
c     ...OFF-AXIS PLY Compliance [S]c (5x5)
      CALL MULT(R2TR, CL(1:5,1:5,IMAT), W2, 5, 5, 5)
      CALL MULT(W2, R2, SC, 5, 5, 5)
           
      CALL dsifa(SC, 5, 5, kpvt, info)
      CALL dsidi(SC, 5, 5, kpvt, det1, inert, wk, 001)
c
c     expand lower triangle
      DO i = 2, 5
         DO j = 1, i - 1
            SC(i, j) = SC(j, i)
         END DO
      END DO
c     equivalent off-axis stiffness [C]c (3x3)
      cc(1, 1, ip) =  sc(1 , 1)
      cc(1, 2, ip) =  sc(1 , 4)
      cc(1, 3, ip) =  sc(1 , 5)
      cc(2, 2, ip) =  sc(4 , 4)
      cc(2, 3, ip) =  sc(4 , 5)
      cc(3, 3, ip) =  sc(5 , 5)
      
      cc(2, 1, ip) =   cc(1, 2, ip)
      cc(3, 1, ip) =   cc(1, 3, ip)
      cc(3, 2, ip) =   cc(2, 3, ip)
c
      END SUBROUTINE OFFPLY1
c ------------------------------------------------------------------
 
c     SCOPE   : 3-D TRANSFORMATION MATRICES
c
      SUBROUTINE TMTR2D(THETA, T2, T1, T3, T4, R2, R1)
      IMPLICIT NONE
      REAL*8 C, CC, CS, CS2, PI, S, SS
      REAL*8 THETA
      REAL*8 T1(5, 5), T2(5, 5), T3(5, 5), T4(5, 5)
      REAL*8 R1(3, 3), R2(3, 3)
      INTEGER II, J, JJ, K
c
      PI = 4.D0*datan(1.D0)
      C = COS(THETA*PI/180.0)
      S = SIN(THETA*PI/180.0)
c
      CC =  C*C
      SS =  S*S
      CS =  C*S
      CS2 = 2.*CS
c
      T1(1, 1) =  CC
      T1(1, 2) =  SS
      T1(1, 5) = -CS2
      T1(2, 1) =  SS
      T1(2, 2) =  CC
      T1(2, 5) =  CS2
      T1(3, 3) =  C
      T1(3, 4) =  S
      T1(4, 3) = -S
      T1(4, 4) =  C
      T1(5, 1) =  CS
      T1(5, 2) = -CS
      T1(5, 5) =  CC - SS
c
      R1(1, 1) =  C
      R1(2, 2) =  C
      R1(1, 2) =  S
      R1(2, 1) = -S
      R1(3, 3) =  1.
c
      DO II = 1, 5
      DO JJ = 1, 5
      T2(II, JJ) = T1(II, JJ)
      END DO
      END DO
c
      T2(1, 5) =  CS2
      T2(2, 5) = -CS2
      T2(3, 4) = -S
      T2(4, 3) =  S
      T2(5, 1) = -CS
      T2(5, 2) =  CS
 
      R2(1, 1) =  C
      R2(2, 2) =  C
      R2(1, 2) = -S
      R2(2, 1) =  S
      R2(3, 3) =  1.
c     T3 IS THE TRANSPOSE OF T2 (Transp Rsigma)
c     T4 IS THE TRANSPOSE OF T1 (Transp Rsigma**-1)
      DO J = 1, 5
      DO K = 1, 5
      T3(J, K) = T2(K, J)
      T4(J, K) = T1(K, J)
      END DO
      END DO
c
      END SUBROUTINE TMTR2D
c
c  ----------------------------------------------------------------
c
      SUBROUTINE MULT(AAA, BCB, CCC, M, K, N)
      IMPLICIT NONE
      INTEGER I, J, K, L, M, N
      REAL*8 AAA(M, K), BCB(K, N), CCC(M, N)
c
c     MULTIPLICATION OF MATRICES
c     CCC = AAA * BCB
c
      DO I = 1, M
         DO J = 1, N
            CCC(I, J) = 0.
         END DO
      END DO
c
      DO L = 1, K
         DO I = 1, M
            DO J = 1, N
               CCC(I, J) = CCC(I, J) + AAA(I, L)*BCB(L, J)
            END DO
         END DO
      END DO
      END SUBROUTINE MULT
c
c  ------------------------------------------------
c
c    Scope: Calculation of laminate shell matrices - stiffness
c           Has the structure to include sublaminates as a discrete layer.
c
c           a, b, d:      extension-shear, coupling, flexure matrices
c
      SUBROUTINE lamstf(NL, thk, h, zbar, CC, w4, w5, A, B, D)
      IMPLICIT NONE
      include '/home/seraf/hGAST/pre_hGAST.dir/paramet.for'
      INTEGER I, K, L
      REAL*8  THK
c     
c     NOTE:  ONLY M>= J  PART IS CALCULATED
 
      INTEGER NL
      REAL*8 cc(3,3,nlayer)
      REAL*8 h (nlayer), zbar(nlayer)
      REAL*8 w4(nlayer),   w5(nlayer)
      REAL*8 a(3, 3), D(3, 3), B(3, 3)
c
      DO k = 1, 3
         DO l = 1, 3
            A(k, l) = 0.
            B(k, l) = 0.
            D(k, l) = 0.
         END DO
      END DO
 
      zbar(1) = (h(1)-thk)/2.
      DO   i  =  2, NL
      zbar(i) =  zbar(i-1) + (h(i-1)+h(i))/2.
      END DO
c
c     ---------------------------------------------------------
c     calculation of stiffness matrices starts here
c     ---------------------------------------------------------
c
      DO i = 1, NL
         w4(i) = h(i)*zbar(i)
         w5(i) = h(i)*(h(i)**2/12.+zbar(i)**2)
c        LAMINATE A, B, D MATRICES
c
         A(1, 1) = A(1, 1) + CC(1, 1, I)*h (i)
         B(1, 1) = B(1, 1) + CC(1, 1, I)*w4(i)
         D(1, 1) = D(1, 1) + CC(1, 1, I)*w5(i)
c
         A(1, 2) = A(1, 2) + CC(1, 2, I)*h (i)
         B(1, 2) = B(1, 2) + CC(1, 2, I)*w4(i)
         D(1, 2) = D(1, 2) + CC(1, 2, I)*w5(i)
c        
         A(1, 3) = A(1, 3) + CC(1, 3, I)*h (i)
         B(1, 3) = B(1, 3) + CC(1, 3, I)*w4(i)
         D(1, 3) = D(1, 3) + CC(1, 3, I)*w5(i)
c       
         A(2, 2) = A(2, 2) + CC(2, 2, I)*h (i)
         B(2, 2) = B(2, 2) + CC(2, 2, I)*w4(i)
         D(2, 2) = D(2, 2) + CC(2, 2, I)*w5(i)
c       
         A(2, 3) = A(2, 3) + CC(2, 3, I)*h (i)
         B(2, 3) = B(2, 3) + CC(2, 3, I)*w4(i)
         D(2, 3) = D(2, 3) + CC(2, 3, I)*w5(i)
c       
         A(3, 3) = A(3, 3) + CC(3, 3, I)*h (i)
         B(3, 3) = B(3, 3) + CC(3, 3, I)*w4(i)
         D(3, 3) = D(3, 3) + CC(3, 3, I)*w5(i)
c        ... lower triangular part...
         A(2, 1) = A(2, 1) + CC(2, 1, I)*h (i)
         B(2, 1) = B(2, 1) + CC(2, 1, I)*w4(i)
         D(2, 1) = D(2, 1) + CC(2, 1, I)*w5(i)
c       
         A(3, 1) = A(3, 1) + CC(3, 1, I)*h (i)
         B(3, 1) = B(3, 1) + CC(3, 1, I)*w4(i)
         D(3, 1) = D(3, 1) + CC(3, 1, I)*w5(i)
c      
         A(3, 2) = A(3, 2) + CC(3, 2, I)* h(i)
         B(3, 2) = B(3, 2) + CC(3, 2, I)*w4(i)
         D(3, 2) = D(3, 2) + CC(3, 2, I)*w5(i)
c
      END DO
 
      END SUBROUTINE LAMSTF
c
c  ------------------------------------------------------------------
c     Scope: Calculation of laminate matrices - specific mass
c
c            am, bm, dm - extensional, coupling, and rotational specific
c                         (per unit area) inertia matrices
c
      SUBROUTINE lamdns(NL, thk, h, zbar, rho, W1, W3, AM, BM, DM)
      IMPLICIT NONE
      include '/home/seraf/hGAST/pre_hGAST.dir/paramet.for'
c
      INTEGER NL
      REAL*8  rho(nlayer)
      REAL*8  zbar(nlayer), h(nlayer)
      REAL*8  W1(nlayer)  ,W3(nlayer)
      REAL*8  AM, BM, DM
      REAL*8  AA, THK
      INTEGER I
c
      am = 0.
      bm = 0.
      dm = 0.
c
c     =======================================
c     ... laminate mass
c     =======================================
c
      DO i  = 1, NL
      w1(i) = h(i)*zbar(i)
      w3(i) = h(i)*(h(i)**2/12.+zbar(i)**2)
      aa = 1.
      am = am + rho(i)*h (i)
      bm = bm + rho(i)*w1(i)
      dm = dm + rho(i)*w3(i)
      END DO
c
      END SUBROUTINE LAMDNS
 
c  -----------------------------------------------
 
c     Scope: Calculation of beam section matrices - stiffness
c
c           a0, b0, d0: section extension, coupling, flexure matrices -- see notes
c                       wr. to the y, z axes at center of cross section
c           a, b, d:       extension, coupling, flexure matrices of face (skin) laminates
c                       about local centeline (of skin)
c                       at curvilinear cordinates (x, s, zeta) -- see notes....
c
      SUBROUTINE secstf(isct, thks, h, zbar, a, b, d, as, A0, B0, D0 , 
     &                  xpnt, ypnt, 
     &                  nskseg, lamsg, isgmt, sknseg, 
     &                  nwebs, webseg, nloop, lpseg, nlpsg, 
     &                  aahstr, area0, aa, area, corweb, n_p)
      IMPLICIT NONE
      include '/home/seraf/hGAST/pre_hGAST.dir/paramet.for'
      INTEGER I, IL, ILOOP, INFO, IP, IP1, IP2, IS, ISCT, ISEG, ISG, 
     &        IWEB, J, LPS
c
      INTEGER nskseg , isgmt(2,mseg),lamsg(mseg)
      INTEGER sknseg(mseg, msct)  , webseg(mwebs, msct), nwebs
      INTEGER lpseg (mseg, mloops), nloop, nlpsg(mloops),kpvt(4)
c
      REAL*8 zbar(nlayer), h(nlayer), ploop(mseg, 3)
      REAL*8 A0(3, 3), D0(3, 3), B0(3, 3), thks(mseg)
      REAL*8 A(3,3,mseg), D(3,3,mseg), B(3,3,mseg), as(2,2,mseg)
      REAL*8 G1, G2
      REAL*8 ro(2),ros(2),ross(2),area0(mloops),aa(mseg),aah, 
     &       aahstr(mloops),rzeta,rs,area(mseg),aaloop(4),aaa(4, 4)
      REAL*8 areb(mseg),areb0(mloops)
      REAL*8 xpnt(maxpoint, msct), ypnt(maxpoint, msct)
      REAL*8 x1, x2, y1, y2, dss, thk
      REAL*8 corweb
      CHARACTER*300 n_p
c
      open(5,file=trim(n_p)//"/other.dir/test.txt",status='unknown')
c
      DO i = 1, 3
         DO j = 1, 3
            a0(i, j) = 0.
            b0(i, j) = 0.
            d0(i, j) = 0.
         END DO
      END DO
c
      DO        iseg   = 1, mseg
         aa    (iseg)  = 0.
         area  (iseg)  = 0.
         areb  (iseg)  = 0.
      END DO
      DO        iloop  = 1, nloop
         area0 (iloop) = 0.
         areb0 (iloop) = 0.
         aaloop(iloop) = 0.
      END DO
c
c     ... calculate Ao, A_h for each section & loop
c
      DO isg = 1, nskseg
         thk = thks(isg)
         ip1 = isgmt(1, isg)
         ip2 = isgmt(2, isg)
         DO ip = ip1, ip2 - 1
            x1 = xpnt(ip  , isct)
            x2 = xpnt(ip+1, isct)
            y1 = ypnt(ip  , isct)
            y2 = ypnt(ip+1, isct)
            CALL mdllin(x1, y1, x2, y2, dss, ro, ros, thk, 1)
            rzeta = -ro(1)*ros(2) + ro(2)*ros(1)
            rs    =  ro(1)*ros(1) + ro(2)*ros(2)
            area(isg) = area(isg) + 1.0D0*rzeta*dss 
     >                            + 2.0D0*b(3,3,isg)/a(3,3,isg)*dss
            areb(isg) = areb(isg) + 2.0D0*b(3,3,isg)/a(3,3,isg)*dss
            aa  (isg) = aa  (isg) + 1.0D0/a(3,3,isg)*dss
         END DO
      END DO
c
c     integration along webs
c
      IF ( nwebs/=0 ) THEN
         DO iweb = 1, nwebs
            iseg = webseg(iweb, isct)
            thk = thks(iseg)
            ip1 = isgmt(1, iseg)
            ip2 = isgmt(2, iseg)
            x1 = xpnt(ip1, isct)
            x2 = xpnt(ip2, isct)
            y1 = ypnt(ip1, isct)
            y2 = ypnt(ip2, isct)
            CALL mdllin(x1, y1, x2, y2, dss, ro, ros, thk, 1)
            dss   =  dss  - corweb
            rzeta = -ro(1)*ros(2) + ro(2)*ros(1)
            rs    =  ro(1)*ros(1) + ro(2)*ros(2)
            area(iseg) = area(iseg) + 1.0D0*rzeta*dss
     >                              + 2.0D0*b(3,3,iseg)/a(3,3,iseg)*dss
            aa  (iseg) = aa  (iseg) + 1.0D0/a(3,3,iseg)*dss
         END DO
      END IF
c
      DO il = 1, nloop
         DO is = 1, nlpsg(il)
            lps = abs(lpseg(is,il))
            area0 (il) = area0 (il) + area(lps)*lpseg(is, il)/lps      
            areb0 (il) = areb0 (il) + areb(lps)*lpseg(is, il)/lps
            aaloop(il) = aaloop(il) + aa  (lps)
      write(5,*) il,lpseg(is,il),area(lps)*lpseg(is, il)/lps,
     &                           areb(lps)*lpseg(is, il)/lps
         END DO
      END DO
      write(5,*)
c
c     ... calculate equivalent Ao*/A_h* for each loop
c
      IF ( nwebs==0 ) aahstr(1) = area0(1)/aaloop(1)
      IF ( nwebs==1 ) THEN
         aaa(1, 1) =  aaloop(1)
         aaa(2, 2) =  aaloop(2)
         aaa(1, 2) = -aa(webseg(1,isct))
         aaa(2, 1) =  aaa(1, 2)
         aahstr(1) =  area0(1)
         aahstr(2) =  area0(2)
c        ...calculate equivalent aahstr for each loop
         CALL dsifa(aaa, 4, nloop, kpvt, info)
         CALL dsisl(aaa, 4, nloop, kpvt, aahstr)
      END IF
      IF ( nwebs==2 ) THEN
         aaa(1, 1) =  aaloop(1)
         aaa(2, 2) =  aaloop(2)
         aaa(3, 3) =  aaloop(3)
         aaa(1, 2) = -aa(webseg(1,isct))
         aaa(2, 1) =  aaa(1, 2)
         aaa(2, 3) = -aa(webseg(2,isct))
         aaa(3, 2) =  aaa(2, 3)
         aahstr(1) =  area0(1)
         aahstr(2) =  area0(2)
         aahstr(3) =  area0(3)
c        ...calculate equivalent aahstr for each loop
         write(5,*)' AREA1, AREA2, AREA3 '
         write(5,2224) area0(1),area0(2),area0(3)
         write(5,*)' AREB1, AREB2, AREB3 '
         write(5,2224) areb0(1),areb0(2),areb0(3)
         write(5,*)' Lamda 11, 12, 13 '
         write(5,2224) aaa(1,1),aaa(1,2),aaa(1,3)
         write(5,*)' Lamda 22, 23, 23 '
         write(5,2224) aaa(2,2),aaa(2,3),aaa(3,3)
         CALL dsifa(aaa, 4, nloop, kpvt, info)
         CALL dsisl(aaa, 4, nloop, kpvt, aahstr)
          write(5,*)' Cha1, Cha2, Cha3 '
         write(5,2224) aahstr(1),aahstr(2),aahstr(3)
         write(5,*)
      END IF
c
c==============
c     integration along segments
      DO isg = 1, nskseg
         thk = thks(isg)
c        ... determine the loop where the segment       belongs
         DO il = 1, nloop
            DO is = 1, nlpsg(il)
               IF ( isg==abs(lpseg(is,il)) ) iloop = il
            END DO
         END DO
 
         ip1 = isgmt(1, isg)
         ip2 = isgmt(2, isg)
         DO ip = ip1, ip2 - 1
            x1 = xpnt(ip, isct)
            x2 = xpnt(ip+1, isct)
            y1 = ypnt(ip, isct)
            y2 = ypnt(ip+1, isct)
            CALL mdllin(x1, y1, x2, y2, dss, ro, ros, thk, 1)
            rzeta = -ro(1)*ros(2) + ro(2)*ros(1)
            rs    =  ro(1)*ros(1) + ro(2)*ros(2)
            aah   =  aahstr(iloop)/a(3,3,isg)-2.d0*b(3,3,isg)/a(3,3,isg)
            write(5,2223) isg,ro(1),ro(2),dss,ros(1),ros(2),aah
c           calculate stifness terms
            CALL stfmtr(dss, ro, ros, ross, rzeta, rs, aah, 
     &                  a(1:3,1:3,isg), b(1:3,1:3,isg), d(1:3,1:3,isg), 
     &                  a0, b0, d0)
         END DO
      END DO
c
c==============
c     integration along webs
 
      IF ( nwebs/=0 ) THEN
         DO iweb = 1, nwebs
            iseg = webseg(iweb, isct)
            thk = thks(iseg)
            ip1 = isgmt(1, iseg)
            ip2 = isgmt(2, iseg)
            x1 = xpnt(ip1, isct)
            x2 = xpnt(ip2, isct)
            y1 = ypnt(ip1, isct)
            y2 = ypnt(ip2, isct)
c  TAKIS ATTENTION (iweb 1 sign changed)
            IF ( iweb==1 ) aah =-(aahstr(1)-aahstr(2))/a(3,3,iseg)
     &                               -2.d0*b(3,3,iseg)/a(3,3,iseg)
            IF ( iweb==2 ) aah = (aahstr(2)-aahstr(3))/a(3,3,iseg)
     &                               -2.d0*b(3,3,iseg)/a(3,3,iseg)
 
            CALL mdllin(x1, y1, x2, y2, dss, ro, ros, thk, 1)
            dss   =  dss  - corweb
            rzeta = -ro(1)*ros(2) + ro(2)*ros(1)
            rs    =  ro(1)*ros(1) + ro(2)*ros(2)
            write(5,2223) iseg,ro(1),ro(2),dss,ros(1),ros(2),aah
c           calculate stifness terms
            CALL stfmtr(dss, ro, ros, ross, rzeta, rs, aah, 
     &                a(1:3,1:3,iseg), b(1:3,1:3,iseg), d(1:3,1:3,iseg), 
     &                a0, b0, d0)
         END DO
      END IF
c     ...symmetric terms
      a0(3, 2) = a0(2, 3)
      a0(2, 1) = a0(1, 2)
      a0(3, 1) = a0(1, 3)
      d0(2, 1) = d0(1, 2)
      d0(3, 1) = d0(1, 3)
      d0(3, 2) = d0(2, 3)
c
 2223 format(i4,6(1x,e15.6))
 2224 format(3(1x,e15.6))
      close(5)
      END SUBROUTINE SECSTF
c ------------------------------------------------
c     linear calculation of geometric data, derivatives, length etc
      SUBROUTINE mdllin(x1, y1, x2, y2, ds, ro, ros, thk, mode)
      IMPLICIT NONE
      REAL*8 x1, y1, x2, y2, ds, ro(2), ros(2)
      REAL*8 dx, dy, thk, r
      INTEGER mode
      dx = x2 - x1
      dy = y2 - y1
      r=abs(1.0-thk/(2.0*dsqrt(((x1+x2)/2.0)**2.0+((y1+y2)/2.0)**2.0)))
      ds=dsqrt(dx*dx+dy*dy) !=================================================================================> seraf change
      ro(1) = (x1+x2)*.5D0 + .5D0*thk*dy/ds*mode
      ro(2) = (y1+y2)*.5D0 - .5D0*thk*dx/ds*mode
      ros(1) = dx/ds
      ros(2) = dy/ds
      END SUBROUTINE MDLLIN
c ------------------------------------------------
      SUBROUTINE stfmtr(dss, ro, ros, ross, rzeta, rs, aah, a, b, d, a0, 
     &                  b0, d0)
      IMPLICIT NONE
c
      REAL*8 A0(3, 3), D0(3, 3), B0(3, 3)
      REAL*8 A (3, 3), D (3, 3), B (3, 3)
      REAL*8 ro(2), ros(2), ross(2), aah, rzeta, rs
      REAL*8 dss
c
      a0(1, 1) = a0(1, 1) + dss*a(1, 1)
      a0(2, 2) = a0(2, 2) + dss*(a(3,3)*ros(1)**2+a(2,2)*ros(2)**2)
      a0(3, 3) = a0(3, 3) + dss*(a(3,3)*ros(2)**2+a(2,2)*ros(1)**2)
      a0(2, 3) = a0(2, 3) + dss*(a(3,3)*ros(1)*ros(2)-a(2,2)*ros(1)
     &           *ros(2))
c     ... Ext-Shr coupling
      a0(1, 2) = a0(1, 2) + dss*a(1, 3)*ros(1)
      a0(1, 3) = a0(1, 3) + dss*a(1, 3)*ros(2)
c
      b0(1, 1) = b0(1, 1) + dss*(-a(1,1)*ro(1)+b(1,1)*ros(2))
      b0(1, 2) = b0(1, 2) + dss*( a(1,1)*ro(2)+b(1,1)*ros(1))
      b0(2, 3) = b0(2, 3) + dss*(-a(3,3)*aah-2.*b(3,3))*ros(1)
      b0(3, 3) = b0(3, 3) + dss*(-a(3,3)*aah-2.*b(3,3))*ros(2)
c     ... 16 coupling terms
      b0(1, 3) = b0(1, 3) + dss*(-a(1,3)*aah-2.*b(1,3))
      b0(2, 1) = b0(2, 1) + dss*(-a(1,3)*ro(1)+b(1,3)*ros(2))*ros(1)
      b0(3, 1) = b0(3, 1) + dss*(-a(1,3)*ro(1)+b(1,3)*ros(2))*ros(2)
      b0(2, 2) = b0(2, 2) + dss*( a(1,3)*ro(2)+b(1,3)*ros(1))*ros(1)
      b0(3, 2) = b0(3, 2) + dss*( a(1,3)*ro(2)+b(1,3)*ros(1))*ros(2)
c     ... Flexure
      d0(1, 1) = d0(1, 1) + dss*(a(1,1)*ro(1)**2-2*b(1,1)*ro(1)*ros(2)
     &           +d(1,1)*ros(2)**2)
      d0(2, 2) = d0(2, 2) + dss*(a(1,1)*ro(2)**2+2*b(1,1)*ro(2)*ros(1)
     &           +d(1,1)*ros(1)**2)
      d0(1, 2) = d0(1, 2)
     &           + dss*(-a(1,1)*ro(1)*ro(2)-b(1,1)*(ro(1)*ros(1)-ro(2)
     &           *ros(2))+d(1,1)*ros(1)*ros(2))
c
      d0(3, 3) = d0(3, 3) + dss*(a(3,3)*aah**2+4.*b(3,3)*aah+4.*d(3,3))
c     ... B-Tr coupling
      d0(1, 3) = d0(1, 3)
     &           + dss*(a(1,3)*aah*ro(1)-b(1,3)*(aah*ros(2)-2.*ro(1))
     &           -2.*d(1,3)*ros(2))
      d0(2, 3) = d0(2, 3)
     &           + dss*(-a(1,3)*aah*ro(2)-b(1,3)*(aah*ros(1)+2.*ro(2))
     &           -2.*d(1,3)*ros(1))
 
c
      END SUBROUTINE STFMTR
c  ------------------------------------------------------------------
c
c     Scope: Calculation of hollow beam section matrices - mass density
c
c            am0, bm0, dm0:       section extension, coupling, flexure matrices -- see notes
c                              wr. to the y, z axes at center of cross section
c            am, bm, dm:           extension, coupling, flexure matrices of face (skin) laminates
c                              about local centeline (of skin)
c                              at curvilinear cordinates (x, s, zeta) -- see notes....
c
      SUBROUTINE secdns(isct,thks,h,zbar, ama, bma, dma, am0, bm0, dm0, 
     &                  xpnt, ypnt, 
     &                  nskseg , lamsg , isgmt , sknseg,
     &                  nwebs, webseg, nloop , lpseg , nlpsg , 
     &                  aahstr, area0, aa, corweb)
      IMPLICIT NONE
      include '/home/seraf/hGAST/pre_hGAST.dir/paramet.for'
      INTEGER I, IL, ILOOP, IP, IP1, IP2, IS, ISCT, ISEG, ISG, IWEB, J
      REAL*8 S_SEGA, S_WEBA
c
      REAL*8 zbar(nlayer),h(nlayer)
      REAL*8 AM0(3, 3), DM0(3, 3), BM0(3, 3), thks(mseg) 
      REAL*8 ama(mseg), dma(mseg), bma(mseg), AM, DM, BM
      REAL*8 ro(2), ros(2), ross(2), rzeta, rs
      REAL*8 xpnt(maxpoint, msct), ypnt(maxpoint, msct)
      REAL*8 x1, y1, x2, y2, dss, thk
      REAL*8 area0(mloops),aa(mseg), aahstr(mloops),darea0,daa,psio, 
     &       psios(maxpoint), aah      
      REAL*8 corweb
c
      INTEGER nskseg, isgmt(2, mseg), lamsg(mseg)
      INTEGER sknseg(mseg , msct)  , webseg(mwebs, msct), nwebs
      INTEGER lpseg (mseg , mloops), nloop, nlpsg(mloops)
c
      DO i = 1, 3
         DO j = 1, 3
            am0(i, j) = 0.
            bm0(i, j) = 0.
            dm0(i, j) = 0.
         END DO
      END DO
c     ..integrate over skin mid-line
      psios(1) = 0.
c
      DO isg = 1, nskseg
         thk = thks(isg)
         s_segA = 0.
c        ... determine the loop where the segment belongs
         DO il = 1, nloop
            DO is = 1, nlpsg(il)
               IF ( isg==abs(lpseg(is,il)) ) iloop = il
            END DO
         END DO
c
         am = ama(isg)
         dm = dma(isg)
         bm = bma(isg)
c
         ip1 = isgmt(1, isg)
         ip2 = isgmt(2, isg)
         DO ip = ip1, ip2 - 1
            x1 = xpnt(ip, isct)
            x2 = xpnt(ip+1, isct)
            y1 = ypnt(ip, isct)
            y2 = ypnt(ip+1, isct)
c
            darea0 = 0.
            daa = 0.
            CALL mdllin(x1, y1, x2, y2, dss, ro, ros, thk, 1)
            rzeta = -ro(1)*ros(2) + ro(2)*ros(1)
            rs    =  ro(1)*ros(1) + ro(2)*ros(2)
            s_segA = s_segA + dss
c
            am0(1, 1) = am0(1, 1) + dss*  am
            bm0(1, 1) = bm0(1, 1) + dss*(-am*ro(1)+bm*ros(2)) !                   !bm0(1, 1) = bm0(1, 1) + dss*(-am*ro(1)+bm*ros(2))
            bm0(1, 2) = bm0(1, 2) + dss*( am*ro(2)+bm*ros(1)) !===> seraf changes !bm0(1, 2) = bm0(1, 2) + dss*( am*ro(2)+bm*ros(1))
            bm0(2, 3) = bm0(2, 3) + dss*(-am*ro(2)-bm*ros(1)) !                   !bm0(2, 3) = bm0(2, 3) + dss*(-am*rzeta-bm)
            bm0(3, 3) = bm0(3, 3) + dss*( am*ro(1)-bm*ros(2)) !                   !bm0(3, 3) = bm0(3, 3) + dss*( am*rs) 
c
            dm0(1, 1) = dm0(1, 1)
     &                  + dss*(am*ro(1)**2-2*bm*ro(1)*ros(2)+dm*ros(2)
     &                  **2)
            dm0(2, 2) = dm0(2, 2)
     &                  + dss*(am*ro(2)**2+2*bm*ro(2)*ros(1)+dm*ros(1)
     &                  **2)
            dm0(1, 2) = dm0(1, 2)
     &                  + dss*(-am*ro(1)*ro(2)-bm*(ro(1)*ros(1)-ro(2)
     &                  *ros(2))+dm*ros(1)*ros(2))
            dm0(3, 3) = dm0(3, 3)
     &                  + dss*(am*(ro(1)**2+ro(2)**2)+2.*bm*rzeta+dm)
 
c           ...calculate warping
c           SKIPPED 
c
         END DO
         WRITE (*, *) 'ISEG, DS = ', isg, s_segA
      END DO
c
c     integration along webs
c
      IF ( nwebs/=0 ) THEN
         s_webA = 0.
         DO iweb = 1, nwebs
            iseg = webseg(iweb, isct)
            thk = thks(iseg)
            ip1 = isgmt(1, iseg)
            ip2 = isgmt(2, iseg)
            x1 = xpnt(ip1, isct)
            x2 = xpnt(ip2, isct)
            y1 = ypnt(ip1, isct)
            y2 = ypnt(ip2, isct)
            am = ama(iseg)
            dm = dma(iseg)
            bm = bma(iseg)
c
            darea0 = 0.
            daa = 0.
            CALL mdllin(x1, y1, x2, y2, dss, ro, ros, thk, 1)
            dss   =  dss  - corweb
            rzeta = -ro(1)*ros(2) + ro(2)*ros(1)
            rs    =  ro(1)*ros(1) + ro(2)*ros(2)
            s_webA = s_webA + dss
c
c           calculate density terms
c
            am0(1, 1) = am0(1, 1) + dss*am
                                                                                    !bm0(1, 1) = bm0(1, 1) + dss*(-am*ro(1)+bm*ros(2))
                                                                                    !bm0(1, 2) = bm0(1, 2) + dss*(am*ro(2)+bm*ros(1))
                                                                                    !bm0(2, 3) = bm0(2, 3) + dss*(-am*rzeta-bm)
                                                                                    !bm0(3, 3) = bm0(3, 3) + dss*(am*rs)
            bm0(1, 1) = bm0(1, 1) + dss*(-am*ro(1)+bm*ros(2)) !                     !bm0(1, 1) = bm0(1, 1) + dss*(-am*ro(1)+bm*ros(2))
            bm0(1, 2) = bm0(1, 2) + dss*( am*ro(2)+bm*ros(1)) !===> seraf changes   !bm0(1, 2) = bm0(1, 2) + dss*( am*ro(2)+bm*ros(1))
            bm0(2, 3) = bm0(2, 3) + dss*(-am*ro(2)-bm*ros(1)) !                     !bm0(2, 3) = bm0(2, 3) + dss*(-am*rzeta-bm)
            bm0(3, 3) = bm0(3, 3) + dss*( am*ro(1)-bm*ros(2)) !                     !bm0(3, 3) = bm0(3, 3) + dss*( am*rs) 
c
            dm0(1, 1) = dm0(1, 1)
     &                  + dss*(am*ro(1)**2-2*bm*ro(1)*ros(2)+dm*ros(2)
     &                  **2)
            dm0(2, 2) = dm0(2, 2)
     &                  + dss*(am*ro(2)**2+2*bm*ro(2)*ros(1)+dm*ros(1)
     &                  **2)
            dm0(1, 2) = dm0(1, 2)
     &                  + dss*(-am*ro(1)*ro(2)-bm*(ro(1)*ros(1)-ro(2)
     &                  *ros(2))+dm*ros(1)*ros(2))
            dm0(3, 3) = dm0(3, 3)
     &                  + dss*(am*(ro(1)**2+ro(2)**2)+2.*bm*rzeta+dm)
c
c           ...calculate warping
c        SKIPPED
c
         END DO
         WRITE (*, *) 'WEBS, DS = ', s_webA
      END IF
c     ...symmetric terms
      am0(2, 2) = am0(1, 1)
      am0(3, 3) = am0(1, 1)
      dm0(2, 1) = dm0(1, 2)
      dm0(3, 1) = dm0(1, 3)
      dm0(3, 2) = dm0(2, 3)
c
      END SUBROUTINE SECDNS
c  ------------------------------------------------------------------
c
      SUBROUTINE strains(isct,thks, h, zbar, A0, B0, D0, a,b,d,  
     &                   Loads, EpsCur,  xpnt, ypnt, 
     &                   nskseg, lamsg, isgmt, sknseg, 
     &                   nwebs, webseg, nloop, lpseg, 
     &                   nlpsg, rshr  , nlay ,  hl, 
     &                   angl, aahstr , cl   , cc , 
     &                   ezz, ezs, ezj, szz, szs, szj,
     &                   im, strength,fcrit1,fcrit2,corweb, n_p)
      IMPLICIT NONE
      include '/home/seraf/hGAST/pre_hGAST.dir/paramet.for'
      INTEGER I, II, IL, ILAM, ILOOP, IP, IP1, IP2, IS, ISCT, ISEG, ISG,
     &        IWEB, J, K, jj
c
      REAL*8 zbar(nlayer), h(nlayer) 
      REAL*8 hl  (nlayer,mlamin), angl(nlayer,mlamin)
      REAL*8 thks(mseg)
      REAL*8 cl(5, 5, mater),cc(3,3,nlayer) 
      REAL*8 A0(3, 3),       D0(3, 3), B0(3, 3)
      REAL*8 G1, G2
      REAL*8 ro(2), ros(2), ross(2)
      REAL*8 rzeta, rs, tt
      REAL*8 xpnt(maxpoint, msct), ypnt(maxpoint, msct)
      REAL*8 x1, x2, y1, y2, dss, thk
      REAL*8 gint(2), pint(2), f(2), g(2), qzero(2), sta66, rshr(2)
      REAL*8 scstf(6, 6), seccmp(6, 6), det1(2), wk(6)
      REAL*8 aahstr(nloop)
      REAL*8 Loads(6), EpsCur(6)
      REAL*8 AAH      
      REAL*8 corweb
c
      INTEGER NL, NL2
      INTEGER nskseg, isgmt(2, mseg), lamsg(mseg)
      INTEGER kpvt(6), inert(3), info
      INTEGER sknseg(mseg, msct), webseg(mwebs, msct), nwebs
      INTEGER lpseg(mseg, mloops), nloop, nlpsg(mloops)
      INTEGER im(nlayer, mlamin), nlay(mlamin)
 
      REAL*8 ez0, ezx0, ezy0, kapx, kapy, kapz
      REAL*8 ezz(maxpoint,nlayer), ezs (maxpoint,nlayer),
     $       ezj(maxpoint,nlayer)
      REAL*8 szz(maxpoint,nlayer), szs (maxpoint,nlayer),
     $       szj(maxpoint,nlayer)
      REAL*8 twu(maxpoint,nlayer), twuR(maxpoint, nlayer)
      REAL*8 strength(mater,mater)
      
      REAL*8 a(3,3,mseg), b(3,3,mseg), d(3,3,mseg)
      REAL*8 offset
      REAL*8 stres(6)
      REAL*8 stren(9)
      REAL*8 fcrit1,fcrit2    
       
      REAL*8 qshear(mseg,maxpoint)
      REAL*8 Toshear,qshear0
      REAL*8 szzz(maxpoint,nlayer),ezzz(maxpoint,nlayer)
      REAL*8 szss(maxpoint,nlayer),ezss(maxpoint,nlayer)
      REAL*8 str_twu(1:7), str_nor(1:7), str_per(1:7)
      CHARACTER*300::n_p
      CHARACTER*400::path
c
      fcrit1=-1.d8
      fcrit2= 1.d8
c 
      CALL SHEARFLO(isct, a , b, A0, B0, D0, xpnt, ypnt, 
     &     nskseg,lamsg,isgmt,sknseg,nwebs, webseg, 
     &     nloop, lpseg, nlpsg,corweb,Loads(2),Loads(3),Toshear,qshear)
      write(*,*)' Torque of Shear is ',Toshear
  
c    
      path=trim(adjustl(n_p))//"/other.dir/stresses.dat"
      OPEN (18, FILE=path, STATUS='unknown')
      path=trim(adjustl(n_p))//"/other.dir/stresses_tswu.dat"
      OPEN(UNIT=44, FILE=path)                                !===================================> seraf changes
      path=trim(adjustl(n_p))//"/other.dir/stresses_norm.dat"
      OPEN(UNIT=45, FILE=path)                                !===================================> seraf changes
      path=trim(adjustl(n_p))//"/other.dir/stresses_peri.dat"
      OPEN(UNIT=46, FILE=path)                                !===================================> seraf changes
c     ...invert stiffness to calculate section compliance
      DO i = 1, 3
         DO j = 1, 3
            seccmp(i  , j  ) = a0(i, j)
            seccmp(i  , j+3) = b0(i, j)
            seccmp(i+3, j  ) = b0(j, i)
            seccmp(i+3, j+3) = d0(i, j)
         END DO
      END DO
      CALL dsifa(seccmp, 6, 6, kpvt, info)
      CALL dsidi(seccmp, 6, 6, kpvt, det1, inert, wk, 001)
c
c     expand lower triangle
      DO i = 2, 6
         DO j = 1, i - 1
            seccmp(i, j) = seccmp(j, i)
         END DO
      END DO
c
      DO i = 1, 6
         EpsCur(i) = 0.
         DO j = 1, 6
            EpsCur(i) = EpsCur(i) + seccmp(i, j)*Loads(j)
         END DO
      END DO
c
      DO i = 1, 6
         stres(i) = 0.
      END DO
c
      ez0  = EpsCur(1)
      ezx0 = EpsCur(2)
      ezy0 = EpsCur(3)
      kapy = EpsCur(4)
      kapx = EpsCur(5)
      kapz = EpsCur(6)
c
      DO isg = 1, nskseg
         ilam = lamsg(isg)
         nl = nlay(ilam)
         thk = thks(isg)
c        ... calculate zbar
         zbar(1) = (hl(1,ilam)-thk)/2.
         DO k = 2, nl
            zbar(k) = zbar(k-1) + (hl(k-1,ilam)+hl(k,ilam))/2.
         END DO
c        ... determine the loop where the segment belongs for aah
         DO il = 1, nloop
            DO is = 1, nlpsg(il)
               IF ( isg==abs(lpseg(is,il)) ) iloop = il
            END DO
         END DO
c
         DO k = 1, nl
         CALL offply (im(k,ilam), k, angl(k,ilam),cl,cc)
c         CALL offply1(im(k,ilam), k, angl(k,ilam),cl,cc)
         enddo
c         
         ip1 = isgmt(1, isg)
         ip2 = isgmt(2, isg)
C        WRITE(44,*) (ip2-ip1) !=========================================================> seraf changes 
         DO ip = ip1, ip2 - 1
            x1 = xpnt(ip  , isct)
            x2 = xpnt(ip+1, isct)
            y1 = ypnt(ip  , isct)
            y2 = ypnt(ip+1, isct)
            CALL mdllin(x1, y1, x2, y2, dss, ro, ros, thk, 1)
            rzeta = -ro(1)*ros(2) + ro(2)*ros(1)
            rs    =  ro(1)*ros(1) + ro(2)*ros(2)
            aah = aahstr(iloop)/a(3,3,isg)-2.d0*b(3,3,isg)/a(3,3,isg)

            nl2=nl                              !
            DO k=1,nl                           !
               if (im(k,ilam).eq.4) nl2=nl2-1   !=======================================> seraf changes
            END DO                              !
C           WRITE(44,*) nl2                     !

            str_twu=0.0
            str_nor=0.0
            str_per=0.0
            DO k = 1, nl
               offset = angl(k, ilam)
               ezz (ip, k) =  ez0 - kapy*(ro(1)-ros(2)*zbar(k))
     &                            + kapx*(ro(2)+ros(1)*zbar(k))
               ezs (ip, k) =  ezx0*ros(1) + ezy0*ros(2)
     &                            - kapz*(aah+2.*zbar(k))
               ezss(ip, k) =  (qshear(isg,ip)-2.*kapz*b(3,3,isg))
     &                                               /a(3,3,isg)
     &                            - kapz*(aah+2.*zbar(k))    
               ezj (ip, k) = -ezx0*ros(2) + ezy0*ros(1)
               szz (ip, k) =  ezz(ip, k)*cc(1, 1, k) +
     &                        ezs(ip, k)*cc(1, 3, k)
               szzz(ip, k) =  ezz(ip, k)*cc(1, 1, k) +
     &                       ezss(ip, k)*cc(1, 3, k)
               szs (ip, k) =  ezz(ip, k)*cc(3, 1, k) +
     &                        ezs(ip, k)*cc(3, 3, k)
               szss(ip, k) =  ezz(ip, k)*cc(3, 1, k) +
     &                       ezss(ip, k)*cc(3, 3, k)
               szj (ip, k) =  ezj(ip, k)*cc(2, 2, k)
               if (abs(offset)>0.001) szss(ip,k)=-szss(ip,k)
               stres(1)   =  szzz(ip, k)
               stres(5)   =  szj (ip, k)
               stres(6)   =  szss(ip, k)
               DO ii = 1, 9
               stren(ii) = strength(ii, im(k,ilam))
               END DO
               CALL FailCrit(stres, stren, offset, twu(ip,k),twuR(ip,k))
               if(im(k,ilam).ne.4.and.twu (ip,k).ne.0.d0) 
     &                             fcrit1=max(fcrit1,twu (ip,k))
               if(im(k,ilam).ne.4.and.twuR(ip,k).ne.0.d0) 
     &                             fcrit2=min(fcrit2,twuR(ip,k))
               if (im(k,ilam).ne.4) then
               WRITE (18,99001) isg, ip, k , im(k,ilam), ro(1), ro(2),
     &                           ezz(ip, k), ezs(ip, k), ezj(ip, k),
     &                           szz(ip, k), szs(ip, k), szj(ip, k),
     &                           twu(ip, k),twuR(ip, k),
     &                          szzz(ip, k),szss(ip, k)
               end if
               if (k<=(nl2/2.0)) then
                  if (im(k,ilam).eq.1) then
                     str_twu(1)=twuR(ip,k)
                     str_nor(1)=szzz(ip,k)
                     str_per(1)=szss(ip,k)
                  end if
                  if (im(k,ilam).eq.2) then
                     str_twu(2)=twuR(ip,k)
                     str_nor(2)=szzz(ip,k)
                     str_per(2)=szss(ip,k)
                  end if
                  if (im(k,ilam).eq.3) then
                     str_twu(3)=twuR(ip,k)
                     str_nor(3)=szzz(ip,k)
                     str_per(3)=szss(ip,k)
                  end if
               else
                  if (im(k,ilam).eq.3) then
                     str_twu(5)=twuR(ip,k)
                     str_nor(5)=szzz(ip,k)
                     str_per(5)=szss(ip,k)
                  end if
                  if (im(k,ilam).eq.2) then
                     str_twu(6)=twuR(ip,k)
                     str_nor(6)=szzz(ip,k)
                     str_per(6)=szss(ip,k)
                  end if
                  if (im(k,ilam).eq.1) then
                     str_twu(7)=twuR(ip,k)
                     str_nor(7)=szzz(ip,k)
                     str_per(7)=szss(ip,k)
                  end if
               end if
            END DO
            WRITE(44,97579) (str_twu(jj), jj=1,7) !===================================> seraf changes
            WRITE(45,97579) (str_nor(jj), jj=1,7) !===================================> seraf changes
            WRITE(46,97579) (str_per(jj), jj=1,7) !===================================> seraf changes
 
         END DO
      END DO
c
c     along webs
c
      IF ( nwebs/=0 ) THEN
         DO iweb = 1, nwebs
            iseg = webseg(iweb, isct)
            ilam = lamsg(iseg)
            nl = nlay(ilam)
            thk = thks(iseg)
c           ... calculate zbar
            zbar(1) = (hl(1,ilam)-thk)/2.
            DO k = 2, nl
               zbar(k) = zbar(k-1) + (hl(k-1,ilam)+hl(k,ilam))/2.
            END DO
c
          DO k = 1, nl
          CALL offply (im(k,ilam), k, angl(k,ilam),cl,cc)
c          CALL offply1(im(k,ilam), k, angl(k,ilam),cl,cc)
          enddo
c 
            ip1 = isgmt(1, iseg)
            ip2 = isgmt(2, iseg)
C           WRITE(44,*) "1" !=========================================================> seraf changes 
            x1 = xpnt(ip1, isct)
            x2 = xpnt(ip2, isct)
            y1 = ypnt(ip1, isct)
            y2 = ypnt(ip2, isct)
            CALL mdllin(x1, y1, x2, y2, dss, ro, ros, thk, 1)
            dss   =  dss  - corweb
            rzeta = -ro(1)*ros(2) + ro(2)*ros(1)
            rs    =  ro(1)*ros(1) + ro(2)*ros(2)
            IF ( iweb==1 ) aah =-(aahstr(1)-aahstr(2))/a(3,3,iseg)
     &                               -2.d0*b(3,3,iseg)/a(3,3,iseg)
            IF ( iweb==2 ) aah = (aahstr(2)-aahstr(3))/a(3,3,iseg)
     &                               -2.d0*b(3,3,iseg)/a(3,3,iseg)

            nl2=nl                              !
            DO k=1,nl                           !
               if (im(k,ilam).eq.4) nl2=nl2-1   !=======================================> seraf changes
            END DO                              !
C           WRITE(44,*) nl2                     !

            str_twu=0.0
            DO k = 1, nl
               offset = angl(k, ilam)
               ezz (ip1, k) =  ez0 - kapy*(ro(1)-ros(2)*zbar(k))
     &                             + kapx*(ro(2)+ros(1)*zbar(k))
               ezs (ip1, k) =  ezx0*ros(1) + ezy0*ros(2)
     &                             - kapz*(aah+2.*zbar(k))
               ezss(ip1, k) =  (qshear(iseg,ip1)-2.*kapz*b(3,3,iseg))
     &                                                  /a(3,3,iseg)
     &                             - kapz*(aah+2.*zbar(k))   
               ezj (ip1, k) = -ezx0*ros(2) + ezy0*ros(1)
               szz (ip1, k) =  ezz(ip1, k)*cc(1, 1, k) +
     &                         ezs(ip1, k)*cc(1, 3, k)
               szzz(ip1, k) =  ezz(ip1, k)*cc(1, 1, k) +
     &                        ezss(ip1, k)*cc(1, 3, k)
               szs (ip1, k) =  ezz(ip1, k)*cc(3, 1, k) + 
     &                         ezs(ip1, k)*cc(3, 3, k)
               szss(ip1, k) =  ezz(ip1, k)*cc(3, 1, k) + 
     &                        ezss(ip1, k)*cc(3, 3, k)
               szj (ip1, k) =  ezj(ip1, k)*cc(2, 2, k)
               stres(1)    =  szzz(ip1, k)
               stres(5)    =  szj (ip1, k)
               stres(6)    =  szss(ip1, k)
               DO ii = 1, 9
               stren(ii) = strength(ii, im(k,ilam))
               END DO
               CALL FailCrit(stres,stren,offset,twu(ip1,k),twuR(ip1,k))
               if(im(k,ilam).ne.4.and.twu (ip1,k).ne.0.d0) 
     &                             fcrit1=max(fcrit1,twu (ip1,k))
               if(im(k,ilam).ne.4.and.twuR(ip1,k).ne.0.d0) 
     &                             fcrit2=min(fcrit2,twuR(ip1,k))
               if (im(k,ilam).ne.4) then
               WRITE (18,99001) iseg, ip1, k , im(k,ilam), ro(1), ro(2),
     &                            ezz(ip1, k), ezs(ip1, k), ezj(ip1, k),
     &                            szz(ip1, k), szs(ip1, k), szj(ip1, k),
     &                            twu(ip1, k),twuR(ip1, k),
     &                           szzz(ip1, k),szss(ip1, k)
               end if
               if (k<=(nl2/2.0)) then
                  if (im(k,ilam).eq.1) then
                     str_twu(1)=twuR(ip1,k)
                     str_nor(1)=szzz(ip1,k)
                     str_per(1)=szss(ip1,k)
                  end if
                  if (im(k,ilam).eq.2) then
                     str_twu(2)=twuR(ip1,k)
                     str_nor(2)=szzz(ip1,k)
                     str_per(2)=szss(ip1,k)
                  end if
                  if (im(k,ilam).eq.3) then
                     str_twu(3)=twuR(ip1,k)
                     str_nor(3)=szzz(ip1,k)
                     str_per(3)=szss(ip1,k)
                  end if
               else
                  if (im(k,ilam).eq.3) then
                     str_twu(5)=twuR(ip1,k)
                     str_nor(5)=szzz(ip1,k)
                     str_per(5)=szss(ip1,k)
                  end if
                  if (im(k,ilam).eq.2) then
                     str_twu(6)=twuR(ip1,k)
                     str_nor(6)=szzz(ip1,k)
                     str_per(6)=szss(ip1,k)
                  end if
                  if (im(k,ilam).eq.1) then
                     str_twu(7)=twuR(ip1,k)
                     str_nor(7)=szzz(ip1,k)
                     str_per(7)=szss(ip1,k)
                  end if
               end if
            END DO
            write(44,97579) (str_twu(jj), jj=1,7)
            write(45,97579) (str_nor(jj), jj=1,7)
            write(46,97579) (str_per(jj), jj=1,7)
         END DO
      END IF
      CLOSE (18)
      CLOSE (44)
      CLOSE (45)
      CLOSE (46)
c
99001 FORMAT (4I4, 12(1x,e12.4))
97579 FORMAT(7(F30.7, 3X))
c
      END SUBROUTINE STRAINS
 
      SUBROUTINE FailCrit(stress, strength, offset, twu, twuR)
      IMPLICIT NONE
      INTEGER I, J
c     subroutine for the determination of failure based on Tsai-Wu
c     written for 3D stress case
c     thermal stresses are not taken into account
 
      REAL*8 stress(6)
      REAL*8 strength(9)
      REAL*8 offset
c     offset is the offset angle of the layer
c     stress(1)=sigma(n), stress(2)=sigma(s), stress(3)=sigma(j)=0.
c     stress(4)=sigma(sj)=0., stress(5)=sigma(nj), stress(6)=sigma(ns)
c     strength:

c     strength(1)= Tension strength in direction 1 (along the fibres)   :  XT 
c     strength(2)= Compressive strength in direction 1 (along the fibres): XC 
c     strength(3)= Tension strength in direction 2 (Transverse to the fibres): YT
c     strength(4)= Compressive strength in direction 2 (Transverse to the fibres): YC 
c     strength(5)= Tension strength in direction 3 (through the thickness): ZT can be neglected 
c     strength(6)= Compressive strength in direction 3 (Through the thickness): ZC ' can be neglected 
c     strength(7)= Transverse Shear strength through the thickness: S23 
c     strength(8)= Transverse Shear strength: S13 
c     strength(9)= Shear strengh in plane: S12
c
      REAL*8 twu, twuR                !this will be the output values
c     twu the output value of tsai-wu failure criterion. Failure when
c     twu>0. twuR is the strength Ratio for the criterion. Failure when
c     twuR>1
      REAL*8 p, thita, c, s, T(6, 6)                   ! local variables
      REAL*8 Son(6), Fij(6, 6), Fi(6), ACON, BCON, DCON
                                                       ! local variables
      REAL*8 R(2)
      REAL*8 a, b, sxt, sxc, syt, syc, szt, szc, s12, s13, s23
      INTEGER caseis                        ! criterion selection
 
c     caseis shows the interaction terms Fij (when i.ne.j) now set
c     caseis=1 for regular tsai-hahn
      Caseis = 1
 
c     TRANSFORM ANGLE/ DEFINE COSINUS AND SINUS
      p = 4.d0*datan(1.d0)
      thita = p*offset/180.
      c = cos(thita)
      s = sin(thita)
 
c     ROTATION of APPLIED STRESSES on LAYER MAIN AXIS
c     transformation matrix 6x6
      T = 0.
      T(1, 1) = c**2.
      T(1, 2) = s**2.
      T(1, 6) = 2.*c*s
      T(2, 1) = s**2.
      T(2, 2) = c**2.
      T(2, 6) = -2.*c*s
      T(3, 3) = 1.
      T(4, 4) = c
      T(4, 5) = -s
      T(5, 4) = s
      T(5, 5) = c
      T(6, 1) = -c*s
      T(6, 2) = c*s
      T(6, 6) = c**2. - s**2.
c     the above matrix can be reduced by deleting row 3 to a 5x5 matrix
 
c     (stress should be reduced also)
      Son = 0.
      DO i = 1, 6
         DO j = 1, 6
            Son(i) = Son(i) + T(i, j)*stress(j)
         END DO
      END DO
 
c     Define failure tensor coefficients
      Fij = 0.
      Fij(1, 1) = 1./Strength(1)/Strength(2)
      Fij(2, 2) = 1./Strength(3)/Strength(4)
      Fij(3, 3) = 1./Strength(5)/Strength(6)
      Fij(4, 4) = 1./Strength(7)**2.
      Fij(5, 5) = 1./Strength(8)**2.
      Fij(6, 6) = 1./Strength(9)**2.
 
      Fi(1) = 1./Strength(1) - 1./Strength(2)
      Fi(2) = 1./Strength(3) - 1./Strength(4)
      Fi(3) = 1./Strength(5) - 1./Strength(6)
 
      IF ( CASEIS==1 ) THEN               !Tsai Hahn regular
         Fij(1, 2) = -.5*dSQRT(Fij(1,1)*Fij(2,2))
         Fij(2, 1) = Fij(1, 2)
         Fij(1, 3) = -.5*dSQRT(Fij(1,1)*Fij(3,3))
         Fij(3, 1) = Fij(1, 3)
         Fij(2, 3) = -.5*dSQRT(Fij(2,2)*Fij(3,3))
         Fij(3, 2) = Fij(2, 3)
      ELSE IF ( CASEIS==2 ) THEN          !paraboloid surface
         Fij(1, 2) = .5*(Fij(3,3)-Fij(1,1)-Fij(2,2))
         Fij(2, 1) = Fij(1, 2)
         Fij(1, 3) = .5*(Fij(2,2)-Fij(1,1)-Fij(3,3))
         Fij(3, 1) = Fij(1, 3)
         Fij(2, 3) = .5*(Fij(1,1)-Fij(2,2)-Fij(3,3))
         Fij(3, 2) = Fij(2, 3)
c        for inplane only reduced to
c        Fij(1,2)=-.5*Fij(1,1)
c        Fij(2,1)=Fij(1,2)
      ELSE IF ( CASEIS==3 ) THEN
         Fij(1, 2) = -Fij(1, 1)*Fij(2, 2)/(Fij(1,1)+Fij(2,2))
         Fij(2, 1) = Fij(1, 2)
         Fij(1, 3) = -Fij(1, 1)*Fij(3, 3)/(Fij(1,1)+Fij(3,3))
         Fij(3, 1) = Fij(1, 3)
         Fij(2, 3) = -Fij(2, 2)*Fij(3, 3)/(Fij(2,2)+Fij(3,3))
         Fij(3, 2) = Fij(2, 3)
      ELSE IF ( CASEIS==4 ) THEN         !Neglect coupling terms
         Fij(1, 2) = 0.
         Fij(2, 1) = Fij(1, 2)
         Fij(1, 3) = 0.
         Fij(3, 1) = Fij(1, 3)
         Fij(2, 3) = 0.
         Fij(3, 2) = Fij(2, 3)
      END IF
 
c     CALCULATING STRENGTH RATIOS
      ACON = 0.
      DO I = 1, 6
         DO J = 1, 6
            ACON = ACON + Fij(I, J)*Son(I)*Son(J)
         END DO
      END DO
c
      BCON = 0.
      DO I = 1, 3
         BCON = BCON + Fi(I)*Son(I)
      END DO
 
      DCON = -1.D0
      R(1) = (-BCON+dSQRT(BCON**2.-4.*ACON*DCON))/2.D0/ACON
      R(2) = (-BCON-dSQRT(BCON**2.-4.*ACON*DCON))/2.D0/ACON
 
c ##############     twuR = max(R(1), R(2)) ! -----> *** testing ***
      sxt=strength(1)  ; syt=strength(3)  ; szt=strength(5)
      sxc=-strength(2) ; syc=-strength(4) ; szc=-strength(6)
      s23=strength(7)  ; s13=strength(8)  ; s12=strength(9)
      if (abs(offset)>0.001) then
         Son(6)=-0.50*Son(6)
      else
         Son(6)=-Son(6)
      end if
c     Son(1)=(-0.00005*(thita**3)+0.0017*(thita**2)
c    $        -0.0268*thita+0.9988)*Son(1)
c     Son(6)=(-0.0018*(thita**2)+0.0173*thita+2.0593)*Son(6)
c     end if

      twuR=-(Son(1)**2)/(sxt*sxc)+Son(1)*(1/sxt+1/sxc)+(Son(6)/s23)**2

      twu = ACON + BCON + DCON
 
      END SUBROUTINE FAILCRIT
      
      
      SUBROUTINE TRLFSM(A0, B0, D0, An, Bn, Dn, xc, yc)
      IMPLICIT NONE
      include '/home/seraf/hGAST/pre_hGAST.dir/paramet.for'
      INTEGER i,j
      REAL*8 xc,yc
      REAL*8 A0(3,3)     , B0(3, 3)    , D0(3, 3)
      REAL*8 An(3,3)     , Bn(3, 3)    , Dn(3, 3)
c
      do i    =    1,3
      do j    =    1,3
      An(i,j) = A0(i,j) 
      enddo
      enddo
    
      Bn(1,1) = B0(1,1) -xc*A0(1,1)
      Bn(1,2) = B0(1,2) +yc*A0(1,1)
      Bn(2,1) = B0(2,1) -xc*A0(1,2)
      Bn(2,2) = B0(2,2) +yc*A0(1,2)
      Bn(3,1) = B0(3,1) -xc*A0(1,3)
      Bn(3,2) = B0(3,2) +yc*A0(1,3)
      
      Bn(1,3) = B0(1,3) +xc*A0(1,3)-yc*A0(2,1)
      Bn(2,3) = B0(2,3) +xc*A0(2,3)-yc*A0(2,2)
      Bn(3,3) = B0(3,3) +xc*A0(3,3)-yc*A0(2,3)
       
      Dn(1,1) = D0(1,1) -2.*xc*B0(1,1)           +xc*xc*A0(1,1)
      Dn(2,2) = D0(2,2) +2.*yc*B0(1,2)           +yc*yc*A0(1,1)
      Dn(1,2) = D0(1,2) -   xc*B0(1,2)+yc*B0(1,1)-xc*yc*A0(1,1)
      Dn(1,3) = D0(1,3) -   xc*B0(1,3)+xc*B0(3,1)-xc*xc*A0(1,3)
     &                  -   yc*B0(2,1)           +xc*yc*A0(1,2)
      Dn(2,3) = D0(2,3) +   yc*B0(1,3)+xc*B0(3,2)-yc*yc*A0(1,2)
     &                  -   yc*B0(2,2)           +xc*yc*A0(1,3)
      Dn(3,3) = D0(3,3) +xc*xc*A0(3,3)-2.*xc*yc*A0(2,3)+yc*yc*A0(2,2)
     &                  +2.*xc*B0(3,3)                 -2.*yc*B0(2,3)
      Dn(2,1) = Dn(1,2)
      Dn(3,1) = Dn(1,3)
      Dn(3,2) = Dn(2,3)
      
      END SUBROUTINE TRLFSM
      
            
      SUBROUTINE ROTFSM(A0, B0, D0, An, Bn, Dn, THETA)
      IMPLICIT NONE
      include '/home/seraf/hGAST/pre_hGAST.dir/paramet.for'
      INTEGER i,j
      REAL*8 theta
      REAL*8 A0(3,3)     , B0(3, 3)    , D0(3, 3)
      REAL*8 An(3,3)     , Bn(3, 3)    , Dn(3, 3)
      REAL*8 c,s,c2,sc,cc,ss,pi
c
      PI = 4.D0*datan(1.D0)
      c  = COS(THETA*PI/180.0)
      s  = SIN(THETA*PI/180.0)
      cc = c*c
      ss = s*s
      sc = s*c
      c2 = cc-ss
      
      An(1,1) =    A0(1,1) 
      An(1,2) =  c*A0(1,2)             + s*A0(1,3)
      An(1,3) = -s*A0(1,2)             + c*A0(1,3) 
      An(2,2) = cc*A0(2,2)+2*sc*A0(2,3)+ss*A0(3,3)
      An(2,3) = sc*  (A0(3,3)-A0(2,2)) +c2*A0(2,3) 
      An(3,3) = ss*A0(2,2)-2*sc*A0(2,3)+cc*A0(3,3)
      An(2,1) = An(1,2)
      An(3,1) = An(1,3)
      An(3,2) = An(2,3)
    
      Bn(1,1) =  c*B0(1,1)             - s*B0(1,2)
      Bn(1,2) =  s*B0(1,1)             + c*B0(1,2)
      Bn(1,3) = B0(1,3) 
      Bn(2,1) = cc*B0(2,1)-sc*(B0(2,2)-B0(3,1))-ss*B0(3,2) 
      Bn(2,2) = cc*B0(2,2)+sc*(B0(2,1)+B0(3,2))+ss*B0(3,1)
      Bn(2,3) =  c*B0(2,3)             + s*B0(3,3)
      Bn(3,1) = cc*B0(3,1)-sc*(B0(2,1)+B0(3,2))+ss*B0(2,2)
      Bn(3,2) = cc*B0(3,2)-sc*(B0(2,2)-B0(3,1))-ss*B0(2,1)
      Bn(3,3) = -s*B0(2,3)             + c*B0(3,3)
       
      Dn(1,1) = cc*D0(1,1)-2*sc*D0(1,2)+ss*D0(2,2) 
      Dn(1,2) = sc*  (D0(1,1)-D0(2,2)) +c2*D0(1,2) 
      Dn(1,3) =  c*D0(1,3)             - s*D0(2,3)
      Dn(2,2) = ss*D0(1,1)+2*sc*D0(1,2)+cc*D0(2,2)  
      Dn(2,3) =  s*D0(1,3)             + c*D0(2,3)
      Dn(3,3) = D0(3,3) 
      Dn(2,1) = Dn(1,2)
      Dn(3,1) = Dn(1,3)
      Dn(3,2) = Dn(2,3)
      
      END SUBROUTINE ROTFSM
     
c ------------------------------------------------
c     Scope:  calculates shear center of section
c             see notes
c
c            a0, b0, d0:       section extension, coupling, flexure matrices -- see notes
c                         wr. to the y, z axes at center of cross section
c            a, b, d:           extension, coupling, flexure matrices of face (skin) laminates
c                         about local centeline (of skin)
c                          at curvilinear cordinates (x, s, zeta) -- see notes....
c            yshr, zshr:  position of shear center w.r. to coordinate center
c
c
      SUBROUTINE SHEARFLO(isct, a , b, A0, B0, D0, xpnt, ypnt, 
     &                   nskseg,lamsg,isgmt,sknseg,nwebs, webseg, 
     &                   nloop, lpseg, nlpsg, corweb,QX0,QY0,To0,q2,n_p)
      IMPLICIT NONE
      include '/home/seraf/hGAST/pre_hGAST.dir/paramet.for'
      INTEGER I, IP, IP1, IP2, ISCT, ISEG, ISG, IWEB, J, K, 
     &        MCOU,  NCOU,JP1, JP2 , IL , IS , LPS , ILOOP
c
      REAL*8 A0(3,3)     , B0(3, 3)    , D0(3, 3)
      REAL*8 A (3,3,mseg), B (3,3,mseg)
      REAL*8 ro(2), ros(2), ross(2)
      REAL*8 rzeta, rs
      REAL*8 xpnt(maxpoint, msct), ypnt(maxpoint, msct)
      REAL*8 x1, x2, y1, y2, dss, thk
      REAL*8 scstf(6, 6), seccmp(6, 6), det1(2), wk(6)
      REAL*8 px1,px2,px3
      REAL*8 py1,py2,py3
      REAL*8 x0,y0
      REAL*8 deter      
      REAL*8 corweb
c
      INTEGER nskseg, isgmt(2, mseg), lamsg(mseg)
      INTEGER kpvt(6), inert(3), info
      INTEGER sknseg(mseg, msct), webseg(mwebs, msct), nwebs
      INTEGER lpseg (mseg, mloops), nloop, nlpsg(mloops)
c  
      REAL*8 QX0,QY0,To0
      REAL*8 qmid   
      REAL*8 fintx(mseg,maxpoint),finty(mseg,maxpoint)
      REAL*8 q0   (mseg,maxpoint),q1   (mseg,maxpoint)
      REAL*8 q2   (mseg,maxpoint)
      REAL*8 qsup (mseg), qsuw(mwebs)
      REAL*8 qalf (mloops,mloops), qrhs  (mloops),qres(mloops)
      REAL*8 aa   (mseg)  ,area  (mseg)
      REAL*8 area0(mloops),aahstr(mloops),aaa(4, 4),aaloop(4)
      INTEGER kqvt(4)
      INTEGER NODE(mwebs),ISGL(mwebs)
      INTEGER ISGLMIN, ISGLMAX, IWEBMIN, IWEBMAX
      CHARACTER*300 n_p
c
c     ..initialize parameters
      DO isg = 1,mseg
      DO i   = 1,maxpoint
      fintx(isg,i)  = 0.
      finty(isg,i)  = 0.
      q0   (isg,i)  = 0.
      q1   (isg,i)  = 0.
      q2   (isg,i)  = 0.
      END DO
      END DO
      
      DO     isg    = 1, mseg
      aa    (isg)   = 0.
      area  (isg)   = 0.
      END DO
      DO     iloop  = 1, nloop
      area0 (iloop) = 0.
      aaloop(iloop) = 0.
      aahstr(iloop) = 0.
      END DO
c      
c     ...invert stiffness to calculate section compliance

      open (1,file=trim(n_p)//"/other.dir/seccmp.dat",status='unknown')
      DO i = 1, 3
         DO j = 1, 3
            seccmp(i  , j  ) = A0(i, j)
            seccmp(i  , j+3) = B0(i, j)
            seccmp(i+3, j  ) = B0(j, i)
            seccmp(i+3, j+3) = D0(i, j)
         END DO
      END DO
      do i=1,6
      write(1,771)(seccmp(i,j),j=1,6)
      enddo
      write(1,*)
      
      CALL dsifa(seccmp, 6, 6, kpvt, info)
      CALL dsidi(seccmp, 6, 6, kpvt, det1, inert, wk, 001)
c
c     expand lower triangle
      DO i = 2, 6
         DO j = 1, i - 1
            seccmp(i, j) = seccmp(j, i)
         END DO
      END DO
c
      do i=1,6
      write(1,771)(seccmp(i,j),j=1,6)
      enddo
      write(1,*)
  771 format(6(1x,e15.6))
c  
c     WEBS
c
      IF ( nwebs/=0 ) THEN
         DO iweb = 1, nwebs
            isg = webseg(iweb, isct)
            ip1 = isgmt (1, isg)
            ip2 = isgmt (2, isg)
      fintx(isg,ip1) = 0. 
      finty(isg,ip1) = 0. 
      q0   (isg,ip1) = 0.
            x1  = xpnt(ip1,isct)
            x2  = xpnt(ip2,isct)
            y1  = ypnt(ip1,isct)
            y2  = ypnt(ip2,isct)
            CALL mdllin(x1, y1, x2, y2, dss, ro, ros, thk, 1)
            dss   =  dss  - corweb
      px1=a(1,1,isg)*(seccmp(1,4)-ro (1)*seccmp(4,4)+ro (2)*seccmp(5,4))
      px2=b(1,1,isg)*(            ros(1)*seccmp(5,4)+ros(2)*seccmp(4,4))
      px3=a(1,3,isg)*(            ros(1)*seccmp(2,4)+ros(2)*seccmp(3,4))
      py1=a(1,1,isg)*(seccmp(1,5)-ro (1)*seccmp(4,5)+ro (2)*seccmp(5,5))
      py2=b(1,1,isg)*(            ros(1)*seccmp(5,5)+ros(2)*seccmp(4,5))
      py3=a(1,3,isg)*(            ros(1)*seccmp(2,5)+ros(2)*seccmp(3,5))
      fintx(isg,ip2) = fintx(isg,ip1)+(px1+px2+px3)*dss 
      finty(isg,ip2) = finty(isg,ip1)+(py1+py2+py3)*dss  
      q0   (isg,ip2) = QX0*fintx(isg,ip2)-QY0*finty(isg,ip2) 
      q1   (isg,ip1) =    0.5*q0(isg,ip2)
      qsuw (iweb)    =        q1(isg,ip1)    
         END DO
      END IF
c
c     SKIN
c
      DO isg = 1, nskseg
         ip1 = isgmt(1, isg)
         ip2 = isgmt(2, isg)
         DO ip = ip1, ip2 - 1
            x1 = xpnt(ip  , isct)
            x2 = xpnt(ip+1, isct)
            y1 = ypnt(ip  , isct)
            y2 = ypnt(ip+1, isct)
            CALL mdllin(x1, y1, x2, y2, dss, ro, ros, thk, 1)
      px1=a(1,1,isg)*(seccmp(1,4)-ro (1)*seccmp(4,4)+ro (2)*seccmp(5,4))
      px2=b(1,1,isg)*(            ros(1)*seccmp(5,4)+ros(2)*seccmp(4,4))
      px3=a(1,3,isg)*(            ros(1)*seccmp(2,4)+ros(2)*seccmp(3,4))
      py1=a(1,1,isg)*(seccmp(1,5)-ro (1)*seccmp(4,5)+ro (2)*seccmp(5,5))
      py2=b(1,1,isg)*(            ros(1)*seccmp(5,5)+ros(2)*seccmp(4,5))
      py3=a(1,3,isg)*(            ros(1)*seccmp(2,5)+ros(2)*seccmp(3,5))
      fintx(isg,ip+1) = fintx(isg,ip)+(px1+px2+px3)*dss 
      finty(isg,ip+1) = finty(isg,ip)+(py1+py2+py3)*dss  
      q0   (isg,ip+1) = QX0*fintx(isg,ip+1)-QY0*finty(isg,ip+1)
      q1   (isg,ip  ) = 0.5*(  q0(isg,ip  )+       q0(isg,ip+1))
         END DO
      END DO
c
      DO isg = 1, nskseg
         ip1 = isgmt(1, isg)
         ip2 = isgmt(2, isg)
         IF(isg.eq.1) then
         ELSE
               jp1 = isgmt(1, isg-1)
               jp2 = isgmt(2, isg-1)
               DO ip = ip1, ip2-1 
               q1(isg,ip) = q0(isg,ip) + q1(isg-1,jp2-1)
               END DO
         ENDIF
      END DO
c
c     DETECT ISGL(i),IWEBL(i) 
c   
      IF ( nwebs/=0 ) THEN
c       
        DO   iweb  =  1,nwebs
        node(iweb) =  isgmt(2,webseg(iweb,isct)) 
          DO    isg  =  1,nskseg
          IF( isgmt(1,isg)==node(iweb)) isgl(iweb)=isg
          END DO
          write(*,*)iweb,node(iweb),isgl(iweb)
        END DO 
c      
        IF(nwebs==2) THEN
c        
          IF(isgl(1).LT.isgl(2))THEN
          isglmin=isgl(1)
          iwebmin=1
          isglmax=isgl(2)
          iwebmax=2
          write(*,*)isglmin,isglmax,iwebmin,iwebmax
          ELSE
          isglmin=isgl(2)
          iwebmin=2
          isglmax=isgl(1)
          iwebmax=1
          write(*,*)isglmin,isglmax,iwebmin,iwebmax
          END IF
c      
          DO isg = 1, nskseg
          IF(isg.lt.isglmin) qsup(isg) = 0.
          IF(isg.ge.isglmin) qsup(isg) = qsuw(iwebmin)
          IF(isg.ge.isglmax) qsup(isg) = qsuw(iwebmin)+qsuw(iwebmax)
          ip1 = isgmt(1, isg)
          ip2 = isgmt(2, isg)
                DO ip = ip1, ip2-1 
                q1(isg,ip) = q1(isg,ip)+qsup(isg)
                END DO
          END DO
        END IF  !nwebs=2
      END IF ! nwebs/=0
c
c     LOOPS
c          
      DO isg = 1, nskseg
         ip1 = isgmt(1, isg)
         ip2 = isgmt(2, isg)
         deter = a(3,3,isg)-a(1,3,isg)*a(1,3,isg)/a(1,1,isg)
         DO ip = ip1, ip2 - 1
            x1 = xpnt(ip  , isct)
            x2 = xpnt(ip+1, isct)
            y1 = ypnt(ip  , isct)
            y2 = ypnt(ip+1, isct)
            dss = DSQRT((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))
            area(isg) = area(isg) + dss*q1(isg,ip)/deter
            aa  (isg) = aa  (isg) + dss           /deter
         END DO
      END DO
      
      IF ( nwebs/=0 ) THEN
         DO iweb  = 1, nwebs
            isg   = webseg(iweb, isct)
            deter = a(3,3,isg)-a(1,3,isg)*a(1,3,isg)/a(1,1,isg)
            ip1   = isgmt(1, isg)
            ip2   = isgmt(2, isg)
            x1    = xpnt(ip1,isct)
            x2    = xpnt(ip2,isct)
            y1    = ypnt(ip1,isct)
            y2    = ypnt(ip2,isct)
            dss = DSQRT((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))
            dss   =  dss  - corweb
            area(isg) = dss*q1(isg,ip1)/deter
            aa  (isg) = dss            /deter
         END DO
      END IF
      
      DO il = 1, nloop
         DO is = 1, nlpsg(il)
            lps = abs(lpseg(is,il))
            area0 (il) = area0 (il) - area(lps)*lpseg(is, il)/lps
            aaloop(il) = aaloop(il) + aa  (lps)
         END DO
      END DO
c      
c     ... solve the system for the loops
c
      IF ( nwebs==0 ) aahstr(1) = area0(1)/aaloop(1)
      IF ( nwebs==1 ) THEN
         aaa(1, 1) =  aaloop(1)
         aaa(2, 2) =  aaloop(2)
         aaa(1, 2) = -aa(webseg(1,isct))
         aaa(2, 1) =  aaa(1, 2)
         aahstr(1) =  area0(1)
         aahstr(2) =  area0(2)
c        ...calculate equivalent aahstr for each loop
         CALL dsifa(aaa, 4, nloop, kqvt, info)
         CALL dsisl(aaa, 4, nloop, kqvt, aahstr)
      END IF
      IF ( nwebs==2 ) THEN
         aaa(1, 1) =  aaloop(1)
         aaa(2, 2) =  aaloop(2)
         aaa(3, 3) =  aaloop(3)
         aaa(1, 2) = -aa(webseg(1,isct))
         aaa(2, 1) =  aaa(1, 2)
         aaa(2, 3) = -aa(webseg(2,isct))
         aaa(3, 2) =  aaa(2, 3)
         aahstr(1) =  area0(1)
         aahstr(2) =  area0(2)
         aahstr(3) =  area0(3)
         CALL dsifa(aaa, 4, nloop, kqvt, info)
         CALL dsisl(aaa, 4, nloop, kqvt, aahstr)
         write(1,*)' Q01, Q02, Q033 '
         write(1,*) aahstr(1),aahstr(2),aahstr(3)
         write(1,*)
      END IF
c
c     Q-correction along segments
c
      DO isg = 1, nskseg
         DO il = 1, nloop
            DO is = 1, nlpsg(il)
               IF ( isg==abs(lpseg(is,il)) ) iloop = il
            END DO
         END DO
 
         ip1 = isgmt(1, isg)
         ip2 = isgmt(2, isg)
         DO ip = ip1, ip2-1  
         q2(isg,ip) = q1(isg,ip)+aahstr(iloop)
         END DO
      END DO 
      
      IF ( nwebs/=0 ) THEN
         DO iweb = 1, nwebs
            isg  = webseg(iweb, isct)
            ip1 = isgmt(1, isg)
            ip2 = isgmt(2, isg)
            x1 = xpnt(ip1, isct)
            x2 = xpnt(ip2, isct)
            y1 = ypnt(ip1, isct)
            y2 = ypnt(ip2, isct)
            IF ( iweb==1 ) then
            q2(isg,ip1) = q1(isg,ip1)-(aahstr(1)-aahstr(2))
            endif
            IF ( iweb==2 ) then
            q2(isg,ip1) = q1(isg,ip1)+(aahstr(2)-aahstr(3))
            endif
         END DO
      END IF
c
c     Torque Calculation
c      
      To0=0 
c      
      DO isg = 1, nskseg
         ip1 = isgmt(1, isg)
         ip2 = isgmt(2, isg)
         DO ip = ip1, ip2 - 1
            x1 = xpnt(ip  , isct)
            x2 = xpnt(ip+1, isct)
            y1 = ypnt(ip  , isct)
            y2 = ypnt(ip+1, isct)
            CALL mdllin(x1, y1, x2, y2, dss, ro, ros, thk, 1)
            rzeta = -ro(1)*ros(2) + ro(2)*ros(1)
            To0   =  To0 - rzeta*q2(isg,ip)*dss
         END DO
      END DO     
c
      IF ( nwebs/=0 ) THEN
         DO iweb = 1, nwebs
            isg = webseg(iweb, isct)
            ip1 = isgmt(1, isg)
            ip2 = isgmt(2, isg)
            x1  = xpnt(ip1,isct)
            x2  = xpnt(ip2,isct)
            y1  = ypnt(ip1,isct)
            y2  = ypnt(ip2,isct)
            CALL mdllin(x1, y1, x2, y2, dss, ro, ros, thk, 1)
            dss   =  dss  - corweb
            rzeta = -ro(1)*ros(2) + ro(2)*ros(1)
            To0   =  To0 - rzeta*q2(isg,ip1)*dss
         END DO
      END IF
c   
c    NEXT TWO LOOPS ONLY FOR PRINTOUT      
           
      DO isg = 1, nskseg
         ip1 = isgmt(1, isg)
         ip2 = isgmt(2, isg)
               DO ip = ip1, ip2-1 
            x1 = xpnt(ip  , isct)
            x2 = xpnt(ip+1, isct)
            y1 = ypnt(ip  , isct)
            y2 = ypnt(ip+1, isct)
            CALL mdllin(x1, y1, x2, y2, dss, ro, ros, thk, 1)
      WRITE (1,99002)isg,ip,ro(1),ro(2),q0(isg,ip),q1(isg,ip),q2(isg,ip)
               END DO
      END DO
      
      IF ( nwebs/=0 ) THEN
         DO iweb = 1, nwebs
            isg = webseg(iweb, isct)
            ip1 = isgmt(1, isg)
            ip2 = isgmt(2, isg)
            x1    = xpnt(ip1,isct)
            y1    = ypnt(ip1,isct)
            x2    = xpnt(ip2,isct)
            y2    = ypnt(ip2,isct)
            CALL mdllin(x1, y1, x2, y2, dss, ro, ros, thk, 1)
      WRITE (1,99002) isg,ip1,ro(1),ro(2),
     &             q0(isg,ip1),q1(isg,ip1),q2(isg,ip1)
         END DO
      END IF
      
      close(1)
99002 format(2i5,5(1x,e12.4))    
      
      END SUBROUTINE SHEARFLO
