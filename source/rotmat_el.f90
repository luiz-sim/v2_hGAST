!--------------------------------------------------------------------------------
!
!  Subrouine : ROTMAT_el ----------------
!
!  Communication of the q's and their derivatives in TIME
!  Derivation of rotation matrices
!
!--------------------------------------------------------------------------------
 Subroutine ROTMAT_el
!--------------------------------------------------------------------------------

 Use Cbeam
 Use Rotate

   implicit none

   real(8) :: ANG, Hlen
   integer :: nn_el, nbod, nnsub, nn, kdof, NDF, nb0_el
   integer :: i, I1, I2, NQ, NQ0, NQR, m1, m2, NNODESUBB, nnodT


   if ((NBODBTWT_el == 0).or.(NQS == 0)) return

!--- Zero transformation Matrices
   if (ICASE_el/=4) nnodT = NNODTWT_el
   if (ICASE_el==4) nnodT = NNODT_el

   do nn_el = 1, nnodT
      nq    = QSnode(nn_el)%Tot_el(2)

      transf_mat(nn_el)%A0_el     (1:nq,1:3,1:3) = 0.d0;
      transf_mat(nn_el)%AT0_el    (1:nq,1:3,1:3) = 0.d0;
      transf_mat(nn_el)%ATDA0_el  (1:nq,1:3,1:3) = 0.d0;
      transf_mat(nn_el)%ATDA1_el  (1:nq,1:3,1:3) = 0.d0;
      transf_mat(nn_el)%ATDDA0_el (1:nq,1:3,1:3) = 0.d0;
      transf_mat(nn_el)%ATDDA1_el (1:nq,1:3,1:3) = 0.d0;
      transf_mat(nn_el)%ATDDA2_el (1:nq,1:3,1:3) = 0.d0;
      transf_mat(nn_el)%R0_el     (1:nq,    1:3) = 0.d0;
      transf_mat(nn_el)%ATDDR0_el (1:nq,    1:3) = 0.d0;
      transf_mat(nn_el)%ATDDR1_el (1:nq,    1:3) = 0.d0;
      transf_mat(nn_el)%ATDDR2_el (1:nq,    1:3) = 0.d0;
      transf_mat(nn_el)%ATDR0_el  (1:nq,    1:3) = 0.d0;
      transf_mat(nn_el)%ATDR1_el  (1:nq,    1:3) = 0.d0;
   enddo

!-- only for blades --- (WT here)
   do nn_el = 1, NNODTWT_el
      transf_mat(nn_el)%ATPA_el  (      1:3,1:3) = 0.d0;
!!    transf_mat(nn_el)%DATPA_el (      1:3,1:3) = 0.d0;
!!    transf_mat(nn_el)%DDATPA_el(      1:3,1:3) = 0.d0;
   enddo



!---- For the Tower
      nbod      = NBLADE_el+2
      NNODESUBB = nsub_tow+1

      do i  = 1, IQfl_rot
         I2 = NDFBT_el + IQfl_tr + i
         call SET_QBr ( i, i, UT_el (I2), UT1_el (I2), UT2_el (I2) ) !rotations: floater x3
      enddo

         call SET_QBr ( IQfl_rot + 1, 1, PIhalf, 0.d0, 0.d0 )        !rotation : Ax(pi/2)

      NDFROT0 = IQfl_rot + 1
      NDFROT  = NDFROT0  + 3*NNODESUBB

      do i  = 1, IQfl_tr
         I2 = NDFBT_el + i
         call SET_QBt ( i, i, UT_el (I2), UT1_el (I2), UT2_el (I2) ) !translations: floater x3
      enddo

         call SET_QBt ( IQfl_tr + 1 , 2, Htow0 , 0.d0, 0.d0 )        !translation : Ry(Htow0)

      NDFTRA0 = IQfl_tr + 1
      NDFTRA  = NDFTRA0 + 3*NNODESUBB

      do nn   = 1, NNODESUBB
         kdof = (nn-1)*3
         NDF  =  ACCU%NQSACC_el (nbod-1) + ACCU%NQSBACC_el(nbod, nn-1)

         I1   = NDFROT0  + kdof
         I2   = NDFBT_el + NQSP + NDF

         call SET_QBr ( I1+1, 1, UT_el (I2+4), UT1_el (I2+4), UT2_el (I2+4) )
         call SET_QBr ( I1+2, 3, UT_el (I2+6), UT1_el (I2+6), UT2_el (I2+6) )
         call SET_QBr ( I1+3, 2, UT_el (I2+5), UT1_el (I2+5), UT2_el (I2+5) )

         if (nn > 1) then
             nb0_el = body(nbod)%NBODTGLB_el(nn-1)
             Hlen   = subbody(nb0_el)%ALENGB_el
         else !nn=1
             Hlen   = 0.d0
         endif

         I1   = NDFTRA0 + kdof

         call SET_QBt ( I1+1, 1, UT_el (I2+1)       , UT1_el (I2+1), UT2_el (I2+1) )
         call SET_QBt ( I1+2, 2, UT_el (I2+2) + Hlen, UT1_el (I2+2), UT2_el (I2+2) )
         call SET_QBt ( I1+3, 3, UT_el (I2+3)       , UT1_el (I2+3), UT2_el (I2+3) )
      enddo !nn

      NDFROT0tow = NDFROT0
      NDFTRA0tow = NDFTRA0


      call ROT_TRAN_BODYB (0,     1,NDFROT, 1,NDFTRA)


      NQ  = NDFROT
      NQ0 = 1

      call AA0B0C0       ( NQ0, NQ )
      call R_DR_DDR_INIT ( nbod    )
!old  call A_DA_DDA ( NDFROT0 ) !for 1st R done to R_DR_DDR_INIT


      do nnsub  = 1, NNODESUBB
         nn_el  = body(nbod)%NNODTGLB_el(nnsub)

         NDFROT = NDFROT0 + nnsub*3
         NDFTRA = NDFTRA0 + nnsub*3

         NQR    = NDFROT0 + (nnsub-1)*3
         m1     = NDFTRA0 + (nnsub-1)*3+1
         m2     = NDFTRA0 +  nnsub   *3
         NQ     = NDFROT

         call R_DR_DDR   ( NQR, m1, m2 ) !calculate R using previous A
         call A_DA_DDA   ( NQ          )
         call ATDA_ATDDA ( NQ ,nbod    )
         call ATDR_ATDDR

         call ROTMAT_TRANS_tower     ( nn_el, nbod, nnsub )
         call ROTMAT_TRANS_tmptoglob ( nn_el )
      enddo !nnsub


!---- For the Shaft
      nbod      = NBLADE_el+1
      NNODESUBB = nsub_sh+1
      I1        = NDFBT_el + NQSY
      I2        = NDFBT_el + NQSW

     if      (IAPPL_el  == 0 ) then
      call SET_QBr ( IQtow_rot_T + 1, 3, PIhalf   , 0.d0      , 0.d0       ) !rotation: Az(pi/2)
     elseif  (IAPPL_el  == 1.or. IAPPL_el == 2 ) then
      call SET_QBr ( IQtow_rot_T + 1, 3, 0.d0     , 0.d0      , 0.d0       ) !rotation: Az(   ) !dummy, for VAWT and heli the shaft is the extension of the tower
     endif
      call SET_QBr ( IQtow_rot_T + 2, 1, YAW      , 0.d0      , 0.d0       ) !rotation: Ax(YAW)
      call SET_QBr ( IQtow_rot_T + 3, 1, UT_el(I1), UT1_el(I1), UT2_el(I1) ) !rotation: Ax(Yaw dof)
      call SET_QBr ( IQtow_rot_T + 4, 3, -TILT    , 0.d0      , 0.d0       ) !rotation: Az(-TILT)
      call SET_QBr ( IQtow_rot_T + 5, 2, UT_el(I2), UT1_el(I2), UT2_el(I2) ) !rotation: Ay(omega dof)
!qq   call SET_QBr ( IQtow_rot_T + 5, 2, UT_el(I2)-PHI0(1), UT1_el(I2), UT2_el(I2) ) !omega dof including initial azimuthal offset

      NDFROT0 = IQtow_rot_T + 5
      NDFROT  = NDFROT0 + 3*NNODESUBB

      call SET_QBt ( IQtow_tr_T  + 1, 1, Hshof    , 0.d0      , 0.d0       ) !translation: Rx(Hshof)
      call SET_QBt ( IQtow_tr_T  + 2, 2, Hsh      , 0.d0      , 0.d0       ) !translation: Ry(Hsh  )

      NDFTRA0 = IQtow_tr_T + 2
      NDFTRA  = NDFTRA0 + 3*NNODESUBB

      do nn   = 1, NNODESUBB
         kdof = (nn-1)*3
         NDF  =  ACCU%NQSACC_el (nbod-1) + ACCU%NQSBACC_el(nbod, nn-1)

         I1   = NDFROT0  + kdof
         I2   = NDFBT_el + NQSP + NDF

         call SET_QBr ( I1+1, 1, UT_el (I2+4), UT1_el (I2+4), UT2_el (I2+4) )
         call SET_QBr ( I1+2, 3, UT_el (I2+6), UT1_el (I2+6), UT2_el (I2+6) )
         call SET_QBr ( I1+3, 2, UT_el (I2+5), UT1_el (I2+5), UT2_el (I2+5) )

         if (nn > 1) then
             nb0_el = body(nbod)%NBODTGLB_el(nn-1)
             Hlen   = subbody(nb0_el)%ALENGB_el
         else !nn=1
             Hlen   = 0d0
         endif

         I1   = NDFTRA0 + kdof

         call SET_QBt ( I1+1, 1, UT_el (I2+1)       , UT1_el (I2+1), UT2_el (I2+1) )
         call SET_QBt ( I1+2, 2, UT_el (I2+2) + Hlen, UT1_el (I2+2), UT2_el (I2+2) )
         call SET_QBt ( I1+3, 3, UT_el (I2+3)       , UT1_el (I2+3), UT2_el (I2+3) )
      enddo !nn

      NDFROT0sha = NDFROT0
      NDFTRA0sha = NDFTRA0

      call ROT_TRAN_BODYB (0,  IQtow_rot_T+1,NDFROT, IQtow_tr_T+1,NDFTRA)

      NQ  = NDFROT
      NQ0 = IQtow_rot_T + 1

      call AA0B0C0       ( NQ0, NQ )
      call ROTMAT_YAW_el

!------ set non-rotating shaft base
      NQ  = IQtow_rot_T + 4
      call A_DA_DDA ( NQ )

      ATsftNoR_el (1:3,1:3) = AT (1:3,1:3)

      call R_DR_DDR_INIT ( nbod    )
!     call A_DA_DDA ( NDFROT0 ) !for 1st R not need, done already to R_DR_DDR_INIT


      do nnsub  = 1, NNODESUBB
         nn_el  = body(nbod)%NNODTGLB_el(nnsub)

         NDFROT = NDFROT0 + nnsub*3
         NDFTRA = NDFTRA0 + nnsub*3

         NQR    = NDFROT0 + (nnsub-1)*3
         m1     = NDFTRA0 + (nnsub-1)*3+1
         m2     = NDFTRA0 +  nnsub   *3
         NQ     = NDFROT

         call R_DR_DDR   ( NQR, m1, m2 ) !calculate R using previous A
         call A_DA_DDA   ( NQ          )
         call ATDA_ATDDA ( NQ ,nbod    )
         call ATDR_ATDDR

         call ROTMAT_TRANS_tower     ( nn_el, nbod+1, nsub_tow+1 )
         call ROTMAT_TRANS_shaft     ( nn_el, nbod  , nnsub      )
         call ROTMAT_TRANS_tmptoglob ( nn_el )
      enddo !nnsub


!---- For the Blades
   do nbod      = 1,NBLADE_el
      NNODESUBB = nsub_bl+1

     if      (IAPPL_el  == 0 ) then
      call SET_QBr ( IQsh_rot_T + 1, 3, -PIhalf    , 0.d0      , 0.d0       ) !rotation: Az(-pi/2)
      call SET_QBr ( IQsh_rot_T + 2, 2,  PIhalf    , 0.d0      , 0.d0       ) !rotation: Ay( pi/2)
      call SET_QBr ( IQsh_rot_T + 3, 3,  PHI0(nbod), 0.d0      , 0.d0       ) !rotation: Az( PHI0)
      call SET_QBr ( IQsh_rot_T + 4, 1,  CONE      , 0.d0      , 0.d0       ) !rotation: Ax( CONE)
     elseif  (IAPPL_el  == 1 ) then
!     call SET_QBr ( IQsh_rot_T + 1, 3,  0.d0      , 0.d0      , 0.d0       )
!     call SET_QBr ( IQsh_rot_T + 2, 2,  PIhalf    , 0.d0      , 0.d0       )
!     call SET_QBr ( IQsh_rot_T + 3, 2,  PHI0(nbod), 0.d0      , 0.d0       )
!     call SET_QBr ( IQsh_rot_T + 4, 1,  CONE      , 0.d0      , 0.d0       )
!------ in this way the strut will be placed horizonltaly...
      call SET_QBr ( IQsh_rot_T + 1, 3,  PIhalf    , 0.d0      , 0.d0       ) !rotation: Az( pi/2)
      call SET_QBr ( IQsh_rot_T + 2, 2,  PIhalf    , 0.d0      , 0.d0       ) !rotation: Ay( pi/2) !rotate in order to be similar to the HAWT orientation
      call SET_QBr ( IQsh_rot_T + 3, 3,  PHI0(nbod), 0.d0      , 0.d0       ) !rotation: Az( PHI0)
      call SET_QBr ( IQsh_rot_T + 4, 1,  CONE      , 0.d0      , 0.d0       ) !rotation: Ax( CONE)
     elseif  (IAPPL_el  == 2 ) then
      call SET_QBr ( IQsh_rot_T + 1, 1, -PIhalf    , 0.d0      , 0.d0       ) !rotation: Ax(-pi/2) !same to global
      call SET_QBr ( IQsh_rot_T + 2, 2,  0.d0      , 0.d0      , 0.d0       ) !rotation: Ay(     ) !dummy
      call SET_QBr ( IQsh_rot_T + 3, 3,  PHI0(nbod), 0.d0      , 0.d0       ) !rotation: Az( PHI0)
      call SET_QBr ( IQsh_rot_T + 4, 1,  CONE      , 0.d0      , 0.d0       ) !rotation: Ax( CONE)
     endif

      NDFROT0 = IQsh_rot_T + 4
      NDFROT  = NDFROT0 + 4*NNODESUBB

!!    call SET_QBt ( IQsh_tr_T  + 1, 2,  Hhub      , 0.d0      , 0.d0       )

      NDFTRA0 = IQsh_tr_T      !!+ 1
      NDFTRA  = NDFTRA0 +3*NNODESUBB

      do nn   = 1, NNODESUBB
         kdof = (nn-1)*4
         NDF  =  ACCU%NQSACC_el (nbod-1) + ACCU%NQSBACC_el(nbod, nn-1)

         I1   = NDFROT0+kdof 
         I2   = NDFBT_el + NQSP + NDF
         I    = IAX_Ben_Swe

         call SET_QBr ( I1+1, 1, UT_el (I2+4), UT1_el (I2+4), UT2_el (I2+4) )
         call SET_QBr ( I1+2, 3, UT_el (I2+6), UT1_el (I2+6), UT2_el (I2+6) )
         call SET_QBr ( I1+3, 2, UT_el (I2+5), UT1_el (I2+5), UT2_el (I2+5) )

!--------- Introduce Precurved angles
        if (I>0) then
         ANG  =  body(nbod)%PRECURV(nn) - body(nbod)%PRECURV(nn-1)
         call SET_QBr ( I1+4, I, ANG         , 0.d0         , 0.d0          )
        else
!--------- VAWT with helix: set directly the rotation matrix
         call SET_QBr ( I1+4, I, 0.d0        , 0.d0         , 0.d0          )
         AQ ( I1+4,1:3,1:3) = body(nbod)%PRECURV_mat(nn,1:3,1:3)
         AQ1( I1+4,1:3,1:3) = 0.d0;
         AQ2( I1+4,1:3,1:3) = 0.d0;
         AQ3( I1+4,1:3,1:3) = 0.d0;
        endif

         kdof = (nn-1)*3

         if (nn > 1) then
             nb0_el = body(nbod)%NBODTGLB_el(nn-1)
             Hlen   = subbody(nb0_el)%ALENGB_el
         else !nn=1
             Hlen   = Hhub
         endif

         I1   = NDFTRA0 + kdof

         call SET_QBt ( I1+1, 1, UT_el (I2+1)       , UT1_el (I2+1), UT2_el (I2+1) )
         call SET_QBt ( I1+2, 2, UT_el (I2+2) + Hlen, UT1_el (I2+2), UT2_el (I2+2) )
         call SET_QBt ( I1+3, 3, UT_el (I2+3)       , UT1_el (I2+3), UT2_el (I2+3) )
      enddo !nn


         NQ  = NDFROT
         NQ0 = IQsh_rot_T + 1
      if (nbod /= 1) &
         NQ0 = NDFROT0-1

      call ROT_TRAN_BODYB (0, NQ0,   NQ , IQsh_tr_T+1, NDFTRA)
      call AA0B0C0        ( NQ0, NQ )
      call R_DR_DDR_INIT  ( nbod    )

!------ Set ATP_el
                       NQ  = NDFROT0       !after  cone
      if (IAERO_el==1) NQ  = NDFROT0 - 1   !before cone
      call A_DA_DDA       ( NQ      )

      ATP     (     1:3,1:3) = AT   (1:3,1:3)
!c    DATP    (     1:3,1:3) = transpose(DA (1:3,1:3))
!c    DDATP   (     1:3,1:3) = transpose(DDA(1:3,1:3))

      ATP_el  (nbod,1:3,1:3) = ATP  (1:3,1:3)
!c    DATP_el (nbod,1:3,1:3) = DATP (1:3,1:3)
!c    DDATP_el(nbod,1:3,1:3) = DDATP(1:3,1:3)


!------ Set ATPWr_el
      NQ                     = NDFROT0 + 3 !+3 after pitch, +4 after pitch and precurve rotation
      call A_DA_DDA       ( NQ      )
      ATPWr_el(nbod,1:3,1:3) = AT   (1:3,1:3)

      call A_DA_DDA       ( NDFROT0 ) !for 1st R

!------ Set ATPWrG_el
!                         NQ  = NDFROT0       !after  cone
      ATPWrG_el(nbod,1:3,1:3) = AT   (1:3,1:3)

      do nnsub  = 1, NNODESUBB                    !the pitch rotation is added here
         nn_el  = body(nbod)%NNODTGLB_el(nnsub)

         NDFROT = NDFROT0 + nnsub*4
         NDFTRA = NDFTRA0 + nnsub*3

         NQR    = NDFROT0 + (nnsub-1)*4
         m1     = NDFTRA0 + (nnsub-1)*3+1
         m2     = NDFTRA0 +  nnsub   *3 
         NQ     = NDFROT

         call R_DR_DDR   ( NQR, m1, m2 ) !calculate R using previous A
         call A_DA_DDA   ( NQ          )
         call ATDA_ATDDA ( NQ ,nbod    )
         call ATDR_ATDDR

         call ROTMAT_TRANS_tower     ( nn_el, NBLADE_el+2, nsub_tow+1 )
         call ROTMAT_TRANS_shaft     ( nn_el, NBLADE_el+1, nsub_sh +1 )
         call ROTMAT_TRANS_blade     ( nn_el, nbod       , nnsub      )
         call ROTMAT_TRANS_tmptoglob ( nn_el )
      enddo !nnsub
   enddo !nbod


!---- For all jacket's bodies, in case of a flexible floater
   if ( ICASE_el /= 4 ) return


   call SET_QBr          ( IQfl_rot + 1, 1, 1.d0, 0.d0, 0.d0 ) !just set [der,dder]QB=0
   call SET_QBt          ( IQfl_tr  + 1, 1, 1.d0, 0.d0, 0.d0 )
   call R_DR_DDR_INIT_ja (                                 1 ) !Store 3 floater's R,DR,DDR to reduce cost


   do nbod    = NBODBTWT_el+1,NBODBT_el
      nn_el   = body(nbod)%NNODTGLB_el(1)

      NDFROT0 = IQfl_rot
      NDFROT  = IQfl_rot + 1
      NDFTRA0 = IQfl_tr
      NDFTRA  = IQfl_tr  + 1

      NQ0     = NDFROT0
      NQ      = NDFROT
      m1      = NDFTRA

      call ROT_TRAN_BODYB_ja      ( transf_mat(nn_el)%Aja_el(1:3,1:3),        &
                                    transf_mat(nn_el)%Rja_el(1:3    ), NQ, m1   )
      call AA0B0C0                ( NQ0,  NQ )
      call R_DR_DDR_INIT_ja       (        2 )
      call A_DA_DDA               ( NQ       )
      call ATDA_ATDDA             ( NQ ,nbod )
      call ATDR_ATDDR
      call ROTMAT_TRANS_ja        ( nn_el    )
      call ROTMAT_TRANS_tmptoglob ( nn_el    )
   enddo !nbod


 END Subroutine ROTMAT_el
!--------------------------------------------------------------------------------
!
!  Subrouine :ROTMAT_TRANS_tower -------------
!
!  Transfer local matrices to global variables
!--------------------------------------------------------------------------------
 Subroutine ROTMAT_TRANS_tower(nn_el,nbod,nnsub)
!--------------------------------------------------------------------------------

 Use Cbeam
 Use Rotate

   implicit none

   integer :: nn_el, nbod, nnsub, ikp(3), nn, i, kk, nqglob, nqloc


   do i=1,NQS
      rot_tmp(i)%R0     (1:3    ) = 0.d0
      rot_tmp(i)%ATDR0  (1:3    ) = 0.d0
      rot_tmp(i)%ATDR1  (1:3    ) = 0.d0
      rot_tmp(i)%ATDDR0 (1:3    ) = 0.d0
      rot_tmp(i)%ATDDR1 (1:3    ) = 0.d0
      rot_tmp(i)%ATDDR2 (1:3    ) = 0.d0
      rot_tmp(i)%A0     (1:3,1:3) = 0.d0
      rot_tmp(i)%AT0    (1:3,1:3) = 0.d0
      rot_tmp(i)%ATDA0  (1:3,1:3) = 0.d0
      rot_tmp(i)%ATDA1  (1:3,1:3) = 0.d0
      rot_tmp(i)%ATDDA0 (1:3,1:3) = 0.d0
      rot_tmp(i)%ATDDA1 (1:3,1:3) = 0.d0
      rot_tmp(i)%ATDDA2 (1:3,1:3) = 0.d0
   enddo


   ikp(1)=4
   ikp(2)=6
   ikp(3)=5

   transf_mat(nn_el)%R_el      (1:3    ) = R      (1:3    )
   transf_mat(nn_el)%DR_el     (1:3    ) = DR     (1:3    )
   transf_mat(nn_el)%ATDR_el   (1:3    ) = ATDR   (1:3    )
   transf_mat(nn_el)%ATDDR_el  (1:3    ) = ATDDR  (1:3    )
   transf_mat(nn_el)%A_el      (1:3,1:3) = A      (1:3,1:3)
   transf_mat(nn_el)%DA_el     (1:3,1:3) = DA     (1:3,1:3)
   transf_mat(nn_el)%AT_el     (1:3,1:3) = AT     (1:3,1:3)
   transf_mat(nn_el)%ATDA_el   (1:3,1:3) = ATDA   (1:3,1:3)
   transf_mat(nn_el)%ATDDA_el  (1:3,1:3) = ATDDA  (1:3,1:3)

   transf_mat(nn_el)%ATPA_el  (1:3,1:3) = ATPA    (1:3,1:3)
!c transf_mat(nn_el)%DATPA_el (1:3,1:3) = DATPA   (1:3,1:3)
!c transf_mat(nn_el)%DDATPA_el(1:3,1:3) = DDATPA  (1:3,1:3)


!--- Floating Rotations         !   (QT_el) = QBr
   do i = 1, IQfl_rot
      call ROTMAT_TRANS (nn_el, i + IQfl_tr , i        , 0)
   enddo

!--- Floating Translations
   do i = 1, IQfl_tr
      call ROTMAT_TRANS (nn_el, i           , i+IQrot_T, 1)
   enddo

!--- Tower's Q's [from subbodies]
   do nn = 1, nnsub
      do kk  = 1, 3 !subbody translations (go 1st to global Qs)
         nqglob = NQSP + ACCU%NQSACC_el (nbod-1) + ACCU%NQSBACC_el(nbod, nn-1) + kk
         nqloc  = IQrot_T + NDFTRA0tow + (nn-1)*3 + kk

         call ROTMAT_TRANS (nn_el, nqglob, nqloc, 1)
      enddo
      do kk = 1, 3 !subbody rotations (go after translations to global Qs)
         nqglob = NQSP + ACCU%NQSACC_el (nbod-1) + ACCU%NQSBACC_el(nbod, nn-1) + ikp(kk)
         nqloc  = NDFROT0tow + (nn-1)*3 + kk 

         call ROTMAT_TRANS (nn_el, nqglob, nqloc, 0)
      enddo
   enddo !nn
 

 END Subroutine ROTMAT_TRANS_tower
!--------------------------------------------------------------------------------
!
!  Subrouine :ROTMAT_TRANS_shaft -------------
!
!  Transfer local matrices to global variables
!--------------------------------------------------------------------------------
 Subroutine ROTMAT_TRANS_shaft(nn_el,nbod,nnsub)
!--------------------------------------------------------------------------------

 Use Cbeam
 Use Rotate

   implicit none

   integer :: nn_el, nbod, nnsub, ikp(3), nn, kk, nqglob, nqloc


   ikp(1)=4
   ikp(2)=6
   ikp(3)=5

!--- Omega                      (QT_el) -> QBr
   call ROTMAT_TRANS    (nn_el, NQSW, 5+IQtow_rot_T , 0    )

!--- Yaw free Movement equation
   call ROTMAT_TRANS    (nn_el, NQSY, 3+IQtow_rot_T , 0    )

!--- Shaft's Q's [from subbodies]
   do nn = 1, nnsub
      do kk  = 1, 3 !subbody translations (go 1st to global Qs)
         nqglob = NQSP + ACCU%NQSACC_el (nbod-1) + ACCU%NQSBACC_el(nbod, nn-1) + kk
         nqloc  = IQrot_T + NDFTRA0sha + (nn-1)*3 + kk

         call ROTMAT_TRANS ( nn_el,nqglob,nqloc, 1 )
      enddo
      do kk = 1, 3 !subbody rotations (go after translations to global Qs)
         nqglob = NQSP + ACCU%NQSACC_el (nbod-1) + ACCU%NQSBACC_el(nbod, nn-1) + ikp(kk)
         nqloc  = NDFROT0sha + (nn-1)*3 + kk 

         call ROTMAT_TRANS ( nn_el,nqglob,nqloc, 0 )
      enddo
   enddo !nn
 

 END Subroutine ROTMAT_TRANS_shaft
!--------------------------------------------------------------------------------
!
!  Subrouine :ROTMAT_TRANS_blade -------------
!
!  Transfer local matrices to global variables
!--------------------------------------------------------------------------------
 Subroutine ROTMAT_TRANS_blade(nn_el,nbod,nnsub)
!--------------------------------------------------------------------------------

 Use Cbeam
 Use Rotate

   implicit none

   integer :: nn_el, nbod, nnsub, ikp(3), nn, kk, nqglob, nqloc


   ikp(1)=4
   ikp(2)=6
   ikp(3)=5

!--- Blade's Q's [from subbodies]
   do nn = 1, nnsub
      do kk  = 1, 3 !subbody translations (go 1st to global Qs)
         nqglob = NQSP + ACCU%NQSACC_el (nbod-1) + ACCU%NQSBACC_el(nbod, nn-1) + kk
         nqloc  = IQrot_T + NDFTRA0 + (nn-1)*3 + kk

         call ROTMAT_TRANS ( nn_el,nqglob,nqloc, 1 )
      enddo
      do kk = 1, 3 !subbody rotations (go after translations to global Qs)
         nqglob = NQSP + ACCU%NQSACC_el (nbod-1) + ACCU%NQSBACC_el(nbod, nn-1) + ikp(kk)
         nqloc  = NDFROT0 + (nn-1)*4 + kk 

         call ROTMAT_TRANS ( nn_el,nqglob,nqloc, 0 )
      enddo
   enddo !nn
 

 END Subroutine ROTMAT_TRANS_blade
!--------------------------------------------------------------------------------
 Subroutine ROTMAT_TRANS_ja (nn_el)
!--------------------------------------------------------------------------------

 Use Cbeam
 Use Rotate

   implicit none

   integer :: nn_el, i


   transf_mat(nn_el)%R_el      (1:3    ) = R      (1:3    )
   transf_mat(nn_el)%DR_el     (1:3    ) = DR     (1:3    )
   transf_mat(nn_el)%ATDR_el   (1:3    ) = ATDR   (1:3    )
   transf_mat(nn_el)%ATDDR_el  (1:3    ) = ATDDR  (1:3    )
   transf_mat(nn_el)%A_el      (1:3,1:3) = A      (1:3,1:3)
   transf_mat(nn_el)%DA_el     (1:3,1:3) = DA     (1:3,1:3)
   transf_mat(nn_el)%AT_el     (1:3,1:3) = AT     (1:3,1:3)
   transf_mat(nn_el)%ATDA_el   (1:3,1:3) = ATDA   (1:3,1:3)
   transf_mat(nn_el)%ATDDA_el  (1:3,1:3) = ATDDA  (1:3,1:3)
!qqtransf_mat(nn_el)%ATPA_el   (1:3,1:3) = ATPA   (1:3,1:3)
!c transf_mat(nn_el)%DATPA_el  (1:3,1:3) = DATPA  (1:3,1:3)
!c transf_mat(nn_el)%DDATPA_el (1:3,1:3) = DDATPA (1:3,1:3)

!--- Floating Rotations
   do i = 1, IQfl_rot
      call ROTMAT_TRANS (nn_el, i + IQfl_tr , i           , 0)
   enddo

!--- Floating Translations
   do i = 1, IQfl_tr
      call ROTMAT_TRANS (nn_el, i           , i+IQrot_T   , 1)
   enddo


 END Subroutine ROTMAT_TRANS_ja
!--------------------------------------------------------------------------------
!
!  Subroutine :ROTMAT_TRANS     -------------
!
!  Transfer local matrices to global variables
!
!--------------------------------------------------------------------------------
 Subroutine ROTMAT_TRANS_tmptoglob ( nn_el )
!--------------------------------------------------------------------------------

 Use Cbeam
 Use Rotate

   implicit none

   integer :: nn_el, nq, iq


   do iq = 1, QSnode(nn_el)%Tot_el(2)
      nq = QSnode(nn_el)%IORD_el(iq)

      transf_mat(nn_el)%R0_el     (iq,1:3    ) = rot_tmp(nq)%R0     (1:3    )
      transf_mat(nn_el)%ATDR0_el  (iq,1:3    ) = rot_tmp(nq)%ATDR0  (1:3    )
      transf_mat(nn_el)%ATDR1_el  (iq,1:3    ) = rot_tmp(nq)%ATDR1  (1:3    )
      transf_mat(nn_el)%ATDDR0_el (iq,1:3    ) = rot_tmp(nq)%ATDDR0 (1:3    )
      transf_mat(nn_el)%ATDDR1_el (iq,1:3    ) = rot_tmp(nq)%ATDDR1 (1:3    )
      transf_mat(nn_el)%ATDDR2_el (iq,1:3    ) = rot_tmp(nq)%ATDDR2 (1:3    )
      transf_mat(nn_el)%A0_el     (iq,1:3,1:3) = rot_tmp(nq)%A0     (1:3,1:3)
      transf_mat(nn_el)%AT0_el    (iq,1:3,1:3) = rot_tmp(nq)%AT0    (1:3,1:3)
      transf_mat(nn_el)%ATDA0_el  (iq,1:3,1:3) = rot_tmp(nq)%ATDA0  (1:3,1:3)
      transf_mat(nn_el)%ATDA1_el  (iq,1:3,1:3) = rot_tmp(nq)%ATDA1  (1:3,1:3)
      transf_mat(nn_el)%ATDDA0_el (iq,1:3,1:3) = rot_tmp(nq)%ATDDA0 (1:3,1:3)
      transf_mat(nn_el)%ATDDA1_el (iq,1:3,1:3) = rot_tmp(nq)%ATDDA1 (1:3,1:3)
      transf_mat(nn_el)%ATDDA2_el (iq,1:3,1:3) = rot_tmp(nq)%ATDDA2 (1:3,1:3)
   enddo


 END Subroutine ROTMAT_TRANS_tmptoglob
!--------------------------------------------------------------------------------
!
!  Subroutine :ROTMAT_TRANS     -------------
!
!  Transfer local matrices to global variables
!
!--------------------------------------------------------------------------------
 Subroutine ROTMAT_TRANS (nn_el, nqglob, nqloc, IRR)
!--------------------------------------------------------------------------------

 Use Cbeam
 Use Rotate

   implicit none

   integer :: nn_el, nqglob, nqloc, IRR


      rot_tmp(nqglob)%R0     (1:3    ) = R0      (nqloc,1:3    )
      rot_tmp(nqglob)%ATDR0  (1:3    ) = ATDR0   (nqloc,1:3    )
      rot_tmp(nqglob)%ATDR1  (1:3    ) = ATDR1   (nqloc,1:3    )
      rot_tmp(nqglob)%ATDDR0 (1:3    ) = ATDDR0  (nqloc,1:3    )
      rot_tmp(nqglob)%ATDDR1 (1:3    ) = ATDDR1  (nqloc,1:3    )
      rot_tmp(nqglob)%ATDDR2 (1:3    ) = ATDDR2  (nqloc,1:3    )

   if (IRR == 1) then
      rot_tmp(nqglob)%A0     (1:3,1:3) = 0.d0;
      rot_tmp(nqglob)%AT0    (1:3,1:3) = 0.d0;
      rot_tmp(nqglob)%ATDA0  (1:3,1:3) = 0.d0;
      rot_tmp(nqglob)%ATDA1  (1:3,1:3) = 0.d0;
      rot_tmp(nqglob)%ATDDA0 (1:3,1:3) = 0.d0;
      rot_tmp(nqglob)%ATDDA1 (1:3,1:3) = 0.d0;
      rot_tmp(nqglob)%ATDDA2 (1:3,1:3) = 0.d0;
   else
      rot_tmp(nqglob)%A0     (1:3,1:3) = A0      (nqloc,1:3,1:3)
      rot_tmp(nqglob)%AT0    (1:3,1:3) = AT0     (nqloc,1:3,1:3)
      rot_tmp(nqglob)%ATDA0  (1:3,1:3) = ATDA0   (nqloc,1:3,1:3)
      rot_tmp(nqglob)%ATDA1  (1:3,1:3) = ATDA1   (nqloc,1:3,1:3)
      rot_tmp(nqglob)%ATDDA0 (1:3,1:3) = ATDDA0  (nqloc,1:3,1:3)
      rot_tmp(nqglob)%ATDDA1 (1:3,1:3) = ATDDA1  (nqloc,1:3,1:3)
      rot_tmp(nqglob)%ATDDA2 (1:3,1:3) = ATDDA2  (nqloc,1:3,1:3)
   endif


 END Subroutine ROTMAT_TRANS
!--------------------------------------------------------------------------------
 Subroutine SET_QBr ( IQB, IAX, QB, DQB, DDQB )
!--------------------------------------------------------------------------------

 Use Rotate

   implicit none

   real(8) :: QB,DQB,DDQB
   integer :: IQB,IAX


   IAX_QBr(IQB) = IAX
   QBr    (IQB) = QB
   derQBr (IQB) = DQB
   dderQBr(IQB) = DDQB


 END Subroutine SET_QBr
!--------------------------------------------------------------------------------
 Subroutine SET_QBt ( IQB, IAX, QB, DQB, DDQB )
!--------------------------------------------------------------------------------

 Use Rotate

   implicit none

   real(8) :: QB,DQB,DDQB
   integer :: IQB,IAX


   IAX_QBt(IQB) = IAX
   QBt    (IQB) = QB
   derQBt (IQB) = DQB
   dderQBt(IQB) = DDQB


 END Subroutine SET_QBt
!--------------------------------------------------------------------------------
 Subroutine ROTMAT_YAW_el
!--------------------------------------------------------------------------------

 Use Cbeam
 Use Rotate

   implicit none

   integer :: NQ


   AT0yaw_el = 0.0;
!! A0yaw_el  = 0.0;

   if (IACT_YAW == 0) goto 1


!--- 2 rot (Az(pi/2), Ax(YAW fixed))
   NQ = IQtow_rot_T + 2

   call A_DA_DDA ( NQ )

   call ROTMAT_TRANS_tower_yaw ( NBLADE_el+2, nsub_tow+1)


!--- 3 rot (Az(pi/2), Ax(YAW fixed), Ax(YAW elastic))
 1 NQ = IQtow_rot_T + 3

   call A_DA_DDA ( NQ )

   Anacelle_el (1:3,1:3) = A (1:3,1:3)


 END Subroutine ROTMAT_YAW_el
!--------------------------------------------------------------------------------
!
!  Subrouine :ROTMAT_TRANS_tower_yaw --------
!
!  Transfer local matrices to global variables
!--------------------------------------------------------------------------------
 Subroutine ROTMAT_TRANS_tower_yaw(nbod,nnsub)
!--------------------------------------------------------------------------------

 Use Cbeam
 Use Rotate

   implicit none

   integer :: nbod, nnsub, ikp(3), nn, i, kk, nqglob, nqloc


!! Ayaw_el  (1:3,1:3) = A  (1:3,1:3)
   ATyaw_el (1:3,1:3) = AT (1:3,1:3)


!--- Floating Rotations    (QT_el) = QBr
   do i = 1, IQfl_rot
      nqglob = i + IQfl_tr
      nqloc  = i
      AT0yaw_el (nqglob,1:3,1:3) = AT0 (nqloc,1:3,1:3)
   enddo

!--- All subbodies
   ikp(1)=4
   ikp(2)=6
   ikp(3)=5

   do nn = 1, nnsub
      do kk = 1, 3 !subbody rotations (go after translations to global Qs)
         nqglob = NQSP + ACCU%NQSACC_el (nbod-1) + ACCU%NQSBACC_el(nbod, nn-1) + ikp(kk)
         nqloc  = NDFROT0tow + (nn-1)*3 + kk 

!!       A0yaw_el  (nqglob,1:3,1:3) = A0  (nqloc,1:3,1:3)
         AT0yaw_el (nqglob,1:3,1:3) = AT0 (nqloc,1:3,1:3)
      enddo
   enddo !nn
 

 END Subroutine ROTMAT_TRANS_tower_yaw
!--------------------------------------------------------------------------------
 Subroutine ROTMAT_TRANS_floater
!--------------------------------------------------------------------------------

 Use Cbeam
 Use Rotate

   implicit none


      R_float (        1:3) =    R (        1:3)
     DR_float (        1:3) =   DR (        1:3)
    DDR_float (        1:3) =  DDR (        1:3)

     AT_float (    1:3,1:3) =    AT(    1:3,1:3)
      A_float (    1:3,1:3) =    A (    1:3,1:3)
     DA_float (    1:3,1:3) =   DA (    1:3,1:3)
    DDA_float (    1:3,1:3) =  DDA (    1:3,1:3)
     A0_float (1:3,1:3,1:3) =   A0 (1:3,1:3,1:3)
    DA0_float (1:3,1:3,1:3) =  DA0 (1:3,1:3,1:3)
   DDA0_float (1:3,1:3,1:3) = DDA0 (1:3,1:3,1:3)


 END Subroutine ROTMAT_TRANS_floater
