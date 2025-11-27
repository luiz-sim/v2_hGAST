!-----------------------------------------------------------------------
!
!-- Subroutine :BOD2BOD_LOADS
!   Communication of loads resulting from body to body
!   interconnection 
!
!-----------------------------------------------------------------------
 Subroutine BOD2BOD_LOADSJack
!-----------------------------------------------------------------------

 Use Cbeam

   implicit None

   integer, parameter :: nja=6
   real(8) :: AMLOC0   (NDFPEM_el,NDFPEM_el)
   real(8) :: ACLOC0   (NDFPEM_el,NDFPEM_el)
   real(8) :: AKLOC0   (NDFPEM_el,NDFPEM_el)
   real(8) :: AFLOC0   (NDFPEM_el          )
   real(8) :: AMLOCNQ0 (NDFPEM_el,nja      )
   real(8) :: ACLOCNQ0 (NDFPEM_el,nja      )
   real(8) :: AKLOCNQ0 (NDFPEM_el,nja      )
   real(8) :: AMB2B    (6,NDFPEM_el)
   real(8) :: ACB2B    (6,NDFPEM_el)
   real(8) :: AKB2B    (6,NDFPEM_el)
   real(8) :: AMB2BNQ  (6,nja      )
   real(8) :: ACB2BNQ  (6,nja      )
   real(8) :: AKB2BNQ  (6,nja      )
   real(8) :: AFB2B    (6          ) 
   real(8) :: AMB2BR   (6,NDFPEM_el)
   real(8) :: ACB2BR   (6,NDFPEM_el)
   real(8) :: AKB2BR   (6,NDFPEM_el)
   real(8) :: AMB2BNQR (6,nja      )
   real(8) :: ACB2BNQR (6,nja      )
   real(8) :: AKB2BNQR (6,nja      )
   real(8) :: AFB2BR   (6          )
   real(8) :: ATRANS(6), COEFM, COEFK
   integer :: nbod,nbsub,nb_el,nn_el,nbc  ,iel,innb_el,ndof,i,iloc
   integer :: mnod,mn,mdof,j,jj,indx,kdof,iecon,inibcon,nn,jloc
   integer :: NQS_ja,nq,ibcon,ibsubcon,isbcon


                      NQS_ja = 6
   if (ICASE_el /= 4) NQS_ja = 0

   do    nbod  = 1+NBODBTWT_el, NBODBT_el
         COEFM = COEFM_el(nbod)
         COEFK = COEFK_el(nbod)
      do nbsub = 1, body(nbod)%NBODSUB_el
         nb_el = body(nbod)%NBODTGLB_el(nbsub)
         nn_el = body(nbod)%NNODTGLB_el(nbsub)

         do nbc     = 1, boundco (nb_el)%nbcpb
            iel     = boundco (nb_el)%nel (nbc)
            innb_el = boundco (nb_el)%nod (nbc)

            AMB2B   = 0.d0;
            ACB2B   = 0.d0;
            AKB2B   = 0.d0;
            AFB2B   = 0.d0;
            AMB2BNQ = 0.d0;
            ACB2BNQ = 0.d0;
            AKB2BNQ = 0.d0;

            AMLOC0    ( 1:NDFPEM_el, 1:NDFPEM_el ) = subbody(nb_el)%AMLOC_el   ( iel, 1:NDFPEM_el, 1:NDFPEM_el )
            ACLOC0    ( 1:NDFPEM_el, 1:NDFPEM_el ) = subbody(nb_el)%ACLOC_el   ( iel, 1:NDFPEM_el, 1:NDFPEM_el )
            AKLOC0    ( 1:NDFPEM_el, 1:NDFPEM_el ) = subbody(nb_el)%AKLOC_el   ( iel, 1:NDFPEM_el, 1:NDFPEM_el )
            AFLOC0    ( 1:NDFPEM_el              ) = subbody(nb_el)%AFLOC_el   ( iel, 1:NDFPEM_el              )

            do nq=1,NQS_ja
              AMLOCNQ0  ( 1:NDFPEM_el, nq        ) = subbody(nb_el)%AMLOCNQ_el ( iel, 1:NDFPEM_el, nq          )
              ACLOCNQ0  ( 1:NDFPEM_el, nq        ) = subbody(nb_el)%ACLOCNQ_el ( iel, 1:NDFPEM_el, nq          )
              AKLOCNQ0  ( 1:NDFPEM_el, nq        ) = subbody(nb_el)%AKLOCNQ_el ( iel, 1:NDFPEM_el, nq          )
            enddo


!------------ LOC_LOADS_el
            do    ndof = 1, 6
                  i    = ndof
                  iloc = subbody(nb_el)%NDFPNACC_el(innb_el-1) + ndof
               do mnod = 1,NNPE_el
                  mn   = subbody(nb_el)%NODMB(iel,mnod)

                  do mdof = 1,subbody(nb_el)%NDFPN_el(mnod)
                     j    =                          subbody(nb_el)%NDFPNACC_el (mnod-1)+mdof
                     jj   = ACCU%NDFTACC_el(nb_el-1)+subbody(nb_el)%NDFPNBACC_el(mn  -1)+mdof

                     AMB2B (i,j) = AMB2B(i,j) + AMLOC0  (iloc,j)
                     ACB2B (i,j) = ACB2B(i,j) + ACLOC0  (iloc,j)
                     AKB2B (i,j) = AKB2B(i,j) + AKLOC0  (iloc,j)

                     ACLOC0(iloc,j) =   COEFM * AMLOC0  (iloc,j) &
                                       +COEFK * AKLOC0  (iloc,j)

                     ACB2B (i,j) = ACB2B(i,j) + ACLOC0  (iloc,j)
                     AFB2B (i  ) = AFB2B(i  ) - ACLOC0  (iloc,j)*UT1_el(jj)

                  enddo!mdof
               enddo!mnod

               do nq=1,NQS_ja
                  j=nq

                  AMB2BNQ(i,j) = AMB2BNQ(i,j) + AMLOCNQ0 (iloc,j)
                  ACB2BNQ(i,j) = ACB2BNQ(i,j) + ACLOCNQ0 (iloc,j)
                  AKB2BNQ(i,j) = AKB2BNQ(i,j) + AKLOCNQ0 (iloc,j)
               enddo!nq

                  AFB2B  (i  ) = AFB2B  (i  ) + AFLOC0   (iloc  )
            enddo!ndof


!------------ LOADS_TRANS_el
            ibcon    = boundco (nb_el)%nbcon (nbc)
           if (ibcon == 0) cycle
            ibsubcon = 1
            isbcon   = body    (ibcon)%NBODTGLB_el(ibsubcon)
            iecon    = boundco (nb_el)%nelcon(nbc)
            inibcon  = boundco (nb_el)%nodcon(nbc)

            do ndof=1,6
               indx   = boundco (nb_el)%indx (nbc,ndof)

               if (indx  == 0) cycle

               do kdof=1,6

                  ATRANS (1:6) = boundco (nb_el)%Atrans2 (nbc,kdof,1:6) !Atrans2_el (nb_el,nbc,kdof,1:6)

                  do j              = 1, subbody(nb_el)%NDFPE_el
                     AMB2BR(kdof,j) = ATRANS(ndof)*AMB2B(ndof,j)
                     ACB2BR(kdof,j) = ATRANS(ndof)*ACB2B(ndof,j)
                     AKB2BR(kdof,j) = ATRANS(ndof)*AKB2B(ndof,j)
                  enddo

                  do j                = 1, NQS_ja
                     AMB2BNQR(kdof,j) = ATRANS(ndof)*AMB2BNQ(ndof,j)
                     ACB2BNQR(kdof,j) = ATRANS(ndof)*ACB2BNQ(ndof,j)
                     AKB2BNQR(kdof,j) = ATRANS(ndof)*AKB2BNQ(ndof,j)
                  enddo
                     AFB2BR  (kdof  ) = ATRANS(ndof)*AFB2B  (ndof  )

!------------------ ASSEMBLY_MATRIX_el: Assembly of the communication loads into the global Matrices
                  nn   = subbody(isbcon)%NODMB(iecon,inibcon)
                  iloc = kdof
                  i    = ACCU%NDFTACC_el(isbcon-1)+subbody(isbcon)%NDFPNBACC_el(nn-1)+kdof

                  do mnod = 1,NNPE_el ! For all nodes per element
                     mn   =  subbody(nb_el)%NODMB(iel,mnod)

                     do mdof = 1,subbody(nb_el)%NDFPN_el(mnod) ! For all d.f. per node
                         jloc =                          subbody(nb_el)%NDFPNACC_el (mnod-1)+mdof
                         j    = ACCU%NDFTACC_el(nb_el-1)+subbody(nb_el)%NDFPNBACC_el(mn  -1)+mdof

                         AM_el(i,j) =  AM_el(i,j) + AMB2BR (iloc,jloc)
                         AC_el(i,j) =  AC_el(i,j) + ACB2BR (iloc,jloc)
                         AK_el(i,j) =  AK_el(i,j) + AKB2BR (iloc,jloc)
                     enddo!mdof
                  enddo!mnod

                  do nq   = 1, NQS_ja ! For all q's
                     jloc = nq
                     j    = NDFBT_el+nq
                  
                     AM_el(i,j) =  AM_el(i,j) + AMB2BNQR (iloc,jloc)
                     AC_el(i,j) =  AC_el(i,j) + ACB2BNQR (iloc,jloc)
                     AK_el(i,j) =  AK_el(i,j) + AKB2BNQR (iloc,jloc)
                  enddo!nq

                  AQ_el(i) =  AQ_el(i) + AFB2BR (iloc)
               enddo!kdof
            enddo!ndof

         enddo!nbc

      enddo!nbsub
   enddo!nbod


 END Subroutine BOD2BOD_LOADSJack
!--------------------------------------------------------------------------------
!
!  Subrouine : BOD2BOD_LOADS_WT_JA    ------------------
!
!  Loads Communication between WT and Jacket if ICASE_el == 2
!
!----------------------------------------------------------------------
 Subroutine BOD2BOD_LOADS_WT_JA (nbodA)
!----------------------------------------------------------------------

 Use Cbeam

   implicit None

   real(8) :: Aba      (3,3        )
   real(8) :: Aba0     (3,3,NQS    )
   real(8) :: Rba      (3          )
   real(8) :: Rba0     (3  ,NQS    )
   real(8) :: AMB2B    (6,NDFPEM_el)
   real(8) :: ACB2B    (6,NDFPEM_el)
   real(8) :: AKB2B    (6,NDFPEM_el)
   real(8) :: AMB2BNQ  (6,NQS      )
   real(8) :: ACB2BNQ  (6,NQS      )
   real(8) :: AKB2BNQ  (6,NQS      )
   real(8) :: AFB2B    (6          )
   integer :: nbod, i, j, NDFPE, iq, nq
   integer :: nnsubB, nbsubB, nnB_el, nbB_el, nelB
   integer :: nbodA , nbsubA, nnA_el, nbA_el, nelA, nodA, nn


   if ( (NBODBTWT_el == 0).or.(ICASE_el /= 2) ) return


   nbod = NBODBTWT_el

!------ Body (B)
      nbsubB  = 1
      nnsubB  = nbsubB
      nelB    = 1

      nnB_el  = body(nbod)%NNODTGLB_el(nnsubB)
      nbB_el  = body(nbod)%NBODTGLB_el(nbsubB)

      NDFPE   = subbody(nbB_el)%NDFPE_el


!------ Body (A)
      nbsubA = body(nbodA)%NBODSUB_el
      nnA_el = body(nbodA)%NNODTGLB_el(nbsubA)
      nbA_el = body(nbodA)%NBODTGLB_el(nbsubA)

      nelA   = subbody(nbA_el)%NTEB_el
      nodA   = NNPE_el
      nn     = subbody(nbA_el)%NODMB(nelA,nodA)

      i      = ACCU%NDFTACC_el(nbA_el-1)+subbody(nbA_el)%NDFPNBACC_el(nn-1)!+1:6
      j      = ACCU%NDFTACC_el(nbB_el-1)!1:NDFPE !+subbody(nb_el)%NDFPNBACC_el(mn-1)+mdof


!------ Aba, Aba0, Rba, Rba0
      Aba(1:3,1:3) = matmul (transf_mat(nnA_el)%AT_el(1:3,1:3),transf_mat(nnB_el)%A_el(1:3,1:3))

      do iq = 1, QSnode(nnB_el)%Tot_el(2)
         nq = QSnode(nnB_el)%IORD_el(iq)

         Aba0 ( 1:3, 1:3, nq ) = matmul( transf_mat(nnA_el)%AT_el(1:3,1:3), transf_mat(nnB_el)%A0_el(iq,1:3,1:3) )
      enddo!iq 


      call LOC_LOADS_el        ( AMB2B  , ACB2B  , AKB2B  , AFB2B,&
                                 AMB2BNQ, ACB2BNQ, AKB2BNQ       ,&
                                 nbB_el  ,nnB_el  ,nelB   , 1       )!nodB


      call LOADS_TRANS_el      ( AMB2B  , ACB2B  , AKB2B  , AFB2B,&
                                 AMB2BNQ, ACB2BNQ, AKB2BNQ       ,&
                                 Aba0   , Aba    , Rba0   , Rba  ,&
                                          nnB_el , 0             ,& !IcalcR
                                 NDFPE  , NDFPEM_el                 )


      call ASSEMBLY_MATRIX_el  ( AMB2B  , ACB2B  , AKB2B  , AFB2B      ,&
                                 AMB2BNQ, ACB2BNQ, AKB2BNQ, 0.d0 ,0.d0 ,&           !COEFM_el(nbod) , COEFK_el(nbod),&
                                          nnB_el , i      , j    , 2   ,&           !IQRayl:: 1:qs with Rayleigh, 2:no Rayleigh to qs
                                 6      , NDFPE  , 6      , NDFPEM_el       )
                               !IMAX    , JMAX   , Idi    , Jdi



 END Subroutine BOD2BOD_LOADS_WT_JA
!----------------------------------------------------------------------
 Subroutine QS_TOWER_WT_JA (IP, nbodA)
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   integer :: nbod,nbsub,iq,nbodA,nbsubA,nbA_el,nnA_el,nelA,nodA,mnA
   integer :: NQ_el,i,ju1,ju,IP


!-- Set tower's bottom q's equal to jacket's top deflections
!
!-- Local U of the jacket body which gives movements, transfered to Tower c.s before the 6 first q's
!
!-- (only 1 rotation by pi/2 over x axis is done)
!
!-- AT90*A_el(nnA_el)* Uloc(nbA_el,9:15) = Qtower(1:6)
!
   nbod           = NBODBTWT_el
   nbsub          = 1
   iq             = NDFBT_el + NQSP + ACCU%NQSACC_el (nbod-1) !!!+ACCU%NQSBACC_el(nbod,nbsub-1)


   nbsubA         = body(nbodA)%NBODSUB_el
   nbA_el         = body(nbodA)%NBODTGLB_el (nbsubA)
   nnA_el         = body(nbodA)%NNODTGLB_el (nbsubA)
   nelA           = subbody(nbA_el)%NTEB_el
   nodA           = NNPE_el
   mnA            = subbody(nbA_el)%NODMB (nelA, nodA)

   ju             = ACCU%NDFTACC_el(nbA_el-1) + subbody(nbA_el)%NDFPNBACC_el(mnA-1)


!- Assign the first 6 qs to the last 6 elastic dof of the jacket
!- REDUCTION enabled

      NQ_el = iq-NDFBT_el!NQSP + ACCU%NQSACC_el (nbod-1) + ACCU%NQSBACC_el(nbod, nnsub-1)
   do IQ    = 1, 6
      NQ_el = NQ_el + 1
      i     = NDFBT_el + NQ_el
      ju1   = ju + IMDOF_el(IQ)

      INDSYSQ(NQ_el)  = 0
      if (IEIG>0.or.ICASE_el/=2) &
      INDSYSQ(NQ_el)  = ju1

      if (IP.eq.0) then
         AK_el(i,i  ) = 1.d0
         AK_el(i,ju1) =-1.d0
         AQ_el(i    ) = UT_el(ju1)-UT_el(i)
      else
         AM_el(i,i  ) = 1.d0
         AM_el(i,ju1) =-1.d0
         AQ_el(i    ) = UT2_el(ju1)-UT2_el(i)
      endif
   enddo !IQ


 END Subroutine QS_TOWER_WT_JA
!-----------------------------------------------------------------------
!
!-- Subroutine :BOUNDCO_el
!
!   BOUNDary COnditions
!
!
!-----------------------------------------------------------------------
 Subroutine BOUNDCOJack_el(IP)
!-----------------------------------------------------------------------

 Use Cbeam

   implicit None

   integer :: IP,nbod,nbsub,nb_el,nbc,ndof,indx,iel,inodel,inode,NDFB0,ii
   integer :: iecon,incon,inodc,NDFBC,k,kkc,ibcon,ibsubcon,isbcon
     

!--- Boundary Conditions for Movements and Rotations 
   do    nbod  = 1+NBODBTWT_el, NBODBT_el
   do    nbsub = 1, body(nbod)%NBODSUB_el
         nb_el = body(nbod)%NBODTGLB_el(nbsub)

      do nbc      = 1, boundco (nb_el)%nbcpb
         iel      = boundco (nb_el)%nel (nbc)
         inodel   = boundco (nb_el)%nod (nbc)
         ibcon    = boundco (nb_el)%nbcon (nbc)

         do ndof  = 1, 6
            indx  = boundco (nb_el)%indx (nbc,ndof)
            if (indx /= 1) cycle
        
            inode =                            subbody(nb_el)%NODMB(iel,inodel)
            NDFB0 = ACCU%NDFTACC_el(nb_el-1) + subbody(nb_el)%NDFPNBACC_el(inode-1)
            ii    = NDFB0 + ndof

            if (ibcon == 0) then
               call ZERO_DOF ( IP, ii, 0, 0 )
               cycle
            endif

         ibsubcon = 1
         isbcon   = body    (ibcon)%NBODTGLB_el(ibsubcon)
         iecon    = boundco (nb_el)%nelcon(nbc)
         incon    = boundco (nb_el)%nodcon(nbc)
            inodc =                             subbody(isbcon)%NODMB(iecon,incon)
            NDFBC = ACCU%NDFTACC_el(isbcon-1) + subbody(isbcon)%NDFPNBACC_el(inodc-1)

            AM_el (ii, 1:NDFT_el) = 0.d0;
            AC_el (ii, 1:NDFT_el) = 0.d0;
            AK_el (ii, 1:NDFT_el) = 0.d0;
            AQ_el (ii           ) = 0.d0   !not needed, because it is set below
            
            if (IP == 0) then
               AK_el(ii,ii) = 1.d0
               AQ_el(ii   ) = -UT_el(ii)
               
               do k             = 1, 6
                  kkc           = NDFBC + k
                  AK_el(ii,kkc) =           - boundco (nb_el)%Atrans1 (nbc,ndof,k)               !Atrans1_el (nb_el,nbc,ndof,k)
                  AQ_el(ii    ) = AQ_el(ii) + boundco (nb_el)%Atrans1 (nbc,ndof,k) * UT_el (kkc)
               enddo
            else
               AM_el(ii,ii) = 1.d0
               AQ_el(ii   ) = -UT2_el(ii)
               
               do k             = 1, 6
                  kkc           = NDFBC + k
                  AM_el(ii,kkc) =           - boundco (nb_el)%Atrans1 (nbc,ndof,k)               !Atrans1_el (nb_el,nbc,ndof,k)
                  AQ_el(ii    ) = AQ_el(ii) + boundco (nb_el)%Atrans1 (nbc,ndof,k) * UT2_el (kkc)
               enddo
            endif !IP
         enddo !ndof

      enddo !nbc

   enddo !nbsub
   enddo !nbod


 END Subroutine BOUNDCOJack_el
