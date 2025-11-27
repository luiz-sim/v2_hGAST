!-------------------------------------------------------------------------
! Foundation Modeling types supported [ITTYPE_FOUND_el]
!
!   0. AF  Apparently Fixed
!   1. CS  Concentrated Spring
!   2. DS  Distributed  Springs
!  (3. DS with p-y curves)
!-------------------------------------------------------------------------
 Subroutine Foundation_contrib
!-------------------------------------------------------------------------

 Use Cbeam

   implicit None

   real(8) :: AKLOC0  (NDFPEM_el,NDFPEM_el)
   real(8) :: AFLOC0  (NDFPEM_el          )
   integer :: nb_el, NDFPE, iel, nel, i, j, nn, ibod


   if (ITYPE_FOUND_el == 0) return


   do    ibod  = 1, NBODTFND_el
         nb_el = foundation(ibod)%nb
         NDFPE = subbody(nb_el)%NDFPE_el

      do iel   = 1, foundation(ibod)%nteb_fnd
         nel   = foundation(ibod)%nel (iel)
         nn    =                                subbody(nb_el)%NODMB       (nel,1)!nnod)
         i     = ACCU%NDFTACC_el(nb_el-1)   +   subbody(nb_el)%NDFPNBACC_el(nn-1 )
         j     = i

         AKLOC0(1:NDFPE, 1:NDFPE) = foundation(ibod)%Stiff (iel,1:NDFPE, 1:NDFPE)
         AFLOC0(1:NDFPE         ) = - matmul ( AKLOC0(1:NDFPE, 1:NDFPE), UT_el(j+1:j+NDFPE) )

!assembly local matrices to global
         AK_el (i+1:i+NDFPE, j+1:j+NDFPE) = AK_el(i+1:i+NDFPE, j+1:j+NDFPE) + AKLOC0(1:NDFPE, 1:NDFPE)
         AQ_el (i+1:i+NDFPE             ) = AQ_el(i+1:i+NDFPE             ) + AFLOC0(1:NDFPE         )

         if ( ITYPE_FOUND_el == 1) cycle

!store local matrices
         subbody(nb_el)%AKLOC_el(nel, 1:NDFPEM_el, 1:NDFPEM_el) = subbody(nb_el)%AKLOC_el(nel, 1:NDFPEM_el, 1:NDFPEM_el) + AKLOC0(1:NDFPEM_el, 1:NDFPEM_el)
         subbody(nb_el)%AFLOC_el(nel, 1:NDFPEM_el             ) = subbody(nb_el)%AFLOC_el(nel, 1:NDFPEM_el             ) + AFLOC0(1:NDFPEM_el             )
      enddo !iel
   enddo !ibod


 END Subroutine Foundation_contrib
!-------------------------------------------------------------------------
 Subroutine Foundation_Init
!-------------------------------------------------------------------------

 Use Cbeam
 Use Paths

   implicit None

   real(8) :: AKLOC0    (NDFPEM_el,NDFPEM_el)
   real(8) :: SHAPETAK  (NDFPEM_el,NEQPEM_el)
   real(8) :: SHAPET    (NDFPEM_el,NEQPEM_el)
   real(8) :: SHAPE     (NEQPEM_el,NDFPEM_el)
   real(8) :: AK        (NEQPEM_el,NEQPEM_el)
   real(8) :: HTA, ALLOC, HTAL, Stiff(6,6),PHIX, PHIZ
   integer :: ibod,nbod,nbsub,nb_el, NDFPE, NEQPE, iel,nel, NTSprPE, ispr, i, j, nteb_fnd, indx(6)
   integer :: ibcon, ibcapplied, nbc, inod


   if (ITYPE_FOUND_el == 0) return

                            write(10,*)
                            write(10,*) 'FOUNDATION Modeling is enabled'
   if (ITYPE_FOUND_el == 1) write(10,*) 'with Concentrated Spring Model'
   if (ITYPE_FOUND_el == 2) write(10,*) 'with Distributed Springs Model'

   open (21,file=trim(file_foundation))

   read (21,*)                                                                 ! title
   read (21,*) NBODTFND_el                                                     ! total sub-bodies at which foundation acts

   Allocate ( foundation (NBODTFND_el) )

   stiff(:,:) = 0.d0;


   do ibod  = 1, NBODTFND_el

      read(21,*) nbod, nteb_fnd
      
      nbsub = 1                                                                ! by default foundation is applied at the 1st subbody...
      nb_el = body(nbod)%NBODTGLB_el(nbsub)
      NDFPE = subbody(nb_el)%NDFPE_el
      NEQPE = subbody(nb_el)%NEQPE_el

                        foundation(ibod)%nb       = nb_el                      ! sub-body        at which foundation acts
                        foundation(ibod)%nteb_fnd = nteb_fnd                   ! number of elements at which foundation acts
             Allocate ( foundation(ibod)%nel       (nteb_fnd            ) )    ! nel of sub-body    at which foundation acts
             Allocate ( foundation(ibod)%Stiff     (nteb_fnd,NDFPE,NDFPE) )    ! Local element Stiffness matrix from foundation
!            Allocate ( foundation(ibod)%Damp      (nteb_fnd,NDFPE,NDFPE) )    ! Local element Damping   matrix from foundation)

      do iel  = 1, foundation(ibod)%nteb_fnd

         read (21,*) nel, NTSprPE !number of total spring per element

         PHIX = subbody(nb_el)%PHIX_el(nel)
         PHIZ = subbody(nb_el)%PHIZ_el(nel)

         AKLOC0 = 0.d0;

         do ispr = 1, NTSprPE

            if    (ITYPE_FOUND_el==1) then
                  read(21,*) HTA, (Stiff(1,j),j=1,6)
               do i = 2, 6
                  read(21,*)      (Stiff(i,j),j=1,6)
               enddo
            elseif(ITYPE_FOUND_el==2) then
                  read(21,*) HTA, (Stiff(j,j),j=1,6)
            endif

            ALLOC  = subbody(nb_el)%ALENG_el(nel)
            HTAL   = HTA/ALLOC

            call SSHAPEFUNC15 ( HTAL, ALLOC, PHIX, PHIZ,& 
                               SHAPE, NDFPE, NEQPE       )

            SHAPET   = transpose (SHAPE )
            AK(IMDOF_el(1:6),IMDOF_el(1:6)) = Stiff(1:6,1:6)

            SHAPETAK = matmul( SHAPET  , AK    )
            AKLOC0   = matmul( SHAPETAK, SHAPE ) + AKLOC0
         enddo !ispr

        foundation(ibod)%nel   (iel                 ) = nel                       ! nel of sub-body    at which foundation acts
        foundation(ibod)%Stiff (iel,1:NDFPE, 1:NDFPE) = AKLOC0 (1:NDFPE, 1:NDFPE) ! Local element Stiffness matrix from foundation
     enddo !iel

   enddo !ibod


!-- set the corresponding b.c. depending on the Foundation type
      read(21,*)
      read(21,*)!** boundary conditions applied to all sub-bodies at which foundation acts
      read(21,*) (indx(i),i = 1, 6)
   close(21)

   do ibod  = 1, NBODTFND_el
      nb_el = foundation(ibod)%nb

      ibcapplied = 0
      do nbc = 1, boundco (nb_el)%nbcpb

         ibcon = boundco (nb_el)%nbcon(nbc)
         if (ibcon /= 0) cycle

         iel  = boundco (nb_el)%nel (nbc)
         inod = boundco (nb_el)%nod (nbc)

         if ((iel/=1).or.(inod/=1)) cycle
         boundco (nb_el)%indx  (nbc,IMDOF_el(1:6)) = indx(1:6)
         ibcapplied = 1                                         ! the b.c are correctly modified
      enddo !nbc

      if (ibcapplied==0) then
         write(* ,*) 'Foundation error, stoping!',ibod,nb_el
         write(10,*) 'Foundation error, stoping!',ibod,nb_el
         stop
      endif
   enddo !ibod


 END Subroutine Foundation_Init
