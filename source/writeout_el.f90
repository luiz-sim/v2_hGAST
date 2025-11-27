!--- to do ---
! OK - write mode for each body
! OK - write mode: remove R0 and also output r so {r,u,v,w,thx,thy,thz}
! OK - loads loop one
! OK - subtract R0 in deform
!      deform torsion wrt init precurved angle -- no, but wrt to the local deformed state
! OK - switch for local or global output [jacket - prebend blade, rotating shaft]
!          JA,TO: global      , local
!          SH   : non-rotating, rotating
!          BLD  : global      , local 
!      set output order and inform the users and update the ppt!
!----------------------------------------------------------------------
 Subroutine OUTPUT (NTIME)
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   integer, intent(in) :: NTIME

   integer      :: NTIMEWRITE !, NTIMEWRITE0
   integer      :: itap, j
   character    :: CNUM1(3)
   character*80 :: outfil


#  ifndef    ASCII
#     define ASCII 1
#  endif

!--- manage output frequency
                      NTIMEWRITE = 1
   if (DT_el<=0.03d0) NTIMEWRITE = 2
    !  if (     NTIME < NTIMEWRITE0    ) return
       if (mod((NTIME),(NTIMEWRITE))/=0) then
          call ROTMAT_el
          return
       endif
       

   if (IAERO_el == 1) &
   call WRITEOUTA

   call write_loads
   call ROTMAT_el
   call write_deform ( UT_el)
   call write_qsdof
   call WRITEOUT_zita

   if (ITRUSS_el > 0) &
   call WRITELEM_tr(NTIME)

!!!call Writeout_MORISON
!qq                            outfil =  'geometry_'
!qqcall WRITEOUT_GLOB  (UT_el, outfil, NTIME, 5)

   if (ICASE_el == 3) &
   call floater_simulator( 3, 99 ) !iter=99 dummy for output


!------ Writeout Buoyancy taps info
   do itap   = 1, NBuoyaTap_el
      call INT_2_CHAR ( 3, CNUM1, itap )
#    if   ASCII == 1
      outfil = 'buoy_tap'//CNUM1(1)//CNUM1(2)//CNUM1(3)//'.dat'
      open (201, file=outfil,access='append')
      write(201,100)                           &
#    elif ASCII == 0
      outfil = 'buoy_tap'//CNUM1(1)//CNUM1(2)//CNUM1(3)//'.bin'
      open (201, file=outfil,access='append',form='UNFORMATTED')
      write(201    )                           &
#    endif
       sngl(TIME                            ), & ! 1
       sngl(buoyancy_tap(itap)%Output_res(1)), & ! 2  !Fhyd
       sngl(buoyancy_tap(itap)%Output_res(2)), & ! 3  !Fdyn
       sngl(buoyancy_tap(itap)%Output_res(3)), & ! 4  !Frel
       sngl(buoyancy_tap(itap)%Output_res(4)), & ! 5  !Fdrag
       sngl(buoyancy_tap(itap)%Output_res(5)), & ! 6  !Finer
       sngl(buoyancy_tap(itap)%Output_res(6)), & ! 7  !Faddmass
       sngl(buoyancy_tap(itap)%Output_res(7))    ! 8  !Z0g
      close(201)
   enddo !itap


!---- Write Total Buoyancy
   if ((ICASE_el == 2).or.(ICASE_el == 4)) then

      BuoyancyT_el (1:6) = BuoyancyTtap_el (1:6) + &
                           BuoyancyTside_el(1:6)

#    if   ASCII == 1
      open  (1,file='buoy_tot.dat',access='append')
      write (1,100)                       &
#    elif ASCII == 0
      open  (1,file='buoy_tot.bin',access='append',form='UNFORMATTED')
      write (1    )                       &
#    endif
       sngl(TIME)                        ,&
     ( sngl(BuoyancyT_el    (j)), j=1,6 ),&
     ( sngl(BuoyancyTtap_el (j)), j=1,6 ),&
     ( sngl(BuoyancyTside_el(j)), j=1,6 )
      close(1)
   endif

 100  format (150en24.8e3)


 END Subroutine OUTPUT
!----------------------------------------------------------------------
 Subroutine write_deform ( UT )
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   real(8)      :: UT  (*)
   real(8)      :: AL2G(3,3), AL  (3,3)
   real(8)      :: AT0G(3,3), AT0L(3,3)
   real(8)      :: X0  (6)  , Xg  (6)  , Xl  (6)
   real(8)      :: XL0p(6)  , XL0 (6)
   real(8)      :: PSI0   , PSIB
   integer      :: nbod, nbsub, nel, nn0_el, numnod, nb_el, nn_el, nnod, nnod0, nn, i, j
   character    :: CNUM1(3), CNUM2(2)
   character*80 :: outfil
   integer      :: IwriteG, IwriteL  !(1:write global/local loads, 0:not)


   do nbod         = 1, NBODBT_el
      IwriteG      = IWRITE_el (min(nbod,NBODBTWT_el+1),1)  ! flag for writing global/local signals (loads,deform) [NBODBT_el+1,2] 1:global, local
      IwriteL      = IWRITE_el (min(nbod,NBODBTWT_el+1),2)  ! flag for writing global/local signals (loads,deform) [NBODBT_el+1,2] 1:global, local

      if (nbod<=NBLADE_el) then; PSI0 = -UT_el(NDFBT_el+NQSW) - UThub_el + PHI0(nbod)
      else                     ; PSI0 = -UT_el(NDFBT_el+NQSW) + PHI0(1)
      endif
                                 PSIB = dmod (PSI0,PI2)*R2D

      call INT_2_CHAR ( 3, CNUM1, nbod )


!--------- Define global and local c.s. per body
      if     (nbod <= NBLADE_el) then
!--------- wrt blade root after cone and pitch, before pre-curved angle
         AT0G(1:3,1:3) = ATPWrG_el(nbod,1:3,1:3) !Blade c.s. after cone, before pitch 
         AT0L(1:3,1:3) = ATPWr_el (nbod,1:3,1:3) !Blade c.s. after cone, after  pitch, before pre-curved angles
!--------- wrt global c.s. for the VAWT
         if (IAPPL_el == 1) &
         AT0G(1:3,1:3) = ATsftNoR_el (1:3,1:3)
!!       call DIAGO (3, AT0G)
      elseif (nbod == NSHAFT_el) then
!--------- wrt the non-rotating shaft c.s.
         AT0G(1:3,1:3) = ATsftNoR_el (1:3,1:3)
         nn0_el        = body(nbod)%NNODTGLB_el(1)
         AT0L(1:3,1:3) = transf_mat(nn0_el)%AT_el(1:3,1:3)
      else
!--------- wrt global c.s. for jacket and tower
         call DIAGO (3, AT0G)
         nn0_el        = body(nbod)%NNODTGLB_el(1)
         AT0L(1:3,1:3) = transf_mat(nn0_el)%AT_el(1:3,1:3)
      endif


            numnod     = 0
            nnod0      = 1        !1 only for 1st nel of 1st nbsub for each body, otherwise NNPE_el
            XL0 (1:6)  = 0.d0
            XL0p(1:6)  = 0.d0
            AL2G       = matmul (transpose(AT0L),AT0G)
      do    nbsub      = 1, body(nbod)%NBODSUB_el
            nb_el      = body(nbod)%NBODTGLB_el(nbsub)
            nn_el      = body(nbod)%NNODTGLB_el(nbsub)
            AL(1:3,1:3)= matmul (AT0L (1:3,1:3),transf_mat(nn_el)%A_el(1:3,1:3))
            XL0p (1:6) = XL0p(1:6) + XL0(1:6)
         do nel        = 1, subbody(nb_el)%NTEB_el
         do nnod       = nnod0, NNPE_el, NNPE_el-1
            numnod     = numnod + 1
            nn         =                            subbody(nb_el)%NODMB       (nel,nnod)
            i          = ACCU%NDFTACC_el(nb_el-1) + subbody(nb_el)%NDFPNBACC_el(nn-1)

            do j       = 1, 6
              X0(j)    = UT (i+IMDOF_el(j))
            enddo

            if (nbod<=NBLADE_el) X0(2) = X0(2) + subbody(nb_el)%HTA_el(nel+nnod/NNPE_el)

            XL0(1:3)   = matmul( AL(1:3,1:3),X0(1:3) )
            XL0(4:6)   = matmul( AL(1:3,1:3),X0(4:6) )

            call INT_2_CHAR ( 2, CNUM2, numnod )

           if (nbod<=NBLADE_el) then
            Xl   (1:3) =  XL0(1:3) + XL0p(1:3) - subbody(nb_el)%Rinit (nel+nnod/NNPE_el,1:3)
           !Xl   (4:6) = (XL0(4:6) + XL0p(4:6) - subbody(nb_el)%Rinit (nel+nnod/NNPE_el,4:6))*R2D
            Xl   (4:6) = (XL0(4:6) + XL0p(4:6))*R2D
           else
            Xl   (1:3) =  XL0(1:3) + XL0p(1:3)
            Xl   (4:6) = (XL0(4:6) + XL0p(4:6))*R2D
           endif

           if (IwriteG==1) then
            Xg   (1:3) =  matmul(AL2G(1:3,1:3),Xl(1:3))
            Xg   (4:6) =  matmul(AL2G(1:3,1:3),Xl(4:6))
           endif !IwriteG

#       if   ASCII == 1
           if (IwriteG==1) then
            outfil     = 'deformG'//CNUM1(1)//CNUM1(2)//CNUM1(3)//'_'//CNUM2(1)//CNUM2(2)//'.dat'
            open  (1,file=outfil, access='append')
            write (1,100)            &
                  (TIME)            ,& ! 1
                  (PSIB)            ,& ! 2
             (    ( Xg  (j) ),j=1,6)   ! 3-8
            close (1)
           endif !IwriteG

           if (IwriteL==1) then
            outfil     = 'deform'//CNUM1(1)//CNUM1(2)//CNUM1(3)//'_'//CNUM2(1)//CNUM2(2)//'.dat'
            open  (1,file=outfil, access='append')
            write (1,100)            &
                  (TIME)            ,& ! 1
                  (PSIB)            ,& ! 2
             (    ( Xl  (j) ),j=1,6)   ! 3-8
            close (1)
           endif !IwriteL

#       elif ASCII == 0

           if (IwriteG==1) then
            outfil     = 'deformG'//CNUM1(1)//CNUM1(2)//CNUM1(3)//'_'//CNUM2(1)//CNUM2(2)//'.bin'
            open  (1,file=outfil, access='append' ,form='UNFORMATTED')
            write (1    )            &
              sngl(TIME)            ,& ! 1
              sngl(PSIB)            ,& ! 2
             (sngl( Xg  (j) ),j=1,6)   ! 3-8
            close (1)
           endif !IwriteG

           if (IwriteL==1) then
            outfil     = 'deform'//CNUM1(1)//CNUM1(2)//CNUM1(3)//'_'//CNUM2(1)//CNUM2(2)//'.bin'
            open  (1,file=outfil, access='append' ,form='UNFORMATTED')
            write (1    )            &
              sngl(TIME)            ,& ! 1
              sngl(PSIB)            ,& ! 2
             (sngl( Xl  (j) ),j=1,6)   ! 3-8
            close (1)
           endif !IwriteL

#       endif
         enddo !nnod
            nnod0      = NNPE_el
         enddo !nel
      enddo !nbsub
   enddo !nbod

 100  format (150en24.8e3)


 END Subroutine write_deform
!----------------------------------------------------------------------
 Subroutine write_loads
!----------------------------------------------------------------------
!
!   FL
!   FG=[ATbodyG.AG].FL blade
!
!           IwriteG=1      IwriteL=1
!           Global         Local      
!   blade : rotor disk   - local-root 
!   shaft : non-rotating - rotating   
!   tower : global c.s.  - local      
!   jacket: global c.s.  - local      
!                                     
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   real(8)      :: AG     (3,3), AL     (3,3)
   real(8)      :: AT0G   (3,3), AT0L   (3,3)
   real(8)      :: AFB2BG (6  ), AFB2BL (6  )
   real(8)      :: PSI0, PSIB, dir
   integer      :: nbod, IEL, I1, I2
   integer      :: nbsub, nn_el , nb_el , nel, nod, nod_end, nn0_el, j
   character    :: CNUM1(3), CNUM2(2)
   character*80 :: outfil
   real(8)      :: NtokN
   integer      :: IwriteG, IwriteL  !(1:write global/local loads, 0:not)


                      NtokN = 1000.d0
   if (IAPPL_el == 1) NtokN =    1.d0

   do nbod         = 1, NBODBT_el
      IwriteG      = IWRITE_el (min(nbod,NBODBTWT_el+1),1)  ! flag for writing global/local signals (loads,deform) [NBODBT_el+1,2] 1:global, local
      IwriteL      = IWRITE_el (min(nbod,NBODBTWT_el+1),2)  ! flag for writing global/local signals (loads,deform) [NBODBT_el+1,2] 1:global, local

      if (nbod<=NBLADE_el) then; PSI0 = -UT_el(NDFBT_el+NQSW) - UThub_el + PHI0(nbod)
      else                     ; PSI0 = -UT_el(NDFBT_el+NQSW) + PHI0(1)
      endif
                                 PSIB = dmod (PSI0,PI2)*R2D

      call INT_2_CHAR ( 3, CNUM1, nbod )


!--------- Define global and local c.s. per body
      if     (nbod <= NBLADE_el) then
!--------- wrt blade root after cone and pitch, before pre-curved angle
         AT0G(1:3,1:3) = ATPWrG_el(nbod,1:3,1:3) !Blade c.s. after cone, before pitch 
         AT0L(1:3,1:3) = ATPWr_el (nbod,1:3,1:3) !Blade c.s. after cone, after  pitch, before pre-curved angles
!--------- wrt global c.s. for the VAWT
         if (IAPPL_el == 1) &
         AT0G(1:3,1:3) = ATsftNoR_el (1:3,1:3)
!!       call DIAGO (3, AT0G)
      elseif (nbod == NSHAFT_el) then
!--------- wrt the non-rotating shaft c.s.
         AT0G(1:3,1:3) = ATsftNoR_el (1:3,1:3)
         nn0_el        = body(nbod)%NNODTGLB_el(1)
         AT0L(1:3,1:3) = transf_mat(nn0_el)%AT_el(1:3,1:3)
      else
!--------- wrt global c.s. for jacket and tower
         call DIAGO (3, AT0G)
         nn0_el        = body(nbod)%NNODTGLB_el(1)
         AT0L(1:3,1:3) = transf_mat(nn0_el)%AT_el(1:3,1:3)
      endif


            IEL            = 0
      do    nbsub          = 1, body(nbod)%NBODSUB_el
            nb_el          = body(nbod)%NBODTGLB_el(nbsub)
            nn_el          = body(nbod)%NNODTGLB_el(nbsub)
            AG(1:3,1:3)    = matmul (AT0G (1:3,1:3),transf_mat(nn_el)%A_el(1:3,1:3))
            AL(1:3,1:3)    = matmul (AT0L (1:3,1:3),transf_mat(nn_el)%A_el(1:3,1:3))
         do nel            = 1, subbody(nb_el)%NTEB_el
            nod_end        = 1
            if (nbsub ==    body(nbod)%NBODSUB_el .and.   &
                nel   == subbody(nb_el)%NTEB_el         ) &
            nod_end        = NNPE_el        !only for last nel of last nbsub for each body
            dir            = 1.d0           !1. for nod=1, -1 for nod=NNPE_el
         do nod            = 1, nod_end, NNPE_el-1
            IEL            = IEL + 1
            I1             = subbody(nb_el)%NDFPNACC_el ( nod-1      ) + 1
            I2             = subbody(nb_el)%NDFPNACC_el ( nod        )
            AFB2BG ( 1:6 ) = subbody(nb_el)%AFLOC_el    ( nel, I1:I2 )*dir
            AFB2BL ( 1:6 ) = subbody(nb_el)%AFLOC_el    ( nel, I1:I2 )*dir
            dir            =-1.d0

!------------ Rotate LOADS with respect to body c.s. (0)
            call LOADS_TRANS0_NEW ( AFB2BG, AG, 1, 1, 1 )
            call LOADS_TRANS0_NEW ( AFB2BL, AL, 1, 1, 1 )
            call INT_2_CHAR       ( 2, CNUM2, IEL )

#       if   ASCII == 1
           if (IwriteG==1) then
            outfil = 'loadsG'//CNUM1(1)//CNUM1(2)//CNUM1(3)//'_'//CNUM2(1)//CNUM2(2)//'.dat'
            open (3,file=outfil, access='append')
            write(3,100)                                           &
                 (TIME)                                           ,&
                 (PSIB)                                           ,&
            (          (AFB2BG (IMDOF_el(j))/NtokN),j=1,6)        ,&  !wrt chosen global c.s.
                 (sqrt( AFB2BG (IMDOF_el(4))**2 +                  &  !combined bending moment 
                        AFB2BG (IMDOF_el(6))**2   )/NtokN)
            close(3)
           endif !IwriteG

           if (IwriteL==1) then
            outfil = 'loads'//CNUM1(1)//CNUM1(2)//CNUM1(3)//'_'//CNUM2(1)//CNUM2(2)//'.dat'
            open (3,file=outfil, access='append')
            write(3,100)                                           &
                 (TIME)                                           ,&
                 (PSIB)                                           ,&
            (          (AFB2BL (IMDOF_el(j))/NtokN),j=1,6)        ,&  !wrt chosen local c.s.
                 (sqrt( AFB2BL (IMDOF_el(4))**2 +                  &  !combined bending moment 
                        AFB2BL (IMDOF_el(6))**2   )/NtokN)
            close(3)
           endif !IwriteL

#       elif ASCII == 0

           if (IwriteG==1) then
            outfil = 'loadsG'//CNUM1(1)//CNUM1(2)//CNUM1(3)//'_'//CNUM2(1)//CNUM2(2)//'.bin'
            open (3,file=outfil, access='append' ,form='UNFORMATTED')
            write(3    )                                           &
             sngl(TIME)                                           ,&
             sngl(PSIB)                                           ,&
            (sngl      (AFB2BG (IMDOF_el(j))/NtokN),j=1,6)        ,&  !wrt chosen global c.s.
             sngl(sqrt( AFB2BG (IMDOF_el(4))**2 +                  &  !combined bending moment 
                        AFB2BG (IMDOF_el(6))**2   )/NtokN)
            close(3)
           endif !IwriteG

           if (IwriteL==1) then
            outfil = 'loads'//CNUM1(1)//CNUM1(2)//CNUM1(3)//'_'//CNUM2(1)//CNUM2(2)//'.bin'
            open (3,file=outfil, access='append' ,form='UNFORMATTED')
            write(3    )                                           &
             sngl(TIME)                                           ,&
             sngl(PSIB)                                           ,&
            (sngl      (AFB2BL(IMDOF_el(j))/NtokN),j=1,6)         ,&  !wrt chosen local c.s.
             sngl(sqrt( AFB2BL (IMDOF_el(4))**2 +                  &  !combined bending moment 
                        AFB2BL (IMDOF_el(6))**2   )/NtokN)
            close(3)
           endif !IwriteL

#          endif
         enddo !nod
         enddo !nel
      enddo !nbsub
   enddo !nbod

 100  format (150f24.8)
!100  format (150en28.17e3)


 END Subroutine write_loads
!----------------------------------------------------------------------
 Subroutine WRITEOUT_GLOB ( UT, outfil0, idigit, MAXdig )
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   real(8)      :: UT (*)
   real(8)      :: RGlob(3),AGlob(3,3),Y(3),Rdef(3)
   integer      :: nbod, nbsub, nel, nb_el, nn_el, nnod, nn, i, idigit, MAXdig
   character    :: CNUM(MAXdig)
   character*80 :: outfil,outfil0


   call INT_2_CHAR ( MAXdig, CNUM, idigit )

      outfil      = outfil0
   do i           = 1, MAXdig
      outfil      =  trim(outfil)//CNUM(i) !   //CNUM(1:MAXdig)//'.dat'
   enddo
      outfil      =  trim(outfil)//'.dat'

   open (2,file=outfil)


!--- write the instantaneous deflected global geometry
   do       nbod  = 1, NBODBT_el
      do    nbsub = 1, body(nbod)%NBODSUB_el
            nb_el = body(nbod)%NBODTGLB_el(nbsub)
            nn_el = body(nbod)%NNODTGLB_el(nbsub)

            AGlob = transf_mat(nn_el)%A_el(1:3,1:3)
            RGlob = transf_mat(nn_el)%R_el(1:3    )

         do nel   = 1, subbody(nb_el)%NTEB_el
            nnod  = 1
            nn    =                            subbody(nb_el)%NODMB       (nel,nnod)
            i     = ACCU%NDFTACC_el(nb_el-1) + subbody(nb_el)%NDFPNBACC_el(nn-1)

            Y(1)  = UT   (i+IMDOF_el(1))
            Y(2)  = UT   (i+IMDOF_el(2)) + subbody(nb_el)%HTA_el(nel)
            Y(3)  = UT   (i+IMDOF_el(3))

            Rdef  = RGlob+matmul(AGlob,Y)

            write (2,100) Rdef(1),Rdef(2),Rdef(3)
         enddo !nel
      enddo !nbsub

            nbsub = body(nbod)%NBODSUB_el
            nb_el = body(nbod)%NBODTGLB_el(nbsub)
            nn_el = body(nbod)%NNODTGLB_el(nbsub)
            nel   = subbody(nb_el)%NTEB_el
            nnod  = NNPE_el
            nn    =                            subbody(nb_el)%NODMB       (nel,nnod)
            i     = ACCU%NDFTACC_el(nb_el-1) + subbody(nb_el)%NDFPNBACC_el(nn-1)

            Y(1)  = UT   (i+IMDOF_el(1))
            Y(2)  = UT   (i+IMDOF_el(2)) + subbody(nb_el)%HTA_el(nel+1)
            Y(3)  = UT   (i+IMDOF_el(3))

            Rdef  = RGlob+matmul(AGlob,Y)

            write (2,100) Rdef(1),Rdef(2),Rdef(3)
            write (2,*)
            write (2,*)
   enddo !nbod

   close (2)

 100  format (150en24.8e3)


 END Subroutine WRITEOUT_GLOB
!----------------------------------------------------------------------
 Subroutine WRITEOUT_LOC ( UT, outfil0, idigit, MAXdig )
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   real(8)      :: UT  (*)
   real(8)      :: X0  (3)  , T0  (3)  , T0p (3)
   real(8)      :: A   (3,3), AT0 (3,3), X0p (3), AT0A(3,3), Xg(6), y0
   integer      :: nbod, nbsub, nel, nn0_el, nb_el, nn_el, nnod, nnod0, nn, i, j
   integer      :: idigit, MAXdig
   character    :: CNUM(MAXdig), CNUM2(3)
   character*80 :: outfil,outfil0,outfil1


   call INT_2_CHAR ( MAXdig, CNUM, idigit )

      outfil1 =  outfil0
   do i=1,MAXdig
      outfil1 =  trim(outfil1)//CNUM(i) !   //CNUM(1:MAXdig)//'.dat'
   enddo


   do       nbod       = 1, NBODBT_el

      call INT_2_CHAR ( 3, CNUM2, nbod )
            outfil     =  trim(outfil1)//'_'//CNUM2(1)//CNUM2(2)//CNUM2(3)//'.dat'
      open (1,file=outfil)

      if     (nbod <= NBLADE_el) then
!--------- wrt blade root after pitch, before pre-curved angle
            AT0(1:3,1:3) = ATPWr_el(nbod,1:3,1:3)
      elseif (nbod <=  NTOWER_el) then
!--------- wrt body start local c.s. [shaft, tower and floating jacket]
            nn0_el       = body(nbod)%NNODTGLB_el(1)
            AT0(1:3,1:3) = transf_mat(nn0_el)%AT_el(1:3,1:3)
      else
!--------- wrt global c.s. [bottom mounted jacket and its tower]
         call DIAGO (3, AT0)
      endif

!--------- wrt global c.s. for the VAWT
      if (IAPPL_el == 1.and.nbod <= NBLADE_el) call DIAGO (3, AT0)

            nnod0      = 1        !1 only for 1st nel of 1st nbsub for each body, otherwise NNPE_el
            X0 (1:3)   = 0.d0
            X0p(1:3)   = 0.d0
            T0 (1:3)   = 0.d0
            T0p(1:3)   = 0.d0

      do    nbsub      = 1, body(nbod)%NBODSUB_el
            nb_el      = body(nbod)%NBODTGLB_el(nbsub)
            nn_el      = body(nbod)%NNODTGLB_el(nbsub)
            A(1:3,1:3) = transf_mat(nn_el)%A_el(1:3,1:3)
            AT0A       = matmul( AT0, A )
            X0p(1:3)   = X0p(1:3) + X0(1:3)
            T0p(1:3)   = T0p(1:3) + T0(1:3)
         do nel        = 1, subbody(nb_el)%NTEB_el
         do nnod       = nnod0, NNPE_el, NNPE_el-1
            nn         =                            subbody(nb_el)%NODMB       (nel,nnod)
            i          = ACCU%NDFTACC_el(nb_el-1) + subbody(nb_el)%NDFPNBACC_el(nn-1)

            do j       = 1, 3
              X0(j)    = UT (i+IMDOF_el(j  ))
              T0(j)    = UT (i+IMDOF_el(j+3))
            enddo

            if (nbod<=NBLADE_el) X0(2) = X0(2) + subbody(nb_el)%HTA_el(nel+nnod/NNPE_el)

            X0         = matmul( AT0A,X0 )
            T0         = matmul( AT0A,T0 )
           if (nbod<=NBLADE_el) then
            y0         = subbody(nb_el)%Rinit (nel+nnod/NNPE_el,2)
            Xg   (1:3) =  X0(1:3) + X0p(1:3) - subbody(nb_el)%Rinit (nel+nnod/NNPE_el,1:3)
!           Xg   (4:6) = (T0(1:3) + T0p(1:3) - subbody(nb_el)%Rinit (nel+nnod/NNPE_el,4:6))*R2D
            Xg   (4:6) = (T0(1:3) + T0p(1:3))*R2D
           else
            y0         = subbody(nb_el)%HTA_el(nel+nnod/NNPE_el)
            Xg   (1:3) =  X0(1:3) + X0p(1:3)
            Xg   (4:6) = (T0(1:3) + T0p(1:3))*R2D
           endif

            write (1,100)    &
              y0            ,&
             (Xg (j), j=1,6)
         enddo !nnod
            nnod0      = NNPE_el
         enddo !nel
      enddo !nbsub

      write (1,*)
   enddo !nbod

   close (1)

 100  format (150en24.8e3)


 END Subroutine WRITEOUT_LOC
!----------------------------------------------------------------------
 Subroutine write_qsdof
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   integer :: nbod, i,j,k,nb, m(NBLADE_el), nb_el
   integer :: nf, i1,i2,i3,i4
   real(8) :: power
!  real(8), dimension(3) :: RL, RG, DRG, DDRG


!-- WRITE_QS: write NQSP ( Qs hydro and rot )
 
   if (NQS > 0) then


         do nb    = 1, NBLADE_el
            m(nb) = NDFBT_el + NQSP + 5 + ACCU%NQSACC_el(nb-1)
         enddo
            i     = NDFBT_el + NQSW
            k     = NDFBT_el + NQSY

         if (NBODBTWT_el>NBLADE_el) then
            nb_el = body(NSHAFT_el)%NBODTGLB_el(1)                        !nbsub=1
            power = subbody(nb_el) %AFLOC_el   (1,IMDOF_el(5))*UT1_el(i)  !nel  =1
         else
            power = 0.d0
         endif

#if   ASCII == 1
      open  (1,file='qsdof.dat',access='append')
      write (1,100)      (TIME)                                 ,&
                         (-UT_el  (i))                          ,&
                         (-UT1_el (i))                          ,& !shaft speed
                         (-UT2_el (i))                          ,&
                    (    (-UT_el  (m(nb)) * R2D)                ,& !Pitch of Blades       !5, 8,11
                         (-UT1_el (m(nb)) * R2D)                ,&                        !6, 9,12
                         (-UT2_el (m(nb)) * R2D),nb=1,NBLADE_el),&                        !7,10,13
                         ( UT_el  (k)     * R2D)                ,& !Yaw free movement     !14 (if nblad=3)
                         ( UT1_el (k)     * R2D)                ,&
                         ( UT2_el (k)     * R2D)                ,&
                           TGenLSS_el/1000.d0                   ,& !17: GenTorque wrt LSS [kN]
                           power     /1000.d0                      !18: Mechanical power at the gen side [kW]
#elif ASCII == 0
      open  (1,file='qsdof.bin',access='append',form='UNFORMATTED')
      write (1    )  sngl(TIME)                                 ,&
                     sngl(-UT_el  (i))                          ,&
                     sngl(-UT1_el (i))                          ,& !shaft speed
                     sngl(-UT2_el (i))                          ,&                                            
                    (sngl(-UT_el  (m(nb)) * R2D)                ,& !Pitch of Blades       !5, 8,11
                     sngl(-UT1_el (m(nb)) * R2D)                ,&                        !6, 9,12
                     sngl(-UT2_el (m(nb)) * R2D),nb=1,NBLADE_el),&                        !7,10,13
                     sngl( UT_el  (k)     * R2D)                ,& !Yaw free movement     !14 (if nblad=3)
                     sngl( UT1_el (k)     * R2D)                ,&                                            
                     sngl( UT2_el (k)     * R2D)                ,&                                            
                     sngl( TGenLSS_el/1000.d0  )                ,& !17: GenTorque wrt LSS [kN]
                     sngl( power     /1000.d0  )                   !18: Mechanical power at the gen side [kW]
#endif
      close(1)
   endif



   if ( ICASE_el >= 3 ) then

      do nf = 1, NFLOATER_el
         i1 = (nf-1)*6+1
         i2 = (nf-1)*6+IQfl_tr
         i3 = (nf-1)*6+IQfl_tr+1
         i4 = (nf-1)*6+IQfl

#if   ASCII == 1
         if (nf==1) open (1,file='qsdof3.dat' ,access='append')
         if (nf==2) open (1,file='qsdof32.dat',access='append')
         write (1,100)                                 &
               (      TIME)                           ,& ! 1
          (    (      UT_el  (NDFBT_el+i)), i=i1, i2 ),& ! 2- 4
          (    (R2D * UT_el  (NDFBT_el+i)), i=i3, i4 ),& ! 5- 7
          (    (      UT1_el (NDFBT_el+i)), i=i1, i2 ),& ! 8-10
          (    (R2D * UT1_el (NDFBT_el+i)), i=i3, i4 ),& !11-13
          (    (      UT2_el (NDFBT_el+i)), i=i1, i2 ),& !14-16
          (    (R2D * UT2_el (NDFBT_el+i)), i=i3, i4 )   !17-19
#elif ASCII == 0
         if (nf==1) open (1,file='qsdof3.bin' ,access='append',form='UNFORMATTED')
         if (nf==2) open (1,file='qsdof32.bin',access='append',form='UNFORMATTED')
         write (1    )                                 &
           sngl(      TIME)                           ,& ! 1
          (sngl(      UT_el  (NDFBT_el+i)), i=i1, i2 ),& ! 2- 4
          (sngl(R2D * UT_el  (NDFBT_el+i)), i=i3, i4 ),& ! 5- 7
          (sngl(      UT1_el (NDFBT_el+i)), i=i1, i2 ),& ! 8-10
          (sngl(R2D * UT1_el (NDFBT_el+i)), i=i3, i4 ),& !11-13
          (sngl(      UT2_el (NDFBT_el+i)), i=i1, i2 ),& !14-16
          (sngl(R2D * UT2_el (NDFBT_el+i)), i=i3, i4 )   !17-19
#endif
         close(1)
      enddo


!------ Transfer motions to Zcg of the floater (for INNWIND experinet WP4_25).
!      nf = 1
!         i1=(nf-1)*6+1
!         i2=(nf-1)*6+IQfl_tr
!         i3=(nf-1)*6+IQfl_tr+1
!         i4=(nf-1)*6+IQfl
!
!         RL(:)=0.d0; RL(3)=-0.3505d0
!         RG(1:3) =  UT_el(NDFBT_el+i1:NDFBT_el+i2) + matmul(  A_float(1:3,1:3), RL(1:3))
!        DRG(1:3) = UT1_el(NDFBT_el+i1:NDFBT_el+i2) + matmul( DA_float(1:3,1:3), RL(1:3))
!       DDRG(1:3) = UT2_el(NDFBT_el+i1:NDFBT_el+i2) + matmul(DDA_float(1:3,1:3), RL(1:3))
!
!#if   ASCII == 1
!         open  (1,file='qsdof30.dat' ,access='append')
!         write (1,100)                                 &
!               (      TIME)                           ,& ! 1
!          (    (      RG     (         i)), i= 1,  3 ),& ! 2- 4
!          (    (R2D * UT_el  (NDFBT_el+i)), i=i3, i4 ),& ! 5- 7
!          (    (     DRG     (         i)), i= 1,  3 ),& ! 8-10
!          (    (R2D * UT1_el (NDFBT_el+i)), i=i3, i4 ),& !11-13
!          (    (    DDRG     (         i)), i= 1,  3 ),& !14-16
!          (    (R2D * UT2_el (NDFBT_el+i)), i=i3, i4 )   !17-19
!#elif ASCII == 0
!         open  (1,file='qsdof30.bin' ,access='append',form='UNFORMATTED')
!         write (1    )                                 &
!           sngl(      TIME)                           ,& ! 1
!          (sngl(      RG     (         i)), i= 1,  3 ),& ! 2- 4
!          (sngl(R2D * UT_el  (NDFBT_el+i)), i=i3, i4 ),& ! 5- 7
!          (sngl(     DRG     (         i)), i= 1,  3 ),& ! 8-10
!          (sngl(R2D * UT1_el (NDFBT_el+i)), i=i3, i4 ),& !11-13
!          (sngl(    DDRG     (         i)), i= 1,  3 ),& !14-16
!          (sngl(R2D * UT2_el (NDFBT_el+i)), i=i3, i4 )   !17-19
!#endif
!         close(1)

   elseif ( (ICASE_el == 2).and.(NBODBTWT_el > 0) ) then
      nbod     = NBODBTWT_el
 !!   nnsub    = 1
      i        = NDFBT_el + NQSP + ACCU%NQSACC_el (nbod-1)!!! + ACCU%NQSBACC_el(nbod   ,nnsub-1)
#if   ASCII == 1
      open  (1,file='qsdof2.dat',access='append')
      write (1,100)      (TIME)                  ,&
                   (     (UT_el (j)), j=i+1,i+6 ),&
                   (     (UT1_el(j)), j=i+1,i+6 ),&
                   (     (UT2_el(j)), j=i+1,i+6 )
#elif ASCII == 0
      open  (1,file='qsdof2.bin',access='append',form='UNFORMATTED')
      write (1    )  sngl(TIME)                  ,&
                   ( sngl(UT_el (j)), j=i+1,i+6 ),&
                   ( sngl(UT1_el(j)), j=i+1,i+6 ),&
                   ( sngl(UT2_el(j)), j=i+1,i+6 )
#endif
      close(1)
   endif

 100  format (150en24.8e3)


 END Subroutine write_qsdof
!
!
!
!----------------------------------------------------------------------
 Subroutine Set_undeformed_blade
!----------------------------------------------------------------------

 Use Cbeam

   implicit None

   real(8) :: X0  (3)  , T0  (3)  , T0p (3)
   real(8) :: A   (3,3), AT0 (3,3), X0p (3), AT0A(3,3)
   integer :: nbod, nbsub, nel, nb_el, nn_el, nnod, nnod0, nn, i, j !, nn0_el


   do       nbod         = 1, NBLADE_el
!------------ wrt blade root
!!          AT0(1:3,1:3) = ATP_el  (nbod,1:3,1:3)               !after cone, after  pitch, before pre-curved angles
            AT0(1:3,1:3) = ATPWr_el(nbod,1:3,1:3)               !after cone, before pitch, before pre-curved angles
!!          nn0_el       = body(nbod)%NNODTGLB_el(1)
!!          AT0(1:3,1:3) = transf_mat(nn0_el)%AT_el(1:3,1:3)    !after cone, after  pitch, after  pre-curved angles

            nnod0        = 1        !1 only for 1st nel of 1st nbsub for each body, otherwise NNPE_el
            X0 (1:3)     = 0.d0
            X0p(1:3)     = 0.d0
            T0 (1:3)     = 0.d0
            T0p(1:3)     = 0.d0

      do    nbsub        = 1, body(nbod)%NBODSUB_el
            nb_el        = body(nbod)%NBODTGLB_el(nbsub)
            nn_el        = body(nbod)%NNODTGLB_el(nbsub)
            A(1:3,1:3)   = transf_mat(nn_el)%A_el(1:3,1:3)
            AT0A         = matmul( AT0, A )
            X0p(1:3)     = X0p(1:3) + X0(1:3)
            T0p(1:3)     = T0p(1:3) + T0(1:3)
         do nel          = 1, subbody(nb_el)%NTEB_el
         do nnod         = nnod0, NNPE_el, NNPE_el-1
            nn           =                            subbody(nb_el)%NODMB       (nel,nnod)
            i            = ACCU%NDFTACC_el(nb_el-1) + subbody(nb_el)%NDFPNBACC_el(nn-1)

            do j         = 1, 3
              X0(j)      = UT_el(i+IMDOF_el(j  ))
              T0(j)      = UT_el(i+IMDOF_el(j+3))
            enddo !j

            if (nbod<=NBODBTWT_el) X0(2) = X0(2) + subbody(nb_el)%HTA_el(nel+nnod/NNPE_el)

            X0           = matmul( AT0A,X0 )
            T0           = matmul( AT0A,T0 )

            subbody(nb_el)%Rinit (nel+nnod/NNPE_el, 1:3) = X0(1:3) + X0p(1:3)
            subbody(nb_el)%Rinit (nel+nnod/NNPE_el, 4:6) = T0(1:3) + T0p(1:3)

            write (88,100)                                         &
             (subbody(nb_el)%Rinit (nel+nnod/NNPE_el,j)    ,j=1,3),&
             (subbody(nb_el)%Rinit (nel+nnod/NNPE_el,j)*R2D,j=4,6)
         enddo !nnod
            nnod0        = NNPE_el
         enddo !nel
      enddo !nbsub
            write (88,*)
            write (88,*)
   enddo !nbod


 100  format (150en21.8e3)


 END Subroutine Set_undeformed_blade
