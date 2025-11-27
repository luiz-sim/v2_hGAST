!----------------------------------------------------------------------
! Calls:
!
! FLAP_Anemometer
!    tuner             -> CYCLIC_MOM                    -> LOADS_TRANS0_NEW,  ROT_MATRIX0
!                         PID                           
!    lowpasswind       -> NEWMARKN                      -> LUDCMP, LUBKSB
!    [BI,TRI,TETRA]LIN -> LOCATE1                       
!    ACTUATOR_1st      -> TIME_LAG1                     -> FLAP_NEWMARK

! CONTROL_FLAP_CYCL
!    CYCLIC_MOM        -> LOADS_TRANS0_NEW, ROT_MATRIX0
!    filterN           -> NEWMARKN                      -> LUDCMP, LUBKSB
!    PID_CONTROL       -> FLAP_NEWMARK
!    ACTUATOR_1st      -> TIME_LAG1                     -> FLAP_NEWMARK
!
! CONTROL_FLAP_hubaccel
!    beam_kinematics
!    lowpasstower      -> NEWMARKN                      -> LUDCMP, LUBKSB
!    PID_CONTROL       -> FLAP_NEWMARK
!----------------------------------------------------------------------
 Subroutine FLAP_Control (it,ISTEP)
!----------------------------------------------------------------------
 Use Craft
 Use Cbeam
 Use Foilfs_mod

   implicit none

   integer, intent(in) :: it, ISTEP


!qq
!! return
   if (it>1 .or. ISTEP /=1 .or. INDXFlap==0) return

!qqcall CONTROL_FLAP_hubaccel
   call FLAP_Anemometer
   call CONTROL_FLAP_CYCL                       !bld mom ->tow Myaw Mtilt -> Myaw, Mtilt -> bld...


 END Subroutine flap_control
!----------------------------------------------------------------------
 Subroutine FLAP_Anemometer
!----------------------------------------------------------------------
       
 Use Craft
 Use Cbeam
 Use Foilfs_mod
 Use Ctrl_mod

   implicit none

   integer, parameter  :: NH = 2 !Number of harmonics for flap motion

   type      :: type_timelag
     real(8) :: VARP=0.d0                                   ! input  var
     real(8) :: YP  =0.d0 , YP1=0.d0 , YP2=0.d0             ! 1st order filter previous time vars
     real(8) :: TVLAG=0.01d0                                ! time delay
   END type type_timelag

   type (type_timelag)  , save :: timelag(3)                ! NBLADE_el
   real(8)              , save :: Flap_min =-10.d0, Flap_max=10.d0, Flap_vel=20d0 !degrees
   real(8)                     :: VAR, TVLAG,  YP, YP1, YP2

   integer,               save :: initialize=0
   real(8),dimension(NH), save :: Kp  = 0.d0, dpsi  = 0.d0
   real(8)                     :: FlapAng(NBLADE_el), FlapVel(NBLADE_el), psiout(NBLADE_el) , psi

   integer, parameter :: order = 2 , N = 2*order
   real(8),save       :: Y_vel (N), Y1_vel (N) , Y2_vel (N)
   real(8),save       :: Y_yaw (N), Y1_yaw (N) , Y2_yaw (N)
   real(8),save       :: Y_inc (N), Y1_inc (N) , Y2_inc (N)
   real(8)            :: YP_vel(N), YP1_vel(N) , YP2_vel(N)
   real(8)            :: YP_yaw(N), YP1_yaw(N) , YP2_yaw(N)
   real(8)            :: YP_inc(N), YP1_inc(N) , YP2_inc(N)
   real(8)            :: velA     , velF !Actual and Filtered U
   real(8)            :: sheA     , sheF !Actual and Filtered shear exponent
   real(8)            :: yawA     , yawF !Actual and Filtered Yaw
   real(8)            :: incA     , incF !Actual and Filtered Inc
   integer            :: nbod, i , h

   integer,save              :: IwindF   !1:mean, 2:Instant, 3:Filtered
   integer,save              :: nbinT, nshT, nyawT, nincT ! U, Shear, Yaw, Inc
   integer                   :: ibin , ish , iyaw , iinc  ! U, Shear, Yaw, Inc
   real(8),save, allocatable :: spinnerU   (:        )  !(nbinT                        ) )
   real(8),save, allocatable :: spinnerS   (:        )  !(       nshT                  ) )
   real(8),save, allocatable :: spinnerY   (:        )  !(             nyawT           ) )
   real(8),save, allocatable :: spinnerI   (:        )  !(                    nincT    ) )
   real(8),save, allocatable :: spinnerKp  (:,:,:,:,:)  !(nbinT, nshT, nyawT, nincT, NH) )
   real(8),save, allocatable :: spinnerdpsi(:,:,:,:,:)  !(nbinT, nshT, nyawT, nincT, NH) )
   real(8),      allocatable :: tmp4       (:,:,:,:  )  !(nbinT, nshT, nyawT, nincT    ) )
   real(8),      allocatable :: tmp3       (:,  :,:  )  !(nbinT,       nyawT, nincT    ) )
   real(8),      allocatable :: tmp2       (:,  :    )  !(nbinT,       nyawT           ) )


   if ( IFLAP_typ == 30 ) return
   if ( INDXFlap  ==  2 ) return

!--- Initialization
   if (initialize==0) then
!------ filter vars for U, Yaw, Inc
      Y_vel (:)=0.d0; Y1_vel (:)=0.d0; Y2_vel (:)=0.d0;
      Y_yaw (:)=0.d0; Y1_yaw (:)=0.d0; Y2_yaw (:)=0.d0;
      Y_inc (:)=0.d0; Y1_inc (:)=0.d0; Y2_inc (:)=0.d0;


      if (IFLAP_typ == 21 .or. &
          IFLAP_typ == 22        ) then

         if (allocated(spinnerU)) goto 10

         open(101,file='spinner.inp')
            read(101,*) !title
            read(101,*) IwindF !1:mean, 2:Instant, 3:Filtered
            read(101,*) nbinT, nshT, nyawT, nincT ! U, Shear, Yaw, Inc
            if (nbinT>0.and.nshT>0.and.nyawT>0.and.nincT>0) then
               Allocate ( spinnerU   (nbinT                        ) )
               Allocate ( spinnerS   (       nshT                  ) )
               Allocate ( spinnerY   (             nyawT           ) )
               Allocate ( spinnerI   (                    nincT    ) )
               Allocate ( spinnerKp  (nbinT, nshT, nyawT, nincT, NH) )
               Allocate ( spinnerdpsi(nbinT, nshT, nyawT, nincT, NH) )


               do    ish  = 1, nshT
               do    ibin = 1, nbinT
                     read(101,*)                                         &
                      spinnerU    (ibin                     )           ,&
                      spinnerS    (      ish                )          
                  do iinc = 1, nincT
                  do iyaw = 1, nyawT
                     read(101,*)                                         &
                      spinnerY    (           iyaw          )           ,&
                      spinnerI    (                 iinc    )           ,&
                     (spinnerKp   (ibin, ish, iyaw, iinc, h )           ,&
                      spinnerdpsi (ibin, ish, iyaw, iinc, h ) ,h=1, NH)
                  enddo !iyaw
                  enddo !iinc
               enddo !ibin
               enddo !ish

                      spinnerY    (:        ) = spinnerY   (:        )/R2D;
                      spinnerI    (:        ) = spinnerI   (:        )/R2D;
                      spinnerdpsi (:,:,:,:,:) = spinnerdpsi(:,:,:,:,:)/R2D;
            else
                write(*,*)'zero dimensioned vectors in spinner.inp'; stop
            endif

         close(101)
 10   endif !IFLAP_typ

      initialize=1
   endif !initialize


   if     (IFLAP_typ == 1 .or.&
           IFLAP_typ == 2      ) then
!------ Prescribed motion (collective or cyclic)
      Kp  (1) = FLAP_ampl !degrees
      dpsi(1) = 0.d0

   elseif (IFLAP_typ ==20      ) then
!------ Tune the flap feed forward motion
       call tuner (Kp, dpsi)

!------ mean inflow values
       velA    = VELHUBT
       sheA    = SHEXP
       yawA    = WINDYAW 
       incA    = WINDINC

   elseif (IFLAP_typ ==21 .or.&
           IFLAP_typ ==22        ) then
!------ spinner anemometer (cyclic)

!------ set current horizontal velocity, yaw and inclination
     if( IwindF==1) then !1:mean, 2:Instant, 3:Filtered
!------ mean inflow values
       velA    = VELHUBT
       sheA    = SHEXP
       yawA    = WINDYAW 
       incA    = WINDINC
     elseif( IwindF==2.or.IwindF==3) then !1:mean, 2:Instant, 3:Filtered
!------ instantaneous inflow values (caused due to turbulence)
       velA    = dsqrt (dot_product(HUB_VEL(1:3),HUB_VEL(1:3)))
       sheA    = SHEXP
       yawA    =-datan2(HUB_VEL(2)           ,HUB_VEL(1))                          !calculate actual yaw
       incA    = datan2(HUB_VEL(3)*dcos(yawA),HUB_VEL(1))                          !calculate actual inc 
! U={Ux,Uy,Uz}=Az(-yaw).Ay(-inc).{Uo,0,0}= Uo cos(inc) cos(yaw),
!                                         -Uo cos(inc) sin(yaw),
!                                          Uo sin(inc)
! -->  yaw=-atan2(Uy         ,Ux)
! -->  inc= atan2(Uz cos(yaw),Ux)
     endif

       velF    = velA
       sheF    = sheA
       yawF    = yawA
       incF    = incA

!     if (initialize==1) then
!        Y1_vel (2)=velA
!        Y1_yaw (2)=yawA
!        Y1_inc (2)=incA
!        initialize=2
!     endif
   
!------ Filtered velocity, yaw and inclination
     if( IwindF==3) then !1:mean, 2:Instant, 3:Filtered
       YP_vel (:)=Y_vel (:); YP1_vel (:)=Y1_vel (:); YP2_vel (:)=Y2_vel (:);       !update filter solution vector
       YP_yaw (:)=Y_yaw (:); YP1_yaw (:)=Y1_yaw (:); YP2_yaw (:)=Y2_yaw (:);       !update filter solution vector
       YP_inc (:)=Y_inc (:); YP1_inc (:)=Y1_inc (:); YP2_inc (:)=Y2_inc (:);       !update filter solution vector

       call lowpasswind ( YP_vel, YP1_vel, YP2_vel, Y_vel, Y1_vel, Y2_vel, velA, velF, DT_el ) !vel  filter
       call lowpasswind ( YP_yaw, YP1_yaw, YP2_yaw, Y_yaw, Y1_yaw, Y2_yaw, yawA, yawF, DT_el ) !yaw  filter
       call lowpasswind ( YP_inc, YP1_inc, YP2_inc, Y_inc, Y1_inc, Y2_inc, incA, incF, DT_el ) !inc  filter
     endif

!------- Consider outer bounds
       velF      = min(  max( velF , spinnerU(1) )  ,  spinnerU(nbinT)  )
       sheF      = min(  max( sheF , spinnerS(1) )  ,  spinnerS(nshT )  )
       yawF      = min(  max( yawF , spinnerY(1) )  ,  spinnerY(nyawT)  )
       incF      = min(  max( incF , spinnerI(1) )  ,  spinnerI(nincT)  )

!------ set interpolated values based on velocity, shear, yaw and inclination.
      if     (nbinT>1.and.nshT>1.and.nyawT>1.and.nincT>1) then
         Allocate ( tmp4 (nbinT, nshT, nyawT, nincT) )

         do h = 1, NH
            tmp4(1:nbinT, 1:nshT, 1:nyawT, 1:nincT) = spinnerKp  (1:nbinT, 1:nshT, 1:nyawT, 1:nincT, h);
!           call TETRALIN (X1A     ,X2A      ,X3A     ,X4A     ,YA  ,X1 ,X2  ,X3  ,X4  ,Y      ,M    ,N   ,O    ,P    ,MM   ,NM  ,OM   ,PM   )
            call TETRALIN (spinnerU,spinnerS,spinnerY,spinnerI,tmp4,velF,sheF,yawF,incF,Kp  (h),nbinT,nshT,nyawT,nincT,nbinT,nshT,nyawT,nincT)

            tmp4(1:nbinT, 1:nshT, 1:nyawT, 1:nincT) = spinnerdpsi(1:nbinT, 1:nshT, 1:nyawT, 1:nincT, h);
            call TETRALIN (spinnerU,spinnerS,spinnerY,spinnerI,tmp4,velF,sheF,yawF,incF,dpsi(h),nbinT,nshT,nyawT,nincT,nbinT,nshT,nyawT,nincT)
         enddo !h

         Deallocate ( tmp4 )

      elseif (nbinT>1.and.nshT==1.and.nyawT>1.and.nincT>1) then

!--------- No shear interpolation
         Allocate ( tmp3 (nbinT, nyawT, nincT) )

         do h = 1, NH
            tmp3(1:nbinT, 1:nyawT, 1:nincT) = spinnerKp  (1:nbinT, 1, 1:nyawT, 1:nincT, h);
!           call TRILIN   (X1A     ,X2A     ,X3A     ,YA  ,X1  ,X2  ,X3  ,Y      ,M    ,N    ,O    ,MM   ,NM   ,OM   )
            call TRILIN   (spinnerU,spinnerY,spinnerI,tmp3,velF,yawF,incF,Kp  (h),nbinT,nyawT,nincT,nbinT,nyawT,nincT)

            tmp3(1:nbinT, 1:nyawT, 1:nincT) = spinnerdpsi(1:nbinT, 1, 1:nyawT, 1:nincT, h);
            call TRILIN2  (spinnerU,spinnerY,spinnerI,tmp3,velF,yawF,incF,dpsi(h),nbinT,nyawT,nincT,nbinT,nyawT,nincT)
         enddo !h

         Deallocate ( tmp3 )

      elseif (nbinT>1.and.nshT==1.and.nyawT>1.and.nincT==1) then
!--------- No shear interpolation, No inclination interpolation
         Allocate ( tmp2 (nbinT, nyawT) )

         do h = 1, NH
            tmp2(1:nbinT, 1:nyawT) = spinnerKp  (1:nbinT, 1, 1:nyawT, 1, h);
!           call BILIN   (X1A     ,X2A     ,YA  ,X1  ,X2  ,Y      ,M    ,N    ,MM   ,NM   )
            call BILIN   (spinnerU,spinnerY,tmp2,velF,yawF,Kp  (h),nbinT,nyawT,nbinT,nyawT)

            tmp2(1:nbinT, 1:nyawT) = spinnerdpsi(1:nbinT, 1, 1:nyawT, 1, h);
            call BILIN   (spinnerU,spinnerY,tmp2,velF,yawF,dpsi(h),nbinT,nyawT,nbinT,nyawT)
         enddo !h

         Deallocate ( tmp2 )

      else
         write(*,*)'error in spinner interpolation';stop
      endif
   endif !IFLAP_typ

!!!--- Paper flap Torque 2018 jprosp
!!     1p and 3p collective flap motion without delay
!!   do nbod                    = 1, NBLADE_el
!!         psi           = -UT_el(NDFBT_el+NQSW) !-UThub_el             !collective
!!         psiout (nbod) = dmod (psi,PI2)
!!         FlapAng(nbod) = 0.d0
!!      if     (IFLAP_typ == 1) then
!!         FlapAng(nbod) = Kp(1)*cos( 1.d0*psi )
!!         FlapVel(nbod) = Kp(1)*sin( 1.d0*psi ) *        UT1_el(NDFBT_el+NQSW)
!!      elseif (IFLAP_typ == 2) then
!!         FlapAng(nbod) = Kp(1)*cos( 3.d0*psi )
!!         FlapVel(nbod) = Kp(1)*sin( 3.d0*psi ) * 3.d0 * UT1_el(NDFBT_el+NQSW)
!!      endif
!!
!!      if (TIME>=dble(NTIME_INI)*DT_el) then
!!         do i = 1, NSTRIP
!!            if (IFLAPLOC(nbod,i) == 1) FLCONTR (nbod,i) = FlapAng(nbod) !degrees
!!         enddo !i
!!      endif
!!   enddo !nbod
!!!qq 
!!   goto 222


   do nbod                    = 1, NBLADE_el
      if (IFLAP_typ == 1) psi = -UT_el(NDFBT_el+NQSW)-UThub_el             !collective
      if (IFLAP_typ >= 2) psi = -UT_el(NDFBT_el+NQSW)-UThub_el +PHI0(nbod) !cyclic, spinner anemometer
         psiout (nbod) = dmod (psi,PI2)
         FlapAng(nbod) = 0.d0
      do h = 1, NH
         FlapAng(nbod) = FlapAng(nbod) + Kp(h)*cos( h*psi-dpsi(h) )
      enddo
!qq tower top acceleration contribution on the flap action
         FlapAng(nbod) = FlapAng(nbod) + BetaF_accel*R2D
!!       FlapAng(nbod) = min (max(FlapAng(nbod), Flap_min), Flap_max)

!------ time lag (delay, usually ~0.1s)
      VAR         = FlapAng(nbod)
      TVLAG       = timelag(nbod)%TVLAG
      YP          = timelag(nbod)%YP
      YP1         = timelag(nbod)%YP1
      YP2         = timelag(nbod)%YP2

!     call ACTUATOR_1st (VAR, TVLAG, DT   , MinAng  , MaxAng  , MaxVel  , YP, YP1, YP2)
      call ACTUATOR_1st (VAR, TVLAG, DT_el, Flap_min, Flap_max, Flap_vel, YP, YP1, YP2)

      timelag(nbod)%VARP = VAR
      timelag(nbod)%YP   = YP
      timelag(nbod)%YP1  = YP1
      timelag(nbod)%YP2  = YP2
      FlapAng(nbod)      = YP1
      FlapVel(nbod)      = YP2

      if (TIME>=dble(NTIME_INI)*DT_el) then
         do i = 1, NSTRIP
            if (IFLAPLOC(nbod,i) == 1) FLCONTR (nbod,i) = FlapAng(nbod) !degrees
         enddo !i
      endif
   enddo !nbod

!!222 continue

!--- Output
   if (IFLAP_typ ==20) then
#    if   ASCII == 1
      open(23,file='control_IFCS.dat',access='append')
       write(23,'(50f15.5)')                                                            &
#    elif ASCII == 0
      open(23,file='control_IFCS.bin',access='append',form='UNFORMATTED')
       write(23            )                                                            &
#    endif
         sngl(TIME)                                                                   , & ! 1
        (sngl(psiout(i)*R2D), sngl(FlapAng(i)     ), sngl(FlapVel(i)), i=1, NBLADE_el), & ! 2,5,8 / 3,6,9 / 4,7,10
        (sngl(Kp    (h)    ), sngl(dpsi(h)    *R2D)                  , h=1, NH       ), & !11, 12 h=1 13,14 h=2
         sngl(velA         )                                                          , & !15
         sngl(yawA     *R2D)                                                          , & !16
         sngl(incA     *R2D)                                                              !17
      close(23)

   elseif (IFLAP_typ ==21 .or.&
           IFLAP_typ ==22        ) then

#    if   ASCII == 1
      open(23,file='control_IFCS.dat',access='append')
       write(23,'(50f15.5)')                                              &
#    elif ASCII == 0
      open(23,file='control_IFCS.bin',access='append',form='UNFORMATTED')
       write(23            )                                              &
#    endif
         sngl(TIME)                                                                   , & ! 1
        (sngl(psiout(i)*R2D), sngl(FlapAng(i)     ), sngl(FlapVel(i)), i=1, NBLADE_el), & ! 2,5,8 / 3,6,9 / 4,7,10
        (sngl(Kp    (h)    ), sngl(dpsi(h)    *R2D)                  , h=1, NH       ), & !11, 12 h=1 13,14 h=2
         sngl(velA         ), sngl(velF           )                                   , & !15,16
         sngl(yawA     *R2D), sngl(yawF       *R2D)                                   , & !17,18
         sngl(incA     *R2D), sngl(incF       *R2D)                                   , & !19,20
         sngl(Y_vel(1:2)   ), sngl(Y1_vel(1:2)    ), sngl(Y2_vel(1:2))                , & !21-26
         sngl(Y_yaw(1:2)   ), sngl(Y1_yaw(1:2)    ), sngl(Y2_yaw(1:2))                , & !27-32
         sngl(Y_inc(1:2)   ), sngl(Y1_inc(1:2)    ), sngl(Y2_inc(1:2))                    !33-38
      close(23)
   elseif (IFLAP_typ == 1 .or.&
           IFLAP_typ == 2        ) then

#    if   ASCII == 1
      open(23,file='flap.dat',access='append')
       write(23,'(50f15.5)')                                              &
#    elif ASCII == 0
      open(23,file='flap.bin',access='append',form='UNFORMATTED')
       write(23            )                                              &
#    endif
         sngl(TIME)                                                                   , & ! 1
        (sngl(psiout(i)*R2D), sngl(FlapAng(i)     ), sngl(FlapVel(i)), i=1, NBLADE_el), & ! 2,5,8 / 3,6,9 / 4,7,10
        (sngl(Kp    (h)    ), sngl(dpsi(h)    *R2D)                  , h=1, NH       )    !11, 12 h=1 13,14 h=2
      close(23)
   endif


 END Subroutine FLAP_Anemometer
!------------------------------------------------------------------------
 Subroutine tuner (beta_flap_amp, dpsi)
!------------------------------------------------------------------------

 use Cbeam
 Use Craft

   implicit none

   integer, parameter     :: NH = 2 !Number of harmonics for flap motion
   real(8), intent(  out) :: beta_flap_amp(NH) , dpsi(NH)

   integer, save          :: initialize = 0
   real(8), save          :: theta_tilt(NH) = 0.d0, theta_yaw(NH) = 0.d0

   integer                :: init_cycles   = 1
   integer, save          :: current_cycle = 1
   real(8), save          :: psip

   integer :: h , j
   real(8) :: M_NRT_C (NH,2), psi  (NBLADE_el)

   integer,                     save :: length = 0
   real(8), dimension(NH,1000), save :: dummy1 , dummy2
   logical                           :: check
   integer                           :: first , last
   real(8), allocatable    ,    save :: Mtilt(:,:) , Myaw(:,:)
   real(8), dimension(NH  ),    save :: Mtilt0 = 0.d0 , Myaw0 = 0.d0
   real(8), dimension(NH,2)          :: error , output
!--- for GS
   real(8) :: V,GSS !!, Va,Vb, GSa,GSb


!--- Find azimuth angle
   call CYCLIC_MOM ( M_NRT_C, psi, NH)

   M_NRT_C(1:NH,1:2) = -2.d0/3.d0*M_NRT_C(1:NH,1:2);

   psi(:) = mod(psi(:),2*pi)*R2D

   if (initialize==0) then
      psip       = psi(1)
      initialize = 1
   endif

!--- Current cycle counter
   if (psi(1)-psip <0.d0) then
      current_cycle = current_cycle + 1
   endif

!  write(*,*)'tuner:',psi(1),psip,current_cycle,init_cycles

   psip   = psi(1)

   open(372,file='nikos.dat',access='append')
      write(372,'(7f15.5,i5)') TIME , psi(1) , (M_NRT_C(1,j) , j=1,2)
   close(372)


!--- Initial cycles that will not be taken into account
   if ( current_cycle < init_cycles ) return


!--- Penultimate cycle before tuning, figuring out the size of the window
   if ( current_cycle == init_cycles ) then
         length           = length + 1
         dummy1(:,length) = M_NRT_C(:,1)
         dummy2(:,length) = M_NRT_C(:,2)
      return
   endif


!--- Fishing new values of loads
!--- Window allocation--
   check = allocated(Mtilt) .and. allocated(Myaw)
   if ( .not.(check) ) then 
      allocate( Mtilt(NH,length), Myaw(NH,length) )
      do    h          = 1, NH
         do j          = 1, length
            Mtilt(h,j) = dummy1(h,j)
            Myaw (h,j) = dummy2(h,j)
         enddo
      enddo
   endif

!--- Window shift--
   first  = lbound(Myaw,2)
   last   = ubound(Myaw,2)
   Mtilt (:,first:last-1) = Mtilt  (:,first+1:last)
   Myaw  (:,first:last-1) = Myaw   (:,first+1:last)
   Mtilt (:,      last  ) = M_NRT_C(:,1)
   Myaw  (:,      last  ) = M_NRT_C(:,2)


!--- Find new pair of (beta_flap_amp,dpsi)
   do h = 1,NH
      Mtilt0(h)     = sum(Mtilt(h,:))/dble(length)
      Myaw0 (h)     = sum(Myaw (h,:))/dble(length)
   enddo

   error(:,1)       = Mtilt0(:)
   error(:,2)       = Myaw0 (:)

!--- GS
   GSS=1.d0
   V  = HUB_VEL(1)
!  Va= 4.0d0; GSa=1.00d0
!  Vb=11.4d0; GSb=0.25d0

!  if     (HUB_VEL(1)<=Va) then; GSS=GSa
!  elseif (HUB_VEL(1)>=Vb) then; GSS=GSb
!  else                        ; GSS=(Vb-V)/(Vb-Va)*GSa + (V-Va)/(Vb-Va)*GSb
!  endif

!--- PID
   call PID ( error , GSS, DT_el , output )

   theta_tilt(:)    = output(:,1)
   theta_yaw (:)    = output(:,2)


!--- Calculate the pair (beta_flap_amp,dpsi)
   beta_flap_amp(:) = sqrt(theta_tilt(:)**2 + theta_yaw(:)**2)
   dpsi(:)          = atan2( theta_yaw(:) , theta_tilt(:) )

!--- Turn radians to degrees
   beta_flap_amp(:) = beta_flap_amp(:)*R2D
   dpsi(:)          = dpsi(:)              ! We need dpsi in radians
   theta_tilt(:)    = theta_tilt(:)*R2D
   theta_yaw (:)    = theta_yaw (:)*R2D


   open(373,file='nikos1.dat',access='append')
      write(373,'(1000e15.7)') TIME , ( beta_flap_amp(h) , dpsi(h)*R2D , theta_tilt(h) , theta_yaw(h) , Mtilt0(h) , Myaw0(h) , h=1,NH)
      !                        1             2,8            3,9              4,10            5,11        6,12        7,13
   close(373)
   open(373,file='tuner.dat',access='append')
      write(373,'(1000e15.7)') TIME , ( beta_flap_amp(h) , dpsi(h)*R2D , h=1,NH ), V, GSS
      !                        1             2,4            3,5                    6, 7
   close(373)


 END Subroutine tuner
!------------------------------------------------------------------------
 Subroutine PID ( error, GS, dt, output )
!------------------------------------------------------------------------

   implicit none

   integer, parameter                      :: NH = 2 !Number of harmonics for flap motion

   real(8), dimension(NH,2), intent(in   ) :: error
   real(8),                  intent(in   ) :: GS, dt
   real(8), dimension(NH,2), intent(  out) :: output

   integer, save :: initialize = 0
   real(8), save :: GAIN, CONP, CONI, COND !!, GS_NH2

   real(8), dimension(NH,2), save :: errorp = 0.d0, outputp = 0.d0
   real(8), dimension(NH,2)       :: derror
   real(8)                        :: Gain_h
   integer                        :: h

 
   if (initialize==0) then
       open(338,file='pid.inp')
           read(338,*)
           read(338,*) GAIN , CONP , CONI , COND !, GS_NH2
       close(338)
       initialize  = initialize + 1
   endif


   do h = 1,NH
                Gain_h = GS*GAIN
!!!   if (h==1) Gain_h = GS*GAIN
!!!   if (h==2) Gain_h = GS*GAIN /GS_NH2
      derror (h,:) =(error  (h,:) - errorp(h,:))/dt
      output (h,:) = outputp(h,:) + Gain_h * ( CONP*derror(h,:) + CONI*error(h,:) )*dt

      errorp (h,:) = error (h,:)
      outputp(h,:) = output(h,:)
   enddo


 END Subroutine PID
!----------------------------------------------------------------------
!
!
! PID for flap, based on bld root moments tranformed to NR frame.
!
!----------------------------------------------------------------------
 Subroutine CONTROL_FLAP_CYCL
!-----------------------------------------------------------------------
       
 Use Craft
 Use Cbeam
 Use Foilfs_mod
 Use Ctrl_mod

   implicit none

   integer, parameter   :: NH = 2 !Number of harmonics for flap motion

!--- filter and pid vars
   integer, parameter   :: N = 4

   type      :: type_filter
     real(8) :: YP(N  )=0., YP1(N)=0., YP2(N)=0.            ! filter previous time vars
     real(8) :: A (N,N)   , B  (N)   , C  (N)   , D         ! filter parameters
     real(8) :: XX        , XX1                             ! input  var
     real(8) :: XXP=0.    , XXP1=0.                         ! input  var previous time
     real(8) :: YYP=0.    , YYP1=0.                         ! output var previous time
   END type type_filter

   type      :: type_pid
     real(8) :: YP=0.     , YP1=0.   , YP2=0.               ! PID previous time vars
!!!  real(8) :: GAIN      , CONP     , CONI     , COND      ! PID parameters
     real(8) ::             CONP     , CONI     , COND      ! PID parameters
   END type type_pid

   type      :: type_timelag
     real(8) :: VARP=0.                                     ! input  var
     real(8) :: YP  =0.   , YP1=0.   , YP2=0.               ! 1st order filter previous time vars
     real(8) :: TVLAG                                       ! time delay
   END type type_timelag

   type (type_filter ), save :: filter (4*NH)
   type (type_pid    ), save :: pid    (2*NH)
   type (type_timelag), save :: timelag(3)                  ! NBLADE_el
   real(8)            , save :: GAIN_PITCH, MaxPitchAmpl
   real(8)            , save :: GAIN_FLAP , MinFlapAng, MaxFlapAng, MaxFlapVel !, MaxFlapAcc
   integer            , save :: initialize=0
   integer            , save :: nhT

   integer :: nn, nbod, i,j, kf, ki, h
   real(8) :: Y (N  ), Y1 (N), Y2 (N)
   real(8) :: YP(N  ), YP1(N), YP2(N)
   real(8) :: A (N,N), B  (N), C  (N), D
   real(8) :: XX     , XX1
   real(8) :: YY     , YY1
   real(8) :: GAIN   , CONP  , CONI  , COND 
   real(8) :: VAR    , DVAR  , DDVAR
   real(8) :: TVLAG

   real(8) :: theta_yaw(NH), theta_tilt(NH)
   real(8) :: TIME_INI, TD1, TD2
   real(8) :: FlapAng(NBLADE_el), Bang(NBLADE_el), FlapVel(NBLADE_el)
   real(8) :: PVEL   (NBLADE_el), PACC(NBLADE_el)

   real(8) :: M_NRT_C (NH,2), psi (NBLADE_el)
   character :: CNUM(2)
   real(8) :: GSS, aa


   if (IFLAP_typ /= 22 .and.&
       IFLAP_typ /= 30       ) return


!--- Initialization
   if (initialize == 0) then
!------ Filter parameters for two filters (Cyclic load)
     open(101,file='flapPID.inp')

!------ read IFC parameters
!------ read filters
         read (101,*) nhT
   do    h  = 1, nhT
      do nn = 1, 4        !3p Myaw, 6p Myaw, 3p Mtilt, 6p Mtilt
         kf = (h-1)*4 + nn
         read (101,*)
         do i=1,N
            if (i==1) read (101,*)      &
             (filter(kf)%A(i,j),j=1,N) ,&
              filter(kf)%B(i  )        ,&
              filter(kf)%C(i  )        ,&
              filter(kf)%D
            if (i> 1) read (101,*)      &
             (filter(kf)%A(i,j),j=1,N) ,&
              filter(kf)%B(i  )        ,&
              filter(kf)%C(i  )
         enddo !i
      enddo !nn

!------ read PID
      do nn = 1, 2        !PID for Myaw, PID for Mtilt
         ki = (h-1)*2 + nn
         read (101,*)!PID
         read (101,*) pid(ki)%CONP, pid(ki)%CONI, pid(ki)%COND
      enddo !nn
   enddo !h

         read (101,*)!Gain for cyclic flap
         read (101,*) GAIN_FLAP

         if (IFLAP_typ == 22) GAIN_FLAP = 0.d0

         read (101,*)!Delay for flap [sec]
         read (101,*) TVLAG
                      timelag(1:NBLADE_el)%TVLAG = TVLAG

         read (101,*)!Min, Max Flap Angle [deg], max velocity[deg/s], max acceleration[deg/s^2]
         read (101,*) MinFlapAng, MaxFlapAng, MaxFlapVel !, MaxFlapAcc
                      MinFlapAng = MinFlapAng/R2D
                      MaxFlapAng = MaxFlapAng/R2D
                      MaxFlapVel = MaxFlapVel/R2D
                     !MaxFlapAcc = MaxFlapAcc/R2D

!------ read IPC parameters
         read (101,*)!Gain for individual pitch
         read (101,*) GAIN_PITCH

         read (101,*)!Max amplitude of individual pitch [deg]
         read (101,*) MaxPitchAmpl
                      MaxPitchAmpl = MaxPitchAmpl/R2D
     close(101)

     initialize = 1
   endif !initialization

      TIME_INI      = dble(NTIME_INI)*DT_el
      TD1           =  5.d0
      TD2           = 10.d0

   if (TIME<TIME_INI) return

!--- Cyclic sensor (moments)
   call CYCLIC_MOM ( M_NRT_C, psi, NH)

!     filter(1)%XX  =  M_NRT(2)  !Myaw
!     filter(3)%XX  =  M_NRT(1)  !Mtilt
   do h                = 1, nhT
      kf               = (h-1)*4
      filter(kf+1)%XX  =  M_NRT_C(h,2)  !Myaw
      filter(kf+3)%XX  =  M_NRT_C(h,1)  !Mtilt
      filter(kf+1)%XX1 = (filter(kf+1)%XX - filter(kf+1)%XXP ) / DT_el
      filter(kf+3)%XX1 = (filter(kf+3)%XX - filter(kf+3)%XXP ) / DT_el
   enddo !h


!--- 3p/6p Band low pass filters and PID (I) for Myaw and Mtilt to derive byaw and btilt
   do h              = 1, nhT
   do nn             = 1, 2 !1: Myaw, 2: Mtilt
!------ Band Stop 6p
      kf             = (h-1)*4 + (nn-1)*2 + 1
      A  (1:N,1:N)   = filter(kf)%A  (1:N,1:N);
      B  (1:N    )   = filter(kf)%B  (1:N    );
      C  (1:N    )   = filter(kf)%C  (1:N    );
      D              = filter(kf)%D
      YP (1:N    )   = filter(kf)%YP (1:N    );
      YP1(1:N    )   = filter(kf)%YP1(1:N    );
      YP2(1:N    )   = filter(kf)%YP2(1:N    );
      XX             = filter(kf)%XX
      XX1            = filter(kf)%XX1

      call filterN ( A,B,C,D, YP,YP1,YP2 ,Y,Y1,Y2, DT_el, N, XX,XX1, YY,YY1 )

      filter(kf)%YP (1:N) = Y (1:N);
      filter(kf)%YP1(1:N) = Y1(1:N);
      filter(kf)%YP2(1:N) = Y2(1:N);
      filter(kf)%XXP      = XX
      filter(kf)%XXP1     = XX1
      filter(kf)%YYP      = YY
      filter(kf)%YYP1     = YY1

!------ Band Stop 3p
      kf             = (h-1)*4 + (nn-1)*2 + 2
      filter(kf)%XX  = YY
      filter(kf)%XX1 = YY1
      A  (1:N,1:N)   = filter(kf)%A  (1:N,1:N);
      B  (1:N    )   = filter(kf)%B  (1:N    );
      C  (1:N    )   = filter(kf)%C  (1:N    );
      D              = filter(kf)%D
      YP (1:N    )   = filter(kf)%YP (1:N    );
      YP1(1:N    )   = filter(kf)%YP1(1:N    );
      YP2(1:N    )   = filter(kf)%YP2(1:N    );
      XX             = filter(kf)%XX
      XX1            = filter(kf)%XX1

      call filterN ( A,B,C,D, YP,YP1,YP2 ,Y,Y1,Y2, DT_el, N, XX,XX1, YY,YY1 )

      filter(kf)%YP (1:N) = Y (1:N);
      filter(kf)%YP1(1:N) = Y1(1:N);
      filter(kf)%YP2(1:N) = Y2(1:N);
      filter(kf)%XXP      = XX
      filter(kf)%XXP1     = XX1
      filter(kf)%YYP      = YY
      filter(kf)%YYP1     = YY1


!------ PID (I)
!qq GS for both Myaw, Mtilt based on bld_pitch (function of Uwind)
     aa  = 0.3d0
!    GSS = (aa+(1.d0-aa)* min(PITCH_Coll_filt*R2D,08.d0) /08.d0)
     GSS = 1.d0
      ki          = (h-1)*2 + nn
      VAR         = YY
      DVAR        = YY1
      DDVAR       = 0.d0
      GAIN        = GSS
      CONP        = pid(ki)%CONP
      CONI        = pid(ki)%CONI
      COND        = pid(ki)%COND
      YP  (1)     = pid(ki)%YP
      YP1 (1)     = pid(ki)%YP1
      YP2 (1)     = pid(ki)%YP2
      
      call PID2_CONTROL( GAIN,     DT_el           ,&
                         CONP    , CONI    , COND  ,& !Kp ,Ki ,Kd
                         0.d0    , 0.d0    , 0.d0  ,& !Kp2,Ki2,Kd2
                         -1.d30  , 1.d30           ,& !minvar  , maxvar
                         -1.d30  , 1.d30           ,& !mindvar , maxdvar
                         -1.d30  , 1.d30           ,& !minddvar, maxddvar
                         VAR     , DVAR    , DDVAR ,& !VAR , DVAR, DDVAR
                         0.d0    , 0.d0    , 0.d0  ,& !VAR2, DVAR2,DDVAR2
                         YP(1)   , YP1(1)  , YP2(1),&
                         Y (1)   , Y1 (1)  , Y2 (1)   )

      pid(ki)%YP  = Y (1)
      pid(ki)%YP1 = Y1(1)
      pid(ki)%YP2 = Y2(1)
   enddo !nn
      ki            = (h-1)*2
      theta_yaw (h) = pid(ki+1)%YP1
      theta_tilt(h) = pid(ki+2)%YP1
   enddo !h



!--- Inverse Colleman - Set IPC, IFC angles
      Bang(:)    = 0.d0
   do nbod       = 1, NBLADE_el
   do h          = 1, nhT
      Bang(nbod) = Bang(nbod) + theta_tilt(h)*dcos(h*psi(nbod)) + theta_yaw(h)*dsin(h*psi(nbod)) !inv Colleman
   enddo !h

!------ IPC - Cyclic (Individual) Pitch Control
      PITCH_INDIV_PP(nbod) = PITCH_INDIV_P(nbod)
      PITCH_INDIV_P (nbod) = PITCH_INDIV  (nbod)
      PITCH_INDIV   (nbod) = min( max( Bang(nbod)*GAIN_PITCH, -MaxPitchAmpl ), MaxPitchAmpl )
      PVEL          (nbod) = ( 3.d0*PITCH_INDIV(nbod) - 4.d0 * PITCH_INDIV_P(nbod) + PITCH_INDIV_PP(nbod)) / (2*DT_el)
      PACC          (nbod) = (      PITCH_INDIV(nbod) - 2.d0 * PITCH_INDIV_P(nbod) + PITCH_INDIV_PP(nbod)) / (  DT_el**2)

      if (IFLAP_typ == 22) cycle


!------ IFC - Cyclic (Individual) Flap Control
      if     ( TIME <  TIME_INI+TD1        ) then; FlapAng(nbod) =  0.d0
      elseif ( TIME >= TIME_INI+TD1 .and. &
               TIME <  TIME_INI+TD2        ) then; FlapAng(nbod) = (-Bang(nbod)*GAIN_FLAP+BetaF_accel)*(TIME-TIME_INI-TD1)/(TD2-TD1)
      elseif ( TIME >= TIME_INI+TD2        ) then; FlapAng(nbod) = (-Bang(nbod)*GAIN_FLAP+BetaF_accel)
      endif

!------ Flap Actuator (time lag usually ~0.1s)
      VAR         = FlapAng(nbod)
      TVLAG       = timelag(nbod)%TVLAG
      YP  (1)     = timelag(nbod)%YP
      YP1 (1)     = timelag(nbod)%YP1
      YP2 (1)     = timelag(nbod)%YP2

!     call ACTUATOR_1st (VAR, TVLAG, DT   , MinAng    , MaxAng    , MaxVel    , YP   , YP1   , YP2   )
      call ACTUATOR_1st (VAR, TVLAG, DT_el, MinFlapAng, MaxFlapAng, MaxFlapVel, YP(1), YP1(1), YP2(1))

      timelag(nbod)%VARP = VAR
      timelag(nbod)%YP   = YP (1)
      timelag(nbod)%YP1  = YP1(1)
      timelag(nbod)%YP2  = YP2(1)
      FlapAng(nbod)      = YP1(1)
      FlapVel(nbod)      = YP2(1)

      do i = 1,NSTRIP
         if (IFLAPLOC(nbod, i) == 1) then
            FLCONTR(nbod,i) = FlapAng(nbod)*R2D !FLCONTR [deg]
         endif
      enddo ! NSTRIP
   enddo ! nbod


!--- output vars
!--- Flap Angle, Velocity and Acceleration [deg,deg/s,deg/s^2]
   if (GAIN_FLAP > 1.d-10) then
#    if   ASCII == 1
      open(23,file='control_IFC.dat',access='append')
       write(23,'(1000e15.5)')                                 &
#    elif ASCII == 0
      open(23,file='control_IFC.bin',access='append',form='UNFORMATTED')
       write(23              )                                 &
#    endif
         sngl(TIME                 )                          ,& !1
        (sngl(psi        (nbod)*R2D)                          ,& !2,5, 8
         sngl(FlapAng    (nbod)*R2D)                          ,& !3,6, 9
         sngl(FlapVel    (nbod)*R2D)                          ,& !4,7,10     
                                         nbod = 1, NBLADE_el)
      close(23)
   endif


!--- Individual Pitch Angle, Velocity and Acceleration [deg,deg/s,deg/s^2]
   if (GAIN_PITCH > 1.d-10) then
#    if   ASCII == 1
      open(23,file='control_IPC.dat',access='append')
       write(23,'(1000e15.5)')                                 &
#    elif ASCII == 0
      open(23,file='control_IPC.bin',access='append',form='UNFORMATTED')
       write(23              )                                 &
#    endif
         sngl(TIME                 )                          ,& !1
        (sngl(psi        (nbod)*R2D)                          ,& !2, 6, 10
         sngl(PITCH_INDIV(nbod)*R2D)                          ,& !3, 7, 11
         sngl(PVEL       (nbod)*R2D)                          ,& !4, 8, 12
         sngl(PACC       (nbod)*R2D)                          ,& !5, 9, 13
                                         nbod = 1, NBLADE_el)
      close(23)
   endif


!--- BandStop Filter Input/Output , PID Output, Flap Actuator Input/Output
   if (GAIN_FLAP > 1.d-10 .or. GAIN_PITCH > 1.d-10) then
    do    h  = 1, nhT
      call INT_2_CHAR ( 2, CNUM, h )
#    if   ASCII == 1
      open(23,file='control_IPCIFC_'//CNUM(1)//CNUM(2)//'.dat',access='append')
       write(23,'(1000e15.5)')                                 &
#    elif ASCII == 0
      open(23,file='control_IPCIFC_'//CNUM(1)//CNUM(2)//'.bin',access='append',form='UNFORMATTED')
       write(23              )                                 &
#    endif
         sngl(TIME                 )                          ,& !1
        (sngl(filter(kf)%XXP       )                          ,& !2 , 6 , 10, 14   input : 1:Myaw       , 3:Mtilt
         sngl(filter(kf)%XXP1      )                          ,& !3 , 7 , 11, 15   input : 1:d/dt(Myaw) , 3:d/dt(Mtilt)
         sngl(filter(kf)%YYP       )                          ,& !4 , 8 , 12, 16   output: filtered moment
         sngl(filter(kf)%YYP1      )                          ,& !5 , 9 , 13, 17
!        sngl(filter(kf)%YP (1:N)  )                          ,&
!        sngl(filter(kf)%YP1(1:N)  )                          ,&
!        sngl(filter(kf)%YP2(1:N)  )                          ,&
                                    kf = (h-1)*4+1,(h-1)*4+4 ),&
        (sngl(pid   (ki)%YP        )                          ,& !18, 21
         sngl(pid   (ki)%YP1       )                          ,& !19, 22           output: pid(1)%YP1=-theta_yaw, pid(2)%YP1=theta_tilt
         sngl(pid   (ki)%YP2       )                          ,& !20, 23
                                    ki = (h-1)*2+1,(h-1)*2+2 ),&
        (sngl(timelag(nbod)%VARP   )                          ,& !24, 28, 32       input : blade angle after inv Colleman
         sngl(timelag(nbod)%YP     )                          ,& !25, 29, 33  
         sngl(timelag(nbod)%YP1    )                          ,& !26, 30, 34       output:timelag(nbod)%YP1 = FlapAng(nbod)
         sngl(timelag(nbod)%YP2    )                          ,& !27, 31, 35   
                                    nbod = 1, NBLADE_el*(2-h)) !in order to write only for the 1st harmonic
      close(23)
    enddo !h
   endif


 END Subroutine CONTROL_FLAP_CYCL
!----------------------------------------------------------------------
!
!
! PI for flap control based on hub acceleration
!
!----------------------------------------------------------------------
 Subroutine CONTROL_FLAP_hubaccel
       
   Use Craft
   Use Cbeam
   Use Foilfs_mod
   Use Ctrl_mod

   implicit none

   real (8) :: xfl(3), dxfl(3), ddxfl(3), dxlfl(3), ddxlfl(3)
   real (8) :: AKSI0(3)
   integer  :: nbsub, nb_el, nel, nbod, i,j,nn

   real(8)  :: GAIN   , CONP  , CONI  , COND
   real(8)  :: VAR    , DVAR  , DDVAR
   real(8)  :: Y      , Y1    , Y2    
   real(8)  :: YP     , YP1   , YP2   

   integer, save    :: initialize=0

   real (8) :: accel1, accel2, accelF

   integer, parameter :: order = 2 , N = 2*order
   real(8),save       :: Y_acc (N), Y1_acc (N) , Y2_acc (N)
   real(8)            :: YP_acc(N), YP1_acc(N) , YP2_acc(N)
   real(8)            :: A(N,N), B(N), C(N), D

   type      :: type_pid
     real(8) :: YP=0.     , YP1=0.   , YP2=0.               ! PID previous time vars
     real(8) :: GAIN      , CONP     , CONI     , COND      ! PID parameters
   END type type_pid

   type      :: type_filter
     real(8) :: A (N,N)   , B  (N)   , C  (N)   , D         ! filter parameters
   END type type_filter

   type (type_pid    ), save :: pid    (3)
   type (type_filter ), save :: filter (1)                 


!--- Initialization
   if (initialize==0) then

      open (101,file='flap_acc.inp')
       read(101,*)
      do i  = 1, N
       read(101,*) ( filter(1)%A(i,j) , j = 1,N )
      enddo
       read(101,*) ( filter(1)%B(j)   , j = 1,N )
       read(101,*) ( filter(1)%C(j)   , j = 1,N )
       read(101,*)   filter(1)%D
      do nn = 1, 3 !1:flap , 2:gen, 3:pitch
       read(101,*)
       read(101,*) pid(nn)%CONP, pid(nn)%CONI, pid(nn)%COND
       read(101,*)
       read(101,*) pid(nn)%GAIN     
      enddo

      close(101)

!------ filter var for acceleration
      Y_acc (:)=0.d0; Y1_acc (:)=0.d0; Y2_acc (:)=0.d0;

      initialize=1
   endif !initialize


!--- Tower base Acceleration
   nbod             = NBLADE_el + 2
   nbsub            = 1
   nb_el            = body(nbod)%NBODTGLB_el(nbsub)
   nel              = 1
   AKSI0(:)         = 0.d0;

   call beam_kinematics (nbod,nbsub,nel,AKSI0,xfl,dxfl,ddxfl,dxlfl,ddxlfl)

   accel1           = ddxlfl(1)
   if (ICASE_el<3) &
   accel1           = 0.d0


!--- Tower top Acceleration
   nbod             = NBLADE_el + 2
   nbsub            = body   (nbod)%NBODSUB_el
   nb_el            = body   (nbod)%NBODTGLB_el(nbsub)
   nel              = subbody(nb_el)%NTEB_el
   AKSI0(:)         = 0.d0;
   AKSI0(2)         = subbody(nb_el)%HTA_el  (nel+1)

   call beam_kinematics (nbod,nbsub,nel,AKSI0,xfl,dxfl,ddxfl,dxlfl,ddxlfl)

   accel2           = ddxlfl(1)-accel1 !Atop-Abot


!--- Low pass filter acceleration 4th order elliptic filter
   A  (1:N,1:N)     = filter(1)%A  (1:N,1:N);
   B  (1:N    )     = filter(1)%B  (1:N    );
   C  (1:N    )     = filter(1)%C  (1:N    );
   D                = filter(1)%D
      
   YP_acc (:)=Y_acc (:); YP1_acc (:)=Y1_acc (:); YP2_acc (:)=Y2_acc (:);       !update filter solution vector

 !      lowpasstower( A,B,C,D,YP    , YP1    , YP2    , Y    , Y1    , Y2    , XX    , YY    , DT    )
   call lowpasstower( A,B,C,D,YP_acc, YP1_acc, YP2_acc, Y_acc, Y1_acc, Y2_acc, accel2, accelF, DT_el ) !acc  filter


   BetaF_accel      = 0.d0
   GenTrqLSS_accel  = 0.d0
   PITCH_accel      = 0.d0

   if (TIME<dble(NTIME_INI)*DT_el) goto 1

!--- PID 
   do nn            = 1, 3 !1:flap , 2:gen, 3:pitch
      VAR           = accelF
      DVAR          = 0.d0
      DDVAR         = 0.d0
      GAIN          = pid(nn)%GAIN
      CONP          = pid(nn)%CONP
      CONI          = pid(nn)%CONI
      COND          = pid(nn)%COND
      YP            = pid(nn)%YP
      YP1           = pid(nn)%YP1
      YP2           = pid(nn)%YP2

      call PID2_CONTROL( GAIN,     DT_el           ,&
                         CONP    , CONI    , COND  ,& !Kp ,Ki ,Kd
                         0.d0    , 0.d0    , 0.d0  ,& !Kp2,Ki2,Kd2
                         -1.d30  , 1.d30           ,& !minvar  , maxvar
                         -1.d30  , 1.d30           ,& !mindvar , maxdvar
                         -1.d30  , 1.d30           ,& !minddvar, maxddvar
                         VAR     , DVAR    , DDVAR ,& !VAR , DVAR, DDVAR
                         0.d0    , 0.d0    , 0.d0  ,& !VAR2, DVAR2,DDVAR2
                         YP      , YP1     , YP2   ,&
                         Y       , Y1      , Y2       )

      pid(nn)%YP    = Y
      pid(nn)%YP1   = Y1
      pid(nn)%YP2   = Y2
   enddo !nn

   BetaF_accel      = -pid(1)%YP1
   GenTrqLSS_accel  = -pid(2)%YP1
   PITCH_accel      = -pid(3)%YP1

 1 open (63,file='tower_flap.dat',access='append')
      write(63,101) TIME, accel1, accel2, accelF, BetaF_accel*R2D, GenTrqLSS_accel, PITCH_accel*R2D
   close(63)
 101  format (10f15.5)


 END Subroutine CONTROL_FLAP_hubaccel
!----------------------------------------------------------------------
!
! Control library Subroutines
!
!
!----------------------------------------------------------------------
! Performs linear interpolation of a 2variable tabulated function
! YA(MM,NM) given at discrite points X1A(MM),X2A(NM).
! The interpolation is performed @ points X1,X2 and the 
! output value is Y.
! M, N are the number of the given values (M<=MM, N<=NM)
!----------------------------------------------------------------------
 Subroutine BILIN (X1A,X2A,YA,X1,X2,Y,M,N,MM,NM)
!----------------------------------------------------------------------

   implicit none

   integer, intent(in ) :: M,N, MM,NM
   real(8), intent(in ) :: X1A(MM),X2A(NM),YA(MM,NM)
   real(8), intent(in ) :: X1, X2
   real(8), intent(out) :: Y

   real(8) :: GKSI, GHTA !contribution of each nearest point
   integer :: M1  , N1   !result from LOCATE


   call LOCATE1 (MM,M,X1A,X1,M1)
   call LOCATE1 (NM,N,X2A,X2,N1)

   GKSI = (X1-X1A(M1))/(X1A(M1+1)-X1A(M1))
   GHTA = (X2-X2A(N1))/(X2A(N1+1)-X2A(N1))
   Y    = (1.d0-GKSI)*(1.d0-GHTA)*YA(M1  ,N1  ) &
         +      GKSI *(1.d0-GHTA)*YA(M1+1,N1  ) &
         +      GKSI *      GHTA *YA(M1+1,N1+1) &
         +(1.d0-GKSI)*      GHTA *YA(M1  ,N1+1)


 END Subroutine BILIN
!----------------------------------------------------------------------
! Performs linear interpolation of a 3variable tabulated function
! YA(MM,NM,OM) given at discrite points X1A(MM),X2A(NM),X3A(OM).
! The interpolation is performed @ points X1,X2,X3 and the 
! output value is Y.
! M, N, O are the number of the given values (M<=MM, N<=NM, O<=OM)
!----------------------------------------------------------------------
 Subroutine TRILIN (X1A,X2A,X3A,YA,X1,X2,X3,Y,M,N,O,MM,NM,OM)
!----------------------------------------------------------------------

   implicit none

   integer, intent(in ) :: M,N,O, MM,NM,OM
   real(8), intent(in ) :: X1A(MM),X2A(NM),X3A(OM),YA(MM,NM,OM)
   real(8), intent(in ) :: X1, X2, X3
   real(8), intent(out) :: Y

   real(8) :: GKSI, GHTA, GZHTA !contribution of each nearest point
   integer :: M1  , N1  , O1    !result from LOCATE


   call LOCATE1 (MM,M,X1A,X1,M1)
   call LOCATE1 (NM,N,X2A,X2,N1)
   call LOCATE1 (OM,O,X3A,X3,O1)

   GKSI = (X1-X1A(M1))/(X1A(M1+1)-X1A(M1))
   GHTA = (X2-X2A(N1))/(X2A(N1+1)-X2A(N1))
   GZHTA= (X3-X3A(O1))/(X3A(O1+1)-X3A(O1))

   Y    = (1.d0-GKSI) * (1.d0-GHTA) * (1.d0-GZHTA)  *  YA(M1  ,N1  ,O1  ) &
         +      GKSI  * (1.d0-GHTA) * (1.d0-GZHTA)  *  YA(M1+1,N1  ,O1  ) &
         +      GKSI  *       GHTA  * (1.d0-GZHTA)  *  YA(M1+1,N1+1,O1  ) &
         +(1.d0-GKSI) *       GHTA  * (1.d0-GZHTA)  *  YA(M1  ,N1+1,O1  ) &
         +(1.d0-GKSI) * (1.d0-GHTA) *       GZHTA   *  YA(M1  ,N1  ,O1+1) &
         +      GKSI  * (1.d0-GHTA) *       GZHTA   *  YA(M1+1,N1  ,O1+1) &
         +      GKSI  *       GHTA  *       GZHTA   *  YA(M1+1,N1+1,O1+1) &
         +(1.d0-GKSI) *       GHTA  *       GZHTA   *  YA(M1  ,N1+1,O1+1)


 END Subroutine TRILIN
!----------------------------------------------------------------------
! Same as TRILIN, but chooses the shortest azimuthial path
! Proper for interpolation of dpsi(1) and dpsi(2)
!----------------------------------------------------------------------
 Subroutine TRILIN2 (X1A,X2A,X3A,YA,X1,X2,X3,Y,M,N,O,MM,NM,OM)
!----------------------------------------------------------------------

   implicit none

   integer, intent(in   ) :: M,N,O, MM,NM,OM
	real(8), intent(in   ) :: X1A(MM),X2A(NM),X3A(OM)
   real(8), intent(inout) :: YA(MM,NM,OM)
   real(8), intent(in   ) :: X1, X2, X3
   real(8), intent(  out) :: Y

   real(8), parameter :: PI = acos(-1.d0)

   real(8) :: GKSI, GHTA, GZHTA !contribution of each nearest point
   integer :: M1  , N1  , O1    !result from LOCATE

   real(8) :: YA1(2),YA2(2),YA3(2)
   integer :: i


   call LOCATE1 (MM,M,X1A,X1,M1)
   call LOCATE1 (NM,N,X2A,X2,N1)
   call LOCATE1 (OM,O,X3A,X3,O1)

   GKSI = (X1-X1A(M1))/(X1A(M1+1)-X1A(M1))
   GHTA = (X2-X2A(N1))/(X2A(N1+1)-X2A(N1))
   GZHTA= (X3-X3A(O1))/(X3A(O1+1)-X3A(O1))


!--------------------------------------------------------------------------------------
   do i = 0,1
      if (  abs( YA(M1,N1+i,O1) - YA(M1+1,N1+i,O1) )  >  PI  ) then
         if ( YA(M1,N1+i,O1) < YA(M1+1,N1+i,O1) ) then
            YA(M1  ,N1+i,O1) = YA(M1  ,N1+i,O1) + 2.d0*PI
         else
            YA(M1+1,N1+i,O1) = YA(M1+1,N1+i,O1) + 2.d0*PI
         endif
      endif
   enddo

   YA1(1:2) = (1.d0-GKSI)*YA( M1 , N1:N1+1 , O1 ) + GKSI*YA( M1+1 , N1:N1+1 , O1 )

   do i = 1,2
      if (  abs(YA1(i)) > 2.d0*PI )  YA1(i) = mod(abs(YA1(i)),2.d0*PI)*YA1(i)/abs(YA1(i))*2.d0*PI
   enddo
!--------------------------------------------------------------------------------------


!--------------------------------------------------------------------------------------
   do i = 0,1
      if (  abs( YA(M1,N1+i,O1+1) - YA(M1+1,N1+i,O1+1) )  >  PI  ) then
         if ( YA(M1,N1+i,O1+1) < YA(M1+1,N1+i,O1+1) ) then
            YA(M1  ,N1+i,O1+1) = YA(M1  ,N1+i,O1+1) + 2.d0*PI
         else
            YA(M1+1,N1+i,O1+1) = YA(M1+1,N1+i,O1+1) + 2.d0*PI
         endif
      endif
   enddo

   YA2(1:2) = (1.d0-GKSI )*YA( M1 , N1:N1+1 , O1+1 ) + GKSI*YA( M1+1 , N1:N1+1 , O1+1 )

   do i = 1,2
      if (  abs(YA2(i)) > 2.d0*PI )  YA2(i) = mod(abs(YA2(i)),2.d0*PI)*YA2(i)/abs(YA2(i))*2.d0*PI
   enddo
!--------------------------------------------------------------------------------------


!--------------------------------------------------------------------------------------
   do i = 1,2
      if (  abs( YA1(i) - YA2(i) )  >  PI  ) then
         if ( YA1(i) < YA2(i) ) then
            YA1(i) = YA1(i) + 2.d0*PI
         else
            YA2(i) = YA2(i) + 2.d0*PI
         endif
      endif
   enddo

   YA3(1:2) = (1.d0-GZHTA)*YA1(1:2) + GZHTA*YA2(1:2)

   do i = 1,2
      if (  abs(YA3(i)) > 2.d0*PI )  YA3(i) = mod(abs(YA3(i)),2.d0*PI)*YA3(i)/abs(YA3(i))*2.d0*PI
   enddo
!--------------------------------------------------------------------------------------


!--------------------------------------------------------------------------------------
   if (  abs( YA3(1) - YA3(2) )  >  PI  ) then
      if ( YA3(1) < YA3(2) ) then
         YA3(1) = YA3(1) + 2.d0*PI
      else
         YA3(2) = YA3(2) + 2.d0*PI
      endif
   endif

   Y = (1.d0-GHTA)*YA3(1) + GHTA*YA3(2)

   if (  abs(Y) > 2.d0*PI )  Y = mod(abs(Y),2.d0*PI)*Y/abs(Y)*2.d0*PI
!--------------------------------------------------------------------------------------


!   Y    = (1.d0-GKSI) * (1.d0-GHTA) * (1.d0-GZHTA)  *  YA(M1  ,N1  ,O1  ) &
!         +      GKSI  * (1.d0-GHTA) * (1.d0-GZHTA)  *  YA(M1+1,N1  ,O1  ) &
!         +      GKSI  *       GHTA  * (1.d0-GZHTA)  *  YA(M1+1,N1+1,O1  ) &
!         +(1.d0-GKSI) *       GHTA  * (1.d0-GZHTA)  *  YA(M1  ,N1+1,O1  ) &
!         +(1.d0-GKSI) * (1.d0-GHTA) *       GZHTA   *  YA(M1  ,N1  ,O1+1) &
!         +      GKSI  * (1.d0-GHTA) *       GZHTA   *  YA(M1+1,N1  ,O1+1) &
!         +      GKSI  *       GHTA  *       GZHTA   *  YA(M1+1,N1+1,O1+1) &
!         +(1.d0-GKSI) *       GHTA  *       GZHTA   *  YA(M1  ,N1+1,O1+1)


 END Subroutine TRILIN2
!----------------------------------------------------------------------
! Performs linear interpolation of a 4variable tabulated function
! YA(MM,NM,OM,PM) given at discrite points X1A(MM),X2A(NM),X3A(OM),X4A(PM).
! The interpolation is performed @ points X1,X2,X3,X4 and the 
! output value is Y.
! M, N, O, P are the number of the given values (M<=MM, N<=NM, O<=OM, P<=PM)
!----------------------------------------------------------------------
 Subroutine TETRALIN (X1A,X2A,X3A,X4A,YA,X1,X2,X3,X4,Y,M,N,O,P,MM,NM,OM,PM)
!----------------------------------------------------------------------

   implicit none

   integer, intent(in ) :: M,N,O,P, MM,NM,OM,PM
   real(8), intent(in ) :: X1A(MM),X2A(NM),X3A(OM),X4A(PM),YA(MM,NM,OM,PM)
   real(8), intent(in ) :: X1, X2, X3, X4
   real(8), intent(out) :: Y

   real(8) :: GKSI, GHTA, GZHTA, GLAMB !contribution of each nearest point
   integer :: M1  , N1  , O1   , P1    !result from LOCATE


   call LOCATE1 (MM,M,X1A,X1,M1)
   call LOCATE1 (NM,N,X2A,X2,N1)
   call LOCATE1 (OM,O,X3A,X3,O1)
   call LOCATE1 (PM,P,X4A,X4,P1)

   GKSI = (X1-X1A(M1))/(X1A(M1+1)-X1A(M1))
   GHTA = (X2-X2A(N1))/(X2A(N1+1)-X2A(N1))
   GZHTA= (X3-X3A(O1))/(X3A(O1+1)-X3A(O1))
   GLAMB= (X4-X4A(P1))/(X4A(P1+1)-X4A(P1))

   Y    = (1.d0-GKSI) * (1.d0-GHTA) * (1.d0-GZHTA) * (1.d0-GLAMB)   *  YA(M1  ,N1  ,O1  ,P  ) &
         +      GKSI  * (1.d0-GHTA) * (1.d0-GZHTA) * (1.d0-GLAMB)   *  YA(M1+1,N1  ,O1  ,P  ) &
         +      GKSI  *       GHTA  * (1.d0-GZHTA) * (1.d0-GLAMB)   *  YA(M1+1,N1+1,O1  ,P  ) &
         +(1.d0-GKSI) *       GHTA  * (1.d0-GZHTA) * (1.d0-GLAMB)   *  YA(M1  ,N1+1,O1  ,P  ) &
         +(1.d0-GKSI) * (1.d0-GHTA) *       GZHTA  * (1.d0-GLAMB)   *  YA(M1  ,N1  ,O1+1,P  ) &
         +      GKSI  * (1.d0-GHTA) *       GZHTA  * (1.d0-GLAMB)   *  YA(M1+1,N1  ,O1+1,P  ) &
         +      GKSI  *       GHTA  *       GZHTA  * (1.d0-GLAMB)   *  YA(M1+1,N1+1,O1+1,P  ) &
         +(1.d0-GKSI) *       GHTA  *       GZHTA  * (1.d0-GLAMB)   *  YA(M1  ,N1+1,O1+1,P  ) &
         +(1.d0-GKSI) * (1.d0-GHTA) * (1.d0-GZHTA) *       GLAMB    *  YA(M1  ,N1  ,O1  ,P+1) &
         +      GKSI  * (1.d0-GHTA) * (1.d0-GZHTA) *       GLAMB    *  YA(M1+1,N1  ,O1  ,P+1) &
         +      GKSI  *       GHTA  * (1.d0-GZHTA) *       GLAMB    *  YA(M1+1,N1+1,O1  ,P+1) &
         +(1.d0-GKSI) *       GHTA  * (1.d0-GZHTA) *       GLAMB    *  YA(M1  ,N1+1,O1  ,P+1) &
         +(1.d0-GKSI) * (1.d0-GHTA) *       GZHTA  *       GLAMB    *  YA(M1  ,N1  ,O1+1,P+1) &
         +      GKSI  * (1.d0-GHTA) *       GZHTA  *       GLAMB    *  YA(M1+1,N1  ,O1+1,P+1) &
         +      GKSI  *       GHTA  *       GZHTA  *       GLAMB    *  YA(M1+1,N1+1,O1+1,P+1) &
         +(1.d0-GKSI) *       GHTA  *       GZHTA  *       GLAMB    *  YA(M1  ,N1+1,O1+1,P+1)


 END Subroutine TETRALIN
!-----------------------------------------------------------------------
 Subroutine LOCATE1 (nm,n,xa,x,k)
!-----------------------------------------------------------------------

   implicit none

   integer, intent(in ) :: nm,n
   real(8), intent(in ) :: xa(nm), x
   integer, intent(out) :: k

   integer :: khi, klo


   klo = 1
   khi = n
   if (xa(1) < xa(n)) then
 1    if (khi-klo > 1) then
         k=(khi+klo)/2
         if (xa(k) > x) then
            khi=k
         else
            klo=k
         endif
         goto 1
      endif
   else
 2    if (khi-klo > 1) then
         k=(khi+klo)/2
         if (xa(k) > x) then
            klo=k
         else
            khi=k
         endif
         goto 2
      endif
   endif

   k   = klo


 END Subroutine LOCATE1
!----------------------------------------------------------------------
 Subroutine lowpasswind ( YP, YP1, YP2, Y , Y1, Y2, XX, YY, DT )
!----------------------------------------------------------------------

   implicit none

   integer, parameter   :: order = 1 , N = 2*order

   real(8), intent(in ) :: DT
   real(8), intent(in ) :: XX
   real(8), intent(in ) :: YP(N), YP1(N), YP2(N)
   real(8), intent(out) :: Y (N), Y1 (N), Y2 (N)
   real(8), intent(out) :: YY

   real(8) :: AM(N,N), AC(N,N), AK(N,N), AQ(N)
   real(8) :: A (N,N), B (N  ), C (N  ), D

   integer :: i !, j

    real(8) :: damp
 
    damp        = 0.20d0 !0.20 default
    A (1:2,1:2) = 0.d0;   A (1,1) =-damp;    A (2,1) = 1.00d0
    B (1:2    ) = 0.d0;   B (1  ) =-damp**2; B (2  ) = damp
    C (1:2    ) = 0.d0;   C (1  ) = 0.00d0;  C (2  ) = 1.00d0
    D           = 0.d0

!  open(312,file='filters.inp')
!     read(312,*)
!     do i = 1,N
!        read(312,*) ( A(i,j) , j = 1,N )
!     enddo
!     read(312,*)    ( B(j)   , j = 1,N )
!     read(312,*)    ( C(j)   , j = 1,N )
!     read(312,*)      D
!  close(312)

   AK(1:N,1:N) = 0.d0;
   AC(1:N,1:N) =-A(1:N,1:N)
   AM(1:N,1:N) = 0.d0; forall(i = 1:N) AM(i,i) = 1.d0
   AQ(1:N    ) = B(1:N)*XX

!       NEWMARKN (UTP,UTP1,UTP2,UT,UT1,UT2,DT,AM,AC,AK,AQ,N) 
   call NEWMARKN (YP ,YP1 ,YP2 ,Y ,Y1 ,Y2 ,DT,AM,AC,AK,AQ,N)

   YY          = dot_product (C(1:N),Y1(1:N)) + D*XX   !answer output


 END Subroutine lowpasswind
!----------------------------------------------------------------------
 Subroutine lowpasstower ( A,B,C,D, YP, YP1, YP2, Y , Y1, Y2, XX, YY, DT )
!----------------------------------------------------------------------

   implicit none

   integer, parameter   :: order = 2 , N = 2*order

   real(8), intent(in ) :: DT
   real(8), intent(in ) :: XX
   real(8), intent(in ) :: YP(N), YP1(N), YP2(N)
   real(8), intent(out) :: Y (N), Y1 (N), Y2 (N)
   real(8), intent(out) :: YY

   real(8) :: AM(N,N), AC(N,N), AK(N,N), AQ(N)
   real(8) :: A (N,N), B (N  ), C (N  ), D
   integer :: i


   AK(1:N,1:N) = 0.d0;
   AC(1:N,1:N) =-A(1:N,1:N)
   AM(1:N,1:N) = 0.d0; forall(i = 1:N) AM(i,i) = 1.d0
   AQ(1:N    ) = B(1:N)*XX

!       NEWMARKN (UTP,UTP1,UTP2,UT,UT1,UT2,DT,AM,AC,AK,AQ,N) 
   call NEWMARKN (YP ,YP1 ,YP2 ,Y ,Y1 ,Y2 ,DT,AM,AC,AK,AQ,N)

   YY          = dot_product (C(1:N),Y1(1:N)) + D*XX   !answer output


 END Subroutine lowpasstower
!----------------------------------------------------------------------
 Subroutine ACTUATOR_1st (VAR, TVLAG, DT, MinAng, MaxAng, MaxVel, YP, YP1, YP2)
!-----------------------------------------------------------------------

   implicit none

   real(8), intent(in   ) :: VAR   , TVLAG , DT
   real(8), intent(in   ) :: MinAng, MaxAng, MaxVel
   real(8), intent(inout) :: YP   , YP1   ,YP2

   real(8)                :: AM, AC, AK, AQ
   real(8)                :: Y    , Y1    ,Y2
   real(8)                :: BITA, GAMMA, UPRE, U1PRE


!     VAR         = FlapAng(nbod)
!     TVLAG       = timelag(nbod)%TVLAG
!     YP          = timelag(nbod)%YP
!     YP1         = timelag(nbod)%YP1
!     YP2         = timelag(nbod)%YP2

      AM       = TVLAG
      AC       = 1.d0
      AK       = 0.d0
      AQ       = VAR
      
      call Newmark_solve1 (AM, AC , AK , AQ, DT, &
                           YP, YP1, YP2,         &
                           Y , Y1 , Y2            )

!------ Saturate flap velocity and angle based on Newmark
         BITA  = 0.25d0
         GAMMA = 0.50d0
!------ Saturate flap velocity
      if (dabs(Y2) > MaxVel) then
         UPRE  = YP + YP1*DT + YP2*(0.5d0 - BITA )*DT**2
         U1PRE =      YP1    + YP2*(1.0d0 - GAMMA)*DT
         Y2    = MaxVel * dabs(Y2)/Y2
         Y     = UPRE  + BITA  * DT**2  * Y2
         Y1    = U1PRE + GAMMA * DT     * Y2
      endif

!------ Saturate flap angle
      if (Y1 > MaxAng .or. Y1 < MinAng) then
         UPRE  = YP + YP1*DT + YP2*(0.5d0 - BITA )*DT**2
         U1PRE =      YP1    + YP2*(1.0d0 - GAMMA)*DT
         Y1    = max( min( Y1, MaxAng ), MinAng )
         Y2    = (Y1-U1PRE) / (GAMMA*DT)
         Y     = UPRE  + BITA  * DT**2 * Y2
      endif
         YP    = Y
         YP1   = Y1
         YP2   = Y2

!     timelag(nbod)%VARP = VAR
!     timelag(nbod)%YP   = Y
!     timelag(nbod)%YP1  = Y1
!     timelag(nbod)%YP2  = Y2
!     FlapAng(nbod)      = Y1


 END Subroutine ACTUATOR_1st
