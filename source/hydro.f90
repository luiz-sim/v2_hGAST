!----------------------------------------------------------------------
!  1. When Wavemod = 6 and velocities-accelerations are read from
!     file, in order to integrate hydro loads up to current z
!     elevation WheelerStretching is set by default 1. The upper U
!     and A are also interpolated in order to avoid artificial high 
!     frequencies near the free surface.
!
!----------------------------------------------------------------------
 Module Hydro
!----------------------------------------------------------------------
  type :: type_floater
!------------------------
!---- Floater's Properties [+linear Moorings] --
    real(8)              :: AM_floater_structural (6,6)    ! m, m*xc, Ixx+m*xc**2, Izz
    real(8)              :: AK_floater_structural (6,6)    ! Mx,My restoring moments due to gravity
    real(8)              :: AK_hydrostatic        (6,6)    ! pi*rho*grav*R**2, Buoyancy * Hhydro_center
    real(8), allocatable :: AM_added_infinity     (:,:)    ! (6,6*nfloaters_hyd)
    real(8)              :: AK_external           (6,6)    ! mooring
    real(8)              :: AC_external           (6,6)
    real(8), allocatable :: AC_retardation        (:,:)    ! (6,6*nfloaters_hyd)
    real(8)              :: AQ_retardation        (  6)
    real(8)              :: AQ_exciting           (  6)    ! total 1st order excitation force
    real(8)              :: AQ_Newman             (  6)    ! total 2nd order difference-frequency force from Newman's approximation
    real(8)              :: AQ_external           (  6)
    real(8)              :: AQ_damp_quad          (  6)
    real(8)              :: AK_morison            (6,6)
    real(8)              :: AC_morison            (6,6)
    real(8)              :: AQ_morison            (  6)
    real(8)              :: Buoyancy_total
    real(8)              :: Moorings_weight
    real(8)              :: Draft
   
!---- Diffraction --
    real(8), allocatable :: wave_force            (:,:)    ! (no_spectral_comp,6)
    real(8), allocatable :: wave_phase            (:,:)    ! (no_spectral_comp,6)
!---- Drift force for Newman's approximation
    real(8), allocatable :: wave_force_Dft_pos    (:,:)    ! (no_spectral_comp,6)
    real(8), allocatable :: wave_force_Dft_neg    (:,:)    ! (no_spectral_comp,6)
   
!---- Radiation --
    real(8), allocatable :: retardation_function  (:,:,:)  ! (6,6*nfloaters_hyd,retardation_length)
    real(8), allocatable :: velocity_buffer       (:,:  )  ! (6                ,retardation_length)

!---- old Floating module --
    real(8), allocatable :: AM_total              (:,:)    ! (6,6*nfloaters_hyd)
    real(8), allocatable :: AC_total              (:,:)    ! (6,6*nfloaters_hyd)
    real(8)              :: AK_total              (6,6)
    real(8)              ::  Q_total              (6  )
    real(8)              ::  Q_truss              (6  )
    real(8)              :: new_velocity          (6  )
    integer              :: IQfl_active_el        (6  )    ! Switch on-off floater's dofs [0:off, 1:on]
    real(8)              :: UflInit               (6  )    ! Initial floater's 6 dofs

  END type type_floater


!-- Floater structure
 type (type_floater), allocatable, save :: floater   (:)   ! [NFLOATER_el/nfloater_hyd]


 integer,              save :: retardation_length

!-- General Variables --
 integer,              save :: LoadMod         !** LoadMod    (1:Morison, 2:Diffraction/Radiation, 3:Diff/Rad + Viscous Drag,21 & 31: enable Newman's approximation)
 integer,              save :: WaveMod         !** WaveMod    (0:none, 1:regular, 2:JONS, 3:PM, 4:read w-A(w)-fi(w), 5:read w-S(w), 6:read directly z,u,a, 7:white noise)
 integer,              save :: CurrMod         !** CurrMod    (0:none, 1:steady current with z power law )
 integer,              save :: StretchMod      !** StretchMod (0:none, 1:Wheeler)
 integer,              save :: NewmanMod       !** NewmanMod  (0:none, 1:yes    )
 integer,              save :: DragMod         
 integer,              save :: IADDMASSELAST   
 integer,              save :: IBUOYANCY       !** Buoyancy Calculation       (0:no, 1:pressure integration method, 2:Volume method )
 integer,              save :: ITYPE_Damp_ext  !** Floater's External damping (0:no, 1: linear external damping, 2: quadratic damping added to rhs)
 integer,              save :: nfloater_hyd    !** number of floaters

 real(8),              save :: gravity
 real(8),              save :: rho
 real(8),              save :: depth
 real(8),              save :: tide
 real(8),              save :: rhoGrav
 real(8),              save :: depth_tide

 real(8),              save :: Time_ramp
 real(8),              save :: Start_Fact

 real(8),              save :: pi_h
 real(8),              save :: pi2_h
 real(8),              save :: delta_t
 real(8),              save :: time_h

 real(8),              save :: Water_upper_lim
 real(8),              save :: Water_lower_lim
 real(8), parameter         :: kd_deep_lim = 340.d0



!-- Wave, Generated from code --
 integer,              save :: no_spectral_comp
 integer,              save :: no_spectral_comp_deep           !after this wavecomponent --> deep water
 integer,              save :: IRand(2)
 real(8), allocatable, save :: spectral_freq           (:)     !(no_spectral_comp  )  !rad/sec
 real(8), allocatable, save :: spectral_ampl           (:)     !(no_spectral_comp  )  !height/2
 real(8), allocatable, save :: spectral_phase          (:)     !(no_spectral_comp  )
 real(8), allocatable, save :: spectral_wavenum        (:)     !(no_spectral_comp  )
 real(8),              save :: spectral_angle
 real(8),              save :: Hs_sp
 real(8),              save :: Tp_sp
 real(8),              save :: DOmega_sp
 real(8),              save :: Omega_Max_sp


!-- Wave, read from file [Read U, A, p, zeta] --
 integer,              save :: IZMax_wav
 integer,              save :: ITMax_wav
 real(8),              save :: DT_wav
 real(8),              save :: Tperiod_wav
 real(8),              save :: Cvel_wav
 real(8),              save :: T0_wav
 real(8), allocatable, save :: U_wav                 (:,:,:)  !(ITMax_wav, IZMax_wav, 3)
 real(8), allocatable, save :: A_wav                 (:,:,:)  !(ITMax_wav, IZMax_wav, 3)
 real(8), allocatable, save :: P_wav                 (  :,:)  !(ITMax_wav, IZMax_wav   )
 real(8), allocatable, save :: Z_wav                 (    :)  !(           IZMax_wav   )
 real(8), allocatable, save :: T_wav                 (    :)  !(ITMax_wav              )
 real(8), allocatable, save :: Zelev_wav             (    :)  !(ITMax_wav              )


!-- Current Properties --
 real(8),              save :: Current_veloc
 real(8),              save :: Current_polaw
 real(8),              save :: Current_angle


!-- Morison's Equation General Variables --
!qq
 integer,              save :: NX_hyd                         ! number of parts/element for hydrodynamic loads (morison/buoyancy) 
 real(8),              save :: CDrag_mor
 real(8),              save :: CMass_mor


!-- Output results
 real(8), allocatable, save :: RES_hyd  (:,:,:,:)             ! nbod, nel, nx_hyd, vars
 

!---- Special application/ sambrela concept
 real(8),              save :: RC_hyd(3)

!qq - for scaling only the outputs
!real(8), parameter         :: scal_hyd  =   1.d0/45.0d0 !ecn innwind experiment
!real(8), parameter         :: rho_salt  = 1025.0d0
!real(8), parameter         :: rho_sweet =  998.2d0
 real(8), parameter         :: scal_hyd  =   1.d0
 real(8), parameter         :: rho_salt  = 1025.0d0
 real(8), parameter         :: rho_sweet = 1025.0d0


 END Module Hydro
!-----------------------------------------------------------------------
 Subroutine Hydro_main ( it )
!-----------------------------------------------------------------------

 Use Hydro
 Use Cbeam

   implicit none

   integer :: it, nf, j


   if (ICASE_el == 0) return

   time_h = TIME
                         Start_Fact = 1.d0
   if (time_h<Time_ramp) Start_Fact = 0.5d0*(1.d0-dcos(PI*time_h/Time_ramp))


!--- Hydrodynamic Problem for floater's dofs
   if (ICASE_el == 3) then
!------ Set floater's velocity
      do nf = 1, nfloater_hyd
         j  = (nf-1)*6
         floater(nf)% new_velocity (1:6) = UT1_el(NDFBT_el+1+j:NDFBT_el+6+j);
      enddo

      call floater_simulator( 1, it )
      return
   endif


!--- Zero results variables store matrix
   if ( (ICASE_el == 2).or.(ICASE_el == 4) ) &
   Res_hyd (:,:,:,:) = 0.d0;


 END Subroutine Hydro_main
!-----------------------------------------------------------------------
 Subroutine Hydro_Init (rho_out, depth_out)
!-----------------------------------------------------------------------

 Use Hydro
 Use Cbeam
 Use Paths

   implicit none

   real(8) :: rho_out, depth_out
   integer :: nbod, nbsub, nb_el, nel, n1, n2, nelM, nf, j, N


   if (ICASE_el == 0) return


   open (1,file=trim(file_hydro)) !'hydro.inp'

   read(1,*)         !*** ---Hydro Module Input--- ***
   read(1,*)         !** LoadMod (1:Morison, 2:Diffraction/Radiation, 3:Diff/Rad + Viscous Drag)
   read(1,*) LoadMod
   read(1,*)         !** WaveMod (0:none, 1:regular, 2:JONS, 3:PM, 4:read w-A(w)-fi(w), 5:read w-S(w), 6:read directly z,u,a)
   read(1,*) WaveMod
   read(1,*)         !** CurrMod (0:none, 1...    )
   read(1,*) CurrMod
   read(1,*)         !** StretchMod (0:none, 1:Wheeler)
   read(1,*) StretchMod

   read(1,*)         !** Buoyancy Calculation (0:no, 1:yes)
   read(1,*) IBUOYANCY
   read(1,*)         !** Elastic Added Mass in Morison's equaiton (0:no, 1:yes)
   read(1,*) IADDMASSELAST

   read(1,*)         !** Morison's Coeff
   read(1,*) CDrag_mor
   read(1,*) CMass_mor
   read(1,*) NX_hyd

   read(1,*)         !** Environment
   read(1,*) gravity          ![m /s^2]
   read(1,*) rho              ![kg/m^3]
   read(1,*) depth            ![m     ]
   read(1,*) tide             ![m     ]
   read(1,*) Water_upper_lim
   read(1,*) Water_lower_lim

   read(1,*)         !** Spectral Properties [not always needed]
   read(1,*) Time_ramp        ![sec   ]
   read(1,*) spectral_angle   ![deg   ]
   read(1,*) Hs_sp            ![m     ]
   read(1,*) Tp_sp            ![sec   ]
   read(1,*) DOmega_sp        ![rad/s ]
   read(1,*) Omega_Max_sp     ![rad/s ]
   read(1,*) IRand(1)
   read(1,*) IRand(2)

   read(1,*)         !** Current  Properties [not always needed]
   read(1,*) Current_veloc    ![m/sec ]
   read(1,*) Current_polaw    ![ -    ]
   read(1,*) Current_angle    ![deg   ]


   if ( ( WaveMod == 6 ).and.( tide /= 0.d0 ) ) then
      write(* ,*) 'Error, for WaveMod 6 tide must be Zero!'
      write(10,*) 'Error, for WaveMod 6 tide must be Zero!'
      stop
   endif


   delta_t         = DT_el
   depth_tide      = depth + tide
   rho_out         = rho
   depth_out       = depth
   nfloater_hyd    = NFLOATER_el !** number of floaters

   spectral_angle  = spectral_angle / R2D
   Current_angle   = Current_angle  / R2D
   pi_h            = dacos(-1.d0)
   pi2_h           = 2d0*pi_h
   rhoGrav         = rho*gravity

!----------------------- Modify LoadMod and NewmanMod for Newman's approximation
                       NewmanMod = 0
   if ( LoadMod == 21) NewmanMod = 1
   if ( LoadMod == 31) NewmanMod = 1
   if ( LoadMod == 21) LoadMod   = 2
   if ( LoadMod == 31) LoadMod   = 3

                       DragMod = 0
   if ( LoadMod == 3 ) DragMod = 1


!--- set: no_spectral_comp, spectral_freq, spectral_ampl, spectral_phase
   call wave_init

   close(1)

!--- calculate Wavenum or read Non linear waves
   call morison_init


   if (LoadMod==1) then
      if (NBODTWT_el==0) then
          nel =maxval ( subbody(         1:NBODT_el)%NTEB_el )
      else
          nel =maxval ( subbody(NBODTWT_el:NBODT_el)%NTEB_el )
      endif
      write(*,*) 'max nel for "morison" body',nel
      allocate ( RES_hyd  (NBODBTWT_el:NBODBT_el,nel,NX_hyd,31) )            ! nbod, nel, nx_hyd, vars
   endif


   if     (ICASE_el == 1) then
!qq   call Monopile_Diffr_init
      if ( LoadMod /= 1 ) then
         write(* ,*) 'Error, Diffraction for monopile not supported'
         write(10,*) 'Error, Diffraction for monopile not supported'
         stop
      endif

      IBUOYANCY  = 0
   elseif (ICASE_el == 2) then
      if (LoadMod /= 1) then
         write(* ,*) 'Error, Diffraction for tripod/jacket not supported'
         write(10,*) 'Error, Diffraction for tripod/jacket not supported'
         stop
      endif
   elseif (ICASE_el == 3) then
      IBUOYANCY  = 0
      StretchMod = 0

      do nf = 1, nfloater_hyd
         j  = (nf-1)*6
         floater(nf)% new_velocity (1:6) = UT1_el(NDFBT_el+1+j:NDFBT_el+6+j);
      enddo

      call floater_simulator( 0, 0 )
   elseif (ICASE_el == 4) then
!------ Just read the floater.inp file for the main matrices [linear stiffness-damping]
         N  = 6 !*nfloater_hyd

      do nf = 1, 1 !nfloater_hyd
         Allocate ( floater(nf)% AC_retardation   (6, N) )
         Allocate ( floater(nf)% AC_total         (6, N) )
         Allocate ( floater(nf)% AM_added_infinity(6, N) )
         Allocate ( floater(nf)% AM_total         (6, N) )
      enddo

      call floater_init

!qq Valid only for 1 floater
      nf = 1

      floater(nf)% AM_total (1:6,1:6) =  floater(nf)% AM_floater_structural (1:6,1:6)         &
                                        +floater(nf)% AM_added_infinity     (1:6,1:6)
      floater(nf)% AC_total (1:6,1:6) =  floater(nf)% AC_external           (1:6,1:6)
      floater(nf)% AK_total (1:6,1:6) =  floater(nf)% AK_external           (1:6,1:6)         &
                                        +floater(nf)% AK_hydrostatic        (1:6,1:6)         &
                                        +floater(nf)% AK_floater_structural (1:6,1:6)
      floater(nf)% Q_total  (    1:6) =  floater(nf)% AQ_external           (    1:6)
!------ addition of total buoyancy, floater's weight and Moorings_weight
      floater(nf)% Q_total  (      3) =  floater(nf)% Q_total(3)                              &
                                        +floater(nf)% Buoyancy_total                          &
                                        -floater(nf)% AM_floater_structural (3  ,3  )*gravity &
                                        -floater(nf)% Moorings_weight
!------ Add the External damping
      if     (ITYPE_Damp_ext==2) then
         write(* ,*)'type2 external damping not supported for flexible floater atm'
         write(10,*)'type2 external damping not supported for flexible floater atm'
         stop
      elseif (ITYPE_Damp_ext==0) then
         floater(nf)% AC_total (1:6,1:6) = 0.d0
      endif
   endif !ICASE_el


   if     ( StretchMod == 0.or.WaveMod==0 ) then
      Water_upper_lim = 0.d0
      Water_lower_lim = 0.d0
   elseif ( StretchMod == 1 ) then
     if     (WaveMod==1           ) then
      Water_upper_lim = Hs_sp/2.d0
      Water_lower_lim =-Hs_sp/2.d0
     elseif (WaveMod==2.or.WaveMod==3) then
      Water_upper_lim = Hs_sp/2.d0 * 2.d0 !1.86
      Water_lower_lim =-Hs_sp/2.d0 * 2.d0 !1.86
     endif
   endif
!------ Consider the tide level
      Water_upper_lim = Water_upper_lim + tide
      Water_lower_lim = Water_lower_lim + tide

   write(10,*)'Water_upper_lim',Water_upper_lim
   write(10,*)'Water_lower_lim',Water_lower_lim



!------ For Tower: transfer vars from hydro to substr_mor
   if ( ((ICASE_el    == 1).or.      &
         (ICASE_el    == 2).or.      &
         (ICASE_el    == 4)    ).and.&
        ( NBODBTWT_el >  0)            ) then

      nbod        = NBODBTWT_el
      nbsub       = 1
      nb_el       = body(nbod)%NBODTGLB_el(nbsub)

      do nel = 1, subbody(nb_el)%NTEB_el
         n1                       = subbody(nb_el)%INODNEL_el(nel,1)
         n2                       = subbody(nb_el)%INODNEL_el(nel,2)
         substr_mor(n1)%Cm_el     = CMass_mor
         substr_mor(n2)%Cm_el     = CMass_mor
         substr_mor(n1)%Cd_el     = CDrag_mor
         substr_mor(n2)%Cd_el     = CDrag_mor
      enddo
   endif


!--- Allocate results variables restore matrix
   if ( (ICASE_el == 2).or.(ICASE_el == 4) ) then
      nelM=0
      do nbod  = NBODBTWT_el+1, NBODBT_el
         nbsub = 1
         nb_el = body(nbod)%NBODTGLB_el(nbsub)
         nelM  = max (subbody(nb_el)%NTEB_el, nelM)
      enddo
!qq
      Deallocate( RES_hyd )
      Allocate  ( RES_hyd(NBODBTWT_el:NBODBT_el, nelM, NX_hyd, 31 ) )
   endif        !        nbod                 , nel , nhyd  , vars


 END Subroutine Hydro_Init
!----------------------------------------------------------------------
 Subroutine morison_init
!----------------------------------------------------------------------

 Use Hydro

   implicit none

   real (8) :: freq, Wavenum
   integer  :: i, ideep


   if (WaveMod == 6) then

      call read_UAPZ

      return
   endif


   open (1,file='Wavenum.out')

   Wavenum               = 0.1d0    !initial value for 1st wave component
   ideep                 = 0
   no_spectral_comp_deep = no_spectral_comp + 1   !if no deep water wave is present

   do i    = 1, no_spectral_comp
      freq = spectral_freq(i)

      call WAVE_NUM ( freq, depth_tide, gravity, Wavenum, kd_deep_lim )

      spectral_wavenum (i) = Wavenum

      write(1,101) freq,Wavenum, pi2_h / Wavenum, Wavenum*depth_tide

      if ( (ideep == 0).and.(Wavenum*depth_tide >= kd_deep_lim) ) then          !-- 520.d0
         ideep = 1
         no_spectral_comp_deep = i
         write( *,*) 'Deep water components from:',no_spectral_comp_deep
         write(10,*) 'Deep water components from:',no_spectral_comp_deep
      endif
   enddo

   close(1)

101 format ( 150f18.6 )


 END Subroutine morison_init
!----------------------------------------------------------------------
 Subroutine read_UAPZ
!----------------------------------------------------------------------

 Use Hydro
 Use Cbeam
 Use Paths

   implicit none

   real(8)      :: c1
   integer      :: it, iz, I,j
   Character    :: CNUM(5)
   Character*80 :: outfil


   open (2,file=trim(dir_wave)//'kinematics.inp')
   open (3,file=trim(dir_wave)//'z.inp')
   open (4,file=trim(dir_wave)//'surface.inp')

   read (2,*) IZMax_wav     !different Z points
   read (2,*) ITMax_wav     !time steps
   read (2,*) DT_wav        !time step     [sec]
   read (2,*) Tperiod_wav   !wave period   [sec]
   read (2,*)    Cvel_wav   !wave velocity [m/sec]
   read (2,*) T0_wav        !time where wave is located at x=0m [sec]

   write(10,*)
   write(10,*) 'different Z points    ',   IZMax_wav
   write(10,*) 'time steps            ',   ITMax_wav
   write(10,*) 'time step     [sec]   ',      DT_wav
   write(10,*) 'wave period   [sec]   ', Tperiod_wav
   write(10,*) 'wave velocity [m/sec] ',    Cvel_wav
   write(10,*) 'T0 at x=0m    [sec]   ',      T0_wav
   write(10,*)

   Allocate ( U_wav    ( ITMax_wav, IZMax_wav, 3) ,&
              A_wav    ( ITMax_wav, IZMax_wav, 3) ,&
              P_wav    ( ITMax_wav, IZMax_wav   ) ,&
              Z_wav    (            IZMax_wav   ) ,&
              T_wav    ( ITMax_wav              ) ,&
              Zelev_wav( ITMax_wav              )    )


   do it=1,ITMax_wav
   do iz=1,IZMax_wav

      read(2,*) U_wav(it,iz,1),&
                U_wav(it,iz,2),&
                U_wav(it,iz,3),&
                A_wav(it,iz,1),&
                A_wav(it,iz,2),&
                A_wav(it,iz,3),&
                P_wav(it,iz  )

   enddo
   enddo


   do iz=1,IZMax_wav

      read(3,*) Z_wav(iz)

   enddo


   read(4,*) !Time [s]      Sea surface elevation [m]
   do it=1,ITMax_wav

      read(4,*) T_wav(it), Zelev_wav(it)

   enddo


   close(2)
   close(3)
   close(4)


!-- correct the inputs  near free surface, when values after z_elev
!-- are given zero and cause noise after the interpolation

   do it=1,ITMax_wav

     call LOCATE ( IZMax_wav, IZMax_wav, Z_wav, Zelev_wav(it), I)

     do j=I+1,IZMax_wav

     !  sl =(yi(n1)-yi(n1-1))/(xi(n1)-xi(n1-1))
     !  y  = yi(n1-1)+sl*(x-xi(n1-1))

        c1 =(U_wav(it,j-1,1)-U_wav(it,j-2,1)) / (Z_wav(j-1)-Z_wav(j-2))
        U_wav(it,j,1) = U_wav(it,j-2,1) + c1 * (Z_wav(j)-Z_wav(j-2))
        c1 =(U_wav(it,j-1,2)-U_wav(it,j-2,2)) / (Z_wav(j-1)-Z_wav(j-2))
        U_wav(it,j,2) = U_wav(it,j-2,2) + c1 * (Z_wav(j)-Z_wav(j-2))
        c1 =(U_wav(it,j-1,3)-U_wav(it,j-2,3)) / (Z_wav(j-1)-Z_wav(j-2))
        U_wav(it,j,3) = U_wav(it,j-2,3) + c1 * (Z_wav(j)-Z_wav(j-2))

        c1 =(A_wav(it,j-1,1)-A_wav(it,j-2,1)) / (Z_wav(j-1)-Z_wav(j-2))
        A_wav(it,j,1) = A_wav(it,j-2,1) + c1 * (Z_wav(j)-Z_wav(j-2))
        c1 =(A_wav(it,j-1,2)-A_wav(it,j-2,2)) / (Z_wav(j-1)-Z_wav(j-2))
        A_wav(it,j,2) = A_wav(it,j-2,2) + c1 * (Z_wav(j)-Z_wav(j-2))
        c1 =(A_wav(it,j-1,3)-A_wav(it,j-2,3)) / (Z_wav(j-1)-Z_wav(j-2))
        A_wav(it,j,3) = A_wav(it,j-2,3) + c1 * (Z_wav(j)-Z_wav(j-2))

        c1 =(P_wav(it,j-1  )-P_wav(it,j-2  )) / (Z_wav(j-1)-Z_wav(j-2))
        P_wav(it,j  ) = P_wav(it,j-2  ) + c1 * (Z_wav(j)-Z_wav(j-2))
     enddo

   enddo !it

!- uncomment to write read waves for checking

!-- depth
!! do it=1,ITMax_wav
!!   call INT_2_CHAR ( 5, CNUM, it )
!!   outfil = 'zwave'//CNUM(1)//CNUM(2)//CNUM(3)//CNUM(4)//CNUM(5)//'.dat'

!!   open (23,file=outfil, access='append')
!!   do iz=1,IZMax_wav

!!      write(23,8) Z_wav(iz)     ,&
!!                  U_wav(it,iz,1),&
!!                  U_wav(it,iz,2),&
!!                  U_wav(it,iz,3),&
!!                  A_wav(it,iz,1),&
!!                  A_wav(it,iz,2),&
!!                  A_wav(it,iz,3),&
!!                  P_wav(it,iz  ),&
!!                  Zelev_wav(it )

!!   enddo

!!   close(23)

!! enddo !it

!-- time
   do iz=1,IZMax_wav

     call INT_2_CHAR ( 5, CNUM, iz )
     outfil = 'twave'//CNUM(3)//CNUM(4)//CNUM(5)//'.dat'

     open (23,file=outfil, access='append')

     do it=1,ITMax_wav
        write(23,8)     &
         T_wav(it)     ,&
         U_wav(it,iz,1),&
         U_wav(it,iz,2),&
         U_wav(it,iz,3),&
         A_wav(it,iz,1),&
         A_wav(it,iz,2),&
         A_wav(it,iz,3),&
         P_wav(it,iz  ),&
         Zelev_wav(it )
     enddo

     close(23)

   enddo !iz

   8 format (150f21.8)

!!!write(*,*) U_wav(ITMax_wav,IZMax_wav,3)
!!!write(*,*) A_wav(ITMax_wav,IZMax_wav,3)
!!!write(*,*) P_wav(ITMax_wav,IZMax_wav  )
!!!write(*,*) T_wav(ITMax_wav)
!!!write(*,*) Z_wav(IZMax_wav)
!!!write(*,*) Zelev_wav(ITMax_wav)


 END Subroutine read_UAPZ
!----------------------------------------------------------------------
 Subroutine wave_init
!----------------------------------------------------------------------

 Use Hydro

   implicit none

   integer :: i, no_inp_wave_freq, n
   real(8) :: delta_om, omega0, R2D
   real(8), allocatable :: inp_spectral_freq(:)
   real(8), allocatable :: inp_spectral_dens(:)


   open (25,file='wave_sp.out')

   write(25,*) 'Wavemod =',Wavemod



   if ( Wavemod == 0 ) then
!------ No Wave
      no_spectral_comp = 0
      return
   elseif ( Wavemod == 1 ) then
!------ Regular Wave
      no_spectral_comp = 1
   elseif ( Wavemod == 2  .or.  Wavemod == 3 .or. Wavemod == 7) then
!------ Irregular Wave (2: Jonswap, 3: Pirson-Moskovitz, 7: white noise)
      no_spectral_comp = int(Omega_Max_sp/DOmega_sp) !!+1 CHECK
   elseif ( Wavemod == 4 ) then
!------ Read Regular Waves (Freq, Ampl, Phase)
      read(1,*) !** WaveMod 4 or 5 inputs 
      read(1,*)  no_spectral_comp
   elseif ( Wavemod == 5 ) then
!------ Read Irregular Waves (Freq, Spectral Dens)
      read(1,*) !** WaveMod 4 or 5 inputs 
      read(1,*)  no_spectral_comp
   elseif ( Wavemod == 6 ) then
!------ Non Linear Wave  (read u,a and apply to Morison's eq.)
      no_spectral_comp = 0
      StretchMod       = 1
      return
   endif !Wavemod


   write (* ,*) 'no of wave components', no_spectral_comp
   write (10,*) 'no of wave components', no_spectral_comp


   Allocate ( spectral_freq    (no_spectral_comp) )
   Allocate ( spectral_ampl    (no_spectral_comp) )
   Allocate ( spectral_phase   (no_spectral_comp) )
   Allocate ( spectral_wavenum (no_spectral_comp) )


   if     ( Wavemod == 1 ) then
      spectral_freq (1) = pi2_h/Tp_sp
      spectral_ampl (1) = Hs_sp / 2d0
      spectral_phase(1) = 0.d0
   elseif ( Wavemod == 2  .or.  Wavemod == 3 .or. Wavemod == 7) then
      omega0   = 0.d0
      delta_om = DOmega_sp

      call random_wave ( delta_om, omega0 )
      call spectrum
   elseif ( Wavemod == 4 ) then
      R2D = 180.d0/pi_h
      read(1,*) !**N  w [rad/s] H/2 [m]   Phase [deg]   K [1/m] 

      do i=1,no_spectral_comp           !High/2
         read(1,*) n, spectral_freq(i),spectral_ampl(i),spectral_phase(i)
         spectral_phase(i) = spectral_phase(i) / R2D
      enddo
   elseif ( Wavemod == 5 ) then
      read(1,*) no_inp_wave_freq

      Allocate ( inp_spectral_freq(no_inp_wave_freq) ,& 
                 inp_spectral_dens(no_inp_wave_freq)    )

      read(1,*) !**N     w [rad/sec]   Spectral density [m^2*sec/rad]

      do i=1,no_inp_wave_freq
         read(1,*) n, inp_spectral_freq(i), inp_spectral_dens(i)
      enddo

      omega0   =  inp_spectral_freq(1)
      delta_om = (inp_spectral_freq(no_inp_wave_freq) - inp_spectral_freq(1))/ dble (no_spectral_comp)

      call random_wave ( delta_om, omega0 )
      call int_read_spec ( inp_spectral_freq, inp_spectral_dens, no_inp_wave_freq )

      Deallocate ( inp_spectral_freq, inp_spectral_dens )
   endif !Wavemod


   write(25,*)' Individual Wave Components'
   write(25,*)'#         Frequency    Wave amplitude   Wave phase'
   write(25,*)'# No       rad/sec           m             rad'


   do i=1,no_spectral_comp
     write(25,101)i,spectral_freq(i),spectral_ampl(i),spectral_phase(i)
   enddo

   close(25)

101 format ( i5,x,150f14.7 )


 END Subroutine wave_init
!----------------------------------------------------------------------
 Subroutine random_wave ( delta_om, omega0 )
!----------------------------------------------------------------------

 Use Hydro

   implicit none

   real(8) :: delta_om, omega0
   real(4) :: rand1, rand2
   integer :: i, nsize_rand

!!!integer,allocatable :: k(:)
!!!integer             :: m

   call random_seed(SIZE=nsize_rand)
   call random_seed(put=IRand(1:nsize_rand) )   !setting starting seed value
!! call random_seed()
!! call random_seed(put=IRand(1:2) )   !setting starting seed value

!!! call random_seed(size=m)
!!! Allocate(k(m))
!!! call random_seed(get=k(1:m))       !this gets the current integers used as a starting value
!!! call random_seed(put=k(1:m) )      !setting starting seed value

   do i=1,no_spectral_comp
      call random_number(rand1)
      call random_number(rand2)

!     write(*,*) i, rand1, rand2

      spectral_phase(i) = dble(rand1) * pi2_h
      spectral_freq (i) = omega0  +  dble(i-1) * delta_om  +  dble(rand2) * delta_om 
   enddo


 END Subroutine random_wave
!----------------------------------------------------------------------
 Subroutine spectrum
!----------------------------------------------------------------------

 Use Hydro

   implicit none

   real(8) :: delta_om, freq, S
   integer :: i


   open (35,file='spectrum.out')

   write(35,*)'#Calculated and Randomized Wave Spectrum'
   write(35,*)'#         Frequency    Spectral dens.  Wave phase'
   write(35,*)'# No       rad/sec      m^2*sec/rad       rad'

   do i=1,no_spectral_comp

      freq = spectral_freq(i)

      call JONSWAP ( freq, S )

      spectral_ampl (i) = S

      write(35,101)i,spectral_freq(i)            , spectral_ampl(i)                   ,&
                                                   spectral_ampl(i)/spectral_freq(i)  ,&  !oc3 based      (u 1:3)
                     spectral_freq(i)*Tp_sp/pi2_h, spectral_ampl(i)/(Hs_sp**2*Tp_sp)      !Faltinsen book (u 4:5)

   enddo

   close(35)

      delta_om                        = ( spectral_freq(2)-0.d0 )/2.d0
      spectral_ampl(1)                = dsqrt( 2.d0*spectral_ampl(1) * delta_om ) 


      delta_om                        = ( Omega_Max_sp - spectral_freq(no_spectral_comp-1) )/2.d0
      spectral_ampl(no_spectral_comp) = dsqrt( 2.d0*spectral_ampl(no_spectral_comp) * delta_om )


   do i=2,no_spectral_comp-1
      delta_om                        = ( spectral_freq(i+1)-spectral_freq(i-1) )/2.d0
      spectral_ampl(i)                = dsqrt( 2.d0*spectral_ampl(i) * delta_om )
   enddo


101 format ( i5,x,150f14.7 )


 END Subroutine spectrum
!----------------------------------------------------------------------
 Subroutine int_read_spec ( inp_spectral_freq, inp_spectral_dens, no_inp_wave_freq )
!----------------------------------------------------------------------

 Use Hydro

   implicit none

   real(8) :: inp_spectral_freq(no_inp_wave_freq)
   real(8) :: inp_spectral_dens(no_inp_wave_freq)
   real(8) :: delta_om, omega, fy_sp 
   integer :: i, no_inp_wave_freq
   real(8), allocatable :: d2fy (:)


   open (35,file='spectrum.out')

   Allocate ( d2fy (no_inp_wave_freq) )

   call spline1 ( inp_spectral_freq, inp_spectral_dens, d2fy, no_inp_wave_freq )


   write(35,*)'#Interpolated and Randomized Wave Spectrum'
   write(35,*)'#         Frequency    Spectral dens.  Wave phase'
   write(35,*)'# No       rad/sec      m^2*sec/rad       rad'

   do i=1,no_spectral_comp

     omega = spectral_freq(i)
     call splint1  (inp_spectral_freq, inp_spectral_dens, d2fy, no_inp_wave_freq, omega, fy_sp )

     if ( fy_sp < 0.d0 ) then

        write(* ,*) 'inperpolated spectral density negative,',i,fy_sp
        write(10,*) 'inperpolated spectral density negative,',i,fy_sp
        spectral_ampl (i) = 0.d0

     else

        write(* ,*) 'interpolated spectral density ok',i,fy_sp
        write(10,*) 'interpolated spectral density ok',i,fy_sp
        spectral_ampl (i) = fy_sp

     endif

     write(35,101)i,spectral_freq(i)            , spectral_ampl(i)                   ,&
                                                  spectral_ampl(i)/spectral_freq(i)  ,&  !oc3 based      (u 1:3)
                    spectral_freq(i)*Tp_sp/pi2_h, spectral_ampl(i)/(Hs_sp**2*Tp_sp)      !Faltinsen book (u 4:5)
   enddo

   close(35)


     delta_om                        = ( spectral_freq(2)-inp_spectral_freq(1) )/2.d0
     spectral_ampl(1)                = dsqrt( 2.d0*spectral_ampl(1) * delta_om ) 


     delta_om                        = ( inp_spectral_freq(no_inp_wave_freq)-spectral_freq(no_spectral_comp-1) )/2.d0
     spectral_ampl(no_spectral_comp) = dsqrt( 2.d0*spectral_ampl(no_spectral_comp) * delta_om )


   do i=2,no_spectral_comp-1
     delta_om                        = ( spectral_freq(i+1)-spectral_freq(i-1) )/2.d0
     spectral_ampl(i)                = dsqrt( 2.d0*spectral_ampl(i) * delta_om )
   enddo


101 format ( i5,x,150f14.7 )

   Deallocate ( d2fy )


 END Subroutine int_read_spec
!----------------------------------------------------------------------
!
!     Subroutine :WAVE_NUM
!
!     Determine Wave Number of each wave component if morison's eq.
!     is used. Wavenum MUST be INITIALISED as an input in order to 
!     use the previous value as starting guess...
!
!-----------------------------------------------------------------------
 Subroutine WAVE_NUM ( Wavefreq, depth, gravity, Wavenum, kd_lim )
!-----------------------------------------------------------------------

   implicit none

   real(8) :: Wavefreq, depth, gravity, Wavenum
   real(8) :: FUNC, DFUNC, DWavenum, Wavefreq2, kd, kd_lim
   integer :: i


   Wavefreq2 = Wavefreq**2

   do i  = 1, 15 
      kd = Wavenum*depth

      if (kd < kd_lim) then                                   !-- 340.d0
         FUNC  = gravity*Wavenum *dtanh(kd) - Wavefreq2
         DFUNC = gravity         *dtanh(kd) +          &
                 gravity*Wavenum/(dcosh(kd))**2*depth
      else
         FUNC  = gravity*Wavenum            - Wavefreq2
         DFUNC = gravity
      endif

      DWavenum = -FUNC/DFUNC
      Wavenum  = Wavenum + DWavenum
   
!w    write(*,*) i,Wavenum,kd

      if (dabs(DWavenum) < 1.d-13) return
   enddo

   write(*,*)'dispertion relation not found',dabs(DWavenum)
   stop


 END Subroutine WAVE_NUM
!----------------------------------------------------------------------
 Subroutine JONSWAP ( freq, S )
!----------------------------------------------------------------------
!
! Implementation of Jonswap (Wavemod=2) and Pierson Moskowitz (Wavemod=3)
! Spectra from IEC.
!
! Here S(w) is given. S(w)=S(f)/2π
!
! 28.2.2014 jim: added White noise (Wavemod=7)
!----------------------------------------------------------------------

 Use Hydro

   implicit none

   real(8) :: sigma, gama, wTpPI2, S, freq
   real(8) :: freq1, freq2


   if (WaveMod == 7) then
!------ White noise ---------------
!------ freq1-freq2 define the frequency band of the noise.
      freq1 = 0.005d0 * pi2_h !!!pi2_h/40.d0 !0.d0 !0.05d0
      freq2 = Omega_Max_sp
      S     = 0.d0
      if (freq<freq1.or.freq>freq2) return
      S     = Hs_sp !Hs_sp**2/( 8.d0*(freq2-freq1)*dble(no_spectral_comp) ) !1.d0
      return
   endif

!--- Pierson Moskowitz type -------
   wTpPI2 = freq*Tp_sp/pi2_h
   S      = 0.3125d0/pi2_h * Hs_sp**2 * Tp_sp * wTpPI2**(-5d0) * dexp( -1.25d0 * wTpPI2**(-4d0) )

   if (WaveMod == 3) return


!--- Jonswap type -----------------
!--- define sigma coefficient
   if ( freq <= pi2_h/Tp_sp ) then
      sigma = 0.07d0
   else
      sigma = 0.09d0
   endif

!--- define gama coefficient
   if     ( Tp_sp/dsqrt(Hs_sp) <= 3.6d0 ) then
      gama  = 5d0
   elseif ( Tp_sp/dsqrt(Hs_sp) >  5.0d0 ) then
      gama  = 1d0
   else
      gama  = dexp ( 5.75d0 - 1.15d0 * Tp_sp/dsqrt(Hs_sp) )
   endif

!  gama  = 2.87d0
   write(10,*)'Jonswap prectrum gama coeff:',gama

   S = S * (  1 - 0.287d0*dlog(gama)  )  * gama**dexp (  -0.5d0* ( (wTpPI2-1)/sigma )**2  )


 END Subroutine JONSWAP

#include "floater.f90"
#include "substruct.f90"
#include "uinflow_sea.f90"
