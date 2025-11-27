!copy paste the following command before run the code, in order
!to be able to call the dll.
!set my_env_var=C:\Gast
!----------------------------------------------------------------------
 Subroutine DLL_Interface ( Gast2dll, dll2Gast, dllname_inp, dllcfg_inp )
!----------------------------------------------------------------------

#if   OS == 0
 use iso_c_binding
#elif OS == 1
 use ifport
 use ifwin
#endif

   implicit None

   character*256, intent(in   ) :: dllname_inp, dllcfg_inp
   real(kind=4) , intent(in   ) :: Gast2dll(*)              ! communication var vector from hGAST to controller
   real(kind=4) , intent(  out) :: dll2Gast(*)              ! communication var vector to controller from hGAST

!----- interface to call the controller DISCON
   interface
     subroutine sub ( avrSWAP, aviFAIL, accINFILE, avcOUTNAME, avcMSG )

      implicit none

      real   (4), intent(inout) :: avrSWAP   (*)            ! The swap array, used to pass data to, and receive data from, the DLL controller.
      integer(4), intent(  out) :: aviFAIL                  ! A flag used to indicate the success of this DLL call set as follows: 0 if the DLL call was successful, >0 if the DLL call was successful but cMessage should be issued as a warning messsage, <0 if the DLL call was unsuccessful or for any other reason the simulation is to be stopped at this point with cMessage as the error message.
      integer(1), intent(in   ) :: accINFILE (*)            ! The address of the first record of an array of 1-byte characters giving the name of the parameter input file, 'DISCON.IN'.
      integer(1), intent(  out) :: avcMSG    (*)            ! The address of the first record of an array of 1-byte characterS giving the message contained in cMessage, which will be displayed by the calling program if aviFAIL <> 0.
      integer(1), intent(in   ) :: avcOUTNAME(*)            ! The address of the first record of an array of 1-byte characterS giving the simulation run name without extension.

     end subroutine sub
   END interface


#if   OS == 0
!---- interface to linux API
   interface
       function dlopen(filename,mode) bind(c,name="dlopen")
           ! void *dlopen(const char *filename, int mode);
           use iso_c_binding
           implicit none
           type     (c_ptr )             :: dlopen
           character(c_char), intent(in) :: filename(*)
           integer  (c_int ), value      :: mode
       end function

       function dlsym(handle,name) bind(c,name="dlsym")
           ! void *dlsym(void *handle, const char *name);
           use iso_c_binding
           implicit none
           type(c_funptr)                :: dlsym
           type(c_ptr   )   , value      :: handle
           character(c_char), intent(in) :: name(*)
       end function

       function dlclose(handle) bind(c,name="dlclose")
           ! int dlclose(void *handle);
           use iso_c_binding
           implicit none
           integer(c_int)        :: dlclose
           type   (c_ptr), value :: handle
       end function
   END interface

!--- variables for dll loading in Linux
   integer(c_int), parameter :: rtld_lazy=1 ! value extracte from the C header file
!  integer(c_int), parameter :: rtld_now =2 ! value extracte from the C header file
   type(c_ptr   )            :: handle
   type(c_funptr)            :: address
   procedure(sub), pointer   :: subLinux

   character(256)            :: libpath
   character(256)            :: procname

#elif OS == 1

!--- variables for dll loading in Windows
   character(256) :: value
   integer(handle):: lib_handle
   integer        :: result
   pointer        :: (p,sub)
#endif
   character( 80) :: dll_name
   character( 80) :: dll_sub_name

!--- Inputs to dll 
   integer, parameter :: npinf  = 256       !
   integer, parameter :: npoutf = 6         ! the exact size without the extension of the output file
   integer, parameter :: npmsg  = 3000      ! 256
   real(kind=4), save :: avrSWAP(1000)
   integer(4)         :: aviFAIL
   character(npinf )  :: cInFile            ! character string giving the name of the parameter input file, 'DISCON.IN'
   character(npmsg )  :: cMessage           ! character string giving a message that will be displayed by the calling program if aviFAIL <> 0.
   character(npoutf)  :: cOutName           ! character string giving the simulation run name without extension.
   integer(1)         :: iInFile (npinf )
   integer(1)         :: iMessage(npmsg )
   integer(1)         :: iOutName(npoutf)

   EQUIVALENCE (iInFile , cInFile )
   EQUIVALENCE (iMessage, cMessage)
   EQUIVALENCE (iOutName, cOutName)


!--- set communicatation variables
   if (int(Gast2dll(1))==0) avrSWAP = 0.;

   avrSWAP ( 1) = Gast2dll( 1)                     ! istatus (0:init, 1:normal)                                  [-    ]         0 then 1
   avrSWAP ( 2) = Gast2dll( 2)                     ! time                                                        [sec  ]         runtime
   avrSWAP ( 3) = Gast2dll(22)                     ! Communication interval                                      [sec  ]         runtime in current implementation is always equat to DT_el
   avrSWAP (20) = max(Gast2dll( 3),1e-5)           ! genspeed                                                    [rad/s]         runtime
   avrSWAP (21) = max(Gast2dll( 4),1e-5)           ! rotorspeed                                                  [rad/s]         runtime
   avrSWAP (60) = Gast2dll( 5)                     ! rotor azimuth angle                                         [rad  ]         runtime
   avrSWAP ( 4) = Gast2dll( 6)                     ! bld1 pitch                                                  [rad  ]         runtime FJS: blade one pointing upwards at 0deg azimuth
   avrSWAP (33) = Gast2dll( 7)                     ! bld2 pitch                                                  [rad  ]         runtime FJS: blade two is at 120deg azimuth (right side of rotor plane when looking with the wind)
   avrSWAP (34) = Gast2dll( 8)                     ! bld3 pitch                                                  [rad  ]         runtime FJS: blade three is at 240deg azimuth (left side of rotor plane when looking with the wind)
   avrSWAP (61) = Gast2dll( 9)                     ! numbler of blds                                             [-    ]         runtime
   avrSWAP (27) = Gast2dll(10)                     ! Horizontal hub vel                                          [m/s  ]         runtime
   avrSWAP (53) = Gast2dll(11)                     ! tower top fore-aft  acceleration                            [m/s^2]         runtime FJS: positive fore to aft when looking with the wind (not rotating with the nacelle)
   avrSWAP (54) = Gast2dll(12)                     ! tower top side-side acceleration                            [m/s^2]         runtime FJS: positive right to left when looking with the wind
   avrSWAP (30) = Gast2dll(31)                     ! blade 1 root out of plane bending moment                    [Nm   ]         runtime FJS: positive in downwind direction
   avrSWAP (31) = Gast2dll(32)                     ! blade 2 root out of plane bending moment                    [Nm   ]         runtime
   avrSWAP (32) = Gast2dll(33)                     ! blade 3 root out of plane bending moment                    [Nm   ]         runtime
   avrSWAP (69) = Gast2dll(34)                     ! blade 1 root in plane bending moment                        [Nm   ]         runtime FJS: positive in direction of rotation
   avrSWAP (70) = Gast2dll(35)                     ! blade 2 root in plane bending moment                        [Nm   ]         runtime
   avrSWAP (71) = Gast2dll(36)                     ! blade 3 root in plane bending moment                        [Nm   ]         runtime
   avrSWAP (23) = Gast2dll(13)                     ! measured gen torque                                         [Nm   ]         runtime
   avrSWAP (15) = Gast2dll(14)                     ! measured gen power                                          [W    ]         runtime
   avrSWAP (35) = Gast2dll(15)                     ! gen contactor (0:off, 1:HS or varible, 2:LS)                [-    ]
   avrSWAP (36) = Gast2dll(16)                     ! brake flag    (0:off, 1:on]                                 [-    ]
   avrSWAP (14) = max(Gast2dll(17),1e-5)           ! measured shaft power                                        [W    ]         runtime
   avrSWAP ( 6) = Gast2dll(18)                     ! min pitch angle                                             [rad  ]
   avrSWAP ( 7) = Gast2dll(19)                     ! max pitch angle                                             [rad  ]
   avrSWAP ( 8) = Gast2dll(20)                     ! min pitch rate                                              [rad/s]
   avrSWAP ( 9) = Gast2dll(21)                     ! max pitch rate                                              [rad/s]
   avrSWAP (10) = Gast2dll(27)                     ! pitch actuator (0:pitch angle demand, 1:pitch rate demand)  [-    ]
   avrSWAP( 17) = Gast2dll(23)                     ! Minimum Generator Speed                                     [rad/s]         ! =   94.24770 rad/s
   avrSWAP( 18) = Gast2dll(24)                     ! Nominal Generator Speed                                     [rad/s]         ! =   175.9290 rad/s  (used for G10X)
   avrSWAP( 19) = Gast2dll(25)                     ! Demanded generator speed above rated [Maximum]              [rad/s]         ! =   198.9674 rad/s  (used for G10X)
   avrSWAP( 22) = Gast2dll(26)                     ! Demanded generator torque            [Nominal]              [Nm   ]         ! = 11795.00 Nm       (used for G10X)
   avrSWAP( 28) = Gast2dll(28)                     ! Kind of pitch control (0:collective, 1:individual)          [-    ]
   avrSWAP ( 5) = Gast2dll(29)                     ! below rate pitch angle set point                            [rad  ]
   avrSWAP (16) = Gast2dll(30)                     ! optimal mode gain                                           [Nms^2/rad^2]
                                                                                                                                                                                                                                     
   avrSWAP( 62) = 0.                               ! maximum number of logging variables
   avrSWAP( 63) = 0.                               ! record of first logging variable
   avrSWAP (49) = real(npmsg )                     ! maximum chars in MESSAGE for BLADED    !3000
   avrSWAP (50) = real(npinf )                     ! file INFILE length                     ! 256
   avrSWAP (51) = real(npoutf)                     ! file OUTNAME length                    !   6
   avrSWAP (64) = 2048.                            ! maximum chars in outname
   avrSWAP(150) =    3.                            ! initial state of the wind turbine
   avrSWAP(117) =    0.                            ! state of the wind turbine for ecn controller for tri-spar 10MW OFWT
                                                   ! (zero for normal production)
!            !! wind turbine status
!            !! STATE 0 - EMERGENCY STOP
!            !! STATE 1 - CRITICAL STOP
!            !! STATE 2 - START UP
!            !! STATE 3 - RUN
!            !! STATE 4 - PITCH FREEZE
!            !! STATE 5 - NORMAL STOP
!            !! STATE 6 - GRID LOSS


!--- Load the library
!--- On Linux you call dlopen and dlsym to get a pointer to the routine, 
!--- on Windows LoadLibrary and GetProcAddress.

 ! open (15, file='dll.inp')
 !  read(15,*) cInFile
 !  read(15,*) dll_name
 ! close(15)

   dll_name     = dllname_inp
   dll_sub_name = "DISCON"
   cInFile      = dllcfg_inp
   cOutName     = "output"

#if   OS == 0
!--- load the library
   libpath      = trim(dll_name)
   handle       = dlopen(trim(libpath)//c_null_char, rtld_lazy)
   if (.not. c_associated(handle)) then
       print*, 'Unable to load DLL'; stop
   endif

!--- load the procedure
   procname     = trim(dll_sub_name)
   address      = dlsym( handle, trim(procname)//c_null_char )
   if ( .not. c_associated(address) ) then
       write(*,*)'Unable to load the procedure', trim(procname); stop
   endif
   call c_f_procpointer(address, subLinux)

   call subLinux ( avrSWAP, aviFAIL, iInFile, iOutName, iMessage )

#elif OS == 1

   result       = getenvqq("my_env_var", value)                                !set my_env_var=C:\Gast
   lib_handle   = LoadLibrary(trim(value)//'\'//trim(dll_name)//achar(0))
   p            = GetProcAddress(lib_handle, trim(dll_sub_name)//achar(0))    

        ! write(*,*) 'my_env_var = ', trim(value)
        ! write(*,*)'jim',result,value, trim(value)
        ! write(*,'(z16.16)') lib_handle
        ! write(*,'(z16.16)') p

   call sub ( avrSWAP, aviFAIL, iInFile, iOutName, iMessage )
#endif

!  if (aviFAIL/=0) write(*,*)'AviFAIL',aviFAIL
!  if (aviFAIL/=0) write(*,*)cMessage

    
!--- Set communicatation variables from the controller
   dll2Gast(1) = avrSWAP(47) ! Demanded generator torque                                                                  [Nm   ]
   dll2Gast(2) = avrSWAP(42) ! pitch1 demand individual (Use the command angles of all blades if using individual pitch)  [rad or rad/s]
   dll2Gast(3) = avrSWAP(43) ! pitch2 demand individual (Use the command angles of all blades if using individual pitch)  [rad or rad/s]
   dll2Gast(4) = avrSWAP(44) ! pitch3 demand individual (Use the command angles of all blades if using individual pitch)  [rad or rad/s]
   dll2Gast(5) = avrSWAP( 1) ! status flag    [0:init, 1:run, -1:last]                                                    [-    ]
   dll2Gast(6) = avrSWAP(35) ! gen contactor  [0:off, 1:HS or varible, 2:LS]                                              [-    ]
   dll2Gast(7) = avrSWAP(36) ! brake flag     [0:off, 1:on]                                                               [-    ]
   dll2Gast(8) = avrSWAP(45) ! demanded collective pitch angle                                                            [rad  ]
   dll2Gast(9) = avrSWAP(46) ! demanded collective pitch rate                                                             [rad/s]


 END Subroutine DLL_Interface


!   INITIALISATION


!  avrSWAP( 52) =                    ! DLL interface version
!  avrSWAP (18) =                                        ! (omGenMod/omGenRef)*omRef;     (used for G10X)
!  avrSWAP (19) =                                        ! omRefSet;                      (used for G10X)

!  avrSWAP (37) =                                       !nacelle angle from north
!  avrSWAP (24) =                                       !yaw error angle
!  avrSWAP( 29) =                                       ! Kind Yaw Control

!  avrSWAP (25) = 0.                                     !start below rate torque speed lookup table
!  avrSWAP (26) = 0.                                     !num points torque lookup table
!  avrSWAP (77) = 0.                                     !yaw bearing MyGl
!  avrSWAP (78) = 0.                                     !yaw bearing MzGl
!  avrSWAP (79) = 1.                                     !Request for load demanded (Pitch atuator torque model)
!  avrSWAP (81) = 0.                                     !variable slip current demand
!  avrSWAP (82) = 0.                                     !nacelle roll acceleration
!  avrSWAP (83) = 0.                                     !nacelle nodding acceleration
!  avrSWAP (84) = 0.                                     !nacelle yaw acceleration
!  avrSWAP (97) = 0.                                     !safety system number that has been activated
!  avrSWAP (99) = 200.                                   !index (BLADED 1-Based) for load demanded (Pitch actuator torque model)
!  avrSWAP (200)=                                        !Pitch actuator torque blade 1
!  avrSWAP (201)=                                        !Pitch actuator torque blade 2
!  avrSWAP (202)=                                        !Pitch actuator torque blade 3


! ---- data from control
!   !!avrSWAP(19)      ! for (G10X) omRef

! ---- demanded yaw actuator torque
! ---- NOT USED
!   !!avrSWAP(48)      ! demanded nacelle yaw rate
!   !!avrSWAP(55)      ! pitch override
!   !!avrSWAP(56)      ! torque override
!   !!avrSWAP(72)      ! generator startup resistance
!   !!avrSWAP(79)      ! request for loads
!   !!avrSWAP(80)      ! variable slip current demand activate
!   !!avrSWAP(81)      ! variable slip current demand
!   !!avrSWAP(92)      ! mean wind speed increment
!   !!avrSWAP(93)      ! turbulence intensity increment
!   !!avrSWAP(94)      ! wind direction increment
!   !!avrSWAP(98)      ! safety system number to activate
! ---- NOT USED
!!    avrSWAP(xx)      ! copying Emergency Chain flag //MODIF_MFO
