!copy paste the following command before run the code, in order
!to be able to call the dll.
!set my_env_var=C:\Gast
!----------------------------------------------------------------------
 Subroutine DLL_Interface ( Gast2dll, dll2Gast, dllname_inp, dllcfg_inp )
!----------------------------------------------------------------------

 Use ifport
 Use ifwin

   implicit None


!-- variables for communication from/to GAST
   Real(kind=4)      :: Gast2dll(*), dll2Gast(*)

   character*256     :: dllname_inp, dllcfg_inp

!-- variables for dll loading
   character(256) value
   character(80) dll_name
   character(80) dll_sub_name
   integer(handle) lib_handle
   integer result

   interface
     subroutine sub ( avrSWAP, aviFAIL, accINFILE, avcOUTNAME, avcMSG )

      implicit none
      REAL(4),    INTENT(INOUT)    :: avrSWAP   (*)                                   ! The swap array, used to pass data to, and receive data from, the DLL controller.

      INTEGER(4), INTENT(  OUT)    :: aviFAIL                                         ! A flag used to indicate the success of this DLL call set as follows: 0 if the DLL call was successful, >0 if the DLL call was successful but cMessage should be issued as a warning messsage, <0 if the DLL call was unsuccessful or for any other reason the simulation is to be stopped at this point with cMessage as the error message.
      INTEGER(1), INTENT(IN   )    :: accINFILE (*)                                   ! The address of the first record of an array of 1-byte CHARACTERs giving the name of the parameter input file, 'DISCON.IN'.
      INTEGER(1), INTENT(  OUT)    :: avcMSG    (*)                                   ! The address of the first record of an array of 1-byte CHARACTERS giving the message contained in cMessage, which will be displayed by the calling program if aviFAIL <> 0.
      INTEGER(1), INTENT(IN   )    :: avcOUTNAME(*)                                   ! The address of the first record of an array of 1-byte CHARACTERS giving the simulation run name without extension.


     end subroutine sub
   end interface

   pointer (p,sub)

!-- Inputs to dll 
   Integer,Parameter :: npinf  = 256  ,&   !
                        npoutf = 6    ,&   !the exact size without the extension of the output file
                        npmsg  = 3000      !256

   Real(kind=4),save :: avrSWAP(0:1000) !qq 0:1000 added
   Integer(4)        :: aviFAIL

   CHARACTER(npinf ) :: cInFile   ! CHARACTER string giving the name of the parameter input file, 'DISCON.IN'
   CHARACTER(npmsg ) :: cMessage  ! CHARACTER string giving a message that will be displayed by the calling program if aviFAIL <> 0.
   CHARACTER(npoutf) :: cOutName  ! CHARACTER string giving the simulation run name without extension.

   INTEGER(1)        :: iInFile (npinf )
   INTEGER(1)        :: iMessage(npmsg )
   INTEGER(1)        :: iOutName(npoutf)


   EQUIVALENCE (iInFile , cInFile )
   EQUIVALENCE (iMessage, cMessage)
   EQUIVALENCE (iOutName, cOutName)


 ! open (15, file='dll.inp')
 !  read(15,*) cInFile
 !  read(15,*) dll_name
 ! close(15)

   cInFile      = dllcfg_inp
   dll_name     = dllname_inp
   cOutName     = "output"
   dll_sub_name = 'DISCON'


!--set communicatation variables
!-----------------------------------

   if (int(Gast2dll(1))==0) avrSWAP = 0.;

   avrSWAP ( 0) = 67                                    ! Number of data of this array
   avrSWAP (10) = Gast2dll( 1)           ! 1            ! istatus (0:init, 1:normal)
   avrSWAP (11) = Gast2dll( 2)           ! 2            ! time                     [sec]
   avrSWAP (12) = Gast2dll(10)           !27            ! Horizontal hub vel       [m/s]   
   avrSWAP (13) = 0.                                    ! Nacelle angle from North [rad]
   avrSWAP (14) = 0.                                    ! Yaw error angle          [rad]
   avrSWAP (15) = Gast2dll( 5)           !50            ! rotor azimuth angle      [rad]
   avrSWAP (16) = Gast2dll( 6)           ! 4            ! bld1 pitch               [rad]
   avrSWAP (17) = Gast2dll( 7)           !33            ! bld2 pitch               [rad]
   avrSWAP (18) = Gast2dll( 8)           !34            ! bld3 pitch               [rad]
   avrSWAP (19) = max(Gast2dll( 4),1e-5) !21            ! rotorspeed               [rad/s]
   avrSWAP (20) = max(Gast2dll( 3),1e-5) !20            ! genspeed                 [rad/s]
   avrSWAP (21) = Gast2dll(14)           !15            ! measured gen power       [W]
   avrSWAP (22) = Gast2dll(13)           !23            ! measured gen torque      [Nm]
   avrSWAP (23) = Gast2dll(15)           !35            ! gen contactor  [0:off, 1:HS or varible, 2:LS]
   avrSWAP (24) = Gast2dll(16)           !36            ! brake flag     [0:off, 1:on]
   avrSWAP (25) = Gast2dll(31)           !30        ! blade 1 root flap-wise   bending moment [Nm]
   avrSWAP (26) = Gast2dll(34)           !69        ! blade 1 root edge-wise   bending moment [Nm]
   avrSWAP (27) = Gast2dll(32)           !31        ! blade 2 root flap-wise   bending moment [Nm]
   avrSWAP (28) = Gast2dll(35)           !70        ! blade 2 root edge-wise   bending moment [Nm]
   avrSWAP (29) = Gast2dll(33)           !32        ! blade 3 root flap-wise   bending moment [Nm]
   avrSWAP (30) = Gast2dll(36)           !71        ! blade 3 root edge-wise   bending moment [Nm]

   avrSWAP (31) = 0.                     !          ! blade 1 root out-of-plae bending moment [Nm]
   avrSWAP (32) = 0.                     !          ! blade 1 root in-plane    bending moment [Nm]
   avrSWAP (33) = 0.                     !          ! blade 2 root out-of-plae bending moment [Nm]
   avrSWAP (34) = 0.                     !          ! blade 2 root in-plane    bending moment [Nm]
   avrSWAP (35) = 0.                     !          ! blade 3 root out-of-plae bending moment [Nm]
   avrSWAP (36) = 0.                     !          ! blade 3 root in-plane    bending moment [Nm]

   avrSWAP (37) = 0.                     !          ! Mx rotating   at hub center [Nm]
   avrSWAP (38) = 0.                     !          ! My rotating   at hub center [Nm]
   avrSWAP (39) = 0.                     !          ! Mz rotating   at hub center [Nm]
   avrSWAP (40) = 0.                     !          ! My stationary at hub center [Nm]
   avrSWAP (41) = 0.                     !          ! Mz stationary at hub center [Nm]
   avrSWAP (42) = 0.                     !          ! My            at tower top [Nm]
   avrSWAP (43) = 0.                     !          ! Mz            at tower top [Nm]
   avrSWAP (44) = 0.                     !          ! Mx            at tower bottom [Nm]
   avrSWAP (45) = 0.                     !          ! My            at tower bottom [Nm]

   avrSWAP (46) = Gast2dll(11)           !53            ! tower top fore-aft  acceleration [m/s^2]
   avrSWAP (47) = Gast2dll(12)           !54            ! tower top side-side acceleration [m/s^2]
   avrSWAP (48) = Gast2dll(30)           !16        ! optimal mode gain (considering the density and the performance) [Nm(rad/s)^2]

   avrSWAP (49) = 1.                     !          ! Generator         performance for the current working conditions [-]
   avrSWAP (50) = 1.                     !          ! Gearbox           performance for the current working conditions [-]
   avrSWAP (51) = 1.                     !          ! Total power train performance for the current working conditions [-]

   avrSWAP (52) = 0.                     !          ! Pitch actuator torque Blade 1 [Nm]
   avrSWAP (53) = 0.                     !          ! Pitch actuator torque Blade 2 [Nm]
   avrSWAP (54) = 0.                     !          ! Pitch actuator torque Blade 3 [Nm]

   avrSWAP (55) = Gast2dll(29)           ! 5        ! below rate pitch angle set point [??]
   avrSWAP (56) = max(Gast2dll(17),1e-5) !14            ! measured shaft power

   avrSWAP (57) = 0.                     !              ! Optimal Model Maximum Speed
   avrSWAP (58) = 0.                     !              ! Start Below Rate Torque Speed LookUp Table
   avrSWAP (59) = 0.                     !              ! Num Points Torque Lookup Table
   avrSWAP (60) = 0.                     !              ! Yaw Bearing MyGL [Nm]
   avrSWAP (61) = 0.                     !              ! Yaw Bearing MzGL [Nm]
   avrSWAP (62) = 0.                     !              ! Varialbe Slip Current Demand [A]

   avrSWAP (63) = 0.                     !              ! Naccelle Roll    Acceleration [m/s^2]
   avrSWAP (64) = 0.                     !              ! Naccelle Nodding Acceleration [m/s^2]
   avrSWAP (65) = 0.                     !              ! Naccelle Yaw     Acceleration [m/s^2]
   avrSWAP (66) = 0.                     !              ! Safety System Number that Has Been Activated
   avrSWAP (67) = 0.                     !              ! Extra Signals


!!!avrSWAP ( 3) = Gast2dll(22)                          ! Communication interval
!!!avrSWAP (61) = Gast2dll( 9)                          ! numbler of blds

!!!
!!!avrSWAP ( 6) = Gast2dll(18)                          ! min pitch angle
!!!avrSWAP ( 7) = Gast2dll(19)                          ! max pitch angle
!!!avrSWAP ( 8) = Gast2dll(20)                          ! min pitch rate
!!!avrSWAP ( 9) = Gast2dll(21)                          ! max pitch rate
!!!avrSWAP (10) = Gast2dll(27)                      ! pitch actuator           [0:pitch angle demand, 1:pitch rate demand]

!!!avrSWAP( 17) = Gast2dll(23)                          ! Minimum Generator Speed                        =   94.24770 rad/s
!!!avrSWAP( 18) = Gast2dll(24)                          ! Nominal Generator Speed                        =   175.9290 rad/s  (used for G10X)
!!!avrSWAP( 19) = Gast2dll(25)                          ! Demanded generator speed above rated [Maximum] =   198.9674 rad/s  (used for G10X)
!!!avrSWAP( 22) = Gast2dll(26)                          ! Demanded generator torque [Nominal]            = 11795.00 Nm       (used for G10X)

!!!avrSWAP( 28) = Gast2dll(28)                      ! Kind of pitch control    [0:collective, 1:individual]

!!!avrSWAP( 62) = 0.                                    ! maximum number of logging variables
!!!avrSWAP( 63) = 0.                                    ! record of first logging variable


!!!avrSWAP (49) = real(npmsg )                          ! maximum chars in MESSAGE for BLADED    !3000
!!!avrSWAP (50) = real(npinf )                          ! file INFILE length                     ! 256
!!!avrSWAP (51) = real(npoutf)                          ! file OUTNAME length                    !   6
!!!avrSWAP (64) = 2048.                                 ! maximum chars in outname

!!!avrSWAP(150) =    3.                                 ! initial state of the wind turbine
!            !! wind turbine status
!            !! STATE 0 - EMERGENCY STOP
!            !! STATE 1 - CRITICAL STOP
!            !! STATE 2 - START UP
!            !! STATE 3 - RUN
!            !! STATE 4 - PITCH FREEZE
!            !! STATE 5 - NORMAL STOP
!            !! STATE 6 - GRID LOSS



!-- Load the dll
!-----------------

   result = getenvqq("my_env_var", value)

   lib_handle   = LoadLibrary(trim(value)//'\'//trim(dll_name)//achar(0))
   p            = GetProcAddress(lib_handle, trim(dll_sub_name)//achar(0))    

          write(*,*) 'my_env_var = ', trim(value)
          write(*,*)'jim',result,value, trim(value)
          write(*,'(z16.16)') lib_handle
          write(*,'(z16.16)') p

  call sub ( avrSWAP, aviFAIL, iInFile, iOutName, iMessage )

  if (aviFAIL/=0) write(*,*)'AviFAIL',aviFAIL
  if (aviFAIL/=0) write(*,*)cMessage

 
    
!--get communicatation variables from the controller
!----------------------------------------------------

!                avrSWAP( 0)          ! Number of data of this array
   dll2Gast(5) = avrSWAP(10)   !1     ! status flag    [0:init, 1:run, -1:last]
   dll2Gast(6) = avrSWAP(11)   !35    ! gen contactor  [0:off, 1:HS or varible, 2:LS]
   dll2Gast(7) = avrSWAP(12)   !36    ! brake flag     [0:off, 1:on]
   dll2Gast(2) = avrSWAP(13)   !42    ! Use the command angles of all blades if using individual pitch
   dll2Gast(3) = avrSWAP(14)   !43    ! "
   dll2Gast(4) = avrSWAP(15)   !44    ! "
!                0              
   dll2Gast(1) = avrSWAP(17)   !47    ! Demanded generator torque
!                        18          ! Yaw         clock-wise actuation [-]
!                        19          ! Yaw counter clock-wise actuation [-]
!                0
!                        21          ! Demanded Nacelle Yaw Rate [rad/s]
!                0
!                0
!                0
!                        25          ! Request for loads=1 [-]
!                0
!                0
!                0
!                0
!                0
!                0

!!!dll2Gast(8) = avrSWAP(45)           ! demanded collective pitch angle
!!!dll2Gast(9) = avrSWAP(46)           ! demanded collective pitch rate


 END Subroutine DLL_Interface
