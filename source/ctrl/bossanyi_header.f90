include '10MW_controller.f90'
include '10MW_controller_fcns.f90'

 Subroutine DISCON (avrSWAP, aviFAIL, accINFILE, avcOUTNAME, avcMSG)
 implicit none

!Compiler specific: Tell the complier that this routine is the entry point for the DLL
!The next line is for the case of the Digital Visual Fortran compiler
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:'DISCON' :: DISCON

!--- Arguments
   real          :: avrSWAP(*)
   integer       :: aviFAIL
   integer*1     :: accINFILE(*), avcOUTNAME(*), avcMSG(*)

!--- Other declarations
   integer       :: iStatus
!  integer, save :: iDone
   real*8 , save :: rPitchDemand, rTorqueDemand, rPI, rPitDem0
   real*8 , save :: array1(1000), array2(100), rGearRatio, rEff

   call CHAR2INT(' ',avcMsg)

   iStatus = NINT(avrSwap(1))
!--- Initialise
   if (iStatus .EQ. 0) THEN
!------ Initialise demands to current values
      rPitDem0 = avrSWAP(4)
      rPitchDemand = avrSWAP(4)
      rTorqueDemand = avrSWAP(23)
      rPI = 4*ATAN(1.0)
!------ Gear Box
      rGearRatio = 50                   ! Gearbox ratio
      rEff       = 0.94        !0.944   ! Electrical efficiency
!------ call DTU code
      array1( 1) = 10000./rEff !10000   ! Rated power [kW]
      array1( 2) = 0.628                ! Minimum rotor speed [rad/s] --> 6   rpm
      array1( 3) = 1.005                ! Rated   rotor speed [rad/s] --> 9.6 rpm
      array1( 4) = 15.6E6               ! Maximum allowable generator torque [Nm]
      array1( 5) = 0           !100     ! Minimum pitch angle, theta_min [deg],
                                        ! if |theta_min|>90, then a table of <wsp,theta_min> is read from a file named 'wptable.n', where n=int(theta_min)
      array1( 6) = 90                   ! Maximum pitch angle [deg]
      array1( 7) = 10                   ! Maximum pitch velocity operation [deg/s]
      array1( 8) = 0.2                  ! Frequency of generator speed filter [Hz]
      array1( 9) = 0.7                  ! Damping ratio of speed filter [-]
      array1(10) = 1.85        !0.64    ! Frequency of free-free DT torsion mode [Hz], if zero no notch filter used
!------ Partial load control parameters
      array1(11) = 1.0013E7    !9.5e6   ! Optimal Cp tracking K factor [Nm/(rad/s)^2],
      array1(12) = 6.83E7      !7.33e7  ! Proportional gain of torque controller [Nm/(rad/s)]
      array1(13) = 1.53E7      !1.32e7  ! Integral gain of torque controller [Nm/rad]
      array1(14) = 0.0                  ! Differential gain of torque controller [Nm/(rad/s^2)]
!------ Full load control parameters
      array1(15) = 1                    ! Generator control switch [1=constant power, 2=constant torque]
      array1(16) = 0.524       !0.592   ! Proportional gain of pitch controller [rad/(rad/s)]
      array1(17) = 0.141       !0.133   ! Integral gain of pitch controller [rad/rad]
      array1(18) = 0.0                  ! Differential gain of pitch controller [rad/(rad/s^2)]
      array1(19) = 0.4e-8               ! Proportional power error gain [rad/W]
      array1(20) = 0.4e-8               ! Integral power error gain [rad/(Ws)]
      array1(21) = 198         !164.13  ! Coefficient of linear term in aerodynamic gain scheduling, KK1 [deg]
      array1(22) = 693         !702.09  ! Coefficient of quadratic term in aerodynamic gain scheduling, KK2 [deg^2]
                                        ! (if zero, KK1 = pitch angle at double gain)
      array1(23) = 1.3                  ! Relative speed for double nonlinear gain [-]
!------ Cut-in simulation parameters
      array1(24) = 0.0         !0.1     ! Cut-in time [s], if zero no cut-in simulated
      array1(25) = 4.0                  ! Time delay for soft start [1/1P]
!------ Cut-out simulation parameters
      array1(26) = 0.0         !710     ! Cut-out time [s], if zero no cut-out simulated
      array1(27) = 5.0                  ! Time constant for 1st order filter lag of torque cut-out [s]
      array1(28) = 1                    ! Stop type [1=linear two pitch speed stop, 2=exponential pitch speed stop]
      array1(29) = 1.0                  ! Time delay for pitch stop 1 [s]
      array1(30) = 20.0                 ! Maximum pitch velocity during stop 1 [deg/s]
      array1(31) = 1.0                  ! Time delay for pitch stop 2 [s]
      array1(32) = 10.0                 ! Maximum pitch velocity during stop 2 [deg/s]
!------ Expert parameters (keep default values unless otherwise given)
      array1(33) = 0.5                  ! Lower angle above lowest minimum pitch angle for switch [deg]
      array1(34) = 0.5                  ! Upper angle above lowest minimum pitch angle for switch [deg]
      array1(35) = 95.0                 ! Ratio between filtered and reference speed for fully open torque limits [%]
      array1(36) = 5.0                  ! Time constant of 1st order filter on wind speed used for minimum pitch [1/1P]
      array1(37) = 5.0                  ! Time constant of 1st order filter on pitch angle for gain scheduling [1/1P]
!------ Drive train damper
      array1(38) = 0.0                  ! Proportional gain of DT damper [Nm/(rad/s)], requires frequency in input 10
!------ Over speed shutdown
      array1(39) = 10.0                 ! Percent maximum over speed of generator speed before cut-out [%]
!------ Additional non-linear pitch control term
      array1(40) = 0.0                  ! Err0 [rad/s]
      array1(41) = 0.0                  ! ErrDot0 [rad/s^2]
      array1(42) = 0.0                  ! PitNonLin1 [rad/s]

!------ The above default parameters may be modified by setting 'external controller parameters' if running from Bladed, or 
!------ otherwise by placing them in a datafile specified in the argument accINFILE. The format of the additional input is:
!------ Two numbers per line, the first being an integer i corresponding to the index of array1 (1 <= i<= 42), and
!------ the second being the value to be assigned to that element of the array, i.e. array1(i)
!------ Users can of course change or extend this to suit their needs
!jim  call GetModifiedParameters(array1, 42, accINFILE, NINT(avrSWAP(50)))

      call init_regulation(array1,array2,rTorqueDemand,rPitchDemand)
   endif

!---- Return demanded variables
!---- Note: previous values can be used before they updated, to simulate a one-sample delay
   if (NINT(avrSWAP(10)) .EQ. 0) THEN
        avrSWAP(42) = rPitchDemand
        avrSWAP(43) = rPitchDemand
        avrSWAP(44) = rPitchDemand
        avrSWAP(45) = rPitchDemand
        avrSWAP(47) = rTorqueDemand
   else
        call CHAR2INT('This DLL needs pitch position actuator',avcMsg)
        aviFail = -1
   endif

!--- Main calculation
   if (aviFail>= 0) then
!------ call DTU code
      array1(1) = avrSWAP( 2)                                      ! time
      array1(2) = avrSWAP(20) / rGearRatio                         ! Generator speed, ref LSS
      array1(3) = avrSWAP( 4)                                      ! Pitch 1
      array1(4) = avrSWAP(33)                                      ! Pitch 2
      array1(5) = avrSWAP(34)                                      ! Pitch 3
      array1(6) = avrSWAP(27) * COS(avrSWAP(24) + avrSWAP(37))     ! Wind u
      array1(7) =-avrSWAP(27) * SIN(avrSWAP(24) + avrSWAP(37))     ! Wind v
      array1(8) = 0.0                                              ! Wind w
      
      call update_regulation(array1,array2)
      rPitchDemand = array2(2)
      rTorqueDemand = array2(1) / rGearRatio
      
!------ Example logging output
      avrSWAP(65) = 3 !Number of logging variables
      call CHAR2INT('Pitch demand:A;Torquedemand:FL;DamperTorque:FL',avcOUTNAME)
      avrSWAP(NINT(avrSWAP(63))) = array2(2)
      avrSWAP(NINT(avrSWAP(63))+1) = array2(1) / rGearRatio
      avrSWAP(NINT(avrSWAP(63))+2) = array2(21) / rGearRatio
   endif

 END Subroutine DISCON
!==============================================================================================
!==============================================================================================
 Subroutine CHAR2INT(accString,avi1Int)
 implicit none

   character*(*) accString
   integer*1 avi1Int(*)
   integer I


   do I = 1,LEN_TRIM(accString)
      avi1Int(I) = ICHAR(accString(I:I))
   enddo

   avi1Int(LEN_TRIM(accString)+1) = 0


 END Subroutine CHAR2INT
!==============================================================================================
!==============================================================================================
 Subroutine GetModifiedParameters(array1, ArrayLength, INFILE, NameLength)
 implicit none

   integer   :: NameLength, ArrayLength
   real*8    :: array1(ArrayLength)
   integer*1 :: INFILE(NameLength)
   character(LEN=NameLength) :: FileName
   integer   :: I, IERR, IREC
   real*8    :: Val


   do I = 1, NameLength
      FileName(I:I) = CHAR(INFILE(I))
   enddo
  
   open (99, FILE=FileName, FORM='FORMATTED', IOSTAT=IERR)
    do while (IERR.EQ.0)
       read (99, *, IOSTAT=IERR) IREC, Val
       if (IERR.EQ.0.AND.IREC.GT.0.AND.IREC.LE.ArrayLength) array1(IREC) = Val
    enddo
   close(99)


 END Subroutine GetModifiedParameters
