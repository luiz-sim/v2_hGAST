!**************************************************************************************************
subroutine normal_operation(GenSpeed, PitchVect, wsp, Pe, TTfa_acc, GenTorqueRef, PitchColRef, dump_array)
   !
   ! Controller for normal operation.
   !
   real(mk), intent(in)    :: PitchVect(3) ! Measured pitch angles [rad].
   real(mk), intent(in)    :: GenSpeed     ! Measured generator speed [rad/s].
   real(mk), intent(in)    :: wsp          ! Measured wind speed [m/s].
   real(mk), intent(in)    :: Pe           ! Measured electrical power [W].
   real(mk), intent(in)    :: TTfa_acc     ! Measured tower top longitudinal acceleration.
   real(mk), intent(out)   :: GenTorqueRef   ! Generator torque reference [Nm].
   real(mk), intent(out)   :: PitchColRef    ! Reference collective pitch [rad].
   real(mk), intent(inout) :: dump_array(50) ! Array for output.
   real(mk) WSPfilt
   real(mk) GenSpeedFilt, dGenSpeed_dtFilt
   real(mk) PitchMean, PitchMeanFilt, PitchMin
   real(mk) GenSpeedRef_full
   real(mk) Qdamp_ref, theta_dam_ref, P_filt
   real(mk) x, y(2)
   !***********************************************************************************************
   ! Inputs and their filtering
   !***********************************************************************************************
   ! Mean pitch angle
   PitchMean = (PitchVect(1) + PitchVect(2) + PitchVect(3)) / 3.0_mk
   ! Low-pass filtering of the rotor speed
   y = lowpass2orderfilt(deltat, stepno, omega2ordervar, GenSpeed)
   GenSpeedFilt = y(1)
   dGenSpeed_dtFilt = y(2)
   ! Low-pass filtering of the mean pitch angle for gain scheduling
   PitchMeanFilt = lowpass1orderfilt(deltat, stepno, pitchfirstordervar, PitchMean)
   PitchMeanFilt = min(PitchMeanFilt, 30.0_mk*degrad)
   ! Low-pass filtering of the nacelle wind speed
   WSPfilt = lowpass1orderfilt(deltat, stepno, wspfirstordervar, wsp)
   ! Minimum pitch angle may vary with filtered wind speed
   PitchMin = GetOptiPitch(WSPfilt)
   !***********************************************************************************************
   ! Limit reference speed for storm control
   !***********************************************************************************************
   if (Vcutout .gt. Vstorm) then
      GenSpeedRef_full = GenSpeedRefMax - max(0.0_mk, &
                         (WSPfilt - Vstorm)/(Vcutout - Vstorm)*(GenSpeedRefMax - GenSpeedRefMin))
   else
      GenSpeedRef_full = GenSpeedRefMax
   endif
   GenSpeedRef_full = max(min(GenSpeedRef_full, GenSpeedRefMax), GenSpeedRefMin)
   !***********************************************************************************************
   ! PID regulation of generator torque
   !***********************************************************************************************
   call torquecontroller(GenSpeed, GenSpeedFilt, dGenSpeed_dtFilt, PitchMean, WSPfilt, PitchMin, &
                         GenSpeedRef_full, Pe, GenTorqueRef, dump_array)
   !***********************************************************************************************
   ! Active DT damping based on filtered rotor speed
   !***********************************************************************************************
   call drivetraindamper(GenSpeed, Qdamp_ref, dump_array)
   TimerGenCutin = TimerGenCutin + deltat
!qq
   if (CutinVar%time.gt.0.0_mk) then
      x = switch_spline(TimerGenCutin, CutinVar%delay, 2.0_mk*CutinVar%delay)
!old  x = switch_spline(time         ,t_generator_cutin+t_cutin_delay,t_generator_cutin+2.d0*t_cutin_delay)
   else
      x = 1.0_mk
   endif
!qqGenTorqueRef = min(max(GenTorqueRef + Qdamp_ref*x, 0.0_mk), GenTorqueMax)
   GenTorqueRef = min(max(GenTorqueRef - Qdamp_ref*x, 0.0_mk), GenTorqueMax)
   !***********************************************************************************************
   ! PID regulation of collective pitch angle
   !***********************************************************************************************
   call pitchcontroller(GenSpeedFilt, dGenSpeed_dtFilt, PitchMeanFilt, Pe, PitchMin, &
                        GenSpeedRef_full, PitchColRef, dump_array)
   !***********************************************************************************************
   ! Active Tower damping based on filtered tower top aceleration
   !***********************************************************************************************
   P_filt = lowpass1orderfilt(deltat, stepno, TTfa_PWRfirstordervar, GenTorqueRef*GenSpeedFilt)
   call towerdamper(TTfa_acc, theta_dam_ref, dump_array)
   x = switch_spline(P_filt, TTfa_PWR_lower*PeRated, TTfa_PWR_upper*PeRated)
   PitchColRef = min(max(PitchColRef + theta_dam_ref*x, PID_pit_var%outmin), PID_pit_var%outmax)
   ! Write into dump array
   dump_array(1) = GenTorqueRef*GenSpeed
   dump_array(2) = WSPfilt
   dump_array(3) = GenSpeedFilt
   dump_array(20) = PitchMeanFilt
end subroutine
!**************************************************************************************************
subroutine torquecontroller(GenSpeed, GenSpeedFilt, dGenSpeed_dtFilt, PitchMean, WSPfilt, &
                            PitchMin, GenSpeedRef_full, Pe, GenTorqueRef, dump_array)
   !
   ! Generator torque controller. Controller that computes the generator torque reference.
   !
   real(mk), intent(in) :: GenSpeed          ! Measured generator speed [rad/s].
   real(mk), intent(in) :: GenSpeedFilt      ! Filtered generator speed [rad/s].
   real(mk), intent(in) :: dGenSpeed_dtFilt  ! Filtered generator acceleration [rad/s**2].
   real(mk), intent(in) :: PitchMean         ! Mean pitch angle [rad].
   real(mk), intent(in) :: PitchMin          ! Minimum pitch angle [rad].
   real(mk), intent(in) :: WSPfilt           ! Filtered wind speed [m/s].
   real(mk), intent(in) :: GenSpeedRef_full  ! Reference generator speed [rad/s].
   real(mk), intent(in) :: Pe                ! Measured electrical power [W].
   real(mk), intent(out) :: GenTorqueRef     ! Generator torque reference [Nm].
   real(mk), intent(inout) :: dump_array(50) ! Array for output.
   real(mk) GenTorqueMin_full, GenTorqueMax_full, GenTorqueMin_partial, GenTorqueMax_partial
   real(mk) GenSpeed_min1, GenSpeed_min2, GenSpeed_max1, GenSpeed_max2, GenSpeedRef
   real(mk) x, switch, switch_pitang_lower, switch_pitang_upper
   real(mk) kgain(3), GenSpeedFiltErr, outmin, outmax
!qq jim T=Kopt*w^pow
   real(mk) :: pow 

   pow = 2.0_mk !default
!  pow = 3.0_mk
   !***********************************************************************************************
   ! Speed ref. changes max. <-> min. for torque contr. and remains at rated for pitch contr.
   !***********************************************************************************************
   select case (PartialLoadControlMode)
   case (1)
      if (GenSpeedFilt .gt. 0.5_mk*(GenSpeedRefMax + GenSpeedRefMin)) then
         GenSpeedRef = GenSpeedRefMax
      else
         GenSpeedRef = GenSpeedRefMin
      endif
   case (2)
      GenSpeedRef = WSPfilt*TSR_opt/R
      GenSpeedRef = min(max(GenSpeedRef, GenSpeedRefMin), GenSpeedRefMax)
   end select
   ! Rotor speed error
   GenSpeedFiltErr = GenSpeedFilt - GenSpeedRef
   !-----------------------------------------------------------------------------------------------
   ! Limits for full load
   !-----------------------------------------------------------------------------------------------
   if (const_power) then 
!qq  GenTorqueMin_full = min((GenTorqueRated*GenSpeedRef_full)/max(GenSpeed    , GenSpeedRefMin), GenTorqueMax)
     GenTorqueMin_full = min((GenTorqueRated*GenSpeedRef_full)/max(GenSpeedFilt, GenSpeedRefMin), GenTorqueMax)
     GenTorqueMax_full = GenTorqueMin_full
   else 
     GenTorqueMin_full = GenTorqueRated
     GenTorqueMax_full = GenTorqueMin_full
   endif
!qq
   dump_array(40) = GenTorqueMin_full 
   dump_array(41) = GenTorqueMax_full 
   !-----------------------------------------------------------------------------------------------
   ! Limits for partial load that opens in both ends
   !-----------------------------------------------------------------------------------------------
   select case (PartialLoadControlMode)
     ! Torque limits for K Omega^2 control of torque
     case (1)
       ! Calculate the constant limits for opening and closing of torque limits
       GenSpeed_min1 = GenSpeedRefMin
       GenSpeed_min2 = GenSpeedRefMin/SwitchVar%rel_sp_open_Qg
       GenSpeed_max1 = (2.0_mk*SwitchVar%rel_sp_open_Qg - 1.0_mk)*GenSpeedRefMax
       GenSpeed_max2 = SwitchVar%rel_sp_open_Qg*GenSpeedRefMax
       ! Compute lower torque limits
       x = switch_spline(GenSpeedFilt, GenSpeed_min1, GenSpeed_min2)
       GenTorqueMin_partial = (Kopt*GenSpeedFilt**pow - Kopt_dot*dGenSpeed_dtFilt)*x
       GenTorqueMin_partial = min(GenTorqueMin_partial, Kopt*GenSpeed_max1**pow)
       x = switch_spline(GenSpeedFilt, GenSpeed_max1, GenSpeed_max2)
       ! Compute upper torque limits
       GenTorqueMax_partial = (Kopt*GenSpeedFilt**pow - Kopt_dot*dGenSpeed_dtFilt)*(1.0_mk - x) + GenTorqueMax_full*x
       GenTorqueMax_partial = max(GenTorqueMax_partial, Kopt*GenSpeed_min2**pow)
     ! Torque limits for PID control of torque
     case (2)
       GenTorqueMin_partial = 0.0_mk
       GenTorqueMax_partial = GenTorqueMax_full
   end select
   ! Interpolation between partial and full load torque limits based on pitch
!qqSwitchVar%pitang_lower   = array1(33)*degrad
!  SwitchVar%pitang_upper   = array1(34)*degrad
!  SwitchVar%rel_sp_open_Qg = array1(35)*0.01_mk
   ! Switch based on pitch
   switch_pitang_lower = SwitchVar%pitang_lower + PitchMin
   switch_pitang_upper = SwitchVar%pitang_upper + PitchMin
   switch = switch_spline(PitchMean, switch_pitang_lower, switch_pitang_upper)
   switch = lowpass1orderfilt(deltat, stepno, switchfirstordervar, switch)
   ! Interpolation between partial and full load torque limits based on switch 1
   outmin = (1.0_mk - switch)*GenTorqueMin_partial + switch*GenTorqueMin_full
   outmax = (1.0_mk - switch)*GenTorqueMax_partial + switch*GenTorqueMax_full
   !***********************************************************************************************
   ! Rotor speed exclusion zone
   !***********************************************************************************************
   call rotorspeedexcl(GenSpeedFilt, Pe/GenSpeed, GenTorqueMin_partial, GenTorqueMax_partial, GenSpeedFiltErr, &
                       outmax, outmin, dump_array)
   PID_gen_var%outmin = outmin
   PID_gen_var%outmax = outmax
   if (PID_gen_var%outmin .gt. PID_gen_var%outmax) PID_gen_var%outmin = PID_gen_var%outmax
   !-----------------------------------------------------------------------------------------------
   ! Compute PID feedback to generator torque demand
   !-----------------------------------------------------------------------------------------------
   kgain = 1.0_mk
   GenTorqueRef = PID(stepno, deltat, kgain, PID_gen_var, GenSpeedFiltErr)
   ! Write into dump array
   dump_array(4) = GenSpeedFiltErr
   dump_array(6) = PID_gen_var%outpro
   dump_array(7) = PID_gen_var%outset
   dump_array(8) = PID_gen_var%outmin
   dump_array(9) = PID_gen_var%outmax
   dump_array(10) = switch
end subroutine torquecontroller
!**************************************************************************************************
subroutine pitchcontroller(GenSpeedFilt, dGenSpeed_dtFilt, PitchMeanFilt, Pe, PitchMin, &
                           GenSpeedRef_full, PitchColRef, dump_array)
   !
   ! Pitch controller. Controller that computes the reference collective pitch angle.
   !
   real(mk), intent(in) :: GenSpeedFilt     ! Filtered generator speed [rad/s].
   real(mk), intent(in) :: dGenSpeed_dtFilt ! Filtered generator acceleration [rad/s**2].
   real(mk), intent(in) :: PitchMeanFilt    ! Filtered mean pitch angle [rad].
   real(mk), intent(in) :: PitchMin         ! Minimum pitch angle [rad].
   real(mk), intent(in) :: GenSpeedRef_full ! Reference generator speed [rad/s].
   real(mk), intent(in) :: Pe               ! Measured electrical power [W].
   real(mk), intent(out) :: PitchColRef     ! Reference collective pitch angle [rad].
   real(mk), intent(inout) :: dump_array(50) ! Array for output.
   real(mk) GenSpeedFiltErr, added_term, aero_gain, aero_damp, kgain(3, 2), err_pitch(2)
   ! Rotor speed error
   GenSpeedFiltErr = GenSpeedFilt - GenSpeedRef_full
   ! Additional nonlinear pitch control term
!qq  Err0       = array1(40)
!    ErrDot0    = array1(41)
!    PitNonLin1 = array1(42)
   if ((PitNonLin1 .gt. 0.0_mk).and.(Err0 .gt. 0.0_mk).and.(ErrDot0.gt.0.0_mk)) then
     added_term = GenSpeedFiltErr/Err0 + dGenSpeed_dtFilt/ErrDot0
     if (added_term .gt. 1.0_mk) then
       AddedPitchRate = PitNonLin1*added_term + AddedPitchRate
     endif
   endif
   ! Limits
   PID_pit_var%outmin = PitchMin
   PID_pit_var%outmax = PitchStopAng
!qq PitchGSVar%invkk1 = 1.0_mk/(array1(21)*degrad)
!   PitchGSVar%invkk2 = 1.0_mk/(array1(22)*degrad*degrad)
   ! Aerodynamic gain scheduling dQ/dtheta
   aero_gain = 1.0_mk + PitchGSVar%invkk1*PitchMeanFilt + PitchGSVar%invkk2*PitchMeanFilt**2
   kgain = 1.0_mk/aero_gain
   ! Nonlinear gain to avoid large rotor speed excursion
!qq rel_limit = array1(23)
   if (rel_limit .ne. 0.0_mk) then
     kgain = kgain*(GenSpeedFiltErr**2 / (GenSpeedRef_full*(rel_limit - 1.0_mk))**2 + 1.0_mk)
   endif
!qq PitchGSVar%kp_speed     = array1(50)
!   PitchGSVar%invkk1_speed = 1.0_mk/(array1(51)*degrad)
!   PitchGSVar%invkk2_speed = 1.0_mk/(array1(52)*degrad*degrad)
   ! Gainscheduling according to dQaero/dOmega
   aero_damp = 1.0_mk + PitchGSVar%invkk1_speed*PitchMeanFilt + &
               PitchGSVar%invkk2_speed*PitchMeanFilt**2
   PID_pit_var%kpro(1) = PID_pit_var%kpro(1) + PitchGSVar%kp_speed*aero_damp

write(3332,'(10e15.5)') Kgain(1,1),Kgain(2,1),PID_pit_var%kpro(1),&
                        (GenSpeedFiltErr**2 / (GenSpeedRef_full*(rel_limit - 1.0_mk))**2 + 1.0_mk)


   !-----------------------------------------------------------------------------------------------
   ! Compute PID feedback to pitch demand
   !-----------------------------------------------------------------------------------------------
   if (DT_mode_filt%f0 .gt. 0.0_mk) then
     err_pitch(1) = notch2orderfilt(deltat, stepno, DT_mode_filt, GenSpeedFiltErr)
     err_pitch(2) = notch2orderfilt(deltat, stepno, pwr_DT_mode_filt, Pe - PeRated)
   else
     err_pitch(1) = GenSpeedFiltErr
     err_pitch(2) = Pe - PeRated
   endif
   PitchColRef = PID2(stepno, deltat, kgain, PID_pit_var, err_pitch, AddedPitchRate)
   ! Write into dump array
   dump_array(11) = GenSpeedFiltErr
   dump_array(12) = err_pitch(2)
   dump_array(13) = PID_pit_var%outpro
   dump_array(14) = PID_pit_var%outset
   dump_array(15) = PID_pit_var%outmin
   dump_array(16) = PID_pit_var%outmax
   dump_array(19) = AddedPitchRate
end subroutine pitchcontroller
