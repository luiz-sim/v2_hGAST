!--------------------------------------!
*** hGAST short description/history
!--------------------------------------!
...
...
...


!--------------------------------------!
*** hGAST built on Linux
!--------------------------------------!
built including GenUVP & vpmfree:
   make -f Makefile_gvp cleanall
   make -f Makefile_gvp all
built only with BEM support:
   make cleanall
   make all


!--------------------------------------!
*** hGAST input  files
!--------------------------------------!
   paths.inp     : Input filenames/paths
   dfile_el.inp  : Main parameters
   machi.inp     : Structural beam properties for main bodies (blade, shaft, tower)
   jacket.inp    : Jacket's parameters (nodes, bodies, properties, b.c., buoyancy)
   foundation.inp: Foundation definition (springs)
   modal.inp     : Parameters for ROM version
   rotorin.inp   : Main aerodynamic parameters
 ./airfoils/     : Directory with look-up tables for CL, CD, CM (i.e. cylind*, du_*, naca* from NREL5MW)
   geomb.inp     : Geometry of the blade
   profilb.inp   : Airfoil profiles of the blade
   drag.inp      : Tower and nacelle drag parameters
   wind.inp      : Defined wind scenario input file
   hydro.inp     : Main hydrodynamic parameters
   floater.inp   : Floater’s parameters
   diffract.inp  : Wave excitation loads (solution of hydrodynamic problem in the frequency domain)            !omega [rad/s], Fx,y,z, Mx,y,z [N or Nm/(rho g A)], PhaseFx,y,z, Mx,y,z [rad], |Drift[Fx,y,z,Mx,y,z| [N or Nm /(rho g A^2)]
   retard.inp    : Retardation functions (solution of hydrodynamic problem in the frequency domain+Fourier)
   morison.inp   : Morison’s equation parameters (viscous hydrodynamic drag)
   truss.inp     : Dynamic mooring lines module main parameters
   eigen.inp     : Parameters for modal (eigenvalue) analysis


!--------------------------------------!
*** hGAST output files ***
!--------------------------------------!

!--- output filenames and corresponding module
** Structural/Dynamic module
      loadsXXX_yy.dat     ,XXX:body, yy :node
      deformXXX_yy.dat    ,XXX:body, yy :node
      qsdof.dat
      qsdof3.dat
** BEM module
      LOADS_aer.dat
      striplXX_ss.dat     ,XX :body, ss :strip
** Environmental excitation
      velhub.dat
      z_elevation.dat
** Modal analysis
      eigen.dat
      modeRZZZ.dat        ,ZZZ:mode,
      modeZZZ_XXX.dat     ,ZZZ:mode, XXX:body


!--- output file description per column
** loadsXXX_yy: Internal loads of body XXX at node yy,
                wrt the local c.s. of each body
                blades       [bodies 1-3 ] [node 1:root - nnt:tip]
                shaft        [body   4   ] [node 1:gen  - nnt:hub]
                tower        [body   5   ] [node 1:base - nnt:top]
                jacket       [bodies 6-..] [orientation is defined in jacket.inp]
       1: Time     [s  ]
       2: Azimuth  [deg]
       3: Fx       [kN ]     [blades: edge-wise, tower: fore-aft ]
       4: Fy       [kN ]     [blades: radial   , tower: vertical ]
       5: Fz       [kN ]     [blades: flap-wise, tower: side-side]
       6: Mx       [kNm]     [blades: flap-wise, tower: side-side]
       7: My       [kNm]     [blades: pitching , tower: yawing   ]
       8: Mz       [kNm]     [blades: edge-wise, tower: fore-aft ]
       9: Mc       [kNm]     [Mc=(Mx^2+Mz^2)^0.5, combinded bending moments]


** deformXXX_yy: Elastic deformations of body XXX at node yy,
                 wrt the local c.s. of each body
                 blades      [bodies 1-3 ] [node 1:root - nnt:tip]
                 shaft       [body   4   ] [node 1:gen  - nnt:hub]
                 tower       [body   5   ] [node 1:base - nnt:top]
                 jacket      [bodies 6-..] [orientation is defined in jacket.inp]
       1: Time     [s  ]
       2: Azimuth  [deg]
       3: Ux       [m  ]     [blades: edge-wise, tower: fore-aft ] ** wrt body local c.s.
       4: Uy       [m  ]     [blades: radial   , tower: vertical ]
       5: Uz       [m  ]     [blades: flap-wise, tower: side-side]
       6: Thx      [deg]     [blades: flap-wise, tower: side-side]
       7: Thy      [deg]     [blades: pitching , tower: yawing   ]
       8: Thz      [deg]     [blades: edge-wise, tower: fore-aft ]
       9: Ux       [m  ]     [blades: edge-wise, tower: fore-aft ] ** wrt station local c.s.
      10: Uy       [m  ]     [blades: radial   , tower: vertical ]
      11: Uz       [m  ]     [blades: flap-wise, tower: side-side]
      12: Thx      [deg]     [blades: flap-wise, tower: side-side]
      13: Thy      [deg]     [blades: pitching , tower: yawing   ]
      14: Thz      [deg]     [blades: edge-wise, tower: fore-aft ]


** qsdof.dat : dofs of additional equations
       1: Time                         [s  ]
       2: Azimuth                      [rad]
       3: Omega                        [rad/s]
       4: Domega/dt                    [rad/s^2]
       5: Blade1 pitch angle           [deg]
       6: Blade1 pitch velocity        [deg/s]
       7: Blade1 pitch acceleration    [deg/s^2]
       8: Blade2 pitch angle           [deg]
       9: Blade2 pitch velocity        [deg/s]
      10: Blade2 pitch acceleration    [deg/s^2]
      11: Blade3 pitch angle           [deg]
      12: Blade3 pitch velocity        [deg/s]
      13: Blade3 pitch acceleration    [deg/s^2]
          ...                          
 3nbld+5: Yaw-spring angle             [deg]
 3nbld+6: Yaw-spring velocity          [deg/s]
 3nbld+7: Yaw-spring acceleration      [deg/s^2]
 3nbld+8: Gen torque demand LSS        [kNm]
 3nbld+9: Mechanical power at gen side [kW ]


** qsdof3.dat : dofs of floating body (floater or helicopter) wrt the ref point (i.e. MWL)
                {surge(x),sway(y),heave(z),roll(x),pitch(y),yawz)}
       1: Time (s)
    2- 7: Motions                      [m    , deg    ]
    8-13: Velocities                   [m/s  , deg/s  ]
   14-10: Accelerations                [m/s^2, deg/s^2]


** LOADS_aer.dat: Integrated aerodynamic outputs
              1: Time                  [s  ]
              2: Azimuth               [deg]
              3: Ct                    [ - ]
              4: Cp                    [ - ]
              5: Lamda                 [ - ]
              6: Thrust                [kN ]
              7: Torque                [kNm]
              8: Power                 [kW ]
              9: Omega                 [rpm]
   10   -10+ nb: Thrust per blade      [kN ] (nb:number of blades)
   11+nb-10+2nb: Torque per blade      [kNm] 


** striplXX_ss: aerodynamic output of blade XX at strip ss
   RD : Rotor Disk c.s
   BLD: Blade local c.s.
       1: Time                         [s  ]
       2: Azimuth                      [deg]
       3: Radius_RD                    [m  ]
       4: Angle of Attack              [deg]
       5: Phi_BLD                      [deg]
       6: Weff_BLD                     [m/s]
       7: axial ind. factor a          [ - ]
       8: circ. ind. factor a'         [ - ]
       9: Cl                           [ - ]
      10: Cd                           [ - ]
      11: Cn                           [ - ]
      12: Ct                           [ - ]
      13: Cm                           [ - ]
      14: Fx_BLD/m  or  Ft             [kN/m]
      15: Fz_BLD/m  or  Fn             [kN/m]
      16: My_BLD/m  or  My             [kN ]
      17: Twist                        [deg]
      18:      Pitch_BLD               [deg]
      19:    Torsion_BLD               [deg]
      20:  d(Torsion+Pitch)_BLD        [deg/s]
      21: dd(Torsion+Pitch)_BLD        [deg/s^2]
      22: Mach                         [ - ]
      23: Reynolds                     [ - ]
      24: Loss Coeff                   [ - ]
      26: Flow skewed angle            [deg]
      25: Wake skewed angle            [deg]
      27: Yaw correct Coeff            [ - ]
      28: Axial  ind. vel.             [m/s]
      29: Periph ind. vel.             [m/s]
   30-32: Uwind_Glob x,y,z             [m/s]
   33-35: Ubody_RD   x,y,z             [m/s]
      36: Weffx_BLD                    [m/s]
      37: Weffz_BLD                    [m/s]
      38: Skewed wake induced velocity [m/s]
   39-44: Loads/m RD Fx,Fy,Fz,Mx,My,Mz [kN/m,N]
   45-47: Weffx,y,z_RD no induction    [m/s] 


** velhub.dat : inflow velocity at hub height
       1: Time                         [s  ]
       2: Ux wind velocity x           [m/s]
       3: Uy wind velocity y           [m/s]
       4: Uz wind velocity z           [m/s]


** z_elevation.dat : free surface elevation
       1: Time                         [s  ]
       2: Elevation                    [m  ]
       3: Azimuth (rotor)              [deg]



** eigen.dat
       1: Frequency                    [Hz]
       2: Damping                      [% ]


** modeRZZZ.dat: global mode-shape of the structure of mode ZZZ (sorted)
     1-3: x,y,z coordinates wrt the global c.s.


** modeZZZ_XXX.dat: local mode-shape ZZZ (sorted) per body XXX
       1: Beamwise direction                            [m]
     2-4: x,y,z modal deflections  wrt body local c.s.
     5-7: x,y,z modal rotations    wrt body local c.s.  
