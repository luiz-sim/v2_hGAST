import os,math,smalt,smalt15MW
import numpy as np

stime=650.0


with open('./steps.txt','w') as f:
	f.write('START\n')

# read data
with open('./blade_construction.txt','r') as f:
    data = [i.split() for i in f]
R_tip,r_hup = float(data[0][1]),float(data[1][1])
triax,biax,uniax,balsa,gelcoat,triax_out,carbon,uniax_out,foam,triax_in = [float(data[i][1]) for i in range(4,17)],[float(data[i][2]) for i in range(4,17)],[float(data[i][3]) for i in range(4,17)],[float(data[i][4]) for i in range(4,17)], [float(data[i][5]) for i in range(4,17)], [float(data[i][6]) for i in range(4,17)], [float(data[i][7]) for i in range(4,17)], [float(data[i][8]) for i in range(4,17)], [float(data[i][9]) for i in range(4,17)], [float(data[i][10]) for i in range(4,17)]
profil,Rthickness = [],[]
for i in range(19,29):
    with open('./'+str(data[i][3]),'r') as f:
        ukn = f.readline()
        data2 = [j.split() for j in f]
    profil.append((float(data[i][1]),[data2[j][1] for j in range(0,len(data2))],[data2[j][2] for j in range(0,len(data2))],[data2[j][3] for j in range(0,len(data2))],[data2[j][4] for j in range(0,len(data2))]))
    Rthickness.append(float(data[i][2]))
ranges=[str(data[i][1]) for i in range(31,43)]
data2=[]
for i in range(0,len(ranges)):
    data2.append([])
    with open('./'+str(ranges[i])) as f:
        ukn=f.readline()
        data2[i]=[j.split() for j in f]
radius,x_aer,z_aer,x_sweep,z_prebent,chord,twist,pit_ax,r_thick,ply,kp01,kp02,kp03,kp04,kp05,kp06,kp07,kp08,kp09,kp10,num_webs = [float(data[i][1]) for i in range(45,66)],[float(data[i][2]) for i in range(45,66)],[float(data[i][3]) for i in range(45,66)],[float(data[i][4]) for i in range(45,66)],[float(data[i][5]) for i in range(45,66)],[float(data[i][6]) for i in range(45,66)],[float(data[i][7]) for i in range(45,66)],[float(data[i][8]) for i in range(45,66)],[float(data[i][9]) for i in range(45,66)],[float(data[i][11]) for i in range(45,66)],[float(data[i][12]) for i in range(45,66)],[float(data[i][13]) for i in range(45,66)],[float(data[i][14]) for i in range(45,66)],[float(data[i][15]) for i in range(45,66)],[float(data[i][16]) for i in range(45,66)],[float(data[i][17]) for i in range(45,66)],[float(data[i][18]) for i in range(45,66)],[float(data[i][19]) for i in range(45,66)],[float(data[i][20]) for i in range(45,66)],[float(data[i][21]) for i in range(45,66)], [int(data[i][22]) for i in range(45,66)]
laminates=[]
for i in range(0,len(radius)):
    laminates.append([[float(data2[j][i][1]),float(data2[j][i][2]),float(data2[j][i][3]),float(data2[j][i][4]),float(data2[j][i][3]),float(data2[j][i][2]),float(data2[j][i][1]),float(data2[j][i][5]),float(data2[j][i][6]),float(data2[j][i][7]), float(data2[j][i][8]), float(data2[j][i][9]) , float(data2[j][i][10])] for j in range(0,len(data2))])
CosMod_params = (float(data[67][1]),float(data[68][1]),float(data[72][1]),float(data[70][1]),float(data[71][1]),[float(data[72][1]),float(data[72][2]),float(data[72][3]),float(data[72][4]),float(data[72][5])],[float(data[73][1]),float(data[73][2]),float(data[73][3]),float(data[73][4])],[float(data[74][1]),float(data[74][2]),float(data[74][3]),float(data[74][4]),float(data[74][5])],[float(data[75][1]),float(data[75][2])])

TsaiHill_aim =[0.3201,0.3599,0.4570,0.5282,0.7608,0.8400,0.9276,0.9756,0.9747,0.9392,0.9048,0.9116,0.8021,0.7527,0.7391,0.6230,0.4856,0.3514,0.1830,1.0739,0.0001]
coeff = [1.0 for j in range(0,21)]
############################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
############################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
# TOWER

##INPUTS
z_elevation=10 #gia onshore 0

D_bottom = 5*2
D_up = 3.25*2

Diameter1 = np.linspace(D_bottom, D_up , 21 )
Diameter1 = [round(Diameter1[i],3) for i in range(0,len(Diameter1))]

Height_up = 150.00-2.75 #τοέκανα έτσι και με 21 αντι για 22 διατομες ώστε να υπαρxει και η 22η τελευταια
Height_down = z_elevation

Height1 = np.linspace(Height_down, Height_up , 21)
Height1 = [round(Height1[i],3) for i in range(0,len(Height1))]
Height1_new = [(Height1[i]-z_elevation) for i in range(0,len(Height1))]


thick1=[0 for i in range(0,len(Diameter1))]

for i in range(0,len(Diameter1)):
    if (Diameter1[i]<=10) and (Diameter1[i]>9.926):
        thick1[i]=0.036456
    elif (Diameter1[i]<=9.926) and (Diameter1[i]>=9.443):
        thick1[i]=0.033779
    elif (Diameter1[i]<=9.443) and (Diameter1[i]>8.833):
        thick1[i]=0.032192
    elif (Diameter1[i]<=8.833) and (Diameter1[i]>8.151):
        thick1[i]=0.030708
    elif (Diameter1[i]<=8.151) and (Diameter1[i]>7.390):
        thick1[i]=0.02910
    elif (Diameter1[i]<=7.390) and (Diameter1[i]>=6.909):
        thick1[i]=0.027213
    elif (Diameter1[i]<=6.909) and (Diameter1[i]>6.748):
        thick1[i]=0.024009
    elif (Diameter1[i]<=6.748) and (Diameter1[i]>=6.5):
        thick1[i]=0.0208026
    #elif (Diameter1[i]<=6.572) and (Diameter1[i]>6.5):
    #    thick1[i]=0.022

with open('./steel.txt','w') as f:
	 f.write('num   S355[m]     [-]       [-]    [-]\n')
	 f.write('  1    %.4f     0.0000    0.0000  0.0000\n'%(0.5*thick1[0])) 
	 f.write('  2    %.4f     0.0000    0.0000  0.0000\n'%(0.5*thick1[1])) 
	 f.write('  3    %.4f     0.0000    0.0000  0.0000\n'%(0.5*thick1[2])) 
	 f.write('  4    %.4f     0.0000    0.0000  0.0000\n'%(0.5*thick1[3])) 
	 f.write('  5    %.4f     0.0000    0.0000  0.0000\n'%(0.5*thick1[4])) 
	 f.write('  6    %.4f     0.0000    0.0000  0.0000\n'%(0.5*thick1[5])) 
	 f.write('  7    %.4f     0.0000    0.0000  0.0000\n'%(0.5*thick1[6])) 
	 f.write('  8    %.4f     0.0000    0.0000  0.0000\n'%(0.5*thick1[7])) 
	 f.write('  9    %.4f     0.0000    0.0000  0.0000\n'%(0.5*thick1[8])) 
	 f.write(' 10    %.4f     0.0000    0.0000  0.0000\n'%(0.5*thick1[9])) 
	 f.write(' 11    %.4f     0.0000    0.0000  0.0000\n'%(0.5*thick1[10])) 
	 f.write(' 12    %.4f     0.0000    0.0000  0.0000\n'%(0.5*thick1[11])) 
	 f.write(' 13    %.4f     0.0000    0.0000  0.0000\n'%(0.5*thick1[12])) 
	 f.write(' 14    %.4f     0.0000    0.0000  0.0000\n'%(0.5*thick1[13])) 
	 f.write(' 15    %.4f     0.0000    0.0000  0.0000\n'%(0.5*thick1[14])) 
	 f.write(' 16    %.4f     0.0000    0.0000  0.0000\n'%(0.5*thick1[15])) 
	 f.write(' 17    %.4f     0.0000    0.0000  0.0000\n'%(0.5*thick1[16])) 
	 f.write(' 18    %.4f     0.0000    0.0000  0.0000\n'%(0.5*thick1[17])) 
	 f.write(' 19    %.4f     0.0000    0.0000  0.0000\n'%(0.5*thick1[18])) 
	 f.write(' 20    %.4f     0.0000    0.0000  0.0000\n'%(0.5*thick1[19])) 
	 f.write(' 21    %.4f     0.0000    0.0000  0.0000\n'%(0.5*thick1[20])) 

TsaiHill_tower_aim = [0.8526, 0.9084, 0.8663, 0.9210, 0.8701, 0.9228, 0.8625, 0.9092, 0.8359, 0.8715, 0.7824, 0.6872, 0.6906, 0.5855, 0.5646, 0.4459, 0.3985, 0.3003, 0.2619, 0.1748, 0.1279]
coeff_tower = [1.0 for i in range(0,21)]

# read data
with open('./tower_construction.txt','r') as f:
    data_tower = [i.split() for i in f]
R_tip_tower,r_hup_tower = float(data_tower[0][1]),float(data_tower[1][1])
triax_tower,biax_tower,uniax_tower,balsa_tower = [float(data_tower[i][1]) for i in range(4,17)],[float(data_tower[i][2]) for i in range(4,17)],[float(data_tower[i][3]) for i in range(4,17)],[float(data_tower[i][4]) for i in range(4,17)]
profil_tower,Rthickness_tower = [],[]
for i in range(19,26):
    with open('./'+str(data_tower[i][3]),'r') as f:
        ukn_tower = f.readline()
        data_tower2 = [j.split() for j in f]
    profil_tower.append((float(data_tower[i][1]),[data_tower2[j][1] for j in range(0,len(data_tower2))],[data_tower2[j][2] for j in range(0,len(data_tower2))],[data_tower2[j][3] for j in range(0,len(data_tower2))],[data_tower2[j][4] for j in range(0,len(data_tower2))]))
    Rthickness_tower.append(float(data_tower[i][2]))
ranges_tower=[str(data_tower[i][1]) for i in range(28,40)]
data_tower2=[]
for i in range(0,len(ranges_tower)):
    data_tower2.append([])
    with open('./'+str(ranges_tower[i])) as f:
        ukn_tower=f.readline()
        data_tower2[i]=[j.split() for j in f]
radius_tower,x_aer_tower,z_aer_tower,x_sweep_tower,z_prebent_tower,chord_tower,twist_tower,pit_ax_tower,r_thick_tower,ply_tower,kp_tower01,kp_tower02,kp_tower03,kp_tower04,kp_tower05,kp_tower06,kp_tower07,kp_tower08,kp_tower09,kp_tower10,num_webs_tower = [float(data_tower[i][1]) for i in range(42,63)],[float(data_tower[i][2]) for i in range(42,63)],[float(data_tower[i][3]) for i in range(42,63)],[float(data_tower[i][4]) for i in range(42,63)],[float(data_tower[i][5]) for i in range(42,63)],[float(data_tower[i][6]) for i in range(42,63)],[float(data_tower[i][7]) for i in range(42,63)],[float(data_tower[i][8]) for i in range(42,63)],[float(data_tower[i][9]) for i in range(42,63)],[float(data_tower[i][11]) for i in range(42,63)],[float(data_tower[i][12]) for i in range(42,63)],[float(data_tower[i][13]) for i in range(42,63)],[float(data_tower[i][14]) for i in range(42,63)],[float(data_tower[i][15]) for i in range(42,63)],[float(data_tower[i][16]) for i in range(42,63)],[float(data_tower[i][17]) for i in range(42,63)],[float(data_tower[i][18]) for i in range(42,63)],[float(data_tower[i][19]) for i in range(42,63)],[float(data_tower[i][20]) for i in range(42,63)],[float(data_tower[i][21]) for i in range(42,63)],[int(data_tower[i][22]) for i in range(42,63)]
laminates_tower=[]
for i in range(0,len(radius_tower)):
    laminates_tower.append([[float(data_tower2[j][i][1]),float(data_tower2[j][i][2]),float(data_tower2[j][i][3]),float(data_tower2[j][i][4]),float(data_tower2[j][i][3]),float(data_tower2[j][i][2]),float(data_tower2[j][i][1])] for j in range(0,len(data_tower2))])

xcm_tower, zcm_tower, EAx_tower, EAz_tower, EIxz_tower, GAxx_tower, GAzz_tower = [0 for i in range(0,len(radius_tower))], [0 for i in range(0,len(radius_tower))], [0 for i in range(0,len(radius_tower))], [0 for i in range(0,len(radius_tower))], [0 for i in range(0,len(radius_tower))], [0 for i in range(0,len(radius_tower))], [0 for i in range(0,len(radius_tower))]
GzAz = [0 for i in range(0,len(radius))]

#for ITER in [0,1]:
args = (profil,radius,chord,twist,pit_ax,kp01,kp02,kp03,kp04,kp05,kp06,kp07,kp08,kp09,kp10,[[[coeff[i]*laminates[i][j][k] for k in range(0,len(laminates[i][j]))] for j in range(0,len(laminates[i]))] for i in range(0,len(laminates))],num_webs,ply,triax,biax,uniax,balsa, gelcoat,triax_out,carbon,uniax_out,foam,triax_in)
K6x6,ukn,M6x6,MC = smalt15MW.STIFFMASS(args)

args_tower = (profil_tower,radius_tower,chord_tower,twist_tower,pit_ax_tower,kp_tower01,kp_tower02,kp_tower03,kp_tower04,kp_tower05,kp_tower06,kp_tower07,kp_tower08,kp_tower09,kp_tower10,[[[coeff_tower[i]*laminates_tower[i][j][k] for k in range(0,len(laminates_tower[i][j]))] for j in range(0,len(laminates_tower[i]))] for i in range(0,len(laminates_tower))],num_webs_tower,ply_tower,triax_tower,biax_tower,uniax_tower,balsa_tower)
K6x6_tower,ukn_tower,M6x6_tower,MC_tower = smalt.STIFFMASS(args_tower)

for (case,vel,pitch) in [('ult13R1_15MW_floating.dir',13.0, 4.13),]:
#for (case,vel,pitch) in [('',13.0, 8.25),]:
    with open('./'+str(case)+'/dfile_el.inp','w') as f:
        f.write('** 10MW DTU (FLOATING) Reference WT input file\n')
        f.write('** Application Parameters\n')
        f.write(' 3   ! Substructure   [0:onshore,1:monopile ,2:jacket    ,3:floating(rig),4:floating(flex)]\n')
        f.write(' 0   ! Application    [0:HAWT   ,1:VAWT     ,2:Helicopter]\n')
        f.write(' 0   ! Elasticity     [0:FEM    ,1:modal                 ]\n')
        f.write(' 1   ! Moorings Truss [0:no     ,1:dicoupled,2:coupled   ]\n')
        f.write(' 0   ! Modal Analysis [0:no(td) ,1:no grav  ,2:gravity   ]\n')
        f.write(' 1   ! Modal Damping  [0:no     ,1:calculate,2:read      ]\n')
        f.write(' 1   ! Aerodynamics   [0:no     ,1:BEM      ,2:Vortex    ]\n')
        f.write(' 1   ! Beam Type      [          1:Timo 1st ,2:Euler 2nd ]\n')
        f.write(' 0   ! Foundation     [0:AF     ,1:CS       ,2:DS        ]\n')
        f.write('10   ! Controller     [0:no     ,1:enabled  ,2:iddling   ]\n')
        f.write('** Bodies/Sub-Bodies/Elements Parameters\n')
        f.write(' 5            ! Total Number of Bodies [jacket not icluded]\n')
        f.write(' 3            ! Total Number of blades\n')
        f.write(' 7  1  0  1   ! Blades--> subbodies number    AND IDOF_el 0:no elasticity, 1:elasticiy enabled, 2:monopile modelling\n')
        f.write(' 1  1  0  1   ! Shaft -->    >>      >>       AND    >>          >>            >>                      >>\n')
        f.write(' 1  1  0  1   ! Tower -->    >>      >>       AND    >>          >>            >>                      >>\n')
        f.write(' 1  3         ! Type of Loads 1:aero, 2:hydro AND Total number of elements of SubBody 1 (blades)\n')
        f.write(' 1  3         !           >>                  AND        >>         >>        SubBody 2 (blades)\n')
        f.write(' 1  3         !           >>                  AND        >>         >>        SubBody 3 (blades)\n')
        f.write(' 1  3         !           >>                  AND        >>         >>        SubBody 4 (blades)\n')
        f.write(' 1  3         !           >>                  AND        >>         >>        SubBody 5 (blades)\n')
        f.write(' 1  3         !           >>                  AND        >>         >>        SubBody 6 (blades)\n')
        f.write(' 1  2         !           >>                  AND        >>         >>        SubBody 7 (blades)\n')
        f.write(' 0  2         !           >>                  AND        >>         >>        SubBody 1 (shaft )\n')
        f.write(' 0 21         !           >>                  AND        >>         >>        SubBody 1 (tower )\n')
        f.write('** Initial Conditions\n')
        f.write('  1              ! Initial static solution [0:no,1:yes]\n')
        f.write('000.00           ! Blade 1: Azimuthal position                         [deg]\n')
        f.write(' %.2f   0.0000   ! Blade 1: Pitch [for heli also qc and qs], Imbalance [deg]\n'%(pitch))
        f.write(' %.2f            ! Rotors Speed [LSS]                                  [RPM]\n'%(7.50*121.119/R_tip))
        f.write('  0.0            ! HSS Torque  [not always needed]                     [Nm ]\n')	
        f.write('0.0       1       ! Initial value surge   ux   AND   active dof [1] or not [0] [ m ]\n')
        f.write('0.0       1       ! Initial value sway    uy           >>     >>         >>    [ m ]\n') 
        f.write('0.0       1       ! Initial value heave   uz           >>     >>         >>    [ m ]\n') 
        f.write('0.0       1       ! Initial value roll    θx           >>     >>         >>    [deg]\n') 
        f.write('0.0       1       ! Initial value pitch   θy           >>     >>         >>    [deg]\n') 
        f.write('0.0       1       ! Initial value yaw     θz           >>     >>         >>    [deg]\n') 
        f.write('** Geometry parameters\n')
        f.write('  %.3f   ! Tower top    Height                    [ m ]\n'%((z_elevation+max(Height1_new)+2.75)))
        f.write('  %.3f    ! Tower bottom Height                    [ m ]\n'%(z_elevation))
        f.write('  0.0       ! Shaft offset from tower top            [ m ]\n')
        f.write('  11.320     ! Shaft length                           [ m ]\n')
        f.write('  %.3f      ! Distance from hub to blade root        [ m ]\n'%(r_hup))
        f.write('  0.00      ! Yaw  angle of the nacelle              [deg]\n')
        f.write('  6.000     ! Tilt angle of the rotor                [deg]\n')
        f.write(' -4.000     ! Cone angle of the blade                [deg]\n')
        f.write('  %.1f      ! Rotors Speed (nominal), defines period [RPM]\n'%(7.50*121.118/R_tip))
        f.write('  1.000     ! Gearbox Ratio                          [ - ]\n')
        f.write('  0         ! Yaw mechanism [0:no,1:yes]             [ - ]\n')
        f.write('99999999999 ! Stiffness of yaw mechanism             [Nm/rad]\n')
        f.write('99999999999 ! Damping   of yaw mechanism             [Nm/(rad/s)]\n')
        f.write('** Damping coefficients\n')
        f.write('  0.0E-03      ! Blade Mass      Rayleigh coefficient\n')
        f.write('  0.0E-05      ! Blade Stiffness Rayleigh coefficient\n')
        f.write('  0.00477465   ! Blade modal damping ratio\n')
        f.write('  0.05  20.0   ! Modal damping ratio and frequency [Hz] for higher modes\n')
        f.write('** Simulation/Solution parameters\n')
        f.write('%.2f       ! Simulation time [sec]\n'%(stime))
        f.write('  0.05     ! Time-step [sec]\n')
        f.write('  0.26     ! Bita coefficient for Newmark Integration\n')
        f.write('  0.50     ! Gama coefficient for Newmark Integration\n')
        f.write('  0        ! Rerun option [0:no,1:continue from backup file]\n')
        f.write('999999     ! Timesteps to keep Backup\n')
        f.write('  1.e-05   ! Maximun error convergance criterion\n')
        f.write(' 15        ! Maximun number of Iterations per time-step\n')
        f.write('** Concentrated masses\n')
        f.write('  0   ! Number of Concentrated masses\n')
        f.write('** nb    Hta      Xcm   Ycm    Zcm    Mass[kg]       IPx      IPy         IPz   Ixy   Ixz   Iyz  [local body c.s.]\n')
        f.write(' 5   %.3f  6.210  0.00   0.00   1070000.00   7.900629e6    7.558513e6  1.013478e7   0.    0.    0.  ! Nacelle [3.45]\n'%(max(Height1_new)+2.75-0.35))
        f.write(' 4    11.3200  0.00   0.00   0.00   178832.00         0.       7.277949e5      0.       0.    0.    0.  ! Hub\n')
        f.write(' 4     0.0000  0.00   6.32   0.00        0.           0.       1.715930e7      0.       0.    0.    0.  ! Generator\n')
        f.write('** Mechanical losses\n')
        f.write('  0                 ! Number of entries\n')
        f.write('**INPUT_TRQ [Nm]   LOSSES_TRQ [Nm]\n')
    with open('./'+str(case)+'/rotorin.inp','w') as f:
        f.write('** AERO input file for INNWIND 10MW RWT\n')
        f.write(' %.2f    ! Hub Velocity     [m/s]\n'%(vel))
        f.write('000.00   ! Wind yaw         [deg]\n')
        f.write('  8.00   ! Wind inclination [deg]\n')
        f.write('  1.225  ! Air Density      [kg/m^3]\n')
        f.write('340.0    ! Speed of sound   [m/s]\n')
        f.write('  0.14   ! Shear Exponent   [ - ]\n')
        f.write('  0.00   ! Veer             [deg/m] slope of yaw angle / m wrt hub height\n')
        f.write('  7     ! wind case        [0:uniform,1:defined wind scenario,2:EOG,3x:ECD,4x:EDC,5x:EWS,6:NTM,7:ETM,8:EWM turb]\n')
        f.write(' 50.     ! Vref             [m/s] I:50.00,II:42.50,III:37.50\n')
        f.write('  0.16   ! TIref            [ - ] A:0.16,B:0.14,C:0.12\n')
        f.write(' 55.     ! Time_gust        [sec]\n')
        f.write('  1      ! Turbulent inflow [0:no turb,1:disk,2:rectangular]\n')
        f.write(' 0.0     ! Time shift of turbulent box [sec]\n\n')
        f.write('  0      ! Parked            [0:N,1:Y]                        --- default 0---\n')
        f.write('  0      ! Disable BEM       [0:N,1:Y]                        --- default 0---\n')
        f.write('  1      ! Dynamic Wake      [0:N,1:Y]                        --- default 1---\n')
        f.write('  1      ! Dynamic Stall     [0:N,1:ONERA,2:Beddoes-Leishman] --- default 1---\n')
        f.write('  7      ! Strip to start DS\n\n')
        f.write(' %.3f    ! BEM Radius        [m] [radius after which BEM is solved]\n'%(10.0*R_tip/89.166))
        f.write(' %.3f    ! Rotor Radius      [m]\n'%(121.119*R_tip/121.119))
        f.write(' 30      ! Strips number (blade elements)\n')
        f.write(' 1       ! Tip-losses        [0:N,1:Y] --- default 1 ---\n')
        f.write(' 0       ! Hub-losses        [0:N,1:Y] --- default 0 ---\n')
        f.write(' 1       ! Tower Shadow      [0:N,1:Y] --- default 1 ---\n')
        f.write(' %.3f    ! Tower bottom R    [m]\n'%(Diameter1[0]/2))
        f.write(' %.3f    ! Tower top    R    [m]\n\n'%(Diameter1[-1]/2))
        f.write(' 0     # ! INDXFlap  [0:Without flaps,1:IFC or IFC+IPC,2:IPC only]\n')
    with open('./'+str(case)+'/profilb.inp','w') as f:
        f.write('** Input geometry for INNWIND 10MW wind turbine\n')
        f.write('10   ! NSPANB\n')
        f.write('r         profil\n')
        f.write('0.00   "CYLINDER.inp"\n')
        f.write('%.3f   "CYLINDER.inp"\n'%(r_hup))
        f.write('%.3f   "FFA-W3-500.inp"\n'%(18.168*R_tip/121.118))
        f.write('%.3f   "FFA-W3-360.inp"\n'%(29.698*R_tip/121.118))
        f.write('%.3f   "FFA-W3-330.inp"\n'%(39.824*R_tip/121.118))
        f.write('%.3f   "FFA-W3-301.inp"\n'%(53.195*R_tip/121.118))
        f.write('%.3f   "FFA-W3-270.inp"\n'%(65.126*R_tip/121.118))
        f.write('%.3f   "FFA-W3-241.inp"\n'%(77.298 *R_tip/121.118))
        f.write('%.3f   "FFA-W3-211.inp"\n'%(93.468 *R_tip/121.118))
        f.write('%.3f   "FFA-W3-211.inp"\n'%(123.00*R_tip/121.118))
    with open('./'+str(case)+'/controls.inp','w') as f:
        lam = 7.0+(8.0-7.0)*(R_tip-0.90*121.119)/(1.10*121.119-0.90*121.119)
        f.write('**Input for Generic controller of the DTU 10MW RWT **\n')
        f.write('15957.4                !  1 ; Rated power, including electrical losses [kW]\n')
        f.write('    1                  !  2 ; LSS(1) or HSS(2) reference [-]\n')
        f.write('    0.5236 !=5.0 rpm   !  3 ; Minimum rotor (LSS) speed [rad/s]\n')
        f.write('    %.4f   !=7.5 rpm   !  4 ; Rated rotor (LSS) speed [rad/s]\n'%(0.7854*121.119/R_tip))
        f.write('   2.1765444e7         !  5 ; Maximum allowable generator torque [Nm]\n')
        f.write('   00.                 !  6 ; Minimum pitch angle, PitchMin [deg], if |PitchMin|>90, then a table of <wsp,PitchMin> is read from a file named wptable.n, where n=int(PitchMin)\n')
        f.write('   90.      #82.00     !  7 ; Maximum pitch angle [deg]\n')
        f.write('   2.                 !  8 ; Maximum pitch velocity operation [deg/s]\n')
        f.write('**Partial load control parameters\n')
        f.write('  %.2f        !  9 ; Optimal Cp tracking K factor [Nm/(rad/s)^2], Qg=K*Omega^2, K=eta*0.5*rho*A*Cp_opt*R^3/lambda_opt^3\n'%(0.50*1.225*math.pi*0.45*(R_tip**5.0)/(lam**3.0)))
        f.write(' 1.64438e8    ! 10 ; Proportional gain of torque controller [Nm/(rad/s)]\n')
        f.write(' 3.68998e7    ! 11 ; Integral gain of torque controller [Nm/rad]\n')
        f.write('  0.          ! 12 ; Differential gain of torque controller [Nm/(rad/s^2)]\n')
        f.write('**Full load control parameters\n')
        f.write('  1                 ! 13 ; Generator control switch [1=constant power, 2=constant torque]\n')
        f.write('  1.177141   !/ 3   ! 14 ; Proportional gain of pitch controller [rad/(rad/s)]\n')
        f.write('  0.421192   !/10   ! 15 ; Integral gain of pitch controller [rad/rad]\n')
        f.write('  0.                ! 16 ; Differential gain of pitch controller [rad/(rad/s^2)]\n')
        f.write('  0.4e-8            ! 17 ; Proportional power error gain [rad/W]\n')
        f.write('  0.4e-8            ! 18 ; Integral power error gain [rad/(Ws)]\n')
        f.write('12.03781            ! 19 ; Coefficient of linear term in aerodynamic gain scheduling, KK1 [deg]\n')
        f.write('683.03643           ! 20 ; Coefficient of quadratic term in aerodynamic gain scheduling, KK2 [deg^2], (if zero, KK1=pitch angle at double gain)\n')
        f.write('  1.2               ! 21 ; Coefficinet of nonlinear gain scheduling [-]\n')
        f.write('**Filters\n')
        f.write('  1.25   10P     ! 22 ; lowpass2 freq of speed & power error used in pitch PID [Hz]\n')
        f.write('  0.7           ! 23 ; lowpass2 damp of speed & power error used in pitch PID [-]\n')
        f.write('  3.33  free2   ! 28 ; notch2 freq1 of speed error used in pitch PID [Hz]\n')
        f.write('  0.2           ! 29 ; notch2 damp1 of speed error used in pitch PID [-]\n')
        f.write('  0.375  3P      ! 24 ; notch2 freq2 of speed error used in pitch PID [Hz]\n')
        f.write('  0.2           ! 25 ; notch2 damp2 of speed error used in pitch PID [-]\n')
        f.write('  1.21  free1   ! 26 ; notch2 freq3 of speed error used in pitch PID [Hz]\n')
        f.write('  0.2           ! 27 ; notch2 damp3 of speed error used in pitch PID [-]\n')
        f.write('  2.5   2.5P    ! 30 ; lowpass1 time lag of omega used in GS_NL, Qmax=P/w_LP1 [s]\n')
        f.write('  8.00  1P      ! 31 ; lowpass1 time lag of pitch used in GS, (derate) [s]\n')
        f.write('  8.00  1P      ! 32 ; lowpass1 time lag of power used in GS_IPC, (min pitch) [s]\n')
        f.write('  2.5   2.5P    ! 33 ; lowpass1 time lag of torque used in Torque limits [s]\n')
        f.write('**Drive Train Damper\n')
        f.write('  1.0e7   ! 34 ; DTD1_gain [Nm/(rad/s^2)] 400*113.2^2=5125696.\n')
        f.write('  1.82    ! 35 ; DTD1_freq 1st collective edge [Hz]\n')
        f.write('  0.2     ! 36 ; DTD1_damp [-]\n')
        f.write('  2.0e7   ! 37 ; DTD2_gain [Nm/(rad/s^2)]\n')
        f.write('  3.35    ! 38 ; DTD2_freq 2nd collective edge [Hz]\n')
        f.write('  0.2     ! 39 ; DTD2_damp [-]\n')
        f.write('  1.0e8   ! 40 ; DTD maximum Torque limit LSS [Nm] ~138kWatt\n')
        f.write('**Out-of-plane moments (IPC, IFC)\n')
        f.write('  2             ! 41 ; Number of harmonic (0:no,1:1P,2:2P) [-]\n')
        f.write('  1.            ! 42 ; GAIN_Pitch (default=1) [-]\n')
        f.write('  0.            ! 43 ; GAIN_Flap (default=1) [-]\n')
        f.write('  0.01          ! 44 ; GS (<=1) at 80% of Prated LP1 [-]\n')
        f.write('  0.  15.00e-11 ! 45 ; Ki (1P)\n')
        f.write('  0.  15.00e-10 ! 46 ; Ki (2P)\n')
        f.write('  0.  37.50e-11 ! 47 ; Kp (1P)\n')
        f.write('  0.  37.50e-10 ! 48 ; Kp (2P)\n')
        f.write('  5.            ! 49 ; Maximum Pitch Amplitude [deg]\n')
        f.write(' 10.            ! 50 ; Maximum Flap Amplitude [deg]\n')
        f.write(' 20             ! 51 ; time step to start IPC [-]\n')
        f.write('**Tower Acceleration (Gen,CPC,CFC)\n')
        f.write('  0.          ! 52 ; GAIN_Gen (default=1) [-]\n')
        f.write('  0.          ! 53 ; GAIN_Pitch (default=1) [-]\n')
        f.write('  0.          ! 54 ; GAIN_Flap (default=1) [-]\n')
        f.write('  1.   0.01   ! 55 ; GS (<=1) at 80% of Prated LP1 [-]\n')
        f.write(' 20.e-10      ! 56 ; Ki\n')
        f.write(' 50.e-10      ! 57 ; Kp\n')
        f.write('  1.0e5       ! 58 ; Maximum GenTorque Amplitude LSS [Nm] ~138kWatt\n')
        f.write('  5.0         ! 59 ; Maximum Pitch Amplitude [deg]\n')
        f.write(' 10.0         ! 60 ; Maximum Flap Amplitude [deg]\n')
        f.write(' 20           ! 61 ; time step to start TTa control [-]\n')
        f.write('**Generator Lag\n')
        f.write('  0.02          ! 62 ; generator lag [sec]\n')
        f.write('**Pitch Actuator\n')
        f.write('  2             ! 63 ; itype ## NOT USED\n')
        f.write('  1.0   **1.6   ! 64 ; frequency of pitch actuator [Hz]\n')
        f.write('  0.7   **0.8   ! 65 ; damping factor of pitch actuator [-]\n')
        f.write(' -5.            ! 66 ; lower angle limit of pitch actuator [deg]\n')
        f.write(' 90.0           ! 67 ; upper angle limit of pitch actuator [deg]\n')
        f.write(' 10.0           ! 68 ; speed limit of pitch actuator [deg/s]\n')
        f.write(' 15.0           ! 69 ; acceleration limit of pitch actuator [deg/s^2]\n')
        f.write('**Flap Actuator\n')
        f.write('  1      ! 70 ; itype ## NOT USED\n')
        f.write('  0.1    ! 71 ; time lag of flap  actuator [Hz]\n')
        f.write('-10.     ! 72 ; lower angle limit of flap actuator [deg]\n')
        f.write(' 10.0    ! 73 ; upper angle limit of flap actuator [deg]\n')
        f.write(' 20.0    ! 74 ; speed limit of flap actuator [deg/s]\n\n\n')
        f.write('**Expert\n')
        f.write(' 1.02    ! vs_Reg1_lim\n')
        f.write(' 0.5     ! vs_pit_switch\n\n')
        f.write('** Pending\n')
        f.write('Start-up\n')
        f.write('Shut down\n')
        f.write('Brake\n')
        f.write('MLoss\n')
        f.write('Faults (blade stuck, emergency shut-down, network loss, short circuit)\n\n')
        f.write('monitoring/alarms\n')
        f.write('IPC 2P\n')
        f.write('IPC: (3P, 4P,...)??\n')
        f.write('check TTa, IFC, generator\n')
    with open('./'+str(case)+'/machi.inp','w') as f:
        f.write('***Beam Characteristics of INNWIND 10MW wind turbine ***\n***BLADE\n1.00000, 2 !mass imbalance, iread\n**num   Y     X[swe]     Z[ben]        m          x_cg      z_cg       ρΙxx          ρIzz          ρIxz          GxA             K12             GxzA          k14          GxAx             k16            EA              k23          -EAx            EAy            EAz               GzA             k34           -GzAz            k36            EIxx            -EIxy          -EIxz               GI              -EIyz           EIzz          D\n')
        for i in range(1,len(radius)):
            f.write(' %d   %.3f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f    %.5f   %.5f   %.5f   %.5f   %.5f   0.0\n'%(i,radius[i-1]*(R_tip-r_hup),x_sweep[i-1],z_prebent[i-1],M6x6[i-1][0][0],MC[i-1][0],MC[i-1][1],M6x6[i-1][3][3],M6x6[i-1][5][5],-M6x6[i-1][3][5],K6x6[i-1][0][0],K6x6[i-1][0][1],-K6x6[i-1][0][2],K6x6[i-1][0][3],K6x6[i-1][0][4],-K6x6[i-1][0][5],K6x6[i-1][1][1],-K6x6[i-1][1][2],K6x6[i-1][1][3],K6x6[i-1][1][4],-K6x6[i-1][1][5],K6x6[i-1][2][2],-K6x6[i-1][2][3],GzAz[i-1] ,K6x6[i-1][2][5],K6x6[i-1][3][3],K6x6[i-1][3][4],-K6x6[i-1][3][5],K6x6[i-1][4][4],-K6x6[i-1][4][5],K6x6[i-1][5][5]))
            f.write('      %.3f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f    %.5f   %.5f   %.5f   %.5f   %.5f   0.0\n'%(  radius[i]*(R_tip-r_hup),  x_sweep[i],  z_prebent[i],  M6x6[i][0][0],  MC[i][0],  MC[i][1],  M6x6[i][3][3],  M6x6[i][5][5],  -M6x6[i][3][5],  K6x6[i][0][0],  K6x6[i][0][1],  -K6x6[i][0][2],  K6x6[i][0][3],  K6x6[i][0][4],  -K6x6[i][0][5],  K6x6[i][1][1],  -K6x6[i][1][2],  K6x6[i][1][3],  K6x6[i][1][4],  -K6x6[i][1][5],  K6x6[i][2][2],  -K6x6[i][2][3],  GzAz[i] ,  K6x6[i][2][5],  K6x6[i][3][3],  K6x6[i][3][4],  -K6x6[i][3][5],  K6x6[i][4][4],  -K6x6[i][4][5],  K6x6[i][5][5]))
        f.write('***BLADE\n1.00000, 2 !mass imbalance, iread\n**num   Y     X[swe]     Z[ben]        m          x_cg      z_cg       ρΙxx          ρIzz          ρIxz          GxA             K12             GxzA          k14          GxAx             k16            EA              k23          -EAx            EAy            EAz               GzA             k34           -GzAz            k36            EIxx            -EIxy          -EIxz               GI              -EIyz           EIzz          D\n')
        for i in range(1,len(radius)):
            f.write(' %d   %.3f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f    %.5f   %.5f   %.5f   %.5f   %.5f   0.0\n'%(i,radius[i-1]*(R_tip-r_hup),x_sweep[i-1],z_prebent[i-1],M6x6[i-1][0][0],MC[i-1][0],MC[i-1][1],M6x6[i-1][3][3],M6x6[i-1][5][5],-M6x6[i-1][3][5],K6x6[i-1][0][0],K6x6[i-1][0][1],-K6x6[i-1][0][2],K6x6[i-1][0][3],K6x6[i-1][0][4],-K6x6[i-1][0][5],K6x6[i-1][1][1],-K6x6[i-1][1][2],K6x6[i-1][1][3],K6x6[i-1][1][4],-K6x6[i-1][1][5],K6x6[i-1][2][2],-K6x6[i-1][2][3],GzAz[i-1] ,K6x6[i-1][2][5],K6x6[i-1][3][3],K6x6[i-1][3][4],-K6x6[i-1][3][5],K6x6[i-1][4][4],-K6x6[i-1][4][5],K6x6[i-1][5][5]))
            f.write('      %.3f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f    %.5f   %.5f   %.5f   %.5f   %.5f   0.0\n'%(  radius[i]*(R_tip-r_hup),  x_sweep[i],  z_prebent[i],  M6x6[i][0][0],  MC[i][0],  MC[i][1],  M6x6[i][3][3],  M6x6[i][5][5],  -M6x6[i][3][5],  K6x6[i][0][0],  K6x6[i][0][1],  -K6x6[i][0][2],  K6x6[i][0][3],  K6x6[i][0][4],  -K6x6[i][0][5],  K6x6[i][1][1],  -K6x6[i][1][2],  K6x6[i][1][3],  K6x6[i][1][4],  -K6x6[i][1][5],  K6x6[i][2][2],  -K6x6[i][2][3],  GzAz[i] ,  K6x6[i][2][5],  K6x6[i][3][3],  K6x6[i][3][4],  -K6x6[i][3][5],  K6x6[i][4][4],  -K6x6[i][4][5],  K6x6[i][5][5]))
        f.write('***BLADE\n1.00000, 2 !mass imbalance, iread\n**num   Y     X[swe]     Z[ben]        m          x_cg      z_cg       ρΙxx          ρIzz          ρIxz          GxA             K12             GxzA          k14          GxAx             k16            EA              k23          -EAx            EAy            EAz               GzA             k34           -GzAz            k36            EIxx            -EIxy          -EIxz               GI              -EIyz           EIzz          D\n')
        for i in range(1,len(radius)):
            f.write(' %d   %.3f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f    %.5f   %.5f   %.5f   %.5f   %.5f   0.0\n'%(i,radius[i-1]*(R_tip-r_hup),x_sweep[i-1],z_prebent[i-1],M6x6[i-1][0][0],MC[i-1][0],MC[i-1][1],M6x6[i-1][3][3],M6x6[i-1][5][5],-M6x6[i-1][3][5],K6x6[i-1][0][0],K6x6[i-1][0][1],-K6x6[i-1][0][2],K6x6[i-1][0][3],K6x6[i-1][0][4],-K6x6[i-1][0][5],K6x6[i-1][1][1],-K6x6[i-1][1][2],K6x6[i-1][1][3],K6x6[i-1][1][4],-K6x6[i-1][1][5],K6x6[i-1][2][2],-K6x6[i-1][2][3],GzAz[i-1] ,K6x6[i-1][2][5],K6x6[i-1][3][3],K6x6[i-1][3][4],-K6x6[i-1][3][5],K6x6[i-1][4][4],-K6x6[i-1][4][5],K6x6[i-1][5][5]))
            f.write('      %.3f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f    %.5f   %.5f   %.5f   %.5f   %.5f   0.0\n'%(  radius[i]*(R_tip-r_hup),  x_sweep[i],  z_prebent[i],  M6x6[i][0][0],  MC[i][0],  MC[i][1],  M6x6[i][3][3],  M6x6[i][5][5],  -M6x6[i][3][5],  K6x6[i][0][0],  K6x6[i][0][1],  -K6x6[i][0][2],  K6x6[i][0][3],  K6x6[i][0][4],  -K6x6[i][0][5],  K6x6[i][1][1],  -K6x6[i][1][2],  K6x6[i][1][3],  K6x6[i][1][4],  -K6x6[i][1][5],  K6x6[i][2][2],  -K6x6[i][2][3],  GzAz[i] ,  K6x6[i][2][5],  K6x6[i][3][3],  K6x6[i][3][4],  -K6x6[i][3][5],  K6x6[i][4][4],  -K6x6[i][4][5],  K6x6[i][5][5]))
        f.write('** SHAFT\n 1.0000, 0 !mass imbalanc, iread\n')
        f.write('**num  y        x        z       m        xcm      zcm      mIxx        mIzz        mIxz          EA         EAx        EAz         EIxx        EIzz         EIxz        GIt         GAx        GAz        GAxx        GAzz       D[m]\n')
        f.write(' 1  0.0000   0.0000   0.0000   1.000   0.0000   0.0000  0.1000E+00  0.1000E+00  0.0000E+00  0.1000E+12  0.0000E+00  0.0000E+00  0.1000E+12  0.1000E+12  0.0000E+00  1.7482E+10  1.0000E+16  1.0000E+16  0.00000000  0.00000000  0.000\n    7.0200   0.0000   0.0000   1.000   0.0000   0.0000  0.1000E+00  0.1000E+00  0.0000E+00  0.1000E+12  0.0000E+00  0.0000E+00  0.1000E+12  0.1000E+12  0.0000E+00  1.7482E+10  1.0000E+16  1.0000E+16  0.00000000  0.00000000  0.000\n')
        f.write(' 2  7.0200   0.0000   0.0000   1.000   0.0000   0.0000  0.1000E+00  0.1000E+00  0.0000E+00  0.1000E+12  0.0000E+00  0.0000E+00  0.1000E+12  0.1000E+12  0.0000E+00  1.7482E+10  1.0000E+16  1.0000E+16  0.00000000  0.00000000  0.000\n    11.320   0.0000   0.0000   1.000   0.0000   0.0000  0.1000E+00  0.1000E+00  0.0000E+00  0.1000E+12  0.0000E+00  0.0000E+00  0.1000E+12  0.1000E+12  0.0000E+00  1.7482E+10  1.0000E+16  1.0000E+16  0.00000000  0.00000000  0.000\n')
        f.write('** Tower  Innwind onshore Tower (starts at z = 10m)\n 1.0000, 0 !mass imbalance, iread\n')
        f.write('**num   y       x      z         m            xcm    zcm          mIzz                   mIxx               mIxz         EA                EAx     EAz           EIxx              EIzz           EIxz        GIt                      GAx                 GAz              GAxx    GAzz    D[m]\n')
        for i in range(1,(len(radius_tower))): #													     							  num                 y                                 x                     z                 m                   xcm              zcm                 mIxx                   mIzz                   mIxz              EA		     EAx		    EAz			      EIxx		EIzz			EIxz			GIt			GAx			GAz	      GAxx	GAzz	D
            f.write(' %d   %.3f   %.3f   %.3f   %.5f     %.3f   %.3f   %.8e         %.8e         %.3f   %.8e        %.3f   %.3f   %.3f   %.3f   %.3f   %.8e        %.8e        %.8e        %.3f   %.3f   %.3f\n'%(i,radius_tower[i-1]*(R_tip_tower-r_hup_tower),x_sweep_tower[i-1],z_prebent_tower[i-1],M6x6_tower[i-1][0][0],xcm_tower[i-1],zcm_tower[i-1],M6x6_tower[i-1][3][3],M6x6_tower[i-1][5][5],-M6x6_tower[i-1][3][5],K6x6_tower[i-1][1][1],EAx_tower[i-1],EAz_tower[i-1],K6x6_tower[i-1][3][3],K6x6_tower[i-1][5][5],EIxz_tower[i-1],K6x6_tower[i-1][4][4], K6x6_tower[i-1][0][0],K6x6_tower[i-1][2][2], GAxx_tower[i-1], GAzz_tower[i-1], Diameter1[i-1] ))
            f.write('      %.3f   %.3f   %.3f   %.5f     %.3f   %.3f   %.8e         %.8e         %.3f   %.8e        %.3f   %.3f   %.3f   %.3f   %.3f   %.8e        %.8e        %.8e        %.3f   %.3f   %.3f\n'%(  radius_tower[i]*(R_tip_tower-r_hup_tower),  x_sweep_tower[i],  z_prebent_tower[i],  M6x6_tower[i][0][0],  xcm_tower[i],  zcm_tower[i],  M6x6_tower[i][3][3],  M6x6_tower[i][5][5],  -M6x6_tower[i][3][5],  K6x6_tower[i][1][1],  EAx_tower[i],  EAz_tower[i],  K6x6_tower[i][3][3],  K6x6_tower[i][5][5],  EIxz_tower[i],  K6x6_tower[i][4][4],   K6x6_tower[i][0][0],  K6x6_tower[i][2][2],   GAxx_tower[i] ,  GAzz_tower[i] ,  Diameter1[i] ))
        f.write('21   %.3f  0.000  0.000     1.000    0.000  0.000 0.1000E+00     0.1000E+00      0.000 0.1000E+13      0.000  0.000 1000000000000.000 1000000000000.000  0.000 0.1000E+13     1.0000E+16     1.0000E+16      0.000  0.000  0.000\n'%(Height1_new[-1]))
        f.write('     %.3f  0.000  0.000     1.000    0.000  0.000 0.1000E+00     0.1000E+00      0.000 0.1000E+13      0.000  0.000 1000000000000.000 1000000000000.000  0.000 0.1000E+13     1.0000E+16     1.0000E+16      0.000  0.000  0.000\n'%((Height1_new[-1]+2.75)))
    with open('./'+str(case)+'/geomb.inp','w') as f:
        f.write('** Input geometry for INNWIND 10MW wind turbine\n41 ! NSPANB\n** r[m] chord[m] twist[deg] x_aer[m] z_aer[m]\n')
        f.write('%.5f   %.5f   %.5f   %.5f   %.5f\n'%(0.0,chord[0],-twist[0],x_aer[0]*np.cos(twist[0]*math.pi/180.0)-z_aer[0]*np.sin(twist[0]*math.pi/180.0),x_aer[0]*np.sin(twist[0]*math.pi/180.0)+z_aer[0]*np.cos(twist[0]*math.pi/180.0)))
        for i in range(1,len(radius)):
            f.write('%.5f   %.5f   %.5f   %.5f   %.5f\n'%(r_hup+(radius[i-1]+radius[i])*(R_tip-r_hup)/2.0,(chord[i-1]+chord[i])/2.0,-(twist[i-1]+twist[i])/2.0,(x_aer[i-1]*np.cos(twist[i-1]*math.pi/180.0)-z_aer[i-1]*np.sin(twist[i-1]*math.pi/180.0)+x_aer[i]*np.cos(twist[i]*math.pi/180.0)-z_aer[i]*np.sin(twist[i]*math.pi/180.0))/2.0,(x_aer[i-1]*np.sin(twist[i-1]*math.pi/180.0)+z_aer[i-1]*np.cos(twist[i-1]*math.pi/180.0)+x_aer[i]*np.sin(twist[i]*math.pi/180.0)+z_aer[i]*np.cos(twist[i]*math.pi/180.0))/2.0))
            f.write('%.5f   %.5f   %.5f   %.5f   %.5f\n'%(r_hup+radius[i]*(R_tip-r_hup),chord[i],-twist[i],x_aer[i]*np.cos(twist[i]*math.pi/180.0)-z_aer[i]*np.sin(twist[i]*math.pi/180.0),x_aer[i]*np.sin(twist[i]*math.pi/180.0)+z_aer[i]*np.cos(twist[i]*math.pi/180.0)))

    #for case in ['ult13R2.dir']:
       # os.system('cd ./'+str(case)+'/ ; ./datclean.inp')
       # os.system('cd ./'+str(case)+'/ ; export OMP_NUM_THREADS=1 ; ~/MANOLO_hGAST/hGAST/hGAST/hGAST')
       # os.system('cd ./'+str(case)+'/ ; ~/hGAST/bin2ascii.out<bin2ascii.inp')

############################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
# BLADE
#    Fx,Fy,Fz,Mx,My,Mz = [],[],[],[],[],[]
#    for (num,cs) in [(0,'01'),(1,'02'),(2,'03'),(3,'04'),(4,'05'),(5,'06'),(6,'07'),(7,'08'),(8,'09'),(9,'10'),(10,'11'),(11,'12'),(12,'13'),(13,'14'),(14,'15'),(15,'16'),(16,'17'),(17,'18'),(18,'19'),(19,'20'),(20,'21')]:
#        fx,fy,fz,mx,my,mz,mc = [],[],[],[],[],[],[]
#        for bla in ['001','002','003']:
#            for case in ['ult13R2.dir']:
#                with open('./'+str(case)+'/loads'+str(bla)+'_'+str(cs)+'.dat','r') as f:
#                    data = [i.split() for i in f]
#                fx.extend([float(data[i][2]) for i in range(501,len(data))]) ; fy.extend([float(data[i][3]) for i in range(501,len(data))]) ; fz.extend([float(data[i][4]) for i in range(501,len(data))])
#                mx.extend([float(data[i][5]) for i in range(501,len(data))]) ; my.extend([float(data[i][6]) for i in range(501,len(data))]) ; mz.extend([float(data[i][7]) for i in range(501,len(data))]) ; mc.extend([float(data[i][8]) for i in range(501,len(data))])
#        Fx.append(fx[mc.index(max(mc))]) ; Fy.append(fy[mc.index(max(mc))]) ; Fz.append(fz[mc.index(max(mc))])
#        Mx.append(mx[mc.index(max(mc))]) ; My.append(my[mc.index(max(mc))]) ; Mz.append(mz[mc.index(max(mc))])
#
#    Loads = [(1350.0*Fx[i],1350.0*Fy[i],1350.0*Fz[i],1350.0*Mx[i],1350.0*My[i],1350.0*Mz[i]) for i in range(0,21)]
#    check,num = True,0
#    while check==True:
#        args1 = (R_tip,r_hup,profil,radius,x_sweep,z_prebent,chord,twist,pit_ax,kp01,kp02,kp03,kp04,kp05,kp06,kp07,kp08,kp09,kp10,[[[(coeff[i]-0.01)*laminates[i][j][k] for k in range(0,len(laminates[i][j]))] for j in range(0,len(laminates[i]))] for i in range(0,len(laminates))],num_webs,ply,triax,biax,uniax,balsa,Loads)
#        TsaiHill1 = smalt.STRESSES(args1)
#        args2 = (R_tip,r_hup,profil,radius,x_sweep,z_prebent,chord,twist,pit_ax,kp01,kp02,kp03,kp04,kp05,kp06,kp07,kp08,kp09,kp10,[[[       coeff[i]*laminates[i][j][k] for k in range(0,len(laminates[i][j]))] for j in range(0,len(laminates[i]))] for i in range(0,len(laminates))],num_webs,ply,triax,biax,uniax,balsa,Loads)
#        TsaiHill2 = smalt.STRESSES(args2)
#        args3 = (R_tip,r_hup,profil,radius,x_sweep,z_prebent,chord,twist,pit_ax,kp01,kp02,kp03,kp04,kp05,kp06,kp07,kp08,kp09,kp10,[[[(coeff[i]+0.01)*laminates[i][j][k] for k in range(0,len(laminates[i][j]))] for j in range(0,len(laminates[i]))] for i in range(0,len(laminates))],num_webs,ply,triax,biax,uniax,balsa,Loads)
#        TsaiHill3 = smalt.STRESSES(args3)
#        coeff_new = [coeff[i]-0.02*(TsaiHill2[i]-TsaiHill_aim[i])/(TsaiHill3[i]-TsaiHill1[i]) if coeff[i]-0.02*(TsaiHill2[i]-TsaiHill_aim[i])/(TsaiHill3[i]-TsaiHill1[i])>0.50 else 0.50 for i in range(0,19)]+[1.0,1.0]
#        if max([abs(coeff[i]-coeff_new[i]) for i in range(0,21)])<0.005 or num>15: check = False
#        coeff,num = [0.25*coeff[i]+0.75*coeff_new[i] for i in range(0,21)],num+1
#
#############################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
## TOWER
#    Fx_tower,Fy_tower,Fz_tower,Mx_tower,My_tower,Mz_tower = [],[],[],[],[],[]
#    for (num,cs) in [(0,'01'),(1,'02'),(2,'03'),(3,'04'),(4,'05'),(5,'06'),(6,'07'),(7,'08'),(8,'09'),(9,'10'),(10,'11'),(11,'12'),(12,'13'),(13,'14'),(14,'15'),(15,'16'),(16,'17'),(17,'18'),(18,'19'),(19,'20'),(20,'21')]:
#        fx_tower,fy_tower,fz_tower,mx_tower,my_tower,mz_tower,mc_tower = [],[],[],[],[],[],[]
#        for tow in ['005']:
#            for case in ['ult13R2.dir']:
#                with open('./'+str(case)+'/loads'+str(tow)+'_'+str(cs)+'.dat','r') as f:
#                    data_tower = [i.split() for i in f]
#                fx_tower.extend([float(data_tower[i][2]) for i in range(501,len(data_tower))]) ; fy_tower.extend([float(data_tower[i][3]) for i in range(501,len(data_tower))]) ; fz_tower.extend([float(data_tower[i][4]) for i in range(501,len(data_tower))])
#                mx_tower.extend([float(data_tower[i][5]) for i in range(501,len(data_tower))]) ; my_tower.extend([float(data_tower[i][6]) for i in range(501,len(data_tower))]) ; mz_tower.extend([float(data_tower[i][7]) for i in range(501,len(data_tower))]) ; mc_tower.extend([float(data_tower[i][8]) for i in range(501,len(data_tower))])
#        Fx_tower.append(fx_tower[mc_tower.index(max(mc_tower))]) ; Fy_tower.append(fy_tower[mc_tower.index(max(mc_tower))]) ; Fz_tower.append(fz_tower[mc_tower.index(max(mc_tower))])
#        Mx_tower.append(mx_tower[mc_tower.index(max(mc_tower))]) ; My_tower.append(my_tower[mc_tower.index(max(mc_tower))]) ; Mz_tower.append(mz_tower[mc_tower.index(max(mc_tower))])
#
#    Loads_tower = [(1350.0*Fx_tower[i],1350.0*Fy_tower[i],1350.0*Fz_tower[i],1350.0*Mx_tower[i],1350.0*My_tower[i],1350.0*Mz_tower[i]) for i in range(0,21)]
#    check,num = True,0
#    while check==True:
#        args_tower1 = (R_tip_tower,r_hup_tower,profil_tower,radius_tower,x_sweep_tower,z_prebent_tower,chord_tower,twist_tower,pit_ax_tower,kp_tower01,kp_tower02,kp_tower03,kp_tower04,kp_tower05,kp_tower06,kp_tower07,kp_tower08,kp_tower09,kp_tower10,[[[(coeff_tower[i]-0.01)*laminates_tower[i][j][k] for k in range(0,len(laminates_tower[i][j]))] for j in range(0,len(laminates_tower[i]))] for i in range(0,len(laminates_tower))],num_webs_tower,ply_tower,triax_tower,biax_tower,uniax_tower,balsa_tower,Loads_tower)
#        TsaiHill_tower1 = smalt.STRESSES(args_tower1)
#        args_tower2 = (R_tip_tower,r_hup_tower,profil_tower,radius_tower,x_sweep_tower,z_prebent_tower,chord_tower,twist_tower,pit_ax_tower,kp_tower01,kp_tower02,kp_tower03,kp_tower04,kp_tower05,kp_tower06,kp_tower07,kp_tower08,kp_tower09,kp_tower10,[[[       coeff_tower[i]*laminates_tower[i][j][k] for k in range(0,len(laminates_tower[i][j]))] for j in range(0,len(laminates_tower[i]))] for i in range(0,len(laminates_tower))],num_webs_tower,ply_tower,triax_tower,biax_tower,uniax_tower,balsa_tower,Loads_tower)
#        TsaiHill_tower2 = smalt.STRESSES(args_tower2)
#        args_tower3 = (R_tip_tower,r_hup_tower,profil_tower,radius_tower,x_sweep_tower,z_prebent_tower,chord_tower,twist_tower,pit_ax_tower,kp_tower01,kp_tower02,kp_tower03,kp_tower04,kp_tower05,kp_tower06,kp_tower07,kp_tower08,kp_tower09,kp_tower10,[[[(coeff_tower[i]+0.01)*laminates_tower[i][j][k] for k in range(0,len(laminates_tower[i][j]))] for j in range(0,len(laminates_tower[i]))] for i in range(0,len(laminates_tower))],num_webs_tower,ply_tower,triax_tower,biax_tower,uniax_tower,balsa_tower,Loads_tower)
#        TsaiHill_tower3 = smalt.STRESSES(args_tower3)
#        coeff_tower_new = [coeff_tower[i]-0.02*(TsaiHill_tower2[i]-TsaiHill_tower_aim[i])/(TsaiHill_tower3[i]-TsaiHill_tower1[i]) if coeff_tower[i]-0.02*(TsaiHill_tower2[i]-TsaiHill_tower_aim[i])/(TsaiHill_tower3[i]-TsaiHill_tower1[i])>0.50 else 0.50 for i in range(0,21)]
#        if max([abs(coeff_tower[i]-coeff_tower_new[i]) for i in range(0,21)])<0.005 or num>15: check = False
#        coeff_tower,num = [0.25*coeff_tower[i]+0.75*coeff_tower_new[i] for i in range(0,21)],num+1
#    with open('./steps.txt','a+') as f:
#        f.write('successfull iteration number %.0f  ...\n'%(ITER))
#
#############################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
## BLADE
#Args = (R_tip,r_hup,profil,radius,x_sweep,z_prebent,chord,twist,pit_ax,kp01,kp02,kp03,kp04,kp05,kp06,kp07,kp08,kp09,kp10,[[[coeff[i]*laminates[i][j][k] for k in range(0,len(laminates[i][j]))] for j in range(0,len(laminates[i]))] for i in range(0,len(laminates))],num_webs,ply,triax,biax,uniax,balsa,CosMod_params)
#Rcost = smalt.COST(Args)
#############################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
#
#with open('./results.txt','w') as f:
#    f.write('#Rcost[$]: %.2f\n'%(Rcost))
#    f.write('#radius[-] coefficients[-] TsaiHill[-] TsaiHill-aim[-]\n')
#    for i in range(0,len(radius)):
#        f.write('%.3f   %.3f   %.4f   %.4f\n'%(radius[i],coeff[i],TsaiHill2[i],TsaiHill_aim[i]))
#
#with open('./ALL_results.txt','a') as f:
#    f.write('#Rcost[$]: %.2f\n'%(Rcost))
#    f.write('#R-Rref[-]: %.6f\n'%(R_tip/89.167))
#    f.write('#H_tower[m]: %.5f\n'%((max(Height1_new)+2.75)))
#    f.write('#radius[-] coefficients[-] TsaiHill[-] TsaiHill-aim[-]\n')
#    for i in range(0,len(radius)):
#        f.write('%.3f   %.3f   %.4f   %.4f\n'%(radius[i],coeff[i],TsaiHill2[i],TsaiHill_aim[i]))
#    f.write('#########################################################################################################\n')
#
#with open('./results_tower.txt','w') as f:
#    f.write('#D_bottom/D_bottom_ref[-]: %.6f\n'%(Diameter1[0]/8.058))
#    f.write('#H_tower[m]: %.5f\n'%((max(Height1_new)+2.75)))
#    f.write('#radius[-] coefficients[-] TsaiHill_Tower[-] TsaiHill_Tower-aim[-] Real-Thickness[m]\n')
#    for i in range(0,len(radius_tower)):
#        f.write('%.3f   %.3f   %.4f   %.4f   %.4f\n'%(radius_tower[i],coeff_tower[i],TsaiHill_tower2[i],TsaiHill_tower_aim[i],coeff_tower[i]*thick1[i]))
#
#with open('./ALL_results_tower.txt','a') as f:
#    f.write('#D_bottom/D_bottom_ref[-]: %.6f\n'%(Diameter1[0]/8.058))
#    f.write('#H_tower[m]: %.5f\n'%((max(Height1_new)+2.75)))
#    f.write('#radius[-] coefficients[-] TsaiHill_Tower[-] TsaiHill_Tower-aim[-] Real-Thickness[m]\n')
#    for i in range(0,len(radius_tower)):
#        f.write('%.3f   %.3f   %.4f   %.4f   %.4f\n'%(radius_tower[i],coeff_tower[i],TsaiHill_tower2[i],TsaiHill_tower_aim[i],coeff_tower[i]*thick1[i]))
#    f.write('#########################################################################################################\n')
#
#with open('./twist_bezier_ALL.gnu','a') as f:
#    f.write('#r   ybez\n')
#    for i in range(0,len(radius)):
#        f.write('%.3f   %.3f\n'%(radius[i],twist[i]))
#    f.write('#########################################################################################################\n')
#
#with open('./chord_ALL.gnu','a') as f:
#    f.write('#r[m]   chord[m]\n')
#    for i in range(0,len(radius)):
#        f.write('%.3f   %.3f\n'%(radius[i],chord[i]))
#    f.write('#########################################################################################################\n')
#
#with open('../twist_bezier.gnu','w') as f:
#    f.write('#r   ybez\n')
#    for i in range(0,len(radius)):
#        f.write('%.3f   %.3f\n'%(radius[i],twist[i]))
#
#with open('../chord.gnu','w') as f:
#    f.write('#r[m]   chord[m]\n')
#    for i in range(0,len(radius)):
#        f.write('%.3f   %.3f\n'%(radius[i],chord[i]))
###########################################################################################################################################################################################################################
#
#with open('../POWER.dir/data.dir/blade_construction.txt','w') as f:
#    f.write('R_tip[m]: %.3f\n'%(R_tip))
#    f.write('R_hup[-]: %.3f\n'%(r_hup))
#    f.write('#################################################\n')
#    f.write('                 TRIAX        BIAX       UNIAX      BALSA\n')
#    f.write('cost[$/Kgr]       %.3f        %.3f       %.3f        %.3f\n'%(triax[0],biax[0],uniax[0],balsa[0]))
#    f.write('density[m3/Kgr]   %.3f        %.3f       %.3f        %.3f\n'%(triax[1],biax[1],uniax[1],balsa[1]))
#    f.write('E11[N/m2]         %.3f        %.3f       %.3f        %.3f\n'%(triax[2],biax[2],uniax[2],balsa[2]))
#    f.write('E22[N/m2]         %.3f        %.3f       %.3f        %.3f\n'%(triax[3],biax[3],uniax[3],balsa[3]))
#    f.write('v12[-]            %.2f        %.2f       %.2f        %.2f\n'%(triax[4],biax[4],uniax[4],balsa[4]))
#    f.write('G12[Nm/m2]        %.3f        %.3f       %.3f        %.3f\n'%(triax[5],biax[5],uniax[5],balsa[5]))
#    f.write('G13[Nm/m2]        %.3f        %.3f       %.3f        %.3f\n'%(triax[6],biax[6],uniax[6],balsa[6]))
#    f.write('G23[Nm/m2]        %.3f        %.3f       %.3f        %.3f\n'%(triax[7],biax[7],uniax[7],balsa[7]))
#    f.write('xc[N/m2]          %.3f        %.3f       %.3f        %.3f\n'%(triax[8],biax[8],uniax[8],balsa[8]))
#    f.write('xt[N/m2]          %.3f        %.3f       %.3f        %.3f\n'%(triax[9],biax[9],uniax[9],balsa[9]))
#    f.write('yc[N/m2]          %.3f        %.3f       %.3f        %.3f\n'%(triax[10],biax[10],uniax[10],balsa[10]))
#    f.write('yt[N/m2]          %.3f        %.3f       %.3f        %.3f\n'%(triax[11],biax[11],uniax[11],balsa[11]))
#    f.write('st[N/m2]          %.3f        %.3f       %.3f        %.3f\n'%(triax[12],biax[12],uniax[12],balsa[12]))
#    f.write('#################################################\n')
#    f.write('num radius[-]   c/t[-]      airfoil\n')
#    f.write(' 1   0.00000   1.00000    cylinder.txt\n')
#    f.write(' 2   0.16650   0.60000   ffa_w3_600.txt\n')
#    f.write(' 3   0.20939   0.48000   ffa_w3_480.txt\n')
#    f.write(' 4   0.29369   0.36000   ffa_w3_360.txt\n')
#    f.write(' 5   0.38831   0.30100   ffa_w3_301.txt\n')
#    f.write(' 6   0.70228   0.24099   ffa_w3_241.txt\n')
#    f.write(' 7   1.00001   0.24099   ffa_w3_241.txt\n')
#    f.write('#################################################\n')
#    f.write('num   range\n')
#    f.write(' 1    tail_v.txt\n')
#    f.write(' 2    tail.txt\n')
#    f.write(' 3    trailing.txt\n')
#    f.write(' 4    caps.txt\n')
#    f.write(' 5    leading.txt\n')
#    f.write(' 6    nose.txt\n')
#    f.write(' 7    leading.txt\n')
#    f.write(' 8    caps.txt\n')
#    f.write(' 9    trailing.txt\n')
#    f.write('10    tail.txt\n')
#    f.write('11    tail_v.txt\n')
#    f.write('12    webs.txt\n')
#    f.write('#################################################\n')
#    f.write('num   radius[-] x-aer[-] y-aer[-] x-sweep[m] z-prebent[m] chord[m] twist[deg]  pit-ax[-]  R.thickness[-]  coefficient[-] ply[deg]  kp-01[-]   kp-02[-]   kp-03[-]   kp-04[-]   kp-05[-]   kp-06[-]   kp-07[-]   kp-08[-]   kp-09[-]   kp-10[-]   num_webs[-]\n')
#    f.write(' 1       %.3f     %.3f     %.3f     %.3f       %.3f      %.3f        %.3f       %.3f         %.3f           %.3f         0.0      %.2f       %.2f       %.2f       %.2f       %.2f      %.2f       %.2f       %.2f        %.2f       %.2f           %d \n'%(radius[0],x_aer[0],z_aer[0],x_sweep[0],z_prebent[0],chord[0],twist[0],pit_ax[0] ,r_thick[0], coeff[0], kp01[0],kp02[0], kp03[0], kp04[0], kp05[0], kp06[0], kp07[0],kp08[0], kp09[0],kp10[0],num_webs[0]))
#    f.write(' 2       %.3f     %.3f     %.3f     %.3f       %.3f      %.3f        %.3f       %.3f         %.3f           %.3f         %.2f      %.2f       %.2f       %.2f       %.2f       %.2f      %.2f       %.2f       %.2f        %.2f       %.2f           %d \n'%(radius[1],x_aer[1],z_aer[1],x_sweep[1],z_prebent[1],chord[1],twist[1],pit_ax[1] ,r_thick[1], coeff[1], ply[1],kp01[1],kp02[1], kp03[1], kp04[1], kp05[1], kp06[1], kp07[1],kp08[1], kp09[1],kp10[1],num_webs[1]))
#    f.write(' 3       %.3f     %.3f     %.3f     %.3f       %.3f      %.3f        %.3f       %.3f         %.3f           %.3f         %.2f      %.2f       %.2f       %.2f       %.2f       %.2f      %.2f       %.2f       %.2f        %.2f       %.2f           %d \n'%(radius[2],x_aer[2],z_aer[2],x_sweep[2],z_prebent[2],chord[2],twist[2],pit_ax[2] ,r_thick[2], coeff[2], ply[2],kp01[2],kp02[2], kp03[2], kp04[2], kp05[2], kp06[2], kp07[2],kp08[2], kp09[2],kp10[2],num_webs[2]))
#    f.write(' 4       %.3f     %.3f     %.3f     %.3f       %.3f      %.3f        %.3f       %.3f         %.3f           %.3f         %.2f      %.2f       %.2f       %.2f       %.2f       %.2f      %.2f       %.2f       %.2f        %.2f       %.2f           %d \n'%(radius[3],x_aer[3],z_aer[3],x_sweep[3],z_prebent[3],chord[3],twist[3],pit_ax[3] ,r_thick[3], coeff[3], ply[3],kp01[3],kp02[3], kp03[3], kp04[3], kp05[3], kp06[3], kp07[3],kp08[3], kp09[3],kp10[3],num_webs[3]))
#    f.write(' 5       %.3f     %.3f     %.3f     %.3f       %.3f      %.3f        %.3f       %.3f         %.3f           %.3f         %.2f      %.2f       %.2f       %.2f       %.2f       %.2f      %.2f       %.2f       %.2f        %.2f       %.2f           %d \n'%(radius[4],x_aer[4],z_aer[4],x_sweep[4],z_prebent[4],chord[4],twist[4],pit_ax[4] ,r_thick[4], coeff[4], ply[4],kp01[4],kp02[4], kp03[4], kp04[4], kp05[4], kp06[4], kp07[4],kp08[4], kp09[4],kp10[4],num_webs[4]))
#    f.write(' 6       %.3f     %.3f     %.3f     %.3f       %.3f      %.3f        %.3f       %.3f         %.3f           %.3f         %.2f      %.2f       %.2f       %.2f       %.2f       %.2f      %.2f       %.2f       %.2f        %.2f       %.2f           %d \n'%(radius[5],x_aer[5],z_aer[5],x_sweep[5],z_prebent[5],chord[5],twist[5],pit_ax[5] ,r_thick[5], coeff[5], ply[5],kp01[5],kp02[5], kp03[5], kp04[5], kp05[5], kp06[5], kp07[5],kp08[5], kp09[5],kp10[5],num_webs[5]))
#    f.write(' 7       %.3f     %.3f     %.3f     %.3f       %.3f      %.3f        %.3f       %.3f         %.3f           %.3f         %.2f      %.2f       %.2f       %.2f       %.2f       %.2f      %.2f       %.2f       %.2f        %.2f       %.2f           %d \n'%(radius[6],x_aer[6],z_aer[6],x_sweep[6],z_prebent[6],chord[6],twist[6],pit_ax[6] ,r_thick[6], coeff[6], ply[6],kp01[6],kp02[6], kp03[6], kp04[6], kp05[6], kp06[6], kp07[6],kp08[6], kp09[6],kp10[6],num_webs[6]))
#    f.write(' 8       %.3f     %.3f     %.3f     %.3f       %.3f      %.3f        %.3f       %.3f         %.3f           %.3f         %.2f      %.2f       %.2f       %.2f       %.2f       %.2f      %.2f       %.2f       %.2f        %.2f       %.2f           %d \n'%(radius[7],x_aer[7],z_aer[7],x_sweep[7],z_prebent[7],chord[7],twist[7],pit_ax[7] ,r_thick[7], coeff[7], ply[7],kp01[7],kp02[7], kp03[7], kp04[7], kp05[7], kp06[7], kp07[7],kp08[7], kp09[7],kp10[7],num_webs[7]))
#    f.write(' 9       %.3f     %.3f     %.3f     %.3f       %.3f      %.3f        %.3f       %.3f         %.3f           %.3f         %.2f      %.2f       %.2f       %.2f       %.2f       %.2f      %.2f       %.2f       %.2f        %.2f       %.2f           %d \n'%(radius[8],x_aer[8],z_aer[8],x_sweep[8],z_prebent[8],chord[8],twist[8],pit_ax[8] ,r_thick[8], coeff[8], ply[8],kp01[8],kp02[8], kp03[8], kp04[8], kp05[8], kp06[8], kp07[8],kp08[8], kp09[8],kp10[8],num_webs[8]))
#    f.write('10       %.3f     %.3f     %.3f     %.3f       %.3f      %.3f        %.3f       %.3f         %.3f           %.3f         %.2f      %.2f       %.2f       %.2f       %.2f       %.2f      %.2f       %.2f       %.2f        %.2f       %.2f           %d \n'%(radius[9],x_aer[9],z_aer[9],x_sweep[9],z_prebent[9],chord[9],twist[9],pit_ax[9] ,r_thick[9], coeff[9], ply[9],kp01[9],kp02[9], kp03[9], kp04[9], kp05[9], kp06[9], kp07[9],kp08[9], kp09[9],kp10[9],num_webs[9]))
#    f.write('11       %.3f     %.3f     %.3f     %.3f       %.3f      %.3f        %.3f       %.3f         %.3f           %.3f         %.2f      %.2f       %.2f       %.2f       %.2f       %.2f      %.2f       %.2f       %.2f        %.2f       %.2f           %d \n'%(radius[10],x_aer[10],z_aer[10],x_sweep[10],z_prebent[10],chord[10],twist[10],pit_ax[10] ,r_thick[10],coeff[10], ply[10],kp01[10],kp02[10],kp03[10],kp04[10], kp05[10], kp06[10],kp07[10],kp08[10],kp09[10],kp10[10],num_webs[10]))
#    f.write('12       %.3f     %.3f     %.3f     %.3f       %.3f      %.3f        %.3f       %.3f         %.3f           %.3f         %.2f      %.2f       %.2f       %.2f       %.2f       %.2f      %.2f       %.2f       %.2f        %.2f       %.2f           %d \n'%(radius[11],x_aer[11],z_aer[11],x_sweep[11],z_prebent[11],chord[11],twist[11],pit_ax[11] ,r_thick[11],coeff[11], ply[11],kp01[11],kp02[11],kp03[11],kp04[11], kp05[11], kp06[11],kp07[11],kp08[11],kp09[11],kp10[11],num_webs[11]))
#    f.write('13       %.3f     %.3f     %.3f     %.3f       %.3f      %.3f        %.3f       %.3f         %.3f           %.3f         %.2f      %.2f       %.2f       %.2f       %.2f       %.2f      %.2f       %.2f       %.2f        %.2f       %.2f           %d \n'%(radius[12],x_aer[12],z_aer[12],x_sweep[12],z_prebent[12],chord[12],twist[12],pit_ax[12] ,r_thick[12],coeff[12], ply[12],kp01[12],kp02[12],kp03[12],kp04[12], kp05[12], kp06[12],kp07[12],kp08[12],kp09[12],kp10[12],num_webs[12]))
#    f.write('14       %.3f     %.3f     %.3f     %.3f       %.3f      %.3f        %.3f       %.3f         %.3f           %.3f         %.2f      %.2f       %.2f       %.2f       %.2f       %.2f      %.2f       %.2f       %.2f        %.2f       %.2f           %d \n'%(radius[13],x_aer[13],z_aer[13],x_sweep[13],z_prebent[13],chord[13],twist[13],pit_ax[13] ,r_thick[13],coeff[13], ply[13],kp01[13],kp02[13],kp03[13],kp04[13], kp05[13], kp06[13],kp07[13],kp08[13],kp09[13],kp10[13],num_webs[13]))
#    f.write('15       %.3f     %.3f     %.3f     %.3f       %.3f      %.3f        %.3f       %.3f         %.3f           %.3f         %.2f      %.2f       %.2f       %.2f       %.2f       %.2f      %.2f       %.2f       %.2f        %.2f       %.2f           %d \n'%(radius[14],x_aer[14],z_aer[14],x_sweep[14],z_prebent[14],chord[14],twist[14],pit_ax[14] ,r_thick[14],coeff[14], ply[14],kp01[14],kp02[14],kp03[14],kp04[14], kp05[14], kp06[14],kp07[14],kp08[14],kp09[14],kp10[14],num_webs[14]))
#    f.write('16       %.3f     %.3f     %.3f     %.3f       %.3f      %.3f        %.3f       %.3f         %.3f           %.3f         %.2f      %.2f       %.2f       %.2f       %.2f       %.2f      %.2f       %.2f       %.2f        %.2f       %.2f           %d \n'%(radius[15],x_aer[15],z_aer[15],x_sweep[15],z_prebent[15],chord[15],twist[15],pit_ax[15] ,r_thick[15],coeff[15], ply[15],kp01[15],kp02[15],kp03[15],kp04[15], kp05[15], kp06[15],kp07[15],kp08[15],kp09[15],kp10[15],num_webs[15]))
#    f.write('17       %.3f     %.3f     %.3f     %.3f       %.3f      %.3f        %.3f       %.3f         %.3f           %.3f         %.2f      %.2f       %.2f       %.2f       %.2f       %.2f      %.2f       %.2f       %.2f        %.2f       %.2f           %d \n'%(radius[16],x_aer[16],z_aer[16],x_sweep[16],z_prebent[16],chord[16],twist[16],pit_ax[16] ,r_thick[16],coeff[16], ply[16],kp01[16],kp02[16],kp03[16],kp04[16], kp05[16], kp06[16],kp07[16],kp08[16],kp09[16],kp10[16],num_webs[16]))
#    f.write('18       %.3f     %.3f     %.3f     %.3f       %.3f      %.3f        %.3f       %.3f         %.3f           %.3f         %.2f      %.2f       %.2f       %.2f       %.2f       %.2f      %.2f       %.2f       %.2f        %.2f       %.2f           %d \n'%(radius[17],x_aer[17],z_aer[17],x_sweep[17],z_prebent[17],chord[17],twist[17],pit_ax[17] ,r_thick[17],coeff[17], ply[17],kp01[17],kp02[17],kp03[17],kp04[17], kp05[17], kp06[17],kp07[17],kp08[17],kp09[17],kp10[17],num_webs[17]))
#    f.write('19       %.3f     %.3f     %.3f     %.3f       %.3f      %.3f        %.3f       %.3f         %.3f           %.3f         %.2f      %.2f       %.2f       %.2f       %.2f       %.2f      %.2f       %.2f       %.2f        %.2f       %.2f           %d \n'%(radius[18],x_aer[18],z_aer[18],x_sweep[18],z_prebent[18],chord[18],twist[18],pit_ax[18] ,r_thick[18],coeff[18], ply[18],kp01[18],kp02[18],kp03[18],kp04[18], kp05[18], kp06[18],kp07[18],kp08[18],kp09[18],kp10[18],num_webs[18]))
#    f.write('20       %.3f     %.3f     %.3f     %.3f       %.3f      %.3f        %.3f       %.3f         %.3f           %.3f         %.2f      %.2f       %.2f       %.2f       %.2f       %.2f      %.2f       %.2f       %.2f        %.2f       %.2f           %d \n'%(radius[19],x_aer[19],z_aer[19],x_sweep[19],z_prebent[19],chord[19],twist[19],pit_ax[19] ,r_thick[19],coeff[19], ply[19],kp01[19],kp02[19],kp03[19],kp04[19], kp05[19], kp06[19],kp07[19],kp08[19],kp09[19],kp10[19],num_webs[19]))
#    f.write('21       %.3f     %.3f     %.3f     %.3f       %.3f      %.3f        %.3f       %.3f         %.3f           %.3f         %.2f      %.2f       %.2f       %.2f       %.2f       %.2f      %.2f       %.2f       %.2f        %.2f       %.2f           %d \n'%(radius[20],x_aer[20],z_aer[20],x_sweep[20],z_prebent[20],chord[20],twist[20],pit_ax[20] ,r_thick[20],coeff[20], ply[20],kp01[20],kp02[20],kp03[20],kp04[20], kp05[20], kp06[20],kp07[20],kp08[20],kp09[20],kp10[20],num_webs[20]))
#    f.write('#################################################\n')  
#    f.write('Resin&Hardener[$/kgr]: 3.63\n')
#    f.write('Adhesive[$/kgr]: 9.00\n')
#    f.write('Paint[$/kgr]: 7.23\n')
#    f.write('T-bolts&Nuts[$/#]: 37.0\n')
#    f.write('Lightning[$/m]: 40.0\n')
#    f.write('C_Rtip[$/m]: 2.12090 0.08413 1.27778 63.48872 2.00444\n')
#    f.write('C_Amolds[$/m2]: 2.28235 0.21783 0.42451 0.53200\n')
#    f.write('C_Aout[$/m2]: 0.06713 0.06814 0.02200 0.16267 0.04278\n')
#    f.write('Labor[-]: 3.89790 2.5225\n')
#
#with open('../POWER.dir/data.dir/tower_construction.txt','w') as f:
#    f.write('R_tip_tower[m]: %.3f\n'%(max(Height1_new)))
#    f.write('R_hup_tower[-]: 0.000 \n')
#    f.write('#################################################\n')
#    f.write('                     S355      	  [-]         [-]       [-]\n')
#    f.write('cost[$/Kgr]        1413.140         0.000       0.000     0.000\n')
#    f.write('density[m3/Kgr]    8500.000         0.000       0.000     0.000\n')
#    f.write('E11[N/m2]          210.00E+9        0.000       0.000     0.000\n')
#    f.write('E22[N/m2]          210.00E+9        0.000       0.000     0.000\n')
#    f.write('v12[-]             0.300            0.000       0.000     0.000\n')
#    f.write('G12[Nm/m2]         80.77E+9         0.000       0.000     0.000\n')
#    f.write('G13[Nm/m2]         80.77E+9         0.000       0.000     0.000\n')
#    f.write('G23[Nm/m2]         80.77E+9         0.000       0.000     0.000\n')
#    f.write('xc[N/m2]           355.00E+6        0.000       0.000     0.000\n')
#    f.write('xt[N/m2]           355.00E+6        0.000       0.000     0.000\n')
#    f.write('yc[N/m2]           355.00E+6        0.000       0.000     0.000\n')
#    f.write('yt[N/m2]           355.00E+6        0.000       0.000     0.000\n')
#    f.write('st[N/m2]           355.00E+6        0.000       0.000     0.000\n')
#    f.write('#################################################\n')
#    f.write('num radius_tower[-]   c/t[-]      airfoil\n')
#    f.write(' 1   0.00000   1.00000    cylinder.txt\n')
#    f.write(' 2   0.16650   1.00000    cylinder.txt\n')
#    f.write(' 3   0.20939   1.00000    cylinder.txt\n')
#    f.write(' 4   0.29369   1.00000    cylinder.txt\n')
#    f.write(' 5   0.38831   1.00000    cylinder.txt\n')
#    f.write(' 6   0.70228   1.00000    cylinder.txt\n')
#    f.write(' 7   1.00001   1.00000    cylinder.txt\n')
#    f.write('#################################################\n')
#    f.write('num   range\n')
#    f.write(' 1    steel.txt\n')
#    f.write(' 2    steel.txt\n')
#    f.write(' 3    steel.txt\n')
#    f.write(' 4    steel.txt\n')
#    f.write(' 5    steel.txt\n')
#    f.write(' 6    steel.txt\n')
#    f.write(' 7    steel.txt\n')
#    f.write(' 8    steel.txt\n')
#    f.write(' 9    steel.txt\n')
#    f.write('10    steel.txt\n')
#    f.write('11    steel.txt\n')
#    f.write('12    steel.txt\n')
#    f.write('#################################################\n')
#    f.write('num   radius_tower[-] x-aer_tower[-] y-aer_tower[-] x-sweep_tower[m] z-prebent_tower[m] chord_tower[m] twist_tower[deg]  pit-ax[-]  R.thickness[-]  coefficient[-] ply_tower[deg]  kp_tower-01[-]   kp_tower-02[-]   kp_tower-03[-]   kp_tower-04[-]   kp_tower-05[-]   kp_tower-06[-]   kp_tower-07[-]   kp_tower-08[-]   kp_tower-09[-]   kp_tower-10[-]   num_webs_tower[-]\n')
#    f.write(' 1       %.3f             0.000          0.000           0.000              0.000           %.3f           0.000           0.500        1.000           %.3f           0.000           0.87              0.81             0.73             0.53             0.07             0.08             0.29             0.49             0.82             0.88               0  \n'%(Height1_new[0]/max(Height1_new),Diameter1[0],coeff_tower[0]))
#    f.write(' 2       %.3f             0.000          0.000           0.000              0.000           %.3f           0.000           0.500        1.000           %.3f           0.000           0.87              0.81             0.72             0.52             0.07             0.08             0.29             0.49             0.82             0.88               0  \n'%(Height1_new[1]/max(Height1_new),Diameter1[1],coeff_tower[1]))
#    f.write(' 3       %.3f             0.000          0.000           0.000              0.000           %.3f           0.000           0.500        1.000           %.3f           0.000           0.87              0.81             0.68             0.51             0.07             0.08             0.29             0.48             0.82             0.88               0  \n'%(Height1_new[2]/max(Height1_new),Diameter1[2],coeff_tower[2]))
#    f.write(' 4       %.3f             0.000          0.000           0.000              0.000           %.3f           0.000           0.500        1.000           %.3f           0.000           0.87              0.81             0.60             0.44             0.07             0.08             0.28             0.46             0.82             0.88               0  \n'%(Height1_new[3]/max(Height1_new),Diameter1[3],coeff_tower[3]))
#    f.write(' 5       %.3f             0.000          0.000           0.000              0.000           %.3f           0.000           0.500        1.000           %.3f           0.000           0.87              0.81             0.53             0.37             0.07             0.08             0.27             0.44             0.82             0.88               0  \n'%(Height1_new[4]/max(Height1_new),Diameter1[4],coeff_tower[4]))
#    f.write(' 6       %.3f             0.000          0.000           0.000              0.000           %.3f           0.000           0.500        1.000           %.3f           0.000           0.87              0.81             0.47             0.32             0.07             0.08             0.25             0.41             0.82             0.88               0  \n'%(Height1_new[5]/max(Height1_new),Diameter1[5],coeff_tower[5]))
#    f.write(' 7       %.3f             0.000          0.000           0.000              0.000           %.3f           0.000           0.500        1.000           %.3f           0.000           0.87              0.81             0.40             0.26             0.07             0.08             0.22             0.37             0.82             0.88               0  \n'%(Height1_new[6]/max(Height1_new),Diameter1[6],coeff_tower[6]))
#    f.write(' 8       %.3f             0.000          0.000           0.000              0.000           %.3f           0.000           0.500        1.000           %.3f           0.000           0.87              0.81             0.39             0.25             0.07             0.08             0.22             0.36             0.82             0.88               0  \n'%(Height1_new[7]/max(Height1_new),Diameter1[7],coeff_tower[7]))
#    f.write(' 9       %.3f             0.000          0.000           0.000              0.000           %.3f           0.000           0.500        1.000           %.3f           0.000           0.87              0.81             0.37             0.23             0.07             0.08             0.21             0.35             0.82             0.88               0  \n'%(Height1_new[8]/max(Height1_new),Diameter1[8],coeff_tower[8]))
#    f.write('10       %.3f             0.000          0.000           0.000              0.000           %.3f           0.000           0.500        1.000           %.3f           0.000           0.87              0.81             0.37             0.23             0.07             0.08             0.22             0.36             0.82             0.88               0  \n'%(Height1_new[9]/max(Height1_new),Diameter1[9],coeff_tower[9]))
#    f.write('11       %.3f             0.000          0.000           0.000              0.000           %.3f           0.000           0.500        1.000           %.3f           0.000           0.87              0.81             0.37             0.23             0.07             0.08             0.22             0.36             0.82             0.88               0  \n'%(Height1_new[10]/max(Height1_new),Diameter1[10],coeff_tower[10]))
#    f.write('12       %.3f             0.000          0.000           0.000              0.000           %.3f           0.000           0.500        1.000           %.3f           0.000           0.87              0.81             0.38             0.23             0.07             0.08             0.23             0.37             0.82             0.88               0  \n'%(Height1_new[11]/max(Height1_new),Diameter1[11],coeff_tower[11]))
#    f.write('13       %.3f             0.000          0.000           0.000              0.000           %.3f           0.000           0.500        1.000           %.3f           0.000           0.87              0.81             0.38             0.23             0.07             0.08             0.23             0.38             0.82             0.88               0  \n'%(Height1_new[12]/max(Height1_new),Diameter1[12],coeff_tower[12]))
#    f.write('14       %.3f             0.000          0.000           0.000              0.000           %.3f           0.000           0.500        1.000           %.3f           0.000           0.87              0.81             0.38             0.23             0.07             0.08             0.23             0.38             0.82             0.88               0  \n'%(Height1_new[13]/max(Height1_new),Diameter1[13],coeff_tower[13]))
#    f.write('15       %.3f             0.000          0.000           0.000              0.000           %.3f           0.000           0.500        1.000           %.3f           0.000           0.87              0.81             0.38             0.23             0.07             0.08             0.23             0.38             0.82             0.88               0  \n'%(Height1_new[14]/max(Height1_new),Diameter1[14],coeff_tower[14]))
#    f.write('16       %.3f             0.000          0.000           0.000              0.000           %.3f           0.000           0.500        1.000           %.3f           0.000           0.87              0.81             0.38             0.23             0.07             0.08             0.23             0.38             0.82             0.88               0  \n'%(Height1_new[15]/max(Height1_new),Diameter1[15],coeff_tower[15]))
#    f.write('17       %.3f             0.000          0.000           0.000              0.000           %.3f           0.000           0.500        1.000           %.3f           0.000           0.87              0.81             0.38             0.23             0.07             0.08             0.23             0.38             0.82             0.88               0  \n'%(Height1_new[16]/max(Height1_new),Diameter1[16],coeff_tower[16]))
#    f.write('18       %.3f             0.000          0.000           0.000              0.000           %.3f           0.000           0.500        1.000           %.3f           0.000           0.87              0.81             0.38             0.23             0.07             0.08             0.23             0.38             0.82             0.88               0  \n'%(Height1_new[17]/max(Height1_new),Diameter1[17],coeff_tower[17]))
#    f.write('19       %.3f             0.000          0.000           0.000              0.000           %.3f           0.000           0.500        1.000           %.3f           0.000           0.87              0.81             0.38             0.23             0.07             0.08             0.23             0.38             0.82             0.88               0  \n'%(Height1_new[18]/max(Height1_new),Diameter1[18],coeff_tower[18]))
#    f.write('20       %.3f             0.000          0.000           0.000              0.000           %.3f           0.000           0.500        1.000           %.3f           0.000           0.87              0.81             0.38             0.23             0.07             0.08             0.23             0.38             0.82             0.88               0  \n'%(Height1_new[19]/max(Height1_new),Diameter1[19],coeff_tower[19]))
#    f.write('21       %.3f             0.000          0.000           0.000              0.000           %.3f           0.000           0.500        1.000           %.3f           0.000           0.87              0.81             0.38             0.23             0.07             0.08             0.23             0.38             0.82             0.88               0  \n'%(Height1_new[20]/max(Height1_new),Diameter1[20],coeff_tower[20]))
#    f.write('#################################################\n')  
#    f.write('Resin&Hardener[$/kgr]: 3.63\n')
#    f.write('Adhesive[$/kgr]: 9.00\n')
#    f.write('Paint[$/kgr]: 7.23\n')
#    f.write('T-bolts&Nuts[$/#]: 37.0\n')
#    f.write('Lightning[$/m]: 40.0\n')
#    f.write('C_Rtip[$/m]: 2.12090 0.08413 1.27778 63.48872 2.00444\n')
#    f.write('C_Amolds[$/m2]: 2.28235 0.21783 0.42451 0.53200\n')
#    f.write('C_Aout[$/m2]: 0.06713 0.06814 0.02200 0.16267 0.04278\n')
#    f.write('Labor[-]: 3.89790 2.5225\n')
