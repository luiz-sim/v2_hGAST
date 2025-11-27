#import matplotlib.pyplot as plt
import os, math
import numpy as np
from scipy.interpolate import interp1d
from numpy.linalg import inv
from scipy.optimize import fsolve

# Angle between vectors
#################################################
angle = lambda v1,v2: np.arccos(np.dot(v1/np.linalg.norm(v1),v2/np.linalg.norm(v2))) if np.cross(v1/np.linalg.norm(v1),v2/np.linalg.norm(v2))>=0.0 else 2.0*math.pi-np.arccos(np.dot(v1/np.linalg.norm(v1),v2/np.linalg.norm(v2)))

# Bezier- Curve
#################################################
def bezier_curve(control, r):
    x_mat=np.array([i[0] for i in control]) ; y_mat=np.array([i[1] for i in control])
    c_mat=np.zeros([len(control),len(control)],dtype=float)
    for i in range(0,len(c_mat)):
       for j in range(0,len(c_mat[0])):
           if i<=j: c_mat[i][j]=((-1)**(j-i-2))*(math.factorial(len(control)-1)/(math.factorial(j)*math.factorial(len(control)-j-1)))*(math.factorial(j)/(math.factorial(i)*math.factorial(j-i)))
    x_100=[] ; y_100=[]
    for i in range(0,101):
        t=(i-1.0)/99.0 ; t_mat=np.array([t**i for i in range(0,len(control))])
        x_100.append(np.matmul(np.matmul(x_mat.transpose(),c_mat),t_mat)) ; y_100.append(np.matmul(np.matmul(y_mat.transpose(),c_mat),t_mat))

    y_bez = []
    for i in range(0,len(r)):
        for j in range(1,len(x_100)):
            if r[i]>=x_100[j-1] and r[i]<x_100[j]+0.0001: y_bez.append(y_100[j-1]+(y_100[j]-y_100[j-1])*(r[i]-x_100[j-1])/(x_100[j]-x_100[j-1]))
    return [(r[i],y_bez[i]) for i in range(0,len(r))]

# Airfoil-geometry
#################################################
def airfoil_geometry(profil, radius, chord, twist, pit_ax, kp01, kp02, kp03, kp04, kp05, kp06, kp07, kp08, kp09, kp10):
    airfoils=[]
    rad_pro=[float(profil[i][0]) for i in range(0,len(profil))]
    for i in range(1,len(rad_pro)):
        if radius>=rad_pro[i-1] and radius<rad_pro[i]:
            x_air=list(np.arange(0.0, 1.01, 0.01))
            con_lp1, con_hp1 = [[float(profil[i-1][1][j]), float(profil[i-1][2][j])] for j in range(0,len(profil[i-1][1])) if profil[i-1][1][j]!='none'], [[float(profil[i-1][3][j]), float(profil[i-1][4][j])] for j in range(0,len(profil[i-1][3]))]
            y_lp1, y_hp1 = [float(bezier_curve(con_lp1, x_air)[j][1]) for j in range(0,len(x_air))], [float(bezier_curve(con_hp1, x_air)[j][1]) for j in range(0,len(x_air))]

            con_lp2, con_hp2 = [[float(profil[i][1][j]), float(profil[i][2][j])] for j in range(0,len(profil[i][1])) if profil[i][1][j]!='none'], [[float(profil[i][3][j]), float(profil[i][4][j])] for j in range(0,len(profil[i][3]))]
            y_lp2, y_hp2 = [float(bezier_curve(con_lp2, x_air)[j][1]) for j in range(0,len(x_air))], [float(bezier_curve(con_hp2, x_air)[j][1]) for j in range(0,len(x_air))]

            y_lp, y_hp = [y_lp1[j]+(y_lp2[j]-y_lp1[j])*(radius-rad_pro[i-1])/(rad_pro[i]-rad_pro[i-1]) for j in range(0,len(y_lp1))], [y_hp1[j]+(y_hp2[j]-y_hp1[j])*(radius-rad_pro[i-1])/(rad_pro[i]-rad_pro[i-1]) for j in range(0,len(y_hp1))]

    for j in range(0,len(x_air)):
        airfoils.append((chord*((x_air[len(x_air)-1-j]-pit_ax)*np.cos(twist*math.pi/180.0)-y_hp[len(x_air)-1-j]*np.sin(twist*math.pi/180.0)), chord*((x_air[len(x_air)-1-j]-pit_ax)*np.sin(twist*math.pi/180.0)+y_hp[len(x_air)-1-j]*np.cos(twist*math.pi/180.0))))
    for j in range(1,len(x_air)):
        airfoils.append((chord*((x_air[j]-pit_ax)*np.cos(twist*math.pi/180.0)-y_lp[j]*np.sin(twist*math.pi/180.0)), chord*((x_air[j]-pit_ax)*np.sin(twist*math.pi/180.0)+y_lp[j]*np.cos(twist*math.pi/180.0))))

    num=0 ; kp=[]
    for j in range(len(x_air)-1,0,-1):
        if abs(x_air[j]-kp01)<0.001 or abs(x_air[j]-kp02)<0.001 or abs(x_air[j]-kp03)<0.001 or abs(x_air[j]-kp04)<0.001 or abs(x_air[j]-kp05)<0.001: kp.append(num)
        num=num+1
    for j in range(1,len(x_air)):
        if abs(x_air[j]-kp06)<0.001 or abs(x_air[j]-kp07)<0.001 or abs(x_air[j]-kp08)<0.001 or abs(x_air[j]-kp09)<0.001 or abs(x_air[j]-kp10)<0.001: kp.append(num)
        num=num+1
    return airfoils, kp

# Grid
#################################################
def grid(air, kp, laminates, num_webs, ply):
    x_geo, y_geo = [air[i][0] for i in range(0,len(air))], [air[i][1] for i in range(0,len(air))]

    # SKIN
    kp.insert(0,0) ; kp.append(len(air))
    k_num=[kp[i]-kp[i-1] for i in range(1,len(kp))]
    X_Grid=[] ; Y_Grid=[] ; C_Grid=[] ; Thick_lay=[] ; Length_lay=[] ; Ply_lay=[] ; Thick_ele=[] ; Length_ele=[] ; Theta_ele=[] ; Center_ele=[]
    num=-1
    for i in range(0,len(k_num)):
        X_Grid.append([]) ; Y_Grid.append([]) ; C_Grid.append([]) ; Thick_lay.append([]) ; Length_lay.append([]) ; Ply_lay.append([]) ; Thick_ele.append([]) ; Length_ele.append([]) ; Theta_ele.append([]) ; Center_ele.append([])
        for j in range(0,k_num[i]):
            num=num+1
            X_Grid[i].append([]) ; Y_Grid[i].append([]) ; C_Grid[i].append([]) ; Thick_lay[i].append([]) ; Length_lay[i].append([]) ; Ply_lay[i].append([])
            if num==0:
                n1=len(air)-1 ; n2=num ; n3=num+1 ; n4=num+2
            elif num==len(air)-2:
                n1=num-1 ; n2=num ; n3=len(air)-1 ; n4=0
            elif num==len(air)-1:
                n1=num-1 ; n2=num ; n3=0 ; n4=1
            else:
                n1=num-1 ; n2=num ; n3=num+1 ; n4=num+2

            theta=angle([1.0, 0.0], [x_geo[n3]-x_geo[n2],y_geo[n3]-y_geo[n2]])
            xloc=[(x_geo[n1]-x_geo[n2])*np.cos(-theta)-(y_geo[n1]-y_geo[n2])*np.sin(-theta), (x_geo[n2]-x_geo[n2])*np.cos(-theta)-(y_geo[n2]-y_geo[n2])*np.sin(-theta), (x_geo[n3]-x_geo[n2])*np.cos(-theta)-(y_geo[n3]-y_geo[n2])*np.sin(-theta), (x_geo[n4]-x_geo[n2])*np.cos(-theta)-(y_geo[n4]-y_geo[n2])*np.sin(-theta)]
            yloc=[(x_geo[n1]-x_geo[n2])*np.sin(-theta)+(y_geo[n1]-y_geo[n2])*np.cos(-theta), (x_geo[n2]-x_geo[n2])*np.sin(-theta)+(y_geo[n2]-y_geo[n2])*np.cos(-theta), (x_geo[n3]-x_geo[n2])*np.sin(-theta)+(y_geo[n3]-y_geo[n2])*np.cos(-theta), (x_geo[n4]-x_geo[n2])*np.sin(-theta)+(y_geo[n4]-y_geo[n2])*np.cos(-theta)]

            delta=0.0 ; xlocele=[xloc[1], xloc[2]] ; ylocele=[yloc[1], yloc[2]]
            for k in range(0,len(laminates[i])):
                if laminates[i][k]>0.0:
                    delta=laminates[i][k]+delta
                    xlocele.append(xloc[2]+delta*np.tan(-(math.pi-angle([xloc[1]-xloc[2], yloc[1]-yloc[2]], [xloc[3]-xloc[2], yloc[3]-yloc[2]]))/2.0)) ; ylocele.append(-delta)
                    xlocele.append(xloc[1]+delta*np.tan((math.pi-angle([xloc[0]-xloc[1], yloc[0]-yloc[1]], [xloc[2]-xloc[1], yloc[2]-yloc[1]]))/2.0))  ; ylocele.append(-delta)
                    Length_lay[i][j].append(float((abs(xlocele[1]-xlocele[0])+abs(xlocele[3]-xlocele[2]))/2.0))
                    Thick_lay[i][j].append(float(laminates[i][k]))
                    if i==3 and (k==2 or k==4):
                        Ply_lay[i][j].append(-float(ply))
                    elif i==7 and (k==2 or k==4):
                        Ply_lay[i][j].append(float(ply))
                    else:
                        Ply_lay[i][j].append(float(0.0))
                    X_Grid[i][j].append([xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+x_geo[n2], xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+x_geo[n2], xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+x_geo[n2], xlocele[3]*np.cos(theta)-ylocele[3]*np.sin(theta)+x_geo[n2], xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+x_geo[n2], None])
                    Y_Grid[i][j].append([xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+y_geo[n2], xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+y_geo[n2], xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+y_geo[n2], xlocele[3]*np.sin(theta)+ylocele[3]*np.cos(theta)+y_geo[n2], xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+y_geo[n2], None])
                    if k==0 or k==6:
                        C_Grid[i][j].append('green')
                    elif k==1 or k==5:
                        C_Grid[i][j].append('yellow')
                    elif k==2 or k==4:
                        C_Grid[i][j].append('red')
                    else:
                        C_Grid[i][j].append('gray')
                    xlocele=[xlocele[3], xlocele[2]] ; ylocele=[ylocele[3], ylocele[2]]
                else:
                    Length_lay[i][j].append(0.0) ; Thick_lay[i][j].append(0.0) ; Ply_lay[i][j].append(0.0)

            Thick_ele[i].append(sum(laminates[i][:]))
            Length_ele[i].append((math.sqrt((X_Grid[i][j][0][1]-X_Grid[i][j][0][0])**2.0+(Y_Grid[i][j][0][1]-Y_Grid[i][j][0][0])**2.0)+math.sqrt((X_Grid[i][j][-1][3]-X_Grid[i][j][-1][2])**2.0+(Y_Grid[i][j][-1][3]-Y_Grid[i][j][-1][2])**2.0))/2.0)
            Theta_ele[i].append(((x_geo[n3]-x_geo[n2])/math.sqrt((x_geo[n3]-x_geo[n2])**2.0+(y_geo[n3]-y_geo[n2])**2.0), (y_geo[n3]-y_geo[n2])/math.sqrt((x_geo[n3]-x_geo[n2])**2.0+(y_geo[n3]-y_geo[n2])**2.0)))
            Center_ele[i].append(((X_Grid[i][j][0][0]+X_Grid[i][j][0][1]+X_Grid[i][j][-1][2]+X_Grid[i][j][-1][3])/4.0, (Y_Grid[i][j][0][0]+Y_Grid[i][j][0][1]+Y_Grid[i][j][-1][2]+Y_Grid[i][j][-1][3])/4.0))
    # WEB-01
    if num_webs==1 or num_webs==2 or num_webs==3:
        X_Grid.append([])                ; Y_Grid.append([])                ; C_Grid.append([])                ; Thick_lay.append([])                ; Length_lay.append([])                ; Ply_lay.append([])             ; Thick_ele.append([])         ; Length_ele.append([])         ; Theta_ele.append([])         ; Center_ele.append([])
        X_Grid[11].extend([[],[],[],[]]) ; Y_Grid[11].extend([[],[],[],[]]) ; C_Grid[11].extend([[],[],[],[]]) ; Thick_lay[11].extend([[],[],[],[]]) ; Length_lay[11].extend([[],[],[],[]]) ; Ply_lay[11].extend([[],[],[],[]])

        Xa=X_Grid[7][0][-1][3] ; Xb=X_Grid[7][0][-1][2] ; Xc=X_Grid[3][-1][-1][3] ; Xd=X_Grid[3][-1][-1][2]
        Ya=Y_Grid[7][0][-1][3] ; Yb=Y_Grid[7][0][-1][2] ; Yc=Y_Grid[3][-1][-1][3] ; Yd=Y_Grid[3][-1][-1][2]

        theta=angle([1.0, 0.0], [Xc-Xb,Yc-Yb])
        xloc=[(Xa-Xb)*np.cos(-theta)-(Ya-Yb)*np.sin(-theta), (Xb-Xb)*np.cos(-theta)-(Yb-Yb)*np.sin(-theta), (Xc-Xb)*np.cos(-theta)-(Yc-Yb)*np.sin(-theta), (Xd-Xb)*np.cos(-theta)-(Yd-Yb)*np.sin(-theta)]
        yloc=[(Xa-Xb)*np.sin(-theta)+(Ya-Yb)*np.cos(-theta), (Xb-Xb)*np.sin(-theta)+(Yb-Yb)*np.cos(-theta), (Xc-Xb)*np.sin(-theta)+(Yc-Yb)*np.cos(-theta), (Xd-Xb)*np.sin(-theta)+(Yd-Yb)*np.cos(-theta)]
        delta=0.0 ; xlocele=[xloc[1], xloc[2]] ; ylocele=[yloc[1], yloc[2]]
        for k in range(0,len(laminates[11])):
            if laminates[11][k]>0.0:
                delta=laminates[11][k]+delta
                xlocele.append(-delta*(xloc[1]-xloc[0])/(yloc[1]-yloc[0]))        ; ylocele.append(-delta)
                xlocele.append(xloc[2]-delta*(xloc[3]-xloc[2])/(yloc[3]-yloc[2])) ; ylocele.append(-delta)

                Length_lay[11][0].append(float((abs(xlocele[1]-xlocele[0])+abs(xlocele[3]-xlocele[2]))/8.0)) ; Length_lay[11][1].append(float((abs(xlocele[1]-xlocele[0])+abs(xlocele[3]-xlocele[2]))/8.0)) ; Length_lay[11][2].append(float((abs(xlocele[1]-xlocele[0])+abs(xlocele[3]-xlocele[2]))/8.0)) ; Length_lay[11][3].append(float((abs(xlocele[1]-xlocele[0])+abs(xlocele[3]-xlocele[2]))/8.0))
                Thick_lay[11][0].append(float(laminates[11][k]))                                             ; Thick_lay[11][1].append(float(laminates[11][k]))                                             ; Thick_lay[11][2].append(float(laminates[11][k]))                                             ; Thick_lay[11][3].append(float(laminates[11][k]))
                Ply_lay[11][0].append(float(0.0))                                                            ; Ply_lay[11][1].append(float(0.0))                                                            ; Ply_lay[11][2].append(float(0.0))                                                            ; Ply_lay[11][3].append(float(0.0))

                X_Grid[11][0].append([xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+0.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+1.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb+1.0*((xlocele[3]*np.cos(theta)-ylocele[3]*np.sin(theta)+Xb)-(xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb))/4.0, xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb+0.0*((xlocele[3]*np.cos(theta)-ylocele[3]*np.sin(theta)+Xb)-(xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb))/4.0, xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+0.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, None])
                Y_Grid[11][0].append([xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+0.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+1.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb+1.0*((xlocele[3]*np.sin(theta)+ylocele[3]*np.cos(theta)+Yb)-(xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb))/4.0, xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb+0.0*((xlocele[3]*np.sin(theta)+ylocele[3]*np.cos(theta)+Yb)-(xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb))/4.0, xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+0.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, None])

                X_Grid[11][1].append([xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+1.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+2.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb+2.0*((xlocele[3]*np.cos(theta)-ylocele[3]*np.sin(theta)+Xb)-(xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb))/4.0, xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb+1.0*((xlocele[3]*np.cos(theta)-ylocele[3]*np.sin(theta)+Xb)-(xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb))/4.0, xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+1.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, None])
                Y_Grid[11][1].append([xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+1.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+2.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb+2.0*((xlocele[3]*np.sin(theta)+ylocele[3]*np.cos(theta)+Yb)-(xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb))/4.0, xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb+1.0*((xlocele[3]*np.sin(theta)+ylocele[3]*np.cos(theta)+Yb)-(xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb))/4.0, xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+1.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, None])

                X_Grid[11][2].append([xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+2.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+3.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb+3.0*((xlocele[3]*np.cos(theta)-ylocele[3]*np.sin(theta)+Xb)-(xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb))/4.0, xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb+2.0*((xlocele[3]*np.cos(theta)-ylocele[3]*np.sin(theta)+Xb)-(xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb))/4.0, xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+2.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, None])
                Y_Grid[11][2].append([xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+2.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+3.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb+3.0*((xlocele[3]*np.sin(theta)+ylocele[3]*np.cos(theta)+Yb)-(xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb))/4.0, xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb+2.0*((xlocele[3]*np.sin(theta)+ylocele[3]*np.cos(theta)+Yb)-(xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb))/4.0, xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+2.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, None])

                X_Grid[11][3].append([xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+3.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+4.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb+4.0*((xlocele[3]*np.cos(theta)-ylocele[3]*np.sin(theta)+Xb)-(xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb))/4.0, xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb+3.0*((xlocele[3]*np.cos(theta)-ylocele[3]*np.sin(theta)+Xb)-(xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb))/4.0, xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+3.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, None])
                Y_Grid[11][3].append([xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+3.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+4.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb+4.0*((xlocele[3]*np.sin(theta)+ylocele[3]*np.cos(theta)+Yb)-(xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb))/4.0, xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb+3.0*((xlocele[3]*np.sin(theta)+ylocele[3]*np.cos(theta)+Yb)-(xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb))/4.0, xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+3.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, None])

                if k==0 or k==6:
                    C_Grid[11][0].append('green')  ; C_Grid[11][1].append('green')  ; C_Grid[11][2].append('green')  ; C_Grid[11][3].append('green')
                elif k==1 or k==5:
                    C_Grid[11][0].append('yellow') ; C_Grid[11][1].append('yellow') ; C_Grid[11][2].append('yellow') ; C_Grid[11][3].append('yellow')
                elif k==2 or k==4:
                    C_Grid[11][0].append('red')    ; C_Grid[11][1].append('red')    ; C_Grid[11][2].append('red')    ; C_Grid[11][3].append('red')
                else:
                    C_Grid[11][0].append('gray')   ; C_Grid[11][1].append('gray')   ; C_Grid[11][2].append('gray')   ; C_Grid[11][3].append('gray')
                xlocele=[xlocele[2], xlocele[3]] ; ylocele=[ylocele[2], ylocele[3]]
            else:
                Length_lay[11][0].append(0.0) ; Thick_lay[11][0].append(0.0) ; Ply_lay[11][0].append(0.0)
                Length_lay[11][1].append(0.0) ; Thick_lay[11][1].append(0.0) ; Ply_lay[11][1].append(0.0)
                Length_lay[11][2].append(0.0) ; Thick_lay[11][2].append(0.0) ; Ply_lay[11][2].append(0.0)
                Length_lay[11][3].append(0.0) ; Thick_lay[11][3].append(0.0) ; Ply_lay[11][3].append(0.0)

        Thick_ele[11].append(sum(laminates[11][:]))
        Thick_ele[11].append(sum(laminates[11][:]))
        Thick_ele[11].append(sum(laminates[11][:]))
        Thick_ele[11].append(sum(laminates[11][:]))

        Length_ele[11].append((math.sqrt((X_Grid[11][0][0][1]-X_Grid[11][0][0][0])**2.0+(Y_Grid[11][0][0][1]-Y_Grid[11][0][0][0])**2.0)+math.sqrt((X_Grid[11][0][-1][3]-X_Grid[11][0][-1][2])**2.0+(Y_Grid[11][0][-1][3]-Y_Grid[11][0][-1][2])**2.0))/2.0)
        Length_ele[11].append((math.sqrt((X_Grid[11][1][0][1]-X_Grid[11][1][0][0])**2.0+(Y_Grid[11][1][0][1]-Y_Grid[11][1][0][0])**2.0)+math.sqrt((X_Grid[11][1][-1][3]-X_Grid[11][1][-1][2])**2.0+(Y_Grid[11][1][-1][3]-Y_Grid[11][1][-1][2])**2.0))/2.0)
        Length_ele[11].append((math.sqrt((X_Grid[11][2][0][1]-X_Grid[11][2][0][0])**2.0+(Y_Grid[11][2][0][1]-Y_Grid[11][2][0][0])**2.0)+math.sqrt((X_Grid[11][2][-1][3]-X_Grid[11][2][-1][2])**2.0+(Y_Grid[11][2][-1][3]-Y_Grid[11][2][-1][2])**2.0))/2.0)
        Length_ele[11].append((math.sqrt((X_Grid[11][3][0][1]-X_Grid[11][3][0][0])**2.0+(Y_Grid[11][3][0][1]-Y_Grid[11][3][0][0])**2.0)+math.sqrt((X_Grid[11][3][-1][3]-X_Grid[11][3][-1][2])**2.0+(Y_Grid[11][3][-1][3]-Y_Grid[11][3][-1][2])**2.0))/2.0)

        Theta_ele[11].append(((X_Grid[11][0][0][1]-X_Grid[11][0][0][0])/math.sqrt((X_Grid[11][0][0][1]-X_Grid[11][0][0][0])**2.0+(Y_Grid[11][0][0][1]-Y_Grid[11][0][0][0])**2.0), (Y_Grid[11][0][0][1]-Y_Grid[11][0][0][0])/math.sqrt((X_Grid[11][0][0][1]-X_Grid[11][0][0][0])**2.0+(Y_Grid[11][0][0][1]-Y_Grid[11][0][0][0])**2.0)))
        Theta_ele[11].append(((X_Grid[11][1][0][1]-X_Grid[11][1][0][0])/math.sqrt((X_Grid[11][1][0][1]-X_Grid[11][1][0][0])**2.0+(Y_Grid[11][1][0][1]-Y_Grid[11][1][0][0])**2.0), (Y_Grid[11][1][0][1]-Y_Grid[11][1][0][0])/math.sqrt((X_Grid[11][1][0][1]-X_Grid[11][1][0][0])**2.0+(Y_Grid[11][1][0][1]-Y_Grid[11][1][0][0])**2.0)))
        Theta_ele[11].append(((X_Grid[11][2][0][1]-X_Grid[11][2][0][0])/math.sqrt((X_Grid[11][2][0][1]-X_Grid[11][2][0][0])**2.0+(Y_Grid[11][2][0][1]-Y_Grid[11][2][0][0])**2.0), (Y_Grid[11][2][0][1]-Y_Grid[11][2][0][0])/math.sqrt((X_Grid[11][2][0][1]-X_Grid[11][2][0][0])**2.0+(Y_Grid[11][2][0][1]-Y_Grid[11][2][0][0])**2.0)))
        Theta_ele[11].append(((X_Grid[11][3][0][1]-X_Grid[11][3][0][0])/math.sqrt((X_Grid[11][3][0][1]-X_Grid[11][3][0][0])**2.0+(Y_Grid[11][3][0][1]-Y_Grid[11][3][0][0])**2.0), (Y_Grid[11][3][0][1]-Y_Grid[11][3][0][0])/math.sqrt((X_Grid[11][3][0][1]-X_Grid[11][3][0][0])**2.0+(Y_Grid[11][3][0][1]-Y_Grid[11][3][0][0])**2.0)))

        Center_ele[11].append(((X_Grid[11][0][0][0]+X_Grid[11][0][0][1]+X_Grid[11][0][-1][2]+X_Grid[11][0][-1][3])/4.0, (Y_Grid[11][0][0][0]+Y_Grid[11][0][0][1]+Y_Grid[11][0][-1][2]+Y_Grid[11][0][-1][3])/4.0))
        Center_ele[11].append(((X_Grid[11][1][0][0]+X_Grid[11][1][0][1]+X_Grid[11][1][-1][2]+X_Grid[11][1][-1][3])/4.0, (Y_Grid[11][1][0][0]+Y_Grid[11][1][0][1]+Y_Grid[11][1][-1][2]+Y_Grid[11][1][-1][3])/4.0))
        Center_ele[11].append(((X_Grid[11][2][0][0]+X_Grid[11][2][0][1]+X_Grid[11][2][-1][2]+X_Grid[11][2][-1][3])/4.0, (Y_Grid[11][2][0][0]+Y_Grid[11][2][0][1]+Y_Grid[11][2][-1][2]+Y_Grid[11][2][-1][3])/4.0))
        Center_ele[11].append(((X_Grid[11][3][0][0]+X_Grid[11][3][0][1]+X_Grid[11][3][-1][2]+X_Grid[11][3][-1][3])/4.0, (Y_Grid[11][3][0][0]+Y_Grid[11][3][0][1]+Y_Grid[11][3][-1][2]+Y_Grid[11][3][-1][3])/4.0))

    # WEB-02
    if num_webs==2 or num_webs==3:
        X_Grid.append([])                ; Y_Grid.append([])                ; C_Grid.append([])                ; Thick_lay.append([])                ; Length_lay.append([])                ; Ply_lay.append([])             ; Thick_ele.append([])         ; Length_ele.append([])         ; Theta_ele.append([])         ; Center_ele.append([])
        X_Grid[12].extend([[],[],[],[]]) ; Y_Grid[12].extend([[],[],[],[]]) ; C_Grid[12].extend([[],[],[],[]]) ; Thick_lay[12].extend([[],[],[],[]]) ; Length_lay[12].extend([[],[],[],[]]) ; Ply_lay[12].extend([[],[],[],[]])

        Xa=X_Grid[3][0][-1][3] ; Xb=X_Grid[3][0][-1][2] ; Xc=X_Grid[7][-1][-1][3] ; Xd=X_Grid[7][-1][-1][2]
        Ya=Y_Grid[3][0][-1][3] ; Yb=Y_Grid[3][0][-1][2] ; Yc=Y_Grid[7][-1][-1][3] ; Yd=Y_Grid[7][-1][-1][2]

        theta=angle([1.0, 0.0], [Xc-Xb,Yc-Yb])
        xloc=[(Xa-Xb)*np.cos(-theta)-(Ya-Yb)*np.sin(-theta), (Xb-Xb)*np.cos(-theta)-(Yb-Yb)*np.sin(-theta), (Xc-Xb)*np.cos(-theta)-(Yc-Yb)*np.sin(-theta), (Xd-Xb)*np.cos(-theta)-(Yd-Yb)*np.sin(-theta)]
        yloc=[(Xa-Xb)*np.sin(-theta)+(Ya-Yb)*np.cos(-theta), (Xb-Xb)*np.sin(-theta)+(Yb-Yb)*np.cos(-theta), (Xc-Xb)*np.sin(-theta)+(Yc-Yb)*np.cos(-theta), (Xd-Xb)*np.sin(-theta)+(Yd-Yb)*np.cos(-theta)]
        delta=0.0 ; xlocele=[xloc[1], xloc[2]] ; ylocele=[yloc[1], yloc[2]]
        for k in range(0,len(laminates[11])):
            if laminates[11][k]>0.0:
                delta=laminates[11][k]+delta
                xlocele.append(-delta*(xloc[1]-xloc[0])/(yloc[1]-yloc[0]))        ; ylocele.append(-delta)
                xlocele.append(xloc[2]-delta*(xloc[3]-xloc[2])/(yloc[3]-yloc[2])) ; ylocele.append(-delta)

                Length_lay[12][0].append(float((abs(xlocele[1]-xlocele[0])+abs(xlocele[3]-xlocele[2]))/8.0)) ; Length_lay[12][1].append(float((abs(xlocele[1]-xlocele[0])+abs(xlocele[3]-xlocele[2]))/8.0)) ; Length_lay[12][2].append(float((abs(xlocele[1]-xlocele[0])+abs(xlocele[3]-xlocele[2]))/8.0)) ; Length_lay[12][3].append(float((abs(xlocele[1]-xlocele[0])+abs(xlocele[3]-xlocele[2]))/8.0))
                Thick_lay[12][0].append(float(laminates[11][k]))                                             ; Thick_lay[12][1].append(float(laminates[11][k]))                                             ; Thick_lay[12][2].append(float(laminates[11][k]))                                             ; Thick_lay[12][3].append(float(laminates[11][k]))
                Ply_lay[12][0].append(float(0.0))                                                            ; Ply_lay[12][1].append(float(0.0))                                                            ; Ply_lay[12][2].append(float(0.0))                                                            ; Ply_lay[12][3].append(float(0.0))

                X_Grid[12][0].append([xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+0.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+1.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb+1.0*((xlocele[3]*np.cos(theta)-ylocele[3]*np.sin(theta)+Xb)-(xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb))/4.0, xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb+0.0*((xlocele[3]*np.cos(theta)-ylocele[3]*np.sin(theta)+Xb)-(xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb))/4.0, xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+0.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, None])
                Y_Grid[12][0].append([xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+0.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+1.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb+1.0*((xlocele[3]*np.sin(theta)+ylocele[3]*np.cos(theta)+Yb)-(xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb))/4.0, xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb+0.0*((xlocele[3]*np.sin(theta)+ylocele[3]*np.cos(theta)+Yb)-(xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb))/4.0, xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+0.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, None])

                X_Grid[12][1].append([xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+1.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+2.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb+2.0*((xlocele[3]*np.cos(theta)-ylocele[3]*np.sin(theta)+Xb)-(xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb))/4.0, xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb+1.0*((xlocele[3]*np.cos(theta)-ylocele[3]*np.sin(theta)+Xb)-(xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb))/4.0, xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+1.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, None])
                Y_Grid[12][1].append([xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+1.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+2.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb+2.0*((xlocele[3]*np.sin(theta)+ylocele[3]*np.cos(theta)+Yb)-(xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb))/4.0, xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb+1.0*((xlocele[3]*np.sin(theta)+ylocele[3]*np.cos(theta)+Yb)-(xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb))/4.0, xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+1.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, None])

                X_Grid[12][2].append([xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+2.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+3.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb+3.0*((xlocele[3]*np.cos(theta)-ylocele[3]*np.sin(theta)+Xb)-(xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb))/4.0, xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb+2.0*((xlocele[3]*np.cos(theta)-ylocele[3]*np.sin(theta)+Xb)-(xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb))/4.0, xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+2.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, None])
                Y_Grid[12][2].append([xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+2.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+3.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb+3.0*((xlocele[3]*np.sin(theta)+ylocele[3]*np.cos(theta)+Yb)-(xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb))/4.0, xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb+2.0*((xlocele[3]*np.sin(theta)+ylocele[3]*np.cos(theta)+Yb)-(xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb))/4.0, xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+2.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, None])

                X_Grid[12][3].append([xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+3.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+4.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb+4.0*((xlocele[3]*np.cos(theta)-ylocele[3]*np.sin(theta)+Xb)-(xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb))/4.0, xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb+3.0*((xlocele[3]*np.cos(theta)-ylocele[3]*np.sin(theta)+Xb)-(xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb))/4.0, xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+3.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, None])
                Y_Grid[12][3].append([xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+3.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+4.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb+4.0*((xlocele[3]*np.sin(theta)+ylocele[3]*np.cos(theta)+Yb)-(xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb))/4.0, xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb+3.0*((xlocele[3]*np.sin(theta)+ylocele[3]*np.cos(theta)+Yb)-(xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb))/4.0, xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+3.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, None])

                if k==0 or k==6:
                    C_Grid[12][0].append('green')  ; C_Grid[12][1].append('green')  ; C_Grid[12][2].append('green')  ; C_Grid[12][3].append('green')
                elif k==1 or k==5:
                    C_Grid[12][0].append('yellow') ; C_Grid[12][1].append('yellow') ; C_Grid[12][2].append('yellow') ; C_Grid[12][3].append('yellow')
                elif k==2 or k==4:
                    C_Grid[12][0].append('red')    ; C_Grid[12][1].append('red')    ; C_Grid[12][2].append('red')    ; C_Grid[12][3].append('red')
                else:
                    C_Grid[12][0].append('gray')   ; C_Grid[12][1].append('gray')   ; C_Grid[12][2].append('gray')   ; C_Grid[12][3].append('gray')
                xlocele=[xlocele[2], xlocele[3]] ; ylocele=[ylocele[2], ylocele[3]]
            else:
                Length_lay[12][0].append(0.0) ; Thick_lay[12][0].append(0.0) ; Ply_lay[12][0].append(0.0)
                Length_lay[12][1].append(0.0) ; Thick_lay[12][1].append(0.0) ; Ply_lay[12][1].append(0.0)
                Length_lay[12][2].append(0.0) ; Thick_lay[12][2].append(0.0) ; Ply_lay[12][2].append(0.0)
                Length_lay[12][3].append(0.0) ; Thick_lay[12][3].append(0.0) ; Ply_lay[12][3].append(0.0)

        Thick_ele[12].append(sum(laminates[11][:]))
        Thick_ele[12].append(sum(laminates[11][:]))
        Thick_ele[12].append(sum(laminates[11][:]))
        Thick_ele[12].append(sum(laminates[11][:]))

        Length_ele[12].append((math.sqrt((X_Grid[12][0][0][1]-X_Grid[12][0][0][0])**2.0+(Y_Grid[12][0][0][1]-Y_Grid[12][0][0][0])**2.0)+math.sqrt((X_Grid[12][0][-1][3]-X_Grid[12][0][-1][2])**2.0+(Y_Grid[12][0][-1][3]-Y_Grid[12][0][-1][2])**2.0))/2.0)
        Length_ele[12].append((math.sqrt((X_Grid[12][1][0][1]-X_Grid[12][1][0][0])**2.0+(Y_Grid[12][1][0][1]-Y_Grid[12][1][0][0])**2.0)+math.sqrt((X_Grid[12][1][-1][3]-X_Grid[12][1][-1][2])**2.0+(Y_Grid[12][1][-1][3]-Y_Grid[12][1][-1][2])**2.0))/2.0)
        Length_ele[12].append((math.sqrt((X_Grid[12][2][0][1]-X_Grid[12][2][0][0])**2.0+(Y_Grid[12][2][0][1]-Y_Grid[12][2][0][0])**2.0)+math.sqrt((X_Grid[12][2][-1][3]-X_Grid[12][2][-1][2])**2.0+(Y_Grid[12][2][-1][3]-Y_Grid[12][2][-1][2])**2.0))/2.0)
        Length_ele[12].append((math.sqrt((X_Grid[12][3][0][1]-X_Grid[12][3][0][0])**2.0+(Y_Grid[12][3][0][1]-Y_Grid[12][3][0][0])**2.0)+math.sqrt((X_Grid[12][3][-1][3]-X_Grid[12][3][-1][2])**2.0+(Y_Grid[12][3][-1][3]-Y_Grid[12][3][-1][2])**2.0))/2.0)

        Theta_ele[12].append(((X_Grid[12][0][0][1]-X_Grid[12][0][0][0])/math.sqrt((X_Grid[12][0][0][1]-X_Grid[12][0][0][0])**2.0+(Y_Grid[12][0][0][1]-Y_Grid[12][0][0][0])**2.0), (Y_Grid[12][0][0][1]-Y_Grid[12][0][0][0])/math.sqrt((X_Grid[12][0][0][1]-X_Grid[12][0][0][0])**2.0+(Y_Grid[12][0][0][1]-Y_Grid[12][0][0][0])**2.0)))
        Theta_ele[12].append(((X_Grid[12][1][0][1]-X_Grid[12][1][0][0])/math.sqrt((X_Grid[12][1][0][1]-X_Grid[12][1][0][0])**2.0+(Y_Grid[12][1][0][1]-Y_Grid[12][1][0][0])**2.0), (Y_Grid[12][1][0][1]-Y_Grid[12][1][0][0])/math.sqrt((X_Grid[12][1][0][1]-X_Grid[12][1][0][0])**2.0+(Y_Grid[12][1][0][1]-Y_Grid[12][1][0][0])**2.0)))
        Theta_ele[12].append(((X_Grid[12][2][0][1]-X_Grid[12][2][0][0])/math.sqrt((X_Grid[12][2][0][1]-X_Grid[12][2][0][0])**2.0+(Y_Grid[12][2][0][1]-Y_Grid[12][2][0][0])**2.0), (Y_Grid[12][2][0][1]-Y_Grid[12][2][0][0])/math.sqrt((X_Grid[12][2][0][1]-X_Grid[12][2][0][0])**2.0+(Y_Grid[12][2][0][1]-Y_Grid[12][2][0][0])**2.0)))
        Theta_ele[12].append(((X_Grid[12][3][0][1]-X_Grid[12][3][0][0])/math.sqrt((X_Grid[12][3][0][1]-X_Grid[12][3][0][0])**2.0+(Y_Grid[12][3][0][1]-Y_Grid[12][3][0][0])**2.0), (Y_Grid[12][3][0][1]-Y_Grid[12][3][0][0])/math.sqrt((X_Grid[12][3][0][1]-X_Grid[12][3][0][0])**2.0+(Y_Grid[12][3][0][1]-Y_Grid[12][3][0][0])**2.0)))

        Center_ele[12].append(((X_Grid[12][0][0][0]+X_Grid[12][0][0][1]+X_Grid[12][0][-1][2]+X_Grid[12][0][-1][3])/4.0, (Y_Grid[12][0][0][0]+Y_Grid[12][0][0][1]+Y_Grid[12][0][-1][2]+Y_Grid[12][0][-1][3])/4.0))
        Center_ele[12].append(((X_Grid[12][1][0][0]+X_Grid[12][1][0][1]+X_Grid[12][1][-1][2]+X_Grid[12][1][-1][3])/4.0, (Y_Grid[12][1][0][0]+Y_Grid[12][1][0][1]+Y_Grid[12][1][-1][2]+Y_Grid[12][1][-1][3])/4.0))
        Center_ele[12].append(((X_Grid[12][2][0][0]+X_Grid[12][2][0][1]+X_Grid[12][2][-1][2]+X_Grid[12][2][-1][3])/4.0, (Y_Grid[12][2][0][0]+Y_Grid[12][2][0][1]+Y_Grid[12][2][-1][2]+Y_Grid[12][2][-1][3])/4.0))
        Center_ele[12].append(((X_Grid[12][3][0][0]+X_Grid[12][3][0][1]+X_Grid[12][3][-1][2]+X_Grid[12][3][-1][3])/4.0, (Y_Grid[12][3][0][0]+Y_Grid[12][3][0][1]+Y_Grid[12][3][-1][2]+Y_Grid[12][3][-1][3])/4.0))

    # WEB-03
    if num_webs==3:
        X_Grid.append([])                ; Y_Grid.append([])                ; C_Grid.append([])                ; Thick_lay.append([])                ; Length_lay.append([])                ; Ply_lay.append([])             ; Thick_ele.append([])         ; Length_ele.append([])         ; Theta_ele.append([])         ; Center_ele.append([])
        X_Grid[13].extend([[],[],[],[]]) ; Y_Grid[13].extend([[],[],[],[]]) ; C_Grid[13].extend([[],[],[],[]]) ; Thick_lay[13].extend([[],[],[],[]]) ; Length_lay[13].extend([[],[],[],[]]) ; Ply_lay[13].extend([[],[],[],[]])

        Xa=X_Grid[1][0][-1][3] ; Xb=X_Grid[1][0][-1][2] ; Xc=X_Grid[9][-1][-1][3] ; Xd=X_Grid[9][-1][-1][2]
        Ya=Y_Grid[1][0][-1][3] ; Yb=Y_Grid[1][0][-1][2] ; Yc=Y_Grid[9][-1][-1][3] ; Yd=Y_Grid[9][-1][-1][2]

        theta=angle([1.0, 0.0], [Xc-Xb,Yc-Yb])
        xloc=[(Xa-Xb)*np.cos(-theta)-(Ya-Yb)*np.sin(-theta), (Xb-Xb)*np.cos(-theta)-(Yb-Yb)*np.sin(-theta), (Xc-Xb)*np.cos(-theta)-(Yc-Yb)*np.sin(-theta), (Xd-Xb)*np.cos(-theta)-(Yd-Yb)*np.sin(-theta)]
        yloc=[(Xa-Xb)*np.sin(-theta)+(Ya-Yb)*np.cos(-theta), (Xb-Xb)*np.sin(-theta)+(Yb-Yb)*np.cos(-theta), (Xc-Xb)*np.sin(-theta)+(Yc-Yb)*np.cos(-theta), (Xd-Xb)*np.sin(-theta)+(Yd-Yb)*np.cos(-theta)]
        delta=0.0 ; xlocele=[xloc[1], xloc[2]] ; ylocele=[yloc[1], yloc[2]]
        for k in range(0,len(laminates[11])):
            if laminates[11][k]>0.0:
                delta=laminates[11][k]+delta
                xlocele.append(-delta*(xloc[1]-xloc[0])/(yloc[1]-yloc[0]))        ; ylocele.append(-delta)
                xlocele.append(xloc[2]-delta*(xloc[3]-xloc[2])/(yloc[3]-yloc[2])) ; ylocele.append(-delta)

                Length_lay[13][0].append(float((abs(xlocele[1]-xlocele[0])+abs(xlocele[3]-xlocele[2]))/8.0)) ; Length_lay[13][1].append(float((abs(xlocele[1]-xlocele[0])+abs(xlocele[3]-xlocele[2]))/8.0)) ; Length_lay[13][2].append(float((abs(xlocele[1]-xlocele[0])+abs(xlocele[3]-xlocele[2]))/8.0)) ; Length_lay[13][3].append(float((abs(xlocele[1]-xlocele[0])+abs(xlocele[3]-xlocele[2]))/8.0))
                Thick_lay[13][0].append(float(laminates[11][k]))                                             ; Thick_lay[13][1].append(float(laminates[11][k]))                                             ; Thick_lay[13][2].append(float(laminates[11][k]))                                             ; Thick_lay[13][3].append(float(laminates[11][k]))
                Ply_lay[13][0].append(float(0.0))                                                            ; Ply_lay[13][1].append(float(0.0))                                                            ; Ply_lay[13][2].append(float(0.0))                                                            ; Ply_lay[13][3].append(float(0.0))

                X_Grid[13][0].append([xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+0.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+1.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb+1.0*((xlocele[3]*np.cos(theta)-ylocele[3]*np.sin(theta)+Xb)-(xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb))/4.0, xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb+0.0*((xlocele[3]*np.cos(theta)-ylocele[3]*np.sin(theta)+Xb)-(xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb))/4.0, xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+0.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, None])
                Y_Grid[13][0].append([xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+0.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+1.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb+1.0*((xlocele[3]*np.sin(theta)+ylocele[3]*np.cos(theta)+Yb)-(xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb))/4.0, xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb+0.0*((xlocele[3]*np.sin(theta)+ylocele[3]*np.cos(theta)+Yb)-(xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb))/4.0, xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+0.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, None])

                X_Grid[13][1].append([xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+1.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+2.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb+2.0*((xlocele[3]*np.cos(theta)-ylocele[3]*np.sin(theta)+Xb)-(xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb))/4.0, xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb+1.0*((xlocele[3]*np.cos(theta)-ylocele[3]*np.sin(theta)+Xb)-(xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb))/4.0, xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+1.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, None])
                Y_Grid[13][1].append([xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+1.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+2.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb+2.0*((xlocele[3]*np.sin(theta)+ylocele[3]*np.cos(theta)+Yb)-(xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb))/4.0, xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb+1.0*((xlocele[3]*np.sin(theta)+ylocele[3]*np.cos(theta)+Yb)-(xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb))/4.0, xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+1.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, None])

                X_Grid[13][2].append([xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+2.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+3.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb+3.0*((xlocele[3]*np.cos(theta)-ylocele[3]*np.sin(theta)+Xb)-(xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb))/4.0, xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb+2.0*((xlocele[3]*np.cos(theta)-ylocele[3]*np.sin(theta)+Xb)-(xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb))/4.0, xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+2.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, None])
                Y_Grid[13][2].append([xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+2.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+3.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb+3.0*((xlocele[3]*np.sin(theta)+ylocele[3]*np.cos(theta)+Yb)-(xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb))/4.0, xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb+2.0*((xlocele[3]*np.sin(theta)+ylocele[3]*np.cos(theta)+Yb)-(xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb))/4.0, xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+2.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, None])

                X_Grid[13][3].append([xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+3.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+4.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb+4.0*((xlocele[3]*np.cos(theta)-ylocele[3]*np.sin(theta)+Xb)-(xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb))/4.0, xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb+3.0*((xlocele[3]*np.cos(theta)-ylocele[3]*np.sin(theta)+Xb)-(xlocele[2]*np.cos(theta)-ylocele[2]*np.sin(theta)+Xb))/4.0, xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb+3.0*((xlocele[1]*np.cos(theta)-ylocele[1]*np.sin(theta)+Xb)-(xlocele[0]*np.cos(theta)-ylocele[0]*np.sin(theta)+Xb))/4.0, None])
                Y_Grid[13][3].append([xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+3.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+4.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb+4.0*((xlocele[3]*np.sin(theta)+ylocele[3]*np.cos(theta)+Yb)-(xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb))/4.0, xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb+3.0*((xlocele[3]*np.sin(theta)+ylocele[3]*np.cos(theta)+Yb)-(xlocele[2]*np.sin(theta)+ylocele[2]*np.cos(theta)+Yb))/4.0, xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb+3.0*((xlocele[1]*np.sin(theta)+ylocele[1]*np.cos(theta)+Yb)-(xlocele[0]*np.sin(theta)+ylocele[0]*np.cos(theta)+Yb))/4.0, None])

                if k==0 or k==6:
                    C_Grid[13][0].append('green')  ; C_Grid[13][1].append('green')  ; C_Grid[13][2].append('green')  ; C_Grid[13][3].append('green')
                elif k==1 or k==5:
                    C_Grid[13][0].append('yellow') ; C_Grid[13][1].append('yellow') ; C_Grid[13][2].append('yellow') ; C_Grid[13][3].append('yellow')
                elif k==2 or k==4:
                    C_Grid[13][0].append('red')    ; C_Grid[13][1].append('red')    ; C_Grid[13][2].append('red')    ; C_Grid[13][3].append('red')
                else:
                    C_Grid[13][0].append('gray')   ; C_Grid[13][1].append('gray')   ; C_Grid[13][2].append('gray')   ; C_Grid[13][3].append('gray')
                xlocele=[xlocele[2], xlocele[3]] ; ylocele=[ylocele[2], ylocele[3]]
            else:
                Length_lay[13][0].append(0.0) ; Thick_lay[13][0].append(0.0) ; Ply_lay[13][0].append(0.0)
                Length_lay[13][1].append(0.0) ; Thick_lay[13][1].append(0.0) ; Ply_lay[13][1].append(0.0)
                Length_lay[13][2].append(0.0) ; Thick_lay[13][2].append(0.0) ; Ply_lay[13][2].append(0.0)
                Length_lay[13][3].append(0.0) ; Thick_lay[13][3].append(0.0) ; Ply_lay[13][3].append(0.0)

        Thick_ele[13].append(sum(laminates[11][:]))
        Thick_ele[13].append(sum(laminates[11][:]))
        Thick_ele[13].append(sum(laminates[11][:]))
        Thick_ele[13].append(sum(laminates[11][:]))

        Length_ele[13].append((math.sqrt((X_Grid[13][0][0][1]-X_Grid[13][0][0][0])**2.0+(Y_Grid[13][0][0][1]-Y_Grid[13][0][0][0])**2.0)+math.sqrt((X_Grid[13][0][-1][3]-X_Grid[13][0][-1][2])**2.0+(Y_Grid[13][0][-1][3]-Y_Grid[13][0][-1][2])**2.0))/2.0)
        Length_ele[13].append((math.sqrt((X_Grid[13][1][0][1]-X_Grid[13][1][0][0])**2.0+(Y_Grid[13][1][0][1]-Y_Grid[13][1][0][0])**2.0)+math.sqrt((X_Grid[13][1][-1][3]-X_Grid[13][1][-1][2])**2.0+(Y_Grid[13][1][-1][3]-Y_Grid[13][1][-1][2])**2.0))/2.0)
        Length_ele[13].append((math.sqrt((X_Grid[13][2][0][1]-X_Grid[13][2][0][0])**2.0+(Y_Grid[13][2][0][1]-Y_Grid[13][2][0][0])**2.0)+math.sqrt((X_Grid[13][2][-1][3]-X_Grid[13][2][-1][2])**2.0+(Y_Grid[13][2][-1][3]-Y_Grid[13][2][-1][2])**2.0))/2.0)
        Length_ele[13].append((math.sqrt((X_Grid[13][3][0][1]-X_Grid[13][3][0][0])**2.0+(Y_Grid[13][3][0][1]-Y_Grid[13][3][0][0])**2.0)+math.sqrt((X_Grid[13][3][-1][3]-X_Grid[13][3][-1][2])**2.0+(Y_Grid[13][3][-1][3]-Y_Grid[13][3][-1][2])**2.0))/2.0)

        Theta_ele[13].append(((X_Grid[13][0][0][1]-X_Grid[13][0][0][0])/math.sqrt((X_Grid[13][0][0][1]-X_Grid[13][0][0][0])**2.0+(Y_Grid[13][0][0][1]-Y_Grid[13][0][0][0])**2.0), (Y_Grid[13][0][0][1]-Y_Grid[13][0][0][0])/math.sqrt((X_Grid[13][0][0][1]-X_Grid[13][0][0][0])**2.0+(Y_Grid[13][0][0][1]-Y_Grid[13][0][0][0])**2.0)))
        Theta_ele[13].append(((X_Grid[13][1][0][1]-X_Grid[13][1][0][0])/math.sqrt((X_Grid[13][1][0][1]-X_Grid[13][1][0][0])**2.0+(Y_Grid[13][1][0][1]-Y_Grid[13][1][0][0])**2.0), (Y_Grid[13][1][0][1]-Y_Grid[13][1][0][0])/math.sqrt((X_Grid[13][1][0][1]-X_Grid[13][1][0][0])**2.0+(Y_Grid[13][1][0][1]-Y_Grid[13][1][0][0])**2.0)))
        Theta_ele[13].append(((X_Grid[13][2][0][1]-X_Grid[13][2][0][0])/math.sqrt((X_Grid[13][2][0][1]-X_Grid[13][2][0][0])**2.0+(Y_Grid[13][2][0][1]-Y_Grid[13][2][0][0])**2.0), (Y_Grid[13][2][0][1]-Y_Grid[13][2][0][0])/math.sqrt((X_Grid[13][2][0][1]-X_Grid[13][2][0][0])**2.0+(Y_Grid[13][2][0][1]-Y_Grid[13][2][0][0])**2.0)))
        Theta_ele[13].append(((X_Grid[13][3][0][1]-X_Grid[13][3][0][0])/math.sqrt((X_Grid[13][3][0][1]-X_Grid[13][3][0][0])**2.0+(Y_Grid[13][3][0][1]-Y_Grid[13][3][0][0])**2.0), (Y_Grid[13][3][0][1]-Y_Grid[13][3][0][0])/math.sqrt((X_Grid[13][3][0][1]-X_Grid[13][3][0][0])**2.0+(Y_Grid[13][3][0][1]-Y_Grid[13][3][0][0])**2.0)))

        Center_ele[13].append(((X_Grid[13][0][0][0]+X_Grid[13][0][0][1]+X_Grid[13][0][-1][2]+X_Grid[13][0][-1][3])/4.0, (Y_Grid[13][0][0][0]+Y_Grid[13][0][0][1]+Y_Grid[13][0][-1][2]+Y_Grid[13][0][-1][3])/4.0))
        Center_ele[13].append(((X_Grid[13][1][0][0]+X_Grid[13][1][0][1]+X_Grid[13][1][-1][2]+X_Grid[13][1][-1][3])/4.0, (Y_Grid[13][1][0][0]+Y_Grid[13][1][0][1]+Y_Grid[13][1][-1][2]+Y_Grid[13][1][-1][3])/4.0))
        Center_ele[13].append(((X_Grid[13][2][0][0]+X_Grid[13][2][0][1]+X_Grid[13][2][-1][2]+X_Grid[13][2][-1][3])/4.0, (Y_Grid[13][2][0][0]+Y_Grid[13][2][0][1]+Y_Grid[13][2][-1][2]+Y_Grid[13][2][-1][3])/4.0))
        Center_ele[13].append(((X_Grid[13][3][0][0]+X_Grid[13][3][0][1]+X_Grid[13][3][-1][2]+X_Grid[13][3][-1][3])/4.0, (Y_Grid[13][3][0][0]+Y_Grid[13][3][0][1]+Y_Grid[13][3][-1][2]+Y_Grid[13][3][-1][3])/4.0))

    return Thick_lay, Length_lay, Ply_lay, Thick_ele, Length_ele, Theta_ele, Center_ele, X_Grid, Y_Grid, C_Grid

# Check Grid
#################################################
def check_grid(X_Grid, Y_Grid, C_Grid, Check):

    if Check==True:
        # Check over-lap
        for i in range(0,len(X_Grid)):
            for j in range(0,len(X_Grid[i])):
                for k in range(0,len(X_Grid[i][j])):
                    v1=[X_Grid[i][j][k][0]-X_Grid[i][j][k][1], Y_Grid[i][j][k][0]-Y_Grid[i][j][k][1]] ; v2=[X_Grid[i][j][k][2]-X_Grid[i][j][k][1], Y_Grid[i][j][k][2]-Y_Grid[i][j][k][1]] ; nor_1=np.cross(v1,v2)
                    v1=[X_Grid[i][j][k][1]-X_Grid[i][j][k][2], Y_Grid[i][j][k][1]-Y_Grid[i][j][k][2]] ; v2=[X_Grid[i][j][k][3]-X_Grid[i][j][k][2], Y_Grid[i][j][k][3]-Y_Grid[i][j][k][2]] ; nor_2=np.cross(v1,v2)
                    if nor_1*nor_2<0.0:
                        Check=False
                        C_Grid[i][j][k]='black'
                        print('Error: can not create grid..... Jacobian is false')

                    for l in range(0,len(X_Grid)):
                        for m in range(0,len(X_Grid[l])):
                            for n in range(0,len(X_Grid[l][m])):
                                for o in range(0,4):
                                    if not(i==l and j==m and k==n):
                                        v1=[X_Grid[i][j][k][0]-X_Grid[l][m][n][o], Y_Grid[i][j][k][0]-Y_Grid[l][m][n][o]] ; v2=[X_Grid[i][j][k][1]-X_Grid[l][m][n][o], Y_Grid[i][j][k][1]-Y_Grid[l][m][n][o]] ; nor_1=np.cross(v1,v2)
                                        v1=[X_Grid[i][j][k][1]-X_Grid[l][m][n][o], Y_Grid[i][j][k][1]-Y_Grid[l][m][n][o]] ; v2=[X_Grid[i][j][k][2]-X_Grid[l][m][n][o], Y_Grid[i][j][k][2]-Y_Grid[l][m][n][o]] ; nor_2=np.cross(v1,v2)
                                        v1=[X_Grid[i][j][k][2]-X_Grid[l][m][n][o], Y_Grid[i][j][k][2]-Y_Grid[l][m][n][o]] ; v2=[X_Grid[i][j][k][3]-X_Grid[l][m][n][o], Y_Grid[i][j][k][3]-Y_Grid[l][m][n][o]] ; nor_3=np.cross(v1,v2)
                                        v1=[X_Grid[i][j][k][3]-X_Grid[l][m][n][o], Y_Grid[i][j][k][3]-Y_Grid[l][m][n][o]] ; v2=[X_Grid[i][j][k][4]-X_Grid[l][m][n][o], Y_Grid[i][j][k][4]-Y_Grid[l][m][n][o]] ; nor_4=np.cross(v1,v2)
                                        if nor_1<-0.0001 and nor_2<-0.0001 and nor_3<-0.0001 and nor_4<-0.0001:
                                            Check=False
                                            C_Grid[i][j][k]='black'
                                            print(nor_1, nor_2, nor_3, nor_4)
                                            print('Error: can not create grid..... Over-lap elements')

    plt.fill_between([], [], [], color='green',  alpha=0.75, label='TRIAX')
    plt.fill_between([], [], [], color='yellow', alpha=0.75, label='BIAX' )
    plt.fill_between([], [], [], color='red',    alpha=0.75, label='UNIAX')
    plt.fill_between([], [], [], color='gray',   alpha=0.75, label='BALSA')
    for i in range(0,len(X_Grid)):
        for j in range(0,len(X_Grid[i])):
            for k in range(0,len(X_Grid[i][j])):
                plt.fill(X_Grid[i][j][k][:], Y_Grid[i][j][k][:], color=C_Grid[i][j][k], alpha=0.75)
                if C_Grid[i][j][k]=='black': plt.fill(X_Grid[i][j][k][:], Y_Grid[i][j][k][:], color=C_Grid[i][j][k], alpha=1.00)
    plt.grid()
    plt.legend()
    plt.show()

  # return Check

# Mass
#################################################
def mass(Thick_lay, Thick_ele, Length_ele, Center_ele, triax, biax, uniax, balsa):
    den=[]
    for i in range(0,len(Thick_lay)):
        den.append([])
        for j in range(0,len(Thick_lay[i])):
            den[i].append((Thick_lay[i][j][0]*triax[1]+Thick_lay[i][j][1]*biax[1]+Thick_lay[i][j][2]*uniax[1]+Thick_lay[i][j][3]*balsa[1]+Thick_lay[i][j][4]*uniax[1]+Thick_lay[i][j][5]*biax[1]+Thick_lay[i][j][6]*triax[1])/sum(Thick_lay[i][j][:]))

    m3x3=np.zeros((3,3),dtype=float) ; x_m=0.0 ; y_m=0.0
    for i in range(0,len(Thick_ele)):
        for j in range(0,len(Thick_ele[i])):
            m3x3=m3x3+Thick_ele[i][j]*Length_ele[i][j]*den[i][j]*np.array([[1.0, 0.0, 0.0], [0.0, Center_ele[i][j][1]**2.0, Center_ele[i][j][0]*Center_ele[i][j][1]], [0.0, 0.0, Center_ele[i][j][0]**2.0]])
            x_m=x_m+Thick_ele[i][j]*Length_ele[i][j]*den[i][j]*Center_ele[i][j][0] ; y_m=y_m+Thick_ele[i][j]*Length_ele[i][j]*den[i][j]*Center_ele[i][j][1]

    Mass=np.array([[m3x3[0][0], 0.0, 0.0, 0.0, m3x3[0][0]*y_m, 0.0], [0.0, m3x3[0][0], 0.0, -m3x3[0][0]*y_m, 0.0, m3x3[0][0]*x_m], [0.0, 0.0, m3x3[0][0], 0.0, -m3x3[0][0]*x_m, 0.0], [0.0, -m3x3[0][0]*y_m, 0.0, m3x3[1][1], 0.0, -m3x3[1][2]], [m3x3[0][0]*y_m, 0.0, -m3x3[0][0]*x_m, 0.0, m3x3[1][1]+m3x3[2][2], 0.0], [0.0, m3x3[0][0]*x_m, 0.0, -m3x3[1][2], 0.0, m3x3[2][2]]])
    return Mass, (x_m/m3x3[0][0], y_m/m3x3[0][0])

# Timoshenko Shear Factor
#################################################
def Tim_fac(Thick_ele, Length_ele, Center_ele, sc, v):
    GA=0.0 ; Ixx=0.0 ; Iyy=0.0
    K_mat=np.zeros((9,9),dtype=float) ; F1_mat=np.zeros((9,1),dtype=float) ; F2_mat=np.zeros((9,1),dtype=float) ; F3_mat=np.zeros((9,1),dtype=float)
    for i in range(0,11):
        for j in range(0,len(Length_ele[i])):
            x=Center_ele[i][j][0]-sc[0] ; y=Center_ele[i][j][1]-sc[1]
            N=np.array([[x, y, x**2.0, x*y, y**2.0, x**3.0, (x**2.0)*y, x*(y**2.0), y**3.0]]) ; Nx=np.array([[1.0, 0.0, 2.0*x, y, 0.0, 3.0*(x**2.0), 2.0*x*y, y**2.0, 0.0]]) ; Ny=np.array([[0.0, 1.0, 0.0, x, 2.0*y, 0.0, x**2.0, 2.0*x*y, 3.0*(y**2.0)]])
            GA=GA+Length_ele[i][j]*Thick_ele[i][j]/(2.0*(v+1.0)) ; Ixx=Ixx+Length_ele[i][j]*Thick_ele[i][j]*(y**2.0) ; Iyy=Iyy+Length_ele[i][j]*Thick_ele[i][j]*(x**2.0)
            K_mat=K_mat+Length_ele[i][j]*Thick_ele[i][j]*(np.matmul(Nx.transpose(),Nx)+np.matmul(Ny.transpose(),Ny))/(2.0*(v+1.0))
            F1_mat=F1_mat+Length_ele[i][j]*Thick_ele[i][j]*(x*N.transpose()+(v*x*y*Ny.transpose()+0.50*v*(x**2.0-y**2.0)*Nx.transpose())/(2.0*(v+1.0)))
            F2_mat=F2_mat+Length_ele[i][j]*Thick_ele[i][j]*(y*N.transpose()+(0.50*v*(y**2.0-x**2.0)*Ny.transpose()+v*x*y*Nx.transpose())/(2.0*(v+1.0)))
            F3_mat=F3_mat+Length_ele[i][j]*Thick_ele[i][j]*(-x*Ny.transpose()+y*Nx.transpose())/(2.0*(v+1.0))
    c1_mat=np.linalg.solve(K_mat,F1_mat) ; c2_mat=np.linalg.solve(K_mat,F2_mat) ; c3_mat=np.linalg.solve(K_mat,F3_mat)

    int_11=0.0 ; int_21=0.0 ; int_31=0.0
    int_12=0.0 ; int_22=0.0 ; int_32=0.0
    for i in range(0,11):
        for j in range(0,len(Length_ele[i])):
            x=Center_ele[i][j][0]-sc[0] ; y=Center_ele[i][j][1]-sc[1]
            N=np.array([[x, y, x**2.0, x*y, y**2.0, x**3.0, (x**2.0)*y, x*(y**2.0), y**3.0]])
            int_11=int_11+Length_ele[i][j]*Thick_ele[i][j]*x*np.matmul(N,c1_mat)/Iyy ; int_12=int_12+Length_ele[i][j]*Thick_ele[i][j]*y*np.matmul(N,c1_mat)/Ixx
            int_21=int_21+Length_ele[i][j]*Thick_ele[i][j]*x*np.matmul(N,c2_mat)/Iyy ; int_22=int_22+Length_ele[i][j]*Thick_ele[i][j]*y*np.matmul(N,c2_mat)/Ixx
            int_31=int_31+Length_ele[i][j]*Thick_ele[i][j]*x*np.matmul(N,c3_mat)/Iyy ; int_32=int_32+Length_ele[i][j]*Thick_ele[i][j]*y*np.matmul(N,c3_mat)/Ixx

    s_mat=np.zeros((3,3),dtype=float) ; r_mat=np.zeros((3,3),dtype=float) ; res=np.zeros((3,3),dtype=float)
    for i in range(0,11):
        for j in range(0,len(Length_ele[i])):
            x=Center_ele[i][j][0]-sc[0] ; y=Center_ele[i][j][1]-sc[1]
            N=np.array([[x, y, x**2.0, x*y, y**2.0, x**3.0, (x**2.0)*y, x*(y**2.0), y**3.0]]) ; Nx=np.array([[1.0, 0.0, 2.0*x, y, 0.0, 3.0*(x**2.0), 2.0*x*y, y**2.0, 0.0]]) ; Ny=np.array([[0.0, 1.0, 0.0, x, 2.0*y, 0.0, x**2.0, 2.0*x*y, 3.0*(y**2.0)]])

            psz1dx=float(np.matmul(Nx,c1_mat)-int_11) ; psz1dy=float(np.matmul(Ny,c1_mat)-int_12)
            psz2dx=float(np.matmul(Nx,c2_mat)-int_21) ; psz2dy=float(np.matmul(Ny,c2_mat)-int_22)
            psz3dx=float(np.matmul(Nx,c3_mat)-int_31) ; psz3dy=float(np.matmul(Ny,c3_mat)-int_32)

            s_mat[0][0]=s_mat[0][0]+Length_ele[i][j]*Thick_ele[i][j]/(2.0*(v+1.0)) ; s_mat[0][1]=0.0                                                        ; s_mat[0][2]=s_mat[0][2]+Length_ele[i][j]*Thick_ele[i][j]*(psz3dx-y)/(2.0*(v+1.0))
            s_mat[1][0]=0.0                                                        ; s_mat[1][1]=s_mat[1][1]+Length_ele[i][j]*Thick_ele[i][j]/(2.0*(v+1.0)) ; s_mat[1][2]=s_mat[1][2]+Length_ele[i][j]*Thick_ele[i][j]*(psz3dy+x)/(2.0*(v+1.0))
            s_mat[2][0]=0.0                                                        ; s_mat[2][1]=0.0                                                        ; s_mat[2][2]=s_mat[2][2]+Length_ele[i][j]*Thick_ele[i][j]*(x*(x+psz3dy)-y*(psz3dx-y))/(2.0*(v+1.0))

            r_mat[0][0]=r_mat[0][0]+Length_ele[i][j]*Thick_ele[i][j]*(psz1dx)/(2.0*(v+1.0)*Iyy)            ; r_mat[0][1]=r_mat[0][1]+Length_ele[i][j]*Thick_ele[i][j]*(psz2dx)/(2.0*(v+1.0)*Ixx)            ; r_mat[0][2]=0.0
            r_mat[1][0]=r_mat[1][0]+Length_ele[i][j]*Thick_ele[i][j]*(psz1dy)/(2.0*(v+1.0)*Iyy)            ; r_mat[1][1]=r_mat[1][1]+Length_ele[i][j]*Thick_ele[i][j]*(psz2dy)/(2.0*(v+1.0)*Ixx)            ; r_mat[1][2]=0.0
            r_mat[2][0]=r_mat[2][0]+Length_ele[i][j]*Thick_ele[i][j]*(x*psz1dy-y*psz1dx)/(2.0*(v+1.0)*Iyy) ; r_mat[2][1]=r_mat[2][1]+Length_ele[i][j]*Thick_ele[i][j]*(x*psz2dy-y*psz2dx)/(2.0*(v+1.0)*Ixx) ; r_mat[2][2]=1.0
    r_mat[0][0]=1.0-r_mat[0][0] ; r_mat[0][1]=-r_mat[0][1]    ; r_mat[0][2]=0.0
    r_mat[1][0]=-r_mat[1][0]    ; r_mat[1][1]=1.0-r_mat[1][1] ; r_mat[1][2]=0.0
    r_mat[2][0]=-r_mat[2][0]    ; r_mat[2][1]=-r_mat[2][1]    ; r_mat[2][2]=1.0
    res=np.matmul(inv(r_mat),s_mat)

   #return (res[1][1]/GA, res[0][0]/GA)
    return (res[0][0]/GA, res[1][1]/GA)

# Stiffness
#################################################
def stiffness(Thick_lay, Length_lay, Ply_lay, Thick_ele, Length_ele, Theta_ele, Center_ele, triax, biax, uniax, balsa, sc, num_webs, TimFac, F):
    sf = 1.125
    E11=[triax[2], biax[2], uniax[2], balsa[2], uniax[2], biax[2], triax[2]] ; E22=[triax[3], biax[3], uniax[3], balsa[3], uniax[3], biax[3], triax[3]] ; v12=[triax[4], biax[4], uniax[4], balsa[4], uniax[4], biax[4], triax[4]]
    G12=[triax[5], biax[5], uniax[5], balsa[5], uniax[5], biax[5], triax[5]] ; G13=[triax[6], biax[6], uniax[6], balsa[6], uniax[6], biax[6], triax[6]] ; G23=[triax[7], biax[7], uniax[7], balsa[7], uniax[7], biax[7], triax[7]]
    xc=[triax[8]/sf,biax[8]/sf,uniax[8]/sf,balsa[8]/sf,uniax[8]/sf,biax[8]/sf,triax[8]/sf]        ; xt=[triax[9]/sf,biax[9]/sf,uniax[9]/sf,balsa[9]/sf,uniax[9]/sf,biax[9]/sf,triax[9]/sf]
    yc=[triax[10]/sf,biax[10]/sf,uniax[10]/sf,balsa[10]/sf,uniax[10]/sf,biax[10]/sf,triax[10]/sf] ; yt=[triax[11]/sf,biax[11]/sf,uniax[11]/sf,balsa[11]/sf,uniax[11]/sf,biax[11]/sf,triax[11]/sf]
    st=[triax[12]/sf,biax[12]/sf,uniax[12]/sf,balsa[12]/sf,uniax[12]/sf,biax[12]/sf,triax[12]/sf]

    TSF = Tim_fac(Thick_ele, Length_ele, Center_ele, (0.0,0.0), (sum(v12[:]))/6.0) if TimFac==True else (1.0, 1.0)

    s5x5=np.zeros((5,5),dtype=float) ; s3x3=np.zeros((3,3),dtype=float) ; r5x5=np.zeros((5,5),dtype=float) ; r3x3=np.zeros((3,3),dtype=float) ; q3x3=[]
    for i in range(0,len(Thick_lay)):
        q3x3.append([])
        for j in range(0,len(Thick_lay[i])):
            q3x3[i].append([])
            for k in range(0,len(Thick_lay[i][j])):
                if Ply_lay[i][j][k]==0.0:
                    q3x3[i][j].append([[E11[k], 0.0, 0.0],[0.0, G13[k], 0.0],[0.0, 0.0, G12[k]]])
                else:
                    theta=Ply_lay[i][j][k]*math.pi/180.0
                    s5x5[0][0]=1.0/E11[k]     ; s5x5[0][1]=-v12[k]/E11[k] ; s5x5[0][2]=0.0        ; s5x5[0][3]=0.0        ; s5x5[0][4]=0.0
                    s5x5[1][0]=-v12[k]/E11[k] ; s5x5[1][1]=1.0/E22[k]     ; s5x5[1][2]=0.0        ; s5x5[1][3]=0.0        ; s5x5[1][4]=0.0
                    s5x5[2][0]=0.0            ; s5x5[2][1]=0.0            ; s5x5[2][2]=1.0/G23[k] ; s5x5[2][3]=0.0        ; s5x5[2][4]=0.0
                    s5x5[3][0]=0.0            ; s5x5[3][1]=0.0            ; s5x5[3][2]=0.0        ; s5x5[3][3]=1.0/G13[k] ; s5x5[3][4]=0.0
                    s5x5[4][0]=0.0            ; s5x5[4][1]=0.0            ; s5x5[4][2]=0.0        ; s5x5[4][3]=0.0        ; s5x5[4][4]=1.0/G12[k]

                    s3x3[0][0]=1.0/E11[k] ; s3x3[0][1]=0.0        ; s3x3[0][2]=0.0
                    s3x3[1][0]=0.0        ; s3x3[1][1]=1.0/G13[k] ; s3x3[1][2]=0.0
                    s3x3[2][0]=0.0        ; s3x3[2][1]=0.0        ; s3x3[2][2]=1.0/G12[k]

                    r5x5[0][0]=(np.cos(theta))**2.0             ; r5x5[0][1]=(np.sin(theta))**2.0            ; r5x5[0][2]=0.0           ; r5x5[0][3]=0.0            ; r5x5[0][4]=2.0*(np.cos(theta))*(np.sin(theta))
                    r5x5[1][0]=(np.sin(theta))**2.0             ; r5x5[1][1]=(np.cos(theta))**2.0            ; r5x5[1][2]=0.0           ; r5x5[1][3]=0.0            ; r5x5[1][4]=-2.0*(np.cos(theta))*(np.sin(theta))
                    r5x5[2][0]=0.0                              ; r5x5[2][1]=0.0                             ; r5x5[2][2]=np.cos(theta) ; r5x5[2][3]=-np.sin(theta) ; r5x5[2][4]=0.0
                    r5x5[3][0]=0.0                              ; r5x5[3][1]=0.0                             ; r5x5[3][2]=np.sin(theta) ; r5x5[3][3]=np.cos(theta)  ; r5x5[3][4]=0.0
                    r5x5[4][0]=-(np.cos(theta))*(np.sin(theta)) ; r5x5[4][1]=(np.cos(theta))*(np.sin(theta)) ; r5x5[4][2]=0.0           ; r5x5[4][3]=0.0            ; r5x5[4][4]=(np.cos(theta))**2.0-(np.sin(theta))**2.0

                    r3x3[0][0]=(np.cos(theta))**2.0             ; r3x3[0][1]=0.0            ; r3x3[0][2]=2.0*(np.cos(theta))*(np.sin(theta))
                    r3x3[1][0]=0.0                              ; r3x3[1][1]=np.cos(theta)  ; r3x3[1][2]=0.0
                    r3x3[2][0]=-(np.cos(theta))*(np.sin(theta)) ; r3x3[2][1]=0.0            ; r3x3[2][2]=(np.cos(theta))**2.0-(np.sin(theta))**2.0


                    q5x5=inv(np.matmul(np.matmul(inv(r5x5).transpose(),s5x5),inv(r5x5)))
                    q5x5t=inv(np.matmul(np.matmul(inv(r3x3).transpose(),s3x3),inv(r3x3)))
                    q3x3[i][j].append([[q5x5t[0][0], 0.0, q5x5[0][4]], [0.0, q5x5[3][3], 0.0], [q5x5[4][0], 0.0, q5x5[4][4]]])

    A3x3=[] ; B3x3=[] ; D3x3=[]
    for i in range(0,len(Thick_lay)):
        A3x3.append([]) ; B3x3.append([]) ; D3x3.append([])
        for j in range(0,len(Thick_lay[i])):
            A=np.zeros((3,3),dtype=float) ; B=np.zeros((3,3),dtype=float) ; D=np.zeros((3,3),dtype=float)
            for k in range(0,len(Thick_lay[i][j])):
                if Thick_lay[i][j][k]>0.0:
                    K=k
                    break
            for k in range(0,len(Thick_lay[i][j])):
                if Thick_lay[i][j][k]>0.0:
                    z= (Thick_lay[i][j][k]-sum(Thick_lay[i][j][:]))/2.0 if k==K else z+(thick_lay+Thick_lay[i][j][k])/2.0
                    thick_lay=Thick_lay[i][j][k]
                    A=A+np.array(q3x3[i][j][k])*Thick_lay[i][j][k]
                    B=B+np.array(q3x3[i][j][k])*Thick_lay[i][j][k]*z
                    D=D+np.array(q3x3[i][j][k])*Thick_lay[i][j][k]*((Thick_lay[i][j][k]**2.0)/12.0+z**2.0)
            A3x3[i].append([[A[0][0],0.0,A[0][2]],[0.0, A[1][1], 0.0],[A[2][0], 0.0, A[2][2]]])
            B3x3[i].append([[B[0][0], B[0][1], B[0][2]],[0.0, 0.0, 0.0],[B[2][0], B[2][1], B[2][2]]])
            D3x3[i].append([[D[0][0], D[0][1], D[0][2]],[D[1][0], D[1][1], D[1][2]],[D[2][0], D[2][1], D[2][2]]])
    #--------------------------------------------------
    if num_webs==0:
        Area=0.0 ; Beta=0.0 ; lamda=0.0
        for i in range(0,len(Length_ele)):
            for j in range(0,len(Length_ele[i])):
                Area, Beta, lamda = Area+Length_ele[i][j]*(-Center_ele[i][j][1]*Theta_ele[i][j][0]+Center_ele[i][j][0]*Theta_ele[i][j][1]), Beta+2.0*Length_ele[i][j]*B3x3[i][j][2][2]/A3x3[i][j][2][2], lamda+Length_ele[i][j]/A3x3[i][j][2][2]
        Cah=-(Area+Beta)/lamda
    #--------------------------------------------------
    if num_webs==1:
        Area=np.array([0.0, 0.0]) ; Beta=np.array([0.0, 0.0]) ; lamda=np.zeros((2,2),dtype=float)
        for i in [6,5,4,11]:
            for j in range(0,len(Length_ele[i])):
                Area[0], Beta[0], lamda[0][0] = Area[0]+Length_ele[i][j]*(-Center_ele[i][j][1]*Theta_ele[i][j][0]+Center_ele[i][j][0]*Theta_ele[i][j][1]), Beta[0]+2.0*Length_ele[i][j]*B3x3[i][j][2][2]/A3x3[i][j][2][2], lamda[0][0]+Length_ele[i][j]/A3x3[i][j][2][2]
        for i in [10,9,8,7,11,3,2,1,0]:
            for j in range(0,len(Length_ele[i])):
                Area[1] = Area[1]-Length_ele[i][j]*(-Center_ele[i][j][1]*Theta_ele[i][j][0]+Center_ele[i][j][0]*Theta_ele[i][j][1]) if i==11 else Area[1]+Length_ele[i][j]*(-Center_ele[i][j][1]*Theta_ele[i][j][0]+Center_ele[i][j][0]*Theta_ele[i][j][1])
                Beta[1], lamda[1][1] = Beta[1]+2.0*Length_ele[i][j]*B3x3[i][j][2][2]/A3x3[i][j][2][2], lamda[1][1]+Length_ele[i][j]/A3x3[i][j][2][2]
        for j in range(0,len(Length_ele[11])):
            lamda[0][1], lamda[1][0] = lamda[0][1]-Length_ele[11][j]/A3x3[11][j][2][2], lamda[1][0]-Length_ele[11][j]/A3x3[11][j][2][2]
        Cah=np.linalg.solve(lamda,-(Area+Beta))
    #--------------------------------------------------
    if num_webs==2:
        Area=np.array([0.0, 0.0, 0.0]) ; Beta=np.array([0.0, 0.0, 0.0]) ; lamda=np.zeros((3,3),dtype=float)
        for i in [6,5,4,11]:
            for j in range(0,len(Length_ele[i])):
                Area[0], Beta[0], lamda[0][0] = Area[0]+Length_ele[i][j]*(-Center_ele[i][j][1]*Theta_ele[i][j][0]+Center_ele[i][j][0]*Theta_ele[i][j][1]), Beta[0]+2.0*Length_ele[i][j]*B3x3[i][j][2][2]/A3x3[i][j][2][2], lamda[0][0]+Length_ele[i][j]/A3x3[i][j][2][2]
        for i in [7, 11, 3, 12]:
            for j in range(0,len(Length_ele[i])):
                Area[1] = Area[1]-Length_ele[i][j]*(-Center_ele[i][j][1]*Theta_ele[i][j][0]+Center_ele[i][j][0]*Theta_ele[i][j][1]) if i==11 or i==12 else Area[1]+Length_ele[i][j]*(-Center_ele[i][j][1]*Theta_ele[i][j][0]+Center_ele[i][j][0]*Theta_ele[i][j][1])
                Beta[1], lamda[1][1] = Beta[1]+2.0*Length_ele[i][j]*B3x3[i][j][2][2]/A3x3[i][j][2][2], lamda[1][1]+Length_ele[i][j]/A3x3[i][j][2][2]
        for i in [10, 9, 8, 12, 2, 1, 0]:
            for j in range(0,len(Length_ele[i])):
                Area[2], Beta[2], lamda[2][2] = Area[2]+Length_ele[i][j]*(-Center_ele[i][j][1]*Theta_ele[i][j][0]+Center_ele[i][j][0]*Theta_ele[i][j][1]), Beta[2]+2.0*Length_ele[i][j]*B3x3[i][j][2][2]/A3x3[i][j][2][2], lamda[2][2]+Length_ele[i][j]/A3x3[i][j][2][2]
        for j in range(0,len(Length_ele[11])):
            lamda[0][1], lamda[1][0] = lamda[0][1]-Length_ele[11][j]/A3x3[11][j][2][2], lamda[1][0]-Length_ele[11][j]/A3x3[11][j][2][2]
        for j in range(0,len(Length_ele[12])):
            lamda[1][2], lamda[2][1] = lamda[1][2]-Length_ele[12][j]/A3x3[12][j][2][2], lamda[2][1]-Length_ele[12][j]/A3x3[12][j][2][2]
        Cah=np.linalg.solve(lamda,-(Area+Beta))
    #--------------------------------------------------
    if num_webs==3:
        Area=np.array([0.0, 0.0, 0.0, 0.0]) ; Beta=np.array([0.0, 0.0, 0.0, 0.0]) ; lamda=np.zeros((4,4),dtype=float)
        for i in [6,5,4,11]:
            for j in range(0,len(Length_ele[i])):
                Area[0], Beta[0], lamda[0][0] = Area[0]+Length_ele[i][j]*(-Center_ele[i][j][1]*Theta_ele[i][j][0]+Center_ele[i][j][0]*Theta_ele[i][j][1]), Beta[0]+2.0*Length_ele[i][j]*B3x3[i][j][2][2]/A3x3[i][j][2][2], lamda[0][0]+Length_ele[i][j]/A3x3[i][j][2][2]
        for i in [7, 11, 3, 12]:
            for j in range(0,len(Length_ele[i])):
                Area[1] = Area[1]-Length_ele[i][j]*(-Center_ele[i][j][1]*Theta_ele[i][j][0]+Center_ele[i][j][0]*Theta_ele[i][j][1]) if i==11 or i==12 else Area[1]+Length_ele[i][j]*(-Center_ele[i][j][1]*Theta_ele[i][j][0]+Center_ele[i][j][0]*Theta_ele[i][j][1])
                Beta[1], lamda[1][1] = Beta[1]+2.0*Length_ele[i][j]*B3x3[i][j][2][2]/A3x3[i][j][2][2], lamda[1][1]+Length_ele[i][j]/A3x3[i][j][2][2]
        for i in [9, 8, 12, 2, 1, 13]:
            for j in range(0,len(Length_ele[i])):
                Area[2]=Area[2]-Length_ele[i][j]*(-Center_ele[i][j][1]*Theta_ele[i][j][0]+Center_ele[i][j][0]*Theta_ele[i][j][1]) if i==13 else Area[2]+Length_ele[i][j]*(-Center_ele[i][j][1]*Theta_ele[i][j][0]+Center_ele[i][j][0]*Theta_ele[i][j][1])
                Beta[2], lamda[2][2] = Beta[2]+2.0*Length_ele[i][j]*B3x3[i][j][2][2]/A3x3[i][j][2][2], lamda[2][2]+Length_ele[i][j]/A3x3[i][j][2][2]
        for i in [10, 13, 0]:
            for j in range(0,len(Length_ele[i])):
                Area[3], Beta[3], lamda[3][3] = Area[3]+Length_ele[i][j]*(-Center_ele[i][j][1]*Theta_ele[i][j][0]+Center_ele[i][j][0]*Theta_ele[i][j][1]), Beta[3]+2.0*Length_ele[i][j]*B3x3[i][j][2][2]/A3x3[i][j][2][2], lamda[3][3]+Length_ele[i][j]/A3x3[i][j][2][2]

        for j in range(0,len(Length_ele[11])):
            lamda[0][1], lamda[1][0] = lamda[0][1]-Length_ele[11][j]/A3x3[11][j][2][2], lamda[1][0]-Length_ele[11][j]/A3x3[11][j][2][2]
        for j in range(0,len(Length_ele[12])):
            lamda[1][2], lamda[2][1] = lamda[1][2]-Length_ele[12][j]/A3x3[12][j][2][2], lamda[2][1]-Length_ele[12][j]/A3x3[12][j][2][2]
        for j in range(0,len(Length_ele[13])):
            lamda[2][3], lamda[3][2] = lamda[2][3]-Length_ele[13][j]/A3x3[13][j][2][2], lamda[3][2]-Length_ele[13][j]/A3x3[13][j][2][2]
        Cah=np.linalg.solve(lamda,-(Area+Beta))

    Stif=np.zeros((6,6),dtype=float)
    for i in range(0,len(Length_ele)):
        for j in range(0,len(Length_ele[i])):
            #-----------------
            if num_webs==0: Ah=(Cah-2.0*B3x3[i][j][2][2])/A3x3[i][j][2][2]
            #-----------------
            if num_webs==1:
                if i in [6,5,4]:
                    Ah=(Cah[0]-2.0*B3x3[i][j][2][2])/A3x3[i][j][2][2]
                elif i in [10,9,8,7,3,2,1,0]:
                    Ah=(Cah[1]-2.0*B3x3[i][j][2][2])/A3x3[i][j][2][2]
                else:
                    Ah=(Cah[0]-Cah[1]-2.0*B3x3[i][j][2][2])/A3x3[i][j][2][2]
            #----------------
            if num_webs==2:
                if i in [6,5,4]:
                    Ah=(Cah[0]-2.0*B3x3[i][j][2][2])/A3x3[i][j][2][2]
                elif i in [7,3]:
                    Ah=(Cah[1]-2.0*B3x3[i][j][2][2])/A3x3[i][j][2][2]
                elif i in [10,9,8,2,1,0]:
                    Ah=(Cah[2]-2.0*B3x3[i][j][2][2])/A3x3[i][j][2][2]
                elif i==11:
                    Ah=(Cah[0]-Cah[1]-2.0*B3x3[i][j][2][2])/A3x3[i][j][2][2]
                else:
                    Ah=(Cah[1]-Cah[2]-2.0*B3x3[i][j][2][2])/A3x3[i][j][2][2]
            #----------------
            if num_webs==3:
                if i in [6,5,4]:
                    Ah=(Cah[0]-2.0*B3x3[i][j][2][2])/A3x3[i][j][2][2]
                elif i in [7,3]:
                    Ah=(Cah[1]-2.0*B3x3[i][j][2][2])/A3x3[i][j][2][2]
                elif i in [9,8,2,1]:
                    Ah=(Cah[2]-2.0*B3x3[i][j][2][2])/A3x3[i][j][2][2]
                elif i in [10,0]:
                    Ah=(Cah[3]-2.0*B3x3[i][j][2][2])/A3x3[i][j][2][2]
                elif i==11:
                    Ah=(Cah[0]-Cah[1]-2.0*B3x3[i][j][2][2])/A3x3[i][j][2][2]
                elif i==12:
                    Ah=(Cah[1]-Cah[2]-2.0*B3x3[i][j][2][2])/A3x3[i][j][2][2]
                else:
                    Ah=(Cah[2]-Cah[3]-2.0*B3x3[i][j][2][2])/A3x3[i][j][2][2]

            sng=float(0.0 if (i==11 or i==12 or i==13) else 1.0)
            # matrix A:
            Stif[0][0]=Stif[0][0]+Length_ele[i][j]*float(TSF[0])*(A3x3[i][j][2][2]*(Theta_ele[i][j][0]**2.0)+A3x3[i][j][1][1]*(Theta_ele[i][j][1]**2.0))                   ; Stif[0][1]=Stif[0][1]+Length_ele[i][j]*A3x3[i][j][0][2]*Theta_ele[i][j][0] ; Stif[0][2]=Stif[0][2]+sng*3.0*Length_ele[i][j]*(A3x3[i][j][2][2]*Theta_ele[i][j][1]*Theta_ele[i][j][0]-A3x3[i][j][1][1]*Theta_ele[i][j][0]*Theta_ele[i][j][1])
            Stif[1][0]=Stif[1][0]+Length_ele[i][j]*A3x3[i][j][0][2]*Theta_ele[i][j][0]                                                                                     ; Stif[1][1]=Stif[1][1]+Length_ele[i][j]*A3x3[i][j][0][0]                    ; Stif[1][2]=Stif[1][2]+sng*Length_ele[i][j]*A3x3[i][j][0][2]*Theta_ele[i][j][1]
            Stif[2][0]=Stif[2][0]+sng*3.0*Length_ele[i][j]*(A3x3[i][j][2][2]*Theta_ele[i][j][1]*Theta_ele[i][j][0]-A3x3[i][j][1][1]*Theta_ele[i][j][0]*Theta_ele[i][j][1]) ; Stif[2][1]=Stif[2][1]+Length_ele[i][j]*A3x3[i][j][0][2]*Theta_ele[i][j][1] ; Stif[2][2]=Stif[2][2]+Length_ele[i][j]*float(TSF[1])*(A3x3[i][j][2][2]*(Theta_ele[i][j][1]**2.0)+A3x3[i][j][1][1]*(Theta_ele[i][j][0]**2.0))
            # matrix B:
            Stif[0][3]=Stif[0][3]+sng*Length_ele[i][j]*(-A3x3[i][j][0][2]*Center_ele[i][j][1]+B3x3[i][j][0][2]*Theta_ele[i][j][0])*Theta_ele[i][j][0] ; Stif[0][4]=Stif[0][4]+0.0                                                              ; Stif[0][5]=Stif[0][5]+sng*Length_ele[i][j]*(A3x3[i][j][0][2]*Center_ele[i][j][0]+B3x3[i][j][0][2]*Theta_ele[i][j][1])*Theta_ele[i][j][0]
            Stif[1][3]=Stif[1][3]+sng*Length_ele[i][j]*(-A3x3[i][j][0][0]*Center_ele[i][j][1]+B3x3[i][j][0][0]*Theta_ele[i][j][0])                    ; Stif[1][4]=Stif[1][4]+sng*Length_ele[i][j]*(-A3x3[i][j][0][2]*Ah-2.0*B3x3[i][j][0][2]) ; Stif[1][5]=Stif[1][5]+sng*Length_ele[i][j]*(A3x3[i][j][0][0]*Center_ele[i][j][0]+B3x3[i][j][0][0]*Theta_ele[i][j][1])
            Stif[2][3]=Stif[2][3]+sng*Length_ele[i][j]*(-A3x3[i][j][0][2]*Center_ele[i][j][1]+B3x3[i][j][0][2]*Theta_ele[i][j][0])*Theta_ele[i][j][1] ; Stif[2][4]=Stif[2][4]+0.0                                                              ; Stif[2][5]=Stif[2][5]+sng*Length_ele[i][j]*(A3x3[i][j][0][2]*Center_ele[i][j][0]+B3x3[i][j][0][2]*Theta_ele[i][j][1])*Theta_ele[i][j][1]
            # matrix BT:
            Stif[3][0]=Stif[3][0]+sng*Length_ele[i][j]*(-A3x3[i][j][0][2]*Center_ele[i][j][1]+B3x3[i][j][0][2]*Theta_ele[i][j][0])*Theta_ele[i][j][0] ; Stif[3][1]=Stif[3][1]+sng*Length_ele[i][j]*(-A3x3[i][j][0][0]*Center_ele[i][j][1]+B3x3[i][j][0][0]*Theta_ele[i][j][0]) ; Stif[3][2]=Stif[3][2]+sng*Length_ele[i][j]*(-A3x3[i][j][0][2]*Center_ele[i][j][1]+B3x3[i][j][0][2]*Theta_ele[i][j][0])*Theta_ele[i][j][1]
            Stif[4][0]=Stif[4][0]+0.0                                                                                                                 ; Stif[4][1]=Stif[4][1]+sng*Length_ele[i][j]*(-A3x3[i][j][0][2]*Ah-2.0*B3x3[i][j][0][2])                                 ; Stif[4][2]=Stif[4][2]+0.0
            Stif[5][0]=Stif[5][0]+sng*Length_ele[i][j]*(A3x3[i][j][0][2]*Center_ele[i][j][0]+B3x3[i][j][0][2]*Theta_ele[i][j][1])*Theta_ele[i][j][0]  ; Stif[5][1]=Stif[5][1]+sng*Length_ele[i][j]*(A3x3[i][j][0][0]*Center_ele[i][j][0]+B3x3[i][j][0][0]*Theta_ele[i][j][1])  ; Stif[5][2]=Stif[5][2]+sng*Length_ele[i][j]*(A3x3[i][j][0][2]*Center_ele[i][j][0]+B3x3[i][j][0][2]*Theta_ele[i][j][1])*Theta_ele[i][j][1]
            # matrix D:
            Stif[3][3]=Stif[3][3]+Length_ele[i][j]*(A3x3[i][j][0][0]*(Center_ele[i][j][1]**2.0)-2.0*B3x3[i][j][0][0]*Center_ele[i][j][1]*Theta_ele[i][j][0]+D3x3[i][j][0][0]*(Theta_ele[i][j][0]**2.0))                                                                ; Stif[3][4]=Stif[3][4]-Length_ele[i][j]*(A3x3[i][j][0][2]*Ah*Center_ele[i][j][1]-B3x3[i][j][0][2]*(Ah*Theta_ele[i][j][0]-2.0*Center_ele[i][j][0])-2.0*D3x3[i][j][0][2]*Theta_ele[i][j][0])  ; Stif[3][5]=Stif[3][5]+Length_ele[i][j]*(-A3x3[i][j][0][0]*Center_ele[i][j][1]*Center_ele[i][j][0]-B3x3[i][j][0][0]*(Center_ele[i][j][1]*Theta_ele[i][j][1]-Center_ele[i][j][0]*Theta_ele[i][j][0])+D3x3[i][j][0][0]*Theta_ele[i][j][1]*Theta_ele[i][j][0])
            Stif[4][3]=Stif[4][3]-Length_ele[i][j]*(A3x3[i][j][0][2]*Ah*Center_ele[i][j][1]-B3x3[i][j][0][2]*(Ah*Theta_ele[i][j][0]-2.0*Center_ele[i][j][0])-2.0*D3x3[i][j][0][2]*Theta_ele[i][j][0])                                                                  ; Stif[4][4]=Stif[4][4]+Length_ele[i][j]*(A3x3[i][j][2][2]*(Ah**2.0)+4.0*B3x3[i][j][2][2]*Ah+4.0*D3x3[i][j][2][2])                                                                           ; Stif[4][5]=Stif[4][5]-Length_ele[i][j]*(-A3x3[i][j][0][2]*Ah*Center_ele[i][j][0]-B3x3[i][j][0][2]*(Ah*Theta_ele[i][j][1]+2.0*Center_ele[i][j][0])-2.0*D3x3[i][j][0][2]*Theta_ele[i][j][1])
            Stif[5][3]=Stif[5][3]+Length_ele[i][j]*(-A3x3[i][j][0][0]*Center_ele[i][j][1]*Center_ele[i][j][0]-B3x3[i][j][0][0]*(Center_ele[i][j][1]*Theta_ele[i][j][1]-Center_ele[i][j][0]*Theta_ele[i][j][0])+D3x3[i][j][0][0]*Theta_ele[i][j][1]*Theta_ele[i][j][0]) ; Stif[5][4]=Stif[5][4]-Length_ele[i][j]*(-A3x3[i][j][0][2]*Ah*Center_ele[i][j][0]-B3x3[i][j][0][2]*(Ah*Theta_ele[i][j][1]+2.0*Center_ele[i][j][0])-2.0*D3x3[i][j][0][2]*Theta_ele[i][j][1]) ; Stif[5][5]=Stif[5][5]+Length_ele[i][j]*(A3x3[i][j][0][0]*(Center_ele[i][j][0]**2.0)+2.0*B3x3[i][j][0][0]*Center_ele[i][j][0]*Theta_ele[i][j][1]+D3x3[i][j][0][0]*(Theta_ele[i][j][1]**2.0))

    Stif[0][4]=sc[1]*Stif[0][0]-sc[0]*Stif[0][2]            ; Stif[4][0]=sc[1]*Stif[0][0]-sc[0]*Stif[0][2]
    Stif[1][4]=sc[1]*Stif[0][1]-sc[0]*Stif[2][1]+Stif[1][4] ; Stif[4][1]=sc[1]*Stif[0][1]-sc[0]*Stif[2][1]+Stif[4][1] 
    Stif[2][4]=sc[1]*Stif[0][2]-sc[0]*Stif[2][2]            ; Stif[4][2]=sc[1]*Stif[0][2]-sc[0]*Stif[2][2]
    Stif[3][4]=sc[1]*Stif[0][3]-sc[0]*Stif[2][3]+Stif[3][4] ; Stif[4][3]=sc[1]*Stif[0][3]-sc[0]*Stif[2][3]+Stif[4][3]
    Stif[4][5]=sc[1]*Stif[0][5]-sc[0]*Stif[2][5]+Stif[4][5] ; Stif[5][4]=sc[1]*Stif[0][5]-sc[0]*Stif[2][5]+Stif[5][4]
    Stif[4][4]=sc[1]*Stif[0][3]-sc[0]*Stif[2][3]+Stif[4][4]

    Q6x6=inv(Stif)
    if len(F)==0:
       Gx=1.0 ; Gy=1.0
    else:
       Gx=-F[2] ; Gy=F[0]
    qx=[] ; qy=[] ; q=[]
    for i in range(0,len(Length_ele)):
        qx.append([0.0]) if i==0 else qx.append([qx[i-1][-1]])
        qy.append([0.0]) if i==0 else qy.append([qy[i-1][-1]])
        q.append([0.0]) if i==0 else q.append([q[i-1][-1]])
        for j in range(1,len(Length_ele[i])):
            qx[i].append(qx[i][j-1]+Length_ele[i][j]*(A3x3[i][j][0][0]*(Q6x6[1][3]-Q6x6[3][3]*Center_ele[i][j][1]+Q6x6[5][3]*Center_ele[i][j][0])+B3x3[i][j][0][0]*(Q6x6[3][3]*Theta_ele[i][j][0]+Q6x6[5][3]*Theta_ele[i][j][1])+A3x3[i][j][0][2]*(Q6x6[2][3]*Theta_ele[i][j][1]+Q6x6[0][3]*Theta_ele[i][j][0])-2.0*B3x3[i][j][0][2]*Q6x6[4][3]))
            qy[i].append(qy[i][j-1]+Length_ele[i][j]*(A3x3[i][j][0][0]*(Q6x6[1][5]-Q6x6[3][5]*Center_ele[i][j][1]+Q6x6[5][5]*Center_ele[i][j][0])+B3x3[i][j][0][0]*(Q6x6[3][5]*Theta_ele[i][j][0]+Q6x6[5][5]*Theta_ele[i][j][1])+A3x3[i][j][0][2]*(Q6x6[2][5]*Theta_ele[i][j][1]+Q6x6[0][5]*Theta_ele[i][j][0])-2.0*B3x3[i][j][0][2]*Q6x6[5][4]))
            q[i].append(q[i][j-1]+Gx*Length_ele[i][j]*(A3x3[i][j][0][0]*(Q6x6[1][3]-Q6x6[3][3]*Center_ele[i][j][1]+Q6x6[5][3]*Center_ele[i][j][0])+B3x3[i][j][0][0]*(Q6x6[3][3]*Theta_ele[i][j][0]+Q6x6[5][3]*Theta_ele[i][j][1])+A3x3[i][j][0][2]*(Q6x6[2][3]*Theta_ele[i][j][1]+Q6x6[0][3]*Theta_ele[i][j][0])-2.0*B3x3[i][j][0][2]*Q6x6[4][3])+Gy*Length_ele[i][j]*(A3x3[i][j][0][0]*(Q6x6[1][5]-Q6x6[3][5]*Center_ele[i][j][1]+Q6x6[5][5]*Center_ele[i][j][0])+B3x3[i][j][0][0]*(Q6x6[3][5]*Theta_ele[i][j][0]+Q6x6[5][5]*Theta_ele[i][j][1])+A3x3[i][j][0][2]*(Q6x6[2][5]*Theta_ele[i][j][1]+Q6x6[0][5]*Theta_ele[i][j][0])-2.0*B3x3[i][j][0][2]*Q6x6[5][4]))
    #----------------------------
    if num_webs==0:
        Alphax=0.0 ; Alphay=0.0 ; Alpha=0.0 ; Beta=0.0
        for i in range(0,len(Length_ele)):
            for j in range(0,len(Length_ele[i])):
                Alphax, Alphay, Alpha, Beta = Alphax+Length_ele[i][j]*qx[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0]), Alphay+Length_ele[i][j]*qy[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0]), Alpha+Length_ele[i][j]*q[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0]), Beta+Length_ele[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0])
        ratex=Alphax/Beta ; ratey=Alphay/Beta ; rate=Alpha/Beta
        Qx=[] ; Qy=[] ; Q=[]
        for i in range(0,len(qx)):
            Qx.append([]) ; Qy.append([]) ; Q.append([])
            for j in range(0,len(qx[i])):
                Qx[i].append(qx[i][j]-ratex) ; Qy[i].append(qy[i][j]-ratey) ; Q[i].append(q[i][j]-rate)
    #----------------------------
    if num_webs==1:
        Alphax=np.array([0.0, 0.0]) ; Alphay=np.array([0.0, 0.0]) ; Alpha=np.array([0.0, 0.0]) ; Beta=np.zeros((2,2),dtype=float)
        for i in [6,5,4,11]:
            for j in range(0,len(Length_ele[i])):
                Alphax[0], Alphay[0], Alpha[0], Beta[0][0] = Alphax[0]+Length_ele[i][j]*qx[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0]), Alphay[0]+Length_ele[i][j]*qy[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0]), Alpha[0]+Length_ele[i][j]*q[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0]), Beta[0][0]+Length_ele[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0])
        for i in [10,9,8,7,11,3,2,1,0]:
            for j in range(0,len(Length_ele[i])):
                if i==11:
                    Alphax[1], Alphay[1], Alpha[1] = Alphax[1]-Length_ele[i][j]*qx[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0]), Alphay[1]-Length_ele[i][j]*qy[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0]), Alpha[1]-Length_ele[i][j]*q[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0])
                else:
                    Alphax[1], Alphay[1], Alpha[1] = Alphax[1]+Length_ele[i][j]*qx[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0]), Alphay[1]+Length_ele[i][j]*qy[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0]), Alpha[1]+Length_ele[i][j]*q[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0])
                Beta[1][1]=Beta[1][1]+Length_ele[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0])
        for j in range(0,len(Length_ele[11])):
            Beta[0][1], Beta[1][0] = Beta[0][1]-Length_ele[11][j]/(A3x3[11][j][2][2]-(A3x3[11][j][0][2]**2.0)/A3x3[11][j][0][0]), Beta[1][0]-Length_ele[11][j]/(A3x3[11][j][2][2]-(A3x3[11][j][0][2]**2.0)/A3x3[11][j][0][0])
        ratex, ratey, rate = np.linalg.solve(Beta,Alphax), np.linalg.solve(Beta,Alphay), np.linalg.solve(Beta,Alpha)
        Qx=[] ; Qy=[] ; Q=[]
        for i in range(0,len(qx)):
            Qx.append([]) ; Qy.append([]) ; Q.append([])
            for j in range(0,len(qx[i])):
                if i in [6,5,4]:
                    Qx[i].append(qx[i][j]-ratex[0]) ; Qy[i].append(qy[i][j]-ratey[0]) ; Q[i].append(q[i][j]-rate[0])
                elif i in [10,9,8,7,3,2,1,0]:
                    Qx[i].append(qx[i][j]-ratex[1]) ; Qy[i].append(qy[i][j]-ratey[1]) ; Q[i].append(q[i][j]-rate[1])
                else:
                    Qx[i].append(qx[i][j]-(ratex[0]-ratex[1])) ; Qy[i].append(qy[i][j]-(ratey[0]-ratey[1])) ; Q[i].append(q[i][j]-(rate[0]-rate[1]))
    #-----------------------------
    if num_webs==2:
        Alphax=np.array([0.0, 0.0, 0.0]) ; Alphay=np.array([0.0, 0.0, 0.0]) ; Alpha=np.array([0.0, 0.0, 0.0]) ; Beta=np.zeros((3,3),dtype=float)
        for i in [6,5,4,11]:
            for j in range(0,len(Length_ele[i])):
                Alphax[0], Alphay[0], Alpha[0], Beta[0][0] = Alphax[0]+Length_ele[i][j]*qx[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0]), Alphay[0]+Length_ele[i][j]*qy[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0]), Alpha[0]+Length_ele[i][j]*q[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0]), Beta[0][0]+Length_ele[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0])
        for i in [7, 11, 3, 12]:
            for j in range(0,len(Length_ele[i])):
                if i==11 or i==12:
                    Alphax[1], Alphay[1], Alpha[1] = Alphax[1]-Length_ele[i][j]*qx[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0]), Alphay[1]-Length_ele[i][j]*qy[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0]), Alpha[1]-Length_ele[i][j]*q[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0])
                else:
                    Alphax[1], Alphay[1], Alpha[1] = Alphax[1]+Length_ele[i][j]*qx[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0]), Alphay[1]+Length_ele[i][j]*qy[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0]), Alpha[1]+Length_ele[i][j]*q[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0])
                Beta[1][1]=Beta[1][1]+Length_ele[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0])
        for i in [10, 9, 8, 12, 2, 1, 0]:
            for j in range(0,len(Length_ele[i])):
                Alphax[2], Alphay[2], Alpha[2], Beta[2][2] = Alphax[2]+Length_ele[i][j]*qx[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0]), Alphay[2]+Length_ele[i][j]*qy[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0]), Alpha[2]+Length_ele[i][j]*q[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0]), Beta[2][2]+Length_ele[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0])
        for j in range(0,len(Length_ele[11])):
            Beta[0][1], Beta[1][0] = Beta[0][1]-Length_ele[11][j]/(A3x3[11][j][2][2]-(A3x3[11][j][0][2]**2.0)/A3x3[11][j][0][0]), Beta[1][0]-Length_ele[11][j]/(A3x3[11][j][2][2]-(A3x3[11][j][0][2]**2.0)/A3x3[11][j][0][0])
        for j in range(0,len(Length_ele[12])):
            Beta[1][2], Beta[2][1] = Beta[1][2]-Length_ele[12][j]/(A3x3[12][j][2][2]-(A3x3[12][j][0][2]**2.0)/A3x3[12][j][0][0]), Beta[2][1]-Length_ele[12][j]/(A3x3[12][j][2][2]-(A3x3[12][j][0][2]**2.0)/A3x3[12][j][0][0])
        ratex, ratey, rate = np.linalg.solve(Beta,Alphax), np.linalg.solve(Beta,Alphay), np.linalg.solve(Beta,Alpha)

        Qx=[] ; Qy=[] ; Q=[]
        for i in range(0,len(qx)):
            Qx.append([]) ; Qy.append([]) ; Q.append([])
            for j in range(0,len(qx[i])):
                if i in [6,5,4]:
                    Qx[i].append(qx[i][j]-ratex[0]) ; Qy[i].append(qy[i][j]-ratey[0]) ; Q[i].append(q[i][j]-rate[0])
                elif i in [7,3]:
                    Qx[i].append(qx[i][j]-ratex[1]) ; Qy[i].append(qy[i][j]-ratey[1]) ; Q[i].append(q[i][j]-rate[1])
                elif i in [10,9,8,2,1,0]:
                    Qx[i].append(qx[i][j]-ratex[2]) ; Qy[i].append(qy[i][j]-ratey[2]) ; Q[i].append(q[i][j]-rate[2])
                elif i==11:
                    Qx[i].append(qx[i][j]-(ratex[0]-ratex[1])) ; Qy[i].append(qy[i][j]-(ratey[0]-ratey[1])) ; Q[i].append(q[i][j]-(rate[0]-rate[1]))
                else:
                    Qx[i].append(qx[i][j]+(ratex[1]-ratex[2])) ; Qy[i].append(qy[i][j]+(ratey[1]-ratey[2])) ; Q[i].append(q[i][j]+(rate[1]-rate[2]))
    #-----------------------------
    if num_webs==3:
        Alphax=np.array([0.0, 0.0, 0.0, 0.0]) ; Alphay=np.array([0.0, 0.0, 0.0, 0.0]) ; Alpha=np.array([0.0, 0.0, 0.0, 0.0]) ; Beta=np.zeros((4,4),dtype=float)
        for i in [6,5,4,11]:
            for j in range(0,len(Length_ele[i])):
                Alphax[0], Alphay[0], Alpha[0], Beta[0][0] = Alphax[0]+Length_ele[i][j]*qx[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0]), Alphay[0]+Length_ele[i][j]*qy[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0]), Alpha[0]+Length_ele[i][j]*q[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0]), Beta[0][0]+Length_ele[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0])
        for i in [7, 11, 3, 12]:
            for j in range(0,len(Length_ele[i])):
                if i==11 or i==12:
                    Alphax[1], Alphay[1], Alpha[1] = Alphax[1]-Length_ele[i][j]*qx[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0]), Alphay[1]-Length_ele[i][j]*qy[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0]), Alpha[1]-Length_ele[i][j]*q[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0])
                else:
                    Alphax[1], Alphay[1], Alpha[1] = Alphax[1]+Length_ele[i][j]*qx[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0]), Alphay[1]+Length_ele[i][j]*qy[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0]), Alpha[1]+Length_ele[i][j]*q[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0])
                Beta[1][1]=Beta[1][1]+Length_ele[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0])
        for i in [9, 8, 12, 2, 1, 13]:
            for j in range(0,len(Length_ele[i])):
                if i==13:
                    Alphax[2], Alphay[2], Alpha[2] = Alphax[2]-Length_ele[i][j]*qx[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0]), Alphay[2]-Length_ele[i][j]*qy[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0]), Alpha[2]-Length_ele[i][j]*q[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0])
                else:
                    Alphax[2], Alphay[2], Alpha[2] = Alphax[2]+Length_ele[i][j]*qx[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0]), Alphay[2]+Length_ele[i][j]*qy[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0]), Alpha[2]+Length_ele[i][j]*q[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0])
                Beta[2][2]=Beta[2][2]+Length_ele[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0])
        for i in [10,13,0]:
            for j in range(0,len(Length_ele[i])):
                Alphax[3], Alphay[3], Alpha[3], Beta[3][3] = Alphax[3]+Length_ele[i][j]*qx[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0]), Alphay[3]+Length_ele[i][j]*qy[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0]), Alpha[3]+Length_ele[i][j]*q[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0]), Beta[3][3]+Length_ele[i][j]/(A3x3[i][j][2][2]-(A3x3[i][j][0][2]**2.0)/A3x3[i][j][0][0])
        for j in range(0,len(Length_ele[11])):
            Beta[0][1],Beta[1][0] = Beta[0][1]-Length_ele[11][j]/(A3x3[11][j][2][2]-(A3x3[11][j][0][2]**2.0)/A3x3[11][j][0][0]), Beta[1][0]-Length_ele[11][j]/(A3x3[11][j][2][2]-(A3x3[11][j][0][2]**2.0)/A3x3[11][j][0][0])
        for j in range(0,len(Length_ele[12])):
            Beta[1][2], Beta[2][1] = Beta[1][2]-Length_ele[12][j]/(A3x3[12][j][2][2]-(A3x3[12][j][0][2]**2.0)/A3x3[12][j][0][0]), Beta[2][1]-Length_ele[12][j]/(A3x3[12][j][2][2]-(A3x3[12][j][0][2]**2.0)/A3x3[12][j][0][0])
        for j in range(0,len(Length_ele[13])):
            Beta[2][3], Beta[3][2] = Beta[2][3]-Length_ele[13][j]/(A3x3[13][j][2][2]-(A3x3[13][j][0][2]**2.0)/A3x3[13][j][0][0]), Beta[3][2]-Length_ele[13][j]/(A3x3[13][j][2][2]-(A3x3[13][j][0][2]**2.0)/A3x3[13][j][0][0])
        ratex, ratey, rate = np.linalg.solve(Beta,Alphax), np.linalg.solve(Beta,Alphay), np.linalg.solve(Beta,Alpha)

        Qx=[] ; Qy=[] ; Q=[]
        for i in range(0,len(qx)):
            Qx.append([]) ; Qy.append([]) ; Q.append([])
            for j in range(0,len(qx[i])):
                if i in [6,5,4]:
                    Qx[i].append(qx[i][j]-ratex[0]) ; Qy[i].append(qy[i][j]-ratey[0]) ; Q[i].append(q[i][j]-rate[0])
                elif i in [7,3]:
                    Qx[i].append(qx[i][j]-ratex[1]) ; Qy[i].append(qy[i][j]-ratey[1]) ; Q[i].append(q[i][j]-rate[1])
                elif i in [9,8,2,1]:
                    Qx[i].append(qx[i][j]-ratex[2]) ; Qy[i].append(qy[i][j]-ratey[2]) ; Q[i].append(q[i][j]-rate[2])
                elif i in [10,0]:
                    Qx[i].append(qx[i][j]-ratex[3]) ; Qy[i].append(qy[i][j]-ratey[3]) ; Q[i].append(q[i][j]-rate[3])
                elif i==11:
                    Qx[i].append(qx[i][j]-(ratex[0]-ratex[1])) ; Qy[i].append(qy[i][j]-(ratey[0]-ratey[1])) ; Q[i].append(q[i][j]-(rate[0]-rate[1]))
                elif i==12:
                    Qx[i].append(qx[i][j]+(ratex[1]-ratex[2])) ; Qy[i].append(qy[i][j]+(ratey[1]-ratey[2])) ; Q[i].append(q[i][j]+(rate[1]-rate[2]))
                else:
                    Qx[i].append(qx[i][j]+(ratex[2]-ratex[3])) ; Qy[i].append(qy[i][j]+(ratey[2]-ratey[3])) ; Q[i].append(q[i][j]+(rate[2]-rate[3]))
    if len(F)==0:
        Tx=0.0 ; Ty=0.0
        for i in range(0,len(Length_ele)):
            for j in range(0,len(Length_ele[i])):
                Tx, Ty = Tx+Length_ele[i][j]*Qx[i][j]*(-Center_ele[i][j][1]*Theta_ele[i][j][0]+Center_ele[i][j][0]*Theta_ele[i][j][1]), Ty+Length_ele[i][j]*Qy[i][j]*(-Center_ele[i][j][1]*Theta_ele[i][j][0]+Center_ele[i][j][0]*Theta_ele[i][j][1])

        return Stif, Tx, Ty
    else:
        Load=np.array(F)
        R6x1=np.matmul(Q6x6,Load.transpose())
        s_n=[] ; s_t=[] ; s_p=[]
        for i in range(0,len(Thick_lay)):
            s_n.append([]) ; s_t.append([]) ; s_p.append([])
            for j in range(0,len(Thick_lay[i])):
                s_n[i].append([]) ; s_t[i].append([]) ; s_p[i].append([])
                if num_webs==0: Ah=(Cah-2.0*B3x3[i][j][2][2])/A3x3[i][j][2][2]
                #-----------------
                if num_webs==1:
                    if i in [6,5,4]:
                        Ah=(Cah[0]-2.0*B3x3[i][j][2][2])/A3x3[i][j][2][2]
                    elif i in [10,9,8,7,3,2,1,0]:
                        Ah=(Cah[1]-2.0*B3x3[i][j][2][2])/A3x3[i][j][2][2]
                    else:
                        Ah=(Cah[0]-Cah[1]-2.0*B3x3[i][j][2][2])/A3x3[i][j][2][2]
                #----------------
                if num_webs==2:
                    if i in [6,5,4]:
                        Ah=(Cah[0]-2.0*B3x3[i][j][2][2])/A3x3[i][j][2][2]
                    elif i in [7,3]:
                        Ah=(Cah[1]-2.0*B3x3[i][j][2][2])/A3x3[i][j][2][2]
                    elif i in [10,9,8,2,1,0]:
                        Ah=(Cah[2]-2.0*B3x3[i][j][2][2])/A3x3[i][j][2][2]
                    elif i==11:
                        Ah=(Cah[0]-Cah[1]-2.0*B3x3[i][j][2][2])/A3x3[i][j][2][2]
                if num_webs==3:
                    if i in [6,5,4]:
                        Ah=(Cah[0]-2.0*B3x3[i][j][2][2])/A3x3[i][j][2][2]
                    elif i in [7,3]:
                        Ah=(Cah[1]-2.0*B3x3[i][j][2][2])/A3x3[i][j][2][2]
                    elif i in [9,8,2,1]:
                        Ah=(Cah[2]-2.0*B3x3[i][j][2][2])/A3x3[i][j][2][2]
                    elif i in [10,0]:
                        Ah=(Cah[3]-2.0*B3x3[i][j][2][2])/A3x3[i][j][2][2]
                    elif i==11:
                        Ah=(Cah[0]-Cah[1]-2.0*B3x3[i][j][2][2])/A3x3[i][j][2][2]
                    elif i==12:
                        Ah=(Cah[1]-Cah[2]-2.0*B3x3[i][j][2][2])/A3x3[i][j][2][2]
                    else:
                        Ah=(Cah[2]-Cah[3]-2.0*B3x3[i][j][2][2])/A3x3[i][j][2][2]

                for k in range(0,len(Thick_lay[i][j])):
                    if Thick_lay[i][j][k]>0.0:
                        K=k
                        break
                for k in range(0,len(Thick_lay[i][j])):
                    if Thick_lay[i][j][k]>0.0:
                        z= (Thick_lay[i][j][k]-sum(Thick_lay[i][j][:]))/2.0 if k==K else z+(thick_lay+Thick_lay[i][j][k])/2.0

                    if Thick_lay[i][j][k]>0.0 and k!=3:
                        s_n[i][j].append((R6x1[1]-R6x1[3]*(Center_ele[i][j][1]-Theta_ele[i][j][0]*z)+R6x1[5]*(Center_ele[i][j][0]+Theta_ele[i][j][1]*z))*q3x3[i][j][k][0][0]+((Q[i][j]-2.0*R6x1[4]*B3x3[i][j][2][2])/A3x3[i][j][2][2]-R6x1[4]*(Ah+2.0*z))*q3x3[i][j][k][0][2])
                        s_p[i][j].append(-(R6x1[1]-R6x1[3]*(Center_ele[i][j][1]-Theta_ele[i][j][0]*z)+R6x1[5]*(Center_ele[i][j][0]+Theta_ele[i][j][1]*z))*q3x3[i][j][k][0][2]-((Q[i][j]-2.0*R6x1[4]*B3x3[i][j][2][2])/A3x3[i][j][2][2]-R6x1[4]*(Ah+2.0*z))*q3x3[i][j][k][2][2])
                        s_t[i][j].append((-R6x1[1]*Theta_ele[i][j][0]+R6x1[0]*Theta_ele[i][j][1])*q3x3[i][j][k][1][1])
                    else:
                        s_n[i][j].append(0.0)
                        s_p[i][j].append(0.0)
                        s_t[i][j].append(0.0)
        Str=[]
        for i in range(0,len(s_n)):
            Str.append([])
            for j in range(0,len(s_n[i])):
                Str[i].append([])
                for k in range(0,len(s_n[i][j])):
                    theta =  -abs(Ply_lay[i][j][k])*math.pi/180.0
                    X = xt[k] if s_n[i][j][k]>=0.0 else xc[k]
                    Y = yt[k] if s_p[i][j][k]>=0.0 else yc[k]
                    S = st[k]
                    STR = np.matmul(np.matmul(np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]]), np.array([[X, S], [S, Y]])), np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]]).transpose())
                    Str[i][j].append((s_n[i][j][k]/STR[0][0])**2.0-abs(s_n[i][j][k]*s_p[i][j][k])/(STR[0][0]**2.0)+(s_p[i][j][k]/STR[1][1])**2.0)
        StrRange=[]
        for i in range(0,11):
            Str_range_max=0.0
            for j in range(0,len(Str[i])):
                for k in range(0,len(Str[i][j])):
                    if Str[i][j][k]>Str_range_max: Str_range_max=Str[i][j][k]
            StrRange.append(Str_range_max)
        return Str,StrRange,max(StrRange)

# Shear Center
#################################################
def shear_center(sc,*arg):
    Thick_lay, Length_lay, Ply_lay, Thick_ele, Length_ele, Theta_ele, Center_ele, triax, biax, uniax, balsa, num_webs, TimFac = arg
    Center_Ele=[] ; [Center_Ele.append([(Center_ele[i][j][0]-sc[0], Center_ele[i][j][1]-sc[1]) for j in range(0,len(Center_ele[i]))]) for i in range(0,len(Center_ele))]
    ukn, Tx, Ty = stiffness(Thick_lay, Length_lay, Ply_lay, Thick_ele, Length_ele, Theta_ele, Center_Ele, triax, biax, uniax, balsa, sc, num_webs, TimFac, [])
    return Tx, Ty

def STIFFMASS(args):
    profil,radius,chord,twist,pit_ax,kp01,kp02,kp03,kp04,kp05,kp06,kp07,kp08,kp09,kp10,laminates,num_webs,ply,triax,biax,uniax,balsa = args
    K6x6,SC,M6x6,MC = [],[],[],[]
    for i in range(0,21):
        airfoil,kp = airfoil_geometry(profil,radius[i],chord[i],twist[i],pit_ax[i],kp01[i],kp02[i],kp03[i],kp04[i],kp05[i],kp06[i],kp07[i],kp08[i],kp09[i],kp10[i])
        Thick_lay,Length_lay,Ply_lay,Thick_ele,Length_ele,Theta_ele,Center_ele,X_Grid,Y_Grid,C_Grid = grid(airfoil,kp,laminates[i],num_webs[i],ply[i])
        m6x6,mc = mass(Thick_lay,Thick_ele,Length_ele,Center_ele,triax,biax,uniax,balsa)
        Timo=True ; sc=tuple(fsolve(shear_center,[0.0,0.0],args=(Thick_lay,Length_lay,Ply_lay,Thick_ele,Length_ele,Theta_ele,Center_ele,triax,biax,uniax,balsa,num_webs[i],Timo),xtol=0.1))
        k6x6,ukn,ukn = stiffness(Thick_lay,Length_lay,Ply_lay,Thick_ele,Length_ele,Theta_ele,Center_ele,triax,biax,uniax,balsa,sc,num_webs[i],Timo,[])
        K6x6.append(k6x6) ; SC.append(sc) ; M6x6.append(m6x6) ; MC.append(mc)
    return K6x6,SC,M6x6,MC

def STRESSES(args):
    R_tip,r_hup,profil,radius,x_sweep,z_prebent,chord,twist,pit_ax,kp01,kp02,kp03,kp04,kp05,kp06,kp07,kp08,kp09,kp10,laminates,num_webs,ply,triax,biax,uniax,balsa,Loads = args
    StrMax,Length,COST = [],[],[]
    for i in range(0,21):
        airfoil,kp = airfoil_geometry(profil,radius[i],chord[i],twist[i],pit_ax[i],kp01[i],kp02[i],kp03[i],kp04[i],kp05[i],kp06[i],kp07[i],kp08[i],kp09[i],kp10[i])
        Thick_lay,Length_lay,Ply_lay,Thick_ele,Length_ele,Theta_ele,Center_ele,X_Grid,Y_Grid,C_Grid = grid(airfoil,kp,laminates[i],num_webs[i],ply[i])
        Length.append(math.sqrt(((R_tip-r_hup)*radius[i])**2.0+x_sweep[i]**2.0+z_prebent[i]**2.0))
        COST.append(sum([Thick_lay[i][j][0]*Length_lay[i][j][0]*triax[0]*triax[1]+Thick_lay[i][j][1]*Length_lay[i][j][1]*biax[0]*biax[1]+Thick_lay[i][j][2]*Length_lay[i][j][2]*uniax[0]*uniax[1]+Thick_lay[i][j][3]*Length_lay[i][j][3]*balsa[0]*balsa[1]+Thick_lay[i][j][4]*Length_lay[i][j][4]*uniax[0]*uniax[1]+Thick_lay[i][j][5]*Length_lay[i][j][5]*biax[0]*biax[1]+Thick_lay[i][j][6]*Length_lay[i][j][6]*triax[0]*triax[1] for i in range(0,len(Thick_lay)) for j in range(0,len(Thick_lay[i]))]))
        Timo=True ; sc=tuple(fsolve(shear_center,[0.0,0.0],args=(Thick_lay,Length_lay,Ply_lay,Thick_ele,Length_ele,Theta_ele,Center_ele,triax,biax,uniax,balsa,num_webs[i],Timo),xtol=0.1))
        ukn,ukn,str_max = stiffness(Thick_lay,Length_lay,Ply_lay,Thick_ele,Length_ele,Theta_ele,Center_ele,triax,biax,uniax,balsa,sc,num_webs[i],Timo,Loads[i])
        StrMax.append(str_max)
    return StrMax#,sum([(COST[i-1]+COST[i])*(Length[i]-Length[i-1])/2.0 for i in range(1,21)])

def COST(args):
    R_tip,r_hup,profil,radius,x_sweep,z_prebent,chord,twist,pit_ax,kp01,kp02,kp03,kp04,kp05,kp06,kp07,kp08,kp09,kp10,laminates,num_webs,ply,triax,biax,uniax,balsa,CosMod_params = args
    Cost_Labor = CosMod_params[8][0]*(R_tip**CosMod_params[8][1])
    Cost_Light = CosMod_params[4]*R_tip
    Cost_TbaNu = CosMod_params[3]*round(math.pi*chord[0]/0.15,0)
    VolTRIAX,VolBIAX,VolUNIAX,VolBALSA,VolRes,Cadh,Lskin,Lwebs,Lblade,CSmass = [],[],[],[],[],[],[],[],[],[]
    for i in range(0,21):
        airfoil,kp = airfoil_geometry(profil,radius[i],chord[i],twist[i],pit_ax[i],kp01[i],kp02[i],kp03[i],kp04[i],kp05[i],kp06[i],kp07[i],kp08[i],kp09[i],kp10[i])
        Thick_lay,Length_lay,Ply_lay,Thick_ele,Length_ele,Theta_ele,Center_ele,X_Grid,Y_Grid,C_Grid = grid(airfoil,kp,laminates[i],num_webs[i],ply[i])
        m6x6,ukn = mass(Thick_lay,Thick_ele,Length_ele,Center_ele,triax,biax,uniax,balsa)
        CSmass.append(m6x6[0][0])
        VolTRIAX.append(sum([Thick_lay[j][k][0]*Length_lay[j][k][0]+Thick_lay[j][k][6]*Length_lay[j][k][6] for j in range(0,len(Thick_lay)) for k in range(0,len(Thick_lay[j]))]))
        VolBIAX.append( sum([Thick_lay[j][k][1]*Length_lay[j][k][1]+Thick_lay[j][k][5]*Length_lay[j][k][5] for j in range(0,len(Thick_lay)) for k in range(0,len(Thick_lay[j]))]))
        VolUNIAX.append(sum([Thick_lay[j][k][2]*Length_lay[j][k][2]+Thick_lay[j][k][4]*Length_lay[j][k][4] for j in range(0,len(Thick_lay)) for k in range(0,len(Thick_lay[j]))]))
        VolBALSA.append(sum([Thick_lay[j][k][3]*Length_lay[j][k][3] for j in range(0,len(Thick_lay)) for k in range(0,len(Thick_lay[j]))]))
        VolRes.append(sum([(1.0-0.725)*Thick_lay[j][k][0]*Length_lay[j][k][0]+(1.0-0.725)*Thick_lay[j][k][1]*Length_lay[j][k][1]+(1.0-0.725)*Thick_lay[j][k][2]*Length_lay[j][k][2]+(1.0-0.725)*Thick_lay[j][k][4]*Length_lay[j][k][4]+(1.0-0.725)*Thick_lay[j][k][5]*Length_lay[j][k][5]+(1.0-0.725)*Thick_lay[j][k][6]*Length_lay[j][k][6] for j in range(0,len(Thick_lay)) for k in range(0,len(Thick_lay[j]))]))
        if num_webs[i]==0: Cadh.append(4.0*0.070*0.007*1200.0*CosMod_params[1])
        if num_webs[i]==1: Cadh.append(8.0*0.070*0.007*1200.0*CosMod_params[1])
        if num_webs[i]==2: Cadh.append(12.0*0.070*0.007*1200.0*CosMod_params[1])
        if num_webs[i]==3: Cadh.append(16.0*0.070*0.007*1200.0*CosMod_params[1])
        Lblade.append(math.sqrt(((R_tip-r_hup)*radius[i])**2.0+x_sweep[i]**2.0+z_prebent[i]**2.0))
        Lskin.append(sum([Length_ele[j][k] for j in range(0,11) for k in range(0,len(Length_ele[j]))]))
        Lwebs.append(sum([Length_ele[j][k] for j in range(11,len(Length_ele)) for k in range(0,len(Length_ele[j]))]))
    Mass = sum([(CSmass[i-1]+CSmass[i])*(Lblade[i]-Lblade[i-1])/2.0 for i in range(1,21)])
    Askin = sum([(Lskin[i-1]+Lskin[i])*(Lblade[i]-Lblade[i-1])/2.0 for i in range(1,21)])
    Awebs = sum([(Lwebs[i-1]+Lwebs[i])*(Lblade[i]-Lblade[i-1])/2.0 for i in range(1,21)])
    Cost_TRIAX = 0.350*triax[0]*triax[1]*sum([(VolTRIAX[i-1]+VolTRIAX[i])*(Lblade[i]-Lblade[i-1])/2.0 for i in range(1,21)])/(1.0-0.12)
    Cost_BIAX  = 0.750*biax[0] *biax[1]* sum([(VolBIAX[i-1] +VolBIAX[i] )*(Lblade[i]-Lblade[i-1])/2.0 for i in range(1,21)])/(1.0-0.12)
    Cost_UNIAX = 0.350*uniax[0]*uniax[1]*sum([(VolUNIAX[i-1]+VolUNIAX[i])*(Lblade[i]-Lblade[i-1])/2.0 for i in range(1,21)])/(1.0-0.12)
    Cost_BALSA = 1.000*balsa[0]*balsa[1]*sum([(VolBALSA[i-1]+VolBALSA[i])*(Lblade[i]-Lblade[i-1])/2.0 for i in range(1,21)])/(1.0-0.30)
    Cost_Resin = CosMod_params[0]*1200.0*sum([(VolRes[i-1]+VolRes[i])*(Lblade[i]-Lblade[i-1])/2.0 for i in range(1,21)])
    Cost_Adhes = sum([(Cadh[i-1]+Cadh[i])*(Lblade[i]-Lblade[i-1])/2.0 for i in range(1,21)])
    Cost_Rtip = [CosMod_params[5][0]*sum([Lblade[i]-Lblade[i-1] for i in range(1,21)]),CosMod_params[5][1]*sum([Lblade[i]-Lblade[i-1] for i in range(1,21)]),CosMod_params[5][2]*sum([Lblade[i]-Lblade[i-1] for i in range(1,21)]),CosMod_params[5][3]*sum([Lblade[i]-Lblade[i-1] for i in range(1,21)]),CosMod_params[5][4]*sum([Lblade[i]-Lblade[i-1] for i in range(1,21)])]
    Cost_Amolds = [2.0*CosMod_params[6][0]*(Askin+Awebs),2.0*CosMod_params[6][1]*(Askin+Awebs),2.0*CosMod_params[6][2]*(Askin+Awebs),2.0*CosMod_params[6][3]*(Askin+Awebs)]
    Cost_Aout = [CosMod_params[7][0]*Askin,CosMod_params[7][1]*Askin,CosMod_params[7][2]*Askin,CosMod_params[7][3]*Askin,CosMod_params[7][4]*Askin]
    Cost_Paint = CosMod_params[2]*1150.0*0.0006*Askin
    Cost_Blade   = Cost_TbaNu+Cost_Light+Cost_Adhes+Cost_Paint+Cost_Resin+Cost_TRIAX+Cost_BIAX+Cost_UNIAX+Cost_BALSA+sum(Cost_Rtip)+sum(Cost_Amolds)+sum(Cost_Aout)+Cost_Labor
    Cost_Hub     = 5.75739*Mass+34280.61
    Cost_Pitch   = 4.30288*(R_tip**2.6578)
    Cost_Spinner = 292.64780*R_tip-4116.84
    if True==False:
        with open('./cost_blade.txt','w') as f:
            f.write('Blades:  %.2f\n'%(3.0*Cost_Blade))
            f.write('Hub:     %.2f\n'%(Cost_Hub))
            f.write('Pitch:   %.2f\n'%(Cost_Pitch))
            f.write('Spinner: %.2f\n'%(Cost_Spinner))
            f.write('--------------\n')
            f.write('ROTOR:   %.2f\n'%(3.0*Cost_Blade+Cost_Hub+Cost_Pitch+Cost_Spinner))
            f.write('##########################\n')
            f.write('###   blade analysis   ###\n')
            f.write('##########################\n')
            f.write('T-bolts & nuts:   %.2f\n'%(Cost_TbaNu))
            f.write('Lightning:        %.2f\n'%(Cost_Light))
            f.write('Adhesive:         %.2f\n'%(Cost_Adhes))
            f.write('Paint:            %.2f\n'%(Cost_Paint))
            f.write('Resin & Hardener: %.2f\n'%(Cost_Resin))
            f.write('materials: %.2f\n'%(Cost_TRIAX+Cost_BIAX+Cost_UNIAX+Cost_BALSA))
            f.write('   TRIAX: %.2f\n'%(Cost_TRIAX))
            f.write('   BIAX:  %.2f\n'%(Cost_BIAX))
            f.write('   UNIAX: %.2f\n'%(Cost_UNIAX))
            f.write('   BALSA: %.2f\n'%(Cost_BALSA))
            f.write('Consumables:    %.2f\n'%(sum(Cost_Rtip)+sum(Cost_Amolds)+sum(Cost_Aout)))
            f.write('   Nonsand tape:       %.2f\n'%(Cost_Rtip[0]))
            f.write('   Chopped strand:     %.2f\n'%(Cost_Rtip[1]))
            f.write('   Tubing:             %.2f\n'%(Cost_Rtip[2]))
            f.write('   Tacky tape:         %.2f\n'%(Cost_Rtip[3]))
            f.write('   Masking tape:       %.2f\n'%(Cost_Rtip[4]))
            f.write('   Peel ply:           %.2f\n'%(Cost_Amolds[0]))
            f.write('   Tackifier adhesive: %.2f\n'%(Cost_Amolds[1]))
            f.write('   Release agent:      %.2f\n'%(Cost_Amolds[2]))
            f.write('   Flow medium:        %.2f\n'%(Cost_Amolds[3]))
            f.write('   Chop fibers:        %.2f\n'%(Cost_Aout[0]))
            f.write('   White lightning:    %.2f\n'%(Cost_Aout[1]))
            f.write('   Hardener:           %.2f\n'%(Cost_Aout[2]))
            f.write('   Putty:              %.2f\n'%(Cost_Aout[3]))
            f.write('   Putty catalyst:     %.2f\n'%(Cost_Aout[4]))
            f.write('Labor & others: %.2f\n'%(Cost_Labor))
            f.write('----------------------------\n')
            f.write('TOTAL: %.2f'%(Cost_Blade))
    return 3.0*Cost_Blade+Cost_Hub+Cost_Pitch+Cost_Spinner
