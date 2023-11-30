import math as m
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#CONSTANTS----------------------------------------
R = 89.17 #m
B = 3
rho = 1.225 #kg/m3
TSR_opt = 8
P_max = 10.64e6

#READING FILES--------------------------------------------------------------------

files=['FFA-W3-241.txt','FFA-W3-301.txt','FFA-W3-360.txt','FFA-W3-480.txt','FFA-W3-600.txt','cylinder.txt']
#Initializing tables    
cl_tab=np.zeros([105,6])
cd_tab=np.zeros([105,6])
cm_tab=np.zeros([105,6])
aoa_tab=np.zeros([105,])
#Readin of tables. Only do this once at startup of simulation
for i in range(np.size(files)):
    aoa_tab[:],cl_tab[:,i],cd_tab[:,i],cm_tab[:,i] = np.loadtxt(files[i], skiprows=0).T

thick_prof=np.zeros(6)
thick_prof[0]=24.1;
thick_prof[1]=30.1;
thick_prof[2]=36;
thick_prof[3]=48;
thick_prof[4]=60;
thick_prof[5]=100;

bladedat = pd.read_csv('bladedat.txt',sep="\t", header=None)
r_ref = bladedat[0].tolist() #m
c_ref = bladedat[2].tolist() #m
beta_ref = bladedat[1].tolist() #deg
tc_ref = bladedat[3].tolist() #%

r_ref.pop()

#FUNCTIONS------------------------------------------------------------

#PLOTS

def simple_graph(x,y,ref,x_label,y_label):
    plt.figure()
    plt.grid()
    plt.xlim(4,25)
    #plt.ylim(y_low,y_top)
    plt.plot(x, y, label = 'BEM')
    plt.plot(x, ref, label = 'DTU report')
    #plt.title((y_label + '(' + x_label + ')'))
    plt.legend()
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.savefig('plots/'+y_label+'.png')


def contourplots(x,y,a,x_label,y_label,a_label):
    axs, fig = plt.subplot(1, 1, figsize=(5, 5))
    axs.set_title(r'Contour Plot')
    [X, Y] = np.meshgrid(x, y)
    cont1 = axs.contourf(X, Y, a, 60, cmap="turbo")
    axs.set_xlabel(x_label)
    axs.set_ylabel(y_label)
    cbar1 = plt.colorbar(cont1, ax=axs[0])
    cbar1.set_label(a_label)
    plt.tight_layout()


def force_coeffs(localalpha,thick,aoa_tab,cl_tab,cd_tab,cm_tab):
    cl_aoa=np.zeros([1,6])
    cd_aoa=np.zeros([1,6])
    cm_aoa=np.zeros([1,6])
    
    #Interpolate to current angle of attack:
    for i in range(np.size(files)):
        cl_aoa[0,i]=np.interp (localalpha,aoa_tab,cl_tab[:,i])
        cd_aoa[0,i]=np.interp (localalpha,aoa_tab,cd_tab[:,i])
        cm_aoa[0,i]=np.interp (localalpha,aoa_tab,cm_tab[:,i])
    
    #Interpolate to current thickness:
    Cl=np.interp (thick,thick_prof,cl_aoa[0,:])
    Cd=np.interp (thick,thick_prof,cd_aoa[0,:])
    Cm=np.interp (thick,thick_prof,cm_aoa[0,:])
    return Cl, Cd, Cm 


def BEM(Vo,TSR,pitch,r,c,twist,thick,aoa_tab,cl_tab,cd_tab,cm_tab):
    a = 0
    aprime = 0
    convergenceFactor = 1e-10
    delta = 1
    deltaPrime = 1

    relax = 0.1
    solidity = (B*c)/(2*m.pi*r)
    count = 0

    while(delta > convergenceFactor and deltaPrime > convergenceFactor):
        count = count + 1
        if (count > 1e4):
            print("No convergence!")
            break

        flowAngle = m.atan(((1-a)*R)/((1+aprime)*TSR*r))
        localalpha =  m.degrees(flowAngle) - (pitch + twist)

        Cl,Cd,Cm = force_coeffs(localalpha,thick,aoa_tab,cl_tab,cd_tab,cm_tab)
        Ct = Cl*m.sin(flowAngle) - Cd*m.cos(flowAngle)
        Cn = Cl*m.cos(flowAngle) + Cd*m.sin(flowAngle)

        F = 2/m.pi*m.acos(m.exp(-B*(R-r)/(2*r*m.sin(abs(flowAngle)))))

        CT = ((1-a)**2*Cn*solidity)/m.sin(flowAngle)**2

        aold = a

        if(aold < 0.33):
            a = (solidity*Cn*(1-aold))/(4*F*m.sin(flowAngle)**2)
        else:
            aStar = CT/(4*F*(1-1/4*(5-3*aold)*aold))
            a = relax*aStar + (1-relax)*aold

        aprimeOld  = aprime
        aprimeStar = (solidity*Ct*(1+aprimeOld))/(4*F*m.sin(flowAngle)*m.cos(flowAngle))
        aprime = relax*aprimeStar + (1-relax)*aprimeOld

        delta = abs(aprime - aprimeOld)
        deltaPrime = abs(aprime - aprimeOld)

    w = TSR*Vo/R
    Vrel = m.sqrt(Vo**2+(w*r)**2)

    Pn = 0.5*rho*Vrel**2*c*Cn
    Pt = 0.5*rho*Vrel**2*c*Ct

    return(Pn, Pt)


#Variables__________
Vo = np.arange(4,25+1,1)
pitch = 0
pitch_delta = 0.1

#Results_____________
P_lst = np.zeros([len(Vo)])
T_lst = np.zeros([len(Vo)])
Cp_lst = np.zeros([len(Vo)])
Ct_lst = np.zeros([len(Vo)])
pitch_lst = np.zeros([len(Vo)])

#MAIN-----------------------------------------------------------------

for i in range(len(Vo)):
    w = TSR_opt*Vo[i]/R
    if w > 1.023:
        w = 1.023
    TSR = w*R/Vo[i]
    P = P_max+1
    T = 0
    pitch += pitch_delta
    print('Vo =',Vo[i], 'w =', w)

    while P > P_max and pitch < 90:
        pitch -= pitch_delta
        Pn_lst = []
        Pt_lst = []

        for j in range(len(r_ref)):
            Pn, Pt = BEM(Vo[i],TSR,pitch,r_ref[j],c_ref[j],beta_ref[j],tc_ref[j],aoa_tab,cl_tab,cd_tab,cm_tab)
            Pn_lst.append(Pn)
            Pt_lst.append(Pt*r_ref[j])

        T = np.trapz(Pn_lst,r_ref)*B
        P = np.trapz(Pt_lst,r_ref)*w*B
        print('\t',round(pitch,1), round(P/1e6,2))

    pitch_lst[i] = pitch
    P_lst[i] = P
    T_lst[i] = T
    Cp_lst[i] = P/(0.5*rho*Vo[i]**3*m.pi*R**2)
    Ct_lst[i] = T/(0.5*rho*Vo[i]**2*m.pi*R**2)

# for i in range(len(Vo)):
#     print(Vo[i], pitch_lst[i], P_lst[i]/1e6)

for i in range(len(Vo)):
    print('Vo(m/s) =', Vo[i], 'T(kN) =', int(T_lst[i]/1000), 'Pitch =', round(pitch_lst[i],1))

ref_pitch = [2.751, 1.966, 0.896, 0, 0, 0, 0, 0, 4.502, 7.266, 9.292, 10.958, 12.499, 13.896, 15.200, 16.432, 17.618, 18.758, 19.860, 20.927, 21.963, 22.975]
ref_power = [280.2, 799.1, 1532.7, 2506.1, 3730.7, 5311.8, 7286.5, 9698.3, 10639.1, 10648.5, 10639.3, 10683.7, 10642.0, 10640.0, 10639.9, 10652.8, 10646.2, 10644.0, 10641.2, 10639.5, 10643.6, 10635.7]
ref_thrust = [225.9, 351.5, 498.1, 643.4, 797.3, 1009.1, 1245.8, 1507.4, 1325.1, 1082.0, 967.9, 890.8, 824.8, 774.0, 732.5, 698.4, 668.1, 642.1, 619.5, 599.8, 582.7, 567.2]
ref_cp = [0.286, 0.418, 0.464, 0.478, 0.476, 0.476, 0.476, 0.476, 0.423, 0.317, 0.253, 0.207, 0.170, 0.142, 0.119, 0.102, 0.087, 0.075, 0.065, 0.057, 0.050, 0.044]
ref_ct = [0.923, 0.919, 0.904, 0.858, 0.814, 0.814, 0.814, 0.814, 0.602, 0.419, 0.323, 0.259, 0.211, 0.175, 0.148, 0.126, 0.109, 0.095, 0.084, 0.074, 0.066, 0.059]


simple_graph(Vo, pitch_lst, ref_pitch, 'Vo (m/s)', 'Pitch (deg)')
simple_graph(Vo, P_lst/1000, ref_power, 'Vo (m/s)', 'Power (kW)')
simple_graph(Vo, T_lst/1000, ref_thrust, 'Vo (m/s)', 'Thrust (kN)')
simple_graph(Vo, Cp_lst, ref_cp, 'Vo (m/s)', 'Cp')
simple_graph(Vo, Ct_lst, ref_ct, 'Vo (m/s)', 'Ct')
plt.close()