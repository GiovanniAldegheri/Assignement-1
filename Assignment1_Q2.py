import math as m
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

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

#Functions____________
def contourplots(pitch, TSR, Cp):
    # Create a figure with two subplots
    fig, axs = plt.subplots(1, 2, figsize=(12, 5))  # 1 row, 2 columns of subplots

    # Subplot 1
    axs[0].set_title(r'$C_p(Chord,Twist)$ Contour Plot')
    [X, Y] = np.meshgrid(TSR, pitch)
    cont1 = axs[0].contourf(Y, X, Cp, 60, cmap="turbo")
    axs[0].set_ylabel(r'Twist (deg)')
    axs[0].set_xlabel(r'Chord (m)')
    cbar1 = plt.colorbar(cont1, ax=axs[0])
    cbar1.set_label(r'$C_p$')

    # Subplot 2
    # axs[1].set_title(r'$C_t(\lambda,\theta_p)$ Contour Plot')
    # [X, Y] = np.meshgrid(pitch, TSR)
    # cont2 = axs[1].contourf(Y, X, Ct)
    # axs[1].set_ylabel(r'$\theta_p$ (deg)')
    # axs[1].set_xlabel(r'$\lambda$ (-)')
    # cbar2 = plt.colorbar(cont2, ax=axs[1])
    # cbar2.set_label(r'$C_t$')

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

def BEM(TSR,pitch,r,c,twist,thick,aoa_tab,cl_tab,cd_tab,cm_tab):
    a = 0
    aprime = 0
    convergenceFactor = (1e-10)
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

    # Vrel = m.sqrt(Vo**2+(w*r)**2)

    # Pn = 0.5*rho*Vrel**2*c*Cn
    # Pt = 0.5*rho*Vrel**2*c*Ct

    Cp = TSR*B*(1-a)**2*Ct*c/(2*m.pi*(m.sin(flowAngle))**2*R)

    return(Cp)

#Constants______________
R = 89.17 #m
B = 3
rho = 1.225 #kg/m3
Vo = 10
r = 71.97

#Interpolate over r, tip speed ratio and pitch
TSR = 8
pitch_fixed = 0
#Blade characteristics
Cp_max = 0
c_max = 0
beta_max = 0
beta = np.arange(-1.11-7,-1.11+3,0.5)
c = np.arange(1,4,0.2)

Cp=np.zeros([len(c),len(beta)])

w = TSR*Vo/R

for i in range(len(c)):
    for j in range(len(beta)):

        Cp[i, j] = BEM(TSR,pitch_fixed,r_ref[9],c[i],beta[j],tc_ref[9],aoa_tab,cl_tab,cd_tab,cm_tab)
        if Cp[i, j] < 0:
            Cp[i, j] = 0


        if (Cp_max < Cp[i,j]):  
            Cp_max = Cp[i,j]
            c_max = c[i]
            beta_max = beta[j]

        print('Cp =',format(Cp[i, j],'.6f'), '\tc =',c[i], '\tbeta =', beta[j])
print('\nBest values', '\nCp =', format(Cp_max,'.6f'), '\tChord =', round(c_max,1), '\tTwist =', round(beta_max,2), '\n')

#Plot the results in a countour plot
contourplots(c, beta, Cp)
plt.show()
