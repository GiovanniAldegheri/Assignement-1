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

r_ref.pop()

#Functions____________
def contourplots(pitch, TSR, Cp, Ct):
    # Create a figure with two subplots
    fig, axs = plt.subplots(1, 2, figsize=(12, 5))  # 1 row, 2 columns of subplots

    # Subplot 1
    axs[0].set_title(r'$C_p(\lambda,\theta_p)$ Contour Plot')
    [X, Y] = np.meshgrid(pitch, TSR)
    cont1 = axs[0].contourf(Y, X, Cp, 60, cmap="turbo")
    axs[0].set_ylabel('Pitch angle (deg)')
    axs[0].set_xlabel('Tip speed ratio (-)')
    cbar1 = plt.colorbar(cont1, ax=axs[0])
    cbar1.set_label(r'$C_p$')

    # Subplot 2
    axs[1].set_title(r'$C_t(\lambda,\theta_p)$ Contour Plot')
    [X, Y] = np.meshgrid(pitch, TSR)
    cont2 = axs[1].contourf(Y, X, Ct, 60, cmap="turbo")
    axs[1].set_ylabel('Pitch angle (deg)')
    axs[1].set_xlabel('Tip speed ratio (-)')
    cbar2 = plt.colorbar(cont2, ax=axs[1])
    cbar2.set_label(r'$C_t$')

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
        aprimeOld  = aprime

        # if(aold > 0.33):
        #         aStar = CT/(4*F*(1-1/4*(5-3*aold)*aold))
        #         a = relax*aStar + (1-relax)*aold
        if(aold > 0.2):
            K = 4*F*(m.sin(flowAngle))**2/(solidity*Cn)
            a = 1+K/2*(1-2*0.2)-0.5*np.sqrt((K*(1-2*0.2)+2)**2+4*(K*0.2**2-1))  
        else:
            a = solidity*Cn/4/F/(m.sin(flowAngle))**2*(1-aold)
            aprime = solidity*Ct/4/F/m.sin(flowAngle)/m.cos(flowAngle)*(1+aprimeOld)
        
            
        aprimeStar = (solidity*Ct*(1+aprimeOld))/(4*F*m.sin(flowAngle)*m.cos(flowAngle))
        aprime = relax*aprimeStar + (1-relax)*aprimeOld

        delta = abs(aprime - aprimeOld)
        deltaPrime = abs(aprime - aprimeOld)

    #Cp = (B*TSR*Ct*(1-a)**2/(2*m.pi*(m.sin(flowAngle))**2)*c/R)
    Vrel = m.sqrt(Vo**2+(w*r)**2)

    Pn = 0.5*rho*Vrel**2*c*Cn
    Pt = 0.5*rho*Vrel**2*c*Ct
    if (m.isnan(Pt)|(m.isnan(Pn))):
        Pt, Pn = 0,0
    return(Pn, Pt)

#Constants______________
R = 89.17 #m
B = 3
rho = 1.225 #kg/m3

#Interpolate over r, tip speed ratio and pitch
Vo = [5,9,11,20]
TSR = [8,8,8,4.56]
pitch = [0,0,0,17.2]

#Blade characteristics
P_max = 0
Cp_max = 0
TSR_max = 0
pitch_max = 0

Pn_lst=np.zeros([len(TSR),len(r_ref)])
Pt_lst=np.zeros([len(TSR),len(r_ref)])

for i in range(len(TSR)):
    w = TSR[i]*Vo[i]/R

    for k in range(len(r_ref)):
        Pn, Pt = BEM(Vo[i],TSR[i],pitch[i],r_ref[k],c_ref[k],beta_ref[k],tc_ref[k],aoa_tab,cl_tab,cd_tab,cm_tab)
        Pn_lst[i][k] = Pn
        Pt_lst[i][k] = Pt

r_str = '0 2.64302 5.37977 8.20274 11.1031 14.071 17.0953 20.1641 23.2647 26.3837 29.5076 32.6228 35.7156 38.773 41.7824 44.732 47.6111 50.4099 53.1201 55.7344 58.247 60.6534 62.9501 65.1352 67.2076 69.1675 71.0159 72.7545 74.386 75.9133 77.3402 78.6705 79.9085 81.0585 82.1252 83.113 84.0265 84.8703 85.6487 86.366'
r_lst = r_str.split(' ')


for i in range(len(r_lst)):
    r_lst[i] = float(r_lst[i])
    r_lst[i] += 2.8

thrust_str = ['','','','']

thrust_str[0] = '39.181 47.6918 113.102 126.034 136.773 158.224 180.691 394.196 825.71 934.533 920.623 1010.61 1099.24 1193.01 1288.35 1385.36 1485.54 1585.33 1680.78 1774.61 1862.36 1941.78 2011.19 2069.71 2121.95 2163.43 2201.56 2226 2247.97 2257.26 2258.47 2246.17 2218.07 2167.01 2091.1 1993.59 1849.02 1640.07 1324.34 0'
thrust_str[1] = '127.914 156.736 369.42 404.966 439.471 504.316 573.861 1319.74 2679.27 3032.55 2985.57 3276.6 3564.36 3871.21 4181.8 4498.18 4825.05 5150.84 5462.62 5769.13 6055.88 6315.54 6544.31 6736.17 6907.38 7043.63 7168.94 7249.76 7322.43 7353.83 7358.85 7319.77 7229.16 7063.65 6817.09 6499.94 6029.32 5348.7 4319.95 0'
thrust_str[2] = '191.181 234.234 552.101 605.305 656.889 753.905 857.941 1971.07 4003.55 4531.41 4461.24 4896.09 5326.24 5784.67 6248.71 6721.38 7209.7 7696.4 8162.15 8620.03 9048.38 9436.26 9778 10064.5 10320.2 10523.7 10710.9 10831.6 10940 10986.9 10994.3 10935.8 10800.4 10553 10184.6 9710.72 9007.56 7990.68 6453.71 0'
thrust_str[3] = '618.785 722.228 1105.03 1093.73 1090.78 1087.18 1067.1 2816.05 5358 5502.45 4357.69 4145.8 3853.26 3643.58 3360.81 3102.81 2883 2693.75 2525.88 2389.39 2268.5 2153.33 2041.7 1932.71 1832.02 1729.84 1647.07 1570.54 1507.71 1456.8 1409.22 1368.38 1320.99 1270.21 1213.33 1134.24 1034.75 898.359 718.062 0'

thrust_lst = ['','','','']

for i in range(4):
    thrust_str[i] = thrust_str[i].replace(',','')
    thrust_str[i] = thrust_str[i].replace(';','')
    thrust_lst[i] = thrust_str[i].split(' ')
    for k in range(len(thrust_lst[i])):
        thrust_lst[i][k] = float(thrust_lst[i][k])

torque_str = ['','','','']

torque_str[0] = '-11.5077 -25.3595 44.5823 34.8795 24.5286 19.7682 8.92608 66.7586 189.515 188.952 192.985 192.582 191.923 194.591 194.247 193.902 193.537 193.112 192.629 192.081 191.491 190.873 190.166 189.311 188.234 186.86 185.05 182.715 179.705 175.902 171.153 165.302 158.158 149.505 139.066 126.407 110.882 91.2024 63.8843 0'
torque_str[1] = '-43.4309 -90.5253 133.316 101.501 72.2097 54.0678 17.9097 228.222 598.126 597.024 611.764 610.977 609.346 618.904 618.209 617.447 616.543 615.398 614.044 612.446 610.716 608.918 607.015 604.547 601.366 597.278 591.774 584.628 575.274 563.376 548.411 529.883 507.171 479.59 446.24 405.689 355.906 292.734 204.95 0'
torque_str[2] = '-64.886 -135.254 199.318 151.841 108.032 80.9981 27.0183 340.671 894.033 892.388 914.4 913.225 910.851 925.124 924.084 922.942 921.587 919.873 917.846 915.456 912.87 910.18 907.334 903.638 898.88 892.765 884.532 873.843 859.853 842.06 819.683 791.98 758.024 716.791 666.938 606.323 531.913 437.496 306.3 0'
torque_str[3] = '-140.68 -261.428 762.145 617.094 530.782 417.67 278.485 1505.42 3110.04 2924.88 2214.28 1946.43 1676.39 1492.11 1283.59 1107.55 963.603 844.48 744.549 663.995 595.862 535.754 482.051 433.701 391.202 351.62 319.469 291.503 268.348 249.922 233.428 219.865 205.998 192.95 180.065 163.265 144.504 120.825 92.3798 0'

torque_lst = ['','','','']

for i in range(4):
    torque_str[i] = torque_str[i].replace(',','')
    torque_str[i] = torque_str[i].replace(';', '')
    torque_lst[i] = torque_str[i].split(' ')
    for k in range(len(torque_lst[i])):
        torque_lst[i][k] = float(torque_lst[i][k])

wind_speed = ['5m/s','9m/s','11m/s','20m/s']

label = [(0, 0), (0, 1), (1, 0), (1, 1)]

plt.figure(figsize=(12,10))
#plt.suptitle('Thrust comparison BEM / Ashes', size = 20)
plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.4, 
                    hspace=0.4)
for i in range(4):
    plt.subplot(2,2,i+1)
    plt.grid()
    plt.xlabel('Radius [m]')
    plt.ylabel('Thrust force [N]')
    plt.title(wind_speed[i])
    plt.plot(r_ref,Pn_lst[i], label='BEM')
    plt.plot(r_lst,thrust_lst[i], label='Ashes')
    plt.legend()
plt.savefig('plots/thrust_Ashes.png')

plt.figure(figsize=(12,10))
#plt.suptitle('Torque comparison BEM / Ashes', size = 20)
plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.4, 
                    hspace=0.4)
for i in range(4):
    plt.subplot(2,2,i+1)
    plt.grid()
    plt.xlabel('Radius [m]')
    plt.ylabel('Torque force [N]')
    plt.title(wind_speed[i])
    plt.plot(r_ref,Pt_lst[i], label='BEM')
    plt.plot(r_lst,torque_lst[i], label='Ashes')
    plt.legend()
plt.savefig('plots/torque_Ashes.png')

plt.close()

#Plot the results in a countour plot
# contourplots(pitch, TSR, Cp, Ct)
# plt.show()
