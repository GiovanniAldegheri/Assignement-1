import math as m
from matplotlib import pyplot as plt
import numpy as np

Cp_max = 0.469
rou = 1.225
R = 89.17
area = R*m.pi**2
P = 10.64e6
TSR = 8

Vo = (P/(0.5*Cp_max*rou*area))**(1/3)
wmax = TSR*Vo/R

print(round(wmax,3))

V = np.arange(0,26,1)
w = []
for i in range(len(V)):
    w.append(TSR*V[i]/R)

plt.plot(V,w)
plt.xlabel(r'$V_o$ (m/s)')
plt.ylabel(r'$\omega$ (rad/s)')
plt.show()