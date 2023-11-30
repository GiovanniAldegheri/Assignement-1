import math as m
import numpy as np
import matplotlib.pyplot as plt

Cp_max = 0.469
rou = 1.225
R = 89.17
area = m.pi*R**2
P = 10.64e6
TSR = 8
cut_in = 4
cut_out = 25

Vo_rat = (P/(0.5*Cp_max*rou*area))**(1/3)
w_max = TSR*Vo_rat/R

print('Rated wind speed:', round(Vo_rat,3), '\nRotational speed:',round(w_max,3))

Vo = np.arange(0,40,0.1)
w = TSR*Vo/R

for i in range(len(Vo)):
    if Vo[i] > cut_out:
        w[i] = 0
    elif Vo[i] > Vo_rat:
        w[i] = w_max
    elif Vo[i] < cut_in:
        w[i] = 0

plt.plot(Vo, w)
plt.xlabel('Vo (m/s)')
plt.ylabel('w (rad/s)')
plt.xlim(0, 30)
plt.grid()
plt.show()