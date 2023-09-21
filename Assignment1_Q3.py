import math as m

Cp_max = 0.469
rou = 1.225
R = 89.17
area = R*m.pi**2
P = 10.64e6
TSR = 8

Vo = (P/(0.5*Cp_max*rou*area))**(1/3)
w = TSR*Vo/R

print(round(w,3))