import math as m
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from Assignment1_Q1_Wilson_Walker import Cp as Cp_wilson
from Assignment1_Q1_Glauert import Cp as Cp_glauert


TSR = np.arange(0,10+1,1)

Cp_glauert_lst = []
Cp_wilson_lst = []


for i in range(len(TSR)):
    Cp_wilson_lst.append(Cp_wilson[i][1])

for i in range(len(TSR)):
    Cp_glauert_lst.append(Cp_glauert[i][3])


plt.figure()
plt.plot(TSR,Cp_glauert_lst, label ='Glauert (pitch = 0)')
plt.plot(TSR,Cp_wilson_lst, label = 'Wilson & Walker (pitch = -2)')
plt.xlabel('TSR')
plt.title('Cp comparison')
plt.legend()
plt.savefig('plots/glauert_wilson.png')
plt.close()
