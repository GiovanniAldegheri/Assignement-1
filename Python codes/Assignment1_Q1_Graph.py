import math as m
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from Assignment1_Q1_Glauert import Cp as Cp_glauert
from Assignment1_Q1_Wilson_Walker import Cp as Cp_wilson

TSR = np.arange(5,10+1,1)

Cp_glauert_lst = []
Cp_wilson_lst = []

for i in range(len(TSR)):
    Cp_glauert_lst.append(Cp_glauert[i][3])
    Cp_wilson_lst.append(Cp_wilson[i][1])

plt.figure()
plt.plot(TSR,Cp_glauert_lst)
plt.plot(TSR,Cp_wilson_lst)
plt.show()
