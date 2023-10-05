import math as m
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

A = 9 #m/s
k = 1.9
P = 10.64e6

V_cutin = 4
V_low = 20
V_cutout = 25
V = np.arange(0,25+1,1) 

h = []
for i in range(len(V)):
    h.append(k/A*(V[i]/A)**(k-1)*m.exp(-(V[i]/A)**k))   #Weibull probability density function

f = m.exp(-(V_low/A)**k) - m.exp(-(V_cutout/A)**k)      #[days] f(20 < V0 < 25) integrated Weibull function
H = 8760 * f                                            #[hours]
AEO = P * H/1e3                                         #[kWh] Same power output at V=20 and V=25 (Q4)

print(H,'Hours lost \n',AEO,'kWH lost')

plt.plot(V,h, color = 'grey')
plt.title('Weibull probability density function')
plt.xlabel(r'$V_0$')
plt.ylabel(r'$h(V_0)$')
plt.fill_between(V, h, where=(V >= 0) & (V <= V_cutin), hatch = '/', color = 'red', alpha=0.2, label = 'Cut-in speed')
plt.fill_between(V, h, where=(V >= V_cutin) & (V <= V_low), color = 'green', alpha=0.2, label = 'Operating range')
plt.fill_between(V, h, where=(V >= V_low) & (V <= V_cutout), hatch = '/',color = 'red', alpha=0.4, label = 'Reduced cut-out speed')
plt.legend(loc='upper right', fontsize=10)

plt.show()