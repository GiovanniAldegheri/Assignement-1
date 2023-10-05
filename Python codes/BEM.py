import math

a = 0
aprime = 0
convergenceFactor = 1e-100
delta = 1
deltaPrime = 1

R = 31
r = 24.5
B = 3
rou = 1.225
Vo = 8.0
w = 2.61
pitch =  math.radians(-3.0)
twist =  math.radians(2.0)
c =  0.5
Cl =  0.5
Cd =  0.01
relax = 0.1
solidity = (B*c)/(2*math.pi*r)
count = 0

while(delta > convergenceFactor and deltaPrime > convergenceFactor):
    count = count + 1
    if (count > 10000):
        print("No convergence!")
        break

    flowAngle = math.atan(((1-a)*Vo)/((1+aprime)*w*r))
    localalpha =  flowAngle - (pitch + twist)

    Ct = Cl*math.sin(flowAngle) - Cd*math.cos(flowAngle)
    Cn = Cl*math.cos(flowAngle) + Cd*math.sin(flowAngle)


    F = 2/math.pi*math.acos(math.exp(-B*(R-r)/(2*r*math.sin(abs(flowAngle)))))
 

    CT = ((1-a)**2*Cn*solidity)/math.sin(flowAngle)**2

    aold = a

    if(aold < 0.33):
        a = (solidity*Cn*(1-aold))/(4*F*math.sin(flowAngle)**2)
    else:
        aStar = CT/(4*F*(1-1/4*(5-3*aold)*aold))
        a = relax*aStar + (1-relax)*aold

    aprimeOld  = aprime
    aprimeStar = (solidity*Ct*(1+aprimeOld))/(4*F*math.sin(flowAngle)*math.cos(flowAngle))
    aprime = relax*aprimeStar + (1-relax)*aprimeOld

    delta = abs(aprime - aprimeOld)

    deltaPrime = abs(aprime - aprimeOld)

Vrel = math.sqrt(Vo**2+(w*r)**2)
Pn = 0.5*rou*Vrel**2*c*Cn
Pt = 0.5*rou*Vrel**2*c*Ct

print('iterations: ',count)
print('a:', round(a,3))
print('a_prime:', round(aprime,3))
print('Pt:', round(Pt,3))
print('Pn:', round(Pn,3))
print('F:', round(F,3))
print('CT:', round(CT,3))


Cp_max = 0
TSR_max = 0
pitch_max = 0