import math as m

a = 0
aprime = 0
convergenceFactor = 1e-100
delta = 1
deltaPrime = 1

R = 89.17
B = 3
rou = 1.225     #air density
Vo = 8.0        #wind speed
w = 2.61
pitch =  m.radians(-3.0)
twist =  m.radians(2.0)
c =  0.5
Cl =  0.5
Cd =  0.01
relax = 0.1
solidity = (B*c)/(2*m.pi*r)
count = 0

while(delta > convergenceFactor and deltaPrime > convergenceFactor):
    count = count + 1
    if (count > 10000):
        print("No convergence!")
        break

    flowAngle = m.atan(((1-a)*Vo)/((1+aprime)*w*r))
    localalpha =  flowAngle - (pitch + twist)

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

Vrel = m.sqrt(Vo**2+(w*r)**2)
Pn = 0.5*rou*Vrel**2*c*Cn
Pt = 0.5*rou*Vrel**2*c*Ct

print('iterations: ',count)
print('a:', round(a,3))
print('a_prime:', round(aprime,3))
print('Pt:', round(Pt,3))
print('Pn:', round(Pn,3))
print('F:', round(F,3))
print('CT:', round(CT,3))

K - 4*F*m.sin()