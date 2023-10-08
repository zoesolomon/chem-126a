import matplotlib.pyplot as plt
import numpy as np
import math


def maxwell_boltzman_distribution(m,v,kB,T):
    normalizationFactor = math.sqrt(m/(2*math.pi*kB*T))
    return (normalizationFactor**3)*(4*math.pi*v**2)*(math.exp(-(m*v**2)/(2*kB*T)))

T1 = 200 #K
T2 = 500 #K 
T3 = 2000 #K 
Av = 6.022 * pow(10, 23)
kB = 1.38064852*10**(-23) #m^2 kg s^(-2) K^(-1)



plt.figure(dpi=150)

v0 = 0
v1 = 5000
MH2 = 2.016 * pow(10, -3) #kg/mol
mH2 = MH2/Av #1 molecule in kg
nSamples = 100000
dv = (v1-v0)/nSamples
velocities = np.linspace(v0,v1, nSamples) #m
frequency = []
for v in velocities:
    frequency.append(maxwell_boltzman_distribution(mH2,v,kB,T1))

plt.plot(velocities,frequency)

frequency2 = []
for v in velocities:
    frequency2.append(maxwell_boltzman_distribution(mH2,v,kB,T2))
    
plt.plot(velocities,frequency2)
    
frequency3 = []
for v in velocities:
    frequency3.append(maxwell_boltzman_distribution(mH2,v,kB,T3))


plt.plot(velocities,frequency3)
plt.xlabel("velocity [m/s]")
plt.ylabel("Probability density")
plt.legend(["200K", "500K", "2000K"], loc ="upper right") 
plt.title("Boltzmann Distribution of H2 at Varied Temperatures")

print("dv: {}".format(dv))
print("area under the curve: {}".format(np.sum(frequency)*dv))

plt.savefig('Boltzmann_H2.png')




plt.figure(dpi=150)

v0 = 0
v1 = 5000
MCO2 = 44.01 * pow(10, -3) #kg/mol
mCO2 = MCO2/Av #1 molecule in kg
nSamples = 100000
dv = (v1-v0)/nSamples
velocities = np.linspace(v0,v1, nSamples) #m
frequency = []
for v in velocities:
    frequency.append(maxwell_boltzman_distribution(mCO2,v,kB,T1))

plt.plot(velocities,frequency)

frequency2 = []
for v in velocities:
    frequency2.append(maxwell_boltzman_distribution(mCO2,v,kB,T2))
    
plt.plot(velocities,frequency2)
    
frequency3 = []
for v in velocities:
    frequency3.append(maxwell_boltzman_distribution(mCO2,v,kB,T3))


plt.plot(velocities,frequency3)
plt.xlabel("velocity [m/s]")
plt.ylabel("Probability density")
plt.legend(["200K", "500K", "2000K"], loc ="upper right") 
plt.title("Boltzmann Distribution of CO2 at Varied Temperatures")

print("dv: {}".format(dv))
print("area under the curve: {}".format(np.sum(frequency)*dv))

plt.savefig('Boltzmann_CO2.png')





plt.figure(dpi=150)

v0 = 0
v1 = 5000
MRn = 222 * pow(10, -3) #kg/mol
mRn = MRn/Av #1 molecule in kg
nSamples = 100000
dv = (v1-v0)/nSamples
velocities = np.linspace(v0,v1, nSamples) #m
frequency = []
for v in velocities:
    frequency.append(maxwell_boltzman_distribution(mRn,v,kB,T1))

plt.plot(velocities,frequency)

frequency2 = []
for v in velocities:
    frequency2.append(maxwell_boltzman_distribution(mRn,v,kB,T2))
    
plt.plot(velocities,frequency2)
    
frequency3 = []
for v in velocities:
    frequency3.append(maxwell_boltzman_distribution(mRn,v,kB,T3))


plt.plot(velocities,frequency3)
plt.xlabel("velocity [m/s]")
plt.ylabel("Probability density")
plt.legend(["200K", "500K", "2000K"], loc ="upper right") 
plt.title("Boltzmann Distribution of Radon at Varied Temperatures")

print("dv: {}".format(dv))
print("area under the curve: {}".format(np.sum(frequency)*dv))

plt.savefig('Boltzmann_Radon.png')



'''
plt.figure(dpi=150)

v0 = 0
v1 = 5000
MRn = 222 * pow(10, -3) #kg/mol
mRn = MRn/Av #1 molecule in kg
MCO2 = 44.01 * pow(10, -3) #kg/mol
mCO2 = MCO2/Av #1 molecule in kg
MH2 = 2.016 * pow(10, -3) #kg/mol
mH2 = MH2/Av #1 molecule in kg
nSamples = 100000
dv = (v1-v0)/nSamples
velocities = np.linspace(v0,v1, nSamples) #m

frequency = []
for v in velocities:
    frequency.append(maxwell_boltzman_distribution(mH2,v,kB,T1))

plt.plot(velocities,frequency)

frequency2 = []
for v in velocities:
    frequency2.append(maxwell_boltzman_distribution(mH2,v,kB,T2))
    
plt.plot(velocities,frequency2)
    
frequency3 = []
for v in velocities:
    frequency3.append(maxwell_boltzman_distribution(mH2,v,kB,T3))

plt.plot(velocities,frequency3)
    
frequency4 = []
for v in velocities:
    frequency4.append(maxwell_boltzman_distribution(mCO2,v,kB,T1))

plt.plot(velocities,frequency4)

frequency5 = []
for v in velocities:
    frequency5.append(maxwell_boltzman_distribution(mCO2,v,kB,T2))
    
plt.plot(velocities,frequency5)
    
frequency6 = []
for v in velocities:
    frequency6.append(maxwell_boltzman_distribution(mCO2,v,kB,T3))

plt.plot(velocities,frequency6)

frequency7 = []
for v in velocities:
    frequency7.append(maxwell_boltzman_distribution(mRn,v,kB,T1))

plt.plot(velocities,frequency7)

frequency8 = []
for v in velocities:
    frequency8.append(maxwell_boltzman_distribution(mRn,v,kB,T2))
    
plt.plot(velocities,frequency8)
    
frequency9 = []
for v in velocities:
    frequency9.append(maxwell_boltzman_distribution(mRn,v,kB,T3))
    
plt.plot(velocities,frequency9)
    
    

plt.xlabel("velocity [m/s]")
plt.ylabel("Probability density")
plt.legend(["H2, 200K", "H2, 500K", "H2, 2000K", "CO2, 200K", "CO2, 500K", "CO2, 2000K", "Rn, 200K", "Rn, 500K", 
            "Rn, 2000K"], loc ="upper right") 
plt.title("Boltzmann Distribution of H2, CO2, and Radon at Varied Temperatures")
'''



