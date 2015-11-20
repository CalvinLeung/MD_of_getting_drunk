# Brownian Motion on the Unit Square
# Parameters
# wallParameter - wall spacing
# particleM - effective particle mass
# NP - particle number
# Temp - temperature which dictates thermal distribution (set k_B = 1)

import matplotlib.pyplot as plt
import numpy as np
from numpy import ma


def genIC(particleM,NP,Temp):
    particleX = np.random.rand(NP,2)
    particleX[:,0] = 0.5 * particleX[:,0]
    particleV = np.random.exponential(1/(2 * particleM * Temp),(NP,2))
    return (particleX, particleV)

def plotState(particleX, particleV):
    plt.figure()
    Q = plt.quiver(particleX[:,0],particleX[:,1],particleV[:,0],particleV[:,1])
    plt.show()
    return 0


def updateV(particleV, gamma, particleM, D):

    # gamma = 6*pi*eta*a, where
    # eta = fluid viscosity
    # a = effective radius
    # Characteristic time step tau = a^ / D
    # D = k_B T / gamma
    zeta = np.random.normal(0,1,(particleV.shape))
    dV = 1/particleM * (-gamma * particleV + D*zeta)
    return particleV + dV

def handleCollision(particleX,particleV,wallParameter):
    particleV = particleV

    return particleV

def randomWalk(timeSteps,NP,Temp,particleM,particleR,eta,wallParam):
    kB = 1 # Boltzmann's constant
    gamma = (6*np.pi*eta*particleR)
    D = Temp * kB / gamma
    tau = particleR**2 / D
    
    xHistory = np.zeros((NP,2,timeSteps))
    vHistory = np.zeros((NP,2,timeSteps))
    
    particleX,particleV = genIC(particleM,NP,Temp)
    
    for i in range(0,timeSteps):
        xHistory[:,:,i] = particleX[:,:]
        vHistory[:,:,i] = particleV[:,:]

        particleV = updateV(particleV,gamma,particleM,D)
        particleX = particleX + particleV
        
    return xHistory,vHistory


particleX, particleV = genIC(1,10,1)
#plotState(particleX,particleV)
timeSteps = 37
NP = 3
Temp = 0.00001
particleM = 100
particleR = 1
eta = 1
wallParam = 0  
xHistory,vHistory = randomWalk(timeSteps, NP, Temp, particleM, particleR, eta, wallParam)
for i in range(timeSteps):
    plotState(xHistory[:,:,i],vHistory[:,:,i])
    input('t = '+str(i))
    
if __name__ == "__main__":
    main()
