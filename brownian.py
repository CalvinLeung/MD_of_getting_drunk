# Brownian Motion on the Unit Square
# Parameters
# wallParameter - wall spacing
# particleM - effective particle mass
# NP - particle number
# Temp - temperature which dictates thermal distribution (set k_B = 1)

import matplotlib.pyplot as plt
import numpy as np
from numpy import ma


def genIC(NP,kB, Temp,particleM):
    particleX = np.random.rand(NP,2)
    particleX[:,0] = 0.5 * particleX[:,0]
    #speed = np.random.normal(0,particleM/(2 * kB * Temp),(NP,1)) # Boltzmann factor for some speed
    #angle = np.random.rand(NP,1) * 2 * np.pi
    #particleV = np.hstack((np.multiply(speed,np.cos(angle)),np.multiply(speed,np.sin(angle))))

    particleV = np.random.normal(0,2*kB*Temp/particleM,(NP,2))
    return (particleX, particleV)

def plotState(particleX, particleV):
    plt.figure()
    Q = plt.quiver(particleX[:,0],particleX[:,1],particleV[:,0],particleV[:,1])
    plt.show()
    return 0

def totalEnergy(particleV,m,kB,Temp):
    NP = particleV.shape[0]
    T = 0.5 * particleM * sum(particleV[:,0] ** 2 + particleV[:,1] ** 2) / NP
    print('Observed average kinetic energy per particle (NP = ' + str(NP) + ') = ' +str(T))
    print('Equipartition predicts 2 QDOF: ' + str(kB*Temp))
    
def updateV(particleV, particleM, timeStep, gamma, kB, Temp):
    #speed = np.sqrt(np.random.normal(0,2*gamma*kB*Temp,(particleV.shape[0],1))) # Fluctuation Disspation Theorem
    #theta = 2*np.pi*np.random.rand(particleV.shape[0],1) #Random scattering angle
    #stoch = np.hstack((np.multiply(zeta,np.cos(theta)),np.multiply(zeta,np.sin(theta))))
    #vx = np.random.normal(0,2*gamma*kB*Temp,(particleV.shape[0],1))
    #vy = np.random.normal(0,2*gamma*kB*Temp,(particleV.shape[0],1))
    #stoch = np.hstack((vx,vy))
    stoch = 0
    dV = (-gamma * particleV + stoch)/particleM * timeStep
    return particleV + dV

def handleCollision(projectedX,particleV,wallParam):
    particleX = projectedX
    # Reflect off a wall
    for i in range(projectedX.shape[0]):
        x,y = projectedX[i,:]
        xd,yd = particleV[i,:]
        while(x < 0 or x > 1 or y < 0 or y > 1):
            print('wall collision')
            print('(x,y) = ' + str((x,y)) + ' being reflected...')
            if x < 0:
                x = abs(x)
                xd = -xd
            elif x > 1:
                x = 2 - x
                xd = -xd
            elif y < 0:
                y = abs(y)
                yd = -yd
            elif y > 1:
                y = 2 - y
                yd = -yd           
        particleX[i,:] = [x,y]
        particleV[i,:] = [xd,yd]
    return particleX,particleV

def randomWalk(totalSteps,timeStep,NP,Temp,particleM,particleR,eta,kB,wallParam):
    # gamma = 6*pi*eta*a, where
    # eta = fluid viscosity
    # a = effective radius

    gamma = (6*np.pi*eta*particleR)   
    xHistory = np.zeros((NP,2,totalSteps))
    vHistory = np.zeros((NP,2,totalSteps))
    
    particleX,particleV = genIC(NP,kB,Temp,particleM)
    #print(particleX)
    #print(particleV)
    for i in range(0,totalSteps):
        #print(i)
        xHistory[:,:,i] = particleX[:,:]
        vHistory[:,:,i] = particleV[:,:]
        
        particleV = updateV(particleV,particleM,timeStep,gamma,kB,Temp)
        print(totalEnergy(particleV,particleM,kB,Temp))
        input('Press enter')
        projectedX = particleX + particleV * timeStep
        #Handle collisions
        #particleX,particleV = handleCollision(projectedX,particleV, wallParam)
        particleX = projectedX
        #input('Press Enter')
    return xHistory,vHistory


#plotState(particleX,particleV)
kB = 1.4e-23 # Boltzmann's constant

totalTime = 1e-8 # seconds of diffusion
NP = 1
Temp = 300 # 300 Kelvin
particleM = 3 * 10 ** -26 # molecular mass of water in kg
particleR = 1 * 10 ** -10 # water is a one angstrom radius sphere, yolo
eta = 1.0 * 10 ** (-6) # viscosity of water in m^2 / s
wallParam = 0
gamma = 6*np.pi*eta*particleR
D = Temp * kB / gamma
tau = particleR**2 / D
timeStep = tau
totalSteps = 1000 #int(totalTime / timeStep)

particleX, particleV = genIC(particleM,10,kB,Temp)

xHistory,vHistory = randomWalk(totalSteps, timeStep, NP, Temp, particleM, particleR, eta, kB,wallParam)
#for i in range(timeSteps):
#    plotState(xHistory[:,:,i],vHistory[:,:,i])
#    input('t = '+str(i))
    
#if __name__ == "__main__":
#    main()
