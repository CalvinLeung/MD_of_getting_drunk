# Brownian Motion on the Unit Square
# Parameters
# wallParameter - wall spacing
# particleM - effective particle mass
# NP - particle number
# Temp - temperature which dictates thermal distribution (set k_B = 1)

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from numpy import ma


def genIC(boxSize, kB, NP, particleM, Temp):
    particleX = boxSize*np.random.rand(NP,3)
    particleX[:,0] = 0.5 * particleX[:,0]
    #speed = np.random.normal(0,particleM/(2 * kB * Temp),(NP,1)) # Boltzmann factor for some speed
    #angle = np.random.rand(NP,1) * 2 * np.pi
    #particleV = np.hstack((np.multiply(speed,np.cos(angle)),np.multiply(speed,np.sin(angle))))

    particleV = np.random.normal(0,np.sqrt(kB*Temp/particleM),(NP,3))
    #particleV = np.zeros((NP,3))
    return (particleX, particleV)

def getState(boxSize,particleV, particleX,plot=False):
    #3D quiver plots and particle counting
    NP = particleX.shape[0]
    if plot:
        normalizedX = np.copy(particleX) 
        normalizedV = np.copy(particleV)
    
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        plt.axis([0,boxSize,0,boxSize])
        ax.set_zbound(lower=0, upper=boxSize)
        ax.quiver(normalizedX[:,0],normalizedX[:,1],normalizedX[:,2],normalizedV[:,0],normalizedV[:,1],normalizedV[:,2],
                  length = 0.05 * boxSize)
        plt.show()
    left = np.sum(particleX[:,0] < boxSize * 0.5,dtype=int)
    return left, NP-left #left and right halves 

def totalEnergy(kB,particleM,particleV,Temp):
    NP = particleV.shape[0]
    T = 0.5 * particleM * sum(particleV[:,0] ** 2 + particleV[:,1] ** 2 + particleV[:,2] ** 2) / NP
    print('Observed average kinetic energy per particle (NP = ' + str(NP) + ') = ' +str(T))
    print('Equipartition predicts (6-3 = 3) QDOF: ' + str(1.5 * kB *Temp))
    return T
    
def updateV(gamma,kB,particleM,particleV,Temp,timeStep):
    stoch = np.sqrt(24*kB*Temp*gamma/timeStep)*np.random.normal(0,np.sqrt(1.0/12.0),(particleV.shape[0],3))
    drag = -gamma * particleV
    dV = ((drag + stoch)/particleM) * timeStep
    return particleV + dV

# def handleCollision(boxSize, particleV,projectedX,wallParam):
#     # Reflect off a wall
#     b = boxSize
#     for i in range(projectedX.shape[0]):
#         x,y,z = projectedX[i,:]
#         xd,yd,zd = particleV[i,:]
#         while(x < 0 or x > b or y < 0 or y > b or z<0 or z>b):
#             if x < 0:
#                 x = abs(x)
#                 xd = -xd
#             elif x > b:
#                 x = 2*b - x
#                 xd = -xd
#             elif y < 0:
#                 y = abs(y)
#                 yd = -yd
#             elif y > b:
#                 y = 2*b - y
#                 yd = -yd           
#             elif z < 0:
#                 z = abs(z)
#                 zd = -zd
#             elif z > b:
#                 z = 2*b - z
#                 zd = -zd
#         projectedX[i,:] = [x,y,z]
#         particleV[i,:] = [xd,yd,zd]
#     return projectedX,particleV

def handleCollision(boxSize, particleV, projectedX, wallParam, NP):
    toohigh = projectedX > boxSize
    toolow  = projectedX < 0

    particleV[toolow] = -particleV[toolow]
    particleV[toohigh]= -particleV[toohigh]

    projectedX[toolow] = -projectedX[toolow]
    projectedX[toohigh] = 2*boxSize - projectedX[toohigh]
    return projectedX, particleV

def randomWalk(boxSize,eta,kB,NP,particleM,particleR,Temp,timeStep,totalSteps,wallParam):
    # gamma = 6*pi*eta*a, where
    # eta = fluid viscosity
    # a = effective radius
    lefts = np.zeros(totalSteps)
    gamma = (6*np.pi*eta*particleR)   
    #xHistory = np.zeros((NP,2,totalSteps))
    #vHistory = np.zeros((NP,2,totalSteps))
    particleX,particleV = genIC(boxSize,kB,NP,particleM,Temp)
    left,right = getState(boxSize, particleV, particleX, True)
    input('Press Enter')
    totalEnergy(kB,particleM,particleV,Temp)
    
    for i in range(0,totalSteps):
        lefts[i],right = getState(boxSize, particleV, particleX, plot=False)
        if i%10000 == 0:
            print("Timestep: "   +  str(i))
            print("# on left: "  + str(lefts[i]))
            print("# on right: " + str(right))
        particleV = updateV(gamma,kB,particleM,particleV,Temp,timeStep)
        projectedX = particleX + particleV * timeStep
        #Handle collisions
        particleX, particleV = handleCollision(boxSize, particleV,projectedX, wallParam, NP)
    totalEnergy(kB,particleM,particleV,Temp)
    return lefts,particleX,particleV

kB = 1.4*10**(-23) # Boltzmann's constant
boxSize = 2e-6 # Our box is 2 mu m right now 
#totalTime = 1e-8 # seconds of diffusion
NP = 100
Temp = 310 # 310 Kelvin
particleM = (162.0/18.0)*(3 * 10 ** (-26)) # molecular mass of nicotine in kg
particleR = 3 * 10 ** (-10) # nicotine is a three angstrom radius sphere, yolo
eta = 1.0 * 10 ** (-6) # viscosity of water in m^2 / s
wallParam = 0
gamma = 6*np.pi*eta*particleR
#D = Temp * kB / gamma
#tau = particleR**2 / D
#print("tau = " + str(tau))
#print("tv = " + str(particleM/gamma))
timeStep = 1e-13
totalSteps = 10000000 #int(totalTime / timeStep)

lefts,particleX,particleV = randomWalk(boxSize,eta,kB,NP,particleM,particleR,Temp,timeStep,totalSteps,wallParam)
getState(boxSize,particleV,particleX,True)
    
#if __name__ == "__main__":
#    main()
