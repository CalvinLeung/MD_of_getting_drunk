# Brownian Motion on the real line
# Parameters
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from numpy import ma


def genIC(boxSize, kB, NP, particleM, Temp):
    particleX = boxSize*(1 - 0.5 * np.random.rand(NP,1))
    particleV = np.random.normal(0,np.sqrt(kB*Temp/particleM),(NP,1))
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
        ax.quiver(normalizedX[:,0],normalizedX[:,1],normalizedX[:,2],normalizedV[:,0],normalizedV[:,1],normalizedV[:,2], length = 0.05 * boxSize)
        ax.set_xlabel('X axis (normal to membrane)')
        ax.set_ylabel('Y axis (impenetrable walls)')
        ax.set_zlabel('Z axis (periodic boundary )')
        plt.show()
    left = np.sum(particleX[:,0] < boxSize * 0.5,dtype=int)
    if plot:
        return left, NP-left,ax #left and right halves
    return left, NP-left

def totalEnergy(kB,particleM,particleV,Temp):
    NP = particleV.shape[0]
    T = 0.5 * particleM * sum(particleV[:] ** 2) / NP
    print('Observed average kinetic energy per particle (NP = ' + str(NP) + ') = ' +str(T))
    print('Equipartition predicts 1 QDOF: ' + str(0.5 * kB *Temp))
    return T
    
def updateV(gamma,kB,particleM,particleV,Temp,timeStep):
    stoch = np.sqrt(24*kB*Temp*gamma/timeStep)*np.random.normal(0,np.sqrt(1.0/12.0),(particleV.shape[0],1))
    drag = -gamma * particleV
    dV = ((drag + stoch)/particleM) * timeStep
    return particleV + dV

def membraneCheck(holeDensityVol,particleR,pTrans):
    #FIX THIS
    return pTrans

def handleCollision(boxSize, holeDensityVol,particleV, projectedX, wallParam, NP, particleR,pTrans):
    # Check if x coordinate are out of the range (0, boxSize)
    toohigh = projectedX > boxSize - particleR
    toolow  = projectedX < 0 + particleR

    rejected = np.random.rand(projectedX.shape[0],projectedX.shape[1]) > membraneCheck(holeDensityVol, particleR,pTrans)

    hitMembraneRight = (particleV[:] >= 0)* (np.abs(projectedX - boxSize * 0.5) < particleR) * rejected    #incident from left, rejected by membrane
    hitMembraneLeft = (particleV[:] < 0)  * (np.abs(projectedX - boxSize * 0.5) < particleR) * rejected    #incident from right, rejected by membrane
    
    
    projectedX[hitMembraneRight] = boxSize - 2* particleR - projectedX[hitMembraneRight]
    projectedX[hitMembraneLeft]  = boxSize + 2* particleR - projectedX[hitMembraneLeft]

    bounce = toolow + hitMembraneRight + toohigh + hitMembraneLeft
    particleV[bounce] = -particleV[bounce]
    #projectedX[hitMembraneRight,0] = 2*(0.5*boxSize - particleR) - projectedX[hitMembraneRight,0]
    #projectedX[hitMembraneLeft,0] = 2*(0.5*boxSize + particleR) - projectedX[hitMembraneLeft,0] 

    return projectedX, particleV

def randomWalk(boxSize,eta,holeDensityVol,kB,NP,particleM,particleR,particleV0,particleX0,pTrans,Temp,timeStep,totalSteps,wallParam):
    # gamma = 6*pi*eta*a, where
    # eta = fluid viscosity
    # a = effective radius
    lefts = np.zeros(totalSteps/10000)
    gamma = (6*np.pi*eta*particleR)
    (particleX,particleV) = (particleX0,particleV0)
    
    for i in range(0,totalSteps):
        if i%10000 == 0:
            left,right = getState(boxSize, particleV, particleX, plot=False)
            lefts[i/10000] = left
            print("Timestep: "   +  str(i))
            print("# on left: "  + str(left))
            print("# on right: " + str(right))
            np.append(lefts,left)
        #Step forward
        projectedX = particleX + particleV * timeStep
        #Handle collisions
        particleX, particleV = handleCollision(boxSize, holeDensityVol,particleV,projectedX, wallParam, NP, particleR,pTrans)
        #Update the velocities stochastically
        particleV = updateV(gamma,kB,particleM,particleV,Temp,timeStep)
    totalEnergy(kB,particleM,particleV,Temp)
    return lefts,particleX,particleV

def TSPlot(timeSeries,timeStep):
    N = np.arange(0,timeSeries.shape[0],1)
    fig = plt.figure()
    plt.plot(N * timeStep, timeSeries)
    plt.show()
    return plt

pTrans = 0.1
kB = 1.4*10**(-23) # Boltzmann's constant
boxSize = 2e-6 # Our box is 2 mu m right now 
#totalTime = 1e-8 # seconds of diffusion
NP = 1000 # typical ethanol concentration in blood
Temp = 310 # 310 Kelvin
particleM = (7.65 * 10 ** (-26)) # molecular mass of nicotine in kg
particleR = 2.33 * 10 ** (-10) # ethanol is a two angstrom radius sphere, yolo
eta = 1.0 * 10 ** (-6) # viscosity of water in m^2 / s
wallParam = [1,1] #[Hole Size (side), hole to wall ratio]
gamma = 6*np.pi*eta*particleR
timeStep = 1e-13
totalSteps = 10000000 #int(totalTime / timeStep)
holeDensityVol = 1.33e-29 # total volume of a single hole in m^3.
#Generate ICs, get some statistics, run random walk, more statistics
particleX0,particleV0 = genIC(boxSize,kB,NP,particleM,Temp)
totalEnergy(kB,particleM,particleV0,Temp)
#getState(boxSize, particleV0, particleX0, True)
lefts,particleX,particleV = randomWalk(boxSize,eta,holeDensityVol,kB,NP,particleM,particleR,particleV0,particleX0,pTrans,Temp,timeStep,totalSteps,wallParam)

TSPlot(lefts,timeStep)
totalEnergy(kB,particleM,particleV,Temp)
np.savetxt("pTrans = " + str(pTrans) + ".csv", lefts, delimiter='\n')

#getState(boxSize, particleV, particleX, True)


def beforeAfter(particleX0,particleV0,particleX,particleV,boxSize):
        normalizedX0 = np.copy(particleX0) 
        normalizedV0 = np.copy(particleV0)
        normalizedX = np.copy(particleX) 
        normalizedV = np.copy(particleV)
        
        fig = plt.figure(figsize = plt.figaspect(2.))
        ax1 = fig.add_subplot(2,1,1,projection='3d')
        plt.axis([0,boxSize,0,boxSize])
        ax1.set_zbound(lower=0, upper=boxSize)
        ax1.quiver(normalizedX0[:,0],normalizedX0[:,1],normalizedX0[:,2],normalizedV0[:,0],normalizedV0[:,1],normalizedV0[:,2], length = 0.05 * boxSize)
        ax1.set_xlabel('X axis (normal to membrane)')
        ax1.set_ylabel('Y axis (impenetrable walls)')
        ax1.set_zlabel('Z axis (periodic boundary )')

        ax2 = fig.add_subplot(2,1,2,projection='3d')
        plt.axis([0,boxSize,0,boxSize])
        ax2.set_zbound(lower=0, upper=boxSize)
        ax2.quiver(normalizedX[:,0],normalizedX[:,1],normalizedX[:,2],normalizedV[:,0],normalizedV[:,1],normalizedV[:,2], length = 0.05 * boxSize)
        ax2.set_xlabel('X axis (normal to membrane)')
        ax2.set_ylabel('Y axis (impenetrable walls)')
        ax2.set_zlabel('Z axis (periodic boundary )')
        plt.show()
    
#if __name__ == "__main__":
#    main()

