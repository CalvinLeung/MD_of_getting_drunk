# Brownian Motion on the Unit Square
# Parameters
# wallParameter - wall spacing
# particleM - effective particle mass
# NP - particle number
# Temp - temperature which dictates thermal distribution (set k_B = 1)
def genIC(particleM,NP,Temp):
    import numpy.random as npr
    particleX = npr.rand(NP,2)
    particleX[:,0] = 0.5 * particleX[:,0]
    particleV = npr.exponential(1/(2 * particleM * Temp),(NP,2))
    return (particleX, particleV)

def plotState(particleX, particleV):
    import matplotlib.pyplot as plt
    import numpy as np
    from numpy import ma

    plt.figure()
    Q = plt.quiver(particleX[:,0],particleX[:,1],particleV[:,0],particleV[:,1])
    plt.show()
    return 0


def evolveTimeStep(particleX, particleV, wallParameter,Temp):
    return particleV

def main():
    particleX, particleV = genIC(1,10,1)
    plotState(particleX,particleV)

if __name__ == "__main__":
    main()
