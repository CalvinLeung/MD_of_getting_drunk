# Brownian Motion on the Unit Square
# Parameters
# wallParameter - wall spacing
# particleM - effective particle mass
# NP - particle number
# Temp - temperature which dictates thermal distribution (set k_B = 1)
def genIC(particleM,NP,T):
    particleX = np.random([NP,2])
    particleV = np.random([NP,2])
    return (particleX, particleV)

def plotState(particleX, particleV)
    import matplotlib.pyplot as plt
    import numpy as np
    from numpy import ma

    plt.figure()
    Q = plt.quiver(particleX,particleV)
    return 0


def evolveTimeStep(particleX, particleV, wallParameter,Temp)
    return particleV

def main():
    particleX, particleV = genIC(1,10,1)
    plotState(particleX,particleV)

if __name__ == "__main__":
    main()
