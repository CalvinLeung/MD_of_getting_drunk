import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
timeStep = 1e-13 * 10000

def TSPlot(timeSeries,timeStep):
    N = np.arange(0,timeSeries.shape[0],1)
    fig = plt.figure()
    plt.logplot(N * timeStep, timeSeries)
    plt.show()
    return plt

lefts100 = genfromtxt('pTrans = 1.csv', delimiter='/n')
nLefts100 = lefts100[1:] / 1000

lefts050 = genfromtxt('pTrans = 0.50.csv', delimiter='/n')
nLefts050 = lefts050[1:]/ 1000

lefts025 = genfromtxt('pTrans = 0.25.csv', delimiter='/n')
nLefts025 = lefts025[1:] / 1000

lefts010 = genfromtxt('pTrans = 0.10.csv', delimiter='/n')
nLefts010 = lefts010[1:] / 1000

lefts001 = genfromtxt('pTrans = 0.01.csv', delimiter='/n')
nLefts001 = lefts001[1:] / 5000

N = np.arange(0,timeSeries.shape[0]-1,1)
fig = plt.figure()
plt.logplot(N * timeStep, lefts100,'-',
            N * timeStep, lefts050,'-',
            N * timeStep, lefts025,'-',
            N * timeStep, lefts010,'-',
            N * timeStep, lefts001,'-')
plt.show()
