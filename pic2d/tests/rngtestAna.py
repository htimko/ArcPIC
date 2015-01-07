#!/usr/bin/env python

import matplotlib.pyplot as plt
import os
import math as m
import numpy as np

def Gaus(v,mu,sigma):
    "Gaussian distribution"
    return np.exp(-0.5*((v-mu)/sigma)**2)/(sigma*m.sqrt(2*m.pi))

mu = 0.5;
sigma = 1.5;

print "reading file..."
gausFile = open("gausRandom.out", 'r')
gausNums = map(float,gausFile.read().split())
gausFile.close()

gausBins = np.linspace(mu-5*sigma, mu+5*sigma);
gausPoints = gausBins[:-1] - np.diff(gausBins)
plt.hist(gausNums,bins=gausBins,normed=True);
plt.plot(gausPoints, Gaus(gausPoints, mu, sigma))
plt.show()
