# -*- coding: utf-8 -*-
"""
Created on Thu May  15 22:11:14 2025

@author: alankar
"""
import numpy as np
import matplotlib.pyplot as plt
from plot_style import *

data = np.loadtxt("./outputs/analysis.dat")

plt.plot(data[1:,1], data[1:,3], color=colors[1])

plt.xlabel(r"Temperature [K]")
plt.ylabel(r"$\Lambda = -\dot{e}/n_H^2$ [$\rm erg\ cm^3 \ s^{-1}$]")

# plt.xlim(xmin=-0.1, xmax=25)
# plt.ylim(ymin=2.0e-02, ymax=14.0) 

plt.xscale("log")
plt.yscale("log")
plt.savefig("cooling-function.png", transparent=True, bbox_inches="tight")
