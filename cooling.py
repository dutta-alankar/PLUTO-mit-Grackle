# -*- coding: utf-8 -*-
"""
Created on Thu May  15 22:11:14 2025

@author: alankar
"""
import numpy as np
import matplotlib.pyplot as plt
from plot_style import *

metals = [0.3, 0.5, 1.0, 3.0]

for i, metal in enumerate(metals):
    data = np.loadtxt(f"./outputs/analysis-Z{metal:.1f}.dat")
    plt.plot(data[1:,1], data[1:,3], color=colors[i], label=f"{metal:.1f}")

plt.xlabel(r"Temperature [K]")
plt.ylabel(r"$\Lambda = -\dot{e}/n_H^2$ [$\rm erg\ cm^3 \ s^{-1}$]")
plt.legend(loc="best", fancybox=True, title=r"$Z/Z_\odot$", framealpha=0.6, ncol=2)

plt.xlim(xmin=1.0e+04, xmax=2.1e+06)
plt.ylim(ymin=5e-24, ymax=3e-21) 

plt.xscale("log")
plt.yscale("log")
plt.savefig("cooling-function.png", transparent=True, bbox_inches="tight")
