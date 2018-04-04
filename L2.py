#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 13:22:41 2018

@author: bartlomiejkos
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import mpl_toolkits.mplot3d.axes3d as axes3d
import pandas as pd

# the scopes need to be changed accordingly to provided n & f in C++ program
# Ex. 2a
data = pd.read_csv("/Users/bartlomiejkos/Downloads/Programming/MetodyNumeryczneII/data.csv", header=None)
data = data.iloc[0:100, 0:100]  # the scope depends on f - printing frequency e.g. if f = 10 then data.iloc[0:100]...
fig = plt.figure()                # plot data from C++ program
ax = fig.add_subplot(111, projection='3d')
Xs = np.linspace(0, 1, 100)
Ys = np.linspace(0, 1, 100)
Xs, Ys = np.meshgrid(Xs, Ys)
data = pd.DataFrame(data)
ax.plot_surface(Xs, Ys, data, rstride=4, cstride=4, alpha=0.4,cmap=cm.jet)

fig = plt.figure()                # calculate and plot explicitly
ax = fig.add_subplot(111, projection='3d')
Xs = np.arange(0, 1, 0.01)
Ys = np.arange(0, 1, 0.01)
Xs, Ys = np.meshgrid(Xs, Ys)
s = (100, 100)
u = np.zeros(s)
for m in range(0, 100):
    for x in range(0, 100):
        u[m][x] = 4 * m * x * (m**2 - x**2)
u = pd.DataFrame(u)
ax.plot_surface(Xs, Ys, u, rstride=4, cstride=4, alpha=0.4,cmap=cm.jet)

# Ex. 2b
data = pd.read_csv("/Users/bartlomiejkos/Downloads/Programming/MetodyNumeryczneII/data2.csv", header=None)
data = data.iloc[0:100, 0:100]  # the scope depends on f - printing frequency e.g. if f = 10 then data.iloc[0:100]...
fig = plt.figure()                # plot data from C++ program
ax = fig.add_subplot(111, projection='3d')
Xs = np.arange(0, 1, 0.01)
Ys = np.arange(0, 1, 0.01)
Xs, Ys = np.meshgrid(Xs, Ys)
data = pd.DataFrame(data)
ax.plot_surface(Xs, Ys, data, rstride=4, cstride=4, alpha=0.4,cmap=cm.jet)

fig = plt.figure()                  # calculate and plot explicitly
ax = fig.add_subplot(111, projection='3d')
Xs = np.arange(0, 1, 0.01)
Ys = np.arange(0, 1, 0.01)
Xs, Ys = np.meshgrid(Xs, Ys)
s = (100, 100)
u = np.zeros(s)
for m in range(0, 100):
    for x in range(0, 100):
        u[m][x] = 1e-4 * np.sin(np.pi * 3 * x) * np.sin(np.pi * 3 * m)
u = pd.DataFrame(u)
ax.plot_surface(Xs, Ys, u, rstride=4, cstride=4, alpha=0.4,cmap=cm.jet)
