# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 20:21:18 2022

@author: Amalia Karalis
"""

import numpy as np
import matplotlib.pyplot as plt

dt = 0.1 # choose an appropriate timestep
Nsteps = 100 # choose an appropriate number of total timesteps

## Setting up initial conditions (vortex centers and circulations)
# Vortex rings
y_v = np.array([10, 30, 50, 70]) # the y positions of the 4 vortices
x_v = np.array([10, 30, 50, 70]) # the x positions of the 4 vortices
k_v = np.array([1, 1, 1, 1]) # the line vortex constant k of the 4 vortices

# Setting up the plot
plt.ion()
fig, ax = plt.subplots(1,1)
# mark the initial positions of vortices
p, = ax.plot(x_v, y_v, 'k+', markersize=10)

#draw the initial velocity streamline
n = 100
ngrid = n # the dimensions of the simulation grid
Y, X = np.mgrid[-ngrid:ngrid:360j, -ngrid:ngrid:360j] #360j sets the resolution of the cartesian grid
vel_x = np.zeros(np.shape(X)) # this holds the x-velocity
vel_y = np.zeros(np.shape(Y)) # this holds the y-velocity

# masking radius for better visualization of the vortex centers
r_mask = 10 # the radius for the mask around the vortex centers
# within this mask, the streamlines are not plotted
# so that the movement of the vortex centers can be seen more clearly

for i in range(len(x_v)): # looping over each vortex
    # computing the total velocity field
    # u = k / r phi-hat 
    # loop through each point in the grid
    for j in range(n):
        for k in range(n):
            # j is the row, k is the column
            r = np.sqrt((x_v[i] - j)**2 + (y_v[i] - k)**2) # the distance from the center of the vortex to the point
            u = k_v[i] / r # the magnitude of the velocity at that location
            v_x, v_y = u*np.array([j - x_v[i], k - y_v[i]]) / r # the vector velocity
            vel_x[j, k] += v_x # the x-coordinate of the velocity at that location
            vel_y[j, k] += v_y  # the y-coordinate of the velocity at that location
    # set the masking area to NaN
    
# set up the boundaries for of the simulation box
ax.set_xlim([-ngrid, ngrid])
ax.set_ylim([-ngrid, ngrid])

# initial plot of the streamlines
ax.streamplot(X, Y, vel_x, vel_y, density=[1,1])

fig.canvas.draw()

# Evolution
count = 0
while count < Nsteps:
    ## Compute the advection velocity due to each vortex
    for i in range(len(x_v)):
        x = x_v[i] # current position of the vortex
        y = y_v[i]
        x_v[i] += vel_x[x, y] * dt # update the position of the vortices
        y_v[i] += vel_y[x, y] * dt 
        
    # re-initialize the total velocity field
    vel_x = np.zeros_like(vel_x)
    vel_y = np.zeros_like(vel_y)
    








