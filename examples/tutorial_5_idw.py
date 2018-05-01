#!/usr/bin/env python
# coding: utf-8

r"""Tutorial 5 Inverse Distance Weighting interpolation technique on a cube"""

import numpy as np
import matplotlib.pyplot as plt

from pygem.params.idwparams import import IDWParameters
from pygem.idw import IDW


parameters_file = './tutorial_5_idw/parameters_idw_cube.prm'

params = IDWParameters()
params.read_parameters(filename=parameters_file)

nx, ny, nz = (10, 10, 10)
mesh = np.zeros((nx * ny * nz, 3))

xv = np.linspace(0, 1, nx)
yv = np.linspace(0, 1, ny)
zv = np.linspace(0, 1, nz)
z, y, x = np.meshgrid(zv, yv, xv)

mesh = np.array([x.ravel(), y.ravel(), z.ravel()])
mesh = mesh.T

fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')
ax.scatter(mesh[:, 0], mesh[:, 1], mesh[:, 2], c='blue', marker='o')
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')
plt.show()

idw = IDW(params, mesh)
idw.perform()
new_mesh_points = idw.modified_mesh_points

fig = plt.figure(2)
ax = fig.add_subplot(111, projection='3d')
ax.scatter(new_mesh_points[:, 0],
           new_mesh_points[:, 1],
           new_mesh_points[:, 2],
           c='red',
           marker='o')
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')
plt.show()
