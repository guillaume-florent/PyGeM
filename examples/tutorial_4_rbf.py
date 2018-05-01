#!/usr/bin/env python
# coding: utf-8

r"""Tutorial 4 RBF interpolation technique on a cube"""

import pygem as pg
import numpy as np
import matplotlib.pyplot as plt

from pygem.utils import plot_rbf_control_points, write_points_in_vtp


params = pg.params.RBFParameters()
params.read_parameters('./tutorial_4_rbf/parameters_rbf_cube.prm')

# create a 10-by-10-by-10 lattice to mimic a cube
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
ax.set_zlabel('Z Axis')

plt.show()

plot_rbf_control_points(params, save_fig=False)

write_points_in_vtp(params.deformed_control_points,
                    outfile='./tutorial_4_rbf/points.vtp',
                    color=(0, 255, 255))

rbf = pg.radial.RBF(params, mesh)
rbf.perform()
new_mesh_points = rbf.modified_mesh_points

fig = plt.figure(2)

ax = fig.add_subplot(111, projection='3d')
ax.scatter(new_mesh_points[:, 0],
           new_mesh_points[:, 1],
           new_mesh_points[:, 2],
           c='red', marker='o')
ax.set_xlabel('X Axis')
ax.set_ylabel('Y Axis')
ax.set_zlabel('Z Axis')

plt.show()
