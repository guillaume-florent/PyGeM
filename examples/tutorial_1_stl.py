#!/usr/bin/env python
# coding: utf-8

r"""Tutorial 1 STL example"""

import pygem as pg
from pygem.utils import write_bounding_box

# Parameters that DO modify the shape
params = pg.params.FFDParameters()
params.read_parameters(filename='./tutorial_1_stl/parameters_test_ffd_sphere.prm')

# Create VTK files to compare the undeformed and deformed lattice in Paraview
write_bounding_box(params,
                   './tutorial_1_stl/params.vtk',
                   write_deformed=False)
write_bounding_box(params,
                   './tutorial_1_stl/params_deformed.vtk',
                   write_deformed=True)

# The params print themselves nicely
print(params)

# Parameters that DO NOT modify the shape
params_null = pg.params.FFDParameters()
params_null.read_parameters(filename='./tutorial_1_stl/null.prm')

# Read the STL file
stl_handler = pg.stlhandler.StlHandler()
mesh_points = stl_handler.parse('./tutorial_1_stl/test_sphere.stl')

# Display the unmodified sphere in 2 possible ways
stl_handler.plot(plot_file='./tutorial_1_stl/test_sphere.stl')
stl_handler.show(show_file='./tutorial_1_stl/test_sphere.stl')

# Apply the freeform transformation that does not change anything
free_form = pg.freeform.FFD(params_null, mesh_points)
free_form.perform()
new_mesh_points = free_form.modified_mesh_points
stl_handler.write(new_mesh_points, './tutorial_1_stl/test_sphere_null.stl')
stl_handler.plot(plot_file='./tutorial_1_stl/test_sphere_null.stl')

# Apply the freeform transformation that modifies the shape
free_form = pg.freeform.FFD(params, mesh_points)
free_form.perform()
new_mesh_points = free_form.modified_mesh_points
stl_handler.write(new_mesh_points, './tutorial_1_stl/test_sphere_mod.stl')

# Display the modified sphere in 2 possible ways
stl_handler.plot(plot_file='./tutorial_1_stl/test_sphere_mod.stl')
stl_handler.show(show_file='./tutorial_1_stl/test_sphere_mod.stl')
