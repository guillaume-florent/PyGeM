#!/usr/bin/env python
# coding: utf-8

r"""Hull FFD example with STL"""

# TODO
#   More complex lattice + find a way to visualize vertices
#   Hydrostatics
#     Since STEP handling seems buggy/touchy, use STL hydrostatics ?
#     Test with STEP first, would make things way easier (waterline package)
#     if it worked

from aocutils.analyze.bounds import stl_bounding_box

import pygem as pg
from pygem.utils import write_bounding_box

stl_start = "./example_hull_stl_syrf/syrf_wide_light_canoe_body.stl"
stl_mod = './example_hull_stl_syrf/syrf_wide_light_canoe_body_mod.stl'

# Bounding box
(xmin, xmax), (ymin, ymax), (zmin, zmax) = stl_bounding_box(stl_start)
print("X min : %f" % xmin)
print("X max : %f" % xmax)
print("Y min : %f" % ymin)
print("Y max : %f" % ymax)
print("Z min : %f" % zmin)
print("Z max : %f" % zmax)

# assert ymin == -ymax

# Initial file loading and visualization
stl_handler = pg.stlhandler.StlHandler()
mesh_points = stl_handler.parse(stl_start)

# Display the unmodified hull in 2 possible ways
# stl_handler.plot(plot_file='./example_hull_stl/SYSSER01_Z0WL.stl')
# stl_handler.show(show_file='./example_hull_stl/SYSSER01_Z0WL.stl')

# Create FFDParameters and move a point
params = pg.params.FFDParameters.null_morphing_box(
    n_control_points=[5, 2, 3],
    length_box=[xmax-xmin, ymax-ymin, zmax - zmin],
    origin_box=[xmin, ymin, zmin])

# Move point symmetrically, manual symmetry enforcement
# params.move_point(i=1, j=1, k=1, direction="Y", displacement=1)
# params.move_point(i=1, j=0, k=1, direction="Y", displacement=-1)

# Move point symmetrically, automatic symmetry enforcement
params.move_point(i=0, j=1, k=1, direction="Y", displacement=0.5, symmetry="XZ")

write_bounding_box(params, './example_hull_stl_syrf/params.vtk', write_deformed=False)
write_bounding_box(params, './example_hull_stl_syrf/params_deformed.vtk', write_deformed=True)

# Deform
free_form = pg.freeform.FFD(params, mesh_points)
free_form.perform()
new_mesh_points = free_form.modified_mesh_points
stl_handler.write(new_mesh_points, stl_mod)

# Visualize deformation
# stl_handler.plot(plot_file='./example_hull_stl/SYSSER01_Z0WL_mod.stl')
stl_handler.show(show_file=stl_mod)

# Hull hydrostatics changes
