#!/usr/bin/env python
# coding: utf-8

r"""Hull FFD example with IGES"""

from aocxchange.iges import IgesImporter
from aocutils.analyze.bounds import BoundingBox

import pygem as pg
from pygem.utils import write_bounding_box

# Bounding box
hull_shape_compound = \
    IgesImporter(filename="./example_hull_iges/SYSSER01_Z0WL.igs").compound
bb = BoundingBox(hull_shape_compound)

print("X min : %f" % bb.x_min)
print("X max : %f" % bb.x_max)
print("Y min : %f" % bb.y_min)
print("Y max : %f" % bb.y_max)
print("Z min : %f" % bb.z_min)
print("Z max : %f" % bb.z_max)

assert abs(abs(bb.y_min) - abs(bb.y_max)) <= 1e-6

# Initial file loading and visualization
iges_handler = pg.igeshandler.IgesHandler()
mesh_points = iges_handler.parse('./example_hull_iges/SYSSER01_Z0WL.igs')
iges_handler.check_topology()

# Display the unmodified hull in 2 possible ways
# iges_handler.plot(plot_file='./example_hull_iges/SYSSER01_Z0WL.igs')
iges_handler.show(show_file='./example_hull_iges/SYSSER01_Z0WL.igs')

# Create FFDParameters and move a point
params = pg.params.FFDParameters.null_morphing_box(
    n_control_points=[3, 2, 2],
    length_box=[bb.x_max-bb.x_min, bb.y_max-bb.y_min, bb.z_max - bb.z_min],
    origin_box=[bb.x_min, bb.y_min, bb.z_min])

# Move point symmetrically, manual symmetry enforcement
# params.move_point(i=1, j=1, k=1, direction="Y", displacement=1)
# params.move_point(i=1, j=0, k=1, direction="Y", displacement=-1)

# Move point symmetrically, automatic symmetry enforcement
params.move_point(i=1, j=1, k=1, direction="Y", displacement=1, symmetry="XZ")

write_bounding_box(params,
                   './example_hull_iges/params.vtk',
                   write_deformed=False)
write_bounding_box(params,
                   './example_hull_iges/params_deformed.vtk',
                   write_deformed=True)

# Deform
free_form = pg.freeform.FFD(params, mesh_points)
free_form.perform()
new_mesh_points = free_form.modified_mesh_points
iges_handler.write(new_mesh_points, './example_hull_iges/SYSSER01_Z0WL_mod.igs')

# Visualize deformation
# iges_handler.plot(plot_file='./example_hull_iges/SYSSER01_Z0WL_mod.igs')
iges_handler.show(show_file='./example_hull_iges/SYSSER01_Z0WL_mod.igs')
