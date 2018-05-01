#!/usr/bin/env python
# coding: utf-8

r"""Tutorial 2 IGES example"""

import pygem as pg

params = pg.params.FFDParameters()
params.read_parameters(filename='./tutorial_2_iges/parameters_test_ffd_iges.prm')

iges_handler = pg.igeshandler.IgesHandler()
mesh_points = iges_handler.parse('./tutorial_2_iges/test_pipe.iges')
iges_handler.check_topology()

iges_handler.show('./tutorial_2_iges/test_pipe.iges')

free_form = pg.freeform.FFD(params, mesh_points)

free_form.perform()

new_mesh_points = free_form.modified_mesh_points

# iges_handler.write(new_mesh_points, './tutorial_2_iges/test_pipe_mod.iges')
iges_handler.write(new_mesh_points,
                   './tutorial_2_iges/test_pipe_mod.iges',
                   1e-3)

iges_handler.show('./tutorial_2_iges/test_pipe_mod.iges')
