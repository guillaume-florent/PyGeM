#!/usr/bin/env python
# coding: utf-8

r"""Tutorial 3 UNV example"""

import pygem as pg

params = pg.params.FFDParameters()
params.read_parameters(filename='./tutorial_3_unv/parameters_test_ffd_pipe_unv_C0.prm')

unv_handler = pg.unvhandler.UnvHandler()
mesh_points = unv_handler.parse('./tutorial_3_unv/test_pipe.unv')

free_form = pg.freeform.FFD(params, mesh_points)
free_form.perform()
new_mesh_points = free_form.original_mesh_points

unv_handler.write(new_mesh_points, './tutorial_3_unv/test_pipe_mod_C0.unv')

# ----

params = pg.params.FFDParameters()
params.read_parameters(filename='./tutorial_3_unv/parameters_test_ffd_pipe_unv_C1.prm')

free_form = pg.freeform.FFD(params, mesh_points)
free_form.perform()
new_mesh_points = free_form.modified_mesh_points

unv_handler.write(new_mesh_points, './tutorial_3_unv/test_pipe_mod_C1.unv')
