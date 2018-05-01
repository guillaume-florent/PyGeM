#!/usr/bin/env python
# coding: utf-8

r"""STEP morphing example"""

import logging
import pygem as pg

logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s :: %(levelname)6s :: '
                              '%(module)20s :: %(lineno)3d :: %(message)s')
handler.setFormatter(formatter)
handler.setLevel(logging.DEBUG)
logger.addHandler(handler)

logger.info("reading parameters")
params = pg.params.FFDParameters()
params.read_parameters(filename='./example_ffd_step/parameters_ffd_step.prm')

step_handler = pg.stephandler.StepHandler()

mesh_points = step_handler.parse('./example_ffd_step/test_pipe.step')

logger.info("show original")
step_handler.show('./example_ffd_step/test_pipe.step')

logger.info("applying transformation")
free_form = pg.freeform.FFD(params, mesh_points)
free_form.perform()

new_mesh_points = free_form.modified_mesh_points

# iges_handler.write(new_mesh_points, './tutorial_2_iges/test_pipe_mod.iges')
logger.info("writing modified file")
step_handler.write(new_mesh_points,
                   './example_ffd_step/test_pipe_mod.step',
                   1e-3)

logger.info("showing modified file")
mesh_points = step_handler.parse('./example_ffd_step/test_pipe_mod.step')
step_handler.show()
