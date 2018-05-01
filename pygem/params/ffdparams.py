"""
Utilities for reading and writing parameters files to perform FFD
geometrical morphing.
"""
try:
    import configparser as configparser
except ImportError:
    import ConfigParser as configparser
import os
import numpy as np
from OCC.Bnd import Bnd_Box
from OCC.BRepBndLib import brepbndlib_Add
from OCC.BRepMesh import BRepMesh_IncrementalMesh
import vtk
import pygem.affine as at


class FFDParameters(object):
    """
    Class that handles the Free Form Deformation parameters in terms of FFD
    bounding box and weight of the FFD control points.

    :param list n_control_points: number of control points in the x, y, and z
        direction. If not provided it is set to [2, 2, 2].

    :cvar numpy.ndarray length_box: dimension of the FFD bounding box, in the
        x, y and z direction (local coordinate system).
    :cvar numpy.ndarray origin_box: the x, y and z coordinates of the origin of
        the FFD bounding box.
    :cvar numpy.ndarray rot_angle: rotation angle around x, y and z axis of the
        FFD bounding box.
    :cvar numpy.ndarray n_control_points: the number of control points in the
        x, y, and z direction.
    :cvar numpy.ndarray array_mu_x: collects the displacements (weights) along
        x, normalized with the box lenght x.
    :cvar numpy.ndarray array_mu_y: collects the displacements (weights) along
        y, normalized with the box lenght y.
    :cvar numpy.ndarray array_mu_z: collects the displacements (weights) along
        z, normalized with the box lenght z.

    :Example: from file

    >>> import pygem.params as ffdp
    >>>
    >>> # Reading an existing file
    >>> params1 = ffdp.FFDParameters()
    >>> params1.read_parameters(
    >>>     filename='tests/test_datasets/parameters_test_ffd_identity.prm')
    >>> 
    >>> # Creating a default parameters file with the right dimensions (if the
    >>> # file does not exists it is created with that name). So it is possible
    >>> # to manually edit it and read it again.
    >>> params2 = ffdp.FFDParameters(n_control_points=[2, 3, 2])
    >>> params2.read_parameters(filename='parameters_test.prm')
    >>> 
    >>> # Creating bounding box of the given shape
    >>> from OCC.IGESControl import IGESControl_Reader
    >>> params3 = ffdp.FFDParameters()
    >>> reader = IGESControl_Reader()
    >>> reader.ReadFile('tests/test_datasets/test_pipe.igs')
    >>> reader.TransferRoots()
    >>> shape = reader.Shape()
    >>> params3.build_bounding_box(shape)

    .. note::
        Four vertex (non coplanar) are sufficient to uniquely identify a
        parallelepiped.
        If the four vertex are coplanar, an assert is thrown when
        affine_points_fit is used.

    """

    def __init__(self, n_control_points=None):
        self.conversion_unit = 1.

        self.lenght_box = np.array([1., 1., 1.])
        self.origin_box = np.array([0., 0., 0.])
        self.rot_angle = np.array([0., 0., 0.])

        if n_control_points is None:
            n_control_points = [2, 2, 2]
        self.n_control_points = np.array(n_control_points)

        self.array_mu_x = np.zeros(self.n_control_points)
        self.array_mu_y = np.zeros(self.n_control_points)
        self.array_mu_z = np.zeros(self.n_control_points)

    @property
    def psi_mapping(self):
        """
        Map from the physical domain to the reference domain.

        :rtype: numpy.ndarray
        """
        return np.diag(np.reciprocal(self.lenght_box))

    @property
    def inv_psi_mapping(self):
        """
        Map from the reference domain to the physical domain.

        :rtype: numpy.ndarray
        """
        return np.diag(self.lenght_box)

    @property
    def rotation_matrix(self):
        """
        The rotation matrix (according to rot_angle_x, rot_angle_y,
        rot_angle_z).

        :rtype: numpy.ndarray
        """
        return at.angles2matrix(
            np.radians(self.rot_angle[2]), np.radians(self.rot_angle[1]),
            np.radians(self.rot_angle[0]))

    @property
    def position_vertices(self):
        """
        The position of the vertices of the FFD bounding box.

        :rtype: numpy.ndarray
        """
        return self.origin_box + np.vstack([
            np.zeros((1, 3)),
            self.rotation_matrix.dot(np.diag(self.lenght_box)).T
        ])

    def read_parameters(self, filename='parameters.prm'):
        """
        Reads in the parameters file and fill the self structure.

        :param string filename: parameters file to be read in.
        """
        if not isinstance(filename, str):
            raise TypeError("filename must be a string")

        # Checks if the parameters file exists. If not it writes the default
        # class into filename.
        if not os.path.isfile(filename):
            self.write_parameters(filename)
            return

        config = configparser.RawConfigParser()
        config.read(filename)

        self.n_control_points[0] = config.getint('Box info',
                                                 'n control points x')
        self.n_control_points[1] = config.getint('Box info',
                                                 'n control points y')
        self.n_control_points[2] = config.getint('Box info',
                                                 'n control points z')

        self.lenght_box[0] = config.getfloat('Box info', 'box lenght x')
        self.lenght_box[1] = config.getfloat('Box info', 'box lenght y')
        self.lenght_box[2] = config.getfloat('Box info', 'box lenght z')

        self.origin_box[0] = config.getfloat('Box info', 'box origin x')
        self.origin_box[1] = config.getfloat('Box info', 'box origin y')
        self.origin_box[2] = config.getfloat('Box info', 'box origin z')

        self.rot_angle[0] = config.getfloat('Box info', 'rotation angle x')
        self.rot_angle[1] = config.getfloat('Box info', 'rotation angle y')
        self.rot_angle[2] = config.getfloat('Box info', 'rotation angle z')

        self.array_mu_x = np.zeros(self.n_control_points)
        self.array_mu_y = np.zeros(self.n_control_points)
        self.array_mu_z = np.zeros(self.n_control_points)

        mux = config.get('Parameters weights', 'parameter x')
        muy = config.get('Parameters weights', 'parameter y')
        muz = config.get('Parameters weights', 'parameter z')

        for line in mux.split('\n'):
            values = np.array(line.split())
            self.array_mu_x[tuple(map(int, values[0:3]))] = float(values[3])

        for line in muy.split('\n'):
            values = line.split()
            self.array_mu_y[tuple(map(int, values[0:3]))] = float(values[3])

        for line in muz.split('\n'):
            values = line.split()
            self.array_mu_z[tuple(map(int, values[0:3]))] = float(values[3])

    def write_parameters(self, filename='parameters.prm'):
        """
        This method writes a parameters file (.prm) called `filename` and fills
        it with all the parameters class members.

        :param string filename: parameters file to be written out.
        """
        if not isinstance(filename, str):
            raise TypeError("filename must be a string")

        output_string = ""
        output_string += '\n[Box info]\n'
        output_string += '# This section collects all the properties of the'
        output_string += ' FFD bounding box.\n'

        output_string += '\n# n control points indicates the number of control'
        output_string += ' points in each direction (x, y, z).\n'
        output_string += '# For example, to create a 2 x 3 x 2 grid, use the'
        output_string += ' following: n control points: 2, 3, 2\n'
        output_string += 'n control points x: ' + str(
            self.n_control_points[0]) + '\n'
        output_string += 'n control points y: ' + str(
            self.n_control_points[1]) + '\n'
        output_string += 'n control points z: ' + str(
            self.n_control_points[2]) + '\n'

        output_string += '\n# box lenght indicates the length of the FFD '
        output_string += 'bounding box along the three canonical directions '
        output_string += '(x, y, z).\n'

        output_string += '# It uses the local coordinate system.\n'
        output_string += '# For example to create a 2 x 1.5 x 3 meters box '
        output_string += 'use the following: lenght box: 2.0, 1.5, 3.0\n'

        output_string += 'box lenght x: ' + str(self.lenght_box[0]) + '\n'
        output_string += 'box lenght y: ' + str(self.lenght_box[1]) + '\n'
        output_string += 'box lenght z: ' + str(self.lenght_box[2]) + '\n'

        output_string += '\n# box origin indicates the x, y, and z coordinates '
        output_string += 'of the origin of the FFD bounding box. That is '
        output_string += 'center of\n'

        output_string += '# rotation of the bounding box. It corresponds to '
        output_string += 'the point coordinates with position [0][0][0].\n'

        output_string += '# See section "Parameters weights" for more '
        output_string += 'details.\n'
        output_string += '# For example, if the origin is equal to 0., 0., 0., '
        output_string += 'use the following: origin box: 0., 0., 0.\n'

        output_string += 'box origin x: ' + str(self.origin_box[0]) + '\n'
        output_string += 'box origin y: ' + str(self.origin_box[1]) + '\n'
        output_string += 'box origin z: ' + str(self.origin_box[2]) + '\n'

        output_string += '\n# rotation angle indicates the rotation angle '
        output_string += 'around the x, y, and z axis of the FFD bounding box '
        output_string += 'in degrees.\n'

        output_string += '# The rotation is done with respect to the box '
        output_string += 'origin.\n'
        output_string += '# For example, to rotate the box by 2 deg along '
        output_string += 'the z '
        output_string += 'direction, use the following: rotation angle: '
        output_string += '0., 0., 2.\n'

        output_string += 'rotation angle x: ' + str(self.rot_angle[0]) + '\n'
        output_string += 'rotation angle y: ' + str(self.rot_angle[1]) + '\n'
        output_string += 'rotation angle z: ' + str(self.rot_angle[2]) + '\n'

        output_string += '\n\n[Parameters weights]\n'
        output_string += '# This section describes the weights of the FFD '
        output_string += 'control points.\n'

        output_string += '# We adopt the following convention:\n'
        output_string += '# For example with a 2x2x2 grid of control points we '
        output_string += 'have to fill a 2x2x2 matrix of weights.\n'

        output_string += '# If a weight is equal to zero you can discard the '
        output_string += 'line since the default is zero.\n'

        output_string += '#\n'
        output_string += '# | x index | y index | z index | weight |\n'
        output_string += '#  --------------------------------------\n'
        output_string += '# |    0    |    0    |    0    |  1.0   |\n'
        output_string += '# |    0    |    1    |    1    |  0.0   | --> you '
        output_string += 'can erase this line without effects\n'
        output_string += '# |    0    |    1    |    0    | -2.1   |\n'
        output_string += '# |    0    |    0    |    1    |  3.4   |\n'

        output_string += '\n# parameter x collects the displacements along x, '
        output_string += 'normalized with the box lenght x.'

        output_string += '\nparameter x:'
        offset = 1
        for i in range(0, self.n_control_points[0]):
            for j in range(0, self.n_control_points[1]):
                for k in range(0, self.n_control_points[2]):
                    output_string += offset * ' ' + str(i) + '   ' + str(
                        j) + '   ' + str(k) + '   ' + str(
                            self.array_mu_x[i][j][k]) + '\n'
                    offset = 13

        output_string += '\n# parameter y collects the displacements along y, '
        output_string += 'normalized with the box lenght y.'

        output_string += '\nparameter y:'
        offset = 1
        for i in range(0, self.n_control_points[0]):
            for j in range(0, self.n_control_points[1]):
                for k in range(0, self.n_control_points[2]):
                    output_string += offset * ' ' + str(i) + '   ' + str(
                        j) + '   ' + str(k) + '   ' + str(
                            self.array_mu_y[i][j][k]) + '\n'
                    offset = 13

        output_string += '\n# parameter z collects the displacements along z, '
        output_string += 'normalized with the box lenght z.'

        output_string += '\nparameter z:'
        offset = 1
        for i in range(0, self.n_control_points[0]):
            for j in range(0, self.n_control_points[1]):
                for k in range(0, self.n_control_points[2]):
                    output_string += offset * ' ' + str(i) + '   ' + str(
                        j) + '   ' + str(k) + '   ' + str(
                            self.array_mu_z[i][j][k]) + '\n'
                    offset = 13

        with open(filename, 'w') as f:
            f.write(output_string)

    def __str__(self):
        """
        This method prints all the FFD parameters on the screen. Its purpose is
        for debugging.
        """
        string = ""
        string += 'conversion_unit = {}\n'.format(self.conversion_unit)
        string += 'n_control_points = {}\n\n'.format(self.n_control_points)
        string += 'lenght_box = {}\n'.format(self.lenght_box)
        string += 'origin_box = {}\n'.format(self.origin_box)
        string += 'rot_angle  = {}\n'.format(self.rot_angle)
        string += '\narray_mu_x =\n{}\n'.format(self.array_mu_x)
        string += '\narray_mu_y =\n{}\n'.format(self.array_mu_y)
        string += '\narray_mu_z =\n{}\n'.format(self.array_mu_z)
        string += '\npsi_mapping = \n{}\n'.format(self.psi_mapping)
        string += '\nrotation_matrix = \n{}\n'.format(self.rotation_matrix)
        string += '\nposition_vertices = {}\n'.format(self.position_vertices)
        return string

    def save_points(self, filename, write_deformed=True):
        """
        Method that writes a vtk file containing the FFD lattice. This method
        allows to visualize where the FFD control points are located before the
        geometrical morphing. If the `write_deformed` flag is set to True the
        method writes out the deformed lattice, otherwise it writes one the
        original undeformed lattice.

        :param str filename: name of the output file.
        :param bool write_deformed: flag to write the original or modified FFD
            control lattice.  The default is set to True.

        :Example:

        >>> from pygem.params import FFDParameters
        >>> 
        >>> params = FFDParameters()
        >>> params.read_parameters(
        >>>     filename='tests/test_datasets/parameters_test_ffd_sphere.prm')
        >>> params.save('tests/test_datasets/box_test_sphere.vtk')

        .. note::
            In order to visualize the points in Paraview, please select the
            **Point Gaussian** representation.

        """
        x = np.linspace(0, self.lenght_box[0], self.n_control_points[0])
        y = np.linspace(0, self.lenght_box[1], self.n_control_points[1])
        z = np.linspace(0, self.lenght_box[2], self.n_control_points[2])

        lattice_y_coords, lattice_x_coords, lattice_z_coords = np.meshgrid(y, x,
                                                                           z)

        if write_deformed:
            box_points = np.array([
                lattice_x_coords.ravel() + self.array_mu_x.ravel() *
                self.lenght_box[0], lattice_y_coords.ravel() +
                self.array_mu_y.ravel() * self.lenght_box[1],
                lattice_z_coords.ravel() + self.array_mu_z.ravel() *
                self.lenght_box[2]
            ])
        else:
            box_points = np.array([
                lattice_x_coords.ravel(), lattice_y_coords.ravel(),
                lattice_z_coords.ravel()
            ])

        n_rows = box_points.shape[1]

        box_points = np.dot(
            self.rotation_matrix,
            box_points) + np.transpose(np.tile(self.origin_box, (n_rows, 1)))

        points = vtk.vtkPoints()

        for box_point in box_points.T:
            points.InsertNextPoint(box_point[0], box_point[1], box_point[2])

        data = vtk.vtkPolyData()
        data.SetPoints(points)

        writer = vtk.vtkPolyDataWriter()
        writer.SetFileName(filename)
        writer.SetInputData(data)
        writer.Write()

    def build_bounding_box(self,
                           shape,
                           tol=1e-6,
                           triangulate=False,
                           triangulate_tol=1e-1):
        """
        Builds a bounding box around the given shape. All parameters are set to
        match the computed box, the deformed FFD points are reset.

        :param shape: the shape to compute the bounding box.
        :type shape: TopoDS_Shape or its subclass
        :param float tol: tolerance of the computed bounding box.
        :param bool triangulate: if True, shape is triangulated before the
            bouning box creation.
        :param float triangulate_tol: tolerance of triangulation (size of
            created triangles).

        .. note::

            Every UV-Surface has to be rectangular. When a solid is created
            surfaces are trimmed. The trimmed part, however, is still saved
            inside a file. It is just *invisible* when drawn in a program.
        """
        bbox = Bnd_Box()
        bbox.SetGap(tol)
        if triangulate:
            BRepMesh_IncrementalMesh(shape, triangulate_tol)
        brepbndlib_Add(shape, bbox, triangulate)
        xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()
        min_xyz = np.array([xmin, ymin, zmin])
        max_xyz = np.array([xmax, ymax, zmax])

        self.origin_box = min_xyz
        self.lenght_box = max_xyz - min_xyz
        self.reset_deformation()

    def reset_deformation(self):
        """
        Set transformation parameters to arrays of zeros.
        """
        self.array_mu_x.fill(0.0)
        self.array_mu_y.fill(0.0)
        self.array_mu_z.fill(0.0)

    @classmethod
    def null_morphing_box(cls,
                          n_control_points,
                          length_box,
                          origin_box,
                          rot_angle=None):
        r"""Constructor of an FFD morphing box with zero displacements

        :param list n_control_points: number of control points
            in the x, y, and z direction.
        :param list length_box: length of the box
            in the x, y, and z direction.
        :param list origin_box: position of the box origin
            in the x, y, and z direction.
        :param list rot_angle: rotation angle
            around the x, y, and z axis. [0., 0., 0.] by default

        :return: FFDParameters with specified geometry and 0 displacements
        :rtype: FFDParameters

        """
        if rot_angle is None:
            rot_angle = [0., 0., 0.]
        p = cls(n_control_points)
        p.lenght_box = np.array(length_box)
        p.origin_box = np.array(origin_box)
        p.rot_angle = np.array(rot_angle)
        p.array_mu_x = np.zeros(n_control_points)
        p.array_mu_y = np.zeros(n_control_points)
        p.array_mu_z = np.zeros(n_control_points)

        # p.rotation_matrix = at.angles2matrix(radians(p.rot_angle[2]),
        #                                      radians(p.rot_angle[1]),
        #                                      radians(p.rot_angle[0]))

        p.position_vertex_0 = p.origin_box
        p.position_vertex_1 = p.position_vertex_0 + \
                              np.dot(p.rotation_matrix,
                                     [p.lenght_box[0], 0, 0])
        p.position_vertex_2 = p.position_vertex_0 + \
                              np.dot(p.rotation_matrix,
                                     [0, p.lenght_box[1], 0])
        p.position_vertex_3 = p.position_vertex_0 + \
                              np.dot(p.rotation_matrix,
                                     [0, 0, p.lenght_box[2]])

        # p.psi_mapping = np.diag(1. / p.lenght_box)
        # p.inv_psi_mapping = np.diag(p.lenght_box)

        return p

    def move_point(self, i, j, k, direction, displacement, symmetry=None):
        r"""Move point of index [i][j][k] by displacement in the specified
        direction

        :param int i: x index
        :param int j: y index
        :param int k: z index
        :param char direction: direction of movement (X Y or Z)
        :param float displacement: value of displacement (normed against box
            dimensions)
        :param str symmetry: symmetry enforcement parameter (None means no
            symmetry enforcement)

        """

        def symmetric_index(index, list_length):
            mid = (list_length - 1) / 2
            delta = index - mid
            return int(mid - delta)

        mu_dir = {"X": self.array_mu_x,
                  "Y": self.array_mu_y,
                  "Z": self.array_mu_z}

        mu_dir[direction][i][j][k] = displacement

        if symmetry is not None:
            if symmetry in ["XY", "YX"]:
                if direction == "Z" and k == symmetric_index(k,
                                                             self.n_control_points[
                                                                 2]):
                    raise ValueError("Cannot move the center point in "
                                     "opposite Z directions at the same time")
                mu_dir[direction][i][j][symmetric_index(k,
                                                        self.n_control_points[
                                                            2])] = -displacement if direction == "Z" else displacement
            elif symmetry in ["XZ", "ZX"]:
                if direction == "Y" and j == symmetric_index(j,
                                                             self.n_control_points[
                                                                 1]):
                    raise ValueError("Cannot move the center point in "
                                     "opposite Y directions at the same time")
                mu_dir[direction][i][
                    symmetric_index(j, self.n_control_points[1])][
                    k] = -displacement if direction == "Y" else displacement
            elif symmetry in ["YZ", "ZY"]:
                if direction == "X" and i == symmetric_index(i,
                                                             self.n_control_points[
                                                                 0]):
                    raise ValueError("Cannot move the center point in "
                                     "opposite X directions at the same time")
                mu_dir[direction][symmetric_index(i, self.n_control_points[0])][
                    j][k] = -displacement if direction == "X" else displacement
            else:
                raise ValueError("Unknown symmetry")
