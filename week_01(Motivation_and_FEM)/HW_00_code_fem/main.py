# ---------------------------------------------------------------------
#
# Copyright (C) 2026 by Peter Munch
#
# This file is part of the class
# Advanced Topics in Numerics of PDEs: FEM for the Navier-Stokes Equations
# given at TU Berlin 2026.
#
# This is free software; you can use it, redistribute
# it, and/or modify it under the terms of the GNU Lesser General
# Public License as published by the Free Software Foundation; either
# version 3.0 of the License, or (at your option) any later version.
#
# ---------------------------------------------------------------------

from meshpy.triangle import MeshInfo, build
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
import scipy


# quadrature
def create_quadrature(dim, cell_type, n_q_points):
    if dim == 1:
        if n_q_points == 2:
            return [[[(-1/math.sqrt(3)+1)/2], 0.5], [[(+1/math.sqrt(3)+1)/2], 0.5]]
        elif n_q_points == 3:
            return [[[(-math.sqrt(3/5)+1)/2], 5/18], [[0.5], 4/9], [[(+math.sqrt(3/5)+1)/2], 5/18]]
        raise Exception('Not implemented!')
    elif dim == 2:
        if cell_type == "hyper_cube":
            quadrature = []

            raise Exception('TODO: Problem 1.4')

            return quadrature

        else:
            if n_q_points == 2:
                return [
                    [[0.17855872826361643, 0.1550510257216822],
                        0.31804138174397717 / 2],
                    [[0.07503111022260812, 0.6449489742783178],
                        0.18195861825602283 / 2],
                    [[0.6663902460147014, 0.1550510257216822],
                        0.31804138174397717 / 2],
                    [[0.28001991549907407, 0.6449489742783178],
                        0.18195861825602283 / 2]
                ]
            raise Exception('Not implemented!')

    raise Exception('Not implemented!')


# shape functions
def reference_shape_value(dim, fe_degree, cell_type, x_ref):
    if dim == 1:
        if fe_degree == 1:
            return np.array([1 - x_ref[0], x_ref[0]])
        elif fe_degree == 2:
            return np.array([2 * x_ref[0] ** 2 - 3 * x_ref[0] + 1, -4*x_ref[0]**2 + 4 * x_ref[0], 2 * x_ref[0]**2 - x_ref[0]])
        raise Exception('Not implemented!')
    elif dim == 2:
        if cell_type == "hyper_cube":
            if fe_degree == 1:
                raise Exception('TODO: Problem 1.4')
                return 0
        else:
            if fe_degree == 1:
                return np.array([1 - x_ref[0] - x_ref[1], x_ref[0], x_ref[1]])
        raise Exception('Not implemented!')

    raise Exception('Not implemented!')


def reference_shape_gradient(dim, fe_degree, cell_type, x_ref):
    if dim == 1:
        if fe_degree == 1:
            return np.array([[-1.0, +1.0]])
        elif fe_degree == 2:
            return np.array([[4*x_ref[0]-3, - 8 * x_ref[0] + 4, 4*x_ref[0]-1]])
        raise Exception('Not implemented!')
    elif dim == 2:
        if cell_type == "hyper_cube":
            if fe_degree == 1:
                raise Exception('TODO: Problem 1.4')
                return 0.0
        else:
            if fe_degree == 1:
                return np.array([[-1.0, +1.0, +0.0], [-1.0, +0.0, +1.0]])
        raise Exception('Not implemented!')

    raise Exception('Not implemented!')


# mapping
def compute_Jacobian(dim, ps, x_ref):
    if dim == 1:
        return np.array([[ps[1][0] - ps[0][0]]])
    elif dim == 2:

        if len(ps) == 4:
            raise Exception('TODO: Problem 1.4')
            return 0.0
        else:
            return np.array([[ps[1][0] - ps[0][0], ps[2][0] - ps[0][0]],
                             [ps[1][1] - ps[0][1], ps[2][1] - ps[0][1]]])

    raise Exception('Not implemented!')


def compute_Jacobian_determinant(dim, ps, x_ref):
    if dim == 1:
        return ps[1][0] - ps[0][0]
    elif dim == 2:
        if len(ps) == 4:
            raise Exception('TODO: Problem 1.4')
            return 0.0
        else:
            return (ps[1][0] - ps[0][0]) * (ps[2][1] - ps[0][1]) - (ps[2][0] - ps[0][0]) * (ps[1][1] - ps[0][1])

    raise Exception('Not implemented!')


def compute_quadrature_point(dim, ps, x_ref):
    return reference_shape_value(dim, 1, "simplex" if len(ps) == 3 else "hyper_cube", x_ref) @ ps


class my_build:
    def __init__(self, points, facets):
        self.points = points
        self.elements = facets


def hextet(mesh):
    element_tri = []

    for e in mesh.elements:
        element_tri += [[e[i] for i in [0, 1, 2]]]
        element_tri += [[e[i] for i in [1, 2, 3]]]

    return my_build(mesh.points, element_tri)


def tethex(mesh):
    dim = len(mesh.points[0])

    if dim != 2:
        raise Exception('Not implemented!')

    new_points = []
    for p in mesh.points:
        new_points = new_points + [p]

    for element in mesh.elements:
        new_point = [0] * dim

        for p in element:
            for d in range(dim):
                new_point[d] += mesh.points[p][d]

        for d in range(dim):
            new_point[d] /= len(element)

        new_points = new_points + [new_point]

    new_elements = []

    n_vertices = len(mesh.points)
    n_cells = len(mesh.elements)

    map = {}
    for c, element in enumerate(mesh.elements):
        indices = []

        for e in element:
            indices = indices + [e]

        for ii in ([[0, 1], [1, 2], [2, 0]] if dim == 2 else [[0], [fe_degree]]):
            face_points = tuple(sorted([element[i] for i in ii]))

            if face_points not in map:
                map[face_points] = len(map)

                new_point = [0] * dim

                for p in face_points:
                    for d in range(dim):
                        new_point[d] += mesh.points[p][d]

                for d in range(dim):
                    new_point[d] /= len(face_points)

                new_points = new_points + [new_point]

            indices = indices + \
                [map[face_points] + n_vertices + n_cells]

        indices = indices + [n_vertices + c]

        new_elements += [[indices[i] for i in [0, 3, 5, 6]]]
        new_elements += [[indices[i] for i in [3, 1, 6, 4]]]
        new_elements += [[indices[i] for i in [5, 6, 2, 4]]]

    return my_build(new_points, new_elements)


def plot_solution(dim, fe_degree, mesh, dofs, solution_vector, fu_g):
    if dim == 1:
        plt.figure(0)
        plt.ion()

        nodes_x = [p[0] for p in dofs.points]
        plt.plot(nodes_x, solution_vector, 'r.',
                 markersize=4, label='computed solution')
        plt.plot(nodes_x, [fu_g(p)
                 for p in dofs.points], 'k-', label='true solution')
        plt.xlabel('x')
        plt.ylabel('u')

        plt.legend()
        plt.grid(True, which='both')
        plt.show()
        plt.pause(5)
        plt.clf()

    elif dim == 2:
        if fe_degree != 1:
            raise Exception('Not implemented!')

        mesh_type = "simplex" if len(mesh.elements[0]) == 3 else "hyper_cube"

        plot_mesh = True

        plt.figure(0)
        plt.ion()

        nodes_x = [p[0] for p in mesh.points]
        nodes_y = [p[1] for p in mesh.points]

        if mesh_type == "hyper_cube":
            mesh_tri = hextet(mesh)
        else:
            mesh_tri = mesh

        triangulation = matplotlib.tri.Triangulation(
            nodes_x, nodes_y, mesh_tri.elements)
        plt.tricontourf(triangulation, solution_vector)

        if plot_mesh:
            for element in mesh.elements:
                for ii in ([[0, 1], [2, 3], [0, 2], [1, 3]] if mesh_type == "hyper_cube" else [[0, 1], [1, 2], [2, 0]]):
                    element_points = [mesh.points[i] for i in element]
                    plt.plot([element_points[i][0] for i in ii], [
                             element_points[i][1] for i in ii], 'k', linewidth=0.5)

        plt.gca().set_aspect('equal')
        plt.show()
        plt.pause(5)
        # plt.savefig("results/contour.pdf", format="pdf")
        plt.clf()

    else:
        raise Exception('Not implemented!')


def compute_error(dim, fe_degree, mesh, dofs, quadrature, solution_vector, fu):
    result = 0

    for cell, element in zip(mesh.elements, dofs.elements):

        n_dofs_per_cell = len(element)
        n_q_points = len(quadrature)

        ps = [mesh.points[i] for i in cell]
        detJ = np.zeros(n_q_points)

        shape_value = np.zeros([n_dofs_per_cell, n_q_points])

        for q, [x_ref, w] in enumerate(quadrature):
            detJ[q] = compute_Jacobian_determinant(dim, ps, x_ref)
            shape_value[:, q] = reference_shape_value(
                dim, fe_degree, "simplex" if len(ps) == 3 else "hyper_cube", x_ref)

        values = [solution_vector[i] for i in element]

        for q, [x_ref, w] in enumerate(quadrature):
            x_real = compute_quadrature_point(dim, ps, x_ref)
            value_true = fu(x_real)

            computed_value = 0.0
            for i in range(0, n_dofs_per_cell):
                computed_value += shape_value[i, q] * values[i]

            result += (value_true - computed_value) ** 2 * detJ[q] * w

    return math.sqrt(result)


# mesh handling
def get_boundary_dofs(dim, fe_degree, dofs):
    map = {}

    if dim == 2 and fe_degree > 1:
        raise Exception('Not implemented!')

    for element in dofs.elements:
        for ii in (([[0, 1], [1, 2], [2, 0]] if len(element) else [[0, 1], [2, 3], [0, 2], [1, 3]]) if dim == 2 else [[0], [fe_degree]]):
            face_points = tuple(sorted([element[i] for i in ii]))

            if face_points in map:
                map[face_points] = map[face_points] + 1
            else:
                map[face_points] = 1

    boundary_dofs = []

    for points, counter in map.items():
        if counter == 1:
            boundary_dofs += points

    return sorted(set(boundary_dofs))


def main():
    # simulation name
    simulation_name = "singularity"  # "gaussian|unit|singularity"

    # mesh setting
    dim = 2
    mesh_type = "simplex"
    n_cells = 20  # only used in 1D
    fe_degree = 1
    n_q_points_1D = fe_degree + 1

    mesh_name = "l"
    max_volume = (0.5/2**3)**2  # only used in 2D

    if simulation_name == "gaussian" and dim == 1:
        gaussian_width = 10

        def fu_f(x): return 2 * gaussian_width * (1 - 2 * gaussian_width *
                                                  (x[0] - 0.5)**2) * math.exp(-gaussian_width * (x[0]-0.5)**2)

        def fu_g(x): return math.exp(-gaussian_width * (x[0]-0.5)**2)

    elif simulation_name == "unit" and dim == 2:
        def fu_f(x): return 1.0
        def fu_g(x): return 1.0
    elif simulation_name == "singularity" and dim == 2:
        def fu_f(x): return 0.0

        def fu_g(x): return np.abs(complex(x[0], x[1]))**(2/3) * \
            np.sin(2/3 * (np.angle(complex(x[0], x[1])) % (2*math.pi)))
    else:
        raise Exception('Not implemented!')

    quadrature = create_quadrature(dim, mesh_type, n_q_points_1D)
    n_q_points = len(quadrature)

    # create mesh
    if dim == 1:
        mesh = my_build([[i/n_cells] for i in range(0, n_cells + 1)],
                        [[i, i + 1] for i in range(0, n_cells)])

        dofs = my_build([[i/n_cells/fe_degree] for i in range(0, n_cells*fe_degree + 1)],
                        [[i*fe_degree + j for j in range(0, fe_degree + 1)] for i in range(0, n_cells)])

    elif dim == 2:
        if mesh_name == "cube" or mesh_name == "l":
            if mesh_name == "cube":
                points = [(-1, -1), (1, -1), (1, 1), (-1, 1)]
                facets = [[0, 1], [1, 2], [2, 3], [3, 0]]
            elif mesh_name == "l":
                points = [(-1, -1), (0, -1), (0, 0), (1, 0), (1, 1), (-1, 1)]
                facets = [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]
            else:
                Exception('Not implemented!')

            mesh_info = MeshInfo()
            mesh_info.set_points(points)
            mesh_info.set_facets(facets)
            mesh = build(mesh_info, max_volume=max_volume)

            if mesh_type == "hyper_cube":
                mesh = tethex(mesh)

        elif mesh_name == "structured cube" or mesh_name == "structured l":

            points = []
            elements = []

            for iy in range(n_cells + 1):
                for ix in range(n_cells + 1):
                    points += [[-1.0 + 2.0 / n_cells *
                                ix, -1.0 + 2.0 / n_cells * iy]]

            for iy in range(n_cells):
                for ix in range(n_cells):
                    elements += [[(ix+0)+(iy+0)*(n_cells+1),
                                  (ix+1)+(iy+0)*(n_cells+1),
                                  (ix+0)+(iy+1)*(n_cells+1),
                                  (ix+1)+(iy+1)*(n_cells+1)]]

            mesh = my_build(points, elements)

            if mesh_name == "structured l":

                points = []
                elements = []

                raise Exception('TODO: Problem 1.4')

                mesh = my_build(points, elements)

        if fe_degree != 1:
            Exception('Not implemented!')

        dofs = my_build(mesh.points, mesh.elements)

    boundary_dofs = get_boundary_dofs(dim, fe_degree, dofs)

    # allocate memory
    rhs_vector = np.zeros(len(dofs.points))

    col_indices = []
    row_indices = []

    for element in dofs.elements:
        for j in element:
            for i in element:
                col_indices += [i]
                row_indices += [j]

    A = scipy.sparse.csr_matrix(
        ([0.0]*len(col_indices), (col_indices, row_indices)))

    # compute system matrix and right-hand-side vector
    for cell, element in zip(mesh.elements, dofs.elements):

        n_dofs_per_cell = len(element)

        # compute quadrature points, Jacobian and its determinant
        ps = [mesh.points[i] for i in cell]
        detJ = np.zeros(n_q_points)

        # evaluate the value and gradient of the
        # shape functions at the quadrature points
        shape_value = np.zeros([n_dofs_per_cell, n_q_points])
        shape_gradient = np.zeros([n_dofs_per_cell, n_q_points, dim])

        for q, [x_ref, w] in enumerate(quadrature):
            J = np.transpose(np.linalg.inv(compute_Jacobian(dim, ps, x_ref)))
            detJ[q] = compute_Jacobian_determinant(dim, ps, x_ref)

            shape_value[:, q] = reference_shape_value(
                dim, fe_degree, mesh_type, x_ref)
            shape_gradient[:, q, :] = np.transpose(
                J @ reference_shape_gradient(dim, fe_degree, mesh_type, x_ref))

        A_cell = np.zeros([n_dofs_per_cell, n_dofs_per_cell])
        rhs_cell = np.zeros([n_dofs_per_cell])

        for q, [x_ref, w] in enumerate(quadrature):
            value = fu_f(compute_quadrature_point(dim, ps, x_ref))

            for i in range(0, n_dofs_per_cell):
                rhs_cell[i] += shape_value[i, q] * value * detJ[q] * w

            for i in range(0, n_dofs_per_cell):
                for j in range(0, n_dofs_per_cell):
                    A_cell[i, j] += np.dot(shape_gradient[i, q, :],
                                           shape_gradient[j, q, :]) * detJ[q] * w

        # apply boundary condition
        for e, i in enumerate(element):
            if i in boundary_dofs:
                rhs_cell[e] = 0

                for j in range(0, len(element)):
                    A_cell[e, j] = 0

                    if not element[j] in boundary_dofs:
                        rhs_cell[j] -= A_cell[j, e] * fu_g(dofs.points[i])

                    A_cell[j, e] = 0
                A_cell[e, e] = 1

        # assemble vector and matrix
        for e, i in enumerate(element):
            rhs_vector[i] += rhs_cell[e]

        for ej, j in enumerate(element):
            for ei, i in enumerate(element):
                A[i, j] += A_cell[ei, ej]

    solution_vector = scipy.sparse.linalg.spsolve(A, rhs_vector)

    # set boundary conditions
    for e, p in enumerate(dofs.points):
        if e in boundary_dofs:
            solution_vector[e] = fu_g(p)

    # computer error
    print(" - error: %f %f" %
          (math.sqrt(max_volume) if dim == 2 else (1.0/n_cells), compute_error(dim, fe_degree, mesh, dofs, quadrature, solution_vector, fu_g)))

    # plot solution
    plot_solution(dim, fe_degree, mesh, dofs, solution_vector, fu_g)


if __name__ == '__main__':
    main()
