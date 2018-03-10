"""
Copyright 2014-2018, David Bernstein

This file is part of Graph.

Graph is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Graph is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Graph.  If not, see <http://www.gnu.org/licenses/>.
"""

from dolfin import *
from math import *

class Graph:
    def __init__(self):
        self.edge_list = []
        self.num_vertices = 0
        self.num_fem_nodes_per_edge = 11
        self.name = "graph"
        self.num_eigenvalues = 7

    def add_edge(self, v0, v1):
        if v0 > self.num_vertices + 1 or v1 > self.num_vertices + 1:
            print "bad vertex index for added edge"
            return

        if (v0 == v1):
            print "identical vertex index in add_edge"
            return

        i0 = min(v0, v1)
        i1 = max(v0, v1)

        if ((i0, i1)) in self.edge_list:
            return

        self.edge_list.append((i0, i1))

        edge_index = self.edge_list.__len__() - 1

        self.num_vertices = max(self.num_vertices - 1, v0, v1) + 1

    def solve_helmholtz(self):
        self.check_dolfin()

        mesh = Mesh(self.get_mesh_name())

        V = FunctionSpace(mesh, "CG", 2)

        u = TrialFunction(V)
        v = TestFunction(V)
        a = dot(grad(u), grad(v))*dx
        b = v * dx

        # mass matrix
        m = u * v * dx
        M = PETScMatrix()
        assemble(m, tensor=M)

        # stiffness matrix
        A = PETScMatrix()
        assemble(a, tensor=A)

        eigensolver = SLEPcEigenSolver(A, M)
        eigensolver.parameters["spectrum"] = "smallest magnitude"
        eigensolver.parameters["solver"] = "lapack"
        eigensolver.parameters["tolerance"] = 1e-10

        print "starting eigensolver"
        eigensolver.solve(self.num_eigenvalues)
        print "finished eigensolver"

        for i in range(0, self.num_eigenvalues):
            lambda_r, lambda_c, x_r, x_c = eigensolver.get_eigenpair(i)
            print i, "&", sqrt(abs(lambda_r))

    def make_mesh(self):
        if self.num_fem_nodes_per_edge < 2:
            print "too few fem nodes per edge"
            return

        if self.num_fem_nodes_per_edge % 2 == 0:
            print "fem nodes per edge must be odd"
            return

        num_edges = self.edge_list.__len__()
        num_mesh_vertices = self.num_vertices + num_edges * (self.num_fem_nodes_per_edge - 2)
        num_mesh_cells = num_edges * (self.num_fem_nodes_per_edge - 1)

        mesh_file = open(self.get_mesh_name(), 'w')
        mesh_file.write('<?xml version="1.0"?>\n')
        mesh_file.write('<dolfin xmlns:dolfin="http://fenicsproject.org">\n')
        mesh_file.write('  <mesh celltype="interval" dim="1">\n')

        # vertices
        mesh_file.write('    <vertices size="' + str(num_mesh_vertices) + '">\n')
        count = 0

        for i in range(0, self.num_vertices):
            self.write_mesh_vertex(count, 0.0, mesh_file)
            count += 1

        edge_length = pi
        dx = edge_length / (self.num_fem_nodes_per_edge - 1)
        for e in range(0, num_edges):
            for i in range(1, self.num_fem_nodes_per_edge - 1):
                if i < (self.num_fem_nodes_per_edge + 1) / 2:
                    x = i * dx
                else:
                    x = edge_length - i * dx
                self.write_mesh_vertex(count, x, mesh_file)
                count += 1

        mesh_file.write('    </vertices>\n')

        # cells
        mesh_file.write('    <cells size="'+str(num_mesh_cells)+'">\n')


        v_start = self.num_vertices
        cell_count = 0
        for e in self.edge_list:
            vertex_list = [e[0]]
            for i in range(0, self.num_fem_nodes_per_edge - 2):
                vertex_list.append(v_start + i)

            vertex_list.append(e[1])
            self.write_mesh_cells(cell_count, vertex_list, mesh_file)
            cell_count += self.num_fem_nodes_per_edge - 1
            v_start += self.num_fem_nodes_per_edge - 2

        mesh_file.write('    </cells>\n')
        mesh_file.write('  </mesh>\n')
        mesh_file.write('</dolfin>')

        mesh_file.close()

    def write_mesh_vertex(self, count, x, mesh_file):
        mesh_file.write('      <vertex index="' + str(count) + '" x="' + str(x) + '"/>\n')
        return

    def write_mesh_cells(self, cell_count, vertex_list, mesh_file):
        for i in range(0, vertex_list.__len__() - 1):
            mesh_file.write('      <interval index="'+ str(cell_count) + '" v0="'+str(vertex_list[i])+'" v1="'+str(vertex_list[i+1])+'"/>\n')
            cell_count += 1

        return

    def print_eigenvalues(self):
        for i in range(0, self.num_eigenvalues):
            lambda_r, lambda_c, x_r, x_c = self.eigensolver.get_eigenpair(i)
            print i, "&", sqrt(abs(lambda_r)) / pi

    def check_dolfin(self):
        if not has_linear_algebra_backend("PETSc"):
            print "Dolfin has not been configured with PETSc"
            exit()

        if not has_slepc():
	        print "Dolfin has not been configured with SLEPc"
	        exit()

    def get_mesh_name(self):
        return self.name + "_mesh.xml"

    def __print__(self):
        print "Number of vertices ", self.num_vertices
        print "Number of edges ", self.edge_list.__len__()

    def print_edge(self, edge_index):
        print self.edge_list[edge_index]

    def print_adjaceny_list(self, vertex_index):
        print "adjacency list of vertex ", vertex_index, ": "
        for i in self.adjacency_list[vertex_index]:
            print i
