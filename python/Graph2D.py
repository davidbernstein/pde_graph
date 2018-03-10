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

import random

class Graph2D:
    def __init__(self):
        self.edge_list = []
        self.vertex_list = []
        self.element_size = 0.05
        self.name = "graph"
        self.mesh = Mesh()
        self.f = []

    def erase(self):
        self.edge_list = []
        self.vertex_list = []
        self.mesh = Mesh()
        return;

    def solve_something(self):
        # function space
        V = FunctionSpace(self.mesh, "CG", 2)

        # edge length scaling coefficient
        g = MeshFunctionDouble(self.mesh,1)
        for c in range(0, self.mesh.num_cells()):
            g.set_value(c, self.f[c])

        DG0 = FunctionSpace(self.mesh, "DG", 0)

        class MyExpr1D(Expression):
            def __init__(self, cell_fun):
                assert(cell_fun.dim()==1)
                self.cell_fun = cell_fun
            def eval_cell(self, values, x, cell):
                values[0] = self.cell_fun[cell.index]

        E = MyExpr1D(g)
        f = Function(DG0)
        f.interpolate(E)

        x0 = self.vertex_list[0][0]
        y0 = self.vertex_list[0][1]
        def boundary(x, on_boundary):
            d = sqrt(pow((x[0]-x0), 2) + pow((x[1]-y0), 2))
            return d < 0.00001

        bc = DirichletBC(V, 1.0, boundary)

        u = Function(V)
        v = TestFunction(V)

        F = inner(grad(u), grad(v))*dx + f*(u**2)*v*dx

        #parameters['newton_solver']['relative_tolerance']=1.0e-6
        u.interpolate(Constant(1.0))
        solve(F == 0, u, bc)

        u.rename('u', 'u')
        file = File(self.name + '_u.pvd')
        file << u

        return

    def generate_random_unit_square(self, num_vertices):
        self.erase()

        random.seed(1)

        #self.add_vertex(0.5, 0.5)
        self.add_vertex(0.0, 1.0)

        for i in range(0, num_vertices-1):
            self.add_vertex(random.random(), random.random())

        #max_degree = num_vertices - 1
        max_degree = num_vertices / 3
        print 'max_degree = ' + str(max_degree)

        for i in range(0, num_vertices):
            out_degree = random.randint(2, max_degree)
            v0 = i
            for j in range(0, out_degree):
                v1 = random.randint(0, num_vertices-1)
                self.add_edge(v0, v1)

        return

    def make_triangle(self):
        self.erase()

        r = 1.0
        start_angle = pi
        angle = 2.0 * pi / 3.0
        
        for i in range(0, 3):
            self.add_vertex(r*cos(start_angle+i*angle), r*sin(start_angle+i*angle))

        self.add_edge(0, 1)
        self.add_edge(1, 2)
        self.add_edge(2, 0)

        return

    def generate_random_unit_circle(self, num_vertices):
        self.erase()
        random.seed(1)

        start_angle = pi
        angle = 2.0*pi/(num_vertices-1)

        r = 1.0
        for i in range(0, num_vertices):
            self.add_vertex(r*cos(start_angle+i*angle), r*sin(start_angle+i*angle))

        max_degree = 3
        print 'max_degree = ' + str(max_degree)

        for i in range(0, num_vertices):
            out_degree = random.randint(2, max_degree)
            degree = 0
            while degree < out_degree:
                v1 = random.randint(0, num_vertices-1)
                if self.add_edge(i, v1):
                    degree+= 1
                    #print str(i) + str(", ") + str(v1)

        return

    def add_vertex(self, x, y):
        if ((x, y)) in self.vertex_list:
            print 'identical vertex already in list'
            return

        self.vertex_list.append((x, y))
        return

    def add_edge(self, v0, v1):
        num_vertices = self.vertex_list.__len__()

        if v0 < 0 or v0 > num_vertices-1:
            print "bad vertex index for added edge"
            return False

        if v1 < 0 or v1 > num_vertices-1:
            print "bad vertex index for added edge"
            return False

        if (v0 == v1):
            return False

        i0 = min(v0, v1)
        i1 = max(v0, v1)

        if ((i0, i1)) in self.edge_list:
            return False

        self.edge_list.append((i0, i1))
        return True

    def make_mesh(self):
        editor = MeshEditor()
        editor.open(self.mesh, 1, 2)

        # edge discretizations
        num_cells = 0
        num_vertices_internal = 0
        n = []
        for i in range(0, self.edge_list.__len__()):
            d = self.get_edge_length(i)

            m = max(int(ceil(d / self.element_size)), 3)
            n.append(m)

            num_cells += m
            num_vertices_internal += m-1

        editor.init_cells(num_cells)
        editor.init_vertices(num_vertices_internal + self.vertex_list.__len__())

        # add vertices which correspond to graph vertices
        v_i = 0
        for v in self.vertex_list:
            editor.add_vertex(v_i, v[0], v[1])
            v_i += 1

        # parallel array to cells
        self.f = []

        cell_index = 0
        # add all other vertices and cells
        for i in range(0, self.edge_list.__len__()):
            v0 = self.edge_list[i][0]
            v1 = self.edge_list[i][1]

            x0 = self.vertex_list[v0][0]
            y0 = self.vertex_list[v0][1]

            x1 = self.vertex_list[v1][0]
            y1 = self.vertex_list[v1][1]

            l2 = (x0-x1)*(x0-x1)+(y0-y1)*(y0-y1)

            # internal vertices
            index = [v0]
            d = self.get_edge_length(i)
            ds = 1.0 / n[i]
            for j in range(0, n[i] - 1):
                s = (j + 1.0) * ds
                x = x0 * (1 - s) + x1 * s
                y = y0 * (1 - s) + y1 * s
                editor.add_vertex(v_i, x, y)
                index.append(v_i)
                v_i += 1

            index.append(v1)

            # add cells
            for j in range(0, n[i]):
                editor.add_cell(cell_index, index[j], index[j+1])
                self.f.append(1/l2)
                cell_index += 1

        editor.close()

        return

    def write_mesh(self):
        file_name = self.name + '_mesh.pvd'
        file = File(file_name)
        file << self.mesh
        return

    def get_edge_length(self, edge_index):
        v0 = self.edge_list[edge_index][0]
        v1 = self.edge_list[edge_index][1]

        x0 = self.vertex_list[v0][0]
        y0 = self.vertex_list[v0][1]

        x1 = self.vertex_list[v1][0]
        y1 = self.vertex_list[v1][1]

        return sqrt(pow((x0-x1),2) + pow((y0-y1),2))

    def check_dolfin(self):
        if not has_linear_algebra_backend("PETSc"):
            print "Dolfin has not been configured with PETSc"
            exit()

    def __print__(self):
        print "Number of vertices ", self.vertex_list.__len__()
        print "Number of edges ", self.edge_list.__len__()

    def print_edge(self, edge_index):
        print self.edge_list[edge_index]

    def print_adjaceny_list(self, vertex_index):
        print "adjacency list of vertex ", vertex_index, ": "
        for i in self.adjacency_list[vertex_index]:
            print i
