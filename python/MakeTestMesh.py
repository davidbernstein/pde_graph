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

from math import *


def make_test_mesh(num_vertices_per_edge):
    num_vertices = num_vertices_per_edge

    mesh_file = open('testMesh.xml', 'w')

    mesh_file.write('<?xml version="1.0"?>\n')
    mesh_file.write('<dolfin xmlns:dolfin="http://fenicsproject.org">\n')
    mesh_file.write('  <mesh celltype="interval" dim="1">\n')

    # vertices
    dx = 1.0 / (num_vertices_per_edge - 1)
    mesh_file.write('    <vertices size="' + str(num_vertices) + '">\n')
    for i in range(0, num_vertices_per_edge):
        if i < (num_vertices_per_edge + 1) / 2:
            x = i * dx
        else:
            x = 1 - i * dx
        mesh_file.write('      <vertex index="' + str(i) + '" x="' + str(x) + '"/>\n')

    mesh_file.write('    </vertices>\n')

    # cells
    num_cells = num_vertices_per_edge - 1
    mesh_file.write('    <cells size="'+str(num_cells)+'">\n')
    for i in range(0, num_cells):
        mesh_file.write('      <interval index="'+ str(i) + '" v0="'+str(i)+'" v1="'+str(i + 1)+'"/>\n')


    mesh_file.write('    </cells>\n')
    mesh_file.write('  </mesh>\n')
    mesh_file.write('</dolfin>')

    mesh_file.close()

    return
