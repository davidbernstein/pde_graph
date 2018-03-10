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


def MakeKinkedMesh(numVerticesPerEdge):
    global vertexCount, nVpE

    nVpE = numVerticesPerEdge

    if numVerticesPerEdge % 2 == 0:
        print "number of vertices must be odd"

    numEdges = 1
    numVertices = numEdges * numVerticesPerEdge

    numCellsPerEdge = numVerticesPerEdge - 1
    numCells = numEdges * numCellsPerEdge

    file = open('kinkedMesh.xml', 'w')

    file.write('<?xml version="1.0"?>\n')
    file.write('<dolfin xmlns:dolfin="http://fenicsproject.org">\n')
    file.write('  <mesh celltype="interval" dim="2">\n')

    angle = 2.0 * pi / 4
    vertexCount = 0

    file.write('    <vertices size="' + str(numVertices) + '">\n')

    n_half = (numVerticesPerEdge - 1) / 2 + 1
    x0 = 0.0
    y0 = 0.0
    x1 = 0.5 * cos(angle)
    y1 = 0.5 * sin(angle)
    WriteVertices(x0, y0, x1, y1, 0, n_half, file)

    x0 = x1
    y0 = y1
    x1 = 2 * x0
    y1 = 0
    WriteVertices(x0, y0, x1, y1, 1, n_half, file)

    file.write('    </vertices>\n')

    # cells
    file.write('    <cells size="' + str(numCells) + '">\n')
    for i in range(0, numCells):
        file.write('      <interval index="' + str(i) + '" v0="' + str(i) + '" v1="' + str(i + 1) + '"/>\n')

    file.write('    </cells>\n')
    file.write('  </mesh>\n')
    file.write('</dolfin>')

    file.close()

    return


def WriteVertices(x0, y0, x1, y1, i0, i1, file):
    global vertexCount, nVpE
    print sqrt(pow(x0 - x1, 2) + pow(y0 - y1, 2))
    for i in range(i0, i1):
        [x, y] = GetVertex(x0, y0, x1, y1, i)
        file.write('      <vertex index="' + str(vertexCount) + '" x="' + str(x) + '" y="' + str(y) + '"/>\n')
        vertexCount = vertexCount + 1

    return


def GetVertex(x0, y0, x1, y1, i):
    global nVpE
    n = (nVpE - 1) / 2 + 1
    s = float(i) / (n - 1)
    x = x0 * (1.0 - s) + x1 * s
    y = y0 * (1.0 - s) + y1 * s
    return [x, y]
