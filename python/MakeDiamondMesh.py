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

def MakeDiamondMesh(numVerticesPerEdge):
	global vertexCount, nVpE

	nVpE = numVerticesPerEdge

	numEdges = 5
	numVertices = numEdges * (numVerticesPerEdge - 1) - 1

	numCellsPerEdge = numVerticesPerEdge - 1
	numCells = numEdges * numCellsPerEdge

	angle = 2.0 * pi / 6

	file = open('diamondMesh.xml', 'w')

	file.write('<?xml version="1.0"?>\n')
	file.write('<dolfin xmlns:dolfin="http://fenicsproject.org">\n')

	file.write('  <mesh celltype="interval" dim="2">\n')

	# vertices
	vertexCount = 0
	file.write('    <vertices size="' + str(numVertices) + '">\n')

	# edge 0
	x0 = 0.5
	y0 = 0
	x1 = 0
	y1 = sin(angle)
	WriteVertices(x0, y0, x1, y1, 0, numVerticesPerEdge - 1, file)

	# edge 1
	x0 = x1
	y0 = y1
	x1 = -0.5
	y1 = 0
	WriteVertices(x0, y0, x1, y1, 0, numVerticesPerEdge - 1, file)

	# edge 2
	x0 = x1
	y0 = y1
	x1 = 0
	y1 = -sin(angle)
	WriteVertices(x0, y0, x1, y1, 0, numVerticesPerEdge - 1, file)

	# edge 3
	x0 = x1
	y0 = y1
	x1 = 0.5
	y1 = 0
	WriteVertices(x0, y0, x1, y1, 0, numVerticesPerEdge - 1, file)

	# edge 4
	x0 = 0.5
	y0 = 0
	x1 = -0.5
	y1 = 0
	WriteVertices(x0, y0, x1, y1, 1, numVerticesPerEdge - 1, file)


	file.write('    </vertices>\n')

	# cells
	file.write('    <cells size="'+str(numCells)+'">\n')
	count = 0
	for e in range(0, 4 * numCellsPerEdge - 1):
		file.write('      <interval index="'+ str(count) + '" v0="'+str(e)+'" v1="'+str(e+1)+'"/>\n')
		count = count + 1

	file.write('      <interval index="'+ str(count) + '" v0="'+str(count)+'" v1="'+str(0)+'"/>\n')
	count = count + 1

	nV = numVertices - numVerticesPerEdge + 2
	file.write('      <interval index="'+ str(count) + '" v0="'+str(0)+'" v1="'+str(nV)+'"/>\n')
	count = count + 1
	
	for i in range(nV, numVertices - 1):
		file.write('      <interval index="'+ str(count) + '" v0="'+str(i)+'" v1="'+str(i + 1)+'"/>\n')
		count = count + 1

	file.write('      <interval index="'+ str(count) + '" v0="'+str(numVertices-1)+'" v1="'+str(2*numVerticesPerEdge-2)+'"/>\n')


	file.write('    </cells>\n')
	file.write('  </mesh>\n')
	file.write('</dolfin>')

	file.close()

	return

def WriteVertices(x0, y0, x1, y1, i0, i1, file):
	global vertexCount, nVpE
	print sqrt(pow(x0-x1, 2) + pow(y0-y1, 2))
	for i in range(i0, i1):
		[x, y] = GetVertex(x0, y0, x1, y1, i)
		file.write('      <vertex index="' + str(vertexCount) + '" x="' + str(x) + '" y="' + str(y) + '"/>\n')
		vertexCount = vertexCount + 1

	return

def GetVertex(x0, y0, x1, y1, i):
	global nVpE
	s = float(i) / (nVpE - 1)
	x = x0 * (1.0 - s) + x1 * s
	y = y0 * (1.0 - s) + y1 * s
	return [x, y]
