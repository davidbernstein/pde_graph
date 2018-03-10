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

def MakeSquareGridMesh(numVerticesPerEdge):
	global vertexCount, nVpE

	nVpE = numVerticesPerEdge

	numEdges = 12
	numVertices = numEdges * numVerticesPerEdge - 4 - 8 - 3

	numCellsPerEdge = numVerticesPerEdge - 1
	numCells = numEdges * numCellsPerEdge

	file = open('squareGridMesh.xml', 'w')

	file.write('<?xml version="1.0"?>\n')
	file.write('<dolfin xmlns:dolfin="http://fenicsproject.org">\n')
	file.write('  <mesh celltype="interval" dim="2">\n')

	# vertices
	vertexCount = 0
	file.write('    <vertices size="' + str(numVertices) + '">\n')

	# bottom row
	WriteVertices(0, 0, 1, 0, 0, numVerticesPerEdge - 1, file)
	WriteVertices(1, 0, 2, 0, 0, numVerticesPerEdge, file)

	# middle row
	WriteVertices(0, 1, 1, 1, 0, numVerticesPerEdge - 1, file)
	WriteVertices(1, 1, 2, 1, 0, numVerticesPerEdge, file)

	# top row
	WriteVertices(0, 2, 1, 2, 0, numVerticesPerEdge - 1, file)
	WriteVertices(1, 2, 2, 2, 0, numVerticesPerEdge, file)

	# left column
	WriteVertices(0, 0, 0, 1, 1, numVerticesPerEdge - 1, file)
	WriteVertices(0, 1, 0, 2, 1, numVerticesPerEdge - 1, file)

	# middle column
	WriteVertices(1, 0, 1, 1, 1, numVerticesPerEdge - 1, file)
	WriteVertices(1, 1, 1, 2, 1, numVerticesPerEdge - 1, file)

	# right column
	WriteVertices(2, 0, 2, 1, 1, numVerticesPerEdge - 1, file)
	WriteVertices(2, 1, 2, 2, 1, numVerticesPerEdge - 1, file)

	file.write('    </vertices>\n')

	# cells
	file.write('    <cells size="'+str(numCells)+'">\n')
	count = 0

	# bottom row
	v0 = 0
	count = 0
	for i in range(0, numCellsPerEdge):
		WriteCell(file, count, v0, v0+1)
		count = count + 1
		v0 = v0 + 1

	for i in range(0, numCellsPerEdge):
		WriteCell(file, count, v0, v0+1)
		count = count + 1
		v0 = v0 + 1

	# middle row
	v0 = v0 + 1
	for i in range(0, numCellsPerEdge):
		WriteCell(file, count, v0, v0+1)
		count = count + 1
		v0 = v0 + 1

	for i in range(0, numCellsPerEdge):
		WriteCell(file, count, v0, v0+1)
		count = count + 1
		v0 = v0 + 1

	# top row
	v0 = v0 + 1
	for i in range(0, numCellsPerEdge):
		WriteCell(file, count, v0, v0+1)
		count = count + 1
		v0 = v0 + 1

	for i in range(0, numCellsPerEdge):
		WriteCell(file, count, v0, v0+1)
		count = count + 1
		v0 = v0 + 1

	# left column
	v0 = v0 + 1
	for i in range(1, numCellsPerEdge-1):
		WriteCell(file, count, v0, v0+1)
		count = count + 1
		v0 = v0 + 1

	v0 = v0 + 1
	for i in range(1, numCellsPerEdge-1):
		WriteCell(file, count, v0, v0+1)
		count = count + 1
		v0 = v0 + 1

	# middle column
	v0 = v0 + 1
	for i in range(1, numCellsPerEdge-1):
		WriteCell(file, count, v0, v0+1)
		count = count + 1
		v0 = v0 + 1

	v0 = v0 + 1
	for i in range(1, numCellsPerEdge-1):
		WriteCell(file, count, v0, v0+1)
		count = count + 1
		v0 = v0 + 1

	# right column
	v0 = v0 + 1
	for i in range(1, numCellsPerEdge-1):
		WriteCell(file, count, v0, v0+1)
		count = count + 1
		v0 = v0 + 1

	v0 = v0 + 1
	for i in range(1, numCellsPerEdge-1):
		WriteCell(file, count, v0, v0+1)
		count = count + 1
		v0 = v0 + 1

	# left bits
	v0 = 3 * (2 * numVerticesPerEdge - 1) 
	WriteCell(file, count, 0, v0)
	count = count + 1

	v0 = v0 + numVerticesPerEdge - 3
	WriteCell(file, count, v0, 2*numVerticesPerEdge-1)
	count = count + 1

	WriteCell(file, count, v0 + 1, 2*numVerticesPerEdge-1)
	count = count + 1

	v0 = v0 + numVerticesPerEdge - 2
	WriteCell(file, count, v0, 4*numVerticesPerEdge-2)
	count = count + 1

	# middle bits
	v0 = v0 + 1
	WriteCell(file, count, v0, numVerticesPerEdge-1)
	count = count + 1

	v0 = v0 + numVerticesPerEdge - 3
	WriteCell(file, count, v0, 3*numVerticesPerEdge-2)
	count = count + 1

	WriteCell(file, count, v0 + 1, 3*numVerticesPerEdge-2)
	count = count + 1

	v0 = v0 + numVerticesPerEdge - 2
	WriteCell(file, count, v0, 5*numVerticesPerEdge-3)
	count = count + 1

	# right bits
	v0 = v0 + 1
	WriteCell(file, count, v0, 2*numVerticesPerEdge-2)
	count = count + 1

	v0 = v0 + numVerticesPerEdge - 3
	WriteCell(file, count, v0, 4*numVerticesPerEdge-3)
	count = count + 1

	WriteCell(file, count, v0 + 1, 4*numVerticesPerEdge-3)
	count = count + 1

	v0 = v0 + numVerticesPerEdge - 2
	WriteCell(file, count, v0, 6*numVerticesPerEdge-4)
	count = count + 1

	file.write('    </cells>\n')
	file.write('  </mesh>\n')
	file.write('</dolfin>')

	file.close()

	return

def WriteCell(file, count, v0, v1):
	file.write('      <interval index="'+ str(count) + '" v0="'+str(v0)+'" v1="'+str(v1)+'"/>\n')
	return

def WriteVertices(x0, y0, x1, y1, i0, i1, file):
	global vertexCount, nVpE
	#print sqrt(pow(x0-x1, 2) + pow(y0-y1, 2))
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
