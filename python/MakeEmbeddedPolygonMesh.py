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

def MakeEmbeddedPolygonMesh(numEdges, numVerticesPerEdge):
	numVertices = numEdges * (numVerticesPerEdge - 1)

	numCellsPerEdge = numVerticesPerEdge - 1
	numCells = numEdges * numCellsPerEdge

	polyAngle = 2.0 * pi / numEdges

	edgeLength = 1.0
	radius = sqrt(0.5 / (1.0 - cos(polyAngle)))

	file = open('polygonMesh.xml', 'w')

	file.write('<?xml version="1.0"?>\n')
	file.write('<dolfin xmlns:dolfin="http://fenicsproject.org">\n')

	file.write('  <mesh celltype="interval" dim="2">\n')

	# vertices
	file.write('    <vertices size="' + str(numVertices) + '">\n')

	angle = 0
	count = 0
	for e in range(0, numEdges):
		x0 = radius * cos(angle)
		y0 = radius * sin(angle)

		angle = angle + polyAngle
		x1 = radius * cos(angle)
		y1 = radius * sin(angle)

		# check for unit distance
		#print "mesh distance ", sqrt(pow(x1-x0, 2) + pow(y0-y1, 2))

		for i in range(0, numVerticesPerEdge - 1):
			s = float(i) / (numVerticesPerEdge - 1)
			x = x0 * (1.0 - s) + x1 * s
			y = y0 * (1.0 - s) + y1 * s
			print s, x, y
			file.write('      <vertex index="' + str(count) + '" x="' + str(x) + '" y="' + str(y) + '"/>\n')
			count = count + 1

	file.write('    </vertices>\n')

	# cells
	file.write('    <cells size="'+str(numCells)+'">\n')
	count = 0
	for e in range(0, numCells - 1):
		file.write('      <interval index="'+ str(count) + '" v0="'+str(e)+'" v1="'+str(e+1)+'"/>\n')
		count = count + 1

	file.write('      <interval index="'+ str(count) + '" v0="'+str(numVertices-1)+'" v1="'+str(0)+'"/>\n')

	file.write('    </cells>\n')
	file.write('  </mesh>\n')
	file.write('</dolfin>')

	file.close()
