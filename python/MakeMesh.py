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

def MakeMesh(numEdges, numVerticesPerEdge):
	numVertices = numEdges * numVerticesPerEdge

	numCellsPerEdge = numVerticesPerEdge - 1
	numCells = numEdges * numCellsPerEdge

	file = open('mesh.xml', 'w')

	file.write('<?xml version="1.0"?>\n')
	file.write('<dolfin xmlns:dolfin="http://fenicsproject.org">\n')

	file.write('  <mesh celltype="interval" dim="1">\n')

	# vertices
	file.write('    <vertices size="' + str(numVertices) + '">\n')
	count = 0
	for e in range(0, numEdges):
		for i in range(0, numVerticesPerEdge):
			x = 2 * e + i / (numVerticesPerEdge - 1.0)
			file.write('      <vertex index="' + str(count) + '" x="' + str(x) + '"/>\n')
			count = count + 1

	file.write('    </vertices>\n')

	# cells
	file.write('    <cells size="'+str(numCells)+'">\n')
	n = 0
	for e in range(0, numEdges):
		vStart = e * numVerticesPerEdge;
		for i in range(vStart, vStart + numCellsPerEdge):
			file.write('      <interval index="'+ str(n) + '" v0="'+str(i)+'" v1="'+str(i+1)+'"/>\n')
			n = n + 1

	file.write('    </cells>\n')
	file.write('  </mesh>\n')
	file.write('</dolfin>')

	file.close()
