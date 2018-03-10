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

#from MakeEmbeddedPolygonMesh import *
#from MakeDiamondMesh import *
from MakeSquareGridMesh import *
#from MakeKinkedMesh import *
#from MakeTestMesh import *

if not has_linear_algebra_backend("PETSc"):
	print "Dolfin has not been configured with PETSc"
	exit()

if not has_slepc():
	print "Dolfin has not been configured with SLEPc"
	exit()

set_log_level(PROGRESS)
#set_log_level(DEBUG)

numEdges = 4
numVerticesPerEdge = 21

#MakeEmbeddedPolygonMesh(numEdges, numVerticesPerEdge)
#mesh = Mesh("polygonMesh.xml")

# MakeDiamondMesh(10)
# mesh = Mesh("diamondMesh.xml")

# MakeSquareGridMesh(numVerticesPerEdge)
# mesh = Mesh("squareGridMesh.xml")

# MakeKinkedMesh(numVerticesPerEdge)
# mesh = Mesh("kinkedMesh.xml")

MakeSquareGridMesh(numVerticesPerEdge)
mesh = Mesh("squareGridMesh.xml")

#plot(mesh, interactive=True)
meshFile = File("mesh.pvd")
meshFile << mesh

V = FunctionSpace(mesh, "CG", 2)

u = TrialFunction(V)
v = TestFunction(V)
a = dot(grad(u), grad(v))*dx
b = v * dx

# mass matrix
m = u * v * dx
M = PETScMatrix()
assemble(m, tensor=M)

# # stiffness matrix
A = PETScMatrix()
assemble(a, tensor=A)

eigensolver = SLEPcEigenSolver(A, M)
eigensolver.parameters["spectrum"] = "smallest magnitude"
eigensolver.parameters["solver"] = "lapack"
eigensolver.parameters["tolerance"] = 1e-10

numEigenvalues = 60

print "starting eigensolver"
eigensolver.solve(numEigenvalues)
print "finished eigensolver"

eigenFunction = Function(V)
for i in range(0, numEigenvalues):
	lambda_r, lambda_c, x_r, x_c = eigensolver.get_eigenpair(i)
	print i, "&", sqrt(abs(lambda_r)) / pi

	fileName = "EigenFunction_" + str(i) + ".pvd"
	file = File(fileName)

	eigenFunction.vector()[:] = x_r
	file << eigenFunction


