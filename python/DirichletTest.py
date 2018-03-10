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
from MakeMesh import *

if not has_linear_algebra_backend("PETSc"):
	print "Dolfin has not been configured with PETSc"
	exit()

if not has_slepc():
	print "Dolfin has not been configured with SLEPc"
	exit()

#set_log_level(PROGRESS)
set_log_level(DEBUG)

numEdges = 1
numVerticesPerEdge = 10

print "making mesh"
MakeMesh(numEdges, numVerticesPerEdge)
mesh = Mesh("mesh.xml")
print "finished mesh"

V = FunctionSpace(mesh, "CG", 1)

u = TrialFunction(V)
v = TestFunction(V)
a = dot(grad(u), grad(v))*dx
b = v * dx

def boundary0(x, on_boundary):
    return x[0] < 0.01

def boundary1(x, on_boundary):
    return abs(1-x[0]) < 0.01
    
bc0 = DirichletBC(V, 0, boundary0)
bc1 = DirichletBC(V, 0, boundary1)

bcs = [bc0, bc1];

# mass matrix
m = u * v * dx
M = PETScMatrix()
assemble(m, tensor=M)
for i in range(0, 2):
    bcs[i].apply(M)

# stiffness matrix
A = PETScMatrix()
assemble_system(a, b, bcs, A_tensor = A)

eigensolver = SLEPcEigenSolver(A, M)
eigensolver.parameters["spectrum"] = "smallest magnitude"
eigensolver.parameters["solver"] = "lapack"
eigensolver.parameters["tolerance"] = 1e-10

numEigenvalues = 9

print "starting eigensolver"
eigensolver.solve(numEigenvalues)
print "finished eigensolver"

for i in range(0, numEigenvalues):
	lambda_r, lambda_c, x_r, x_c = eigensolver.get_eigenpair(i)
	print i, lambda_r, lambda_c, sqrt(abs(lambda_r))

# lambda_r, lambda_c, x_r, x_c = eigensolver.get_eigenpair(0)
# f = Function(V)
# f.vector()[:] = x_r
# plot(f, interactive=True)

