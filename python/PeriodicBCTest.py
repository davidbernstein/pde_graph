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
#set_log_level(DEBUG)

numEdges = 3
numVerticesPerEdge = 10

print "making mesh"
MakeMesh(numEdges, numVerticesPerEdge)
mesh = Mesh("mesh.xml")
print "finished mesh"

class PeriodicBoundary(SubDomain):

    # Left boundary is "target domain" G
    def inside(self, x, on_boundary):
    	return bool(near(x[0], 1) and on_boundary)

    # Map right boundary (H) to left boundary (G)
    def map(self, x, y):
    	y[0] = -1
    	if (x[0] >= 2  and x[0] <= 3):
        	y[0] = 1.0 - x[0] 
        if (x[0] >=4 and x[0] <= 5):
        	y[0] = 3.0 - x[0]

V = FunctionSpace(mesh, "CG", 2, constrained_domain=PeriodicBoundary())
#V = FunctionSpace(mesh, "CG", 1)
u = TrialFunction(V)
v = TestFunction(V)
a = dot(grad(u), grad(v))*dx

# mass matrix
m = u * v * dx
M = assemble(m)
M = as_backend_type(M)

b = v*dx
A,_= assemble_system(a, b)
A = as_backend_type(A)

eigensolver = SLEPcEigenSolver(A, M)
eigensolver.parameters["spectrum"] = "smallest real"
eigensolver.parameters["problem_type"] = "gen_hermitian"
eigensolver.parameters["solver"] = "krylov-schur"
#eigensolver.parameters["solver"] = "lanczos"
eigensolver.parameters["solver"] = "arnoldi"
eigensolver.parameters["tolerance"] = 1e-10
#eigensolver.parameters["verbose"] = True

numEigenvalues = 8

print "starting eigensolver"
eigensolver.solve(numEigenvalues)
print "finished eigensolver"

for i in range(0, numEigenvalues):
	lambda_r, lambda_c, x_r, x_c = eigensolver.get_eigenpair(i)
	print i, lambda_r, sqrt(abs(lambda_r))

lambda_r, lambda_c, x_r, x_c = eigensolver.get_eigenpair(4)
f = Function(V)
f.vector()[:] = x_r
plot(f, interactive=True)

