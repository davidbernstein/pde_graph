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

from Graph2D import *

g = Graph2D()

#g.generate_random_unit_circle(6)
g.make_triangle()

g.make_mesh()
g.write_mesh()
g.solve_something()

g.__print__()