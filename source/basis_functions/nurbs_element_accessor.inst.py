#-+--------------------------------------------------------------------
# Igatools a general purpose Isogeometric analysis library.
# Copyright (C) 2012-2014  by the igatools authors (see authors.txt).
#
# This file is part of the igatools library.
#
# The igatools library is free software: you can use it, redistribute
# it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#-+--------------------------------------------------------------------

# QA (pauletti, Mar 4, 2014 ):
from init_instantiation_data import *
file_output, inst = intialize_instantiation()

include_files = ['#include <igatools/geometry/cartesian_grid.h>\n',
                 '#include <igatools/geometry/cartesian_grid_element_accessor.h>\n'
                 '#include <igatools/basis_functions/bspline_space.h>\n',
                 '#include <../source/geometry/grid_forward_iterator.cpp>\n']
for include in include_files:
    file_output.write(include)
file_output.write('IGA_NAMESPACE_OPEN\n')

spaces = ['NURBSElementAccessor<%d, %d, %d>' %(x.dim, x.range, x.rank)  
          for x in inst.all_ref_sp_dims ]
for sp in spaces:
   file_output.write('template class %s ;\n' %sp)
   file_output.write('template class GridForwardIterator<%s> ;\n' %sp)

file_output.write('IGA_NAMESPACE_CLOSE\n')
file_output.close()




