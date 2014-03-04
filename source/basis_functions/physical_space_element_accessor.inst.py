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

include_files = ['#include <igatools/basis_functions/bspline_space.h>\n',
                 '#include <igatools/basis_functions/bspline_element_accessor.h>\n',
                 '#include <igatools/basis_functions/nurbs_space.h>\n',
                 '#include <igatools/basis_functions/nurbs_element_accessor.h>\n',
                 '#include <igatools/basis_functions/physical_space.h>\n',
                 '#include <igatools/geometry/cartesian_grid_element_accessor.h>\n'
                 '#include <igatools/geometry/mapping_element_accessor.h>\n',
                 '#include <igatools/geometry/push_forward_element_accessor.h>\n',
                 '#include <igatools/basis_functions/physical_space_element_accessor.h>\n',
                 '#include <../source/geometry/grid_forward_iterator.cpp>\n']
for file in include_files:
    file_output.write(file)
file_output.write('IGA_NAMESPACE_OPEN\n')

for space in inst.PhysSpaces:
   accessor = ('PhysicalSpaceElementAccessor<%s>' %space)
   file_output.write('template class %s;\n' %accessor)
   file_output.write('template class GridForwardIterator<%s> ;\n' %accessor)

file_output.write('IGA_NAMESPACE_CLOSE\n')
file_output.close()
    
	

