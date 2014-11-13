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

# QA (pauletti, Mar 19, 2014):
from init_instantiation_data import *

include_files = ['basis_functions/new_bspline_space.h',
                 'basis_functions/nurbs_space.h',
                 'geometry/new_push_forward.h',
                 'geometry/cartesian_grid_element.h',
                 'geometry/mapping_element.h',
                 'geometry/push_forward_element.h',
                 'basis_functions/bspline_element.h',
                 '../../source/basis_functions/physical_space_element.cpp',
                 '../../source/basis_functions/bspline_element.cpp',
                 'basis_functions/nurbs_element_accessor.h',
                 'basis_functions/physical_space_element.h',
                 '../../source/geometry/grid_forward_iterator.cpp']
data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)

for space in inst.AllPhysSpaces:
      f.write( 'template class %s;\n' %space.name)
