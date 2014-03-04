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

###############################################################################
# Common header for instantiation files 
from init_instantiation_data import *
file_output, inst = intialize_instantiation()
###############################################################################


include_files = ['#include <igatools/basis_functions/bspline_space.h>\n',
                 '#include <igatools/basis_functions/bspline_element_accessor.h>\n',
                 '#include <igatools/basis_functions/nurbs_space.h>\n',
                 '#include <igatools/basis_functions/nurbs_element_accessor.h>\n',
                 '#include <igatools/basis_functions/physical_space.h>\n',
                 '#include <igatools/geometry/cartesian_grid_element_accessor.h>\n'
                 '#include <igatools/geometry/mapping_element_accessor.h>\n',
                 '#include <igatools/geometry/push_forward_element_accessor.h>\n',
                 '#include <igatools/basis_functions/physical_space_element_accessor.h>\n']
for include in include_files:
    file_output.write(include)
file_output.write('IGA_NAMESPACE_OPEN\n')

strings = []
spaces = ['BSplineSpace', 'NURBSSpace']

for row in inst.user_table:
   writer = 'Writer<%d, %d, double>' %(row.dim, row.space_dim)
   strings.append('template class %s ;\n' % (writer))
   for name in spaces:
      space_ref  = '%s<%d,%d,%d>' % (name, row.dim, row.range, row.rank)
      PushForward = 'PushForward<Transformation::%s,%d,%d>' %(row.trans_type, row.dim, row.codim)
      space_phys = 'PhysicalSpace<%s,%s>' %(space_ref,PushForward)
      func = 'add_field<%s>(shared_ptr<%s>, const Vector &, const string & )' % (space_phys,space_phys)
      strings.append('template void %s::%s ;\n' % (writer,func))
      func = 'add_field<%s>(shared_ptr<%s>, const Vector &, const string & )' % (space_ref,space_ref)
      strings.append('template void %s::%s ;\n' % (writer,func))


for s in unique(strings): # Removing repeated entries.
    file_output.write(s)


file_output.write('IGA_NAMESPACE_CLOSE\n')
file_output.close()

	
