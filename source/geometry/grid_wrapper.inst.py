#-+--------------------------------------------------------------------
# Igatools a general purpose Isogeometric analysis library.
# Copyright (C) 2012-2015  by the igatools authors (see authors.txt).
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

include_files = ['geometry/cartesian_grid.h']
data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)

grids = ['CartesianGrid<%d>' %(dim) for dim in inst.all_domain_dims]

grid_wrappers = []

for grid in grids:
    grid_wrapper = 'GridWrapper<%s>' % (grid)
    grid_wrappers.append(grid_wrapper)
    f.write('template class %s; \n' % (grid_wrapper))



#---------------------------------------------------
f.write('IGA_NAMESPACE_CLOSE\n')
 
f.write('#ifdef SERIALIZATION\n')
id = 0 
for grid in unique(grids):
    alias = 'GridWrapperAlias%d' %(id)
    f.write('using %s = iga::GridWrapper<iga::%s>; \n' % (alias, grid))
    f.write('BOOST_CLASS_EXPORT(%s) \n' %alias)
    id += 1 
f.write('#endif // SERIALIZATION\n')
     
f.write('IGA_NAMESPACE_OPEN\n')
#---------------------------------------------------
