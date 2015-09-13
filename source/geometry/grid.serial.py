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

from init_instantiation_data import *

data = Instantiation()
(f, inst) = (data.file_output, data.inst)
grids = [] 
for dim in inst.sub_domain_dims:
    grid = 'Grid<%d>' %(dim)
    grids.append(grid)
   
    
for dim in inst.domain_dims:
    grid = 'Grid<%d>' %(dim)   
    grids.append(grid)
         
#---------------------------------------------------
f.write('IGA_NAMESPACE_CLOSE\n')



f.write('#ifdef SERIALIZATION\n')
archives = ['OArchive','IArchive']

id = 0 
for grid in unique(grids):
    alias = 'GridAlias%d' %(id)
    f.write('using %s = iga::%s; \n' % (alias, grid))
#    f.write('CEREAL_REGISTER_TYPE(%s);\n' %(alias))
    for ar in archives:
        f.write('template void %s::serialize(%s&);\n' %(alias,ar))
#    
#    f.write('ALLOW_SHARED_THIS(%s)\n' %alias )
#    
#    f.write('BOOST_CLASS_EXPORT_IMPLEMENT(%s) \n' %alias)
#    f.write('template void %s::serialize(OArchive &, const unsigned int);\n' % alias)
#    f.write('template void %s::serialize(IArchive &, const unsigned int);\n' % alias)
    id += 1 
f.write('#endif // SERIALIZATION\n')
#   
f.write('IGA_NAMESPACE_OPEN\n')
#---------------------------------------------------