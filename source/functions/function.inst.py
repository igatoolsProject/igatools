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

include_files = []

data = Instantiation(include_files)

(f, inst) = (data.file_output, data.inst)

functions = []

for row in inst.all_function_dims:
    dims = '<%d, %d, %d, %d>' %(row.dim, row.codim, row.range, row.rank)
    func = 'Function%s' %(dims) 
    functions.append(func)
    f.write('template class %s ;\n' %(func))



#---------------------------------------------------
#f.write('IGA_NAMESPACE_CLOSE\n')
#  
#f.write('#ifdef SERIALIZATION\n')
#id = 0 
#for func in unique(functions):
#    alias = 'FunctionAlias%d' %(id)
#    f.write('using %s = iga::%s; \n' % (alias, func))
#    
#    f.write('ALLOW_SHARED_THIS(%s)\n' %alias )
#
#    f.write('BOOST_CLASS_EXPORT_IMPLEMENT(%s) \n' %alias)
#    f.write('template void %s::serialize(OArchive &, const unsigned int);\n' % alias)
#    f.write('template void %s::serialize(IArchive &, const unsigned int);\n' % alias)
#    id += 1 
#f.write('#endif // SERIALIZATION\n')
#      
#f.write('IGA_NAMESPACE_OPEN\n')
#---------------------------------------------------
