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
data = Instantiation()
(f, inst) = (data.file_output, data.inst)

class_list = []

class_name = 'TensorIndex'

for dim in inst.all_domain_dims:
    c = '%s<%d>' %(class_name,dim)
    class_list.append(c)




#---------------------------------------------------
f.write('IGA_NAMESPACE_CLOSE\n')

f.write('#ifdef SERIALIZATION\n')

id = 0 
for cls in unique(class_list):
    alias = '%sAlias%d' %(class_name,id)
    f.write('using %s = iga::%s; \n' % (alias, cls))
    f.write('CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(%s,cereal::specialization::member_serialize);\n' %(alias))
    id += 1 
f.write('#endif // SERIALIZATION\n')

f.write('IGA_NAMESPACE_OPEN\n')
#---------------------------------------------------
