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

from init_instantiation_data import *

include_files = ['../../source/base/function_element.cpp',
                 '../../source/geometry/mapping_element.cpp',
                 '../../source/geometry/grid_forward_iterator.cpp']
data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)

for row in inst.all_mapping_dims:
    dims = '<%d, %d>' %(row.dim, row.codim)
    s = 'template class NewMapping%s ;\n' %(dims)
    f.write(s)
    s = 'template class MappingElement%s ;\n' %(dims)
    f.write(s)
    s = 'template class GridForwardIterator<MappingElement%s> ;\n' %(dims)
    f.write(s)

k_members = ['void NewMapping<dim,cod>::fill_cache<k>(ElementAccessor &elem, const int j);',
             'void NewMapping<dim,cod>::init_cache<k>(ElementAccessor &elem);',
             'void NewMapping<dim,cod>::reset<k>(const NewValueFlags flag, const Quadrature<k> &quad);']

for row in inst.all_mapping_dims:
    for fun in k_members:
        for k in range(max(0,row.dim-inst.n_sub_element), row.dim+1):
            s = fun.replace('dim','%d' %row.dim).replace('k','%d' %(k)).replace('cod','%d' %row.codim);
            f.write('template ' + s + '\n')
 
s = 'template class NewMapping<0,0> ;\n' 
f.write(s)      
s = 'template class MappingElement<0,0> ;\n'
f.write(s)
s = 'template class GridForwardIterator<MappingElement<0,0>> ;\n' 
f.write(s) 

