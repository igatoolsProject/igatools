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

include_files = ['utils/array.h',
                 '../../source/base/function_element.cpp',
                 '../../source/geometry/mapping_element.cpp',
                 '../../source/geometry/cartesian_grid_iterator.cpp']
data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)



sub_dim_members = \
['void Mapping<dim,cod>::fill_cache<k>(ElementAccessor &elem, const int j);',
 'void Mapping<dim,cod>::init_cache<k>(ElementAccessor &elem);',
 'void Mapping<dim,cod>::reset<k>(const ValueFlags flag, const Quadrature<k> &quad);',
 'const ValueVector<Points<dim+cod>> & MappingElement<dim,cod>::get_boundary_normals<k>(const int s_id) const;']

for x in inst.sub_mapping_dims:
    dims = '<%d, %d>' %(x.dim, x.codim)
    f.write('template class Mapping%s ;\n' %(dims))
    acc = 'MappingElement%s' % (dims)
    f.write('template class %s ;\n' %(acc))
    for it in inst.iterators:
        iterator = it.replace('Accessor','%s' % (acc) )
        f.write('template class %s; \n' %iterator)
    for fun in sub_dim_members:
        k = x.dim
        s = fun.replace('cod', '%d' % (x.codim)).replace('dim', '%d' % (x.dim)).replace('k', '%d' % (k));
        f.write('template ' + s + '\n')

for x in inst.mapping_dims:
    dims = '<%d, %d>' %(x.dim, x.codim)
    f.write('template class Mapping%s ;\n' %(dims))
    acc = 'MappingElement%s' % (dims)
    f.write('template class %s ;\n' %(acc))
    for it in inst.iterators:
        iterator = it.replace('Accessor','%s' % (acc) )
        f.write('template class %s; \n' %iterator)
    for fun in sub_dim_members:
        for k in inst.sub_dims(x.dim):
            s = fun.replace('dim','%d' %x.dim).replace('k','%d' %(k)).replace('cod','%d' %x.codim);
            f.write('template ' + s + '\n')
 
