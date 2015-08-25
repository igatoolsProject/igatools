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

include_files = ['../../source/geometry/grid_iterator.cpp']
data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)



sub_dim_members = \
['const ValueVector<Points<dim+cod>> & PhysicalDomainElement<dim,cod>::get_boundary_normals<k>(const int s_id) const;']

elements = []

for x in inst.sub_mapping_dims:
    elem = 'PhysicalDomainElement<%d,%d>' %(x.dim, x.codim)
    elements.append(elem)
    for fun in sub_dim_members:
        k = x.dim
        s = fun.replace('cod', '%d' % (x.codim)).replace('dim', '%d' % (x.dim)).replace('k', '%d' % (k));
        f.write('template ' + s + '\n')

for x in inst.mapping_dims:
    elem = 'PhysicalDomainElement<%d,%d>' %(x.dim, x.codim)
    elements.append(elem)
    for fun in sub_dim_members:
        for k in inst.sub_dims(x.dim):
            s = fun.replace('dim','%d' %x.dim).replace('k','%d' %(k)).replace('cod','%d' %x.codim);
            f.write('template ' + s + '\n')
 

for elem in unique(elements):
    f.write('template class %s ;\n' %(elem))
    for it in inst.iterators:
        iterator = it.replace('Accessor','%s' % (elem) )
        f.write('template class %s; \n' %iterator)
