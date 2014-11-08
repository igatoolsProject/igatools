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

include_files = ['geometry/cartesian_grid_element.h']
data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)

cartesian_grids = ['GridElementHandler<%d>' % (dim) for dim in inst.domain_dims]
for row in cartesian_grids:
   f.write('template class %s; \n' % (row))


k_members = ['void GridElementHandler<dim>::fill_cache<k>(ElementAccessor &elem, const int j);',
             'void GridElementHandler<dim>::init_cache<k>(ElementAccessor &elem);',
             'void GridElementHandler<dim>::reset<k>(const NewValueFlags flag, const Quadrature<k> &quad);',
             'Size GridElementHandler<dim>::get_num_points<k>() const;']
for dim in inst.domain_dims:
    for fun in k_members:
        for k in range(dim, max(0,dim-inst.n_sub_element) - 1, -1):
            s = fun.replace('dim','%d' %dim).replace('k','%d' %(k));
            f.write('template ' + s + '\n')
