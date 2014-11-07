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

include_files = ['../../source/basis_functions/new_bspline_space.cpp',
                 '../../source/basis_functions/bspline_element.cpp',
                 '../../source/geometry/grid_forward_iterator.cpp'
                 ]
data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)


for x in inst.really_all_ref_sp_dims:
#     f.write('template class NewBSplineSpace<%d, %d, %d>; \n'
#             %(x.dim, x.range, x.rank))
    f.write('template class SpaceElement<NewBSplineSpace<%d, %d, %d>>; \n'
            %(x.dim, x.range, x.rank))
    f.write('template class BSplineElement<%d, %d, %d>; \n' 
            %(x.dim, x.range, x.rank))
    f.write('template class GridForwardIterator<BSplineElement<%d, %d, %d>>; \n' 
            %(x.dim, x.range, x.rank))
    elemhandler = 'BSplineElementHandler<%d, %d, %d>' %(x.dim, x.range, x.rank)
    f.write('template class BSplineElementHandler<%d, %d, %d>; \n' 
            %(x.dim, x.range, x.rank))
    k_members = ['void %s::fill_cache<k>(ElementAccessor &elem, const int j);' %elemhandler,
             'void %s::init_cache<k>(ElementAccessor &elem);' %elemhandler,
             'void %s::reset<k>(const NewValueFlags flag, const Quadrature<k> &quad);' %elemhandler]
    for fun in k_members:
        for k in range(x.dim, max(0, x.dim - inst.n_sub_element) - 1, -1):
            s = fun.replace('k', '%d' % (k));
            f.write('template ' + s + '\n')           
   
