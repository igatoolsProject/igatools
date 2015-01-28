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

include_files = ['basis_functions/bspline_space.h',
                 '../../source/basis_functions/bspline_element.cpp',
                 '../../source/geometry/cartesian_grid_iterator.cpp'
                 ]
data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)


sub_dim_members = \
[]
#'void elhandler::fill_cache<k>(ElementAccessor &elem, const int j);']
#,
#'void elhandler::init_cache<k>(ElementAccessor &elem);']
# ,
#'void elhandler::reset<k>(const ValueFlags flag, const QuadratureTensorProduct<k> &quad);']


for x in inst.sub_ref_sp_dims:
    space = 'ReferenceSpace<%d, %d, %d>' %(x.dim, x.range, x.rank)
#    f.write('template class SpaceElement<%s>; \n' %space)
    acc = 'ReferenceElement<%d, %d, %d>' %(x.dim, x.range, x.rank)
#    f.write('template class %s; \n' %acc)
    for it in inst.iterators:
        iterator = it.replace('Accessor','%s' % (acc) )
#        f.write('template class %s; \n' %iterator)
    elemhandler = 'ReferenceElementHandler<%d, %d, %d>' %(x.dim, x.range, x.rank)
    f.write('template class %s; \n'  %elemhandler)
#    for fun in sub_dim_members:
#        k = x.dim
#        s = fun.replace('elhandler', elemhandler).replace('k', '%d' % (k));
#        f.write('template ' + s + '\n')


for x in inst.ref_sp_dims:
    space = 'ReferenceSpace<%d, %d, %d>' %(x.dim, x.range, x.rank)
#    f.write('template class SpaceElement<%s>;' %space)
    acc = 'ReferenceElement<%d, %d, %d>' %(x.dim, x.range, x.rank)
#    f.write('template class %s; \n' %acc)
    for it in inst.iterators:
        iterator = it.replace('Accessor','%s' % (acc) )
#        f.write('template class %s; \n' %iterator)
    elemhandler = 'ReferenceElementHandler<%d, %d, %d>' %(x.dim, x.range, x.rank)
    f.write('template class %s; \n'  %elemhandler)
#    for fun in sub_dim_members:
#        for k in inst.sub_dims(x.dim):
#            s = fun.replace('elhandler', elemhandler).replace('k', '%d' % (k));
#            f.write('template ' + s + '\n')


