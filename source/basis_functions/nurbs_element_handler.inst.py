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

include_files = ['basis_functions/nurbs_space.h',
                 '../../source/basis_functions/nurbs_element.cpp',
                 '../../source/geometry/cartesian_grid_iterator.cpp'
                 ]
data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)


sub_dim_members = []
#\
#['void elhandler::fill_cache<k>(ElementAccessor &elem, const int j);',
#'void elhandler::init_cache<k>(ElementAccessor &elem);']
# ,
#'void elhandler::reset<k>(const ValueFlags flag, const Quadrature<k> &quad);']

for x in inst.sub_ref_sp_dims:
#    space = 'NURBSSpace<%d, %d, %d>' %(x.dim, x.range, x.rank)
#    space_elem = 'SpaceElement<%s,%d,0,%d,%d>' %(space,x.dim, x.range, x.rank)
#    f.write('template class %s; \n' %space_elem)
    elem = 'NURBSElement<%d, %d, %d>' %(x.dim, x.range, x.rank)
    f.write('template class %s; \n' %elem)
#    acc = 'ReferenceElement<%d, %d, %d>' %(x.dim, x.range, x.rank)
#    for it in inst.iterators:
#        iterator = it.replace('Accessor','%s' % (acc) )
#        f.write('template class %s; \n' %iterator)
#    f.write('template class CartesianGridIterator<%s>; \n' %acc)
    elemhandler = 'NURBSElementHandler<%d, %d, %d>' %(x.dim, x.range, x.rank)
    f.write('template class %s; \n'  %elemhandler)
    for fun in sub_dim_members:
        k = x.dim
        s = fun.replace('elhandler', elemhandler).replace('k', '%d' % (k));
        f.write('template ' + s + '\n')



for x in inst.ref_sp_dims:
#    space = 'NURBSSpace<%d, %d, %d>' %(x.dim, x.range, x.rank)
#    space_elem = 'SpaceElement<%s,%d,0,%d,%d>' %(space,x.dim, x.range, x.rank)
#    f.write('template class %s; \n' %space_elem)
    elem = 'NURBSElement<%d, %d, %d>' %(x.dim, x.range, x.rank)
    f.write('template class %s; \n' %elem)
#    acc = 'ReferenceElement<%d, %d, %d>' %(x.dim, x.range, x.rank)
#    for it in inst.iterators:
#        iterator = it.replace('Accessor','%s' % (acc) )
#        f.write('template class %s; \n' %iterator)
#    f.write('template class CartesianGridIterator<%s>; \n' %acc)
    elemhandler = 'NURBSElementHandler<%d, %d, %d>' %(x.dim, x.range, x.rank)
    f.write('template class %s; \n'  %elemhandler)
    for fun in sub_dim_members:
        for k in inst.sub_dims(x.dim):
            s = fun.replace('elhandler', elemhandler).replace('k', '%d' % (k));
            f.write('template ' + s + '\n')
        
   
