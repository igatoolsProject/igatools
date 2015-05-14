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
#'void elhandler::reset<k>(const ValueFlags flag, const Quadrature<k> &quad);']


handlers = []

for x in inst.sub_ref_sp_dims:
    space = 'ReferenceSpace<%d, %d, %d>' %(x.dim, x.range, x.rank)
    acc = 'ReferenceElement<%d, %d, %d>' %(x.dim, x.range, x.rank)
    for it in inst.iterators:
        iterator = it.replace('Accessor','%s' % (acc) )
    elemhandler = 'ReferenceElementHandler<%d, %d, %d>' %(x.dim, x.range, x.rank)
    handlers.append(elemhandler)
    f.write('template class %s; \n'  %elemhandler)


for x in inst.ref_sp_dims:
    space = 'ReferenceSpace<%d, %d, %d>' %(x.dim, x.range, x.rank)
    acc = 'ReferenceElement<%d, %d, %d>' %(x.dim, x.range, x.rank)
    for it in inst.iterators:
        iterator = it.replace('Accessor','%s' % (acc) )
    elemhandler = 'ReferenceElementHandler<%d, %d, %d>' %(x.dim, x.range, x.rank)
    handlers.append(elemhandler)
    f.write('template class %s; \n'  %elemhandler)
#    for fun in sub_dim_members:
#        for k in inst.sub_dims(x.dim):
#            s = fun.replace('elhandler', elemhandler).replace('k', '%d' % (k));
#            f.write('template ' + s + '\n')

#---------------------------------------------------
f.write('IGA_NAMESPACE_CLOSE\n')
 
f.write('#ifdef SERIALIZATION\n')
id = 0 
for handler in unique(handlers):
    alias = 'ReferenceElementHandlerAlias%d' %(id)
    f.write('using %s = iga::%s; \n' % (alias, handler))
    f.write('BOOST_CLASS_EXPORT_IMPLEMENT(%s) \n' %alias)
    f.write('template void %s::serialize(OArchive &, const unsigned int);\n' % alias)
    f.write('template void %s::serialize(IArchive &, const unsigned int);\n' % alias)
    id += 1 
f.write('#endif // SERIALIZATION\n')
     
f.write('IGA_NAMESPACE_OPEN\n')
#---------------------------------------------------

