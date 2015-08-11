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

include_files = ['../../source/geometry/cartesian_grid_iterator.cpp']
data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)


for x in inst.sub_mapping_dims:
    dims = '<Transformation::h_grad, %d, %d>' %(x.dim, x.codim)
    acc = 'PushForwardElement%s' % (dims)
    f.write('template class %s ;\n' %(acc))
#    for it in inst.iterators:
#        iterator = it.replace('Accessor','%s' % (acc) )
#        f.write('template class %s; \n' %iterator)

for x in inst.mapping_dims:
    dims = '<Transformation::h_grad, %d, %d>' %(x.dim, x.codim)
    acc = 'PushForwardElement%s' % (dims)
    f.write('template class %s ;\n' %(acc))
#    for it in inst.iterators:
#        iterator = it.replace('Accessor','%s' % (acc) )
#        f.write('template class %s; \n' %iterator)
