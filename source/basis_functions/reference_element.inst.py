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

# QA (pauletti, Jun 27, 2014):
from init_instantiation_data import *

include_files = ['geometry/grid_iterator.h']

data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)


member_functions = set()


elements = ['ReferenceElement<0,0,1>']#, 'ReferenceElement<0,1,1>']

member_func = 'ValueVector<Real> ReferenceElement<0,0,1>::get_w_measures<0>(const int s_id) const'
member_functions.add(member_func)



for x in inst.sub_ref_sp_dims + inst.ref_sp_dims:
    elem = 'ReferenceElement<%d,%d,%d>' %(x.dim, x.range, x.rank)
    elements.append(elem)
    for k in range(0,x.dim+1):
        member_func = 'ValueVector<Real> %s::get_w_measures<%d>(const int s_id) const' % (elem,k)
        member_functions.add(member_func)


for elem in unique(elements):
    f.write('template class %s; \n' %elem)
    f.write('template class GridIterator<%s>; \n' %elem)


for func in member_functions:
    f.write('template %s;\n' %(func))



#---------------------------------------------------


