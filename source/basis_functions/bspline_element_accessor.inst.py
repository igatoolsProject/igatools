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

# QA (pauletti, Jun 27, 2014):
from init_instantiation_data import *

include_files = ['geometry/cartesian_grid.h',
                 'geometry/cartesian_grid_element_accessor.h',
                 'basis_functions/bspline_space.h',
                 '../../source/geometry/grid_forward_iterator.cpp']
data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)

accessors = [('BSplineElementAccessor<%d, %d, %d>' %(x.dim, x.range, x.rank), x.dim)
             for x in inst.really_all_ref_sp_dims]

for acc in accessors:
    f.write('template class %s ;\n' %acc[0])
    f.write('template class GridForwardIterator<%s> ;\n' %acc[0])

#     fun = ('template void ' + acc[0] + '::evaluate_bspline_derivatives<order>' +
#            '(const ComponentContainer<std::array<const BasisValues1d *, %d' %acc[1] 
#            + '>> &,' +
#            'const ' + acc[0] + '::ValuesCache &,' +
#            'ValueTable<Conditional<order==0,Value,Derivative<order>>> &) const; \n')
#     fun_list = [fun.replace('order', str(d)) for d in inst.deriv_order]
#     for s in fun_list:
#         f.write(s)

    fun_2 = ('template ValueTable< Conditional< order==0,'+
               acc[0] + '::Value,' +
               acc[0] + '::Derivative<order> > > ' + 
               acc[0] + '::evaluate_basis_derivatives_at_points<order>' +
               '(const ValueVector<Point>&) const; \n')
    fun_2_list = [fun_2.replace('order', str(d)) for d in inst.deriv_order]
    for s in fun_2_list:
        f.write(s)

