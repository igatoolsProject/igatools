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

###############################################################################
# Common header for instantiation files 
from init_instantiation_data import *
file_output, inst = intialize_instantiation()
###############################################################################

include_files = ['#include <igatools/geometry/cartesian_grid.h>\n',
                 '#include <igatools/geometry/cartesian_grid_element_accessor.h>\n'
                 '#include <igatools/basis_functions/bspline_space.h>\n']

for include in include_files:
    file_output.write(include)
    
file_output.write('IGA_NAMESPACE_OPEN\n')


strings = []   
for row in inst.all_table:
    CA = 'BSplineElementAccessor< %d, %d, %d >' % (row.dim,
                row.range, row.rank)
    strings.append('template class %s ;\n' % (CA))

    master = ('template  void ' + CA + '::evaluate_bspline_derivatives<deriv_order>' +
              '(const FuncPointSize &,' +
              'StaticMultiArray<std::array<const BasisValues1d*, dim_domain>, dim_range, rank> &,'+
              'ValueTable< Derivative<deriv_order> > &) const; \n')
    
    master1 = master.replace('dim_domain', str(row.dim) ).replace('dim_range', str(row.range) ).replace('rank', str(row.rank) );
    for d in inst.deriv_order:
        strings.append(master1.replace('deriv_order', str(d)))
#         VT = 'ValueTable< Derivatives< %d, %d, %d, %d > >' % (row.dim,
#                 row.range, row.rank, d)
#         strings.append('template void %s::evaluate_bspline_derivatives<%d>(const FuncPointSize &,
#             StaticMultiArray<std::array<BasisValues1d*, dim_domain>, dim_range, rank> &elem_values,
#             ValueTable< DerivativeRef_t<deriv_order> > &derivatives_phi_hat%s &, const int) const ;\n' %
#                 (CA,d,VT))


for s in unique(strings): # Removing repeated entries.
    file_output.write(s)


file_output.write('IGA_NAMESPACE_CLOSE\n')

file_output.close()


