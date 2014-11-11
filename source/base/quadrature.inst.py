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

# QA (pauletti, Mar 19, 2014):
from init_instantiation_data import *
data = Instantiation()
(f, inst) = (data.file_output, data.inst)

for dim in inst.all_domain_dims:
    f.write('template class Quadrature<%d>; \n' %dim)
        
for dim in inst.domain_dims:
    for k in inst.sub_dims(dim):
        f.write('template Quadrature<%d> Quadrature<%d>::' %(dim, dim) +
                'collapse_to_sub_element<%d>(const int id) const; \n' %(k) )

ext_members = [ 'Quadrature<dim> extend_sub_elem_quad<k,dim>(const Quadrature<k> &quad, const int sub_elem_id);']   
for dim in inst.domain_dims:
    for fun in ext_members:
        for k in range(0, dim+1):
            s = fun.replace('dim','%d' %dim).replace('k','%d' %k);
            f.write('template ' + s + '\n')

