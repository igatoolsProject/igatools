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

for dim in inst.ref_dom_dims:
    f.write('template class Quadrature<%d>; \n' %dim)
    for k in range(dim+1, max(0, dim - inst.n_sub_element), -1):
        f.write('template Quadrature<%d> Quadrature<%d>::' %(dim, dim) +
                'collapse_to_sub_element<%d>(const int id) const; \n' %(k-1) )
          
# for dim in inst.face_ref_dom_dims:
#     f.write('template Quadrature<%d> extend_face_quad<%d>' %(dim+1, dim) +
#             '(const Quadrature <%d> &, const int);\n' %(dim))
