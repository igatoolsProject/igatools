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

file_output, inst = intialize_instantiation()

file_output.write('IGA_NAMESPACE_OPEN\n')

# instantiating Quadrature
for row in inst.quadratures:
    file_output.write('template class %s; \n' % (row))

file_output.write('\n')

file_output.write('IGA_NAMESPACE_CLOSE\n')

file_output.close()


# ###############################################################################
# # Common header for instantiation files 
# from igatools_instantiation import *
# file_output, inst = intialize_instantiation()
# ###############################################################################
# 
# for dim in inst.ref_dom_dims:
#     file_output.write('template class Quadrature< %d > ;\n' % (dim))
#     file_output.write('template Quadrature< %d > restricted_quad' %(dim) +
#                       '(const Quadrature< %d > &, const int) ;\n' % (dim))
# 
# for dim in inst.face_ref_dom_dims:
#     file_output.write('template Quadrature< %d > extend_quad_dim< %d >' % (dim+1, dim) +
#                       '(const Quadrature < %d > &, const int);\n' %(dim))
# 
# file_output.close()

