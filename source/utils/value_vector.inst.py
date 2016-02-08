#-+--------------------------------------------------------------------
# Igatools a general purpose Isogeometric analysis library.
# Copyright (C) 2012-2016  by the igatools authors (see authors.txt).
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
include_files = ['utils/safe_stl_array.h','/base/tensor.h']
data = Instantiation(include_files)
f = data.file_output
inst = data.inst

normals = ['SafeSTLArray<Points<%d>, %d>' %(x.space_dim, x.codim) for x in inst.sub_mapping_dims + inst.mapping_dims]
curvatures = ['SafeSTLVector<Real>']    
value_vectors=['ValueVector<Real>']

for deriv in inst.derivatives + inst.values + inst.divs:
    value_vectors.append('ValueVector<%s>' % (deriv))

for row in set (value_vectors):
    f.write('template class %s; \n' % (row))
    f.write("template %s operator*(const Real, const %s &) ;\n" % (row,row))
    f.write("template %s operator*(const %s &, const Real) ;\n" % (row,row))
    
for row in set (normals + curvatures):
      f.write('template class ValueVector<%s>; \n' % (row))
