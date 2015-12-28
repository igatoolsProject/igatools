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
include_files = ['utils/element_index.h',
                 'utils/safe_stl_array.h',
                 'utils/safe_stl_vector.h',
                 'basis_functions/bernstein_extraction.h']
data = Instantiation(include_files)
f = data.file_output
inst = data.inst

containers = []

containers.append('std::set<int>')
containers.append('std::vector<int>')
containers.append('std::vector<Real>')
containers.append('std::vector<BernsteinOperator>')
containers.append('std::vector<const BernsteinOperator *>')
containers.append('std::vector<SafeSTLVector<BernsteinOperator>>')
containers.append('std::vector<SafeSTLVector<const BernsteinOperator *>>')
#containers.append('SafeSTLContainer<std::vector<SafeSTLVector<const BernsteinOperator*>>>')
containers.append('std::vector<const SafeSTLVector<BernsteinOperator> *>');

for dim in inst.all_domain_dims:
    t1 = 'ElementIndex<%d>' % (dim)
    containers.append('std::vector<%s>' % (t1))
    t2 = 'TensorIndex<%d>' % (dim)
    containers.append('std::vector<%s>' % (t2))
    t3 = 'SafeSTLArray<Real,2>'
    containers.append('std::array<%s,%d>' % (t3,dim))
    t4 = 'int'
    containers.append('std::array<%s,%d>' % (t4,dim))
#    t5 = 'std::pair<Real,Real>'
#    containers.append('std::array<%s,%d>' % (t5,dim))



for c in unique(containers):
    f.write('template class SafeSTLContainer<%s>;\n' %(c))
#for deriv in inst.derivatives + inst.values + inst.divs:
#    value_vectors.append('ValueVector<%s>' % (deriv))

#for row in set (value_vectors):
#    f.write('template class %s; \n' % (row))
#    f.write("template %s operator*(const Real, const %s &) ;\n" % (row,row))
#    f.write("template %s operator*(const %s &, const Real) ;\n" % (row,row))
    
#for row in set (normals + curvatures):
#      f.write('template class ValueVector<%s>; \n' % (row))
