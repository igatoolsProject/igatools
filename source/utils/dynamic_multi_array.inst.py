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

include_files = ['base/tensor.h','basis_functions/bspline_element_accessor.h']
data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)

ma_list = ['DynamicMultiArray<TensorIndex<%s>,%s>' %(dim,dim) 
           for dim in inst.domain_dims]
ma_list = ma_list + ['DynamicMultiArray<%s,%s>' % (t,dim)
                     for  dim in inst.domain_dims for t in ('Real','Index')]
ma_list = ma_list + ['DynamicMultiArray<%s,2>' %(deriv)
           for deriv in inst.derivatives + inst.values]

for row in ma_list:
    f.write('template class %s; \n' % (row))
    f.write('template LogStream &operator<<(LogStream &, const %s &); \n' % (row))
    

      
      
evaluators = set(['DynamicMultiArray<std::shared_ptr<BSplineElementScalarEvaluator<%d>>,%d>' %(x.dim,x.dim)
              for x in inst.all_ref_sp_dims])
for eval in evaluators:
   f.write('template class %s ;\n' %eval)
