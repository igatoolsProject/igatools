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

# QA needs reviewing
from init_instantiation_data import *

include_files = ['base/tensor.h',
                 'utils/container_view.h',
                 'utils/concatenated_iterator.h',
                 'utils/array.h']
data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)

ma_list = ['DynamicMultiArray<TensorIndex<%s>,%s>' %(dim,dim) 
           for dim in inst.all_domain_dims]
ma_list = ma_list + ['DynamicMultiArray<%s,%s>' % (t,dim)
                     for  dim in inst.all_domain_dims for t in ('Real','Index','bool','vector<Index>')]
normals = ['special_array<Points<%d>, %d>' %(x.dim, x.codim) for x in inst.all_mapping_dims]
curvatures = ['vector<Real>']

ma_list = ma_list + ['DynamicMultiArray<%s,2>' %(deriv)
           for deriv in inst.derivatives + inst.values + inst.divs + normals + curvatures]

    		


for row in ma_list:
    f.write('template class %s; \n' % (row))


# needed by DofDistribution
Vec = 'vector<Index>'
ContView = ('ContainerView<%s>' % Vec)
ConstContView = ('ConstContainerView<%s>' % Vec)
It = ('ConcatenatedIterator<%s>' % ContView)
ConstIt = ('ConcatenatedConstIterator<%s,%s>' %(ContView,ConstContView))
ConstView = ('ConstView<%s,%s>' %(It,ConstIt))
dof_views = set(['DynamicMultiArray<%s,%d>' %(ConstView,dim)
                     for dim in inst.domain_dims])
for dof_view in dof_views:
   f.write('template class %s ;\n' %dof_view)

