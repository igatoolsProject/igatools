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

# QA (pauletti, Mar 19, 2014): review the transform value how is handle
from init_instantiation_data import *
include_files = ['geometry/push_forward.h',
                 'geometry/cartesian_grid_element_accessor.h',
                 'geometry/mapping_element_accessor.h',
                 '../../source/geometry/grid_forward_iterator.cpp']
data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)

output = []
containers = ['ValueTable', 'ValueVector']

# Trasformation for values
for row in unique(inst.all_table + inst.extended_table):
    PF = 'PushForward<Transformation::%s,%d,%d>' %(row.trans_type, row.dim, row.codim)
    push_fwd_elem_acc = 'PushForwardElementAccessor<%s>' %(PF)
    output.append('template class %s ;\n' %(push_fwd_elem_acc) )
    output.append('template class GridForwardIterator<%s> ;\n' %(push_fwd_elem_acc) )
    value_ref  = ("Values<%d,%d,%d>" %(row.dim, row.range, row.rank)) 
    gradient_ref  = ("Derivatives<%d,%d,%d,1>" %(row.dim, row.range, row.rank)) 
    hessian_ref  = ("Derivatives<%d,%d,%d,2>" %(row.dim, row.range, row.rank)) 
    value_phys = ("Values<%d,%d,%d>" %(row.space_dim, row.phys_range, row.phys_rank)) 
    gradient_phys = ("Derivatives<%d,%d,%d,1>" %(row.space_dim, row.phys_range, row.phys_rank)) 
    hessian_phys = ("Derivatives<%d,%d,%d,2>" %(row.space_dim, row.phys_range, row.phys_rank)) 
    topology_id = ("const TopologyId<%d> &" %(row.dim))
    pts = 'const vector<Point<%d>> &' %(row.dim)
    for container in containers:
        v_ref  = 'const %s<%s> &' %(container, value_ref)
        d1v_ref  = 'const %s<%s> &' %(container, gradient_ref)
        d2v_ref  = 'const %s<%s> &' %(container, hessian_ref)
        v_phys = '%s<%s> &' %(container, value_phys)
        d1v_phys = '%s<%s> &' %(container, gradient_phys)
        d2v_phys = '%s<%s> &' %(container, hessian_phys)
        output.append(
                    'template void %s::' %(push_fwd_elem_acc) +
                    'transform_values<%d,%d,%s,Transformation::%s>' %(row.range, row.rank, container,row.trans_type) +
                    '(%s,%s,%s,void *) const ;\n' %(v_ref, v_phys, topology_id)
                    )
        order = 1
        output.append(
            'template void %s::' %(push_fwd_elem_acc) +
            'transform_gradients<%d,%d,%s,Transformation::%s>' %(row.range, row.rank, container,row.trans_type) +
            '(%s,%s,%s,%s,void *) const;\n' %(v_ref,d1v_ref,d1v_phys,topology_id)
        )
        order = 2
        output.append(
            'template void %s::' %(push_fwd_elem_acc) +
            'transform_hessians<%d,%d,%s,Transformation::%s>' %(row.range, row.rank, container,row.trans_type) +
            '(%s,%s,%s,%s,%s,void *) const;\n' %(v_ref,d1v_ref,d2v_ref,d2v_phys,topology_id)
        )

        order = 0
        output.append(
            'template void %s::' %(push_fwd_elem_acc) +
            'transform_basis_derivatives_at_points<%d,%d,%s,Transformation::%s>' %(row.range,row.rank,container,row.trans_type) +
            '(%s,%s,%s,%s,%s,void *) const;\n' %(pts,v_ref,d1v_ref,d2v_ref,v_phys)
        )
        order = 1
        output.append(
            'template void %s::' %(push_fwd_elem_acc) +
            'transform_basis_derivatives_at_points<%d,%d,%s,Transformation::%s>' %(row.range,row.rank,container,row.trans_type) +
            '(%s,%s,%s,%s,%s,void *) const;\n' %(pts,v_ref,d1v_ref,d2v_ref,d1v_phys)
        )
        order = 2
        output.append(
            'template void %s::' %(push_fwd_elem_acc) +
            'transform_basis_derivatives_at_points<%d,%d,%s,Transformation::%s>' %(row.range,row.rank,container,row.trans_type) +
            '(%s,%s,%s,%s,%s,void *) const;\n' %(pts,v_ref,d1v_ref,d2v_ref,d2v_phys)
        )

 


for s in unique(output):
    f.write(s)


