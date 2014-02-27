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


include_files = [ '#include <igatools/geometry/push_forward.h>\n',
                  '#include <igatools/geometry/cartesian_grid_element_accessor.h>\n'
                  '#include <igatools/geometry/mapping_element_accessor.h>\n']

for include in include_files:
    file_output.write(include)


file_output.write('IGA_NAMESPACE_OPEN\n')

output = []
#todo should be read form table
transformation_types = ['Transformation::h_grad']#, 'Transformation::h_div', 'Transformation::h_curl', 'Transformation::l_2']

containers = ['ValueTable', 'ValueVector']

for transformation_name in transformation_types:

# Trasformation for values
    for row in inst.table:
        PF = 'PushForward<%s,%d,%d>' %(transformation_name,row.dim_ref,row.dim_phys-row.dim_ref)
        push_fwd_elem_acc = 'PushForwardElementAccessor<%s>' %(PF)
        output.append('template class %s ;\n' %(push_fwd_elem_acc) )
        value_ref = ("Values<%d,%d,%d>" %(row.dim_ref, row.range_ref, row.rank_ref)) 
        value_phys = ("Values<%d,%d,%d>" %(row.dim_phys, row.range_ref, row.rank_ref)) 
        deriv_ref = ("Derivatives<%d,%d,%d,1>" %(row.dim_ref, row.range_ref, row.rank_ref)) 
        deriv_phys = ("Derivatives<%d,%d,%d,1>" %(row.dim_phys, row.range_ref, row.rank_ref)) 
        for container in containers:
            v_ref  = 'const %s<%s> &' %(container, value_ref)
            v_phys = '%s<%s> &' %(container, value_phys)
            output.append(
                    'template void %s::' %(push_fwd_elem_acc) +
                    'transform_values<%d,%d,%s,%s>' %(row.range_ref, row.rank_ref,container,transformation_name) +
                    '( %s, %s, void *) const ;\n' %(v_ref, v_phys)
                    )
            output.append(
                    'template void %s::' %(push_fwd_elem_acc) +
                    'transform_face_values<%d,%d,%s,%s>' %(row.range_ref, row.rank_ref,container,transformation_name) +
                    '( const Index, %s, %s, void *) const ;\n' %(v_ref, v_phys)
                    )

# Trasformation for gradients
            dv_ref  = 'const %s<%s> &' %(container, deriv_ref)
            dv_phys = '%s<%s> &' %(container, deriv_phys)
            output.append(
                    'template void %s::' %(push_fwd_elem_acc) +
                    'transform_gradients<%d,%d,%s,%s>' %(row.range_ref, row.rank_ref,container, transformation_name) +
                    '(%s,%s,%s, void * ) const;\n' %(v_ref,dv_ref,dv_phys)
                    )
            output.append(
                    'template void %s::' %(push_fwd_elem_acc) +
                    'transform_face_gradients<%d,%d,%s,%s>' %(row.range_ref, row.rank_ref,container, transformation_name) +
                    '(const Index, %s,%s,%s, void * ) const;\n' %(v_ref,dv_ref,dv_phys)
                    )
 


for s in unique(output):
    file_output.write(s)


file_output.write('IGA_NAMESPACE_CLOSE\n')


file_output.close()
