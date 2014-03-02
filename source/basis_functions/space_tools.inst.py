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

include_files = ['#include <igatools/basis_functions/bspline_space.h>\n',
                 '#include <igatools/basis_functions/bspline_element_accessor.h>\n',
                 '#include <igatools/basis_functions/nurbs_space.h>\n',
                 '#include <igatools/basis_functions/nurbs_element_accessor.h>\n',
                 '#include <igatools/basis_functions/physical_space.h>\n',
                 '#include <igatools/geometry/cartesian_grid_element_accessor.h>\n'
                 '#include <igatools/geometry/mapping_element_accessor.h>\n',
                 '#include <igatools/geometry/push_forward_element_accessor.h>\n',
                 '#include <igatools/basis_functions/physical_space_element_accessor.h>\n']

for include in include_files:
    file_output.write(include)

file_output.write('IGA_NAMESPACE_OPEN\n')

master=('template Vector space_tools::projection_l2('
        'const Function<Space::space_dim,Space::dim_range,Space::rank> &f,'
        'std::shared_ptr<const Space> phys_space,'
        'const Quadrature<Space::dim> & );\n')
for sp in inst.PhysSpaces + inst.RefSpaces:
    file_output.write(master.replace('Space', sp))


master=('template void space_tools::project_boundary_values('
        'const Function<Space::space_dim,Space::dim_range,Space::rank> &,'
        'std::shared_ptr<const Space> ,'
        'const Quadrature<Space::dim-1> &,'
        'const std::set<boundary_id>  &,'
        'std::map<Index, Real>  &);\n')
for sp in inst.UserFilteredRefSpaces + inst.UserPhysSpaces:
    file_output.write(master.replace('Space', sp))


master=('template void space_tools::project_boundary_values('
        'const Func<Space> &,'
        'std::shared_ptr<const Space> ,'
        'const Quadrature<Space::dim-1> &,'
        'const boundary_id ,'
        'std::map<Index, Real>  &);\n')
for sp in inst.UserFilteredRefSpaces + inst.UserPhysSpaces:
    file_output.write(master.replace('Space', sp))



master=('template Real space_tools::integrate_difference('
        'std::shared_ptr<const Func<Space> > ,'
        'std::shared_ptr<const Space> ,'
        'const Quadrature< Space::dim > &,'
        'const Norm &,'
        'const Vector &,'
        'std::vector< Real > &);\n');
for sp in inst.UserPhysSpaces + inst.UserRefSpaces:
    file_output.write(master.replace('Space', sp))
    
master=('template std::shared_ptr< FaceSpace<SP> >'
        'space_tools::get_face_space(const std::shared_ptr<const SP> ,'
        'const int ,'
        'std::vector<Index> &);\n');
for sp in inst.UserPhysSpaces + inst.UserFilteredRefSpaces:
    file_output.write(master.replace('SP', sp))   
    
file_output.write('IGA_NAMESPACE_CLOSE\n')

file_output.close()
