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
include_files = ['basis_functions/bspline_space.h',
                 'basis_functions/bspline_element_accessor.h',
                 'basis_functions/nurbs_space.h',
                 'basis_functions/nurbs_element_accessor.h',
                 'basis_functions/physical_space.h',
                 'geometry/cartesian_grid_element_accessor.h',
                 'geometry/mapping_element_accessor.h',
                 'geometry/push_forward_element_accessor.h',
                 'basis_functions/physical_space_element_accessor.h']
data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)


master=('template Vector<LinearAlgebraPackage::trilinos> space_tools::projection_l2('
        'const Function<Space::space_dim,Space::range,Space::rank> &f,'
        'std::shared_ptr<const Space> phys_space,'
        'const Quadrature<Space::dim> & );\n')
for sp in inst.PhysSpaces + inst.RefSpaces:
    f.write(master.replace('Space', sp))


master=('template void space_tools::project_boundary_values('
        'const Function<Space::space_dim,Space::range,Space::rank> &,'
        'std::shared_ptr<const Space> ,'
        'const Quadrature<Space::dim-1> &,'
        'const std::set<boundary_id>  &,'
        'std::map<Index, Real>  &);\n')
for sp in inst.UserFilteredRefSpaces + inst.UserPhysSpaces:
    f.write(master.replace('Space', sp))


master=('template void space_tools::project_boundary_values('
        'const Func<Space> &,'
        'std::shared_ptr<const Space> ,'
        'const Quadrature<Space::dim-1> &,'
        'const boundary_id ,'
        'std::map<Index, Real>  &);\n')
for sp in inst.UserFilteredRefSpaces + inst.UserPhysSpaces:
    f.write(master.replace('Space', sp))



master=('template Real space_tools::integrate_difference('
        'std::shared_ptr<const Func<Space> > ,'
        'std::shared_ptr<const Space> ,'
        'const Quadrature< Space::dim > &,'
        'const Norm &,'
        'const Vector<LinearAlgebraPackage::trilinos> &,'
        'std::vector< Real > &);\n');
for sp in inst.UserPhysSpaces + inst.UserRefSpaces:
    f.write(master.replace('Space', sp))
    
master=('template std::shared_ptr< FaceSpace<SP> >'
        'space_tools::get_face_space(const std::shared_ptr<const SP> ,'
        'const int ,'
        'std::vector<Index> &);\n');
for sp in inst.UserPhysSpaces + inst.UserFilteredRefSpaces:
    f.write(master.replace('SP', sp))   
