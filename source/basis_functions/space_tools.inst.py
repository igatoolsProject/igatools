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
include_files = ['basis_functions/bspline_space.h',
#                  'basis_functions/bspline_element_accessor.h',
#                  'basis_functions/bspline_uniform_quad_cache.h',
#                  'basis_functions/nurbs_space.h',
#                  'basis_functions/nurbs_element_accessor.h',
#                  'basis_functions/nurbs_uniform_quad_cache.h',
#                  'basis_functions/physical_space.h',
#                  'geometry/cartesian_grid_element.h',
#                  'geometry/mapping_element_accessor.h',
#                  'geometry/push_forward_element_accessor.h',
#                  'geometry/push_forward_uniform_quad_cache.h',
#                  'basis_functions/physical_space_element_accessor.h',
#                  'basis_functions/space_uniform_quad_cache.h'
]
data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)

 
# projection_l2=('template Vector<LinAlgebra> '
#         'space_tools::projection_l2<Space,LinAlgebra>('
#         'const typename Space::Func &f,'
#         'std::shared_ptr<const Space> phys_space,'
#         'const QuadratureTensorProduct<Space::dim> & );\n')
# 
# project_boundary_values_1=('template void space_tools::project_boundary_values<Space,LinAlgebra>('
#         'const typename Space::Func &,'
#         'std::shared_ptr<const Space> ,'
#         'const QuadratureTensorProduct<Space::dim-1> &,'
#         'const std::set<boundary_id>  &,'
#         'std::map<Index, Real>  &);\n')
# 
# project_boundary_values_2=('template void space_tools::project_boundary_values<Space,LinAlgebra>('
#         'const typename Space::Func &,'
#         'std::shared_ptr<const Space> ,'
#         'const QuadratureTensorProduct<Space::dim-1> &,'
#         'const boundary_id ,'
#         'std::map<Index, Real>  &);\n')
# 
# integrate_difference=('template Real space_tools::integrate_difference('
#         'const typename Space::Func & ,'
#         'std::shared_ptr<const Space> ,'
#         'const QuadratureTensorProduct< Space::dim > &,'
#         'const Norm &,'
#         'const Vector<LinAlgebra> &,'
#         'vector< Real > &);\n');
# 
# 
# 
# ############################################
# # TRILINOS specific instantiations -- begin
# f.write('#ifdef USE_TRILINOS\n')
# for sp in inst.PhysSpaces + inst.UserRefSpaces:
#     f.write(projection_l2.replace('Space',sp).replace('LinAlgebra','LAPack::trilinos'))
# 
# 
# for sp in inst.UserPhysSpaces + inst.UserRefSpaces:
#    f.write(project_boundary_values_1.replace('Space',sp).replace('LinAlgebra','LAPack::trilinos'))
#   
#  
# for sp in inst.UserPhysSpaces + inst.UserRefSpaces:
#     f.write(project_boundary_values_2.replace('Space',sp).replace('LinAlgebra','LAPack::trilinos'))
#  
# for sp in inst.UserPhysSpaces + inst.UserRefSpaces:
#     f.write(integrate_difference.replace('Space',sp).replace('LinAlgebra','LAPack::trilinos'))
# 
# f.write('#endif\n')
# # TRILINOS specific instantiations -- end
# ############################################
# 
# 
# 
# ############################################
# # PETSC specific instantiations -- begin
# f.write('#ifdef USE_PETSC\n')
# for sp in inst.PhysSpaces + inst.RefSpaces:
#     f.write(projection_l2.replace('Space',sp).replace('LinAlgebra','LAPack::petsc'))
# 
# 
# # for sp in inst.UserFilteredRefSpaces + inst.UserPhysSpaces:
# #     f.write(project_boundary_values_1.replace('Space',sp).replace('LinAlgebra','LAPack::petsc'))
# # 
# # 
# # for sp in inst.UserFilteredRefSpaces + inst.UserPhysSpaces:
# #     f.write(project_boundary_values_2.replace('Space',sp).replace('LinAlgebra','LAPack::petsc'))
# # 
# # 
# # for sp in inst.UserPhysSpaces + inst.UserRefSpaces:
# #     f.write(integrate_difference.replace('Space',sp).replace('LinAlgebra','LAPack::petsc'))
# f.write('#endif\n')
# # PETSC specific instantiations -- end
# ############################################
