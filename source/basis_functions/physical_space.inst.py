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

# QA (pauletti, Mar 19, 2014):
from init_instantiation_data import *

include_files = ['basis_functions/bspline_space.h',
                 'basis_functions/nurbs_space.h',
                 'geometry/push_forward.h',
                 'geometry/cartesian_grid_element.h',
                 'geometry/mapping_element.h',
                 'geometry/push_forward_element.h',
                 'basis_functions/bspline_element.h',
                 'basis_functions/nurbs_element.h',
                 '../../source/basis_functions/physical_space_element.cpp',
                 '../../source/basis_functions/bspline_element.cpp',
                 '../../source/basis_functions/nurbs_element.cpp',
                 'basis_functions/physical_space_element.h']
#                 '../../source/geometry/cartesian_grid_iterator.cpp']
data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)


sub_dim_members = \
 ['std::shared_ptr<typename class::SubSpace<k> > ' +
   'class::get_sub_space(const int s_id, InterSpaceMap<k> &dof_map, ' +
   'std::shared_ptr<CartesianGrid<k>> sub_grid, ' + 
   'std::shared_ptr<InterGridMap<k>> elem_map) const;']


# for x in inst.sub_ref_sp_dims:
#     space = 'BSplineSpace<%d, %d, %d>' %(x.dim, x.range, x.rank)
#     f.write('template class %s ;\n' %space)
#     for fun in sub_dim_members:
#         k = x.dim
#         s = fun.replace('class', space).replace('k', '%d' % (k));
#         f.write('template ' + s + '\n')


# for x in inst.ref_sp_dims:
#     space = 'BSplineSpace<%d, %d, %d>' %(x.dim, x.range, x.rank)
#     f.write('template class %s ;\n' %space)
#     for fun in sub_dim_members:
#         for k in inst.sub_dims(x.dim):
#             s = fun.replace('class', space).replace('k', '%d' % (k));
#             f.write('template ' + s + '\n')
                        
 
for space in inst.SubPhysSpaces:
    x = space.spec
    f.write( 'template class %s;\n' %space.name)
    for fun in sub_dim_members:
        k = x.dim
        s = fun.replace('class', space.name).replace('k', '%d' % (k));
        f.write('template ' + s + '\n')


for space in inst.PhysSpaces:
    x = space.spec
    f.write( 'template class %s;\n' %space.name)
    for fun in sub_dim_members:
        for k in inst.sub_dims(x.dim):
            s = fun.replace('class', space.name).replace('k', '%d' % (k));
            f.write('template ' + s + '\n')
