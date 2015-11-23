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

# include_files = ['geometry/grid_element.h']
include_files = []
data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)
    
# Grids
dims = []
for dim in inst.domain_dims:
    dims.append(dim)
for dim in inst.sub_domain_dims:
    dims.append(dim)

for dim in unique(dims):
    grid_ptr = 'std::shared_ptr<Grid<%d>>' % (dim)
    f.write('template bool ObjectsContainer::is_grid<%d>(const Index&) const;\n' % (dim))
    f.write('template %s ObjectsContainer::get_grid<%d>(const Index&) const;\n' % (grid_ptr, dim))
    f.write('template void ObjectsContainer::add_grid<%d>(const %s, const Index&);\n' % (dim, grid_ptr))
f.write('\n')


space_dims = [[0, 0, 1]]
for x in inst.sub_ref_sp_dims:
    space_dims.append([x.dim, x.range, x.rank])
for x in inst.ref_sp_dims:
    space_dims.append([x.dim, x.range, x.rank])
    
space_dims = unique(space_dims);

# Spline spaces
for dims in space_dims:
    dims_str = '%d, %d, %d' % (dims[0], dims[1], dims[2])
    ss_ptr = 'std::shared_ptr<SplineSpace<%s>>' % (dims_str)
    f.write('template bool ObjectsContainer::is_spline_space<%s>(const Index&) const;\n' % (dims_str))
    f.write('template %s ObjectsContainer::get_spline_space<%s>(const Index&) const;\n' % (ss_ptr, dims_str))
    f.write('template void ObjectsContainer::add_spline_space<%s>(const %s, const Index&);\n' % (dims_str, ss_ptr))
f.write('\n')

# Reference spaces
for dims in space_dims:
    dims_str = '%d, %d, %d' % (dims[0], dims[1], dims[2])
    rs_ptr = 'std::shared_ptr<ReferenceSpaceBasis<%s>>' % (dims_str)
    f.write('template bool ObjectsContainer::is_ref_space<%s>(const Index&) const;\n' % (dims_str))
    f.write('template %s ObjectsContainer::get_ref_space<%s>(const Index&) const;\n' % (rs_ptr, dims_str))
    f.write('template void ObjectsContainer::add_ref_space<%s>(const %s, const Index&);\n' % (dims_str, rs_ptr))
f.write('\n')


# BSpline spaces
for dims in space_dims:
    dims_str = '%d, %d, %d' % (dims[0], dims[1], dims[2])
    bs_ptr = 'std::shared_ptr<BSpline<%s>>' % (dims_str)
    f.write('template bool ObjectsContainer::is_bspline<%s>(const Index&) const;\n' % (dims_str))
    f.write('template %s ObjectsContainer::get_bspline<%s>(const Index&) const;\n' % (bs_ptr, dims_str))
f.write('\n')


# NURBS spaces
for dims in space_dims:
    dims_str = '%d, %d, %d' % (dims[0], dims[1], dims[2])
    nr_ptr = 'std::shared_ptr<NURBS<%s>>' % (dims_str)
    f.write('template bool ObjectsContainer::is_nurbs<%s>(const Index&) const;\n' % (dims_str))
    f.write('template %s ObjectsContainer::get_nurbs<%s>(const Index&) const;\n' % (nr_ptr, dims_str))
f.write('\n')

gf_dims = []
for x in inst.sub_mapping_dims:
    gf_dims.append([x.dim, x.space_dim])
for x in inst.mapping_dims:
    gf_dims.append([x.dim, x.space_dim])
    # The next dimensions are needed by NURBS
    gf_dims.append([x.dim, 1])

# Grid functions
for dims in unique(gf_dims):
    dims_str = '%d, %d' % (dims[0], dims[1])
    gf_ptr = 'std::shared_ptr<GridFunction<%s>>' % (dims_str)
    f.write('template bool ObjectsContainer::is_grid_function<%s>(const Index&) const;\n' % (dims_str))
    f.write('template %s ObjectsContainer::get_grid_function<%s>(const Index&) const;\n' % (gf_ptr, dims_str))
    f.write('template void ObjectsContainer::add_grid_function<%s>(const %s, const Index&);\n' % (dims_str, gf_ptr))
f.write('\n')


# Domain
dm_dims = []
for x in inst.sub_mapping_dims:
    dm_dims.append([x.dim, x.codim])
for x in inst.mapping_dims:
    dm_dims.append([x.dim, x.codim])

for dims in unique(dm_dims):
    dims_str = '%d, %d' % (dims[0], dims[1])
    dm_ptr = 'std::shared_ptr<Domain<%s>>' % (dims_str)
    f.write('template bool ObjectsContainer::is_domain<%s>(const Index&) const;\n' % (dims_str))
    f.write('template %s ObjectsContainer::get_domain<%s>(const Index&) const;\n' % (dm_ptr, dims_str))
    f.write('template void ObjectsContainer::add_domain<%s>(const %s, const Index&);\n' % (dims_str, dm_ptr))
f.write('\n')


# Physical spaces
sp_specs = [[0, 0, 1, 0]]
for sp in inst.SubPhysSpaces:
    sp_specs.append([sp.spec.dim, sp.spec.range, sp.spec.rank, sp.spec.codim])
for sp in inst.PhysSpaces:
    sp_specs.append([sp.spec.dim, sp.spec.range, sp.spec.rank, sp.spec.codim])

for dims in unique(sp_specs):
    dims_str = '%d, %d, %d, %d' % (dims[0], dims[1], dims[2], dims[3])
    ps_ptr = 'std::shared_ptr<PhysicalSpaceBasis<%s>>' % (dims_str)
    f.write('template bool ObjectsContainer::is_phys_space_basis<%s>(const Index&) const;\n' % (dims_str))
    f.write('template %s ObjectsContainer::get_phys_space_basis<%s>(const Index&) const;\n' % (ps_ptr, dims_str))
    f.write('template void ObjectsContainer::add_phys_space_basis<%s>(const %s, const Index&);\n' % (dims_str, ps_ptr))
f.write('\n')


# Functions
for dims in inst.all_function_dims:
    dims_str = '%d, %d, %d, %d' % (dims.dim, dims.codim, dims.range, dims.rank)
    fn_ptr = 'std::shared_ptr<Function<%s>>' % (dims_str)
    f.write('template bool ObjectsContainer::is_function<%s>(const Index&) const;\n' % (dims_str))
    f.write('template %s ObjectsContainer::get_function<%s>(const Index&) const;\n' % (fn_ptr, dims_str))
    f.write('template void ObjectsContainer::add_function<%s>(const %s, const Index&);\n' % (dims_str, fn_ptr))
f.write('\n')