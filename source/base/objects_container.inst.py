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

include_files = []
data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)


# Creating all possible types in the container.
valid_types = []

# Grids
for dim in inst.domain_dims:
    valid_types.append('Grid<%d>' % (dim))
for dim in inst.sub_domain_dims:
    valid_types.append('Grid<%d>' % (dim))

# Spline spaces
valid_types.append('SplineSpace<%d, %d, %d>' % (0, 0, 1))
for x in inst.sub_ref_sp_dims:
    valid_types.append('SplineSpace<%d, %d, %d>' % (x.dim, x.range, x.rank))
for x in inst.ref_sp_dims:
    valid_types.append('SplineSpace<%d, %d, %d>' % (x.dim, x.range, x.rank))

# Reference spaces
valid_types.append('ReferenceBasis<%d, %d, %d>' % (0, 0, 1))
for x in inst.sub_ref_sp_dims:
    valid_types.append('ReferenceBasis<%d, %d, %d>' % (x.dim, x.range, x.rank))
for x in inst.ref_sp_dims:
    valid_types.append('ReferenceBasis<%d, %d, %d>' % (x.dim, x.range, x.rank))

# Grid functions
for x in inst.sub_mapping_dims:
    valid_types.append('GridFunction<%d, %d>' % (x.dim, x.space_dim))
for x in inst.mapping_dims:
    valid_types.append('GridFunction<%d, %d>' % (x.dim, x.space_dim))
    # The next dimensions are needed by NURBS
    valid_types.append('GridFunction<%d, %d>' % (x.dim, 1))

# Domains
for x in inst.sub_mapping_dims:
    valid_types.append('Domain<%d, %d>' % (x.dim, x.codim))
for x in inst.mapping_dims:
    valid_types.append('Domain<%d, %d>' % (x.dim, x.codim))

# Physical spaces
valid_types.append('PhysicalBasis<%d, %d, %d, %d>' % (0, 0, 1, 0))
for sp in inst.SubPhysSpaces:
    valid_types.append('PhysicalBasis<%d, %d, %d, %d>' % (sp.spec.dim, sp.spec.range, sp.spec.rank, sp.spec.codim))
for sp in inst.PhysSpaces:
    valid_types.append('PhysicalBasis<%d, %d, %d, %d>' % (sp.spec.dim, sp.spec.range, sp.spec.rank, sp.spec.codim))

# Functions
for dims in inst.all_function_dims:
    valid_types.append('Function<%d, %d, %d, %d>' % (dims.dim, dims.codim, dims.range, dims.rank))


valid_types = unique(valid_types)

# Writing to file

for tp in valid_types:
    f.write('template bool ObjectsContainer::is_object_present<%s>(const Index&) const;\n' % (tp))
    f.write('template bool ObjectsContainer::is_const_object_present<%s>(const Index&) const;\n' % (tp))
f.write('\n')

for tp in valid_types:
    f.write('template std::shared_ptr<%s> ObjectsContainer::get_object<%s>(const Index&) const;\n' % (tp, tp))
    f.write('template std::shared_ptr<const %s> ObjectsContainer::get_const_object<%s>(const Index&) const;\n' % (tp, tp))
f.write('\n')

for tp in valid_types:
    f.write('template SafeSTLSet<Index> ObjectsContainer::get_object_ids<%s>() const;\n' % (tp))
    f.write('template SafeSTLSet<Index> ObjectsContainer::get_const_object_ids<%s>() const;\n' % (tp))
f.write('\n')

for tp in valid_types:
    f.write('template void ObjectsContainer::insert_object<%s>(const std::shared_ptr<%s>, const bool);\n' % (tp, tp))
    f.write('template void ObjectsContainer::insert_const_object<%s>(const std::shared_ptr<const %s>, const bool);\n' % (tp, tp))
