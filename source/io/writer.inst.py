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

strings = []
spaces = ['BSplineSpace', 'NURBSSpace']
writer_real_types = ['double','float']



for row in inst.user_phy_sp_dims:
    for writer_real_t in writer_real_types:
        writer = 'Writer<%d, %d, %s>' %(row.dim, row.codim, writer_real_t)
        strings.append('template class %s ;\n' % (writer))
for s in unique(strings): # Removing repeated entries.
    f.write(s)



strings = []
for row in inst.user_phy_sp_dims:
    for writer_real_t in writer_real_types:
        writer = 'Writer<%d, %d, %s>' %(row.dim, row.codim, writer_real_t)
        for name in spaces:
            space_ref  = '%s<%d,%d,%d>' % (name, row.dim, row.range, row.rank)
            PushForward = 'PushForward<Transformation::%s,%d,%d>' %(row.trans_type, row.dim, row.codim)
            space_phys = 'PhysicalSpace<%s,%s>' %(space_ref,PushForward)
            func = 'add_field<%s>(shared_ptr<%s>, const Vector<LinAlgebra> &, const string & )' % (space_phys,space_phys)
            strings.append('template void %s::%s ;\n' % (writer,func))
            func = 'add_field<%s>(shared_ptr<%s>, const Vector<LinAlgebra> &, const string & )' % (space_ref,space_ref)
            strings.append('template void %s::%s ;\n' % (writer,func))


############################################
# TRILINOS specific instantiations -- begin
f.write('#ifdef USE_TRILINOS\n')
for s in unique(strings): # Removing repeated entries.
    f.write(s.replace('LinAlgebra','LAPack::trilinos'))
f.write('#endif\n')
# TRILINOS specific instantiations -- end
############################################


############################################
# PETSC specific instantiations -- begin
f.write('#ifdef USE_PETSC\n')
for s in unique(strings): # Removing repeated entries.
    f.write(s.replace('LinAlgebra','LAPack::petsc'))
f.write('#endif\n')
# PETSC specific instantiations -- end
############################################

