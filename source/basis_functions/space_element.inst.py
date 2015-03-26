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

include_files = ['basis_functions/space_element.h']
data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)


sub_dim_members = \
[]


space_elem_list = []

#--------------------------------------------------------------------------------------
# SpaceElement used by ReferenceSpaceElement 
for x in inst.sub_ref_sp_dims:
    space_elem = 'SpaceElement<%d,0,%d,%d>' %(x.dim, x.range, x.rank)
    space_elem_list.append(space_elem)



for x in inst.ref_sp_dims:
    space_elem = 'SpaceElement<%d,0,%d,%d>' %(x.dim, x.range, x.rank)
    space_elem_list.append(space_elem)
#--------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------
# SpaceElement used by PhysicalSpaceElement 
for space in inst.SubPhysSpaces:
    x = space.spec
    space_elem = 'SpaceElement<%d,%d,%d,%d>' %(x.dim,x.codim,x.range, x.rank)
    space_elem_list.append(space_elem)


for space in inst.PhysSpaces:
    x = space.spec
    space_elem = 'SpaceElement<%d,%d,%d,%d>' %(x.dim,x.codim,x.range, x.rank)
    space_elem_list.append(space_elem)
#--------------------------------------------------------------------------------------



for space_elem in unique(space_elem_list):
    f.write('template class %s ;\n' %space_elem)


