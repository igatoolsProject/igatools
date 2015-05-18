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

include_files = ['../../source/basis_functions/space_element.cpp']
data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)


sub_dim_members = \
[]


handlers = ['SpaceElementHandler<0,0,0,1>']

#--------------------------------------------------------------------------------------
# SpaceElement used by ReferenceSpaceElement 
for x in inst.sub_ref_sp_dims + inst.ref_sp_dims:
    handler = 'SpaceElementHandler<%d,0,%d,%d>' %(x.dim, x.range, x.rank)
    handlers.append(handler)

#--------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------
# SpaceElement used by PhysicalSpaceElement 
for space in inst.SubPhysSpaces + inst.PhysSpaces:
    x = space.spec
    handler = 'SpaceElementHandler<%d,%d,%d,%d>' %(x.dim,x.codim,x.range, x.rank)
    handlers.append(handler)

#--------------------------------------------------------------------------------------



#---------------------------------------------------
for handler in unique(handlers):
    f.write('template class %s ;\n' %handler)



#---------------------------------------------------
# f.write('IGA_NAMESPACE_CLOSE\n')
# 
# f.write('#ifdef SERIALIZATION\n')
# id = 0 
# for elem in unique(handlers):
#     alias = 'SpaceElementHandlerAlias%d' %(id)
#     f.write('using %s = iga::%s; \n' % (alias, elem))
#     f.write('BOOST_CLASS_EXPORT_IMPLEMENT(%s) \n' %alias)
#     f.write('template void %s::serialize(OArchive &, const unsigned int);\n' % alias)
#     f.write('template void %s::serialize(IArchive &, const unsigned int);\n' % alias)
#     id += 1 
# f.write('#endif // SERIALIZATION\n')
#     
# f.write('IGA_NAMESPACE_OPEN\n')
#---------------------------------------------------


