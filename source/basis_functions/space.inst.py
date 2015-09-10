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

spaces = []

#for dim in inst.all_domain_dims:
#    space = 'SpaceBase<%d>' %(dim)
#    spaces. append(space)
#    f.write("template class %s ;\n" %(space))


#---------------------------------------------------
#f.write('IGA_NAMESPACE_CLOSE\n')
 
#f.write('#ifdef SERIALIZATION\n')
#id = 0 
#for space in unique(spaces):
#    alias = 'SpaceBaseAlias%d' %(id)
#    f.write('using %s = iga::%s; \n' % (alias, space))
#    f.write('BOOST_CLASS_EXPORT_IMPLEMENT(%s) \n' %alias)
#    f.write('template void %s::serialize(OArchive &, const unsigned int);\n' % alias)
#    f.write('template void %s::serialize(IArchive &, const unsigned int);\n' % alias)
#    id += 1 
#f.write('#endif // SERIALIZATION\n')
    
#f.write('IGA_NAMESPACE_OPEN\n')
#---------------------------------------------------







transformations = ['Transformation::h_grad']

spaces = ['Space<0,0,0,1,Transformation::h_grad>']

for x in inst.all_phy_sp_dims:
    for t in transformations:
        space = 'Space<%d,%d,%d,%d,%s>' %(x.dim,x.codim,x.range,x.rank,t)
        spaces.append(space)

for x in inst.all_ref_sp_dims:
    for t in transformations:
        space = 'Space<%d,0,%d,%d,%s>' %(x.dim,x.range,x.rank,t)
        spaces.append(space)
    
for space in unique(spaces):
    f.write("template class %s ;\n" %(space))


#---------------------------------------------------
#f.write('IGA_NAMESPACE_CLOSE\n')
 
#f.write('#ifdef SERIALIZATION\n')
#id = 0 
#for space in unique(spaces):
#    alias = 'SpaceAlias%d' %(id)
    
#    f.write('using %s = iga::%s; \n' % (alias, space))
#    f.write('using %s = iga::%s; \n' % (alias, space.replace('Transformation','iga::Transformation')))

#    f.write('ALLOW_SHARED_THIS(%s)\n' %alias )
    
#    f.write('BOOST_CLASS_EXPORT_IMPLEMENT(%s) \n' %alias)
#    f.write('template void %s::serialize(OArchive &, const unsigned int);\n' % alias)
#    f.write('template void %s::serialize(IArchive &, const unsigned int);\n' % alias)
#    id += 1 
#f.write('#endif // SERIALIZATION\n')
    
#f.write('IGA_NAMESPACE_OPEN\n')
#---------------------------------------------------
