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



transformations = ['Transformation::h_grad']

sub_dim_members = \
[]

handlers = []
handler_methods = []

handler = 'SpaceElementHandler<0,0,0,1,Transformation::h_grad>'
handlers.append(handler)

handler_method = 'void %s::init_cache<0>(ElementAccessor &)' % (handler)
handler_methods.append(handler_method)
handler_method = 'void %s::fill_cache<0>(ElementAccessor &, const int)' % (handler)
handler_methods.append(handler_method)

#--------------------------------------------------------------------------------------
# SpaceElement used by ReferenceSpaceElement 
for x in inst.sub_ref_sp_dims + inst.ref_sp_dims:
    for t in transformations:
        handler = 'SpaceElementHandler<%d,0,%d,%d,%s>' %(x.dim, x.range, x.rank,t)
        handlers.append(handler)
        for k in inst.sub_dims(x.dim):
#      handler_method = 'void %s::init_cache<%d>(SpaceElement<%d,0,%d,%d> &)' % (handler, k, x.dim, x.range, x.rank)
            handler_method = 'void %s::init_cache<%d>(ElementAccessor &)' % (handler, k)
            handler_methods.append(handler_method)
            handler_method = 'void %s::fill_cache<%d>(ElementAccessor &, const int)' % (handler, k)
            handler_methods.append(handler_method)

for x in inst.ref_sp_dims:
    for t in transformations:
        handler = 'SpaceElementHandler<%d,0,%d,%d,%s>' %(x.dim, x.range, x.rank,t)
        handlers.append(handler)
        for k in inst.sub_dims(x.dim):
#        handler_method = 'void %s::init_cache<%d>(SpaceElement<%d,0,%d,%d> &)' % (handler, k, x.dim, x.range, x.rank)
            handler_method = 'void %s::init_cache<%d>(ElementAccessor &)' % (handler, k)
            handler_methods.append(handler_method)
            handler_method = 'void %s::fill_cache<%d>(ElementAccessor &, const int)' % (handler, k)
            handler_methods.append(handler_method)
#--------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------
# SpaceElement used by PhysicalSpaceElement 
for space in inst.SubPhysSpaces + inst.PhysSpaces:
    x = space.spec
    for t in transformations:
        handler = 'SpaceElementHandler<%d,%d,%d,%d,%s>' %(x.dim,x.codim,x.range, x.rank,t)
        handlers.append(handler)
        for k in inst.sub_dims(x.dim):
#      handler_method = 'void %s::init_cache<%d>(SpaceElement<%d,%d,%d,%d> &)' % (handler, k, x.dim, x.codim, x.range, x.rank)
            handler_method = 'void %s::init_cache<%d>(ElementAccessor &)' % (handler, k)
            handler_methods.append(handler_method)
            handler_method = 'void %s::fill_cache<%d>(ElementAccessor &, const int)' % (handler, k)
            handler_methods.append(handler_method)

for space in inst.PhysSpaces:
    x = space.spec
    for t in transformations:
        handler = 'SpaceElementHandler<%d,%d,%d,%d,%s>' %(x.dim,x.codim,x.range, x.rank,t)
        handlers.append(handler)
        for k in inst.sub_dims(x.dim):
#        handler_method = 'void %s::init_cache<%d>(SpaceElement<%d,%d,%d,%d> &)' % (handler, k, x.dim, x.codim,x.range, x.rank)
            handler_method = 'void %s::init_cache<%d>(ElementAccessor &)' % (handler, k)
            handler_methods.append(handler_method)
            handler_method = 'void %s::fill_cache<%d>(ElementAccessor &, const int)' % (handler, k)
            handler_methods.append(handler_method)
#--------------------------------------------------------------------------------------



#---------------------------------------------------
for handler in unique(handlers):
    f.write('template class %s ;\n' %handler)

for handler_method in unique(handler_methods):
    f.write('template %s ;\n' % handler_method)



#---------------------------------------------------
f.write('IGA_NAMESPACE_CLOSE\n')
 
f.write('#ifdef SERIALIZATION\n')
id = 0 
for elem in unique(handlers):
    alias = 'SpaceElementHandlerAlias%d' %(id)
    f.write('using %s = iga::%s; \n' % (alias, elem))
    f.write('BOOST_CLASS_EXPORT_IMPLEMENT(%s) \n' %alias)
    f.write('template void %s::serialize(OArchive &, const unsigned int);\n' % alias)
    f.write('template void %s::serialize(IArchive &, const unsigned int);\n' % alias)
    id += 1 
f.write('#endif // SERIALIZATION\n')
     
f.write('IGA_NAMESPACE_OPEN\n')
#---------------------------------------------------


