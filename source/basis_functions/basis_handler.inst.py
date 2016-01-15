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




sub_dim_members = \
[]

handlers = set()
handler_methods = set()

handler = 'BasisHandler<0,0,0,1>'
handlers.add(handler)

handler_method = 'void %s::set_flags<0>(const typename basis_element::Flags &flag)' % (handler)
handler_methods.add(handler_method)
handler_method = 'void %s::init_cache<0>(ElementAccessor &elem,const std::shared_ptr<const Quadrature<0>> &quad) const' % (handler)
handler_methods.add(handler_method)
handler_method = 'void %s::init_cache<0>(ElementIterator &elem,const std::shared_ptr<const Quadrature<0>> &quad) const' % (handler)
handler_methods.add(handler_method)
handler_method = 'void %s::fill_cache<0>(ElementAccessor &elem, const int s_id) const' % (handler)
handler_methods.add(handler_method)
handler_method = 'void %s::fill_cache<0>(ElementIterator &elem, const int s_id) const' % (handler)
handler_methods.add(handler_method)


#--------------------------------------------------------------------------------------
# BasisElement used by ReferenceBasisElement 
for x in inst.sub_ref_sp_dims + inst.ref_sp_dims:
    handler = 'BasisHandler<%d,0,%d,%d>' %(x.dim, x.range, x.rank)
    handlers.add(handler)
    for k in range(0,x.dim+1):
        handler_method = 'void %s::set_flags<%d>(const typename basis_element::Flags &flag)' % (handler,k)
        handler_methods.add(handler_method)

        handler_method = 'void %s::init_cache<%d>(ElementAccessor &elem,const std::shared_ptr<const Quadrature<%d>> &quad) const' % (handler,k,k)
        handler_methods.add(handler_method)
        handler_method = 'void %s::init_cache<%d>(ElementIterator &elem,const std::shared_ptr<const Quadrature<%d>> &quad) const' % (handler,k,k)
        handler_methods.add(handler_method)
        
        handler_method = 'void %s::fill_cache<%d>(ElementAccessor &elem, const int s_id) const' % (handler,k)
        handler_methods.add(handler_method)
        handler_method = 'void %s::fill_cache<%d>(ElementIterator &elem, const int s_id) const' % (handler,k)
        handler_methods.add(handler_method)
#--------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------
# BasisElement used by PhysicalBasisElement 
for basis in inst.SubPhysBases + inst.PhysBases:
    x = basis.spec
    handler = 'BasisHandler<%d,%d,%d,%d>' %(x.dim,x.codim,x.range, x.rank)
    handlers.add(handler)
    for k in range(0,x.dim+1):
        handler_method = 'void %s::set_flags<%d>(const typename basis_element::Flags &flag)' % (handler,k)
        handler_methods.add(handler_method)

        handler_method = 'void %s::init_cache<%d>(ElementAccessor &elem,const std::shared_ptr<const Quadrature<%d>> &quad) const' % (handler,k,k)
        handler_methods.add(handler_method)
        handler_method = 'void %s::init_cache<%d>(ElementIterator &elem,const std::shared_ptr<const Quadrature<%d>> &quad) const' % (handler,k,k)
        handler_methods.add(handler_method)
        
        handler_method = 'void %s::fill_cache<%d>(ElementAccessor &elem, const int s_id) const' % (handler,k)
        handler_methods.add(handler_method)
        handler_method = 'void %s::fill_cache<%d>(ElementIterator &elem, const int s_id) const' % (handler,k)
        handler_methods.add(handler_method)
#--------------------------------------------------------------------------------------



#---------------------------------------------------
for handler in handlers:
    f.write('template class %s;\n' %handler)

for handler_method in handler_methods:
    f.write('template %s;\n' % handler_method)



