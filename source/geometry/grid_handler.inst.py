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

include_files = ['geometry/grid_element.h']
data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)




handlers = set()
handler_funcs = set()

for dim in inst.sub_domain_dims + inst.domain_dims:
    handler = 'GridHandler<%d>' %(dim)
    handlers.add(handler)
    for k in range(0,dim+1):
        func = 'void %s::SetFlagsDispatcher::operator()(const Topology<%d> &)' % (handler,k)
        handler_funcs.add(func)
        func = 'void %s::InitCacheDispatcher::operator()(const std::shared_ptr<const Quadrature<%d>> &)' % (handler,k)
        handler_funcs.add(func)
        func = 'void %s::FillCacheDispatcher::operator()(const Topology<%d> &)' % (handler,k)
        handler_funcs.add(func)
        func = 'void %s::set_flags<%d>(const Flags &)' % (handler,k)
        handler_funcs.add(func)
        func = 'void %s::init_cache<%d>(ElementAccessor &,const std::shared_ptr<const Quadrature<%d>> &) const' % (handler,k,k)
        handler_funcs.add(func)
        func = 'void %s::init_cache<%d>(ElementIterator &,const std::shared_ptr<const Quadrature<%d>> &) const' % (handler,k,k)
        handler_funcs.add(func)
        func = 'void %s::fill_cache<%d>(ElementAccessor &,const int) const' % (handler,k)
        handler_funcs.add(func)
        func = 'void %s::fill_cache<%d>(ElementIterator &,const int) const' % (handler,k)
        handler_funcs.add(func)



for handler in handlers:
    f.write('template class %s;\n' %(handler))
    
for func in handler_funcs:
    f.write('template %s;\n' %(func))
