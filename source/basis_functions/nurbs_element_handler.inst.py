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

include_files = ['basis_functions/nurbs_space.h',
                 '../../source/basis_functions/nurbs_element.cpp',
                 '../../source/geometry/cartesian_grid_iterator.cpp'
                 ]
data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)


sub_dim_members = []

        
        
        
handlers = ['NURBSElementHandler<0,0,1>']
elements = ['NURBSElement<0,0,1>']
handler_templated_funcs = []


for x in inst.sub_ref_sp_dims:
    elem = 'NURBSElement<%d,%d,%d>' %(x.dim, x.range, x.rank)
    elements.append(elem)
    handler = 'NURBSElementHandler<%d,%d,%d>' %(x.dim, x.range, x.rank)
    handlers.append(handler)
    for fun in sub_dim_members:
        k = x.dim
        s = fun.replace('elhandler', handler).replace('k', '%d' % (k));
        handler_templated_funcs.append(s)


for x in inst.ref_sp_dims:
    elem = 'NURBSElement<%d,%d,%d>' %(x.dim, x.range, x.rank)
    elements.append(elem)
    handler = 'NURBSElementHandler<%d,%d,%d>' %(x.dim, x.range, x.rank)
    handlers.append(handler)
    for fun in sub_dim_members:
        for k in inst.sub_dims(x.dim):
            s = fun.replace('elhandler', handler).replace('k', '%d' % (k));
            handler_templated_funcs.append(s)
        
        
        
for elem in unique(elements):
    f.write('template class %s; \n' %elem)
    
for handler in unique(handlers):
    f.write('template class %s; \n' %handler)

for func in unique(handler_templated_funcs):        
    f.write('template %s; \n' %func)        
   
