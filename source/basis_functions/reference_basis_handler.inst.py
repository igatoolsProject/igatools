#-+--------------------------------------------------------------------
# Igatools a general purpose Isogeometric analysis library.
# Copyright (C) 2012-2016  by the igatools authors (see authors.txt).
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


handlers = ['ReferenceBasisHandler<0,0,1>']
handler_funcs = set()


for x in inst.sub_ref_sp_dims + inst.ref_sp_dims:
    handler = 'ReferenceBasisHandler<%d,%d,%d>' %(x.dim, x.range, x.rank)
    handlers.append(handler)
#    for k in range(0,x.dim+1):
#        func = 'int %s::get_num_points<%d>() const;' % (handler,k)
#        handler_funcs.add(func)

        
                
    
for handler in unique(handlers):
    f.write('template class %s;\n' %handler)

for func in handler_funcs:        
    f.write('template %s;\n' %func)        
   



