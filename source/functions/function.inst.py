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

funcs = set()

for row in inst.all_function_dims:
    dims = '<%d,%d,%d,%d>' %(row.dim, row.codim, row.range, row.rank)
    func = 'Function%s' %(dims) 
    funcs.add(func)
    f.write('template class %s ;\n' %(func))



#---------------------------------------------------
f.write('#ifdef SERIALIZATION\n')
archives = ['OArchive','IArchive']

for func in funcs:
    for ar in archives:
        f.write('template void %s::serialize(%s&);\n' %(func,ar))
f.write('#endif // SERIALIZATION\n')
#---------------------------------------------------