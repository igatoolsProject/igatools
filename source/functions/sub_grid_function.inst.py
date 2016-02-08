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


classes = set()


for x in inst.sub_mapping_dims + inst.mapping_dims:
    for sdim in range(0,x.dim):
        cl = 'SubGridFunction<%d,%d,%d>' %(sdim,x.dim,x.space_dim)
        classes.add(cl)

#for x in inst.mapping_dims:
#    for sdim in range(0,x.dim):
#        cl = 'SubGridFunction<%d,%d,%d>' %(sdim,x.dim,x.space_dim)
#        classes.append(cl)

    #the next classes are needed by NURBS
        cl = 'SubGridFunction<%d,%d,1>' %(sdim,x.dim)
        classes.add(cl)

 


for cl in classes:
    f.write('template class %s;\n' %(cl))

