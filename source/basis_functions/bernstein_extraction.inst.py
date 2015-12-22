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

# QA (pauletti, Jun 6, 2014):
from init_instantiation_data import *

data = Instantiation()
(f, inst) = (data.file_output, data.inst)

classes = ['BernsteinExtraction<%d,%d,%d>' %(x.dim, x.range, x.rank)  
          for x in inst.all_ref_sp_dims]

classes.append('BernsteinExtraction<0,0,1>')

for c in classes:
   f.write('template class %s ;\n' %c)


#---------------------------------------------------
f.write('#ifdef SERIALIZATION\n')

archives = ['OArchive','IArchive']

classes.append('BernsteinOperator')
for c in unique(classes):
    for ar in archives:
        f.write('template void %s::serialize(%s&);\n' %(c,ar))
f.write('#endif // SERIALIZATION\n')
#---------------------------------------------------
