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

bases = []





bases = ['Basis<0,0,0,1>']

for x in inst.all_phy_sp_dims:
    basis = 'Basis<%d,%d,%d,%d>' %(x.dim,x.codim,x.range,x.rank)
    bases.append(basis)

for x in inst.all_ref_sp_dims:
    basis = 'Basis<%d,0,%d,%d>' %(x.dim,x.range,x.rank)
    bases.append(basis)
    
for basis in unique(bases):
    f.write("template class %s ;\n" %(basis))


#---------------------------------------------------
f.write('#ifdef IGATOOLS_WITH_SERIALIZATION\n')
archives = ['OArchive','IArchive']

for basis in unique(bases):
    for ar in archives:
        f.write('template void %s::serialize(%s&);\n' %(basis,ar))
f.write('#endif // IGATOOLS_WITH_SERIALIZATION\n')
#---------------------------------------------------

