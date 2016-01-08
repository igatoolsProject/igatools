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
data = Instantiation()
(f, inst) = (data.file_output, data.inst)

funcs = set ()
for row in inst.all_function_dims:
    funcs.add ('functions::LinearFunction<%d,%d,%d>' 
         %(row.dim, row.codim, row.range) )

for row in inst.all_function_dims:
    funcs.add ('functions::ConstantFunction<%d,%d,%d,%d>' 
         % (row.dim, row.codim, row.range, row.rank))

for func in funcs:
  f.write('template class %s;\n' %(func))




#s = ('template class functions::SphereFunction<%d>;\n' %1 )
#f.write(s)
#s = ('template class functions::SphereFunction<%d>;\n' %2 )
#f.write(s)

#s = ('template class functions::CylindricalAnnulus<%d>;\n' %3)
#f.write(s)



#---------------------------------------------------
#f.write('#ifdef SERIALIZATION\n')
#archives = ['OArchive','IArchive']
#
#for func in unique(funcs):
#    for ar in archives:
#        f.write('template void %s::serialize(%s&);\n' %(func,ar))
#f.write('#endif // SERIALIZATION\n')
#---------------------------------------------------
