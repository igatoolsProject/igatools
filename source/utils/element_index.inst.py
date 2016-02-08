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

# QA (pauletti, Mar 19, 2014):
from init_instantiation_data import *
data = Instantiation()
(f, inst) = (data.file_output, data.inst)

elem_ids = set()

for dim in inst.all_domain_dims:
    elem_id = 'ElementIndex<%d>' % (dim)
    elem_ids.add(elem_id)
    
    
for elem_id in elem_ids:
   f.write('template class %s; \n' %elem_id)
   f.write('template LogStream & operator<<(LogStream &,const %s &); \n' %(elem_id) )


#---------------------------------------------------
f.write('#ifdef SERIALIZATION\n')
archives = ['OArchive','IArchive']

for elem_id in elem_ids:
    for ar in archives:
        f.write('template void %s::serialize(%s&);\n' %(elem_id,ar))
f.write('#endif // SERIALIZATION\n')
#---------------------------------------------------
