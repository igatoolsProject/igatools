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


#---------------------------------------------------
f.write('IGA_NAMESPACE_CLOSE\n')

archives = ['OArchive','IArchive']

id = 0 
for dim in unique(inst.all_domain_dims):
    alias = 'TensorSizeAlias%d' %(id)
    f.write('using %s = iga::TensorSize<%d>;\n' % (alias,dim))
    for ar in archives:
        f.write('CEREAL_SPECIALIZE_FOR_ARCHIVE(%s,%s,cereal::specialization::member_serialize)\n' % (ar,alias));

    id += 1 
#   
f.write('IGA_NAMESPACE_OPEN\n')
#---------------------------------------------------
