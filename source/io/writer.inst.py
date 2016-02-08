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

strings = []

writer_real_types = ['double']

for x in inst.mapping_dims:
    for writer_real_t in writer_real_types:
        writer = 'Writer<%d, %d, %s>' %(x.dim, x.codim, writer_real_t)
        strings.append('template class %s ;\n' % (writer))

for s in unique(strings): # Removing repeated entries.
    f.write(s)

