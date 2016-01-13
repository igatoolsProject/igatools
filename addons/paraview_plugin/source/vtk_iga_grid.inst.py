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

include_files = []

data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)


domains = set()

for x in inst.sub_mapping_dims:
    if x.dim > 0 and x.dim < 4 and (x.dim + x.codim) > 0 and (x.dim + x.codim) < 4:
      domain = 'Domain<%d,%d>' %(x.dim, x.codim)
      domains.add(domain)

for x in inst.mapping_dims:
    if x.dim > 0 and x.dim < 4 and (x.dim + x.codim) > 0 and (x.dim + x.codim) < 4:
      domain = 'Domain<%d,%d>' %(x.dim, x.codim)
      domains.add(domain)

f.write("namespace paraview_plugin {\n")
 
for domain in domains:
    f.write('template class VtkIgaGrid<%s>;\n' %(domain))

f.write("};")