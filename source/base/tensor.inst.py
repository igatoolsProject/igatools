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

# todo: this class is inline for completeness many things would need to be inst
from init_instantiation_data import *
data = Instantiation()
(f, inst) = (data.file_output, data.inst)

for row in inst.derivatives + inst.values:
    f.write('template class %s ;\n' %row)
    f.write('template EnableIf<%s::is_tensor,%s> operator+ < %s > (const %s &A, const %s &B) noexcept ;\n'
            % (row, row, row, row, row))
    f.write('template EnableIf<%s::is_tensor,%s> operator- < %s > (const %s &A, const %s &B) noexcept ;\n'
            % (row, row, row, row, row))
