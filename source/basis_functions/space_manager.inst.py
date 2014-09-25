#-+--------------------------------------------------------------------
# Igatools a general purpose Isogeometric analysis library.
# Copyright (C) 2012-2014  by the igatools authors (see authors.txt).
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

#spaces = []
#for space in inst.AllRefSpaces_v2:
#    spaces.append( '%s' %space.name)
#for space in inst.PhysSpaces_v2:
#    spaces.append( '%s' %space.name)

f.write( 'using SpacePtrVariant = Variant<\n')
for space in inst.AllRefSpaces_v2:
    f.write( 'std::shared_ptr<%s>,\n' %space.name)
for space in inst.PhysSpaces_v2:
    f.write( 'std::shared_ptr<%s>,\n' %space.name)
f.seek(-2,2);
f.write( '>;\n')

