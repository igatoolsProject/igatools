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

include_files = ['../../source/geometry/cartesian_grid_iterator.cpp',
		         '../../source/base/function.cpp']

data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)

accessors = ['FunctionElement<%d, %d, %d>' %(row.dim, row.range, row.rank) 
			for row in inst.all_function_dims]
for row in accessors:
	f.write('template class %s; \n' %(row))
	for it in inst.iterators:
		iterator = it.replace('Accessor','%s' % (row) )
		f.write('template class %s; \n' %iterator)
