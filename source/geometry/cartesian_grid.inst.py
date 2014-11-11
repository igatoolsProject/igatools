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

include_files = ['geometry/cartesian_grid_element.h']
data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)

for dim in inst.all_domain_dims:
    grid = 'CartesianGrid<%d>' %(dim)
    f.write('template class %s; \n' % (grid))
    
for dim in inst.domain_dims:
    grid = 'CartesianGrid<%d>' %(dim)   
    for k in inst.sub_dims(dim):
        f.write('template std::shared_ptr<CartesianGrid<%d>> %s::get_sub_grid(const int sub_elem_id, InterGridMap<%d> &elem_map) const; \n' % (k,grid,k))
        f.write('template std::array<Points<%d>, %d> ' %(dim,dim-k) +
                '%s::get_normal_space<%d>(const int j) const; \n' % (grid,dim-k))