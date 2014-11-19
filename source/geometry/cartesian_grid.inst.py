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

sub_dim_members = ['CartesianGrid<dim>::template BoundaryNormal<k> ' +
                   'CartesianGrid<dim>::get_boundary_normals<k>(const int s_id) const;',
                   'std::shared_ptr<CartesianGrid<k>> CartesianGrid<dim>::'+
                   'get_sub_grid<k>(const int sub_elem_id, InterGridMap<k> &elem_map) const;']
 
for dim in inst.sub_domain_dims:
    grid = 'CartesianGrid<%d>' %(dim)
    f.write('template class %s; \n' % (grid))
    for fun in sub_dim_members:
        k = dim
        s = fun.replace('k', '%d' % (k)).replace('dim', '%d' % (dim));
        f.write('template ' + s + '\n')  
    
for dim in inst.domain_dims:
    grid = 'CartesianGrid<%d>' %(dim)   
    f.write('template class %s; \n' % (grid))
    for fun in sub_dim_members:
        for k in inst.sub_dims(dim):
            s = fun.replace('k', '%d' % (k)).replace('dim', '%d' % (dim));
            f.write('template ' + s + '\n')
       
