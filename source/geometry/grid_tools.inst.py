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

include_files = ['geometry/grid_element.h']
data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)

grids = ['Grid<%d>' % (dim) for dim in inst.all_domain_dims]
for grid in grids:
    map = 'SafeSTLMap<typename %s::IndexType,typename %s::IndexType>' %(grid,grid)
    f.write('template std::shared_ptr<%s> grid_tools::build_grid_union('
            'const %s &,const %s &,%s &,%s &); \n' % (grid,grid,grid, map,map))
    
for dim in inst.all_domain_dims:
    grid = 'Grid<%d>' % (dim)
    index = 'typename %s::IndexType' % (grid)
    f.write('template iga::SafeSTLMap<%s, %s>' % (index, index) + 
            ' grid_tools::build_map_elements_id_between_grids(' +
            'const %s &,const %s &); \n' % (grid, grid))
#    f.write('template SafeSTLArray<SafeSTLVector<Index>,%d>' % (dim)+
#            ' grid_tools::build_map_intervals_id_between_grids('
#            'const %s &,const %s &); \n' % (grid,grid))
#    f.write('template SafeSTLArray<Index,%d>' % (dim)+
#            ' grid_tools::get_max_num_fine_intervals_in_coarse_interval('
#            'const %s &,const %s &); \n' % (grid,grid))
