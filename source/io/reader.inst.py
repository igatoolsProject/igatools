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

# # QA (pauletti, Mar 19, 2014):
from init_instantiation_data import *
include_files = []
#include_files = ['basis_functions/bspline_space.h',
#                 'basis_functions/nurbs_space.h']
# data = Instantiation(include_files)
data = Instantiation()
(f, inst) = (data.file_output, data.inst)
# 
# 
# grids = ['std::shared_ptr<CartesianGrid<%d>>' %(dim) for dim in inst.user_domain_dims]
# ref_spaces = ['BSplineSpace']#todo use general variable, 'NURBSSpace')
# for grid in grids:
#     f.write('template %s get_cartesian_grid_from_xml(const boost::property_tree::ptree &);\n' %grid)
# 
#maps = ['std::shared_ptr<MapFunction<%d,%d>>' %(row.dim,row.dim+row.codim) for row in inst.mapping_dims]
#for map in maps:
#    f.write('template %s get_mapping_from_xml(const boost::property_tree::ptree &);\n' %map)
#    f.write('template %s get_mapping_from_file(const std::string &);\n' %map)
for x in inst.mapping_dims:
    map = 'std::shared_ptr<MapFunction<%d,%d>>' %(x.dim, x.dim + x.codim)
    f.write('template %s get_mapping_from_file<%d,%d>(const std::string &);\n' %(map,x.dim,x.codim))

# 
# bsp_spaces = ['std::shared_ptr<BSplineSpace<%d,%d,%d>>' %(x.dim, x.range, x.rank)  
#           for x in inst.all_ref_sp_dims ] 
# # nrb_spaces = ['std::shared_ptr<NURBSSpace<%d,%d,%d>>' %(x.dim, x.range, x.rank)  
# #           for x in inst.all_ref_sp_dims ]
# 
# for sp in bsp_spaces:
#     f.write('template %s get_bspline_space_from_xml(const boost::property_tree::ptree &);\n' % (sp))
# # for sp in nrb_spaces:
# #     f.write('template %s get_nurbs_space_from_xml(const boost::property_tree::ptree &);\n' % (sp))
