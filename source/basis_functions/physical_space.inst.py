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

# QA (pauletti, Mar 19, 2014):
from init_instantiation_data import *

include_files = ['basis_functions/bspline_space.h',
                 'basis_functions/nurbs_space.h',
                 'geometry/push_forward.h',
                 'geometry/cartesian_grid_element.h',
                 'geometry/mapping_element.h',
                 'geometry/push_forward_element.h',
                 'basis_functions/bspline_element.h',
                 'basis_functions/physical_space_element.h']

data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)


sub_dim_members = \
 ['std::shared_ptr<typename class::SubSpace<k> > ' +
   'class::get_sub_space(const int s_id, InterSpaceMap<k> &dof_map, ' +
   'std::shared_ptr<CartesianGrid<k>> sub_grid, ' + 
   'InterGridMap &elem_map) const;']


spaces = ['PhysicalSpace<0,0,1,0, Transformation::h_grad>']


for sp in inst.SubPhysSpaces:
    x = sp.spec
#    f.write( 'template class %s;\n' %sp.name)
    spaces.append(sp.name)
    for fun in sub_dim_members:
        k = x.dim
        s = fun.replace('class', sp.name).replace('k', '%d' % (k));
        f.write('template ' + s + '\n')


for sp in inst.PhysSpaces:
    x = sp.spec
#    f.write( 'template class %s;\n' %sp.name)
    spaces.append(sp.name)
    for fun in sub_dim_members:
        for k in inst.sub_dims(x.dim):
            s = fun.replace('class', sp.name).replace('k', '%d' % (k));
            f.write('template ' + s + '\n')



for space in unique(spaces):
    f.write( 'template class %s;\n' %space)


#---------------------------------------------------
f.write('IGA_NAMESPACE_CLOSE\n')
 
f.write('#ifdef SERIALIZATION\n')
id = 0 
for space in unique(spaces):
    alias = 'PhysicalSpaceAlias%d' %(id)
    f.write('using %s = iga::%s; \n' % (alias, space.replace('Transformation','iga::Transformation')))
    f.write('BOOST_CLASS_EXPORT_IMPLEMENT(%s) \n' %alias)
    f.write('template void %s::serialize(OArchive &, const unsigned int);\n' % alias)
    f.write('template void %s::serialize(IArchive &, const unsigned int);\n' % alias)
    id += 1 
f.write('#endif // SERIALIZATION\n')
     
f.write('IGA_NAMESPACE_OPEN\n')
#---------------------------------------------------

