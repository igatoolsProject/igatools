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

# QA (pauletti, Mar 19, 2014):
from init_instantiation_data import *
data = Instantiation()
f = data.file_output
inst = data.inst

sub_dim_members = \
 ['std::shared_ptr<typename class::template SubRefSpace<k>> ' + 
  'class::get_ref_sub_space<k>(const int sub_elem_id, ' + 
  'InterSpaceMap<k> &dof_map, ' + 
  'std::shared_ptr<CartesianGrid<k>> sub_grid) const;',
  'std::shared_ptr<typename class::template SubSpace<k>> ' + 
  'class::get_sub_space<k>(const int sub_elem_id, ' + 
  'InterSpaceMap<k> &dof_map, ' + 
  'std::shared_ptr<CartesianGrid<k>> sub_grid, ' +
  'SubGridMap<k> &elem_map) const;']

spaces = ['ReferenceSpace<0,0,1>']
templated_funcs = []

for x in inst.sub_ref_sp_dims:
    space = 'ReferenceSpace<%d, %d, %d>' %(x.dim, x.range, x.rank)
    spaces.append(space)
    for fun in sub_dim_members:
        k = x.dim
        s = fun.replace('class', space).replace('k', '%d' % (k));
        templated_funcs.append(s)

for x in inst.ref_sp_dims:
    space = 'ReferenceSpace<%d, %d, %d>' %(x.dim, x.range, x.rank)
    spaces.append(space)
    for fun in sub_dim_members:
        for k in inst.sub_dims(x.dim):
            s = fun.replace('class', space).replace('k', '%d' % (k));
            templated_funcs.append(s)


for space in unique(spaces):
    f.write('template class %s ;\n' %space)

for func in unique(templated_funcs):
    f.write('template %s ;\n' %func)


#---------------------------------------------------
f.write('IGA_NAMESPACE_CLOSE\n')
 
f.write('#ifdef SERIALIZATION\n')
id = 0 
for space in unique(spaces):
    alias = 'ReferenceSpaceAlias%d' %(id)
    f.write('using %s = iga::%s; \n' % (alias, space))
    f.write('BOOST_CLASS_EXPORT_IMPLEMENT(%s) \n' %alias)
    f.write('template void %s::serialize(OArchive &, const unsigned int);\n' % alias)
    f.write('template void %s::serialize(IArchive &, const unsigned int);\n' % alias)
    id += 1 
f.write('#endif // SERIALIZATION\n')
    
f.write('IGA_NAMESPACE_OPEN\n')
#---------------------------------------------------

