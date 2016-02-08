#-+--------------------------------------------------------------------
# Igatools a general purpose Isogeometric analysis library.
# Copyright (C) 2012-2016  by the igatools authors (see authors.txt).
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

include_files = ['geometry/grid_element.h',
                 'basis_functions/bspline_element.h']
data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)

sub_dim_members = \
 ['std::shared_ptr<const typename class::template BSpline<sdim,range,rank>> ' + 
  'class::get_sub_bspline_basis<sdim>(const int s_id, ' + 
  'InterBasisMap<sdim> &dof_map, ' +
  'const std::shared_ptr<const Grid<sdim>> &sub_grid) const'
#  ,
#  'std::shared_ptr<typename class::template SubBasis<k>> ' + 
#  'class::get_sub_basis<k>(const int sub_elem_id, ' + 
#  'InterBasisMap<k> &dof_map, SubGridMap<k> &elem_map) const;'
  ]
 


bases = ['BSpline<0,0,1>']
templated_funcs = []

for x in inst.sub_ref_sp_dims:
    basis = 'BSpline<%d,%d,%d>' %(x.dim, x.range, x.rank)
    bases.append(basis)
    for fun in sub_dim_members:
        for k in range(0,max(x.dim-1,0)+1):
#        k = max(x.dim-1,0)
            s = fun.replace('class', basis).replace('sdim','%d' % (k)).replace('range','%d' % (x.range)).replace('rank','%d' % (x.rank));
            templated_funcs.append(s)


for x in inst.ref_sp_dims:
    basis = 'BSpline<%d,%d,%d>' %(x.dim, x.range, x.rank)
    bases.append(basis)
    for fun in sub_dim_members:
        for k in range(0,max(x.dim-1,0)+1):
#        for k in inst.sub_dims(x.dim):
#            if (k < x.dim):
                s = fun.replace('class', basis).replace('sdim','%d' % (k)).replace('range','%d' % (x.range)).replace('rank','%d' % (x.rank));
                templated_funcs.append(s)
            
            
            
for basis in unique(bases):
    f.write('template class %s ;\n' %basis)

for func in unique(templated_funcs):
    f.write('template %s ;\n' %func)



#---------------------------------------------------
f.write('#ifdef SERIALIZATION\n')

archives = ['OArchive','IArchive']

for basis in unique(bases):
    for ar in archives:
        f.write('template void %s::serialize(%s&);\n' %(basis,ar))
f.write('#endif // SERIALIZATION\n')
#---------------------------------------------------


    
