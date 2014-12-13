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
 ['std::shared_ptr<typename class::template SubSpace<k>::MultiplicityTable> class::get_sub_space_mult<k>(const Index s_id) const;', 
  'typename class::template SubSpace<k>::DegreeTable class::get_sub_space_degree<k>(const Index s_id) const;',
  'typename class::template SubSpace<k>::EndBehaviourTable class::get_sub_space_end_b<k>(const Index s_id) const;']         

for x in inst.sub_ref_sp_dims:
    space = 'SplineSpace<%d, %d, %d>' %(x.dim, x.range, x.rank)
    f.write('template class %s ;\n' %space)
    for fun in sub_dim_members:
        k = x.dim
        s = fun.replace('class', space).replace('k', '%d' % (k));
        f.write('template ' + s + '\n')
    

for x in inst.ref_sp_dims:
    space = 'SplineSpace<%d, %d, %d>' %(x.dim, x.range, x.rank)
    f.write('template class %s ;\n' %space)
    for fun in sub_dim_members:
        for k in inst.sub_dims(x.dim):
            s = fun.replace('class', space).replace('k', '%d' % (k));
            f.write('template ' + s + '\n')
            
