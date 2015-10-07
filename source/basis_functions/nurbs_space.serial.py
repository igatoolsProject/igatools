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

include_files = []
data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)


spaces = []


dim = 0
range = 0
rank = 1
space = 'NURBSSpace<%d,%d,%d>' %(dim,range,rank)
spaces.append(space)


for x in inst.sub_ref_sp_dims + inst.ref_sp_dims:
    space = 'NURBSSpace<%d,%d,%d>' %(x.dim, x.range, x.rank)
    spaces.append(space)


         
#---------------------------------------------------
f.write('IGA_NAMESPACE_CLOSE\n')


f.write('#ifdef SERIALIZATION\n')
#archives = ['OArchive','IArchive']


#f.write('using VecBernstOp = iga::SafeSTLVector<iga::BernsteinOperator>;\n');
#f.write('CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(VecBernstOp,cereal::specialization::member_serialize);\n');

id = 0 
for space in unique(spaces):
    sp_alias = 'NURBSSpaceAlias%d' %(id)
    f.write('using %s = iga::%s;\n' % (sp_alias, space));
    f.write('CEREAL_REGISTER_TYPE(%s);\n' %sp_alias);
    id += 1 


f.write('#endif // SERIALIZATION\n')

#   
f.write('IGA_NAMESPACE_OPEN\n')
#---------------------------------------------------
