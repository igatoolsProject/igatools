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

types = ['BasisEndBehaviour']

spaces = []
arrays = [] 
operators = []

dim = 0
range = 0
rank = 1
n_components = range ** rank
space = 'BSpline<%d,%d,%d>' %(dim,range,rank)
spaces.append(space)
for t in types:
    arr = 'SafeSTLArray<%s,%d>' %(t,dim)
    arrays.append(arr)
    arr_arr = 'SafeSTLArray<%s,%d>' %(arr,n_components)
    arrays.append(arr_arr)
#arrays.append('SafeSTLArray<VecBernstOp,%d>' %(dim))
#arrays.append('SafeSTLArray<CartesianProductArray<BernsteinOperator,%d>,%d>' %(dim,n_components))


for x in inst.sub_ref_sp_dims + inst.ref_sp_dims:
    space = 'BSpline<%d,%d,%d>' %(x.dim, x.range, x.rank)
    spaces.append(space)
    n_components = x.range ** x.rank
    for t in types:
        arr = 'SafeSTLArray<%s,%d>' %(t,x.dim)
        arrays.append(arr)
        arr_arr = 'SafeSTLArray<%s,%d>' %(arr,n_components)
        arrays.append(arr_arr)
#    arrays.append('SafeSTLArray<VecBernstOp,%d>' %(x.dim))
#    arrays.append('SafeSTLArray<CartesianProductArray<BernsteinOperator,%d>,%d>' %(x.dim,n_components))


         
#---------------------------------------------------
f.write('IGA_NAMESPACE_CLOSE\n')


archives = ['OArchive','IArchive']


f.write('using VecBernstOp = iga::SafeSTLVector<iga::BernsteinOperator>;\n');
# Already done in spline_space.serial.py
# for ar in archives:
#    f.write('CEREAL_SPECIALIZE_FOR_ARCHIVE(%s,VecBernstOp,cereal::specialization::member_serialize)\n' %(ar));

id = 0 
for space in unique(spaces):
    sp_alias = 'BSplineAlias%d' %(id)
    f.write('using %s = iga::%s;\n' % (sp_alias, space));
    f.write('CEREAL_REGISTER_TYPE(%s)\n' %sp_alias);
    id += 1 


id = 0 
for arr in unique(arrays):
    alias = 'Array_Alias%d' %(id)
    f.write('using %s = %s;\n' % (alias, 
                                  arr.replace('BasisEndBehaviour','iga::BasisEndBehaviour')
                                     .replace('SafeSTLArray','iga::SafeSTLArray')
                                     .replace('CartesianProductArray','iga::CartesianProductArray')
                                     .replace('BernsteinOperator','iga::BernsteinOperator')
                                     .replace('Real','iga::Real')));
    for ar in archives:
        f.write('CEREAL_SPECIALIZE_FOR_ARCHIVE(%s,%s,cereal::specialization::member_serialize)\n' %(ar,alias));
    id += 1


#   
f.write('IGA_NAMESPACE_OPEN\n')
#---------------------------------------------------
