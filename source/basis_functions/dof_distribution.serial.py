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

#types = ['BasisEndBehaviour','std::pair<Real,Real>']

#arr_t_ids = []
arr_dmas = []

dim = 0
range = 0
rank = 1
n_components = range ** rank
#arr_t_ids.append('SafeSTLArray<TensorIndex<%d>,%d>' %(dim,n_components))
arr_dmas.append('SafeSTLArray<DynamicMultiArray<int,%d>,%d>' %(dim,n_components))
arr_dmas.append('SafeSTLArray<CartesianProductArray<int,%d>,%d>' %(dim,n_components))

for x in inst.sub_ref_sp_dims:
    n_components = x.range ** x.rank
#    arr_t_ids.append('SafeSTLArray<TensorIndex<%d>,%d>' %(x.dim,n_components))
    arr_dmas.append('SafeSTLArray<DynamicMultiArray<int,%d>,%d>' %(x.dim,n_components))
    arr_dmas.append('SafeSTLArray<CartesianProductArray<int,%d>,%d>' %(x.dim,n_components))

for x in inst.ref_sp_dims:
    n_components = x.range ** x.rank
#    arr_t_ids.append('SafeSTLArray<TensorIndex<%d>,%d>' %(x.dim,n_components))
    arr_dmas.append('SafeSTLArray<DynamicMultiArray<int,%d>,%d>' %(x.dim,n_components))
    arr_dmas.append('SafeSTLArray<CartesianProductArray<int,%d>,%d>' %(x.dim,n_components))
            


         
#---------------------------------------------------
f.write('IGA_NAMESPACE_CLOSE\n')


f.write('#ifdef SERIALIZATION\n')
archives = ['OArchive','IArchive']

id = 0 
for arr_dma in unique(arr_dmas):
    alias = 'Array_DMA_int_Alias%d' %(id)
    f.write('using %s = %s;\n' % (alias, arr_dma.replace('DynamicMultiArray','iga::DynamicMultiArray')
                                                .replace('CartesianProductArray','iga::CartesianProductArray')
                                                .replace('SafeSTLArray','iga::SafeSTLArray')));
    f.write('CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(%s,cereal::specialization::member_serialize);\n' %alias);
    id += 1 


#id = 0 
#for arr_t_id in unique(arr_t_ids):
#    alias = 'Array_TIndex%d' %(id)
#    f.write('using %s = iga::%s;\n' % (alias, arr_t_id.replace('TensorIndex','iga::TensorIndex')
#                                                      .replace('SafeSTLArray','iga::SafeSTLArray')));
#    f.write('CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(%s,cereal::specialization::member_serialize);\n' %alias);
#    id += 1 


f.write('#endif // SERIALIZATION\n')

#   
f.write('IGA_NAMESPACE_OPEN\n')
#---------------------------------------------------
