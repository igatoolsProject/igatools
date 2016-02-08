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

include_files = []

data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)

         
#---------------------------------------------------
f.write('class BernsteinOperator;\n')

f.write('IGA_NAMESPACE_CLOSE\n')


archives = ['OArchive','IArchive']

aliases = set()

arrays = set()




#-----------------------------------------------
for dim in unique(inst.sub_domain_dims + inst.domain_dims):
    arrays.add('iga::SafeSTLArray<int,%d>' %(dim))

    arrays.add('iga::SafeSTLArray<int,%d>' %(2*dim)) #this is for Grid::boundary_id_
    
    arrays.add('iga::SafeSTLArray<iga::SafeSTLVector<iga::Real>,%d>' % (dim))
#-----------------------------------------------


#-----------------------------------------------
arrays.add('iga::SafeSTLArray<iga::Real,2>')

arrays.add('iga::SafeSTLArray<iga::SafeSTLArray<iga::Real,2>,0>')
arrays.add('iga::SafeSTLArray<iga::SafeSTLArray<iga::SafeSTLArray<iga::Real,2>,0>,0>')
arrays.add('iga::SafeSTLArray<iga::SafeSTLArray<iga::SafeSTLVector<int>,0>,0>')
arrays.add('iga::SafeSTLArray<iga::SafeSTLArray<iga::SafeSTLVector<iga::BernsteinOperator>,0>,0>')
for x in inst.all_ref_sp_dims:
#    space = 'SplineSpace<%d,%d,%d>' %(x.dim, x.range, x.rank)
    n_components = x.range**x.rank
    arrays.add('iga::SafeSTLArray<iga::SafeSTLArray<iga::Real,2>,%d>' %(x.dim))
    arrays.add('iga::SafeSTLArray<iga::TensorSize<%d>,%d>' % (x.dim,n_components))
    arrays.add('iga::SafeSTLArray<iga::SafeSTLArray<iga::SafeSTLArray<iga::Real,2>,%d>,%d>' %(x.dim,n_components))
    arrays.add('iga::SafeSTLArray<iga::SafeSTLArray<iga::SafeSTLVector<int>,%d>,%d>' %(x.dim,n_components))
    arrays.add('iga::SafeSTLArray<iga::SafeSTLArray<iga::SafeSTLVector<iga::BernsteinOperator>,%d>,%d>' %(x.dim,n_components))
#-----------------------------------------------


#-----------------------------------------------
id = 0
for array in arrays:
    alias = 'SafeSTLArrayAlias%d' %(id)
    aliases.add(alias)
    f.write('using %s = %s;\n' %(alias,array))
    id += 1 
#-----------------------------------------------
    
    
    
    
for alias in aliases:
    for ar in archives:
        f.write('CEREAL_SPECIALIZE_FOR_ARCHIVE(%s,%s,cereal::specialization::member_serialize)\n' % (ar,alias));
#   
f.write('IGA_NAMESPACE_OPEN\n')
#---------------------------------------------------
