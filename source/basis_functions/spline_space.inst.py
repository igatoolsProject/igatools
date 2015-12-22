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
data = Instantiation()
f = data.file_output
inst = data.inst

sub_dim_members = \
 ['typename class::template SubSpace<k>::MultiplicityTable class::get_sub_space_mult<k>(const Index s_id) const', 
  'typename class::template SubSpace<k>::DegreeTable class::get_sub_space_degree<k>(const Index s_id) const',
  'typename class::template SubSpace<k>::PeriodicityTable class::get_sub_space_periodicity<k>(const Index s_id) const']




spaces = ['SplineSpace<0,0,1>']
templated_funcs = []

for x in inst.sub_ref_sp_dims:
    space = 'SplineSpace<%d,%d,%d>' %(x.dim, x.range, x.rank)
    spaces.append(space)
    for fun in sub_dim_members:
        for k in range(0,max(x.dim-1,0)+1):
#        k = x.dim
            s = fun.replace('class', space).replace('k', '%d' % (k));
            templated_funcs.append(s)
#            f.write('template ' + s + '\n')
    

for x in inst.ref_sp_dims:
    space = 'SplineSpace<%d,%d,%d>' %(x.dim, x.range, x.rank)
    spaces.append(space)
    for fun in sub_dim_members:
        for k in range(0,max(x.dim-1,0)+1):
#        for k in inst.sub_dims(x.dim):
            s = fun.replace('class', space).replace('k', '%d' % (k));
            templated_funcs.append(s)



#member_types = ['DegreeTable',  
#  'MultiplicityTable',
#  'BoundaryKnotsTable',
#  'KnotsTable',
#  'PeriodicityTable',
#  'EndBehaviourTable']


for space in unique(spaces):
    f.write('template class %s ;\n' %space)
    f.write('template class %s::template ComponentContainer<int>;\n' %(space))
    f.write('template class %s::template ComponentContainer<std::pair<Real,Real>>;\n' %(space))
    dim = '%s::dim' %(space)
    n_components = '%s::n_components' %(space)
    t0 = 'TensorSize<%s>' %(dim)
    f.write('template class %s::template ComponentContainer<%s>;\n' %(space,t0))
    t1 = 'CartesianProductArray<int,%s>' %(dim)
    f.write('template class %s::template ComponentContainer<%s>;\n' %(space,t1))
    f.write('template %s::template ComponentContainer<%s>::ComponentContainer(std::initializer_list<%s>,void *);\n' %(space,t1,t1))
    t2 = 'SafeSTLArray<bool,%s>' %(dim)
    f.write('template class %s::template ComponentContainer<%s>;\n' %(space,t2))
    f.write('template %s::template ComponentContainer<%s>::ComponentContainer(const %s &,void *);\n' %(space,t2,t2))
    t3 = 'SafeSTLArray<BasisEndBehaviour,%s>' %(dim)
    f.write('template class %s::template ComponentContainer<%s>;\n' %(space,t3))
    f.write('template %s::template ComponentContainer<%s>::ComponentContainer(const %s &,void *);\n' %(space,t3,t3))
    f.write('template %s::template ComponentContainer<%s>::ComponentContainer(const SafeSTLArray<int,%s> &,const %s &,void *);\n' %(space,t3,n_components,t3))
    f.write('template %s::template ComponentContainer<%s>::ComponentContainer(std::initializer_list<%s>,void *);\n' %(space,t3,t3))
    t4 = 'SafeSTLArray<std::pair<Real,Real>,%s>' %(dim)
    f.write('template class %s::template ComponentContainer<%s>;\n' %(space,t4))
    t5 = 'SafeSTLArray<SafeSTLVector<Real>,%s>' %(dim)
    f.write('template class %s::template ComponentContainer<%s>;\n' %(space,t5))
    t6 = 'SafeSTLArray<CartesianProductArray<double,2>,%s>' % (dim) 
    f.write('template class %s::template ComponentContainer<%s>;\n' %(space,t6))
    f.write('template %s::template ComponentContainer<%s>::ComponentContainer(std::initializer_list<%s>,void *);\n' %(space,t6,t6))
    t7 = 'unique_ptr<const TensorProductFunctionEvaluator<%s>>'%(dim)
    f.write('template class %s::template ComponentContainer<%s>;\n' %(space,t7))
    t8 = 'SafeSTLArray<BasisValues1d,%s>' % (dim)
    f.write('template class %s::template ComponentContainer<%s>;\n' %(space,t8))
    t9 = 'CartesianProductArray<BernsteinOperator,%s>' % (dim)
    f.write('template class %s::template ComponentContainer<%s>;\n' %(space,t9))
    t10 = 'std::shared_ptr<CartesianProductIndexer<%s>>' % (dim)
    f.write('template class %s::template ComponentContainer<%s>;\n' %(space,t10))
    t11 = 'SafeSTLArray<const BernsteinOperator *,%s>' % (dim)
    f.write('template class %s::template ComponentContainer<%s>;\n' %(space,t11))
    t12 = 'TensorIndex<%s>' %(dim)
    f.write('template class %s::template ComponentContainer<%s>;\n' %(space,t12))
    f.write('template %s::template ComponentContainer<%s>::ComponentContainer(const %s &,void *);\n' %(space,t12,t12))
    f.write('template %s::template ComponentContainer<%s>::ComponentContainer(std::initializer_list<%s>,void *);\n' %(space,t12,t12))

for func in unique(templated_funcs):
    f.write('template %s;\n' %func)




#---------------------------------------------------
f.write('#ifdef SERIALIZATION\n')

archives = ['OArchive','IArchive']

for space in unique(spaces):
    for ar in archives:
        f.write('template void %s::serialize(%s&);\n' %(space,ar))
f.write('#endif // SERIALIZATION\n')
#---------------------------------------------------

