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
include_files = ['utils/tensor_index.h',
                 'utils/element_index.h',
                 'basis_functions/bernstein_extraction.h',
                 'basis_functions/values1d_const_view.h',
                 'basis_functions/nurbs.h',
                 'utils/concatenated_iterator.h',
                 'io/writer.h']
data = Instantiation(include_files)
f = data.file_output
inst = data.inst

classes = []



classes.append('int');
classes.append('Real');
classes.append('SafeSTLVector<Real>');
classes.append('ContainerView<SafeSTLVector<int>>');
classes.append('ConstContainerView<SafeSTLVector<int>>');
classes.append('BernsteinOperator');
classes.append('const BernsteinOperator *');
classes.append('SafeSTLVector<BernsteinOperator>');
classes.append('const SafeSTLVector<BernsteinOperator> *');
for dim in inst.all_domain_dims:
    classes.append('TensorIndex<%d>' %(dim))
    classes.append('ElementIndex<%d>' %(dim))
    classes.append('SafeSTLArray<Real,%d>' %(dim))
    classes.append('SafeSTLArray<SafeSTLArray<Real,%d>,%d>' %(dim,dim))
    classes.append('SafeSTLArray<BernsteinOperator,%d>' %(dim))
    classes.append('SafeSTLArray<const BernsteinOperator *,%d>' %(dim))
    classes.append('SafeSTLArray<SafeSTLVector<BernsteinOperator>,%d>' %(dim));
    classes.append('SafeSTLArray<SafeSTLVector<const BernsteinOperator *>,%d>' %(dim));
    classes.append('SafeSTLArray<const SafeSTLVector<BernsteinOperator> *,%d>' %(dim));

vec_int = 'SafeSTLVector<int>'
t = 'ConstView<ConcatenatedIterator<ContainerView<%s>>,ConcatenatedConstIterator<ContainerView<%s>,ConstContainerView<%s>>>' %(vec_int,vec_int,vec_int)
classes.append('%s' % (t))


space = 'SplineSpace<0,0,1>'
classes.append('typename %s::template ComponentContainer<SafeSTLArray<BasisValues1d,0>>' % (space));
for x in inst.sub_ref_sp_dims + inst.ref_sp_dims:
    space = 'SplineSpace<%d,%d,%d>' %(x.dim, x.range, x.rank)
    classes.append('typename %s::template ComponentContainer<SafeSTLArray<BasisValues1d,%d>>' %(space,x.dim));
    

for x in inst.sub_ref_sp_dims + inst.ref_sp_dims:
    bspline = 'std::shared_ptr<BSpline<%d,%d,%d>>' %(x.dim, x.range, x.rank)
    classes.append('%s' %(bspline))
    nurbs = 'std::shared_ptr<NURBS<%d,%d,%d>>' %(x.dim, x.range, x.rank)
    classes.append('%s' %(nurbs))
    vec = 'SafeSTLVector<int>'
    classes.append('ConstView<MultiArrayIterator<MultiArray<%s,%d>>,MultiArrayConstIterator<MultiArray<%s,%d>>>' %(vec,x.dim,vec,x.dim))




classes.append('ConstView<std::vector<int>::iterator,std::vector<int>::const_iterator>')
classes.append('ConstView<std::vector<int>::iterator,std::vector<int>::iterator>')
classes.append('NonConstView<std::vector<int>::iterator,std::vector<int>::const_iterator>')

# Needed for writer --------------------------------------------------------------#

for dims in inst.dims_list:
    classes.append('SafeSTLArray<int,%d>' % (pow(2, dims.dim)))
    classes.append('SafeSTLArray<int,%d>' % (dims.dim))
    classes.append('SafeSTLVector<SafeSTLArray<int,%d>>' % (pow(2, dims.dim)))
    classes.append('SafeSTLVector<SafeSTLArray<int,%d>>' % (dims.dim))
    
classes.append('SafeSTLVector<int>')

for x in inst.mapping_dims:
    writer_types = ['double']
    for writer_t in writer_types:
        classes.append('Writer<%d,%d,%s>::PointData' % (x.dim, x.codim, writer_t))
        for writer_t_2 in writer_types + ['int']:
          classes.append('Writer<%d,%d,%s>::CellData<%s>' % (x.dim, x.codim, writer_t, writer_t_2))

classes.append('SafeSTLVector<SafeSTLArray<double,3>>')
# classes.append('float')
# classes.append('SafeSTLArray<float,3>')
# classes.append('SafeSTLVector<SafeSTLArray<float,3>>')
classes.append('SafeSTLVector<SafeSTLVector<double>>')
#----------------------------------------------------------------------------------#

classes.append('std::string')



for val in inst.values:
    classes.append('%s' % val)

for der in inst.derivatives:
    classes.append('%s' % der)

for div in inst.divs:
    classes.append('%s' % div)


for dims in inst.dims_list:
    value ='Tensor<%d, %d, tensor::contravariant, Tdouble>' % (dims.range,dims.rank)
    for k in range(0,dims.range+1):
        classes.append('SafeSTLArray<%s,%d>' % (value,k))



vectors = []
for c in unique(classes):
    v = 'SafeSTLVector<%s>' %(c)
    vectors.append('%s' %(v))
    f.write('template class %s;\n' % (v))






#---------------------------------------------------
#f.write('IGA_NAMESPACE_CLOSE\n')
#
#f.write('#ifdef SERIALIZATION\n')
#
#archives = ['OArchive','IArchive']
#
#id = 0 
#for v in unique(vectors):
#    alias = 'SafeSTLVectorAlias%d' %(id)
#    f.write('using %s = %s;\n' %(alias,v))
#    for ar in archives:
#        f.write('template void %s::serialize(%s&);\n' %(alias,ar))
#    id += 1 
#f.write('#endif // SERIALIZATION\n')
#    
#f.write('IGA_NAMESPACE_OPEN\n')
#---------------------------------------------------






