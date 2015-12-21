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

include_files = ['base/value_types.h',
                 'utils/safe_stl_vector.h',
                 'basis_functions/bernstein_extraction.h',
                 'basis_functions/spline_space.h',
                 'basis_functions/values1d_const_view.h',
                 'geometry/grid_function_element.h']

data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)

array_list = []

types = ['int',
         'Real',
         'bool',
         'SafeSTLVector<int>',
         'SafeSTLVector<Real>',
         'std::shared_ptr<SafeSTLVector<Real>>',
         'BernsteinOperator',
         'const BernsteinOperator *',         
         'SafeSTLVector<BernsteinOperator>',
         'const SafeSTLVector<BernsteinOperator> *',
         'SafeSTLVector<const BernsteinOperator *>',
         'SafeSTLVector<SafeSTLVector<BernsteinOperator>>',
         'SafeSTLVector<const SafeSTLVector<BernsteinOperator> *>',
         'BasisValues1d',
         'DenseMatrix',
         'SafeSTLArray<Real,2>',
         'std::pair<Real,Real>',
         'BasisEndBehaviour']

flags = ['grid_element::Flags',
         'domain_element::Flags',
         'space_element::Flags',
         'grid_function_element::Flags',
         'function_element::Flags'];




for dim in inst.all_domain_dims:
    for type in types:
        array = 'SafeSTLArray<%s,%d>' %(type,dim)
        array_list.append(array)
    for k in range(0,dim+1):
        sub_elem = 'typename UnitElement<%d>::template SubElement<%d>' % (dim,k)
        n_sub_elems = skel_size(dim,k)
        array = 'SafeSTLArray<%s,%s>' %(sub_elem,n_sub_elems)
        array_list.append(array)
        array = 'SafeSTLArray<int,%s>' %(n_sub_elems)
        array_list.append(array)
        
        space_dims = [1,dim]
        
        for space_dim in space_dims:
            grid_func = 'GridFunction<%d,%d>' %(dim,space_dim)
            grid_func_D0 = 'typename %s::Value' %(grid_func)
            grid_func_D1 = 'typename %s::template Derivative<1>' %(grid_func)
            grid_func_D2 = 'typename %s::template Derivative<2>' %(grid_func)
            cache_D0 = 'boost::fusion::pair<grid_function_element::_D<0>,DataWithFlagStatus<ValueVector<%s>>>' % (grid_func_D0)
            cache_D1 = 'boost::fusion::pair<grid_function_element::_D<1>,DataWithFlagStatus<ValueVector<%s>>>' % (grid_func_D1)
            cache_D2 = 'boost::fusion::pair<grid_function_element::_D<2>,DataWithFlagStatus<ValueVector<%s>>>' % (grid_func_D2)
            array = 'SafeSTLArray<PointValuesCache<%d,boost::fusion::map<%s,%s,%s>>,%d>' %(dim+1,cache_D0,cache_D1,cache_D2,n_sub_elems)
            array_list.append(array)
        
        for flag in flags:
            array = 'SafeSTLArray<%s,%s>' %(flag,dim+1)
            array_list.append(array)
        


for array in unique(array_list):
    f.write('template class %s; \n' % (array))
#    f.write('template LogStream &operator<<(LogStream &, const %s &); \n' % (row))
#    f.write('template %s operator+(const %s &, const %s &); \n' % (row,row,row))
#    f.write('template %s operator+(const %s &, const Index); \n' % (row,row))
#    f.write('template %s operator-(const %s &, const Index); \n' % (row,row))


