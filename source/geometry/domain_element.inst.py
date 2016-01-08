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

include_files = ['../../source/geometry/grid_iterator.cpp']
data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)


elements = set()
element_funcs = set()

for x in inst.sub_mapping_dims + inst.mapping_dims:
    space_dim = x.dim + x.codim
    elem = 'DomainElement<%d,%d>' %(x.dim, x.codim)
    elements.add(elem)
    for k in range(0,x.dim+1):
        func = 'const ValueVector<typename %s::Point> & %s::get_points<%d>(const int) const' % (elem,elem,k)
        element_funcs.add(func)
        
        func = 'const ValueVector<typename %s::Jacobian> & %s::get_jacobians<%d>(const int) const' % (elem,elem,k)
        element_funcs.add(func)

        func = 'const ValueVector<typename %s::Hessian> & %s::get_hessians<%d>(const int) const' % (elem,elem,k)
        element_funcs.add(func)

        func = 'const ValueVector<Real> & %s::get_measures<%d>(const int) const' % (elem,k)
        element_funcs.add(func)

        func = 'const ValueVector<Real> & %s::get_w_measures<%d>(const int) const' % (elem,k)
        element_funcs.add(func)

        if (k+1 == x.dim):
            func = 'const ValueVector<Points<%d> > & %s::get_boundary_normals<%d>(const int, EnableIf<(%d >= 0)> *) const' %(space_dim,elem,k,k)
            element_funcs.add(func)


 

#accs=  ['DomainElement']
#iters =  ['GridIterator']
#for x in inst.sub_mapping_dims+inst.mapping_dims:
#  for i in range(len(accs)):
#    acc = iters[i] + '<' + accs[i]+ '<%d,%d>' %(x.dim, x.codim) + '>' 
#    f.write('template class %s; \n' %(acc))
    

for element in elements:
    f.write('template class %s;\n' %(element))
    f.write('template class GridIterator<%s>;\n' %(element))

for func in element_funcs:
    f.write('template %s;\n' %(func))
