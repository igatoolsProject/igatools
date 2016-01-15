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

include_files = ['geometry/grid.h',
                 'geometry/grid_element.h',
                 'basis_functions/basis_element.h',
                 '../../source/geometry/grid_iterator.cpp']
data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)


sub_dim_members = \
 ['ValueVector<Real> elem::get_w_measures<k>(const int) const']

elements = set()
elements.add('BasisElement<0,0,0,1>')

#templated_funcs = ['ValueVector<Real> BasisElement<0,0,0,1>::get_w_measures<0>(const int) const']

element_funcs = set()


VTypes = ['basis_element::_Value','basis_element::_Gradient','basis_element::_Hessian','basis_element::_Divergence']



#--------------------------------------------------------------------------------------
# BasisElement used by ReferenceBasisElement 
for x in inst.sub_ref_sp_dims + inst.ref_sp_dims:
    elem = 'BasisElement<%d,0,%d,%d>' %(x.dim, x.range, x.rank)
    elements.add(elem)
    for k in range(0,x.dim+1):
        func = 'ValueVector<Real> %s::get_w_measures<%d>(const int) const' % (elem,k)
        element_funcs.add(func)
            
#--------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------
# BasisElement used by PhysicalBasisElement 

for space in inst.SubPhysBases + inst.PhysBases:
    x = space.spec
    elem = 'BasisElement<%d,%d,%d,%d>' %(x.dim,x.codim,x.range, x.rank)
    elements.add(elem)
    for k in range(0,x.dim+1):
        func = 'ValueVector<Real> %s::get_w_measures<%d>(const int) const' % (elem,k)
        element_funcs.add(func)
#--------------------------------------------------------------------------------------



#---------------------------------------------------
iters =  ['GridIterator']

for elem in elements:
    f.write('template class %s;\n' %elem)
    f.write('template class GridIterator<%s>;\n' %elem)


for func in element_funcs:
    f.write('template %s;\n' %func)
#---------------------------------------------------

