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


 
 

elements = set()
element_funcs = set()

elem = 'PhysicalBasisElement<0,0,1,0>'
elements.add(elem)


func = 'const ValueVector<Real> %s::get_w_measures<0>(const int) const' % (elem)
element_funcs.add(func)
func = 'const ValueVector<Real> & %s::get_measures<0>(const int) const' % (elem)
element_funcs.add(func)



for basis in inst.SubPhysBases + inst.PhysBases:
    x = basis.spec
    space_dim = x.dim + x.codim
    elem = 'PhysicalBasisElement<%d,%d,%d,%d>' %(x.dim,x.range,x.rank,x.codim)
    elements.add(elem)
    for k in range(0,x.dim+1):
        func = 'const ValueVector<Real> %s::get_w_measures<%d>(const int) const' % (elem,k)
        element_funcs.add(func)
        func = 'const ValueVector<Real> & %s::get_measures<%d>(const int) const' % (elem,k)
        element_funcs.add(func)
        if (k+1 == x.dim):
            func = 'const ValueVector<Points<%d>> & %s::get_boundary_normals<%d>(const int) const' % (space_dim,elem,k)
            element_funcs.add(func)





for elem in elements:
    f.write('template class %s;\n' %elem)


for func in element_funcs:
    f.write('template %s;\n' %func)

      
