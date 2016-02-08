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

include_files = ['../../source/geometry/grid_iterator.cpp']
data = Instantiation(include_files)
(f, inst) = (data.file_output, data.inst)

sub_dim_members = [
             'Real Element::get_measure<k>(const int j) const;',
             'const ValueVector<Real> & Element::get_weights<k>(const int j) const;',
             'const Points<k> Element::get_side_lengths<k>(const int j) const;',
             'const ValueVector<typename Element::Point> & Element::get_points<k>(const int j) const;',
             'bool Element::is_boundary<k>(const Index sub_elem_id) const;',
             'bool Element::is_boundary<k>() const;',
             'std::shared_ptr<const Quadrature<k>> Element::get_quad<k>() const']

elements = set()
element_funcs = set()

for dim in inst.domain_dims + inst.sub_domain_dims:
    elem = 'GridElement<%d>' %(dim)
    elements.add(elem)
    for fun in sub_dim_members:
        for k in range(0,dim+1):
          s = fun.replace('k', '%d' % (k)).replace('Element', '%s' % (elem));
          element_funcs.add(s)

 

for elem in elements:
    f.write('template class %s;\n' %(elem))
    f.write('template class GridIterator<%s>;\n' %(elem))

for func in element_funcs:
    f.write('template %s;\n' %(func))
  
