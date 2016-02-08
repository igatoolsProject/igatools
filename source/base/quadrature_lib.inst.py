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
###############################################################################

# QA (pauletti, Mar 19, 2014):
from init_instantiation_data import *
data = Instantiation()
f = data.file_output
inst = data.inst


quad_types = ['QGauss','QGaussLobatto','QUniform','QTrapez']

quad_classes = []

for dim in inst.all_domain_dims:
    for quad_type in quad_types:
        quad = '%s<%d>' % (quad_type,dim)
        quad_classes.append(quad)
   
for quad in unique(quad_classes):   
    f.write('template class %s;\n' % (quad))
    

#f.write('IGA_NAMESPACE_CLOSE\n')
#
#for quad in unique(quad_classes):   
#    f.write('extern template class std::shared_ptr<iga::%s>;\n' % (quad))
#
#f.write('IGA_NAMESPACE_OPEN\n')

