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
data = Instantiation()
(f, inst) = (data.file_output, data.inst)

funcs = set ()

for dim in inst.domain_dims:
    funcs.add ('grid_functions::BallGridFunction<%d>' % (dim) )    
    funcs.add ('grid_functions::ConstantGridFunction<%d,%d>' % (dim,dim) )
    funcs.add ('grid_functions::ConstantGridFunction<%d,1>' % (dim) )
    funcs.add ('grid_functions::IdentityGridFunction<%d>' % (dim) )

for dim in inst.sub_domain_dims:
    funcs.add ('grid_functions::ConstantGridFunction<%d,%d>' % (dim,dim) )
    funcs.add ('grid_functions::IdentityGridFunction<%d>' % (dim) )


funcs.add ('grid_functions::LinearGridFunction<0,0>')
for dim in inst.domain_dims:
    funcs.add ('grid_functions::LinearGridFunction<%d,1>' % (dim) )
    funcs.add ('grid_functions::ConstantGridFunction<%d,1>' % (dim) )
    for sub_dim in range(0,dim+1):
        funcs.add ('grid_functions::LinearGridFunction<%d,%d>' % (sub_dim,dim) )
        funcs.add ('grid_functions::ConstantGridFunction<%d,%d>' % (sub_dim,dim) )


funcs.add('grid_functions::SphereGridFunction<1>')
funcs.add('grid_functions::SphereGridFunction<2>')


funcs.add('grid_functions::TriangleGridFunction<2>')
funcs.add('grid_functions::TriangleGridFunction<3>')



for func in funcs:
    f.write('template class %s;\n' %(func))


#---------------------------------------------------
#f.write('#ifdef IGATOOLS_WITH_SERIALIZATION\n')
#archives = ['OArchive','IArchive']
#
#for func in unique(funcs):
#    for ar in archives:
#        f.write('template void %s::serialize(%s&);\n' %(func,ar))
#f.write('#endif // IGATOOLS_WITH_SERIALIZATION\n')
#---------------------------------------------------



