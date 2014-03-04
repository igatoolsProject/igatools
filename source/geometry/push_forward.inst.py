#-+--------------------------------------------------------------------
# Igatools a general purpose Isogeometric analysis library.
# Copyright (C) 2012-2014  by the igatools authors (see authors.txt).
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

# QA (pauletti, Mar 4, 2014 ):
from init_instantiation_data import *
file_output, inst = intialize_instantiation()
file_output.write('IGA_NAMESPACE_OPEN\n')

pf_dims = unique( [ [x.dim, x.codim, x.trans_type] for x in 
                   inst.all_table + inst.extended_table] )
all_pf_args = [PForwRow(x) for x in pf_dims]
# We need the pushforwards for all physical spaces and
# all reference spaces
pf_args = ['Transformation::%s, %d, %d'
            %(x.trans_type, x.dim, x.codim) for x in all_pf_args]
# pf_args = pf_args + unique(['Transformation::h_grad, %d, 0' %(x.dim)
#                             for x in inst.all_ref_sp_dims])
for pf in unique(pf_args):
    file_output.write('template class PushForward<%s> ;\n' %(pf))

file_output.write('IGA_NAMESPACE_CLOSE\n')
file_output.close()
