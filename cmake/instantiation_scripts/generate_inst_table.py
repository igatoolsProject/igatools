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


## @package generate_inst_table
# Script to create an instantiantion table that...
#
# @author: antolin, 2013
# @author: pauletti, 2013

 # Getting a dictionary or arguments.
from sys import argv as sysargv
from sys import exit
from os import sep as ossep
args = dict([arg.split('=') for arg in sysargv[1:]])
for k in args.keys():
    args[k] = [int(i) for i in args[k].split()]


# Openning the output file.
file_output = open('instantiation_table.txt', 'w')
file_output.write('# Each line describes a physical discrete space that the\n')
file_output.write('# library want to be compile for. \n')
file_output.write('# dim   codim  range  rank trans_type\n')


#if len(args['dim']) == 0 and len(args['range']) == 0 and len(args['rank']) == 0:
file_output.write(
 '''1          0          1         1        h_grad
    1          1          1         1        h_grad
    1          1          2         1        h_grad
    1          1          3         1        h_grad
    1          2          1         1        h_grad
    1          2          3         1        h_grad
    2          0          1         1        h_grad
    2          0          2         1        h_grad
    2          1          1         1        h_grad
    2          1          3         1        h_grad
    3          0          1         1        h_grad
    3          0          3         1        h_grad''')
print('Default instantiation table was created.')

#elif args.has_key('dom_dim') and args.has_key('space_dim') and args.has_key('range_rank') and \
#        len(args['dom_dim']) == len(args['space_dim']) and len(args['dom_dim']) == len(args['range_rank']) \
#        and len(args['dom_dim']) > 0:
#    for i in range(len(args['dom_dim'])):
#        file_output.write('    %d          %d          %d         %d\n' %
#                (args['dom_dim'][i],
#                 args['space_dim'][i],
#                 1 if args['range_rank'][i] == 0 else args['space_dim'][i],
#                 args['range_rank'][i]))
#
#    print('User defined instantiation table was created.')

#else:
#    file_output.write('''    1          1          1         0
#    1          2          1         0
#    1          2          2         1
#    1          3          1         0
#    1          3          3         1
#    2          2          1         0
#    2          2          2         1
#    2          3          1         0
#    2          3          3         1
#    3          3          1         0
#    3          3          3         1''')
#    print('!!!  There is an error in the user defined instantiation table.')
#    print('Default instantiation table was created.')


file_output.close()
