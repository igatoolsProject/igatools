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
###############################################################################

# QA (pauletti, Mar 19, 2014):
from init_instantiation_data import *

include_files = ['base/config.h']

data = Instantiation(include_files)
f = data.file_output
inst = data.inst

f.write("""
template <int dim> class Grid;
template <int dim, int range> class GridFunction;
template <int dim, int codim> class Domain;
template <int dim, int range, int rank> class SplineSpace;
template <int dim, int range, int rank> class BSpline;
template <int dim, int range, int rank> class NURBS;
template <int dim, int range, int rank> class ReferenceSpaceBasis;
template <int dim, int range, int rank, int codim> class PhysicalSpaceBasis;
template <int dim, int codim, int range, int rank> class Function;

""")

f.write("""
/**
 * @brief This is a helper struct for defining containers for all
 * the igatools instantiated types.
 *
 * This is a helper struct for defining containers for all
 * the igatools instantiated types.
 *
 * The types are stored as shared pointers into <tt>boost::fusion</tt>
 * vectors.
 *
 * @author P. Antolin
 * @date 2015
 */
struct InstantiatedTypes
{
public:
""")


# Grids
f.write("""

  /** All grid instantiations. */
  typedef boost::mpl::vector<""")

types = []
for dim in inst.domain_dims:
    types.append('Grid<%d>' % (dim))
for dim in inst.sub_domain_dims:
    types.append('Grid<%d>' % (dim))

types = unique(types)
f.write(" %s,\n" % (types[0]));
for i in range(1, len(types) - 1):
  f.write("                              %s,\n" % (types[i]));
f.write("                              %s >  Grids;\n" % (types[-1]));



# Spline spaces
f.write("""

  /** All spline space instantiations. */
  typedef boost::mpl::vector<""")

types = []
types.append('SplineSpace<%d, %d, %d>' % (0, 0, 1))
for x in inst.sub_ref_sp_dims:
    types.append('SplineSpace<%d, %d, %d>' % (x.dim, x.range, x.rank))
for x in inst.ref_sp_dims:
    types.append('SplineSpace<%d, %d, %d>' % (x.dim, x.range, x.rank))

types = unique(types)
f.write(" %s,\n" % (types[0]));
for i in range(1, len(types) - 1):
  f.write("                              %s,\n" % (types[i]));
f.write("                              %s >  SplineSpaces;\n" % (types[-1]));




# Reference spaces

f.write("""

  /** All reference space basis instantiations. */
  typedef boost::mpl::vector<""")

types = []
types.append('ReferenceSpaceBasis<%d, %d, %d>' % (0, 0, 1))
for x in inst.sub_ref_sp_dims:
    types.append('ReferenceSpaceBasis<%d, %d, %d>' % (x.dim, x.range, x.rank))
for x in inst.ref_sp_dims:
    types.append('ReferenceSpaceBasis<%d, %d, %d>' % (x.dim, x.range, x.rank))

types = unique(types)
f.write(" %s,\n" % (types[0]));
for i in range(1, len(types) - 1):
  f.write("                              %s,\n" % (types[i]));
f.write("                              %s >  RefSpaceBases;\n" % (types[-1]));




# Grid functions

f.write("""

  /** All grid function instantiations. */
  typedef boost::mpl::vector<""")

types = []
for x in inst.sub_mapping_dims:
    types.append('GridFunction<%d, %d>' % (x.dim, x.space_dim))
for x in inst.mapping_dims:
    types.append('GridFunction<%d, %d>' % (x.dim, x.space_dim))
    # The next dimensions are needed by NURBS
    types.append('GridFunction<%d, %d>' % (x.dim, 1))

types = unique(types)
f.write(" %s,\n" % (types[0]));
for i in range(1, len(types) - 1):
  f.write("                              %s,\n" % (types[i]));
f.write("                              %s >  GridFunctions;\n" % (types[-1]));




# Domains

f.write("""

  /** All domain instantiations. */
  typedef boost::mpl::vector<""")

types = []
for x in inst.sub_mapping_dims:
    types.append('Domain<%d, %d>' % (x.dim, x.codim))
for x in inst.mapping_dims:
    types.append('Domain<%d, %d>' % (x.dim, x.codim))

types = unique(types)
f.write(" %s,\n" % (types[0]));
for i in range(1, len(types) - 1):
  f.write("                              %s,\n" % (types[i]));
f.write("                              %s >  Domains;\n" % (types[-1]));




# Physical spaces

f.write("""

  /** All physical space basis instantiations. */
  typedef boost::mpl::vector<""")

types = []
types.append('PhysicalSpaceBasis<%d, %d, %d, %d>' % (0, 0, 1, 0))
for sp in inst.SubPhysSpaces:
    types.append('PhysicalSpaceBasis<%d, %d, %d, %d>' % (sp.spec.dim, sp.spec.range, sp.spec.rank, sp.spec.codim))
for sp in inst.PhysSpaces:
    types.append('PhysicalSpaceBasis<%d, %d, %d, %d>' % (sp.spec.dim, sp.spec.range, sp.spec.rank, sp.spec.codim))

types = unique(types)
f.write(" %s,\n" % (types[0]));
for i in range(1, len(types) - 1):
  f.write("                              %s,\n" % (types[i]));
f.write("                              %s >  PhysSpaces;\n" % (types[-1]));




# Functions

f.write("""

  /** All function instantiations. */
  typedef boost::mpl::vector<""")
types = []
# Functions
for dims in inst.all_function_dims:
    types.append('Function<%d, %d, %d, %d>' % (dims.dim, dims.codim, dims.range, dims.rank))

f.write(" %s,\n" % (types[0]));
for i in range(1, len(types) - 1):
  f.write("                              %s,\n" % (types[i]));
f.write("                              %s >  Functions;\n" % (types[-1]));



f.write("""
};

""")