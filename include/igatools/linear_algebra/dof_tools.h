//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2016  by the igatools authors (see authors.txt).
//
// This file is part of the igatools library.
//
// The igatools library is free software: you can use it, redistribute
// it and/or modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation, either
// version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//-+--------------------------------------------------------------------


#ifndef __DOF_TOOLS_H_
#define __DOF_TOOLS_H_

#include <igatools/base/config.h>
#include <igatools/linear_algebra/epetra.h>

IGA_NAMESPACE_OPEN

#ifdef USE_TRILINOS

/**
 * Collection of routines to handle the relation
 * between the linear algebra and the degrees of freedom of
 * the spaces.
 */
namespace dof_tools
{
using namespace EpetraTools;
/**
 * Modifies the matrix, the unknown and rhs of a linear system
 * to impose dirichlet constraints on the dofs.
 * todo: //TODO: apply_dirichlet_constraint? and document more.
 */
void apply_boundary_values(const std::map<Index,Real> &boundary_values,
                           Matrix &matrix,
                           Vector &rhs,
                           Vector &solution);

} // end of namespace dof_tools

#endif // USE_TRILINOS

IGA_NAMESPACE_CLOSE

#endif /* __DOF_TOOLS_H_ */
