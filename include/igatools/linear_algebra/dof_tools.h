//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2014  by the igatools authors (see authors.txt).
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
#include <igatools/linear_algebra/sparsity_pattern.h>
#include <igatools/basis_functions/function_space.h>


#include <memory>

IGA_NAMESPACE_OPEN

template <LAPack la_pack>
class Vector;

template <LAPack la_pack>
class Matrix;

/**
 * Collection of routines to handle the relation
 * between the linear algebra and the degrees of freedom of
 * the spaces.
 */
namespace dof_tools
{
/**
 * Modifies the matrix, the unknown and rhs of a linear system
 * to impose dirichlet constraints on the dofs.
 * todo: //TODO: apply_dirichlet_constraint? and document more.
 */
template <LAPack la_pack>
void apply_boundary_values(const std::map<Index,Real> &boundary_values,
                           Matrix<la_pack> &matrix,
                           Vector<la_pack> &rhs,
                           Vector<la_pack> &solution);


/**
 * This function returns the connectivity of two sets of dofs index, using an all-to-all
 * strategy, i.e. each dof in the set @p row_dofs is connected with all dofs in the set @p col_dofs.
 */
inline
std::map<Index,std::set<Index> >
get_dofs_connectvity_all_to_all(const std::set<Index> &row_dofs,const std::set<Index> &col_dofs)
{
    std::map<Index,std::set<Index>> dofs_connectivity;
    for (const Index row_dof : row_dofs)
        dofs_connectivity[row_dof].insert(col_dofs.begin(),col_dofs.end());

    return dofs_connectivity;
}


} // end of namespace dof_tools

IGA_NAMESPACE_CLOSE

#endif /* __DOF_TOOLS_H_ */
