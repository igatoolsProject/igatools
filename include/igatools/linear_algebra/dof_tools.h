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

class Vector;
class Matrix;


/**
 * Collection of routines to handle the relation
 * between the linear algebra and the degrees of freedom of
 * the spaces.
 */
namespace dof_tools
{

/**
 * Construct the sparsity pattern associated with the degrees
 * of freedom of the space.
 * @note This function is called when SpaceType is derived from FunctionSpace.
 */
template < class SpaceType >
SparsityPattern
get_sparsity_pattern(const SpaceType &space, Enable_if<Is_function_space<SpaceType>()>* = nullptr);


/**
 * Construct the sparsity pattern associated with the degrees
 * of freedom of the space.
 * @note This function is called when SpaceType is not derived from FunctionSpace.
 */
template < class SpaceType >
SparsityPattern
get_sparsity_pattern(const SpaceType &space, Enable_if<!Is_function_space<SpaceType>()>* = nullptr)
{
    AssertThrow(false,ExcMessage("The function argument is not a function space."));
    SparsityPattern graph;
    return graph;
}


/**
 * @todo document this function
 * @note This function is called when SpaceType is derived from FunctionSpace.
 */
template < class SpaceType >
SparsityPattern
get_sparsity_pattern(const std::vector< std::shared_ptr< SpaceType > > &space,
                     Enable_if<Is_function_space<SpaceType>()>* = nullptr);

/**
 * @todo document this function
 * @note This function is called when SpaceType is derived from FunctionSpace.
 */
template < class SpaceType >
SparsityPattern
get_sparsity_pattern(const std::vector< std::shared_ptr< SpaceType > > &space,
                     Enable_if<!Is_function_space<SpaceType>()>* = nullptr)
{
    AssertThrow(false,ExcMessage("The function argument is not a function space."));
    SparsityPattern graph;
    return graph;
}


/**
 * Returns the dof ids of the space.
 * @note This function is called when SpaceType is derived from FunctionSpace.
 */
template < class SpaceType >
std::vector<Index>
get_dofs(const SpaceType &space, Enable_if<Is_function_space<SpaceType>()>* = nullptr);


/**
 * Returns the dof ids of the space.
 * @note This function is called when SpaceType is not derived from FunctionSpace.
 */
template < class SpaceType >
std::vector<Index>
get_dofs(const SpaceType &space, Enable_if<!Is_function_space<SpaceType>()>* = nullptr)
{
    AssertThrow(false,ExcMessage("The function argument is not a function space."));
    std::vector<Index> vec;
    return vec;
}



/**
 *
 * @warning This function only works when both spaces have the same cartesian grid.
 * @note This function is called when SpaceType is derived from FunctionSpace.
 * @todo document this function
 */
template < class Space1, class Space2 >
inline
SparsityPattern
get_sparsity_pattern(const Space1 &space_rows,
                     const Space2 &space_cols,
                     Enable_if<
                     Is_function_space<Space1>()  &&Is_function_space<Space2>() >* = nullptr)
{
    const auto &row_dofs = get_dofs(space_rows);
    const auto &col_dofs = get_dofs(space_cols);

    SparsityPattern sparsity_pattern(row_dofs, col_dofs);

    // adding the global dof keys to the map representing the dof connectivity
    using set_t = std::set<Index>;
    const set_t empty;
    for (const auto &dof : row_dofs)
        sparsity_pattern.insert(std::pair<Index, set_t>(dof, empty));

    auto element_rows     = space_rows.begin() ;
    const auto element_rows_end = space_rows.end() ;
    auto element_cols     = space_cols.begin() ;

    for (; element_rows != element_rows_end ; ++element_rows, ++element_cols)
    {
        const auto &dofs_element_rows = element_rows->get_local_to_global() ;
        const auto &dofs_element_cols = element_cols->get_local_to_global() ;

        auto dofs_rows_begin = dofs_element_rows.cbegin() ;
        auto dofs_rows_end   = dofs_element_rows.cend() ;

        auto dofs_cols_begin = dofs_element_cols.cbegin() ;
        auto dofs_cols_end   = dofs_element_cols.cend() ;


        for (auto dofs_it = dofs_rows_begin ; dofs_it != dofs_rows_end ; ++dofs_it)
        {
            // get the map element corresponding to the current dof in the
            // current element
            sparsity_pattern[ *dofs_it ].insert(dofs_cols_begin, dofs_cols_end) ;
        }
    }

    return sparsity_pattern;
}

/**
 *
 * @note This function is called when SpaceType is not derived from FunctionSpace.
 * @todo document this function
 */
template < class Space1, class Space2 >
inline
SparsityPattern
get_sparsity_pattern(const Space1 &space_rows,
                     const Space2 &space_cols,
                     Enable_if<
                     !(Is_function_space<Space1>() &&Is_function_space<Space2>())>* = nullptr)
{
    AssertThrow(false,ExcMessage("At least one function argument is not a function space."));
    SparsityPattern graph;
    return graph;
}


/**
 * @warning This function only works when both spaces have the same cartesian grid.
 * @note This function is called when SpaceType is derived from FunctionSpace.
 * @todo document this function
 */
template < class SpaceType >
SparsityPattern get_sparsity_pattern(
    const std::vector< std::shared_ptr< SpaceType > > &space_rows,
    const std::vector< std::shared_ptr< SpaceType > > &space_cols,
    Enable_if<Is_function_space<SpaceType>()>* = nullptr);

/**
 * @warning This function only works when both spaces have the same cartesian grid.
 * @note This function is called when SpaceType is not derived from FunctionSpace.
 * @todo document this function
 */
template < class SpaceType >
SparsityPattern get_sparsity_pattern(
    const std::vector< std::shared_ptr< SpaceType > > &space_rows,
    const std::vector< std::shared_ptr< SpaceType > > &space_cols,
    Enable_if<!Is_function_space<SpaceType>()>* = nullptr)
{
    AssertThrow(false,ExcMessage("The function arguments are not function space."));
    SparsityPattern graph;
    return graph;
}

/**
 * Modifies the matrix, the unknown and rhs of a linear system
 * to impose dirichlet constraints on the dofs.
 * todo: document more.
 */
void apply_boundary_values(const std::map<Index,Real> &boundary_values,
                           Matrix      &matrix,
                           Vector      &rhs,
                           Vector      &solution);
} // end of namespace dof_tools

IGA_NAMESPACE_CLOSE

#endif /* __DOF_TOOLS_H_ */
