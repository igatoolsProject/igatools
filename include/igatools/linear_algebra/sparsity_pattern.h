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

#ifndef SPARSITY_PATTERN_H_
#define SPARSITY_PATTERN_H_

#include <igatools/base/config.h>
#include <igatools/basis_functions/dofs_manager.h>

#include <set>
#include <vector>
#include <map>

IGA_NAMESPACE_OPEN

/** @addtogroup linear_algebra
 *@{
*/

/**
 *
 * @todo SparsityPattern should be a wrapper to the Trilinos sparsity pattern
 * @todo Missing documentation
 * @todo Re-design in order to be used with different linear algebra package
 * (Trilinos, PETSc,etc.)
 *
 */
class SparsityPattern
    : public std::map< Index, std::set< Index > >
{
public:

    SparsityPattern() = delete;

    /**
     * Constructs the SparsityPattern associated with the DofsManager of one space.
     */
    template<int space_dim>
    SparsityPattern(const DofsManager<space_dim> &dofs_manager);

    /**
     * Construct the sparsity pattern associated with the DofsManager of two space.
     *
     * @warning This function only works when both spaces have the same number of elements.
     */
    template<int space_dim>
    SparsityPattern(const DofsManager<space_dim> &dofs_manager_rows,const DofsManager<space_dim> &dofs_manager_cols);
    /*
        SparsityPattern(const std::vector< Index > row_dofs,
                        const std::vector< Index > col_dofs) ;
    //*/
    SparsityPattern(const SparsityPattern &) ;

    /**
     * Returns the number of dofs in rows.
     */
    Size get_num_row_dofs() const ;

    /**
     * Returns the number of dofs in columns.
     */
    Size get_num_col_dofs() const ;


    /**
     * Returns the vector of the column dofs.
     */
    const std::vector< Index > get_col_dofs() const ;

    /**
     * Returns the vector of the row dofs.
     */
    const std::vector< Index > get_row_dofs() const ;


    /**
     * todo: document me.
     */
    std::vector< long unsigned int > get_num_overlapping_funcs() const ;


private:

    /**
     * todo: document me.
     */
    SparsityPattern &operator=(const SparsityPattern &) ;

    std::vector< Index > row_dofs_ ;
    std::vector< Index > col_dofs_ ;


};
/**@}*/



IGA_NAMESPACE_CLOSE

#endif
