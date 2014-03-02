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

#include <igatools/linear_algebra/sparsity_pattern.h>

using std::set;
using std::vector;

IGA_NAMESPACE_OPEN

SparsityPattern::SparsityPattern(const std::vector< Index > row_dofs,
                                 const std::vector< Index > col_dofs)
    :
    map< Index, set< Index > >(),
    row_dofs_(row_dofs),
    col_dofs_(col_dofs)
{}


SparsityPattern::SparsityPattern(const SparsityPattern &sparsity_pattern)
    :
    map< Index, set< Index > >(sparsity_pattern),
    row_dofs_(sparsity_pattern.row_dofs_),
    col_dofs_(sparsity_pattern.col_dofs_)
{}

int
SparsityPattern::get_num_row_dofs() const
{
    return (row_dofs_.size()) ;
}

int
SparsityPattern::get_num_col_dofs() const
{
    return (col_dofs_.size()) ;
}


const vector< Index >
SparsityPattern::get_col_dofs() const
{
    return (col_dofs_) ;
}

const vector< Index >
SparsityPattern::get_row_dofs() const
{
    return (row_dofs_) ;
}


vector< long unsigned int >
SparsityPattern::get_num_overlapping_funcs() const
{

    vector< long unsigned int > num_overlapping_funcs ;

    auto dof     = row_dofs_.cbegin() ;
    auto dof_end = row_dofs_.cend() ;

    for (; dof != dof_end ; ++dof)
        num_overlapping_funcs.emplace_back(this->at(*dof).size()) ;

    return (num_overlapping_funcs) ;
}

IGA_NAMESPACE_CLOSE

