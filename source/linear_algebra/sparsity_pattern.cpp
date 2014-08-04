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
using std::pair;

IGA_NAMESPACE_OPEN

SparsityPattern::
SparsityPattern(const DofsManager &dofs_manager)
{
    Assert(dofs_manager.is_space_insertion_open() == false,ExcInvalidState());

    // build the dofs graph
    const auto &dofs_view = dofs_manager.get_dofs_view();

    for (const auto &dof : dofs_view)
        row_dofs_.push_back(dof);

    Assert(!row_dofs_.empty(),ExcEmptyObject());
    col_dofs_ = row_dofs_;

    using DofsInRow = set<Index>;
    DofsInRow empty_set;

    // adding the global dof keys to the map representing the dof connectivity
    for (const auto &dof : dofs_view)
        this->insert(pair<Index,DofsInRow>(dof,empty_set));

    const auto &spaces_info = dofs_manager.get_spaces_info();
    Assert(!spaces_info.empty(),ExcEmptyObject());
    for (const auto &space : spaces_info)
        for (const auto element_dofs : space.second.get_elements_dofs_view())
            for (const auto &dof : element_dofs)
                (*this)[dof].insert(element_dofs.begin(),element_dofs.end());
}

SparsityPattern::
SparsityPattern(const DofsManager &dofs_manager_rows,const DofsManager &dofs_manager_cols)
{
    Assert(dofs_manager_rows.is_space_insertion_open() == false,ExcInvalidState());
    Assert(dofs_manager_cols.is_space_insertion_open() == false,ExcInvalidState());

    // build the dofs graph
    const auto &dofs_view_rows = dofs_manager_rows.get_dofs_view();
    const auto &dofs_view_cols = dofs_manager_cols.get_dofs_view();

    for (const auto &dof : dofs_view_rows)
        row_dofs_.push_back(dof);
    Assert(!row_dofs_.empty(),ExcEmptyObject());

    for (const auto &dof : dofs_view_cols)
        col_dofs_.push_back(dof);
    Assert(!col_dofs_.empty(),ExcEmptyObject());

    using DofsInRow = set<Index>;
    DofsInRow empty_set;

    for (const auto &dof : row_dofs_)
        this->insert(pair<Index,DofsInRow>(dof,empty_set));



    const auto &spaces_rows_info = dofs_manager_rows.get_spaces_info();
    Assert(!spaces_rows_info.empty(),ExcEmptyObject());

    const auto &spaces_cols_info = dofs_manager_cols.get_spaces_info();
    Assert(!spaces_cols_info.empty(),ExcEmptyObject());

    //check the equality of num. patches on each space
    Assert(spaces_rows_info.size() == spaces_cols_info.size(),
           ExcDimensionMismatch(spaces_rows_info.size(),spaces_cols_info.size()));

    auto space_row_iterator = spaces_rows_info.cbegin();
    auto space_row_iterator_end = spaces_rows_info.cend();
    auto space_col_iterator = spaces_cols_info.cbegin();
    for (; space_row_iterator != space_row_iterator_end ; ++space_row_iterator, ++space_col_iterator)
    {
        const auto &space_row = space_row_iterator->second;
        const auto dofs_elements_view_space_row = space_row.get_elements_dofs_view();

        const auto &space_col = space_col_iterator->second;
        const auto dofs_elements_view_space_col = space_col.get_elements_dofs_view();

        //check the equality of num. elements on each patch
        Assert(dofs_elements_view_space_row.size() == dofs_elements_view_space_col.size(),
               ExcDimensionMismatch(dofs_elements_view_space_row.size(),dofs_elements_view_space_col.size()));

        auto dofs_row_iterator     = dofs_elements_view_space_row.cbegin();
        auto dofs_row_iterator_end = dofs_elements_view_space_row.end();
        auto dofs_col_iterator = dofs_elements_view_space_col.cbegin();

        for (; dofs_row_iterator != dofs_row_iterator_end ; ++dofs_row_iterator, ++dofs_col_iterator)
        {
            for (const auto &dof_row : *dofs_row_iterator)
                (*this)[dof_row].insert(dofs_col_iterator->cbegin(),dofs_col_iterator->cend());

        }
    }

    /*
        const Index n_elements = elements_dofs_rows.get_num_elements();
        for (Index ielem = 0 ; ielem < n_elements ; ++ielem)
        {
            const auto &dofs_rows = elements_dofs_rows[ielem];
            const auto &dofs_cols = elements_dofs_cols[ielem];

            for (const auto &dof_row : dofs_rows)
                (*this)[dof_row].insert(dofs_cols.begin(),dofs_cols.end());
        }
        //*/
}

/*
SparsityPattern::SparsityPattern(const std::vector< Index > row_dofs,
                                 const std::vector< Index > col_dofs)
    :
    map< Index, set< Index > >(),
    row_dofs_(row_dofs),
    col_dofs_(col_dofs)
{}
//*/

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

