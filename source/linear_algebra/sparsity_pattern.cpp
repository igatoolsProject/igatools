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

using std::pair;

IGA_NAMESPACE_OPEN

SparsityPattern::
SparsityPattern(const SpaceManager &space_manager)
{
    Assert(space_manager.is_spaces_insertion_open() == false,ExcInvalidState());

    // build the dofs graph

    //-----------------------------------------------------------
    // adding the dofs id -- begin
    const auto &dofs_view = space_manager.get_dofs_view();

    for (const auto &dof : dofs_view)
        row_dofs_.push_back(dof);
    // adding the dofs id -- end
    //-----------------------------------------------------------


    //-----------------------------------------------------------
    // adding the linear constraints id -- begin
    const auto &linear_constraints = space_manager.get_linear_constraints();

    Index row_id = space_manager.get_num_dofs();
    for (const auto &lc : linear_constraints)
        row_dofs_.push_back(row_id++);
    // adding the linear constraints id -- end
    //-----------------------------------------------------------



    //-----------------------------------------------------------
    // copying the row ids to the col ids -- begin
    Assert(!row_dofs_.empty(),ExcEmptyObject());
    col_dofs_ = row_dofs_;
    // copying the row ids to the col ids -- end
    //-----------------------------------------------------------





    using DofsInRow = set<Index>;
    DofsInRow empty_set;

    // adding the global dof keys to the map representing the dof connectivity
    for (const auto &dof : dofs_view)
        this->insert(pair<Index,DofsInRow>(dof,empty_set));


    //-----------------------------------------------------------
    // adding the DOF-DOF contribution -- begin
    const auto &spaces_info = space_manager.get_spaces_info();
    Assert(!spaces_info.empty(),ExcEmptyObject());
    for (const auto &space : spaces_info)
        for (const auto element_dofs : space.second->get_elements_dofs_view())
            for (const auto &dof : element_dofs.second)
                (*this)[dof].insert(element_dofs.second.begin(),element_dofs.second.end());
    // adding the DOF-DOF contribution -- end
    //-----------------------------------------------------------




    //-----------------------------------------------------------
    // adding the LC-DOF/DOF-LC contributions -- begin
    row_id = space_manager.get_num_dofs();
    for (const auto &lc : linear_constraints)
    {
        const auto lc_dofs = lc->get_dofs_id();

        //-----------------------------------------------------------
        // adding the LC-DOF contribution -- begin
        for (const auto lc_dof : lc_dofs)
            (*this)[row_id].emplace(lc_dof);
        // adding the LC-DOF contribution -- end
        //-----------------------------------------------------------


        //-----------------------------------------------------------
        // adding the DOF-LC contribution -- begin
        for (const auto lc_dof : lc_dofs)
            (*this)[lc_dof].emplace(row_id);
        // adding the DOF-LC contribution -- end
        //-----------------------------------------------------------


        row_id++;
    }
    // adding the LC-DOF/DOF-LC contributions -- end
    //-----------------------------------------------------------


}

SparsityPattern::
SparsityPattern(const SpaceManager &space_manager_rows,const SpaceManager &space_manager_cols)
{
    Assert(space_manager_rows.is_spaces_insertion_open() == false,ExcInvalidState());
    Assert(space_manager_cols.is_spaces_insertion_open() == false,ExcInvalidState());

    // build the dofs graph
    const auto &dofs_view_rows = space_manager_rows.get_dofs_view();
    const auto &dofs_view_cols = space_manager_cols.get_dofs_view();

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



    const auto &spaces_rows_info = space_manager_rows.get_spaces_info();
    Assert(!spaces_rows_info.empty(),ExcEmptyObject());

    const auto &spaces_cols_info = space_manager_cols.get_spaces_info();
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
        const auto dofs_elements_view_space_row = space_row->get_elements_dofs_view();

        const auto &space_col = space_col_iterator->second;
        const auto dofs_elements_view_space_col = space_col->get_elements_dofs_view();

        //check the equality of num. elements on each patch
        Assert(dofs_elements_view_space_row.size() == dofs_elements_view_space_col.size(),
               ExcDimensionMismatch(dofs_elements_view_space_row.size(),dofs_elements_view_space_col.size()));

        auto dofs_row_iterator     = dofs_elements_view_space_row.cbegin();
        auto dofs_row_iterator_end = dofs_elements_view_space_row.end();
        auto dofs_col_iterator = dofs_elements_view_space_col.cbegin();

        for (; dofs_row_iterator != dofs_row_iterator_end ; ++dofs_row_iterator, ++dofs_col_iterator)
            for (const auto &dof_row : dofs_row_iterator->second)
                (*this)[dof_row].insert(dofs_col_iterator->second.cbegin(),dofs_col_iterator->second.cend());
    }
}

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

