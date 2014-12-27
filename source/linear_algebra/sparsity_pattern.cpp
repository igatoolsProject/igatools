//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2015  by the igatools authors (see authors.txt).
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

#ifdef USE_PETSC

using std::set;

using std::pair;

IGA_NAMESPACE_OPEN

SparsityPattern::
SparsityPattern(const SpaceManager &space_manager)
{
    Assert(space_manager.is_spaces_insertion_open() == false,ExcInvalidState());

    // build the dofs graph

    //-----------------------------------------------------------
    // adding the row dofs id -- begin
    {
        const auto set_row_dofs = space_manager.get_row_dofs();

        row_dofs_.insert(row_dofs_.end(),set_row_dofs.begin(),set_row_dofs.end());
    }
    // adding the row dofs id -- end
    //-----------------------------------------------------------

    //-----------------------------------------------------------
    // adding the col dofs id -- begin
    {
        const auto set_col_dofs = space_manager.get_col_dofs();

        col_dofs_.insert(col_dofs_.end(),set_col_dofs.begin(),set_col_dofs.end());
    }
    // adding the col dofs id -- end
    //-----------------------------------------------------------


    //-----------------------------------------------------------
    const auto &spaces_connections = space_manager.get_spaces_connections();
    Assert(!spaces_connections.empty(),ExcEmptyObject());
    for (const auto &sp_conn : spaces_connections)
    {
        if (sp_conn.is_unique_space())
        {
            // adding the contribution of the dofs defined within the space itself-- begin
            const auto &space = sp_conn.get_space_row();
            for (const auto element_dofs : space.get_elements_dofs_view())
                for (const auto &dof : element_dofs.second)
                    (*this)[dof].insert(element_dofs.second.begin(),element_dofs.second.end());
            // adding the contribution of the dofs defined within the space -- end
        }



        // adding the extra contribution to the connectivity defined within the spaces connection -- begin
        const auto &extra_dofs_connectivity = sp_conn.get_extra_dofs_connectivity();
        for (const auto &connectivity_map_entry : extra_dofs_connectivity)
        {
            const auto   row_id = connectivity_map_entry.first;
            const auto &cols_id = connectivity_map_entry.second;

            (*this)[row_id].insert(cols_id.begin(),cols_id.end());
        }
        // adding the extra contribution to the connectivity defined within the spaces connection -- end
    }
    //-----------------------------------------------------------
}


#if 0
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
#endif

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
SparsityPattern::get_num_dof_connections() const
{

    vector< long unsigned int > num_dof_connections ;

    for (const auto &map_entry : (*this))
        num_dof_connections.emplace_back(map_entry.second.size()) ;

    return (num_dof_connections) ;
}

IGA_NAMESPACE_CLOSE

#endif
