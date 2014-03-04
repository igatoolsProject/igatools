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


#include <igatools/linear_algebra/dof_tools.h>
#include <igatools/base/exceptions.h>
#include <igatools/linear_algebra/distributed_matrix.h>


using std::vector ;
using std::map ;
using std::set ;
using std::pair ;
using std::shared_ptr ;


IGA_NAMESPACE_OPEN



namespace dof_tools
{

template < class SpaceType >
SparsityPattern
get_sparsity_pattern(const SpaceType &space, EnableIf<Is_function_space<SpaceType>() > *)
{
    //--------------------------------------------------------------------------
    // build the dofs graph
    const vector< Index > &dofs = get_dofs(space) ;

    SparsityPattern sparsity_pattern(dofs, dofs) ;

    typedef set< Index > Set_t ;
    typedef vector< Index > Vec_t ;

    Set_t empty_set ;

    // adding the global dof keys to the map representing the dof connectivity
    for (const auto & dof : dofs)
        sparsity_pattern.insert(pair< Index, Set_t >(dof, empty_set)) ;


    auto element     = space.begin() ;
    const auto element_end = space.end() ;
    for (; element != element_end ; ++element)
    {
        const Vec_t &dofs_element = element->get_local_to_global() ;

        typename Vec_t::const_iterator dofs_begin = dofs_element.cbegin() ;
        typename Vec_t::const_iterator dofs_end   = dofs_element.cend() ;


        for (auto dofs_it = dofs_begin ; dofs_it != dofs_end ; ++dofs_it)
        {
            // get the map element corresponding to the current dof
            sparsity_pattern[ *dofs_it ].insert(dofs_begin, dofs_end) ;
        }
    }

    return (sparsity_pattern) ;
}


template < class SpaceType >
SparsityPattern
get_sparsity_pattern(const vector< shared_ptr< SpaceType > > &space_multipatch,
                     EnableIf<Is_function_space<SpaceType>()> *)
{

    //--------------------------------------------------------------------------
    // build the dofs graph

    typedef set< Index > Set_t ;
    typedef vector< Index > Vec_t ;

    Set_t empty_set ;

    // adding the global dof keys to the map representing the dof connectivity
    Set_t dofs_set;
    for (const auto & space : space_multipatch)
    {
        const vector< Index > &dofs_space = get_dofs(*space) ;
        for (const auto & dof : dofs_space)
            dofs_set.insert(dof) ;
    }

    Vec_t dofs_vector(dofs_set.begin(), dofs_set.end());

    SparsityPattern sparsity_pattern(dofs_vector, dofs_vector) ;
    for (const auto & dof : dofs_vector)
        sparsity_pattern.insert(pair< Index, Set_t >(dof, empty_set)) ;


    // now the keys are initialized, then fill the set of dofs corresponding to each key
    for (const auto & space : space_multipatch)
    {

        auto element     = space->begin() ;
        const auto element_end = space->end() ;
        for (; element != element_end ; ++element)
        {
            const Vec_t &dofs_element = element->get_local_to_global() ;

            typename Vec_t::const_iterator dofs_begin = dofs_element.cbegin() ;
            typename Vec_t::const_iterator dofs_end   = dofs_element.cend() ;


            for (auto dofs_it = dofs_begin ; dofs_it != dofs_end ; ++dofs_it)
            {
                // get the map element corresponding to the current dof
                sparsity_pattern[ *dofs_it ].insert(dofs_begin, dofs_end) ;
            }
        }
    }

    return (sparsity_pattern) ;
}




template < class SpaceType >
SparsityPattern
get_sparsity_pattern(
    const vector< shared_ptr<SpaceType> > &space_multipatch_rows,
    const vector< shared_ptr<SpaceType> > &space_multipatch_cols,
    EnableIf<Is_function_space<SpaceType>()> *)
{
    Assert(space_multipatch_rows.size() == space_multipatch_cols.size(),
           ExcDimensionMismatch(space_multipatch_rows.size(), space_multipatch_cols.size())) ;


    //--------------------------------------------------------------------------
    // build the dofs graph

    typedef set< Index > Set_t ;
    typedef vector< Index > Vec_t ;

    Set_t empty_set ;

    // adding the global dof keys to the map representing the dof connectivity
    Set_t row_dofs_set;
    for (const auto & space : space_multipatch_rows)
    {
        const vector< Index > &dofs_space = get_dofs(*space) ;
        for (const auto & dof : dofs_space)
            row_dofs_set.insert(dof) ;
    }


    Set_t col_dofs_set;
    for (const auto & space : space_multipatch_cols)
    {
        const vector< Index > &dofs_space = get_dofs(*space) ;
        for (const auto & dof : dofs_space)
            col_dofs_set.insert(dof) ;
    }

    Vec_t row_dofs_vector(row_dofs_set.begin(), row_dofs_set.end());
    Vec_t col_dofs_vector(col_dofs_set.begin(), col_dofs_set.end());

    SparsityPattern sparsity_pattern(row_dofs_vector, col_dofs_vector) ;
    for (const auto & dof : row_dofs_vector)
        sparsity_pattern.insert(pair< Index, Set_t >(dof, empty_set)) ;


    // now the keys are initialized, then fill the set of dofs corresponding to each key

    auto space_rows     = space_multipatch_rows.begin() ;
    const auto space_rows_end = space_multipatch_rows.end() ;
    auto space_cols     = space_multipatch_cols.begin() ;

    for (; space_rows != space_rows_end; ++space_rows, ++space_cols)
    {
        auto element_rows     = (*space_rows)->begin() ;
        auto element_rows_end = (*space_rows)->end() ;
        auto element_cols     = (*space_cols)->begin() ;
        for (; element_rows != element_rows_end ; ++element_rows, ++element_cols)
        {
            const Vec_t &dofs_rows_element = element_rows->get_local_to_global() ;
            const Vec_t &dofs_cols_element = element_cols->get_local_to_global() ;

            typename Vec_t::const_iterator dofs_rows_begin = dofs_rows_element.cbegin() ;
            typename Vec_t::const_iterator dofs_rows_end   = dofs_rows_element.cend() ;
            typename Vec_t::const_iterator dofs_cols_begin = dofs_cols_element.cbegin() ;
            typename Vec_t::const_iterator dofs_cols_end   = dofs_cols_element.cend() ;


            for (auto dofs = dofs_rows_begin ; dofs != dofs_rows_end ; ++dofs)
            {
                // get the map element corresponding to the current dof
                sparsity_pattern[ *dofs ].insert(dofs_cols_begin, dofs_cols_end) ;
            }
        }
    }

    return (sparsity_pattern) ;
}


template < class SpaceType >
vector<Index>
get_dofs(const SpaceType &space, EnableIf<Is_function_space<SpaceType>()> *)
{
    auto element = space.begin() ;
    const auto element_end = space.end() ;

    set<Index> dofs_set ;

    for (; element != element_end ; ++element)
    {
        const vector< Index > &element_dofs = element->get_local_to_global() ;
        for (const Index & dof : element_dofs)
            dofs_set.insert(dof) ;
    }


    vector<Index> space_dofs(dofs_set.begin(), dofs_set.end()) ;

    return (space_dofs) ;
}



void apply_boundary_values(const std::map<Index,Real> &boundary_values,
                           Matrix      &matrix,
                           Vector      &rhs,
                           Vector      &solution)
{
    /*
    LogStream out ;

    out <<"Before" <<std::endl;
    out <<"RHS" <<std::endl;
    rhs.print(out) ;
    out <<"Matrix" <<std::endl;
    matrix.print(out) ;
    //*/
    std::vector<Index> constrained_rows;

    auto dof = boundary_values.begin() ;
    const auto dof_end = boundary_values.end() ;
    for (; dof != dof_end ; ++dof)
    {
        Index row_id = dof->first;
        const Real bc_value  = dof->second;
        const Real mat_value = matrix(row_id,row_id);


        // set the matrix in write mode
        matrix.resume_fill();

        // set the selected row to 0.0
        matrix.clear_row(row_id);

        // set the diagonal element corresponding to the entry
        // (row_id,row_id) to mat_value
        matrix.add_entry(row_id, row_id, mat_value);

        // communicate the matrix values to the different processors
        matrix.fill_complete();

        rhs(row_id) = bc_value * mat_value;
        solution(row_id) = bc_value;
    }

}


} ;
IGA_NAMESPACE_CLOSE

#include <igatools/linear_algebra/dof_tools.inst>
