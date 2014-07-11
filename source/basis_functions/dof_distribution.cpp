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


#include <igatools/basis_functions/dof_distribution.h>

using std::vector;
using std::shared_ptr;

IGA_NAMESPACE_OPEN

template<int dim, int range, int rank>
DofDistribution<dim, range, rank>::
DofDistribution(shared_ptr<CartesianGrid<dim> > grid,
                const MultiplicityTable &accum_mult,
                const SpaceDimensionTable &n_basis,
                const SpaceDimensionTable &n_elem_basis,
                DistributionPolicy pol)
    :
    element_loc_to_global_(grid->get_num_elements_dim()),
    element_loc_to_global_view_(grid->get_num_elements_dim())
{
    Assert(pol == DistributionPolicy::standard, ExcNotImplemented());

    // fills the standard distribution, sorted by component and
    // by direction x moves faster
    for (int comp = 0, j = 0; comp < Space::n_components; ++comp)
    {
        index_distribution_(comp).resize(n_basis(comp));
        for (auto &x : index_distribution_(comp))
            x = j++;
    }
    /*
        for (const auto elem : *grid)
        {
            const auto index = elem.get_tensor_index();
            auto &basis_list = element_loc_to_global_(index);
            auto basis = basis_list.begin();

            for (int comp = 0; comp < Space::n_components; ++comp)
            {
                auto origin = accum_mult(comp).cartesian_product(index);
                auto increment = n_elem_basis(comp);

                auto comp_dofs = index_distribution_(comp).get_sub_array(origin, increment).get_data();
                element_loc_to_global_(index).insert(basis, comp_dofs.begin(), comp_dofs.end());
                for (auto x : element_loc_to_global_)
                    basis = element_loc_to_global_(index).end();
            }
        }
        //*/
    element_loc_to_global_ =
        this->create_element_loc_to_global_from_index_distribution(grid,accum_mult,n_elem_basis,index_distribution_);
}

template<int dim, int range, int rank>
DynamicMultiArray<vector<Index>, dim>
DofDistribution<dim, range, rank>::
create_element_loc_to_global_from_index_distribution(
    shared_ptr<const CartesianGrid<dim> > grid,
    const MultiplicityTable &accum_mult,
    const SpaceDimensionTable &n_elem_basis,
    const IndexDistributionTable &index_distribution) const
{
    DynamicMultiArray<vector<Index>, dim> element_loc_to_global(grid->get_num_elements_dim());
    for (const auto elem : *grid)
    {
        const auto index = elem.get_tensor_index();
        auto &basis_list = element_loc_to_global(index);
        auto basis = basis_list.begin();

        for (int comp = 0; comp < Space::n_components; ++comp)
        {
            auto origin = accum_mult(comp).cartesian_product(index);
            auto increment = n_elem_basis(comp);

            auto comp_dofs = index_distribution(comp).get_sub_array(origin, increment).get_data();
            element_loc_to_global(index).insert(basis, comp_dofs.begin(), comp_dofs.end());
            for (auto x : element_loc_to_global)
                basis = element_loc_to_global(index).end();
        }
    }
    return element_loc_to_global;
}


// TODO (pauletti, May 28, 2014): inline this
template<int dim, int range, int rank>
const vector<Index> &
DofDistribution<dim, range, rank>::
get_loc_to_global_indices(const TensorIndex<dim> &j) const
{
    return element_loc_to_global_(j);
}



// TODO (antolin, Jul 11, 2014): inline this
template<int dim, int range, int rank>
const std::vector<Index> &
DofDistribution<dim, range, rank>::
get_loc_to_global_indices(const Index &j) const
{
    return element_loc_to_global_(j);
}



template<int dim, int range, int rank>
void
DofDistribution<dim, range, rank>::
add_dofs_offset(const Index offset)
{
    for (auto &dofs_component : index_distribution_)
        for (auto &dof_id : dofs_component)
            dof_id += offset;

    for (auto &dofs_element : element_loc_to_global_)
        for (auto &dof_id : dofs_element)
            dof_id += offset;
}


template<int dim, int range, int rank>
void
DofDistribution<dim, range, rank>::
print_info(LogStream &out) const
{
    for (int comp = 0; comp < Space::n_components; ++comp)
        index_distribution_(comp).print_info(out);
    out << std::endl;

    for (auto x : element_loc_to_global_)
        out << x <<  std::endl;

}

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/dof_distribution.inst>

