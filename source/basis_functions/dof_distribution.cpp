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
//    element_loc_to_global_(grid->get_num_elements_dim()),
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
#if 0
    element_loc_to_global_ =
        this->create_element_loc_to_global_from_index_distribution(grid,accum_mult,n_elem_basis,index_distribution_);
#endif

    for (const auto elem : *grid)
    {
        const auto index = elem.get_tensor_index();
        auto &dofs_elem_view = element_loc_to_global_view_(index);

        vector<DofsComponentConstView> dofs_elem_ranges;


        for (int comp = 0; comp < Space::n_components; ++comp)
        {
            const auto &index_distribution_comp = index_distribution_(comp);

            auto origin = accum_mult(comp).cartesian_product(index);
            Index origin_flat_id = index_distribution_comp.tensor_to_flat(origin);

            auto increment = n_elem_basis(comp);

            using VecIt = vector<Index>::const_iterator;
            const VecIt comp_dofs_begin = index_distribution_comp.get_data().begin();

            if (dim == 0)
            {}
            else if (dim == 1)
            {
                const VecIt pos_begin = index_distribution_comp.get_data().begin() + origin_flat_id;
                const VecIt pos_end   = pos_begin+increment[0];

                dofs_elem_ranges.emplace_back(DofsComponentConstView(pos_begin,pos_end));
            }
            else if (dim == 2)
            {
                TensorIndex<dim> incr_t_id;
                for (incr_t_id(1) = 0 ; incr_t_id(1) < increment(1); ++incr_t_id(1))
                {
                    TensorIndex<dim> pos_t_id = origin + incr_t_id;

                    Index pos_flat_id = index_distribution_comp.tensor_to_flat(pos_t_id);

                    const VecIt pos_begin = comp_dofs_begin + pos_flat_id;
                    const VecIt pos_end = pos_begin + increment(0);

                    dofs_elem_ranges.emplace_back(DofsComponentConstView(pos_begin,pos_end));
                } // end loop incr_t_id(1)
            }
            else if (dim == 3)
            {
                TensorIndex<dim> incr_t_id;
                for (incr_t_id(2) = 0 ; incr_t_id(2) < increment(2); ++incr_t_id(2))
                {
                    for (incr_t_id(1) = 0 ; incr_t_id(1) < increment(1); ++incr_t_id(1))
                    {
                        TensorIndex<dim> pos_t_id = origin + incr_t_id;

                        Index pos_flat_id = index_distribution_comp.tensor_to_flat(pos_t_id);

                        const VecIt pos_begin = comp_dofs_begin + pos_flat_id;
                        const VecIt pos_end = pos_begin + increment(0);

                        dofs_elem_ranges.emplace_back(DofsComponentConstView(pos_begin,pos_end));

                    } // end loop incr_t_id(1)
                } // end loop incr_t_id(2)

            }
            else
            {
                Assert(false,ExcNotImplemented());
                AssertThrow(false,ExcNotImplemented());
            }

        }
        //*/

        if (dim != 0)
        {
            dofs_elem_view = DofsView(
                                 DofsConstIterator(dofs_elem_ranges,0),
                                 DofsConstIterator(dofs_elem_ranges,IteratorState::pass_the_end));
        }
//*/
    }
}

#if 0
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
            basis = element_loc_to_global(index).end();
        }
    }
    return element_loc_to_global;
}
#endif

// TODO (pauletti, May 28, 2014): inline this
template<int dim, int range, int rank>
vector<Index>
DofDistribution<dim, range, rank>::
get_loc_to_global_indices(const TensorIndex<dim> &j) const
{
    const auto &dofs_elem_view = element_loc_to_global_view_(j);

    return vector<Index>(dofs_elem_view.begin(),dofs_elem_view.end());
}



// TODO (antolin, Jul 11, 2014): inline this
template<int dim, int range, int rank>
std::vector<Index>
DofDistribution<dim, range, rank>::
get_loc_to_global_indices(const Index &j) const
{
    const auto &dofs_elem_view = element_loc_to_global_view_(j);

    return vector<Index>(dofs_elem_view.begin(),dofs_elem_view.end());
}



template<int dim, int range, int rank>
void
DofDistribution<dim, range, rank>::
add_dofs_offset(const Index offset)
{
    for (auto &dofs_component : index_distribution_)
        for (auto &dof_id : dofs_component)
            dof_id += offset;
#if 0
    for (auto &dofs_element : element_loc_to_global_)
        for (auto &dof_id : dofs_element)
            dof_id += offset;
#endif
}


template<int dim, int range, int rank>
void
DofDistribution<dim, range, rank>::
print_info(LogStream &out) const
{
    for (int comp = 0; comp < Space::n_components; ++comp)
        index_distribution_(comp).print_info(out);
    out << std::endl;

#if 0
    for (auto x : element_loc_to_global_)
        out << x <<  std::endl;
#endif

    int i = 0;
    for (auto dofs_elem : element_loc_to_global_view_)
    {
        out << this->get_loc_to_global_indices(i++) << std::endl;
        /*
        out << "[ ";
        for (auto x : dofs_elem)
            out << x << " ";
        out << "]" << std::endl;
        //*/
    }
}

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/dof_distribution.inst>

