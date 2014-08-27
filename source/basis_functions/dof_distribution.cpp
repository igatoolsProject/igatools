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
using std::map;
using std::shared_ptr;
using std::make_shared;

IGA_NAMESPACE_OPEN

template<int dim, int range, int rank>
DofDistribution<dim, range, rank>::
DofDistribution(shared_ptr<CartesianGrid<dim> > grid,
                const MultiplicityTable &accum_mult,
                const SpaceDimensionTable &n_basis,
                const DegreeTable &degree_table,
                DistributionPolicy pol)
    :
    TensorSizedContainer<dim>(grid->get_num_intervals()),
    elements_loc_to_global_flat_view_(
        make_shared<map<Index,DofsConstView>>(map<Index,DofsConstView>())),
    policy_(pol)
{
    Assert(pol == DistributionPolicy::standard, ExcNotImplemented());

    //-----------------------------------------------------------------------
    // fills the standard distribution, sorted by component and
    // by direction x moves faster
    Index dof_id = 0;
    for (int comp = 0 ; comp < Space::n_components ; ++comp)
    {
        index_table_(comp).resize(n_basis(comp));
        for (auto &x : index_table_(comp))
            x = dof_id++;
    }
    //-----------------------------------------------------------------------



    //-----------------------------------------------------------------------
    // creating the dofs view from the dofs components views -- begin
    vector<DofsComponentView> components_views;
    for (auto &dofs_distribution_comp : index_table_)
        components_views.emplace_back(dofs_distribution_comp.get_flat_view());

    dofs_view_ = DofsView(
                     DofsIterator(components_views,0),
                     DofsIterator(components_views,IteratorState::pass_the_end));
    // creating the dofs view from the dofs components views -- end
    //-----------------------------------------------------------------------


    //-----------------------------------------------------------------------
    SpaceDimensionTable n_elem_basis;
    for (int iComp = 0 ; iComp <  Space::n_components ; ++iComp)
        n_elem_basis(iComp) = TensorSize<dim>(degree_table(iComp)+1);

    this->create_element_loc_to_global_view(grid,accum_mult,n_elem_basis);
    //-----------------------------------------------------------------------
}

template<int dim, int range, int rank>
Index
DofDistribution<dim, range, rank>::
get_min_dof_id() const
{
    return *std::min_element(dofs_view_.cbegin(),dofs_view_.cend());
}

template<int dim, int range, int rank>
Index
DofDistribution<dim, range, rank>::
get_max_dof_id() const
{
    return *std::max_element(dofs_view_.cbegin(),dofs_view_.cend());
}


template<int dim, int range, int rank>
auto
DofDistribution<dim, range, rank>::
basis_flat_to_tensor(const Index index, const Index comp) const -> TensorIndex<dim>
{
    return index_table_(comp).flat_to_tensor(index);
}


template<int dim, int range, int rank>
Index
DofDistribution<dim, range, rank>::
basis_tensor_to_flat(const TensorIndex<dim> &tensor_index,
                     const Index comp) const
{
    return index_table_(comp).tensor_to_flat(tensor_index);
}


template<int dim, int range, int rank>
void
DofDistribution<dim, range, rank>::
create_element_loc_to_global_view(
    std::shared_ptr<const CartesianGrid<dim> > grid,
    const MultiplicityTable &accum_mult,
    const SpaceDimensionTable &n_elem_basis)
{
    for (const auto elem : *grid)
    {
        const auto t_index = elem.get_tensor_index();

        vector<DofsComponentConstView> dofs_elem_ranges;
        for (int comp = 0; comp < Space::n_components; ++comp)
        {
            const auto &index_table_comp = index_table_(comp);

            auto origin = accum_mult(comp).cartesian_product(t_index);
            Index origin_flat_id = index_table_comp.tensor_to_flat(origin);

            auto increment = n_elem_basis(comp);

            using VecIt = vector<Index>::const_iterator;
            const VecIt comp_dofs_begin = index_table_comp.get_data().begin();

            if (dim == 0)
            {
                const VecIt pos_begin = comp_dofs_begin + origin_flat_id;
                const VecIt pos_end   = pos_begin+1; // one dof for spaces with dim==0

                dofs_elem_ranges.emplace_back(DofsComponentConstView(pos_begin,pos_end));

            } // end if (dim == 0)
            else if (dim == 1)
            {
                const VecIt pos_begin = comp_dofs_begin + origin_flat_id;
                const VecIt pos_end   = pos_begin+increment[0];

                dofs_elem_ranges.emplace_back(DofsComponentConstView(pos_begin,pos_end));

            } // end else if (dim == 1)
            else if (dim == 2)
            {
                TensorIndex<dim> incr_t_id;
                for (incr_t_id[1] = 0 ; incr_t_id[1] < increment[1]; ++incr_t_id[1])
                {
                    TensorIndex<dim> pos_t_id = origin + incr_t_id;

                    Index pos_flat_id = index_table_comp.tensor_to_flat(pos_t_id);

                    const VecIt pos_begin = comp_dofs_begin + pos_flat_id;
                    const VecIt pos_end = pos_begin + increment[0];

                    dofs_elem_ranges.emplace_back(DofsComponentConstView(pos_begin,pos_end));
                } // end loop incr_t_id(1)

            } // end else if (dim == 2)
            else if (dim == 3)
            {
                TensorIndex<dim> incr_t_id;
                for (incr_t_id[2] = 0 ; incr_t_id[2] < increment[2]; ++incr_t_id[2])
                {
                    for (incr_t_id[1] = 0 ; incr_t_id[1] < increment[1]; ++incr_t_id[1])
                    {
                        TensorIndex<dim> pos_t_id = origin + incr_t_id;

                        Index pos_flat_id = index_table_comp.tensor_to_flat(pos_t_id);

                        const VecIt pos_begin = comp_dofs_begin + pos_flat_id;
                        const VecIt pos_end = pos_begin + increment[0];

                        dofs_elem_ranges.emplace_back(DofsComponentConstView(pos_begin,pos_end));

                    } // end loop incr_t_id(1)
                } // end loop incr_t_id(2)

            } // end else if (dim == 3)
            else
            {
                Assert(false,ExcNotImplemented());
                AssertThrow(false,ExcNotImplemented());
            }

        } // end loop elem

        const auto f_index = elem.get_flat_index();

        (*elements_loc_to_global_flat_view_)[f_index] =
            DofsConstView(DofsConstIterator(dofs_elem_ranges,0),
                          DofsConstIterator(dofs_elem_ranges,IteratorState::pass_the_end));
    }
}




template<int dim, int range, int rank>
std::vector<Index>
DofDistribution<dim, range, rank>::
get_loc_to_global_indices(const CartesianGridElement<dim> &element) const
{
    const auto &dofs_elem_view = elements_loc_to_global_flat_view_->at(element.get_flat_index());
    return vector<Index>(dofs_elem_view.begin(),dofs_elem_view.end());
}


template<int dim, int range, int rank>
void
DofDistribution<dim, range, rank>::
add_dofs_offset(const Index offset)
{
    for (auto &dofs_component : index_table_)
        for (auto &dof_id : dofs_component)
            dof_id += offset;
}

template<int dim, int range, int rank>
auto
DofDistribution<dim, range, rank>::
get_index_table() const -> const IndexDistributionTable &
{
    return index_table_;
}

template<int dim, int range, int rank>
auto
DofDistribution<dim, range, rank>::
get_elements_view() const -> std::shared_ptr<const std::map<Index,DofsConstView>>
{
    return elements_loc_to_global_flat_view_;
}


template<int dim, int range, int rank>
auto
DofDistribution<dim, range, rank>::
get_dofs_view() -> DofsView &
{
    return dofs_view_;
}

template<int dim, int range, int rank>
auto
DofDistribution<dim, range, rank>::
get_dofs_view() const -> const DofsView &
{
    return dofs_view_;
}


#if 0
template<int dim, int range, int rank>
auto
DofDistribution<dim, range, rank>::
get_num_dofs_element(const Index elem_flat_id) const -> Size
{
    DofsPerElementTable dofs_per_element_table;
    const auto &dofs_element_view = elements_loc_to_global_flat_view_->at(elem_flat_id);

    return dofs_element_view.get_num_entries();
}
#endif

template<int dim, int range, int rank>
auto
DofDistribution<dim, range, rank>::
get_num_dofs_element(const CartesianGridElement<dim> &element) const -> Size
{
//  DofsPerElementTable dofs_per_element_table;
    const auto &dofs_element_view = elements_loc_to_global_flat_view_->at(element.get_flat_index());

    return dofs_element_view.get_num_entries();
}


template<int dim, int range, int rank>
void
DofDistribution<dim, range, rank>::
print_info(LogStream &out) const
{
    for (int comp = 0; comp < Space::n_components; ++comp)
        index_table_(comp).print_info(out);
    out << std::endl;

//    int i = 0;
    for (const auto &dofs_elem : *elements_loc_to_global_flat_view_)
    {
//        out << this->get_loc_to_global_indices(dofs_elem.first) << std::endl;

        out << "[ ";
        for (auto x : dofs_elem.second)
            out << x << " ";
        out << "]" << std::endl;
        //*/
    }
}

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/dof_distribution.inst>

