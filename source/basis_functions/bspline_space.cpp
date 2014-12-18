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

#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/space_manager.h>
#include <igatools/base/sub_function.h>
#include <igatools/base/identity_function.h>

using std::endl;
using std::array;

using std::shared_ptr;
using std::make_shared;
using std::const_pointer_cast;

IGA_NAMESPACE_OPEN

template<int dim_, int range_, int rank_>
BSplineSpace<dim_, range_, rank_>::
BSplineSpace(const int degree, shared_ptr<GridType> grid)
    :
    BSplineSpace(TensorIndex<dim>(degree), grid)
{}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
create(const int degree, shared_ptr< GridType > knots) -> shared_ptr<self_t>
{
    return shared_ptr<self_t>(new self_t(degree, knots));
}



template<int dim_, int range_, int rank_>
BSplineSpace<dim_, range_, rank_>::
BSplineSpace(const TensorIndex<dim> &degree, shared_ptr<GridType> knots)
    :
    BSplineSpace(DegreeTable(degree), knots, true)
{}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
create(const TensorIndex<dim> &degree, shared_ptr<GridType> knots)
-> shared_ptr<self_t>
{
    return shared_ptr<self_t>(new self_t(degree, knots));
}



template<int dim_, int range_, int rank_>
BSplineSpace<dim_, range_, rank_>::
BSplineSpace(const DegreeTable &deg,
             std::shared_ptr<GridType> knots,
             const bool homogeneous_range)
    :
    BaseSpace(deg, knots, BaseSpace::InteriorReg::maximum),
    dof_distribution_global_(knots,BaseSpace::accumulated_interior_multiplicities(),
                             BaseSpace::get_num_basis_table(),BaseSpace::get_degree()),
    dof_distribution_patch_(knots,BaseSpace::accumulated_interior_multiplicities(),
                            BaseSpace::get_num_basis_table(),BaseSpace::get_degree()),
    operators_(knots,
               BaseSpace::compute_knots_with_repetition(this->get_end_behaviour()),
               BaseSpace::accumulated_interior_multiplicities(), deg)
{
    // create a signal and a connection for the grid refinement
    this->connect_refinement_h_function(
        std::bind(&self_t::refine_h_after_grid_refinement, this,
                  std::placeholders::_1,std::placeholders::_2));
}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
create(const DegreeTable &deg,
       std::shared_ptr<GridType> knots,
       const bool homogeneous_range) -> shared_ptr<self_t>
{
    return shared_ptr<self_t>(new self_t(deg, knots, homogeneous_range));
}



template<int dim_, int range_, int rank_>
BSplineSpace<dim_, range_, rank_>::
BSplineSpace(const DegreeTable &deg,
             std::shared_ptr<GridType> knots,
             std::shared_ptr<const MultiplicityTable> interior_mult,
             const EndBehaviourTable &end_b)
    :
    BaseSpace(deg, knots, interior_mult, end_b),
    dof_distribution_global_(knots,BaseSpace::accumulated_interior_multiplicities(),
                             BaseSpace::get_num_basis_table(),BaseSpace::get_degree()),
    dof_distribution_patch_(knots,BaseSpace::accumulated_interior_multiplicities(),
                            BaseSpace::get_num_basis_table(),BaseSpace::get_degree()),
    operators_(knots,
               BaseSpace::compute_knots_with_repetition(this->get_end_behaviour()),
               BaseSpace::accumulated_interior_multiplicities(), deg)
{
    // create a signal and a connection for the grid refinement
    this->connect_refinement_h_function(
        std::bind(&self_t::refine_h_after_grid_refinement, this,
                  std::placeholders::_1,std::placeholders::_2));
}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
create(const DegreeTable &deg,
       std::shared_ptr<GridType> knots,
       std::shared_ptr<const MultiplicityTable> interior_mult,
       const EndBehaviourTable &ends)
-> shared_ptr<self_t>
{
    return shared_ptr<self_t>(new self_t(deg, knots, interior_mult, ends));
}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
get_reference_space() const -> shared_ptr<const self_t>
{
    return this->shared_from_this();
}

template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
create_element(const Index flat_index) const -> std::shared_ptr<ReferenceElement<dim_,range_,rank_> >
{
    using Elem = BSplineElement<dim_,range_,rank_>;
    auto elem = shared_ptr<Elem>(new Elem(this->shared_from_this(),flat_index));
    Assert(elem != nullptr, ExcNullPtr());

    return elem;
}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
begin() const -> ElementIterator
{
    return ElementIterator(this->create_element(0));
}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
last() const -> ElementIterator
{
    Assert(false,ExcNotImplemented());
    //TODO: the index of the last element in the grid is not correct, because we need to use the
    // index of the last ACTIVE element in the grid
    return ElementIterator(this->create_element(this->get_grid()->get_num_active_elems()-1));
}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
end() const -> ElementIterator
{
    return ElementIterator(this->create_element(IteratorState::pass_the_end));
}


#if 0
template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
get_element_handler() const -> ElementHandler
{
    return ElementHandler(this->shared_from_this());
}
#endif


template<int dim_, int range_, int rank_>
template<int k>
auto
BSplineSpace<dim_, range_, rank_>::
get_ref_sub_space(const int s_id,
                  InterSpaceMap<k> &dof_map,
                  std::shared_ptr<CartesianGrid<k>> sub_grid) const
-> std::shared_ptr<SubRefSpace<k> >
{
    if (!(sub_grid))
    {
        typename GridType::template InterGridMap<k>  elem_map;
        sub_grid   = this->get_grid()->template get_sub_grid<k>(s_id, elem_map);
    }
    auto sub_mult   = this->template get_sub_space_mult<k>(s_id);
    auto sub_degree = this->template get_sub_space_degree<k>(s_id);
    auto end_b      = this->template get_sub_space_end_b<k>(s_id);

    auto sub_space = BSplineSpace<k,range,rank>::create(sub_degree, sub_grid, sub_mult,end_b);

    auto &k_elem = UnitElement<dim>::template get_elem<k>(s_id);

    // Crating the mapping between the space degrees of freedom
    const auto &active_dirs = k_elem.active_directions;
    const int n_dir = k_elem.constant_directions.size();

    TensorIndex<dim> tensor_index;
    int comp_i = 0;
    dof_map.resize(sub_space->get_num_basis());
    for (auto comp : components)
    {
        const int n_basis = sub_space->get_num_basis(comp);
        const auto &sub_local_indices = sub_space->get_dof_distribution_patch().get_index_table()[comp];
        const auto &elem_global_indices = dof_distribution_global_.get_index_table()[comp];

        for (Index sub_i = 0; sub_i < n_basis; ++sub_i, ++comp_i)
        {
            const auto sub_base_id = sub_local_indices.flat_to_tensor(sub_i);

            for (int j=0; j<k; ++j)
                tensor_index[active_dirs[j]] =  sub_base_id[j];
            for (int j=0; j<n_dir; ++j)
            {
                auto dir = k_elem.constant_directions[j];
                auto val = k_elem.constant_values[j];
                const int fixed_id = val * (this->get_num_basis(comp, dir) - 1);
                tensor_index[dir] = fixed_id;

            }
            dof_map[comp_i] = elem_global_indices(tensor_index);
        }

    }

    return sub_space;
}



template<int dim_, int range_, int rank_>
template<int k>
auto
BSplineSpace<dim_, range_, rank_>::
get_sub_space(const int s_id, InterSpaceMap<k> &dof_map,
              std::shared_ptr<CartesianGrid<k>> sub_grid,
              std::shared_ptr<typename GridType::template InterGridMap<k>> elem_map) const
-> std::shared_ptr<SubSpace<k> >
{
    using SubMap = SubMapFunction<k, dim, space_dim>;
    auto grid =  this->get_grid();
//    typename GridType::template InterGridMap<k> elem_map;
//    auto sub_grid = this->get_grid()->template get_sub_grid<k>(s_id, elem_map);

    auto sub_ref_space = get_ref_sub_space(s_id, dof_map, sub_grid);
    auto F = IdentityFunction<dim>::create(grid);
    auto sub_map_func = SubMap::create(sub_grid, F, s_id, *elem_map);
    auto sub_space = SubSpace<k>::create(sub_ref_space, sub_map_func);
    return sub_space;
}





template<int dim_, int range_, int rank_>
void
BSplineSpace<dim_, range_, rank_>::
refine_h_after_grid_refinement(
    const std::array<bool,dim> &refinement_directions,
    const GridType &grid_old)
{
    dof_distribution_global_ = DofDistribution<dim, range, rank>(
                                   this->get_grid(),
                                   BaseSpace::accumulated_interior_multiplicities(),
                                   BaseSpace::get_num_basis_table(),
                                   BaseSpace::get_degree());

    dof_distribution_patch_ = DofDistribution<dim, range, rank>(
                                  this->get_grid(),
                                  BaseSpace::accumulated_interior_multiplicities(),
                                  BaseSpace::get_num_basis_table(),
                                  BaseSpace::get_degree());

    operators_ = BernsteinExtraction<dim, range, rank>(
                     this->get_grid(),
                     BaseSpace::compute_knots_with_repetition(this->get_end_behaviour()),
                     BaseSpace::accumulated_interior_multiplicities(),
                     this->get_degree());
}





template<int dim_, int range_, int rank_>
void
BSplineSpace<dim_, range_, rank_>::
add_dofs_offset(const Index offset)
{
    dof_distribution_global_.add_dofs_offset(offset);
}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
get_dof_distribution_global() const -> const DofDistribution<dim, range, rank> &
{
    return dof_distribution_global_;
}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
get_dof_distribution_global() -> DofDistribution<dim, range, rank> &
{
    return dof_distribution_global_;
}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
get_dof_distribution_patch() const -> const DofDistribution<dim, range, rank> &
{
    return dof_distribution_patch_;
}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
get_dof_distribution_patch() -> DofDistribution<dim, range, rank> &
{
    return dof_distribution_patch_;
}


template<int dim_, int range_, int rank_>
Index
BSplineSpace<dim_, range_, rank_>::
get_global_dof_id(const TensorIndex<dim> &tensor_index,
                  const Index comp) const
{
    return dof_distribution_global_.get_index_table()[comp](tensor_index);
}


template<int dim_, int range_, int rank_>
vector<Index>
BSplineSpace<dim_, range_, rank_>::
get_loc_to_global(const CartesianGridElement<dim> &element) const
{
    return dof_distribution_global_.get_loc_to_global_indices(element);
}



template<int dim_, int range_, int rank_>
vector<Index>
BSplineSpace<dim_, range_, rank_>::
get_loc_to_patch(const CartesianGridElement<dim> &element) const
{
    return dof_distribution_patch_.get_loc_to_global_indices(element);
}


template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
get_space_manager() -> shared_ptr<SpaceManager>
{
    auto space_manager = make_shared<SpaceManager>(SpaceManager());

    auto this_space = this->shared_from_this();

    space_manager->spaces_insertion_open();
    space_manager->add_space(this_space);
    space_manager->spaces_insertion_close();


    space_manager->spaces_connectivity_open();
    space_manager->add_spaces_connection(this_space);
    space_manager->spaces_connectivity_close();

    return space_manager;
}

template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
get_space_manager() const -> std::shared_ptr<const SpaceManager>
{
    return const_cast<self_t &>(*this).get_space_manager();
}

template<int dim_, int range_, int rank_>
void
BSplineSpace<dim_, range_, rank_>::
print_info(LogStream &out) const
{
    out.begin_item("Spline Space:");
    BaseSpace::print_info(out);
    out.end_item();

    out.begin_item("Patch Basis Indices:");
    dof_distribution_patch_.print_info(out);
    out.end_item();

    out.begin_item("Global Basis Indices:");
    dof_distribution_global_.print_info(out);
    out.end_item();

    out.begin_item("Bernstein Extraction:");
    operators_.print_info(out);
    out.end_item();
}

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/bspline_space.inst>
