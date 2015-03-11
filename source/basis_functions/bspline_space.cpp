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

#include <igatools/basis_functions/bspline_space.h>
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
BSplineSpace(const int degree,
             std::shared_ptr<GridType> knots,
             const InteriorReg interior_reg,
             const bool periodic,
             const BasisEndBehaviour end_b)
    :
    BSplineSpace(Degrees(degree), knots, interior_reg,
                 Periodicity(filled_array<bool,dim>(periodic)),
                 EndBehaviour(filled_array<BasisEndBehaviour,dim>(end_b)))
{}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
create(const int degree,
       std::shared_ptr<GridType> knots,
       const InteriorReg interior_reg,
       const bool periodic,
       const BasisEndBehaviour end_b) -> shared_ptr<self_t>
{
    auto sp = shared_ptr<self_t>(new self_t(degree,knots,interior_reg,periodic,end_b));
    Assert(sp != nullptr, ExcNullPtr());

    sp->create_connection_for_h_refinement(sp);

    return sp;
}



template<int dim_, int range_, int rank_>
BSplineSpace<dim_, range_, rank_>::
BSplineSpace(const Degrees &deg,
             std::shared_ptr<GridType> knots,
             const InteriorReg interior_reg,
             const Periodicity &periodic,
             const EndBehaviour &end_b)
    :
    BSplineSpace(DegreeTable(true,deg),
                 knots,
                 SpaceData::get_multiplicity_from_regularity(interior_reg,DegreeTable(true,deg),
                                                             knots->get_num_intervals()),
                 PeriodicTable(true,periodic),
                 EndBehaviourTable(true,end_b))
{}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
create(const Degrees &deg,
       std::shared_ptr<GridType> knots,
       const InteriorReg interior_reg,
       const Periodicity &periodic,
       const EndBehaviour &end_b)
-> shared_ptr<self_t>
{
    auto sp = shared_ptr<self_t>(new self_t(deg, knots, interior_reg, periodic, end_b));
    Assert(sp != nullptr, ExcNullPtr());

    sp->create_connection_for_h_refinement(sp);

    return sp;
}



template<int dim_, int range_, int rank_>
BSplineSpace<dim_, range_, rank_>::
BSplineSpace(const DegreeTable &deg,
             std::shared_ptr<GridType> knots,
             std::shared_ptr<const MultiplicityTable> interior_mult,
             const PeriodicTable &periodic,
             const EndBehaviourTable &end_b)
    :
    BaseSpace(SpaceData::create(deg, knots, interior_mult, periodic)),
    end_b_(end_b),
    dof_distribution_global_(
        this->space_data_->get_num_basis_table(),
        this->space_data_->get_degree(),
        this->space_data_->get_periodic_table()),
    operators_(
        this->space_data_->get_grid(),
        this->space_data_->compute_knots_with_repetition(end_b),
        this->space_data_->accumulated_interior_multiplicities(),
        this->space_data_->get_degree()),
    end_interval_(end_b.get_comp_map())
{
    //------------------------------------------------------------------------------
// TODO (pauletti, Dec 24, 2014): after it work it should be recoded properly
    const auto rep_knots =
        this->space_data_->compute_knots_with_repetition(end_b_);
    const auto &degt = this->get_degree();
    for (auto i : end_interval_.get_active_components_id())
        for (int dir=0; dir<dim; ++dir)
        {
            const auto p = deg[i][dir];

            const auto x1 = knots->get_knot_coordinates().get_data_direction(dir)[1];
            const auto a = knots->get_knot_coordinates().get_data_direction(dir)[0];
            const auto x0 = rep_knots[i].get_data_direction(dir)[p];
            end_interval_[i][dir].first = (x1-a) / (x1-x0);

            const auto xk= *(knots->get_knot_coordinates().get_data_direction(dir).end()-2);
            const auto b = *(knots->get_knot_coordinates().get_data_direction(dir).end()-1);
            const auto xk1 = *(rep_knots[i].get_data_direction(dir).end() - (p+1));
            end_interval_[i][dir].second = (b-xk) / (xk1-xk);
        }
    //------------------------------------------------------------------------------



    //------------------------------------------------------------------------------
    dof_distribution_global_.add_dofs_property(property_active_);
    dof_distribution_global_.set_all_dofs_property_status(property_active_,true);
    //------------------------------------------------------------------------------
}




template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
create(const DegreeTable &deg,
       std::shared_ptr<GridType> knots,
       std::shared_ptr<const MultiplicityTable> interior_mult,
       const PeriodicTable &periodic,
       const EndBehaviourTable &end_b)
-> shared_ptr<self_t>
{
    auto sp = shared_ptr<self_t>(new self_t(deg, knots, interior_mult, periodic, end_b));
    Assert(sp != nullptr, ExcNullPtr());

    sp->create_connection_for_h_refinement(sp);

    return sp;
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
    auto elem = shared_ptr<Elem>(new Elem(
        this->shared_from_this(),
        flat_index));
    Assert(elem != nullptr, ExcNullPtr());

    return elem;
}


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
    auto sub_mult   = this->space_data_->template get_sub_space_mult<k>(s_id);
    auto sub_degree = this->space_data_->template get_sub_space_degree<k>(s_id);
    auto sub_periodic = this->space_data_->template get_sub_space_periodicity<k>(s_id);

    using SubRefSp = BSplineSpace<k,range,rank>;

    using SubEndBT = typename SubRefSp::EndBehaviourTable;
    auto &k_elem = UnitElement<dim>::template get_elem<k>(s_id);
    const auto &active_dirs = k_elem.active_directions;

    SubEndBT sub_end_b(end_b_.get_comp_map());
    for (int comp : end_b_.get_active_components_id())
        for (int j=0; j<k; ++j)
            sub_end_b[comp][j] = end_b_[comp][active_dirs[j]];
    auto sub_space =
    SubRefSp::create(sub_degree, sub_grid, sub_mult, sub_periodic, sub_end_b);

    // Creating the mapping between the space degrees of freedom
    const int n_dir = k_elem.constant_directions.size();

    TensorIndex<dim> tensor_index;
    int comp_i = 0;
    dof_map.resize(sub_space->get_num_basis());
    for (auto comp : SpaceData::components)
    {
        const int n_basis = sub_space->get_num_basis(comp);
        const auto &sub_local_indices = sub_space->get_dof_distribution_global().get_index_table()[comp];
        const auto &elem_global_indices = dof_distribution_global_.get_index_table()[comp];

        for (Index sub_i = 0; sub_i < n_basis; ++sub_i, ++comp_i)
        {
            const auto sub_base_id = sub_local_indices.flat_to_tensor(sub_i);

            for (int j=0; j<k; ++j)
                tensor_index[active_dirs[j]] = sub_base_id[j];
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
                                   this->space_data_->get_num_basis_table(),
                                   this->space_data_->get_degree(),
                                   this->space_data_->get_periodic_table());
    operators_ = BernsteinExtraction<dim, range, rank>(
                     this->get_grid(),
                     this->space_data_->compute_knots_with_repetition(end_b_),
                     this->space_data_->accumulated_interior_multiplicities(),
                     this->get_degree());
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
void
BSplineSpace<dim_, range_, rank_>::
create_connection_for_h_refinement(std::shared_ptr<self_t> space)
{
    using SlotType = typename CartesianGrid<dim>::SignalRefineSlot;

    auto refinement_func_bspline_space =
        std::bind(&self_t::refine_h_after_grid_refinement,
                  space.get(),
                  std::placeholders::_1,
                  std::placeholders::_2);
    this->connect_refinement_h_function(
        SlotType(refinement_func_bspline_space).track_foreign(space));
}


template<int dim_, int range_, int rank_>
void
BSplineSpace<dim_, range_, rank_>::
print_info(LogStream &out) const
{
    out.begin_item("Spline Space:");
    this->space_data_->print_info(out);
    out.end_item();


    out.begin_item("DoFs Distribution:");
    dof_distribution_global_.print_info(out);
    out.end_item();

    out.begin_item("Bernstein Extraction:");
    operators_.print_info(out);
    out.end_item();
}

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/bspline_space.inst>
