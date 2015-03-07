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
#include <igatools/basis_functions/space_manager.h>
#include <igatools/base/sub_function.h>
#include <igatools/base/identity_function.h>
#include <igatools/utils/multi_array_utils.h>

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
        this->space_data_->get_grid(),
        this->space_data_->accumulated_interior_multiplicities(),
        this->space_data_->get_num_basis_table(),
        this->space_data_->get_degree(),
        this->space_data_->get_periodic_table()),
#if 0
    dof_distribution_patch_(
        this->space_data_->get_grid(),
        this->space_data_->accumulated_interior_multiplicities(),
        this->space_data_->get_num_basis_table(),
        this->space_data_->get_degree(),
        this->space_data_->get_periodic_table()),
#endif
    operators_(
        this->space_data_->get_grid(),
        this->space_data_->compute_knots_with_repetition(end_b),
        this->space_data_->accumulated_interior_multiplicities(),
        this->space_data_->get_degree()),
    end_interval_(end_b.get_comp_map())
{
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
    auto elem = shared_ptr<Elem>(new Elem(this->shared_from_this(),flat_index));
    Assert(elem != nullptr, ExcNullPtr());

    return elem;
}


template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
begin() const -> ElementIterator
{
    return ElementIterator(this->create_element(0),GridType::elems_property_none);
}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
last() const -> ElementIterator
{
    Assert(false,ExcNotImplemented());
    //TODO: the index of the last element in the grid is not correct, because we need to use the
    // index of the last ACTIVE element in the grid
    return ElementIterator(this->create_element(this->get_grid()->get_num_all_elems()-1),GridType::elems_property_none);
}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
end() const -> ElementIterator
{
    return ElementIterator(this->create_element(IteratorState::pass_the_end),GridType::elems_property_none);
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
                                   this->get_grid(),
                                   this->space_data_->accumulated_interior_multiplicities(),
                                   this->space_data_->get_num_basis_table(),
                                   this->space_data_->get_degree(),
                                   this->space_data_->get_periodic_table());
#if 0
    dof_distribution_patch_ = DofDistribution<dim, range, rank>(
                                  this->get_grid(),
                                  this->space_data_->accumulated_interior_multiplicities(),
                                  this->space_data_->get_num_basis_table(),
                                  this->space_data_->get_degree(),
                                  this->space_data_->get_periodic_table());
#endif
    operators_ = BernsteinExtraction<dim, range, rank>(
                     this->get_grid(),
                     this->space_data_->compute_knots_with_repetition(end_b_),
                     this->space_data_->accumulated_interior_multiplicities(),
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


#if 0
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
#endif


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
get_element_dofs(
    const CartesianGridElement<dim> &element,
    const DofDistribution<dim, range, rank> &dofs_distribution) const
{
    const auto &accum_mult = this->space_data_->accumulated_interior_multiplicities();
    const auto &index_table = dofs_distribution.get_index_table();

    const auto &degree_table = this->space_data_->get_degree();

    vector<Index> element_dofs;
    const auto &elem_tensor_id = element.get_tensor_index();

    using Topology = UnitElement<dim>;

    for (int comp = 0 ; comp < SpaceData::n_components ; ++comp)
    {
        //-----------------------------------------------------------------
        // building the lookup table for the local dof id on the current component of the element --- begin
        // TODO (MM, March 06, 2015): this can be put on the SplineSpace constructor for optimization
        const auto &degree_comp = degree_table[comp];

        TensorSize<dim> dofs_t_size_elem_comp;
        for (const auto dir : Topology::active_directions)
            dofs_t_size_elem_comp[dir] = degree_comp[dir] + 1;

        const auto dofs_f_size_elem_comp = dofs_t_size_elem_comp.flat_size();

        vector<Index> elem_comp_dof_f_id(dofs_f_size_elem_comp);
        std::iota(elem_comp_dof_f_id.begin(),elem_comp_dof_f_id.end(),0);

        vector<TensorIndex<dim>> elem_comp_dof_t_id;
        const auto w_dofs_elem_comp = MultiArrayUtils<dim>::compute_weight(dofs_t_size_elem_comp);
        for (const auto dof_f_id : elem_comp_dof_f_id)
            elem_comp_dof_t_id.emplace_back(MultiArrayUtils<dim>::flat_to_tensor_index(dof_f_id,w_dofs_elem_comp));
        // building the lookup table for the local dof id on the current component of the element --- end
        //-----------------------------------------------------------------



        //-----------------------------------------------------------------
        const auto &index_table_comp = index_table[comp];

        const auto dof_t_origin = accum_mult[comp].cartesian_product(elem_tensor_id);
        for (const auto loc_dof_t_id : elem_comp_dof_t_id)
        {
            const auto dof_t_id = dof_t_origin + loc_dof_t_id;
            element_dofs.emplace_back(index_table_comp(dof_t_id));
        }
        //-----------------------------------------------------------------

    } // end comp loop

    return element_dofs;
}

template<int dim_, int range_, int rank_>
vector<Index>
BSplineSpace<dim_, range_, rank_>::
get_loc_to_global(const CartesianGridElement<dim> &element) const
{
    return this->get_element_dofs(element,dof_distribution_global_);
}



template<int dim_, int range_, int rank_>
vector<Index>
BSplineSpace<dim_, range_, rank_>::
get_loc_to_patch(const CartesianGridElement<dim> &element) const
{
    const auto elem_dofs_global = this->get_loc_to_global(element);
    vector<Index> elem_dofs_local;
    for (const auto dof_global : elem_dofs_global)
    	elem_dofs_local.push_back(dof_distribution_global_.global_to_patch_local(dof_global));

    return elem_dofs_local;
#if 0
    return this->get_element_dofs(element,dof_distribution_patch_);
#endif
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

//#if 0
    out.begin_item("Patch Basis Indices:");
    dof_distribution_global_.print_info(out);
    out.end_item();
//#endif
    out.begin_item("Global Basis Indices:");
    dof_distribution_global_.print_info(out);
    out.end_item();

    out.begin_item("Bernstein Extraction:");
    operators_.print_info(out);
    out.end_item();
}

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/bspline_space.inst>
