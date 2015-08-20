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
#include <igatools/basis_functions/bspline_element_handler.h>
#include <igatools/functions/sub_function.h>
#include <igatools/functions/identity_function.h>

using std::endl;

using std::shared_ptr;
using std::make_shared;
using std::const_pointer_cast;

IGA_NAMESPACE_OPEN

template<int dim_, int range_, int rank_>
BSplineSpace<dim_, range_, rank_>::
BSplineSpace(const int degree,
             const std::shared_ptr<GridType> &grid,
             const InteriorReg interior_reg,
             const bool periodic,
             const BasisEndBehaviour end_b)
    :
    BSplineSpace(Degrees(degree), grid, interior_reg,
                Periodicity(periodic),
                EndBehaviour(end_b))
{}

template<int dim_, int range_, int rank_>
BSplineSpace<dim_, range_, rank_>::
BSplineSpace(const int degree,
             const std::shared_ptr<const GridType> &grid,
             const InteriorReg interior_reg,
             const bool periodic,
             const BasisEndBehaviour end_b)
    :
    BSplineSpace(Degrees(degree), grid, interior_reg,
                Periodicity(periodic),
                EndBehaviour(end_b))
{}


template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
create_nonconst(const int degree,
                const std::shared_ptr<GridType> &grid,
                const InteriorReg interior_reg,
                const bool periodic,
                const BasisEndBehaviour end_b) -> shared_ptr<self_t>
{
    auto sp = shared_ptr<self_t>(new self_t(degree,grid,interior_reg,periodic,end_b));
    Assert(sp != nullptr, ExcNullPtr());

#ifdef MESH_REFINEMENT
    sp->create_connection_for_insert_knots(sp);
#endif

    return sp;
}

template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
create(const int degree,
       const std::shared_ptr<const GridType> &grid,
       const InteriorReg interior_reg,
       const bool periodic,
       const BasisEndBehaviour end_b) -> shared_ptr<const self_t>
{
    shared_ptr<self_t> sp(new self_t(degree,grid,interior_reg,periodic,end_b));
    Assert(sp != nullptr, ExcNullPtr());

    return sp;
}


template<int dim_, int range_, int rank_>
BSplineSpace<dim_, range_, rank_>::
BSplineSpace(const Degrees &deg,
             const std::shared_ptr<GridType> &grid,
             const InteriorReg interior_reg,
             const Periodicity &periodic,
             const EndBehaviour &end_b)
    :
    BSplineSpace(DegreeTable(true,deg),
                grid,
                SpaceData::get_multiplicity_from_regularity(interior_reg,DegreeTable(true,deg),
                                                            grid->get_num_intervals()),
                PeriodicityTable(true,periodic),
                EndBehaviourTable(true,end_b))
{}

template<int dim_, int range_, int rank_>
BSplineSpace<dim_, range_, rank_>::
BSplineSpace(const Degrees &deg,
             const std::shared_ptr<const GridType> &grid,
             const InteriorReg interior_reg,
             const Periodicity &periodic,
             const EndBehaviour &end_b)
    :
    BSplineSpace(DegreeTable(true,deg),
                grid,
                SpaceData::get_multiplicity_from_regularity(interior_reg,DegreeTable(true,deg),
                                                            grid->get_num_intervals()),
                PeriodicityTable(true,periodic),
                EndBehaviourTable(true,end_b))
{}


template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
create_nonconst(const Degrees &deg,
                const std::shared_ptr<GridType> &grid,
                const InteriorReg interior_reg,
                const Periodicity &periodic,
                const EndBehaviour &end_b)
-> shared_ptr<self_t>
{
    auto sp = shared_ptr<self_t>(new self_t(deg, grid, interior_reg, periodic, end_b));
    Assert(sp != nullptr, ExcNullPtr());

#ifdef MESH_REFINEMENT
    sp->create_connection_for_insert_knots(sp);
#endif

    return sp;
}

template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
create(const Degrees &deg,
       const std::shared_ptr<const GridType> &grid,
       const InteriorReg interior_reg,
       const Periodicity &periodic,
       const EndBehaviour &end_b)
-> shared_ptr<const self_t>
{
    auto sp = shared_ptr<self_t>(new self_t(deg, grid, interior_reg, periodic, end_b));
    Assert(sp != nullptr, ExcNullPtr());

    return sp;
}

template<int dim_, int range_, int rank_>
BSplineSpace<dim_, range_, rank_>::
BSplineSpace(const std::shared_ptr<SpaceData> &space_data,
             const EndBehaviourTable &end_b)
    :
    BaseSpace(space_data->get_ptr_grid()),
    space_data_(space_data),
    end_b_(end_b),
    operators_(*space_data_,end_b),
    end_interval_(end_b.get_comp_map()),
    dof_distribution_(std::make_shared<DofDistribution<dim_,range_,rank_>>(
                         space_data->get_num_basis_table(),
                         space_data->get_degree_table(),
                         space_data->get_periodic_table()))
{
//    Assert(space_data_ != nullptr,ExcNullPtr());
    Assert(dof_distribution_ != nullptr,ExcNullPtr());

    //------------------------------------------------------------------------------
// TODO (pauletti, Dec 24, 2014): after it work it should be recoded properly

    const auto &sp_data = *this->space_data_;
    const auto &grid = *sp_data.get_ptr_const_grid();
    const auto &degree_table = sp_data.get_degree_table();
    const auto rep_knots = sp_data.compute_knots_with_repetition(end_b_);

    for (auto i : end_interval_.get_active_components_id())
    {
        for (int dir=0; dir<dim; ++dir)
        {
            const auto p = degree_table[i][dir];

            const auto &knots_coord_dir = grid.get_knot_coordinates().get_data_direction(dir);

            const auto x1 = knots_coord_dir[1];
            const auto a = knots_coord_dir[0];
            const auto x0 = rep_knots[i].get_data_direction(dir)[p];
            end_interval_[i][dir].first = (x1-a) / (x1-x0);

            const auto xk= *(knots_coord_dir.end()-2);
            const auto b = *(knots_coord_dir.end()-1);
            const auto xk1 = *(rep_knots[i].get_data_direction(dir).end() - (p+1));
            end_interval_[i][dir].second = (b-xk) / (xk1-xk);
        } // end loop dir
    } // end loop i
    //------------------------------------------------------------------------------



#if 0
    //------------------------------------------------------------------------------
    this->dof_distribution_->add_dofs_property(this->dofs_property_active_);
    this->dof_distribution_->set_all_dofs_property_status(this->dofs_property_active_,true);
    //------------------------------------------------------------------------------
#endif
}

template<int dim_, int range_, int rank_>
BSplineSpace<dim_, range_, rank_>::
BSplineSpace(const std::shared_ptr<const SpaceData> &space_data,
             const EndBehaviourTable &end_b)
    :
    BaseSpace(space_data->get_ptr_const_grid()),
    space_data_(space_data),
    end_b_(end_b),
    operators_(*space_data_,end_b),
    end_interval_(end_b.get_comp_map()),
    dof_distribution_(std::make_shared<DofDistribution<dim_,range_,rank_>>(
                         space_data->get_num_basis_table(),
                         space_data->get_degree_table(),
                         space_data->get_periodic_table()))
{
//    Assert(space_data_ != nullptr,ExcNullPtr());
    Assert(dof_distribution_ != nullptr,ExcNullPtr());

    //------------------------------------------------------------------------------
// TODO (pauletti, Dec 24, 2014): after it work it should be recoded properly

    const auto &sp_data = *this->space_data_;
    const auto &grid = *sp_data.get_ptr_const_grid();
    const auto &degree_table = sp_data.get_degree_table();
    const auto rep_knots = sp_data.compute_knots_with_repetition(end_b_);

    for (auto i : end_interval_.get_active_components_id())
    {
        for (int dir=0; dir<dim; ++dir)
        {
            const auto p = degree_table[i][dir];

            const auto &knots_coord_dir = grid.get_knot_coordinates().get_data_direction(dir);

            const auto x1 = knots_coord_dir[1];
            const auto a = knots_coord_dir[0];
            const auto x0 = rep_knots[i].get_data_direction(dir)[p];
            end_interval_[i][dir].first = (x1-a) / (x1-x0);

            const auto xk= *(knots_coord_dir.end()-2);
            const auto b = *(knots_coord_dir.end()-1);
            const auto xk1 = *(rep_knots[i].get_data_direction(dir).end() - (p+1));
            end_interval_[i][dir].second = (b-xk) / (xk1-xk);
        } // end loop dir
    } // end loop i
    //------------------------------------------------------------------------------


#if 0
    //------------------------------------------------------------------------------
    this->dof_distribution_->add_dofs_property(this->dofs_property_active_);
    this->dof_distribution_->set_all_dofs_property_status(this->dofs_property_active_,true);
    //------------------------------------------------------------------------------
#endif
}

template<int dim_, int range_, int rank_>
BSplineSpace<dim_, range_, rank_>::
BSplineSpace(const DegreeTable &deg,
             const std::shared_ptr<GridType> &grid,
             const MultiplicityTable &interior_mult,
             const PeriodicityTable &periodic,
             const EndBehaviourTable &end_b)
    :
    BSplineSpace(SpaceData::create(deg, grid, interior_mult, periodic),end_b)
{}

template<int dim_, int range_, int rank_>
BSplineSpace<dim_, range_, rank_>::
BSplineSpace(const DegreeTable &deg,
             const std::shared_ptr<const GridType> &grid,
             const MultiplicityTable &interior_mult,
             const PeriodicityTable &periodic,
             const EndBehaviourTable &end_b)
    :
    BSplineSpace(SpaceData::create(deg, grid, interior_mult, periodic),end_b)
{}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
create_nonconst(const DegreeTable &deg,
                const std::shared_ptr<GridType> &grid,
                const MultiplicityTable &interior_mult,
                const PeriodicityTable &periodic,
                const EndBehaviourTable &end_b)
-> shared_ptr<self_t>
{
    auto sp = shared_ptr<self_t>(new self_t(deg, grid, interior_mult, periodic, end_b));
    Assert(sp != nullptr, ExcNullPtr());

#ifdef MESH_REFINEMENT
    sp->create_connection_for_insert_knots(sp);
#endif

    return sp;
}

template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
create(const DegreeTable &deg,
       const std::shared_ptr<const GridType> &grid,
       const MultiplicityTable &interior_mult,
       const PeriodicityTable &periodic,
       const EndBehaviourTable &end_b)
-> shared_ptr<const self_t>
{
    auto sp = shared_ptr<self_t>(new self_t(deg, grid, interior_mult, periodic, end_b));
    Assert(sp != nullptr, ExcNullPtr());

    return sp;
}


template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
create_nonconst(const std::shared_ptr<SpaceData> &space_data,
                const EndBehaviourTable &end_b)
-> shared_ptr<self_t>
{
    auto sp = shared_ptr<self_t>(new self_t(space_data, end_b));
    Assert(sp != nullptr, ExcNullPtr());

#ifdef MESH_REFINEMENT
    sp->create_connection_for_insert_knots(sp);
#endif

    return sp;
}

template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
create(const std::shared_ptr<const SpaceData> &space_data,
       const EndBehaviourTable &end_b)
-> shared_ptr<const self_t>
{
    auto sp = shared_ptr<self_t>(new self_t(space_data, end_b));
    Assert(sp != nullptr, ExcNullPtr());

    return sp;
}


template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
get_this_space() const -> shared_ptr<const self_t>
{
    auto ref_sp = const_cast<self_t *>(this)->shared_from_this();
    auto bsp_space = std::dynamic_pointer_cast<self_t>(ref_sp);
    Assert(bsp_space != nullptr,ExcNullPtr());

    return bsp_space;
}

template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
create_element(const Index flat_index) const
-> std::shared_ptr<SpaceElement<dim_,0,range_,rank_,Transformation::h_grad> >
{
    using Elem = BSplineElement<dim_,range_,rank_>;

    auto elem = make_shared<Elem>(this->get_this_space(),flat_index);
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
        InterGridMap elem_map;
        sub_grid   = this->get_ptr_const_grid()->template get_sub_grid<k>(s_id, elem_map);
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
    SubRefSp::create_nonconst(sub_degree, sub_grid, sub_mult, sub_periodic, sub_end_b);

    // Creating the mapping between the space degrees of freedom
    const int n_dir = k_elem.constant_directions.size();
#ifndef NDEBUG
    for (int comp : end_b_.get_active_components_id())
        for (int j=0; j<n_dir; ++j)
            Assert(end_b_[comp][k_elem.constant_directions[j]] ==
            BasisEndBehaviour::interpolatory,
            ExcNotImplemented());
#endif


    TensorIndex<dim> tensor_index;
    int comp_i = 0;
    dof_map.resize(sub_space->get_num_basis());
    const auto &sub_space_index_table = sub_space->get_ptr_const_dof_distribution()->get_index_table();
    const auto     &space_index_table = this->get_ptr_const_dof_distribution()->get_index_table();
    for (auto comp : SpaceData::components)
    {
        const auto n_basis = sub_space->get_num_basis(comp);
        const auto &sub_local_indices = sub_space_index_table[comp];
        const auto &elem_global_indices = space_index_table[comp];

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
              InterGridMap &elem_map) const
-> std::shared_ptr<SubSpace<k> >
{
    using SubMap = SubMapFunction<k, dim, space_dim>;
    auto grid = const_pointer_cast<CartesianGrid<dim_> >(this->get_ptr_const_grid());

    auto sub_ref_space = get_ref_sub_space(s_id, dof_map, sub_grid);
    auto F = IdentityFunction<dim>::create(grid);
    auto sub_map_func = SubMap::create(sub_grid, *F, s_id, elem_map);
    auto sub_space = SubSpace<k>::create_nonconst(sub_ref_space, sub_map_func);
    return sub_space;
}








template<int dim_, int range_, int rank_>
void
BSplineSpace<dim_, range_, rank_>::
get_element_dofs(
    const Index element_id,
    SafeSTLVector<Index> &dofs_global,
    SafeSTLVector<Index> &dofs_local_to_patch,
    SafeSTLVector<Index> &dofs_local_to_elem,
    const std::string &dofs_property) const
{
    const auto &sp_data = *space_data_;
    const auto &accum_mult = sp_data.accumulated_interior_multiplicities();

    const auto &dof_distr = *this->dof_distribution_;
    const auto &index_table = dof_distr.get_index_table();

    dofs_global.clear();
    dofs_local_to_patch.clear();
    dofs_local_to_elem.clear();

    const auto &elem_tensor_id = this->get_ptr_const_grid()->flat_to_tensor(element_id);

    Index dof_loc_to_elem = 0;
    for (const auto comp : SpaceData::components)
    {
        const auto &index_table_comp = index_table[comp];

        const auto dof_t_origin = accum_mult[comp].cartesian_product(elem_tensor_id);

        const auto &elem_comp_dof_t_id = sp_data.get_dofs_tensor_id_elem_table()[comp];

//        if (dofs_property == DofProperties::active)
//        {
//            for (const auto loc_dof_t_id : elem_comp_dof_t_id)
//            {
//                const auto dof_global = index_table_comp(dof_t_origin + loc_dof_t_id);
//                dofs_global.emplace_back(dof_global);
//
//                const auto dof_loc_to_patch = this->dof_distribution_->global_to_patch_local(dof_global);
//                dofs_local_to_patch.emplace_back(dof_loc_to_patch);
//
//                dofs_local_to_elem.emplace_back(dof_loc_to_elem);
//
//                ++dof_loc_to_elem;
//            } // end loop loc_dof_t_id
//        }
//        else
        {
            for (const auto loc_dof_t_id : elem_comp_dof_t_id)
            {
                const auto dof_global = index_table_comp(dof_t_origin + loc_dof_t_id);
                if (dof_distr.test_if_dof_has_property(dof_global, dofs_property))
                {
                    dofs_global.emplace_back(dof_global);

                    const auto dof_loc_to_patch = dof_distr.global_to_patch_local(dof_global);
                    dofs_local_to_patch.emplace_back(dof_loc_to_patch);

                    dofs_local_to_elem.emplace_back(dof_loc_to_elem);

                }
                ++dof_loc_to_elem;
            } // end loop loc_dof_t_id
        }

    } // end comp loop
}



#ifdef MESH_REFINEMENT

template<int dim_, int range_, int rank_>
void
BSplineSpace<dim_, range_, rank_>::
create_connection_for_insert_knots(std::shared_ptr<self_t> space)
{
    Assert(space != nullptr, ExcNullPtr());
    Assert(&(*space) == &(*this), ExcMessage("Different objects."));

    auto func_to_connect =
        std::bind(&self_t::rebuild_after_insert_knots,
                  space.get(),
                  std::placeholders::_1,
                  std::placeholders::_2);

    using SlotType = typename CartesianGrid<dim>::SignalInsertKnotsSlot;
    this->get_ptr_grid()->connect_insert_knots(SlotType(func_to_connect).track_foreign(space));
}


template<int dim_, int range_, int rank_>
void
BSplineSpace<dim_, range_, rank_>::
rebuild_after_insert_knots(
    const SafeSTLArray<SafeSTLVector<Real>,dim> &knots_to_insert,
    const CartesianGrid<dim> &old_grid)
{
    this->ref_space_previous_refinement_ =
        BSplineSpace<dim_,range_,rank_>::create(
            this->space_data_->get_spline_space_previous_refinement(),
            this->end_b_);


    this->dof_distribution_ = shared_ptr<DofDistribution<dim_,range_,rank_>>(
                                  new DofDistribution<dim_,range_,rank_>(
                                      this->space_data_->get_num_basis_table(),
                                      this->space_data_->get_degree_table(),
                                      this->space_data_->get_periodic_table()));

    operators_ = BernsteinExtraction<dim,range,rank>(*this->space_data_,end_b_);
    /*
                     *this->get_ptr_grid(),
                     this->space_data_->compute_knots_with_repetition(end_b_),
                     this->space_data_->accumulated_interior_multiplicities(),
                     this->space_data_->get_degree_table());
                     //*/
}
#endif //MESH_REFINEMENT

template<int dim_, int range_, int rank_>
void
BSplineSpace<dim_, range_, rank_>::
print_info(LogStream &out) const
{
    out.begin_item("Spline Space:");
    this->space_data_->print_info(out);
    out.end_item();


    out.begin_item("DoFs Distribution:");
    this->dof_distribution_->print_info(out);
    out.end_item();


    out.begin_item("Bernstein Extraction:");
    operators_.print_info(out);
    out.end_item();
}


template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
get_degree_table() const -> const DegreeTable &
{
    return this->space_data_->get_degree_table();
}

template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
get_periodicity() const -> const PeriodicityTable &
{
    return space_data_->get_periodicity();
}

template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
get_end_behaviour_table() const -> const EndBehaviourTable &
{
    return end_b_;
};


template<int dim_, int range_, int rank_>
bool
BSplineSpace<dim_, range_, rank_>::
is_bspline() const
{
    return true;
}

template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
get_elem_handler() const
-> std::shared_ptr<SpaceElementHandler<dim_,0,range_,rank_,Transformation::h_grad>>
{
    return ElementHandler::create(this->get_this_space());
}


template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
get_ptr_const_dof_distribution() const -> shared_ptr<const DofDistribution<dim,range,rank> >
{
    return dof_distribution_;
}

template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
get_ptr_dof_distribution() -> shared_ptr<DofDistribution<dim,range,rank> >
{
    return dof_distribution_;
}



#ifdef SERIALIZATION

template<int dim_, int range_, int rank_>
template<class Archive>
void
BSplineSpace<dim_, range_, rank_>::
serialize(Archive &ar, const unsigned int version)
{
    ar &boost::serialization::make_nvp("ReferenceSpace",
                                       boost::serialization::base_object<BaseSpace>(*this));

    ar &boost::serialization::make_nvp("space_data_",space_data_);
//    Assert(space_data_ != nullptr,ExcNullPtr());

    ar &boost::serialization::make_nvp("end_b_",end_b_);

    ar &boost::serialization::make_nvp("operators_",operators_);

    ar &boost::serialization::make_nvp("end_interval_",end_interval_);

    ar &boost::serialization::make_nvp("dof_distribution_",dof_distribution_);

//    ar &boost::serialization::make_nvp("dofs_tensor_id_elem_table_",dofs_tensor_id_elem_table_);
}

#endif // SERIALIZATION

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/bspline_space.inst>
