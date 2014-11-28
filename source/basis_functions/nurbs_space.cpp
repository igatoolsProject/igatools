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

#include <igatools/basis_functions/nurbs_space.h>
#include <igatools/basis_functions/space_manager.h>
//#include <igatools/basis_functions/space_tools.h>

//#include <igatools/base/sub_function.h>
//#include <igatools/base/exceptions.h>

#ifdef NURBS
using std::array;

using std::endl;
using std::shared_ptr;
using std::make_shared;
using std::bind;
using std::placeholders::_1;
using std::placeholders::_2;

IGA_NAMESPACE_OPEN


#if 0

template <int dim_, int range_, int rank_>
NURBSSpace<dim_, range_, rank_>::
NURBSSpace(const int degree,
           shared_ptr< GridType > knots,
           const WeightsTable &weights)
    :
    NURBSSpace(DegreeTable(TensorIndex<dim>(degree)), knots, weights)
{}



template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
create(const int degree,
       shared_ptr< GridType > knots,
       const WeightsTable &weights) -> shared_ptr<self_t>
{
    return shared_ptr<self_t>(new self_t(degree, knots,  weights));
}



template <int dim_, int range_, int rank_>
NURBSSpace<dim_, range_, rank_>::
NURBSSpace(const DegreeTable &degree,
           shared_ptr<GridType> knots,
           const WeightsTable &weights)
    :
    BaseSpace(knots),
    sp_space_(SpSpace::create(degree,knots)),
    weights_(weights)
{

    create_refinement_connection();
    perform_post_construction_checks();
}



template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
create(const DegreeTable &degree, shared_ptr<GridType> knots,
       const WeightsTable &weights)
-> shared_ptr<self_t>
{
    return shared_ptr<self_t>(new self_t(degree, knots, weights));
}



template <int dim_, int range_, int rank_>
NURBSSpace<dim_, range_, rank_>::
NURBSSpace(const DegreeTable &deg,
           std::shared_ptr<GridType> knots,
           std::shared_ptr<const MultiplicityTable> interior_mult,
           const EndBehaviourTable &ends,
           const WeightsTable &weights)
    :
    BaseSpace(knots),
    sp_space_(SpSpace::create(deg, knots, interior_mult, ends)),
    weights_(weights)
{
    create_refinement_connection();
    perform_post_construction_checks();
}



template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
create(const DegreeTable &deg,
       std::shared_ptr<GridType> knots,
       std::shared_ptr<const MultiplicityTable> interior_mult,
       const EndBehaviourTable &ends,
       const WeightsTable &weights) -> shared_ptr<self_t>
{
    return shared_ptr<self_t>(new self_t(deg, knots, interior_mult, ends, weights));
}

#endif

template <int dim_, int range_, int rank_>
NURBSSpace<dim_, range_, rank_>::
NURBSSpace(std::shared_ptr<SpSpace> bs_space,
           const WeightsTable &weights,
           const WeightFunctionPtr weight_func)
    :
    BaseSpace(bs_space->get_grid()),
    sp_space_(bs_space),
    weights_(weights),
    weight_func_(weight_func)
{
    Assert(weight_func_ != nullptr, ExcNullPtr());

    Assert(*this->get_grid() == *weight_func_->get_grid(),ExcMessage("Mismatching grids."));

#if 0
    const auto &degree_table = sp_space_->get_degree_table();
    const auto w_grid = CartesianGrid<dim>::create(*bs_space->get_grid());

    for (const auto &comp : weights_.get_active_components_id())
    {
        Assert(false,ExcNotImplemented());
        const auto &interior_mult_tmp = (*sp_space_->get_interior_mult())[comp];


        auto w_sp = WeightSpace::create(degree_table[comp], w_grid, std::shared_ptr< const MultiplicityTable > interior_mult, const EndBehaviourTable &ends=EndBehaviourTable())

                    const vector<Real> coeffs = weights_[comp].get_data();
        weights_func[comp] = WeightFunction::create();
    }
#endif

//    create_refinement_connection();
//    perform_post_construction_checks();
}


template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
create(std::shared_ptr<SpSpace> bs_space,
       const WeightsTable &weights,
       const WeightFunctionPtr weight_func) -> shared_ptr<self_t>
{
    return shared_ptr<self_t>(new self_t(bs_space, weights,weight_func));
}


#if 0
template <int dim_, int range_, int rank_>
void
NURBSSpace<dim_, range_, rank_>::
create_refinement_connection()
{
    // create a signal and a connection for the grid refinement
    this->connect_refinement_h_function(
        bind(&self_t::refine_h_weights, this, std::placeholders::_1, std::placeholders::_2));
}
template <int dim_, int range_, int rank_>
void
NURBSSpace<dim_, range_, rank_>::
perform_post_construction_checks() const
{
#ifndef NDEBUG
    // check that the number of weights is equal to the number of basis functions in the space
    for (auto comp : components)
    {
        Assert(sp_space_->get_num_basis(comp) == weights_[comp].flat_size(),
               ExcDimensionMismatch(sp_space_->get_num_basis(comp),weights_[comp].flat_size()));
    }
#endif
}
#endif


template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
begin() const -> ElementIterator
{
    return ElementIterator(std::enable_shared_from_this<self_t>::shared_from_this(), 0);
}



template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
last() const -> ElementIterator
{
    return ElementIterator(std::enable_shared_from_this<self_t>::shared_from_this(),
                           this->get_grid()->get_num_active_elems() - 1);
}



template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
end() const -> ElementIterator
{
    return ElementIterator(std::enable_shared_from_this<self_t>::shared_from_this(),
                           IteratorState::pass_the_end);
}


#if 0
template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
get_weights() const -> const WeightsTable &
{
    return weights_;
}



template <int dim_, int range_, int rank_>
void
NURBSSpace<dim_, range_, rank_>::
reset_weights(const WeightsTable &weights)
{
    weights_ = weights;
    perform_post_construction_checks();
}


template<int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
get_ref_face_space(const Index face_id,
                   vector<Index> &face_to_element_dofs,
                   typename GridType::FaceGridMap &elem_map) const
-> std::shared_ptr<RefFaceSpace>
{
    auto f_space = sp_space_->get_ref_face_space(face_id, face_to_element_dofs, elem_map);

    // TODO (pauletti, Jun 11, 2014): this should be put and completed in
    // get_face_weigjts()
    const auto &v_weights = weights_;
    //const auto &active_dirs = UnitElement<dim>::face_active_directions[face_id];
    typename RefFaceSpace::WeightsTable f_weights(v_weights.get_comp_map());

    const auto n_basis = f_space->get_num_basis_table();
    for (int comp : f_weights.get_active_components_id())
    {
        f_weights[comp].resize(n_basis[comp],1.0);
        //        for (auto j : RefFaceSpace::dims)
        //            f_weights(comp).copy_data_direction(j, v_weights(comp).get_data_direction(active_dirs[j]));
    }


    return RefFaceSpace::create(f_space, f_weights);
}



template<int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
get_face_space(const Index face_id,
               vector<Index> &face_to_element_dofs) const
-> std::shared_ptr<FaceSpace>
{
    auto elem_map = std::make_shared<typename GridType::FaceGridMap>();
    auto face_ref_sp = get_ref_face_space(face_id, face_to_element_dofs, *elem_map);
    auto map  = get_push_forward()->get_mapping();

    auto fmap = MappingSlice<FaceSpace::PushForwardType::dim, FaceSpace::PushForwardType::codim>::
    create(map, face_id, face_ref_sp->get_grid(), elem_map);
    auto fpf = FaceSpace::PushForwardType::create(fmap);
    auto face_space = FaceSpace::create(face_ref_sp,fpf);

    return face_space;
}

template <int dim_, int range_, int rank_>
void
NURBSSpace<dim_, range_, rank_>::
refine_h_weights(
    const std::array<bool,dim> &refinement_directions,
    const GridType &grid_old1)
{
    auto grid = this->get_grid();
    auto grid_old = this->get_grid()->get_grid_pre_refinement();

    auto knots_with_repetitions_pre_refinement = sp_space_->get_spline_space_previous_refinement()
                                                 ->compute_knots_with_repetition(
                                                     sp_space_->get_end_behaviour());
    auto knots_with_repetitions = sp_space_->compute_knots_with_repetition(
                                      sp_space_->get_end_behaviour());

    for (int direction_id = 0; direction_id < dim; ++direction_id)
    {
        if (refinement_directions[direction_id])
        {
            // knots in the refined grid along the selected direction
            vector<Real> knots_new = grid->get_knot_coordinates(direction_id);

            // knots in the original (unrefined) grid along the selected direction
            vector<Real> knots_old = grid_old->get_knot_coordinates(direction_id);

            vector<Real> knots_added(knots_new.size());

            // find the knots in the refined grid that are not present in the old grid
            auto it = std::set_difference(
                          knots_new.begin(),knots_new.end(),
                          knots_old.begin(),knots_old.end(),
                          knots_added.begin());

            knots_added.resize(it-knots_added.begin());
            /*
                        Assert(false,ExcNotImplemented());
                        AssertThrow(false,ExcNotImplemented());
            //*/
            for (const int comp_id : weights_.get_active_components_id())
            {
                const int p = sp_space_->get_degree()[comp_id][direction_id];
                const auto &U = knots_with_repetitions_pre_refinement[comp_id].get_data_direction(direction_id);
                const auto &X = knots_added;
                const auto &Ubar = knots_with_repetitions[comp_id].get_data_direction(direction_id);


                const int m = U.size()-1;
                const int r = X.size()-1;
                const int a = space_tools::find_span(p,X[0],U);
                const int b = space_tools::find_span(p,X[r],U)+1;

                const int n = m-p-1;

                const auto Pw = weights_[comp_id];
                const auto old_sizes = Pw.tensor_size();
                Assert(old_sizes[direction_id] == n+1,
                       ExcDimensionMismatch(old_sizes[direction_id], n+1));

                auto new_sizes = old_sizes;
                new_sizes[direction_id] += r+1; // r+1 new weights in the refinement direction
                Assert(new_sizes[direction_id] ==
                       sp_space_->get_num_basis(comp_id,direction_id),
                       ExcDimensionMismatch(new_sizes[direction_id],
                                            sp_space_->get_num_basis(comp_id,direction_id)));

                DynamicMultiArray<Real,dim> Qw(new_sizes);

                for (Index j = 0; j <= a-p; ++j)
                    Qw.copy_slice(direction_id,j,Pw.get_slice(direction_id,j));

                for (Index j = b-1; j <= n; ++j)
                    Qw.copy_slice(direction_id,j+r+1,Pw.get_slice(direction_id,j));

                Index i = b + p - 1;
                Index k = b + p + r;
                for (Index j = r; j >= 0; --j)
                {
                    while (X[j] <= U[i] && i > a)
                    {
                        Qw.copy_slice(direction_id,k-p-1,Pw.get_slice(direction_id,i-p-1));
                        k = k-1;
                        i = i-1;
                    }
                    Qw.copy_slice(direction_id,k-p-1,
                                  Qw.get_slice(direction_id,k-p));

                    for (Index l = 1; l <= p; ++l)
                    {
                        Index ind = k-p+l;

                        Real alfa = Ubar[k+l] - X[j];
                        if (fabs(alfa) == 0.0)
                        {
                            Qw.copy_slice(direction_id,ind-1,Qw.get_slice(direction_id,ind));
                        }
                        else
                        {
                            alfa = alfa / (Ubar[k+l] - U[i-p+l]);

                            Qw.copy_slice(direction_id,ind-1,
                                          alfa  * Qw.get_slice(direction_id,ind-1) +
                                          (1.0-alfa) * Qw.get_slice(direction_id,ind));
                        }
                    } // end loop l
                    k = k-1;
                } // end loop j

                weights_[comp_id] = Qw;

            } // end loop comp_id

        } // end if (refinement_directions[direction_id])

    } // end loop direction_id

    this->perform_post_construction_checks();
}

#endif

template <int dim_, int range_, int rank_>
Size
NURBSSpace<dim_, range_, rank_>::
get_num_basis() const
{
    return sp_space_->get_num_basis();
}

template <int dim_, int range_, int rank_>
Size
NURBSSpace<dim_, range_, rank_>::
get_num_basis(const int i) const
{
    return sp_space_->get_num_basis(i);
}

template <int dim_, int range_, int rank_>
Size
NURBSSpace<dim_, range_, rank_>::
get_num_basis(const int comp, const int dir) const
{
    return sp_space_->get_num_basis(comp, dir);
}

template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
get_num_basis_table() const ->
const SpaceDimensionTable &
{
    return sp_space_->get_num_basis_table();
}

template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
get_degree() const -> const DegreeTable &
{
    return sp_space_->get_degree();
}

#if 0
template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
get_interior_mult() const -> std::shared_ptr<const MultiplicityTable>
{
    return sp_space_->get_interior_mult();
}
#endif

template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
get_loc_to_global(const CartesianGridElement<dim> &element) const -> vector<Index>
{
    return sp_space_->get_loc_to_global(element);
}



template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
get_loc_to_patch(const CartesianGridElement<dim> &element) const -> vector<Index>
{
    return sp_space_->get_loc_to_patch(element);
}



template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
get_spline_space() const -> const std::shared_ptr<SpSpace>
{
    return sp_space_;
}



template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
get_dof_distribution_global() const -> const DofDistribution<dim, range, rank> &
{
    return sp_space_->get_dof_distribution_global();
}



template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
get_dof_distribution_global() -> DofDistribution<dim, range, rank> &
{
    return sp_space_->get_dof_distribution_global();
}



template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
get_dof_distribution_patch() const -> const DofDistribution<dim, range, rank> &
{
    return sp_space_->get_dof_distribution_patch();
}



template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
get_dof_distribution_patch() -> DofDistribution<dim, range, rank> &
{
    return sp_space_->get_dof_distribution_patch();
}

#if 0
template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
get_push_forward() -> std::shared_ptr<PushForwardType>
{
    return sp_space_->get_push_forward();
}

template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
get_push_forward() const -> std::shared_ptr<const PushForwardType>
{
    return sp_space_->get_push_forward();
}

template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
get_reference_space() const -> std::shared_ptr<const self_t >
{
    return this->shared_from_this();
}
#endif

template <int dim_, int range_, int rank_>
void
NURBSSpace<dim_, range_, rank_>::
add_dofs_offset(const Index offset)
{
    sp_space_->add_dofs_offset(offset);
}


template<int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
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
NURBSSpace<dim_, range_, rank_>::
get_space_manager() const -> std::shared_ptr<const SpaceManager>
{
    return const_cast<self_t &>(*this).get_space_manager();
}



template<int dim_, int range_, int rank_>
template<int k>
auto
NURBSSpace<dim_, range_, rank_>::
get_ref_sub_space(const int s_id,
                  InterSpaceMap<k> &dof_map,
                  std::shared_ptr<CartesianGrid<k>> sub_grid) const
-> std::shared_ptr<SubRefSpace<k> >
{
    //TODO (martinelli Nov 27,2014): implement this function
#if 0
    if (!(sub_grid))
    {
        typename GridType::template InterGridMap<k>  elem_map;
        sub_grid   = this->get_grid()->template get_sub_grid<k>(s_id, elem_map);
    }
    auto sub_mult   = this->template get_sub_space_mult<k>(s_id);
    auto sub_degree = this->template get_sub_space_degree<k>(s_id);

    auto sub_space = SubRefSpace<k>::create(sub_degree, sub_grid, sub_mult);

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
#endif
    Assert(false,ExcNotImplemented());
    AssertThrow(false,ExcNotImplemented());
    return nullptr;
}



template<int dim_, int range_, int rank_>
template<int k>
auto
NURBSSpace<dim_, range_, rank_>::
get_sub_space(const int s_id, InterSpaceMap<k> &dof_map,
              std::shared_ptr<CartesianGrid<k>> sub_grid,
              std::shared_ptr<typename GridType::template InterGridMap<k>> elem_map) const
-> std::shared_ptr<SubSpace<k> >
{
    //TODO (martinelli Nov 27,2014): implement this function
#if 0
    using SubMap = SubMapFunction<k, dim, space_dim>;
    auto grid =  this->get_grid();
//    typename GridType::template InterGridMap<k> elem_map;
//    auto sub_grid = this->get_grid()->template get_sub_grid<k>(s_id, elem_map);

    auto sub_ref_space = get_ref_sub_space(s_id, dof_map, sub_grid);
    auto F = IdentityFunction<dim>::create(grid);
    auto sub_map_func = SubMap::create(sub_grid, F, s_id, *elem_map);
    auto sub_space = SubSpace<k>::create(sub_ref_space, sub_map_func);
    return sub_space;
#endif
    Assert(false,ExcNotImplemented());
    AssertThrow(false,ExcNotImplemented());
    return nullptr;
}

template <int dim_, int range_, int rank_>
void
NURBSSpace<dim_, range_, rank_>::
print_info(LogStream &out) const
{
    out.begin_item("BSpline Space:");
    sp_space_->print_info(out);
    out.end_item();

    out.begin_item("Weights:");
    for (auto w : weights_)
    {
        w.print_info(out);
    }
    out.end_item();
}


IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/nurbs_space.inst>

#endif // #ifdef NURBS

