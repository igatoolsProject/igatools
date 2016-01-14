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

#include <igatools/basis_functions/nurbs.h>
#include <igatools/basis_functions/nurbs_element_handler.h>
//#include <igatools/basis_functions/space_tools.h>

//#include <igatools/base/sub_function.h>
//#include <igatools/base/exceptions.h>

#ifdef USE_NURBS


using std::to_string;
using std::endl;
using std::shared_ptr;
using std::make_shared;
using std::bind;
using std::placeholders::_1;
using std::placeholders::_2;

IGA_NAMESPACE_OPEN



template <int dim_, int range_, int rank_>
NURBS<dim_, range_, rank_>::
NURBS(const SharedPtrConstnessHandler<BSpBasis> &bsp_basis,
      const SharedPtrConstnessHandler<WeightFunction> &weight_func)
  :
  bsp_basis_(bsp_basis),
  weight_func_(weight_func)
{
#ifndef NDEBUG
  Assert(this->get_grid() == weight_func_->get_grid(),ExcMessage("Mismatching grids."));

  const auto w_as_ig_func =
    std::dynamic_pointer_cast<const IgGridFunction<dim,1>>(weight_func_.get_ptr_const_data());
  Assert(w_as_ig_func != nullptr,ExcNullPtr());

  const auto w_func_basis = w_as_ig_func->get_basis();
  Assert(w_func_basis->is_bspline(),
         ExcMessage("The space for the weight function is not BSpline."));

  const auto &n_basis_table = this->get_spline_space()->get_dof_distribution()->get_num_dofs_table();
  int comp_id = 0;
  for (const auto &n_basis_comp : n_basis_table)
  {
    Assert(n_basis_comp == w_func_basis->get_spline_space()->get_dof_distribution()->get_num_dofs_table()[0],
           ExcMessage("Mismatching number of basis functions and weight "
                      "coefficients for scalar component " + to_string(comp_id)));

    ++comp_id;
  }
#endif
}



template <int dim_, int range_, int rank_>
auto
NURBS<dim_, range_, rank_>::
create(const std::shared_ptr<BSpBasis> &bs_space,
       const std::shared_ptr<WeightFunction> &weight_func) -> shared_ptr<self_t>
{
  auto sp = shared_ptr<self_t>(
    new self_t(SharedPtrConstnessHandler<BSpBasis>(bs_space),
  SharedPtrConstnessHandler<WeightFunction>(weight_func)));
  Assert(sp != nullptr, ExcNullPtr());

#ifdef MESH_REFINEMENT
  sp->create_connection_for_insert_knots(sp);
#endif // MESH_REFINEMENT

  return sp;
}

template <int dim_, int range_, int rank_>
auto
NURBS<dim_, range_, rank_>::
const_create(const std::shared_ptr<const BSpBasis> &bs_space,
             const std::shared_ptr<const WeightFunction> &weight_func) -> shared_ptr<const self_t>
{
  auto sp = shared_ptr<self_t>(
    new self_t(SharedPtrConstnessHandler<BSpBasis>(bs_space),
  SharedPtrConstnessHandler<WeightFunction>(weight_func)));
  Assert(sp != nullptr, ExcNullPtr());

  return sp;
}


template<int dim_, int range_, int rank_>
auto
NURBS<dim_, range_, rank_>::
get_this_basis() const -> std::shared_ptr<const self_t >
{
  auto ref_sp = const_cast<self_t *>(this)->shared_from_this();
  auto nrb_space = std::dynamic_pointer_cast<self_t>(ref_sp);
  Assert(nrb_space != nullptr,ExcNullPtr());

  return nrb_space;
}




template<int dim_, int range_, int rank_>
auto
NURBS<dim_, range_, rank_>::
create_element_begin(const PropId &property) const
-> std::unique_ptr<BasisElement<dim_,0,range_,rank_> >
{
  using Elem = NURBSElement<dim_,range_,rank_>;

//  const auto &id_elems_with_property = this->get_grid()->get_elements_with_property(property);
  return std::make_unique<Elem>(this->get_this_basis(),
  bsp_basis_->create_bspline_element_begin(property),
  weight_func_->create_element_begin(property));
}

template<int dim_, int range_, int rank_>
auto
NURBS<dim_, range_, rank_>::
create_element_end(const PropId &property) const
-> std::unique_ptr<BasisElement<dim_,0,range_,rank_> >
{
  using Elem = NURBSElement<dim_,range_,rank_>;

//  const auto &id_elems_with_property = this->get_grid()->get_elements_with_property(property);
  return std::make_unique<Elem>(this->get_this_basis(),
  bsp_basis_->create_bspline_element_end(property),
  weight_func_->create_element_end(property));
}


template<int dim_, int range_, int rank_>
auto
NURBS<dim_, range_, rank_>::
create_ref_element_begin(const PropId &property) const
-> std::unique_ptr<ReferenceElement<dim_,range_,rank_> >
{
  using Elem = NURBSElement<dim_,range_,rank_>;

  return std::make_unique<Elem>(this->get_this_basis(),
  bsp_basis_->create_bspline_element_begin(property),
  weight_func_->create_element_begin(property));
}

template<int dim_, int range_, int rank_>
auto
NURBS<dim_, range_, rank_>::
create_ref_element_end(const PropId &property) const
-> std::unique_ptr<ReferenceElement<dim_,range_,rank_> >
{
  using Elem = NURBSElement<dim_,range_,rank_>;

  return std::make_unique<Elem>(this->get_this_basis(),
  bsp_basis_->create_bspline_element_end(property),
  weight_func_->create_element_end(property));
}

template <int dim_, int range_, int rank_>
auto
NURBS<dim_, range_, rank_>::
get_weight_func() const -> std::shared_ptr<const WeightFunction>
{
  return weight_func_.get_ptr_const_data();
}

template <int dim_, int range_, int rank_>
auto
NURBS<dim_, range_, rank_>::
get_grid() const -> std::shared_ptr<const Grid<dim_>>
{
  return bsp_basis_->get_grid();
}



#if 0
template <int dim_, int range_, int rank_>
void
NURBS<dim_, range_, rank_>::
reset_weights(const WeightsTable &weights)
{
  weights_ = weights;
  perform_post_construction_checks();
}


template<int dim_, int range_, int rank_>
auto
NURBS<dim_, range_, rank_>::
get_ref_face_space(const Index face_id,
                   SafeSTLVector<Index> &face_to_element_dofs,
                   typename GridType::FaceGridMap &elem_map) const
-> std::shared_ptr<RefFaceSpace>
{
  auto f_space = bsp_basis_->get_ref_face_space(face_id, face_to_element_dofs, elem_map);

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
NURBS<dim_, range_, rank_>::
get_face_space(const Index face_id,
               SafeSTLVector<Index> &face_to_element_dofs) const
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
NURBS<dim_, range_, rank_>::
refine_h_weights(
  const SafeSTLArray<bool,dim> &refinement_directions,
  const GridType &grid_old1)
{
  auto grid = this->get_grid();
  auto grid_old = this->get_grid()->get_grid_pre_refinement();

  auto knots_with_repetitions_pre_refinement = bsp_basis_->get_spline_space_previous_refinement()
                                               ->compute_knots_with_repetition(
                                                 bsp_basis_->get_end_behaviour());
  auto knots_with_repetitions = bsp_basis_->compute_knots_with_repetition(
                                  bsp_basis_->get_end_behaviour());

  for (int direction_id = 0; direction_id < dim; ++direction_id)
  {
    if (refinement_directions[direction_id])
    {
      for (const int comp_id : weights_.get_active_components_id())
      {
        const int p = bsp_basis_->get_degree_table()[comp_id][direction_id];
        const auto &U = knots_with_repetitions_pre_refinement[comp_id].get_data_direction(direction_id);
        const auto &Ubar = knots_with_repetitions[comp_id].get_data_direction(direction_id);


        SafeSTLVector<Real> knots_added(Ubar.size());

        // find the knots in the refined space that are not present in the old space
        auto it = std::set_difference(
                    Ubar.begin(),Ubar.end(),
                    U.begin(),U.end(),
                    knots_added.begin());

        knots_added.resize(it-knots_added.begin());
        const auto &X = knots_added;


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
               bsp_basis_->get_num_basis(comp_id,direction_id),
               ExcDimensionMismatch(new_sizes[direction_id],
                                    bsp_basis_->get_num_basis(comp_id,direction_id)));

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
auto
NURBS<dim_, range_, rank_>::
get_bspline_basis() const -> const std::shared_ptr<const BSpBasis>
{
  return bsp_basis_.get_ptr_const_data();
}





template<int dim_, int range_, int rank_>
template<int sdim>
auto
NURBS<dim_, range_, rank_>::
get_sub_nurbs_space(const int s_id,
                    InterSpaceMap<sdim> &dof_map,
                    const std::shared_ptr<const Grid<sdim>> &sub_grid) const
-> std::shared_ptr<const NURBS<sdim,range_,rank_> >
{
  static_assert(sdim == 0 || (sdim > 0 && sdim < dim_),
  "The dimensionality of the sub_grid is not valid.");

  auto sub_bsp_basis = bsp_basis_->template get_sub_bspline_space<sdim>(s_id,dof_map,sub_grid);
  auto space_sub_grid = sub_bsp_basis->get_grid();

  auto sub_w_func = weight_func_->template get_sub_function<sdim>(s_id,space_sub_grid);

  auto sub_nrb_space = NURBS<sdim,range_,rank_>::const_create(sub_bsp_basis,sub_w_func);

  return sub_nrb_space;
}


#if 0
template<int dim_, int range_, int rank_>
template<int k>
auto
NURBS<dim_, range_, rank_>::
get_sub_space(const int s_id, InterSpaceMap<k> &dof_map,
              SubGridMap<k> &elem_map) const
-> std::shared_ptr<SubSpace<k> >
{
  //TODO (martinelli Nov 27,2014): implement this function
  static_assert(k == 0 || (k > 0 && k < dim_),
  "The dimensionality of the sub_grid is not valid.");

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
#endif

template <int dim_, int range_, int rank_>
void
NURBS<dim_, range_, rank_>::
print_info(LogStream &out) const
{
  out.begin_item("BSpline Space:");
  bsp_basis_->print_info(out);
  out.end_item();

  out.begin_item("Weight function :");
  weight_func_->print_info(out);
  out.end_item();
}



template <int dim_, int range_, int rank_>
bool
NURBS<dim_, range_, rank_>::
is_bspline() const
{
  return false;
}

#if 0
template <int dim_, int range_, int rank_>
auto
NURBS<dim_, range_, rank_>::
get_degree_table() const -> const DegreeTable &
{
  return this->bsp_basis_->get_degree_table();
}
#endif





#if 0
template <int dim_, int range_, int rank_>
auto
NURBS<dim_, range_, rank_>::
get_periodicity() const -> const PeriodicityTable &
{
  return bsp_basis_->get_periodicity();
}
#endif

template <int dim_, int range_, int rank_>
auto
NURBS<dim_, range_, rank_>::
get_end_behaviour_table() const -> const EndBehaviourTable &
{
  return bsp_basis_->get_end_behaviour_table();
}


template <int dim_, int range_, int rank_>
auto
NURBS<dim_, range_, rank_>::
create_cache_handler() const -> std::unique_ptr<BasisElementHandler<dim_,0,range_,rank_>>
{
  return std::unique_ptr<ElementHandler>(new ElementHandler(this->get_this_basis()));
}


template<int dim_, int range_, int rank_>
std::shared_ptr<const SplineSpace<dim_,range_,rank_>>
                                                   NURBS<dim_, range_, rank_>::
                                                   get_spline_space() const
{
  return bsp_basis_->get_spline_space();
}


#ifdef MESH_REFINEMENT


template<int dim_, int range_, int rank_>
void
NURBS<dim_, range_, rank_>::
refine_h(const Size n_subdivisions)
{
  //the refinement of the BSpline also refines the weight_fucntion (they share the same Grid)
  bsp_basis_.get_ptr_data()->refine_h(n_subdivisions);
}

template<int dim_, int range_, int rank_>
void
NURBS<dim_, range_, rank_>::
rebuild_after_insert_knots(
  const SafeSTLArray<SafeSTLVector<Real>,dim_> &knots_to_insert,
  const Grid<dim_> &old_grid)
{
  auto bsp_basis_pre_refinement =
    std::dynamic_pointer_cast<const BSpBasis>(bsp_basis_->get_basis_previous_refinement());
  Assert(bsp_basis_pre_refinement != nullptr,ExcNullPtr());

  auto weight_func_pre_refinement_ =
    std::dynamic_pointer_cast<const WeightFunction>(weight_func_->get_grid_function_previous_refinement());
  Assert(weight_func_pre_refinement_ != nullptr,ExcNullPtr());


  this->ref_basis_previous_refinement_ =
    self_t::const_create(bsp_basis_pre_refinement,weight_func_pre_refinement_);
}




#endif //MESH_REFINEMENT


#ifdef SERIALIZATION

template<int dim_, int range_, int rank_>
template<class Archive>
void
NURBS<dim_, range_, rank_>::
serialize(Archive &ar)
{
  using std::to_string;
  const std::string base_name = "ReferenceBasis_" +
                                to_string(dim_) + "_" +
                                to_string(0) + "_" +
                                to_string(range_) + "_" +
                                to_string(rank_);

  ar &make_nvp(base_name,base_class<BaseSpace>(this));
  ar &make_nvp("bsp_basis_",bsp_basis_);

  ar &make_nvp("weight_func_",weight_func_);
}
///@}

#endif // SERIALIZATION

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/nurbs.inst>

#endif // #ifdef USE_NURBS

