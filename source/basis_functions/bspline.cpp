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

#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/bspline_handler.h>
#include <igatools/functions/sub_function.h>
//#include <igatools/functions/identity_function.h>
//#include <igatools/functions/grid_function_lib.h>


using std::endl;

using std::shared_ptr;
using std::make_shared;
using std::const_pointer_cast;

IGA_NAMESPACE_OPEN


template<int dim_, int range_, int rank_>
BSpline<dim_, range_, rank_>::
BSpline(const SharedPtrConstnessHandler<SpSpace> &spline_space,
        const EndBehaviourTable &end_b)
  :
  spline_space_(spline_space),
  end_b_(end_b),
  operators_(*spline_space_,end_b),
  end_interval_(end_b.get_comp_map())
{
  //------------------------------------------------------------------------------
// TODO (pauletti, Dec 24, 2014): after it work it should be recoded properly

  const auto &sp_space = *this->spline_space_;
  const auto &grid = *sp_space.get_grid();
  const auto &degree_table = sp_space.get_degree_table();
  const auto rep_knots = sp_space.compute_knots_with_repetition(end_b_);

  const auto &knots_coord = grid.get_knots();
  for (auto i : end_interval_.get_active_components_id())
  {
    const auto &rep_knots_i = rep_knots[i];

    auto &end_interval_i = end_interval_[i];

    for (int dir=0; dir<dim; ++dir)
    {
      const auto p = degree_table[i][dir];

      const auto &knots_coord_dir = *knots_coord[dir];

      const auto &rep_knots_i_dir = rep_knots_i[dir];

      auto &end_interval_i_dir = end_interval_i[dir];

      const auto x1 = knots_coord_dir[1];
      const auto a = knots_coord_dir[0];
      const auto x0 = rep_knots_i_dir[p];
      end_interval_i_dir[0] = (x1-a) / (x1-x0);

      const auto xk= *(knots_coord_dir.end()-2);
      const auto b = *(knots_coord_dir.end()-1);
      const auto xk1 = *(rep_knots_i_dir.end() - (p+1));
      end_interval_i_dir[1] = (b-xk) / (xk1-xk);
    } // end loop dir
  } // end loop i
  //------------------------------------------------------------------------------
}



template<int dim_, int range_, int rank_>
auto
BSpline<dim_, range_, rank_>::
create(const std::shared_ptr<SpSpace> &spline_space,
       const EndBehaviourTable &end_b)
-> shared_ptr<self_t>
{
  auto sp = shared_ptr<self_t>(
    new self_t(SharedPtrConstnessHandler<SpSpace>(spline_space), end_b));
  Assert(sp != nullptr, ExcNullPtr());

#ifdef MESH_REFINEMENT
  sp->create_connection_for_insert_knots(sp);
#endif

  return sp;
}

template<int dim_, int range_, int rank_>
auto
BSpline<dim_, range_, rank_>::
const_create(const std::shared_ptr<const SpSpace> &spline_space,
             const EndBehaviourTable &end_b)
-> shared_ptr<const self_t>
{
  auto sp = shared_ptr<self_t>(
    new self_t(SharedPtrConstnessHandler<SpSpace>(spline_space), end_b));
  Assert(sp != nullptr, ExcNullPtr());

  return sp;
}


template<int dim_, int range_, int rank_>
auto
BSpline<dim_, range_, rank_>::
create(const std::shared_ptr<SpSpace> &spline_space,
       const BasisEndBehaviour &end_b)
-> std::shared_ptr<self_t>
{
  return self_t::create(spline_space,EndBehaviour(end_b));
}

template<int dim_, int range_, int rank_>
auto
BSpline<dim_, range_, rank_>::
const_create(const std::shared_ptr<const SpSpace> &spline_space,
             const BasisEndBehaviour &end_b)
-> std::shared_ptr<const self_t>
{
  return self_t::const_create(spline_space,EndBehaviour(end_b));
}



template<int dim_, int range_, int rank_>
auto
BSpline<dim_, range_, rank_>::
get_this_basis() const -> shared_ptr<const self_t>
{
  auto ref_sp = const_cast<self_t *>(this)->shared_from_this();
  auto bsp_basis = std::dynamic_pointer_cast<self_t>(ref_sp);
  Assert(bsp_basis != nullptr,ExcNullPtr());

  return bsp_basis;
}

template<int dim_, int range_, int rank_>
auto
BSpline<dim_, range_, rank_>::
create_element_begin(const PropId &property) const
-> std::unique_ptr<BasisElement<dim_,0,range_,rank_> >
{
  using Elem = BSplineElement<dim_,range_,rank_>;

//  const auto &id_elems_with_property = this->get_grid()->get_elements_with_property(property);
  return std::make_unique<Elem>(
    this->get_this_basis(),
    this->get_grid()->create_element_begin(property));
}

template<int dim_, int range_, int rank_>
auto
BSpline<dim_, range_, rank_>::
create_element_end(const PropId &property) const
-> std::unique_ptr<BasisElement<dim_,0,range_,rank_> >
{
  using Elem = BSplineElement<dim_,range_,rank_>;

//  const auto &id_elems_with_property = this->get_grid()->get_elements_with_property(property);
  return std::make_unique<Elem>(this->get_this_basis(),
  this->get_grid()->create_element_end(property));
}



template<int dim_, int range_, int rank_>
auto
BSpline<dim_, range_, rank_>::
create_ref_element_begin(const PropId &property) const
-> std::unique_ptr<ReferenceBasisElement<dim_,range_,rank_> >
{
  using Elem = BSplineElement<dim_,range_,rank_>;

  return std::make_unique<Elem>(this->get_this_basis(),
  this->get_grid()->create_element_begin(property));
}

template<int dim_, int range_, int rank_>
auto
BSpline<dim_, range_, rank_>::
create_ref_element_end(const PropId &property) const
-> std::unique_ptr<ReferenceBasisElement<dim_,range_,rank_> >
{
  using Elem = BSplineElement<dim_,range_,rank_>;

//  const auto &id_elems_with_property = this->get_grid()->get_elements_with_property(property);
  return std::make_unique<Elem>(this->get_this_basis(),
  this->get_grid()->create_element_end(property));
}


template<int dim_, int range_, int rank_>
auto
BSpline<dim_, range_, rank_>::
create_bspline_element_begin(const PropId &property) const
-> std::unique_ptr<BSplineElement<dim_,range_,rank_> >
{
  using Elem = BSplineElement<dim_,range_,rank_>;

  return std::make_unique<Elem>(this->get_this_basis(),
  this->get_grid()->create_element_begin(property));
}

template<int dim_, int range_, int rank_>
auto
BSpline<dim_, range_, rank_>::
create_bspline_element_end(const PropId &property) const
-> std::unique_ptr<BSplineElement<dim_,range_,rank_> >
{
  using Elem = BSplineElement<dim_,range_,rank_>;

  return std::make_unique<Elem>(this->get_this_basis(),
  this->get_grid()->create_element_end(property));
}

template<int dim_, int range_, int rank_>
template<int sdim>
auto
BSpline<dim_, range_, rank_>::
get_sub_bspline_basis(const int s_id,
                      InterBasisMap<sdim> &dof_map,
                      const std::shared_ptr<const Grid<sdim>> &sub_grid_in) const
-> std::shared_ptr<const BSpline<sdim, range_, rank_> >
{
  static_assert(sdim == 0 || (sdim > 0 && sdim < dim_),
  "The dimensionality of the sub_grid is not valid.");

  std::shared_ptr<const Grid<sdim>> sub_grid;
  if (sub_grid_in != nullptr)
  {
#ifndef NDEBUG
    typename Grid<dim_>::template SubGridMap<sdim> elem_map;
    sub_grid = this->get_grid()->template get_sub_grid<sdim>(s_id, elem_map);
    Assert(*sub_grid_in == *sub_grid,ExcMessage("Invalid input grid."));
#endif
    sub_grid = sub_grid_in;
  }
  else
  {
    typename Grid<dim_>::template SubGridMap<sdim> elem_map;
    sub_grid = this->get_grid()->template get_sub_grid<sdim>(s_id, elem_map);
  }

  auto sub_mult   = this->spline_space_->template get_sub_space_mult<sdim>(s_id);
  auto sub_degree = this->spline_space_->template get_sub_space_degree<sdim>(s_id);
  auto sub_periodic = this->spline_space_->template get_sub_space_periodicity<sdim>(s_id);

  using SubBasis = BSpline<sdim,range,rank>;

  using SubEndBT = typename SubBasis::EndBehaviourTable;
  auto &k_elem = UnitElement<dim>::template get_elem<sdim>(s_id);
  const auto &active_dirs = k_elem.active_directions;

  SubEndBT sub_end_b(end_b_.get_comp_map());
  for (int comp : end_b_.get_active_components_id())
    for (int j=0; j<sdim; ++j)
      sub_end_b[comp][j] = end_b_[comp][active_dirs[j]];

  auto sub_spline_space =
  SplineSpace<sdim,range_,rank_>::const_create(sub_degree,sub_grid,sub_mult,sub_periodic);
  auto sub_basis = SubBasis::const_create(sub_spline_space, sub_end_b);

  // Creating the mapping between the space degrees of freedom
  const int n_dir = k_elem.constant_directions.size();
  TensorIndex<dim> tensor_index;
  int comp_i = 0;
  dof_map.resize(sub_basis->get_num_basis());
  const auto &sub_space_index_table = sub_basis->get_spline_space()->get_dof_distribution()->get_index_table();
  const auto &space_index_table = this->get_spline_space()->get_dof_distribution()->get_index_table();
  const auto &sub_dof_distribution = *sub_basis->get_spline_space()->get_dof_distribution();
  const auto &dof_distribution = *this->get_spline_space()->get_dof_distribution();
  for (auto comp : SpSpace::components)
  {
    const auto n_basis = sub_dof_distribution.get_num_dofs_comp(comp);
    const auto &sub_local_indices = sub_space_index_table[comp];
    const auto &elem_global_indices = space_index_table[comp];

    for (Index sub_i = 0; sub_i < n_basis; ++sub_i, ++comp_i)
    {
      const auto sub_base_id = sub_local_indices.flat_to_tensor(sub_i);

      for (int j=0; j<sdim; ++j)
        tensor_index[active_dirs[j]] = sub_base_id[j];
      for (int j=0; j<n_dir; ++j)
      {
        auto dir = k_elem.constant_directions[j];
        auto val = k_elem.constant_values[j];
        const int fixed_id = val * (dof_distribution.get_num_dofs_comp(comp,dir) - 1);
        tensor_index[dir] = fixed_id;

      }
      dof_map[comp_i] = elem_global_indices(tensor_index);
    }
  }

  return sub_basis;
}






template<int dim_, int range_, int rank_>
auto
BSpline<dim_, range_, rank_>::
get_grid() const -> std::shared_ptr<const Grid<dim_>>
{
  return spline_space_->get_grid();
}








#ifdef MESH_REFINEMENT



template <int dim_,int range_,int rank_>
void
BSpline<dim_,range_,rank_>::
refine_h(const Size n_subdivisions)
{
  spline_space_.get_ptr_data()->refine_h(n_subdivisions);
}


template<int dim_, int range_, int rank_>
void
BSpline<dim_, range_, rank_>::
rebuild_after_insert_knots(
  const SafeSTLArray<SafeSTLVector<Real>,dim> &knots_to_insert,
  const Grid<dim> &old_grid)
{
  this->ref_basis_previous_refinement_ =
    BSpline<dim_,range_,rank_>::const_create(
      this->spline_space_->get_spline_space_previous_refinement(),
      this->end_b_);

  operators_ = BernsteinExtraction<dim,range,rank>(*this->spline_space_,end_b_);
}



#endif //MESH_REFINEMENT

template<int dim_, int range_, int rank_>
void
BSpline<dim_, range_, rank_>::
print_info(LogStream &out) const
{
  out.begin_item("Spline Space:");
  this->spline_space_->print_info(out);
  out.end_item();


  out.begin_item("DoFs Distribution:");
  this->get_spline_space()->get_dof_distribution()->print_info(out);
  out.end_item();


  out.begin_item("Bernstein Extraction:");
  operators_.print_info(out);
  out.end_item();
}

#if 0
template<int dim_, int range_, int rank_>
auto
BSpline<dim_, range_, rank_>::
get_periodicity() const -> const PeriodicityTable &
{
  return spline_space_->get_periodicity();
}
#endif

template<int dim_, int range_, int rank_>
auto
BSpline<dim_, range_, rank_>::
get_end_behaviour_table() const -> const EndBehaviourTable &
{
  return end_b_;
}


template<int dim_, int range_, int rank_>
bool
BSpline<dim_, range_, rank_>::
is_bspline() const
{
  return true;
}

template<int dim_, int range_, int rank_>
auto
BSpline<dim_, range_, rank_>::
create_cache_handler() const
-> std::unique_ptr<BasisHandler<dim_,0,range_,rank_>>
{
  using Handler = BSplineHandler<dim_,range_,rank_>;
  return std::unique_ptr<Handler>(new Handler(this->get_this_basis()));
}



template<int dim_, int range_, int rank_>
std::shared_ptr<const SplineSpace<dim_,range_,rank_> >
BSpline<dim_, range_, rank_>::
get_spline_space() const
{
  return spline_space_.get_ptr_const_data();
}


#ifdef SERIALIZATION

template<int dim_, int range_, int rank_>
template<class Archive>
void
BSpline<dim_, range_, rank_>::
serialize(Archive &ar)
{
  using std::to_string;
  const std::string base_name = "ReferenceBasis_" +
                                to_string(dim_) + "_" +
                                to_string(0) + "_" +
                                to_string(range_) + "_" +
                                to_string(rank_);

  ar &make_nvp(base_name,base_class<RefBasis>(this));

  ar &make_nvp("spline_space_",spline_space_);

  ar &make_nvp("end_b_",end_b_);

  ar &make_nvp("operators_",operators_);

  ar &make_nvp("end_interval_",end_interval_);


//    ar &make_nvp("dofs_tensor_id_elem_table_",dofs_tensor_id_elem_table_);
}
#endif // SERIALIZATION

IGA_NAMESPACE_CLOSE


#include <igatools/basis_functions/bspline.inst>


#ifdef SERIALIZATION

//using BSpAlias0_1_1 = iga::BSpline<0,1,1>;
//CEREAL_REGISTER_DYNAMIC_INIT(BSpAlias0_1_1);

#endif // SERIALIZATION
