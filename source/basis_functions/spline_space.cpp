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

#include <igatools/basis_functions/spline_space.h>
#include <igatools/base/array_utils.h>
#include <igatools/utils/vector_tools.h>
#include <igatools/utils/multi_array_utils.h>
#include <igatools/basis_functions/dof_distribution.h>
#include <igatools/utils/unique_id_generator.h>


using std::unique_ptr;
using std::shared_ptr;
using std::make_shared;

IGA_NAMESPACE_OPEN

template<int dim_, int range_, int rank_>
const Size SplineSpace<dim_, range_, rank_>::n_components;


template<int dim_, int range_, int rank_>
const SafeSTLArray<Size, SplineSpace<dim_, range_, rank_>::n_components>
SplineSpace<dim_, range_, rank_>::components =
  sequence<SplineSpace<dim_, range_, rank_>::n_components>();




template<int dim_,int range_,int rank_>
SplineSpace<dim_,range_,rank_>::
SplineSpace(const int degree,
            const SharedPtrConstnessHandler<GridType> &grid,
            const InteriorReg interior_reg,
            const bool periodic)
  :
  SplineSpace(Degrees(degree), grid, interior_reg,
             Periodicity(periodic))
{}

template<int dim_,int range_,int rank_>
auto
SplineSpace<dim_,range_,rank_>::
create(const int degree,
       const std::shared_ptr<GridType> &grid,
       const InteriorReg interior_reg,
       const bool periodic) -> shared_ptr<self_t>
{
  auto sp = shared_ptr<self_t>(
    new self_t(degree,SharedPtrConstnessHandler<GridType>(grid),interior_reg,periodic));
  Assert(sp != nullptr, ExcNullPtr());

#ifdef MESH_REFINEMENT
  sp->create_connection_for_insert_knots(sp);
#endif

  return sp;
}

template<int dim_,int range_,int rank_>
auto
SplineSpace<dim_,range_,rank_>::
const_create(const int degree,
             const std::shared_ptr<const GridType> &grid,
             const InteriorReg interior_reg,
             const bool periodic) -> shared_ptr<const self_t>
{
  auto sp = shared_ptr<const self_t>(
    new self_t(degree,SharedPtrConstnessHandler<GridType>(grid),interior_reg,periodic));
  Assert(sp != nullptr, ExcNullPtr());

  return sp;
}


template<int dim_,int range_,int rank_>
SplineSpace<dim_,range_,rank_>::
SplineSpace(const Degrees &deg,
            const SharedPtrConstnessHandler<GridType> &grid,
            const InteriorReg interior_reg,
            const Periodicity &periodic)
  :
  SplineSpace(DegreeTable(true,deg),
             grid,
             self_t::get_multiplicity_from_regularity(interior_reg,DegreeTable(true,deg),
                                                      grid->get_num_intervals()),
             PeriodicityTable(true,periodic))
{}


template<int dim_,int range_,int rank_>
auto
SplineSpace<dim_,range_,rank_>::
create(const Degrees &deg,
       const std::shared_ptr<GridType> &grid,
       const InteriorReg interior_reg,
       const Periodicity &periodic)
-> shared_ptr<self_t>
{
  auto sp = shared_ptr<self_t>(
    new self_t(deg, SharedPtrConstnessHandler<GridType>(grid), interior_reg, periodic));
  Assert(sp != nullptr, ExcNullPtr());

#ifdef MESH_REFINEMENT
  sp->create_connection_for_insert_knots(sp);
#endif

  return sp;
}

template<int dim_,int range_,int rank_>
auto
SplineSpace<dim_,range_,rank_>::
const_create(const Degrees &deg,
             const std::shared_ptr<const GridType> &grid,
             const InteriorReg interior_reg,
             const Periodicity &periodic)
-> shared_ptr<const self_t>
{
  auto sp = shared_ptr<const self_t>(
    new self_t(deg, SharedPtrConstnessHandler<GridType>(grid), interior_reg, periodic));
  Assert(sp != nullptr, ExcNullPtr());

  return sp;
}

#if 0
template<int dim_,int range_,int rank_>
auto
BSpline<dim_,range_,rank_>::
create(const DegreeTable &deg,
       const std::shared_ptr<GridType> &grid,
       const MultiplicityTable &interior_mult,
       const PeriodicityTable &periodic,
       const EndBehaviourTable &end_b)
-> shared_ptr<self_t>
{
  auto sp = shared_ptr<self_t>(
    new self_t(deg, SharedPtrConstnessHandler<GridType>(grid), interior_mult, periodic, end_b));
  Assert(sp != nullptr, ExcNullPtr());

#ifdef MESH_REFINEMENT
  sp->create_connection_for_insert_knots(sp);
#endif

  return sp;
}

template<int dim_,int range_,int rank_>
auto
BSpline<dim_,range_,rank_>::
const_create(const DegreeTable &deg,
             const std::shared_ptr<const GridType> &grid,
             const MultiplicityTable &interior_mult,
             const PeriodicityTable &periodic,
             const EndBehaviourTable &end_b)
-> shared_ptr<const self_t>
{
  auto sp = shared_ptr<const self_t>(
    new self_t(deg, SharedPtrConstnessHandler<GridType>(grid), interior_mult, periodic, end_b));
  Assert(sp != nullptr, ExcNullPtr());

  return sp;
}
#endif



template<int dim_, int range_, int rank_>
SplineSpace<dim_, range_, rank_>::
SplineSpace(const DegreeTable &deg,
            const SharedPtrConstnessHandler<GridType> &grid,
            const MultiplicityTable &interior_mult,
            const PeriodicityTable &periodic)
  :
  grid_(grid),
  interior_mult_(interior_mult),
  deg_(deg),
  periodic_(periodic),
  object_id_(UniqueIdGenerator::get_unique_id())
{
  this->init();
}

template<int dim_, int range_, int rank_>
SplineSpace<dim_, range_, rank_>::
SplineSpace()
  :
  object_id_(UniqueIdGenerator::get_unique_id())
{}



template<int dim_, int range_, int rank_>
std::shared_ptr<SplineSpace<dim_,range_,rank_> >
SplineSpace<dim_, range_, rank_>::
create(const DegreeTable &deg,
       const std::shared_ptr<GridType> &grid,
       const MultiplicityTable &interior_mult,
       const PeriodicityTable &periodic)
{
  using SpSpace = SplineSpace<dim_,range_,rank_>;
  auto sp = std::shared_ptr<SpSpace>(new SpSpace(
                                       deg,
                                       SharedPtrConstnessHandler<Grid<dim_>>(grid),
                                       interior_mult,periodic));
  Assert(sp != nullptr, ExcNullPtr());

#ifdef MESH_REFINEMENT
  sp->create_connection_for_insert_knots(sp);
#endif // MESH_REFINEMENT

  return sp;
}

template<int dim_, int range_, int rank_>
std::shared_ptr<const SplineSpace<dim_,range_,rank_> >
SplineSpace<dim_, range_, rank_>::
const_create(const DegreeTable &deg,
             const std::shared_ptr<const GridType> &grid,
             const MultiplicityTable &interior_mult,
             const PeriodicityTable &periodic)
{
  using SpSpace = SplineSpace<dim_,range_,rank_>;
  auto sp = std::shared_ptr<const SpSpace>(new SpSpace(
                                             deg,
                                             SharedPtrConstnessHandler<Grid<dim_>>(grid),
                                             interior_mult,periodic));
  Assert(sp != nullptr, ExcNullPtr());

  return sp;
}


template<int dim_, int range_, int rank_>
void
SplineSpace<dim_, range_, rank_>::
init()
{
//    Assert(grid_ != nullptr,ExcNullPtr());

  //------------------------------------------------------------------------------
  // the default value of a bool variable is undefined, so we need to set
  // set the values of the inactive components of the perodicity table to true or false (we use false)
  using PeriodicAsMArray = StaticMultiArray<Periodicity,range_,rank_>;
  PeriodicAsMArray &periodic_as_static_m_array = static_cast<PeriodicAsMArray &>(periodic_);
  for (const auto inactive_id : periodic_.get_inactive_components_id())
    periodic_as_static_m_array[inactive_id] = Periodicity(false);
  //------------------------------------------------------------------------------



  //------------------------------------------------------------------------------
  // Determine the dimensionality of the spline space --- begin
  typename TensorSizeTable::base_t n_basis;
  for (const auto comp : components)
  {
    const auto &deg_comp =              deg_[comp];
    const auto &mult_comp = interior_mult_[comp];

    const auto &periodic_comp = periodic_[comp];

    for (const auto dir : UnitElement<dim_>::active_directions)
    {
      const auto &deg =  deg_comp[dir];
      const auto &mult = mult_comp.get_data_direction(dir);

      const auto order = deg + 1;

#ifndef NDEBUG
      Assert(mult.size() == grid_->get_num_knots_dim()[dir]-2,
             ExcMessage("Interior multiplicity size does not match the grid"));
      if (!mult.empty())
      {
        auto result = std::minmax_element(mult.begin(), mult.end());
        Assert((*result.first > 0) && (*result.second <= order),
               ExcMessage("multiplicity values not between 0 and p+1"));
      }
#endif

      n_basis[comp][dir] = std::accumulate(
                             mult.begin(),
                             mult.end(),
                             periodic_comp[dir] ? 0 : order);

#ifndef NDEBUG
      if (periodic_comp[dir])
        Assert(n_basis[comp][dir] > order,
               ExcMessage("Not enough basis functions"));
#endif

    } // end loop dir
  } // end loop comp

  space_dim_ = n_basis;
  // Determine the dimensionality of the spline space --- end
  //------------------------------------------------------------------------------







  //------------------------------------------------------------------------------
  // building the lookup table for the local dof id on the current component of an element --- begin
  for (const auto comp : components)
  {
    const auto dofs_t_size_elem_comp = TensorSize<dim_>(deg_[comp]+1);
    const auto dofs_f_size_elem_comp = dofs_t_size_elem_comp.flat_size();

    auto &elem_comp_dof_t_id = dofs_tensor_id_elem_table_[comp];
    elem_comp_dof_t_id.resize(dofs_f_size_elem_comp);

    const auto w_dofs_elem_comp = MultiArrayUtils<dim_>::compute_weight(dofs_t_size_elem_comp);

    for (int dof_f_id = 0 ; dof_f_id < dofs_f_size_elem_comp ; ++dof_f_id)
      elem_comp_dof_t_id[dof_f_id] = MultiArrayUtils<dim_>::flat_to_tensor_index(dof_f_id,w_dofs_elem_comp);
  }
  // building the lookup table for the local dof id on the current component of an element --- end
  //------------------------------------------------------------------------------



  //------------------------------------------------------------------------------
  dof_distribution_ = std::make_shared<DofDistribution<dim_,range_,rank_>>(
                        this->get_num_basis_table(),
                        this->get_degree_table(),
                        this->get_periodic_table());
  Assert(dof_distribution_ != nullptr, ExcNullPtr());
  //------------------------------------------------------------------------------

}



template<int dim_, int range_, int rank_>
std::shared_ptr<const Grid<dim_> >
SplineSpace<dim_, range_, rank_>::
get_grid() const
{
  return grid_.get_ptr_const_data();
}

template<int dim_, int range_, int rank_>
int
SplineSpace<dim_, range_, rank_>::
get_max_degree() const
{
  int max_degree = 0;

  const auto &degree_table = this->get_degree_table();
  for (const auto &degree_comp : degree_table)
    for (const auto &degree_comp_dim : degree_comp)
      max_degree = std::max(max_degree,degree_comp_dim);

  return max_degree;
}


#ifdef MESH_REFINEMENT

template<int dim_, int range_, int rank_>
void
SplineSpace<dim_, range_, rank_>::
refine_h(const Size n_subdivisions)
{
  grid_.get_ptr_data()->refine(n_subdivisions);
//  Assert(false,ExcNotImplemented());
}


template<int dim_, int range_, int rank_>
void
SplineSpace<dim_, range_, rank_>::
rebuild_after_insert_knots(
  const SafeSTLArray<SafeSTLVector<Real>,dim_> &knots_to_insert,
  const Grid<dim_> &old_grid)
{
//    const auto refined_grid = grid_;
  auto grid_pre_refinement = grid_.get_ptr_data()->get_grid_pre_refinement();

  Assert(&(*grid_pre_refinement) == &old_grid,ExcMessage("Different grids."));

  spline_space_previous_refinement_ =
    SplineSpace<dim_,range_,rank_>::const_create(
      deg_,
      grid_pre_refinement,
      interior_mult_,
      periodic_);


  const auto &old_unique_knots = old_grid.get_knots();

#ifndef NDEBUG
  //---------------------------------------------------------------------------------------
  // check that the new knots are internal to the grid --- begin
  for (const auto dir : UnitElement<dim_>::active_directions)
  {
    const auto &knots_dir = *old_unique_knots[dir];

    for (const auto &knt_value : knots_to_insert[dir])
    {
      Assert(knt_value > knots_dir.front() && knt_value < knots_dir.back(),
             ExcMessage("The knot value" + std::to_string(knt_value) +
                        " inserted along the direction " + std::to_string(dir) +
                        " is not contained in the open interval (" +
                        std::to_string(knots_dir.front()) + " , " +
                        std::to_string(knots_dir.back()) + ")"));
    }
  }
  // check that the new knots are internal to the grid --- end
  //---------------------------------------------------------------------------------------
#endif


  //---------------------------------------------------------------------------------------
  // build the new internal knots with repetitions --- begin
  for (const auto comp : components)
  {
    const auto &interior_mult_comp = spline_space_previous_refinement_->interior_mult_[comp];

    for (const auto dir : UnitElement<dim_>::active_directions)
    {
      const auto &old_unique_knots_dir = *old_unique_knots[dir];
      const auto &interior_mult_comp_dir = interior_mult_comp.get_data_direction(dir);

      const int n_internal_knots_old = old_unique_knots_dir.size() - 2;
      Assert(n_internal_knots_old == interior_mult_comp_dir.size(),
             ExcDimensionMismatch(n_internal_knots_old,interior_mult_comp_dir.size()));

      SafeSTLVector<Real> new_internal_knots(knots_to_insert[dir]);
      for (int i = 0 ; i < n_internal_knots_old ; ++i)
      {
        new_internal_knots.insert(new_internal_knots.end(),
                                  interior_mult_comp_dir[i],
                                  old_unique_knots_dir[i+1]);
      }
      std::sort(new_internal_knots.begin(),new_internal_knots.end());

      SafeSTLVector<Real> new_internal_knots_unique;
      SafeSTLVector<int> new_internal_mult_dir;
      vector_tools::count_and_remove_duplicates(
        new_internal_knots,
        new_internal_knots_unique,
        new_internal_mult_dir);

#ifndef NDEBUG
      //---------------------------------------------------------------------------------------
      // check that the new unique knots are the same of the new grid --- begin
      const auto &unique_knots_dir_refined_grid = grid_->get_knot_coordinates(dir);

      const int n_internal_knots_unique_new = unique_knots_dir_refined_grid.size() - 2;
      Assert(n_internal_knots_unique_new == new_internal_knots_unique.size(),
             ExcDimensionMismatch(n_internal_knots_unique_new,new_internal_knots_unique.size()));
      for (int i = 0 ; i < n_internal_knots_unique_new ; ++i)
      {
        Assert(new_internal_knots_unique[i] == unique_knots_dir_refined_grid[i+1],
               ExcMessage(
                 "The knot value " + std::to_string(new_internal_knots_unique[i]) +
                 " is not present along the direction " + std::to_string(dir) +
                 " of the refined grid."));
      }
      // check that the new unique knots are the same of the new grid --- end
      //---------------------------------------------------------------------------------------


      //---------------------------------------------------------------------------------------
      // check that the multiplicity is compatible with the degree --- begin
      for (const auto mult : new_internal_mult_dir)
        Assert(mult <= deg_[comp][dir],ExcIndexRange(mult,1,deg_[comp][dir]+1))

        // check that the multiplicity is compatible with the degree --- end
        //---------------------------------------------------------------------------------------
#endif

        interior_mult_[comp].copy_data_direction(dir,new_internal_mult_dir);
    }
  }
  // build the new internal knots with repetitions --- end
  //---------------------------------------------------------------------------------------

  this->init();
}


template<int dim_, int range_, int rank_>
void
SplineSpace<dim_, range_, rank_>::
create_connection_for_insert_knots(std::shared_ptr<SplineSpace<dim_,range_,rank_>> space)
{
  Assert(space != nullptr, ExcNullPtr());
  Assert(&(*space) == &(*this), ExcMessage("Different objects."));

  auto func_to_connect =
    std::bind(&SplineSpace<dim_,range_,rank_>::rebuild_after_insert_knots,
              space.get(),
              std::placeholders::_1,
              std::placeholders::_2);

  using SlotType = typename Grid<dim_>::SignalInsertKnotsSlot;
  grid_.get_ptr_data()->connect_insert_knots(
    SlotType(func_to_connect).track_foreign(space));
}

#endif // MESH_REFINEMENT



template<int dim_, int range_, int rank_>
auto
SplineSpace<dim_, range_, rank_>::
compute_knots_with_repetition(const EndBehaviourTable &ends,
                              const BoundaryKnotsTable &boundary_knots) const
-> KnotsTable
{

#ifndef NDEBUG
  for (auto iComp : components)
  {
    for (const int j : UnitElement<dim_>::active_directions)
    {
      const auto &l_knots = boundary_knots[iComp][j].get_data_direction(0);
      const auto &r_knots = boundary_knots[iComp][j].get_data_direction(1);

      if (periodic_[iComp][j])
      {
        Assert(ends[iComp][j] == BasisEndBehaviour::periodic,
               ExcMessage("Periodic inconsistency"));
        Assert(l_knots.size()==0,
               ExcMessage("Periodic inconsistency"));
        Assert(r_knots.size()==0,
               ExcMessage("Periodic inconsistency"));
      }
      else
      {
        if (ends[iComp][j] == BasisEndBehaviour::interpolatory)
        {
          Assert(l_knots.size()==0,
                 ExcMessage("Interpolatory inconsistency"));
          Assert(r_knots.size()==0,
                 ExcMessage("Interpolatory inconsistency"));
        }
        if (ends[iComp][j] == BasisEndBehaviour::end_knots)
        {
          const auto &knots = grid_->get_knot_coordinates(j);
          const int m = deg_[iComp][j] + 1;
          Assert(l_knots.size() == m,
                 ExcMessage("Wrong number of boundary knots"));
          Assert(r_knots.size() == m,
                 ExcMessage("Wrong number of boundary knots"));
          Assert(knots.front() >= l_knots.back(),
                 ExcMessage("Boundary knots should be smaller or equal a"));
          Assert(knots.back() <= r_knots.front(),
                 ExcMessage("Boundary knots should be greater or equal b"));
          Assert(std::is_sorted(l_knots.begin(), l_knots.end()),
                 ExcMessage("Boundary knots is not sorted"));
          Assert(std::is_sorted(r_knots.begin(), r_knots.end()),
                 ExcMessage("Boundary knots is not sorted"));
        }
      }
    }
  }
#endif

  KnotsTable result;

  for (int comp = 0; comp < n_components; ++comp)
  {
    const auto   &degree_comp = deg_[comp];
    const auto     &mult_comp = interior_mult_[comp];
    const auto &periodic_comp = periodic_[comp];

    for (const auto dir : UnitElement<dim_>::active_directions)
    {
      const auto deg = degree_comp[dir];
      const auto order = deg + 1;
      const auto &knots = grid_->get_knot_coordinates(dir);
      const auto &mult  = mult_comp.get_data_direction(dir);

      const int m = order;
      const int K = std::accumulate(mult.begin(),mult.end(),0);


      SafeSTLVector<Real> rep_knots(2*order+K);

      auto rep_it = rep_knots.begin() + m;
      auto m_it = mult.begin();
      auto k_it = ++knots.begin();
      auto end = mult.end();
      for (; m_it !=end; ++m_it, ++k_it)
      {
        for (int iMult = 0; iMult < *m_it; ++iMult, ++rep_it)
          *rep_it = *k_it;
      }


      if (periodic_comp[dir])
      {
        const Real a = knots.front();
        const Real b = knots.back();
        const Real L = b - a;
        for (int i=0; i<m; ++i)
        {
          rep_knots[i] = rep_knots[K+i] - L;
          rep_knots[K+m+i] = rep_knots[m+i] + L;
        }
      }
      else
      {
        const auto &endb = ends[comp][dir];
        switch (endb)
        {
          case BasisEndBehaviour::interpolatory:
          {
            const Real a = knots.front();
            const Real b = knots.back();

            for (int i=0; i<m; ++i)
            {
              rep_knots[i] = a;
              rep_knots[m+K+i] = b;
            }
          }
          break;
          case BasisEndBehaviour::end_knots:
          {
            const auto &left_knts = boundary_knots[comp][dir].get_data_direction(0);
            const auto &right_knts = boundary_knots[comp][dir].get_data_direction(1);

            for (int i=0; i<m; ++i)
            {
              rep_knots[i]     = left_knts[i];
              rep_knots[K+m+i] = right_knts[i];
            }
          }
          break;
          case BasisEndBehaviour::periodic:
            Assert(false, ExcMessage("Impossible"));
        }
      }

      result[comp][dir] = rep_knots;
    } // end loop dir
  } // end loop comp

  return result;
}



template<int dim_, int range_, int rank_>
template<int k>
auto
SplineSpace<dim_, range_, rank_>::
get_sub_space_mult(const Index sub_elem_id) const
-> typename SubSpace<k>::MultiplicityTable
{
  using SubMultT = typename SubSpace<k>::MultiplicityTable;
  const auto &v_mult = interior_mult_;

  auto &k_elem = UnitElement<dim_>::template get_elem<k>(sub_elem_id);
  const auto &active_dirs = k_elem.active_directions;

  auto sub_mult = SubMultT(v_mult.get_comp_map());
  for (int comp : sub_mult.get_active_components_id())
  {
    for (int j=0; j<k; ++j)
      sub_mult[comp].copy_data_direction(j, v_mult[comp].get_data_direction(active_dirs[j]));
  }
  return sub_mult;
}



template<int dim_, int range_, int rank_>
template<int k>
auto
SplineSpace<dim_, range_, rank_>::
get_sub_space_degree(const Index sub_elem_id) const
-> typename SubSpace<k>::DegreeTable
{
  using SubDegreeT = typename SubSpace<k>::DegreeTable;
  auto &k_elem = UnitElement<dim_>::template get_elem<k>(sub_elem_id);
  const auto &active_dirs = k_elem.active_directions;

  SubDegreeT sub_degree(deg_.get_comp_map());

  for (int comp : sub_degree.get_active_components_id())
    for (int j=0; j<k; ++j)
      sub_degree[comp][j] = deg_[comp][active_dirs[j]];

  return sub_degree;
}



template<int dim_, int range_, int rank_>
template<int k>
auto
SplineSpace<dim_, range_, rank_>::
get_sub_space_periodicity(const Index sub_elem_id) const
-> typename SubSpace<k>::PeriodicityTable
{
  using SubPeriodicT = typename SubSpace<k>::PeriodicityTable;
  auto &k_elem = UnitElement<dim_>::template get_elem<k>(sub_elem_id);
  const auto &active_dirs = k_elem.active_directions;

  SubPeriodicT sub_periodic(periodic_.get_comp_map());

  for (int comp : periodic_.get_active_components_id())
    for (int j=0; j<k; ++j)
      sub_periodic[comp][j] = periodic_[comp][active_dirs[j]];

  return sub_periodic;
}



template<int dim_, int range_, int rank_>
auto SplineSpace<dim_, range_, rank_>::
accumulated_interior_multiplicities() const -> MultiplicityTable
{
  MultiplicityTable result;
  for (const auto comp : components)
  {
    const auto &mult_comp = interior_mult_[comp];

    for (const auto j : UnitElement<dim_>::active_directions)
    {
      // Assert(!periodic_[iComp][j], ExcMessage("periodic needs to be implemented"));
      const auto &mult = mult_comp.get_data_direction(j);

      SafeSTLVector<Size> accum_mult;
      const int size = mult.size();
      accum_mult.reserve(size + 1);
      accum_mult.push_back(0);
      for (int i = 0; i < size; ++i)
        accum_mult.push_back(accum_mult[i] + mult[i]);

      result[comp].copy_data_direction(j, accum_mult);

      //TODO(pauletti, May 3, 2014): write some post assertions
    }
  }
  return result;
}



template<int dim_, int range_, int rank_>
auto
SplineSpace<dim_, range_, rank_>::
get_multiplicity_from_regularity(const InteriorReg reg,
                                 const DegreeTable &deg,
                                 const TensorSize<dim_> &n_elem)
-> MultiplicityTable
{
  auto res = MultiplicityTable(deg.get_comp_map());


  for (int comp : res.get_active_components_id())
    for (const auto dir : UnitElement<dim_>::active_directions)
    {
      int val;
      switch (reg)
      {
        case InteriorReg::maximum:
          val = 1;
          break;
        case InteriorReg::minimun:
          val = deg[comp][dir] + 1;
          break;
      }

      const auto size = n_elem[dir]-1;

      SafeSTLVector<Size> mult(size);
      if (size>0)
        res[comp].copy_data_direction(dir, SafeSTLVector<Size>(size, val));
    }
  return res;
}



template<int dim_, int range_, int rank_>
void
SplineSpace<dim_, range_, rank_>::
get_element_dofs(
  const typename GridType::IndexType &elem_id,
  SafeSTLVector<Index> &dofs_global,
  SafeSTLVector<Index> &dofs_local_to_patch,
  SafeSTLVector<Index> &dofs_local_to_elem,
  const std::string &dofs_property) const
{
  const auto &accum_mult = this->accumulated_interior_multiplicities();

  const auto &dof_distr = *(this->get_dof_distribution());
  const auto &index_table = dof_distr.get_index_table();

  dofs_global.clear();
  dofs_local_to_patch.clear();
  dofs_local_to_elem.clear();

  const auto &elem_t_id = elem_id.get_tensor_index();

  Index dof_loc_to_elem = 0;
  for (const auto comp : components)
  {
    const auto &index_table_comp = index_table[comp];

    const auto dof_t_origin = accum_mult[comp].cartesian_product(elem_t_id);

    const auto &elem_comp_dof_t_id = this->get_dofs_tensor_id_elem_table()[comp];

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

template<int dim_, int range_, int rank_>
void
SplineSpace<dim_, range_, rank_>::
print_info(LogStream &out) const
{
  out.begin_item("Knots without repetition:");
  grid_->print_info(out);
  out.end_item();

  out.begin_item("Degrees:");
  deg_.print_info(out);
  out.end_item();

  out.begin_item("Interior multiplicities:");
  const MultiplicityTable &interior_mult_ref = interior_mult_;
  for (const auto &v : interior_mult_ref)
    v.print_info(out);
  out.end_item();

  out.begin_item("Dimensionality Table:");
  space_dim_.print_info(out);
  out.end_item();
}



template<int dim_, int range_, int rank_>
Index
SplineSpace<dim_, range_, rank_>::
get_object_id() const
{
  return object_id_;
}


template<int dim_, int range_, int rank_>
auto
SplineSpace<dim_, range_, rank_>::
get_degree_table() const -> const DegreeTable &
{
  return deg_;
}

template<int dim_, int range_, int rank_>
auto
SplineSpace<dim_, range_, rank_>::
get_periodicity() const -> const PeriodicityTable &
{
  return periodic_;
}

template<int dim_, int range_, int rank_>
auto
SplineSpace<dim_, range_, rank_>::
get_components_map() const -> const SafeSTLArray<Index,n_components> &
{
  return interior_mult_.get_comp_map();
}

template<int dim_, int range_, int rank_>
auto
SplineSpace<dim_, range_, rank_>::
get_active_components_id() const -> const SafeSTLVector<Index> &
{
  return interior_mult_.get_active_components_id();
}

template<int dim_, int range_, int rank_>
Size
SplineSpace<dim_, range_, rank_>::
get_num_basis() const
{
  return space_dim_.total_dimension();
}

template<int dim_, int range_, int rank_>
Size
SplineSpace<dim_, range_, rank_>::
get_num_basis(const int comp) const
{
  return space_dim_.get_component_size(comp);
}

template<int dim_, int range_, int rank_>
Size
SplineSpace<dim_, range_, rank_>::
get_num_basis(const int comp, const int dir) const
{
  return space_dim_[comp][dir];
}

template<int dim_, int range_, int rank_>
auto
SplineSpace<dim_, range_, rank_>::
get_num_basis_table() const -> const TensorSizeTable &
{
  return space_dim_;
}


template<int dim_, int range_, int rank_>
auto
SplineSpace<dim_, range_, rank_>::
get_interior_mult() const -> const MultiplicityTable &
{
  return interior_mult_;
}

template<int dim_, int range_, int rank_>
auto
SplineSpace<dim_, range_, rank_>::
get_periodic_table() const -> const PeriodicityTable &
{
  return periodic_;
}


template<int dim_, int range_, int rank_>
auto
SplineSpace<dim_, range_, rank_>::
get_dofs_tensor_id_elem_table() const
-> const ComponentContainer<SafeSTLVector<TensorIndex<dim_> > > &
{
  return dofs_tensor_id_elem_table_;
}


template<int dim_, int range_, int rank_>
auto
SplineSpace<dim_, range_, rank_>::
get_dof_distribution() const ->
std::shared_ptr<const DofDistribution<dim_,range_,rank_> >
{
  return dof_distribution_;
}

template<int dim_, int range_, int rank_>
auto
SplineSpace<dim_, range_, rank_>::
get_dof_distribution() ->
std::shared_ptr<DofDistribution<dim_,range_,rank_> >
{
  return dof_distribution_;
}

#ifdef SERIALIZATION
template<int dim_, int range_, int rank_>
template<class Archive>
void
SplineSpace<dim_, range_, rank_>::
serialize(Archive &ar)
{
  ar &make_nvp("grid_",grid_);

  ar &make_nvp("interior_mult_",interior_mult_);

  ar &make_nvp("deg_", deg_);

  ar &make_nvp("space_dim_", space_dim_);

  ar &make_nvp("periodic_", periodic_);

  ar &make_nvp("dofs_tensor_id_elem_table_",dofs_tensor_id_elem_table_);

  ar &make_nvp("dof_distribution_",dof_distribution_);

#ifdef MESH_REFINEMENT
  using self_t = SplineSpace<dim_,range_,rank_>;
  auto tmp = std::const_pointer_cast<self_t>(spline_space_previous_refinement_);
  ar &make_nvp("spline_space_previous_refinement_",tmp);
  spline_space_previous_refinement_ = tmp;
#endif


}
#endif // SERIALIZATION


template<int dim_, int range_, int rank_>
Size
SplineSpace<dim_, range_, rank_>::
TensorSizeTable::
get_component_size(const int comp) const
{
  return (*this)[comp].flat_size();
}


template<int dim_, int range_, int rank_>
Size
SplineSpace<dim_, range_, rank_>::
TensorSizeTable::
total_dimension() const
{
  Index total_dimension = 0;
  for (const auto comp : components)
    total_dimension += this->get_component_size(comp);

  return total_dimension;
}

template<int dim_, int range_, int rank_>
auto
SplineSpace<dim_, range_, rank_>::
TensorSizeTable::
get_offset() const -> ComponentContainer<Size>
{
  ComponentContainer<Size> offset;
  offset[0] = 0;
  for (int comp = 1; comp < n_components; ++comp)
    offset[comp] = offset[comp-1] + this->get_component_size(comp-1);

  return offset;

}

template<int dim_, int range_, int rank_>
void
SplineSpace<dim_, range_, rank_>::
TensorSizeTable::
print_info(LogStream &out) const
{
  base_t::print_info(out);

  out.begin_item("Scalar components dimensions:");
  out << "[ ";
  for (auto comp : components)
    out << this->get_component_size(comp) << " ";
  out << "]";
  out.end_item();

  out << "Total Dimension: " << total_dimension() << std::endl;
}


template<int dim_, int range_, int rank_>
SplineSpace<dim_, range_, rank_>::
TensorSizeTable::
TensorSizeTable(const base_t &n_basis)
  :
  base_t(n_basis)
{}





IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/spline_space.inst>
