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


#include <igatools/basis_functions/reference_space_basis.h>
#include <igatools/basis_functions/reference_element.h>
#include <igatools/base/array_utils.h>

#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/nurbs.h>
#include <igatools/functions/grid_function_lib.h>


using std::shared_ptr;
using std::make_shared;
using std::const_pointer_cast;

IGA_NAMESPACE_OPEN



template<int dim_,int range_,int rank_>
template<int sdim>
auto
ReferenceSpaceBasis<dim_,range_,rank_>::
get_ref_sub_space(const int sub_elem_id,
                  InterSpaceMap<sdim> &dof_map,
                  const std::shared_ptr<const Grid<sdim>> &sub_grid) const
-> std::shared_ptr<const SubRefSpace<sdim> >
{
  static_assert(sdim == 0 || (sdim > 0 && sdim < dim),
  "The dimensionality of the sub_grid is not valid.");

  std::shared_ptr< const SubRefSpace<sdim> > sub_ref_space;
  if (this->is_bspline())
  {
    const auto bsp_space = dynamic_cast<const BSpline<dim,range,rank> *>(this);
    Assert(bsp_space != nullptr,ExcNullPtr());
    sub_ref_space = bsp_space->template get_sub_bspline_space<sdim>(sub_elem_id,dof_map,sub_grid);
  }
  else
  {
#ifdef USE_NURBS
    const auto nrb_space = dynamic_cast<const NURBS<dim,range,rank> *>(this);
    Assert(nrb_space != nullptr,ExcNullPtr());
    sub_ref_space = nrb_space->template get_sub_nurbs_space<sdim>(sub_elem_id,dof_map,sub_grid);
#else
    Assert(false,ExcMessage("NURBS support disabled from configuration cmake parameters."));
    AssertThrow(false,ExcMessage("NURBS support disabled from configuration cmake parameters."));
#endif
  }

  Assert(sub_ref_space != nullptr, ExcNullPtr());
  return sub_ref_space;
}



template<int dim_,int range_,int rank_>
template<int k>
auto
ReferenceSpaceBasis<dim_,range_,rank_>::
get_sub_space(const int s_id,
              InterSpaceMap<k> &dof_map,
              SubGridMap<k> &elem_map) const
-> std::shared_ptr<const SubSpace<k> >
{
  static_assert(k == 0 || (k > 0 && k < dim),
  "The dimensionality of the sub_grid is not valid.");

  using SubGridFunc = grid_functions::LinearGridFunction<k,dim>;
  using Grad  = typename SubGridFunc::template Derivative<1>;
  using Value = typename SubGridFunc::Value;
  Grad A;
  Value b;

  const auto &sub_elem = UnitElement<dim>::template get_elem<k>(s_id);
  const auto &active_dirs = sub_elem.active_directions;
  const auto &constant_dirs = sub_elem.constant_directions;
  const auto &constant_vals = sub_elem.constant_values;

  int i = 0;
  for (const int active_dir : active_dirs)
    A[i++][active_dir] = 1.0;


  const auto grid = this->get_grid();

  i = 0;
  for (const int constant_dir : constant_dirs)
  {
    const int constant_val = constant_vals[i++];

    const auto &knots_const_direction = grid->get_knot_coordinates(constant_dir);

    b[constant_dir] = (constant_val == 0) ?
    knots_const_direction.front() :
    knots_const_direction.back();
  }

  auto sub_ref_space = this->template get_ref_sub_space<k>(s_id,dof_map,nullptr);

  auto sub_grid = sub_ref_space->get_grid();

  auto sub_grid_func = SubGridFunc::const_create(sub_grid,A,b);

  using SubDomain = Domain<k,dim-k>;
  auto sub_domain = SubDomain::const_create(sub_grid_func);


  auto sub_space = SubSpace<k>::const_create(sub_ref_space, sub_domain);

  Assert(sub_space != nullptr, ExcNullPtr());
  return sub_space;
}







#if 0
template<int dim_,int range_,int rank_>
int
ReferenceSpaceBasis<dim_,range_,rank_>::
get_max_degree() const
{
  return this->get_spline_space()->get_max_degree();
}
#endif


#ifdef MESH_REFINEMENT
template<int dim_,int range_,int rank_>
void
ReferenceSpaceBasis<dim_,range_,rank_>::
create_connection_for_insert_knots(const std::shared_ptr<self_t> &space)
{
  Assert(space != nullptr, ExcNullPtr());
  Assert(&(*space) == &(*this), ExcMessage("Different objects."));

  auto func_to_connect =
    std::bind(&self_t::rebuild_after_insert_knots,
              space.get(),
              std::placeholders::_1,
              std::placeholders::_2);

  using SlotType = typename Grid<dim>::SignalInsertKnotsSlot;
  std::const_pointer_cast<Grid<dim>>(this->get_grid())->connect_insert_knots(SlotType(func_to_connect).track_foreign(space));
}


template<int dim_,int range_,int rank_>
auto
ReferenceSpaceBasis<dim_,range_,rank_>::
get_basis_previous_refinement() const -> std::shared_ptr<const self_t>
{
  return ref_basis_previous_refinement_;
}

#endif // MESH_REFINEMENT



#ifdef SERIALIZATION

template<int dim_,int range_,int rank_>
template<class Archive>
void
ReferenceSpaceBasis<dim_,range_,rank_>::
serialize(Archive &ar)
{
  using std::to_string;
  const std::string base_name = "Space_" +
                                to_string(dim_) + "_" +
                                to_string(0) + "_" +
                                to_string(range_) + "_" +
                                to_string(rank_) + "_hgrad";

  ar &make_nvp(base_name,base_class<base_t>(this));

#ifdef MESH_REFINEMENT
  auto tmp = std::const_pointer_cast<RefBasis>(ref_basis_previous_refinement_);
  ar &make_nvp("ref_basis_previous_refinement_",tmp);
  ref_basis_previous_refinement_ = std::const_pointer_cast<const RefBasis>(tmp);
#endif
}
#endif // SERIALIZATION


IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/reference_space_basis.inst>

