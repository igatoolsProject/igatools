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


#include <igatools/basis_functions/reference_space.h>
#include <igatools/basis_functions/reference_element.h>
#include <igatools/base/array_utils.h>

#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/nurbs_space.h>
#include <igatools/geometry/grid_function_lib.h>


using std::shared_ptr;
using std::make_shared;
using std::const_pointer_cast;

IGA_NAMESPACE_OPEN



template<int dim, int range, int rank>
template<int k>
auto
ReferenceSpace<dim, range, rank>::
get_ref_sub_space(const int sub_elem_id,
                  InterSpaceMap<k> &dof_map) const
-> std::shared_ptr< SubRefSpace<k> >
{
  static_assert(k == 0 || (k > 0 && k < dim),
  "The dimensionality of the sub_grid is not valid.");

  std::shared_ptr< SubRefSpace<k> > sub_ref_space;
  if (this->is_bspline())
  {
    const auto bsp_space = dynamic_cast<const BSplineSpace<dim,range,rank> *>(this);
    Assert(bsp_space != nullptr,ExcNullPtr());
    sub_ref_space = bsp_space->template get_sub_bspline_space<k>(sub_elem_id,dof_map);
  }
  else
  {
#ifdef NURBS
    const auto nrb_space = dynamic_cast<const NURBSSpace<dim,range,rank> *>(this);
    Assert(nrb_space != nullptr,ExcNullPtr());
    sub_ref_space = nrb_space->template get_sub_nurbs_space<k>(sub_elem_id,dof_map);
#else
    Assert(false,ExcMessage("NURBS support disabled from configuration cmake parameters."));
    AssertThrow(false,ExcMessage("NURBS support disabled from configuration cmake parameters."));
#endif
  }

  Assert(sub_ref_space != nullptr, ExcNullPtr());
  return sub_ref_space;
}



template<int dim, int range, int rank>
template<int k>
auto
ReferenceSpace<dim, range, rank>::
get_sub_space(const int s_id,
              InterSpaceMap<k> &dof_map,
              SubGridMap<k> &elem_map) const
-> std::shared_ptr<SubSpace<k> >
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


  const auto grid = this->get_ptr_const_grid();

  i = 0;
  for (const int constant_dir : constant_dirs)
  {
    const int constant_val = constant_vals[i++];

    const auto &knots_const_direction = grid->get_knot_coordinates(constant_dir);

    b[constant_dir] = (constant_val == 0) ?
    knots_const_direction.front() :
    knots_const_direction.back();
  }

  auto sub_ref_space = this->template get_ref_sub_space<k>(s_id, dof_map);

  auto sub_grid = sub_ref_space->get_ptr_grid();

  auto sub_grid_func = SubGridFunc::create(sub_grid,A,b);

  using SubDomain = Domain<k,dim-k>;
  auto sub_domain = SubDomain::create(sub_grid_func);


  auto sub_space = SubSpace<k>::create(sub_ref_space, sub_domain);

  Assert(sub_space != nullptr, ExcNullPtr());
  return sub_space;
}







template<int dim, int range, int rank>
int
ReferenceSpace<dim, range, rank>::
get_max_degree() const
{
  int max_degree = 0;

  const auto &degree_table = this->get_degree_table();
  for (const auto &degree_comp : degree_table)
    for (const auto &degree_comp_dim : degree_comp)
      max_degree = std::max(max_degree,degree_comp_dim);

  return max_degree;
}



#if 0
#ifdef SERIALIZATION
template<int dim, int range, int rank>
template<class Archive>
void
ReferenceSpace<dim, range, rank>::
serialize(Archive &ar, const unsigned int version)
{
  ar &boost::serialization::make_nvp("ReferenceSpace_base_t",
                                     boost::serialization::base_object<base_t>(*this));

  auto tmp = const_pointer_cast<RefSpace>(ref_space_previous_refinement_);
  ar &boost::serialization::make_nvp("ref_space_previous_refinement_",tmp);
  ref_space_previous_refinement_ = const_pointer_cast<const RefSpace>(tmp);
}
#endif // SERIALIZATION
#endif

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/reference_space.inst>

