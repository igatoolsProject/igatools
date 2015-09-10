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
#include <igatools/functions/identity_function.h>


using std::shared_ptr;
using std::make_shared;
using std::const_pointer_cast;

IGA_NAMESPACE_OPEN


template<int dim_, int range_, int rank_>
ReferenceSpace<dim_, range_, rank_>::
ReferenceSpace(
  const std::shared_ptr<Grid<dim_>> &grid)
  :
  base_t(grid,std::make_shared<IdentityFunction<dim,dim>>(grid))
{}

template<int dim_, int range_, int rank_>
ReferenceSpace<dim_, range_, rank_>::
ReferenceSpace(
  const std::shared_ptr<const Grid<dim_>> &grid)
  :
  base_t(grid,std::make_shared<IdentityFunction<dim,dim>>(grid))
{}


template<int dim, int range, int rank>
template<int k>
auto
ReferenceSpace<dim, range, rank>::
get_ref_sub_space(const int sub_elem_id,
                  InterSpaceMap<k> &dof_map,
                  std::shared_ptr<Grid<k>> sub_grid) const
-> std::shared_ptr< SubRefSpace<k> >
{
  std::shared_ptr< SubRefSpace<k> > sub_ref_space;
  if (this->is_bspline())
  {
    const auto bsp_space = dynamic_cast<const BSplineSpace<dim,range,rank> *>(this);
    Assert(bsp_space != nullptr,ExcNullPtr());
    sub_ref_space = bsp_space->get_ref_sub_space(sub_elem_id,dof_map,sub_grid);
  }
  else
  {
#ifdef NURBS
    //TODO (MM, Dec 22, 2014): implement NURBSSpace::get_ref_sub_space()
#ifndef NDEBUG
    const auto nrb_space = dynamic_cast<const NURBSSpace<dim,range,rank> *>(this);
#endif
    Assert(nrb_space != nullptr,ExcNullPtr());
    Assert(false,ExcNotImplemented());
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
get_sub_space(const int s_id, InterSpaceMap<k> &dof_map,
              std::shared_ptr<Grid<k>> sub_grid,
              SubGridMap<k> &elem_map) const
-> std::shared_ptr<SubSpace<k> >
{
  std::shared_ptr<SubSpace<k> > sub_space;
  if (this->is_bspline())
  {
    const auto bsp_space =
    dynamic_cast<const BSplineSpace<dim,range,rank> *>(this);
    Assert(bsp_space != nullptr, ExcNullPtr());
    sub_space = bsp_space->get_sub_space(s_id,dof_map,sub_grid, elem_map);
  }
  else
  {
#ifdef NURBS
    const auto nrb_space =
    dynamic_cast<const NURBSSpace<dim,range,rank> *>(this);
    Assert(nrb_space != nullptr, ExcNullPtr());
    sub_space = nrb_space->get_sub_space(s_id,dof_map,sub_grid, elem_map);
#else
    Assert(false,ExcMessage("NURBS support disabled from configuration cmake parameters."));
    AssertThrow(false,ExcMessage("NURBS support disabled from configuration cmake parameters."));
#endif
  }

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

