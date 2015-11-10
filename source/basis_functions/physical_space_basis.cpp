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

#include <igatools/basis_functions/physical_space_basis.h>
#include <igatools/functions/function.h>
#include <igatools/functions/sub_function.h>
#include <igatools/basis_functions/phys_space_element_handler.h>
#include <igatools/geometry/push_forward.h>


using std::shared_ptr;
using std::make_shared;
using std::unique_ptr;

using std::const_pointer_cast;
using std::endl;

IGA_NAMESPACE_OPEN


template <int dim_, int range_, int rank_, int codim_>
const SafeSTLArray<int, PhysicalSpaceBasis<dim_, range_, rank_, codim_>::n_components>
PhysicalSpaceBasis<dim_, range_, rank_, codim_>::components =
  sequence<PhysicalSpaceBasis<dim_, range_, rank_, codim_>::n_components>();


template <int dim_, int range_, int rank_, int codim_>
PhysicalSpaceBasis<dim_, range_, rank_, codim_>::
PhysicalSpaceBasis(const SharedPtrConstnessHandler<RefBasis> &ref_basis,
                   const SharedPtrConstnessHandler<PhysDomain> &phys_domain,
                   const Transformation &transformation_type)
  :
  ref_basis_(ref_basis),
  phys_domain_(phys_domain),
  transformation_type_(transformation_type)
{
  Assert(this->get_grid() == phys_domain_->get_grid_function()->get_grid(),
         ExcMessage("The space and the physical domain must have the same grid!"));

//TODO(pauletti, Jan 18, 2014): put static assert on h_div, h_curl range and rank
}



template <int dim_, int range_, int rank_, int codim_>
auto
PhysicalSpaceBasis<dim_, range_, rank_, codim_>::
create(const shared_ptr<RefBasis> &ref_basis,
       const shared_ptr<PhysDomain> &phys_domain,
       const Transformation &transformation_type) -> shared_ptr<self_t>
{
  auto sp = shared_ptr<self_t>(
    new self_t(SharedPtrConstnessHandler<RefBasis>(ref_basis),
  SharedPtrConstnessHandler<PhysDomain>(phys_domain),
  transformation_type));
  Assert(sp != nullptr,ExcNullPtr());

#ifdef MESH_REFINEMENT
  sp->create_connection_for_insert_knots(sp);
#endif // MESH_REFINEMENT

  return sp;
}

template <int dim_, int range_, int rank_, int codim_>
auto
PhysicalSpaceBasis<dim_, range_, rank_, codim_>::
const_create(const shared_ptr<const RefBasis> &ref_basis,
             const shared_ptr<const PhysDomain> &phys_domain,
             const Transformation &transformation_type) -> shared_ptr<const self_t>
{
  auto sp = shared_ptr<const self_t>(
    new self_t(SharedPtrConstnessHandler<RefBasis>(ref_basis),
  SharedPtrConstnessHandler<PhysDomain>(phys_domain),transformation_type));
  Assert(sp != nullptr,ExcNullPtr());

  return sp;
}

template <int dim_, int range_, int rank_, int codim_>
auto
PhysicalSpaceBasis<dim_, range_, rank_, codim_>::
get_this_basis() const -> std::shared_ptr<const self_t >
{
  auto sp = const_cast<self_t *>(this)->shared_from_this();
  auto this_space = std::dynamic_pointer_cast<self_t>(sp);
  Assert(this_space != nullptr,ExcNullPtr());

  return this_space;
}


template <int dim_, int range_, int rank_, int codim_>
auto
PhysicalSpaceBasis<dim_, range_, rank_, codim_>::
create_element(const ListIt &index, const PropId &property) const
-> std::unique_ptr<SpaceElement<dim_,codim_,range_,rank_>>
{
  std::unique_ptr<SpaceElement<dim_,codim_,range_,rank_>>
  elem = std::make_unique<ElementAccessor>(this->get_this_basis(),index,property);
  Assert(elem != nullptr, ExcNullPtr());

  return elem;
}



template <int dim_, int range_, int rank_, int codim_>
auto
PhysicalSpaceBasis<dim_, range_, rank_, codim_>::
get_reference_basis() const -> shared_ptr<const RefBasis>
{
  return ref_basis_.get_ptr_const_data();
}

template <int dim_, int range_, int rank_, int codim_>
auto
PhysicalSpaceBasis<dim_, range_, rank_, codim_>::
get_physical_domain() const -> std::shared_ptr<const Domain<dim_,codim_>>
{
  return phys_domain_.get_ptr_const_data();
}

template <int dim_, int range_, int rank_, int codim_>
auto
PhysicalSpaceBasis<dim_, range_, rank_, codim_>::
get_grid() const -> std::shared_ptr<const Grid<dim_>>
{
  return ref_basis_->get_grid();
}


#if 0
template <int dim_, int range_, int rank_, int codim_>
auto
PhysicalSpaceBasis<dim_, range_, rank_, codim_>::
get_reference_basis() -> shared_ptr<RefBasis>
{
  return ref_basis_.get_ptr_data();
}
#endif


template <int dim_, int range_, int rank_, int codim_>
template<int sdim>
auto
PhysicalSpaceBasis<dim_, range_, rank_, codim_>::
get_sub_space(const int s_id, InterSpaceMap<sdim> &dof_map,
              const std::shared_ptr<const Grid<sdim>> &sub_grid,
              SubGridMap<sdim> &elem_map,
              EnableIf<(dim_ != 0) &&(sdim>=0)> *) const
-> std::shared_ptr<const SubSpace<sdim> >
{
  static_assert((dim_ == 0 && sdim == 0) || (dim_ > 0 && sdim < dim_),
  "The dimensionality of the sub_grid is not valid.");


  auto sub_ref_basis = ref_basis_->get_ref_sub_space(s_id, dof_map, sub_grid);

  shared_ptr<const GridFunction<sdim,space_dim>> sub_func;

  const auto F = this->phys_domain_->get_grid_function();
  using IgFunc = IgGridFunction<dim,dim+codim>;
  auto ig_func = std::dynamic_pointer_cast<const IgFunc>(F);
  if (ig_func)
  {
    auto coeffs = ig_func->get_coefficients();
    IgCoefficients sub_coeffs;
    const int n_sub_dofs = dof_map.size();
    for (int sub_dof = 0 ; sub_dof < n_sub_dofs ; ++ sub_dof)
      sub_coeffs[sub_dof] = coeffs[dof_map[sub_dof]];

    sub_func = IgGridFunction<sdim,space_dim>::const_create(sub_ref_basis,sub_coeffs);
  }
  else
  {
    AssertThrow(false,ExcMessage("IgFunction"));
    AssertThrow(false,ExcNotImplemented());
  }
  auto sub_domain = Domain<sdim,space_dim-sdim>::const_create(sub_func);

  auto sub_phys_basis = SubSpace<sdim>::const_create(sub_ref_basis, sub_domain);

  return sub_phys_basis;
}



#if 0
template <int dim_, int range_, int rank_, int codim_>
auto
PhysicalSpaceBasis<dim_, range_, rank_, codim_>::
get_face_space(const Index face_id,
               SafeSTLVector<Index> &face_to_element_dofs) const -> shared_ptr<FaceSpace>
{
  auto elem_map = std::make_shared<typename GridType::FaceGridMap >();
  auto face_ref_sp = ref_basis_->get_ref_face_space(face_id, face_to_element_dofs, *elem_map);
  auto map  = push_forward_->get_mapping();

  auto fmap = MappingSlice<FaceSpace::PushForwardType::dim, FaceSpace::PushForwardType::codim>::
  create(map, face_id, face_ref_sp->get_grid(), elem_map);
  auto fpf = FaceSpace::PushForwardType::create(fmap);
  auto face_space = FaceSpace::create(face_ref_sp,fpf);

  return face_space;
}


template <int dim_, int range_, int rank_, int codim_>
Index
PhysicalSpaceBasis<dim_, range_, rank_, codim_>::
get_id() const
{
  return ref_basis_->get_id();
}
#endif


template <int dim_, int range_, int rank_, int codim_>
void
PhysicalSpaceBasis<dim_, range_, rank_, codim_>::
get_element_dofs(
  const IndexType element_id,
  SafeSTLVector<Index> &dofs_global,
  SafeSTLVector<Index> &dofs_local_to_patch,
  SafeSTLVector<Index> &dofs_local_to_elem,
  const std::string &dofs_property) const
{
  return ref_basis_->get_element_dofs(
           element_id,
           dofs_global,
           dofs_local_to_patch,
           dofs_local_to_elem,dofs_property);
}







template <int dim_, int range_, int rank_, int codim_>
void
PhysicalSpaceBasis<dim_, range_, rank_, codim_>::
print_info(LogStream &out) const
{
  out.begin_item("Reference space:");
  ref_basis_->print_info(out);
  out.end_item();

  out.begin_item("Physical domain:");
  phys_domain_->print_info(out);
  out.end_item();
}



template <int dim_, int range_, int rank_, int codim_>
auto
PhysicalSpaceBasis<dim_, range_, rank_, codim_>::
create_cache_handler() const
-> std::unique_ptr<SpaceElementHandler<dim_,codim_,range_,rank_>>
{
  return std::make_unique<ElementHandler>(this->get_this_basis());
}





template <int dim_, int range_, int rank_, int codim_>
Transformation
PhysicalSpaceBasis<dim_, range_, rank_, codim_>::
get_transformation_type() const
{
  return transformation_type_;
}


template <int dim_, int range_, int rank_, int codim_>
std::shared_ptr<const SplineSpace<dim_,range_,rank_>>
                                                   PhysicalSpaceBasis<dim_, range_, rank_, codim_>::
                                                   get_spline_space() const
{
  return ref_basis_->get_spline_space();
}


#ifdef MESH_REFINEMENT

template <int dim_, int range_, int rank_, int codim_>
void
PhysicalSpaceBasis<dim_, range_, rank_, codim_>::
refine_h(const Size n_subdivisions)
{
  //the refinement of the ReferenceSpaceBasis also refines the Domain (they share the same Grid)
  ref_basis_.get_ptr_data()->refine_h(n_subdivisions);
}


template <int dim_, int range_, int rank_, int codim_>
void
PhysicalSpaceBasis<dim_, range_, rank_, codim_>::
rebuild_after_insert_knots(
  const SafeSTLArray<SafeSTLVector<Real>,dim> &knots_to_insert,
  const Grid<dim> &old_grid)
{
  auto prev_ref_basis = ref_basis_->get_basis_previous_refinement();
  Assert(prev_ref_basis != nullptr, ExcNullPtr());

  auto prev_phys_domain = phys_domain_->get_domain_previous_refinement();
  Assert(prev_phys_domain != nullptr, ExcNullPtr());

  this->phys_basis_previous_refinement_ =
    self_t::const_create(prev_ref_basis,prev_phys_domain);
}

template <int dim_, int range_, int rank_, int codim_>
void
PhysicalSpaceBasis<dim_, range_, rank_, codim_>::
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
  std::const_pointer_cast<Grid<dim>>(ref_basis_->get_grid())->connect_insert_knots(SlotType(func_to_connect).track_foreign(space));
}

template <int dim_, int range_, int rank_, int codim_>
auto
PhysicalSpaceBasis<dim_, range_, rank_, codim_>::
get_basis_previous_refinement() const -> std::shared_ptr<const base_t>
{
  return phys_basis_previous_refinement_;
}

#endif

#if 0
#ifdef SERIALIZATION
template <int dim_, int range_, int rank_, int codim_>
template<class Archive>
void
PhysicalSpaceBasis<dim_, range_, rank_, codim_>::
serialize(Archive &ar, const unsigned int version)
{
  ar &boost::serialization::make_nvp("PhysicalSpaceBasis_base_t",
                                     boost::serialization::base_object<base_t>(*this));

  ar.template register_type<BSpline<dim_,range_,rank_> >();
#ifdef USE_NURBS
  ar.template register_type<NURBS<dim_,range_,rank_> >();
#endif // NURBS

  ar &boost::serialization::make_nvp("ref_basis_",ref_basis_);

  auto tmp = const_pointer_cast<self_t>(phys_space_previous_refinement_);
  ar &boost::serialization::make_nvp("phys_space_previous_refinement_",tmp);
  phys_space_previous_refinement_ = tmp;
}

///@}
#endif // SERIALIZATION
#endif

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/physical_space_basis.inst>


