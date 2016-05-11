//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2016  by the igatools authors (see authors.txt).
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

#include <igatools/basis_functions/physical_basis.h>
#include <igatools/functions/formula_grid_function.h>
#include <igatools/basis_functions/physical_basis_handler.h>
#include <igatools/geometry/push_forward.h>


using std::shared_ptr;
using std::make_shared;
using std::unique_ptr;

using std::const_pointer_cast;
using std::endl;

IGA_NAMESPACE_OPEN


template <int dim_, int range_, int rank_, int codim_>
const SafeSTLArray<int, PhysicalBasis<dim_, range_, rank_, codim_>::n_components>
PhysicalBasis<dim_, range_, rank_, codim_>::components =
  sequence<PhysicalBasis<dim_, range_, rank_, codim_>::n_components>();


template <int dim_, int range_, int rank_, int codim_>
PhysicalBasis<dim_, range_, rank_, codim_>::
PhysicalBasis(const SharedPtrConstnessHandler<RefBasis> &ref_basis,
              const SharedPtrConstnessHandler<PhysDomain> &phys_domain,
              const Transformation &transformation_type)
  :
  ref_basis_(ref_basis),
  phys_domain_(phys_domain),
  transformation_type_(transformation_type)
{
  Assert(this->get_grid() == phys_domain_->get_grid_function()->get_grid(),
         ExcMessage("The basis and the physical domain must have the same grid!"));

//TODO(pauletti, Jan 18, 2014): put static assert on h_div, h_curl range and rank
}



template <int dim_, int range_, int rank_, int codim_>
auto
PhysicalBasis<dim_, range_, rank_, codim_>::
create(const shared_ptr<RefBasis> &ref_basis,
       const shared_ptr<PhysDomain> &phys_domain,
       const Transformation &transformation_type) -> shared_ptr<self_t>
{
  auto sp = shared_ptr<self_t>(
    new self_t(SharedPtrConstnessHandler<RefBasis>(ref_basis),
  SharedPtrConstnessHandler<PhysDomain>(phys_domain),
  transformation_type));
  Assert(sp != nullptr,ExcNullPtr());

#ifdef IGATOOLS_WITH_MESH_REFINEMENT
  sp->create_connection_for_insert_knots(sp);
#endif // IGATOOLS_WITH_MESH_REFINEMENT

  return sp;
}

template <int dim_, int range_, int rank_, int codim_>
auto
PhysicalBasis<dim_, range_, rank_, codim_>::
const_create(const shared_ptr<const RefBasis> &ref_basis,
             const shared_ptr<const PhysDomain> &phys_domain,
             const Transformation &transformation_type) -> shared_ptr<const self_t>
{
  auto sp = shared_ptr<self_t>(
    new self_t(SharedPtrConstnessHandler<RefBasis>(ref_basis),
  SharedPtrConstnessHandler<PhysDomain>(phys_domain),transformation_type));
  Assert(sp != nullptr,ExcNullPtr());

  return sp;
}

template <int dim_, int range_, int rank_, int codim_>
auto
PhysicalBasis<dim_, range_, rank_, codim_>::
get_this_basis() const -> std::shared_ptr<const self_t >
{
  auto sp = const_cast<self_t *>(this)->shared_from_this();
  auto this_basis = std::dynamic_pointer_cast<self_t>(sp);
  Assert(this_basis != nullptr,ExcNullPtr());

  return this_basis;
}

#if 0
template <int dim_, int range_, int rank_, int codim_>
auto
PhysicalBasis<dim_, range_, rank_, codim_>::
create_element(const ListIt &index, const PropId &prop) const
-> std::unique_ptr<BasisElement<dim_,codim_,range_,rank_>>
{
  return std::unique_ptr<ElementAccessor>(
    new ElementAccessor(this->get_this_basis(),
  index,
  this->get_reference_basis()->create_ref_element(index,prop),
  this->get_domain()->create_element(index,prop),
  prop));
}
#endif

template <int dim_, int range_, int rank_, int codim_>
auto
PhysicalBasis<dim_, range_, rank_, codim_>::
create_element_begin(const PropId &prop) const
-> std::unique_ptr<BasisElement<dim_,codim_,range_,rank_>>
{
  return std::unique_ptr<ElementAccessor>(
    new ElementAccessor(this->get_this_basis(),
  this->get_reference_basis()->create_ref_element_begin(prop),
  this->get_domain()->create_element_begin(prop)));
}

template <int dim_, int range_, int rank_, int codim_>
auto
PhysicalBasis<dim_, range_, rank_, codim_>::
create_element_end(const PropId &prop) const
-> std::unique_ptr<BasisElement<dim_,codim_,range_,rank_>>
{
  return std::unique_ptr<ElementAccessor>(
    new ElementAccessor(this->get_this_basis(),
  this->get_reference_basis()->create_ref_element_end(prop),
  this->get_domain()->create_element_end(prop)));
}


template <int dim_, int range_, int rank_, int codim_>
auto
PhysicalBasis<dim_, range_, rank_, codim_>::
get_reference_basis() const -> shared_ptr<const RefBasis>
{
  return ref_basis_.get_ptr_const_data();
}



template <int dim_, int range_, int rank_, int codim_>
auto
PhysicalBasis<dim_, range_, rank_, codim_>::
get_reference_basis() -> shared_ptr<RefBasis>
{
  return ref_basis_.get_ptr_data();
}



template <int dim_, int range_, int rank_, int codim_>
auto
PhysicalBasis<dim_, range_, rank_, codim_>::
get_domain() const -> std::shared_ptr<const Domain<dim_,codim_>>
{
  return phys_domain_.get_ptr_const_data();
}



template <int dim_, int range_, int rank_, int codim_>
auto
PhysicalBasis<dim_, range_, rank_, codim_>::
get_domain() -> std::shared_ptr<Domain<dim_,codim_>>
{
  return phys_domain_.get_ptr_data();
}



template <int dim_, int range_, int rank_, int codim_>
auto
PhysicalBasis<dim_, range_, rank_, codim_>::
get_grid() const -> std::shared_ptr<const Grid<dim_>>
{
  return ref_basis_->get_grid();
}


#if 0
template <int dim_, int range_, int rank_, int codim_>
auto
PhysicalBasis<dim_, range_, rank_, codim_>::
get_reference_basis() -> shared_ptr<RefBasis>
{
  return ref_basis_.get_ptr_data();
}
#endif


template <int dim_, int range_, int rank_, int codim_>
template<int sdim>
auto
PhysicalBasis<dim_, range_, rank_, codim_>::
get_sub_basis(const int s_id, InterBasisMap<sdim> &dof_map,
              const std::shared_ptr<const Grid<sdim>> &sub_grid,
              SubGridMap<sdim> &elem_map,
              EnableIf<(dim_ != 0) &&(sdim>=0)> *) const
-> std::shared_ptr<const SubBasis<sdim> >
{
  static_assert((dim_ == 0 && sdim == 0) || (dim_ > 0 && sdim < dim_),
  "The dimensionality of the sub_grid is not valid.");


  const auto sub_ref_basis = ref_basis_->get_ref_sub_basis(s_id, dof_map, sub_grid);

  const auto sub_domain = this->phys_domain_->get_sub_domain(s_id,elem_map,sub_grid);


  auto sub_phys_basis = SubBasis<sdim>::const_create(sub_ref_basis, sub_domain);

  return sub_phys_basis;
}



#if 0
template <int dim_, int range_, int rank_, int codim_>
auto
PhysicalBasis<dim_, range_, rank_, codim_>::
get_face_basis(const Index face_id,
               SafeSTLVector<Index> &face_to_element_dofs) const -> shared_ptr<FaceBasis>
{
  auto elem_map = std::make_shared<typename GridType::FaceGridMap >();
  auto face_ref_bs = ref_basis_->get_ref_face_basis(face_id, face_to_element_dofs, *elem_map);
  auto map  = push_forward_->get_mapping();

  auto fmap = MappingSlice<FaceBasis::PushForwardType::dim, FaceBasis::PushForwardType::codim>::
  create(map, face_id, face_ref_bs->get_grid(), elem_map);
  auto fpf = FaceBasis::PushForwardType::create(fmap);
  auto face_basis = FaceBasis::create(face_ref_bs,fpf);

  return face_basis;
}


template <int dim_, int range_, int rank_, int codim_>
Index
PhysicalBasis<dim_, range_, rank_, codim_>::
get_id() const
{
  return ref_basis_->get_id();
}
#endif








template <int dim_, int range_, int rank_, int codim_>
void
PhysicalBasis<dim_, range_, rank_, codim_>::
print_info(LogStream &out) const
{
  out.begin_item("Reference basis:");
  ref_basis_->print_info(out);
  out.end_item();

  out.begin_item("Physical domain:");
  phys_domain_->print_info(out);
  out.end_item();
}



template <int dim_, int range_, int rank_, int codim_>
auto
PhysicalBasis<dim_, range_, rank_, codim_>::
create_cache_handler() const
-> std::unique_ptr<BasisHandler<dim_,codim_,range_,rank_>>
{
  return std::unique_ptr<Handler>(
    new Handler(this->get_this_basis()));
}





template <int dim_, int range_, int rank_, int codim_>
Transformation
PhysicalBasis<dim_, range_, rank_, codim_>::
get_transformation_type() const
{
  return transformation_type_;
}


template <int dim_, int range_, int rank_, int codim_>
std::shared_ptr<const SplineSpace<dim_,range_,rank_> >
PhysicalBasis<dim_, range_, rank_, codim_>::get_spline_space() const
{
  return ref_basis_->get_spline_space();
}



template <int dim_, int range_, int rank_, int codim_>
std::shared_ptr<SplineSpace<dim_,range_,rank_> >
PhysicalBasis<dim_, range_, rank_, codim_>::get_spline_space()
{
  return ref_basis_.get_ptr_data()->get_spline_space();
}



#ifdef IGATOOLS_WITH_MESH_REFINEMENT
template <int dim_, int range_, int rank_, int codim_>
void
PhysicalBasis<dim_, range_, rank_, codim_>::refine_h(const Size n_subdivisions)
{
  //the refinement of the ReferenceBasis also refines the Domain (they share the same Grid)
  ref_basis_.get_ptr_data()->refine_h(n_subdivisions);
}


template <int dim_, int range_, int rank_, int codim_>
void
PhysicalBasis<dim_, range_, rank_, codim_>::rebuild_after_insert_knots(
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
PhysicalBasis<dim_, range_, rank_, codim_>::
create_connection_for_insert_knots(const std::shared_ptr<self_t> &basis)
{
  Assert(basis != nullptr, ExcNullPtr());
  Assert(&(*basis) == &(*this), ExcMessage("Different objects."));

  auto func_to_connect =
    std::bind(&self_t::rebuild_after_insert_knots,
              basis.get(),
              std::placeholders::_1,
              std::placeholders::_2);

  using SlotType = typename Grid<dim>::SignalInsertKnotsSlot;
  std::const_pointer_cast<Grid<dim>>(ref_basis_->get_grid())->connect_insert_knots(SlotType(func_to_connect).track_foreign(basis));
}

template <int dim_, int range_, int rank_, int codim_>
auto
PhysicalBasis<dim_, range_, rank_, codim_>::
get_basis_previous_refinement() const -> std::shared_ptr<const base_t>
{
  return phys_basis_previous_refinement_;
}

#endif

#ifdef IGATOOLS_WITH_SERIALIZATION

template <int dim_, int range_, int rank_, int codim_>
template<class Archive>
void
PhysicalBasis<dim_, range_, rank_, codim_>::
serialize(Archive &ar)
{
  using std::to_string;
  const std::string base_name = "Basis_" +
                                to_string(dim_) + "_" +
                                to_string(codim_) + "_" +
                                to_string(range_) + "_" +
                                to_string(rank_);

  ar &make_nvp(base_name,base_class<base_t>(this));


  ar &make_nvp("ref_basis_",ref_basis_);

  ar &make_nvp("phys_domain_",phys_domain_);

  Transformation transformation_type_tmp = transformation_type_;
  ar &make_nvp("transformation_type_",transformation_type_tmp);
  const_cast<Transformation &>(transformation_type_) = transformation_type_tmp;

#ifdef IGATOOLS_WITH_MESH_REFINEMENT
  auto tmp = std::const_pointer_cast<self_t>(phys_basis_previous_refinement_);
  ar &make_nvp("phys_basis_previous_refinement_",tmp);
  phys_basis_previous_refinement_ = tmp;
#endif // IGATOOLS_WITH_MESH_REFINEMENT
}

#endif // IGATOOLS_WITH_SERIALIZATION

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/physical_basis.inst>


