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

#include <igatools/basis_functions/physical_basis_element.h>
#include <igatools/base/exceptions.h>


using std::shared_ptr;
using std::make_shared;
using std::make_unique;

IGA_NAMESPACE_OPEN

template<int dim_,int range_,int rank_,int codim_>
PhysicalBasisElement<dim_,range_,rank_,codim_>::
PhysicalBasisElement(const std::shared_ptr<ContainerType> &phys_basis,
                     GridIterator<RefElemAccessor> &&ref_basis_element,
                     GridIterator<PhysDomainElem> &&phys_domain_element)
  :
  parent_t(phys_basis),
  ref_basis_element_(std::move(ref_basis_element)),
  phys_domain_element_(std::move(phys_domain_element)),
  phys_basis_(phys_basis)
{}

template<int dim_,int range_,int rank_,int codim_>
auto
PhysicalBasisElement<dim_,range_,rank_,codim_>::
get_physical_basis() const -> std::shared_ptr<const PhysBasis>
{
  return phys_basis_;
}

template<int dim_,int range_,int rank_,int codim_>
auto
PhysicalBasisElement<dim_,range_,rank_,codim_>::
get_grid_element() -> GridElement<dim_> &
{
  return get_domain_element().get_grid_function_element().get_grid_element();
}

template<int dim_,int range_,int rank_,int codim_>
auto
PhysicalBasisElement<dim_,range_,rank_,codim_>::
get_grid_element() const -> const GridElement<dim_> &
{
  return get_domain_element().get_grid_function_element().get_grid_element();
}


#if 0
template<int dim_,int range_,int rank_,int codim_>
PhysicalBasisElement<dim_,range_,rank_,codim_>::
PhysicalBasisElement(const PhysicalBasisElement<dim_,range_,rank_,codim_> &in,
                     const CopyPolicy &copy_policy)
  :
  parent_t(in,copy_policy)
{
  if (copy_policy == CopyPolicy::shallow)
  {
    ref_basis_element_ = in.ref_basis_element_;
    phys_domain_element_ = in.phys_domain_element_;
  }
  else
  {
    ref_basis_element_ = std::dynamic_pointer_cast<RefElemAccessor>(in.ref_basis_element_->clone());
    phys_domain_element_ = make_shared<PhysDomainElem>(*in.phys_domain_element_);
  }

  Assert(false,ExcNotTested());
}
#endif



template<int dim_,int range_,int rank_,int codim_>
void
PhysicalBasisElement<dim_,range_,rank_,codim_>::
operator++()
{
  ++(*phys_domain_element_);

  ++(*ref_basis_element_);
}

template<int dim_,int range_,int rank_,int codim_>
void
PhysicalBasisElement<dim_,range_,rank_,codim_>::
move_to(const IndexType &elem_id)
{
  ref_basis_element_->move_to(elem_id);
  phys_domain_element_->move_to(elem_id);
}



template<int dim_,int range_,int rank_,int codim_>
template <int sdim>
auto
PhysicalBasisElement<dim_,range_,rank_,codim_>::
get_points(const int s_id) const -> const ValueVector<PhysPoint>
{
//  using _Point = typename PhysDomainElem::_Point;
  return phys_domain_element_->template get_points<sdim>(s_id);
}

template<int dim_,int range_,int rank_,int codim_>
auto
PhysicalBasisElement<dim_,range_,rank_,codim_>::
get_element_points() const -> const ValueVector<PhysPoint>
{
  return this->template get_points<dim>(0);
}


template<int dim_,int range_,int rank_,int codim_>
auto
PhysicalBasisElement<dim_,range_,rank_,codim_>::
get_element_w_measures() const -> const ValueVector<Real>
{
  return this->template get_w_measures<dim>(0);
}


template<int dim_,int range_,int rank_,int codim_>
auto
PhysicalBasisElement<dim_,range_,rank_,codim_>::
get_index() const -> IndexType
{
  return parent_t::get_index();
}

#if 0
template<int dim_,int range_,int rank_,int codim_>
void
PhysicalBasisElement<dim_,range_,rank_,codim_>::
move_to(const Index flat_index)
{
  this->get_grid_element().move_to(flat_index);
  ref_basis_element_->move_to(flat_index);
  phys_domain_element_->move_to(flat_index);
}
#endif



template<int dim_,int range_,int rank_,int codim_>
auto
PhysicalBasisElement<dim_,range_,rank_,codim_>::
get_ref_basis_element() const -> const RefElemAccessor &
{
  return dynamic_cast<const RefElemAccessor &>(*ref_basis_element_);
}

template<int dim_,int range_,int rank_,int codim_>
auto
PhysicalBasisElement<dim_,range_,rank_,codim_>::
get_ref_basis_element() -> RefElemAccessor &
{
  return dynamic_cast<RefElemAccessor &>(*ref_basis_element_);
}

#if 0
template<int dim_,int range_,int rank_,int codim_>
auto
PhysicalBasisElement<dim_,range_,rank_,codim_>::
get_grid() const -> const std::shared_ptr<const Grid<dim> >
{
  return this->get_ref_basis_element().get_grid();
}
#endif

template<int dim_,int range_,int rank_,int codim_>
auto
PhysicalBasisElement<dim_,range_,rank_,codim_>::
get_domain_element() const -> const PhysDomainElem &
{
  return *phys_domain_element_;
}

template<int dim_,int range_,int rank_,int codim_>
auto
PhysicalBasisElement<dim_,range_,rank_,codim_>::
get_domain_element() -> PhysDomainElem &
{
  return *phys_domain_element_;
}
//*/


template<int dim_,int range_,int rank_,int codim_>
void
PhysicalBasisElement<dim_,range_,rank_,codim_>::
print_info(LogStream &out) const
{
  using std::to_string;
  out.begin_item("PhysicalBasisElement<" +
                 to_string(dim_) + "," +
                 to_string(range_) + "," +
                 to_string(rank_) + "," +
                 to_string(codim_) + ">:");

  out.begin_item("ReferenceBasisElement<" +
                 to_string(RefBasis::dim) + "," +
                 to_string(RefBasis::range) + "," +
                 to_string(RefBasis::rank) + ">");
  ref_basis_element_->print_info(out);
  out.end_item();

  out.begin_item("DomainElement<" +
                 to_string(dim_) + "," +
                 to_string(codim_) + ">");
  phys_domain_element_->print_info(out);
  out.end_item();

  std::string transf_type;
  auto type = phys_basis_->get_transformation_type();
  if (type == Transformation::h_grad)
    transf_type = "h_grad";
  else if (type == Transformation::h_div)
    transf_type = "h_div";
  else if (type == Transformation::h_curl)
    transf_type = "h_curl";
  else if (type == Transformation::l_2)
    transf_type = "l_2";
  out << "Transformation: " << transf_type << std::endl;

  out.end_item();
}

template<int dim_,int range_,int rank_,int codim_>
void
PhysicalBasisElement<dim_,range_,rank_,codim_>::
print_cache_info(LogStream &out) const
{
  using std::to_string;
  out.begin_item("ReferenceBasisElement<" +
                 to_string(RefBasis::dim) + "," +
                 to_string(RefBasis::range) + "," +
                 to_string(RefBasis::rank) + "> cache:");
  ref_basis_element_->print_cache_info(out);
  out.end_item();

  out.begin_item("DomainElement<" +
                 to_string(dim_) + "," +
                 to_string(codim_) + "> cache:");
  phys_domain_element_->print_cache_info(out);
  out.end_item();
}



template<int dim_,int range_,int rank_,int codim_>
bool
PhysicalBasisElement<dim_,range_,rank_,codim_>::
operator==(const parent_t &a) const
{
  return !((*this) != a);
}

template<int dim_,int range_,int rank_,int codim_>
bool
PhysicalBasisElement<dim_,range_,rank_,codim_>::
operator!=(const parent_t &a) const
{
  const self_t &elem = dynamic_cast<const self_t &>(a);
  Assert(this->is_comparable_with(elem),
         ExcMessage("The elements are not comparable"));
  return (ref_basis_element_ != elem.ref_basis_element_) ||
         (phys_domain_element_ != elem.phys_domain_element_);
}



template<int dim_,int range_,int rank_,int codim_>
bool
PhysicalBasisElement<dim_,range_,rank_,codim_>::
is_comparable_with(const self_t &elem) const
{
  return (this->get_basis() == elem.get_basis());
}


template<int dim_,int range_,int rank_,int codim_>
template <int k>
auto
PhysicalBasisElement<dim_,range_,rank_,codim_>::
get_w_measures(const int j) const -> const ValueVector<Real>
{
  return phys_domain_element_->template get_w_measures<k>(j);
}

template<int dim_,int range_,int rank_,int codim_>
template <int k>
auto
PhysicalBasisElement<dim_,range_,rank_,codim_>::
get_measures(const int j) const -> const ValueVector<Real> &
{
  return phys_domain_element_->template get_measures<k>(j);
}


template<int dim_,int range_,int rank_,int codim_>
template <int k>
auto
PhysicalBasisElement<dim_,range_,rank_,codim_>::
get_boundary_normals(const int s_id) const
-> const ValueVector<Points<space_dim> > &
{
  return phys_domain_element_->template get_boundary_normals<k>(s_id);
}


IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/physical_basis_element.inst>
