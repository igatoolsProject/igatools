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

#include <igatools/geometry/physical_domain.h>
#include <igatools/geometry/physical_domain_element.h>
#include <igatools/functions/function.h>

using std::shared_ptr;

IGA_NAMESPACE_OPEN


//namespace
//{
//
//ValueFlags
//mapping_to_function_flags(const ValueFlags &flags)
//{
//    ValueFlags valid_func_flags = ValueFlags::value |
//                                  ValueFlags::gradient |
//                                  ValueFlags::hessian |
//                                  ValueFlags::divergence |
//                                  ValueFlags::point;
//
//    ValueFlags transfer_flags = ValueFlags::measure |
//                                ValueFlags::w_measure |
//                                ValueFlags::boundary_normal |
//                                valid_func_flags;
//
//
//    ValueFlags f_flags = flags & transfer_flags;
//
//    if (contains(flags, ValueFlags::measure) ||
//        contains(flags, ValueFlags::w_measure) ||
//        contains(flags, ValueFlags::inv_gradient) ||
//        contains(flags, ValueFlags::outer_normal))
//        f_flags |=  ValueFlags::gradient;
//
//    if (contains(flags, ValueFlags::inv_hessian) ||
//        contains(flags, ValueFlags::curvature))
//        f_flags |=  ValueFlags::gradient | ValueFlags::hessian;
//
//    return f_flags;
//}
//};



template<int dim_, int codim_>
PhysicalDomain<dim_, codim_>::
PhysicalDomain(std::shared_ptr<const GridType> grid,
               std::shared_ptr<FuncType> F)
  :
  grid_(grid),
  func_(F)
{
  Assert(grid_ != nullptr, ExcNullPtr());
}



template<int dim_, int codim_>
PhysicalDomain<dim_, codim_>::
~PhysicalDomain()
{}



//template<int dim_, int codim_>
//auto
//PhysicalDomain<dim_, codim_>::
//create(std::shared_ptr<FuncType> F)-> std::shared_ptr<self_t>
//{
//    return std::shared_ptr<self_t>(new self_t(F));
//}


template<int dim_, int codim_>
auto
PhysicalDomain<dim_, codim_>::
get_grid() const -> std::shared_ptr<const CartesianGrid<dim_> >
{
  return grid_;
}



template<int dim_, int codim_>
auto
PhysicalDomain<dim_, codim_>::
get_function() const -> std::shared_ptr<FuncType>
{
  return func_;
}



template<int dim_, int codim_>
auto
PhysicalDomain<dim_, codim_>::
create_cache_handler() const
-> std::shared_ptr<ElementHandler>
{
  return std::make_shared<ElementHandler>(this->shared_from_this());
}



template<int dim_, int codim_>
auto
PhysicalDomain<dim_, codim_>::
create_element(const ListIt &index, const PropId &prop) const
-> std::shared_ptr<ConstElementAccessor>
{
  using Elem = ConstElementAccessor;
  auto elem = std::make_shared<Elem>(this->shared_from_this(), index, prop);
  Assert(elem != nullptr,ExcNullPtr());

  return elem;
}



template<int dim_, int codim_>
auto
PhysicalDomain<dim_, codim_>::
begin(const PropId &prop) -> ElementIterator
{
  return ElementIterator(this->shared_from_this(),
  grid_->get_element_property(prop).begin(),
  prop);
}



template<int dim_, int codim_>
auto
PhysicalDomain<dim_, codim_>::
end(const PropId &prop) -> ElementIterator
{
  return ElementIterator(this->shared_from_this(),
  grid_->get_element_property(prop).end(),
  prop);
}



template<int dim_, int codim_>
auto
PhysicalDomain<dim_, codim_>::
begin(const PropId &prop) const -> ElementConstIterator
{
  return this->cbegin(prop);
}



template<int dim_, int codim_>
auto
PhysicalDomain<dim_, codim_>::
end(const PropId &prop) const -> ElementConstIterator
{
  return this->cend(prop);
}



template<int dim_, int codim_>
auto
PhysicalDomain<dim_, codim_>::
cbegin(const PropId &prop) const -> ElementConstIterator
{
  return ElementConstIterator(this->shared_from_this(),
                              grid_->get_element_property(prop).end(),
                              prop);
}



template<int dim_, int codim_>
auto
PhysicalDomain<dim_, codim_>::
cend(const PropId &prop) const -> ElementConstIterator
{
  return ElementConstIterator(this->shared_from_this(),
                              grid_->get_element_property(prop).end(),
                              prop);
}

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/physical_domain.inst>

