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

#include <igatools/basis_functions/physical_basis_element_handler.h>
#include <igatools/basis_functions/physical_basis_element.h>

#include <functional>

using std::shared_ptr;



IGA_NAMESPACE_OPEN

template<int dim, int codim>
using MapFunc= Function<dim, 0, dim + codim, 1>;


template<int dim_,int range_,int rank_,int codim_>
PhysicalBasisElementHandler<dim_,range_,rank_,codim_>::
PhysicalBasisElementHandler(std::shared_ptr<const PhysSpace> space)
  :
  base_t(space),
  ref_space_handler_(space->get_reference_basis()->create_cache_handler()),
  phys_domain_handler_(space->get_domain()->create_cache_handler()),
  phys_space_(space)
{}





template<int dim_,int range_,int rank_,int codim_>
auto
PhysicalBasisElementHandler<dim_,range_,rank_,codim_>::
print_info(LogStream &out) const -> void
{
  ref_space_handler_->print_info(out);
  //  PFCache::print_info(out);
}



template<int dim_,int range_,int rank_,int codim_>
PhysicalBasisElementHandler<dim_,range_,rank_,codim_>::
SetFlagsDispatcher::
SetFlagsDispatcher(const typename space_element::Flags phys_elem_flag,
                  const Transformation &transformation_type,
                  RefElemHandler &ref_space_handler,
                  PhysDomainHandler &phys_domain_handler,
                  SafeSTLArray<typename space_element::Flags, dim+1> &flags)
  :
  phys_elem_flag_(phys_elem_flag),
  transformation_type_(transformation_type),
  ref_space_handler_(ref_space_handler),
  phys_domain_handler_(phys_domain_handler),
  flags_(flags)
{}


template<int dim_,int range_,int rank_,int codim_>
template<int sdim>
void
PhysicalBasisElementHandler<dim_,range_,rank_,codim_>::
SetFlagsDispatcher::
operator()(const Topology<sdim> &topology)
{
  ref_space_handler_.template set_flags<sdim>(
    phys_space_to_reference_space_flag(transformation_type_,phys_elem_flag_));

  phys_domain_handler_.template set_flags<sdim>(
    phys_space_to_domain_flag(transformation_type_,phys_elem_flag_));

  flags_[sdim] = phys_elem_flag_;
}


template<int dim_,int range_,int rank_,int codim_>
PhysicalBasisElementHandler<dim_,range_,rank_,codim_>::
InitCacheDispatcher::
InitCacheDispatcher(const RefElemHandler &ref_space_handler,
                    const PhysDomainHandler &phys_domain_handler,
                    const SafeSTLArray<typename space_element::Flags, dim+1> &flags,
                    BaseElem &elem)
  :
  ref_space_handler_(ref_space_handler),
  phys_domain_handler_(phys_domain_handler),
  flags_(flags),
  elem_(elem)
{}


template<int dim_,int range_,int rank_,int codim_>
template<int sdim>
void
PhysicalBasisElementHandler<dim_,range_,rank_,codim_>::
InitCacheDispatcher::
operator()(const std::shared_ptr<const Quadrature<sdim>> &quad)
{
  using PhysSpaceElem = PhysicalBasisElement<dim_,range_,rank_,codim_>;
  auto &phys_space_elem  = dynamic_cast<PhysSpaceElem &>(elem_);

  ref_space_handler_.template init_cache<sdim>(
    *phys_space_elem.ref_space_element_,quad);

  phys_domain_handler_.init_cache(
    *phys_space_elem.phys_domain_element_,quad);


  auto &cache = phys_space_elem.all_sub_elems_cache_;

  const auto n_basis = phys_space_elem.get_num_basis(DofProperties::active);

  const auto n_pts = quad->get_num_points();

  const auto flag = flags_[sdim];

  for (auto &s_id: UnitElement<dim_>::template elems_ids<sdim>())
  {
    auto &s_cache = cache.template get_sub_elem_cache<sdim>(s_id);
    s_cache.resize(flag, n_pts, n_basis);
  }
}


template<int dim_,int range_,int rank_,int codim_>
PhysicalBasisElementHandler<dim_,range_,rank_,codim_>::
FillCacheDispatcher::
FillCacheDispatcher(const int s_id,
                    const RefElemHandler &ref_space_handler,
                    const PhysDomainHandler &phys_domain_handler,
                    const self_t &phys_space_handler,
                    BaseElem &elem)
  :
  s_id_(s_id),
  ref_space_handler_(ref_space_handler),
  phys_domain_handler_(phys_domain_handler),
  phys_space_handler_(phys_space_handler),
  elem_(elem)
{}


template<int dim_,int range_,int rank_,int codim_>
template<int sdim>
void
PhysicalBasisElementHandler<dim_,range_,rank_,codim_>::
FillCacheDispatcher::
operator()(const Topology<sdim> &topology)
{
  using PhysSpaceElem = PhysicalBasisElement<dim_,range_,rank_,codim_>;
  auto &phys_space_elem  = dynamic_cast<PhysSpaceElem &>(elem_);

  auto &ref_space_elem = *phys_space_elem.ref_space_element_;
  auto &phys_domain_elem = *phys_space_elem.phys_domain_element_;

  ref_space_handler_.template fill_cache<sdim>(ref_space_elem,s_id_);

  phys_domain_handler_.template fill_cache<sdim>(phys_domain_elem,s_id_);

  auto &all_sub_elems_cache = phys_space_elem.all_sub_elems_cache_;
  auto &sub_elem_cache = all_sub_elems_cache.template get_sub_elem_cache<sdim>(s_id_);

  using _Value = typename BaseElem::_Value;
  using _Gradient = typename BaseElem::_Gradient;
  using _Hessian = typename BaseElem::_Hessian;
  using _Divergence = typename BaseElem::_Divergence;

  const auto phys_space = phys_space_elem.get_physical_space();
  const typename PhysSpace::PushFwd push_fwd(phys_space->get_transformation_type());
  if (sub_elem_cache.template status_fill<_Value>())
  {
    push_fwd.template
    transform_0<RefBasis::range,RefBasis::rank,sdim>(
      s_id_,
      ref_space_elem,
      phys_domain_elem,
      sub_elem_cache);
  }
  if (sub_elem_cache.template status_fill<_Gradient>())
  {
    push_fwd.template
    transform_1<RefBasis::range,RefBasis::rank,sdim>(
      s_id_,ref_space_elem,phys_domain_elem,sub_elem_cache);
  }
  if (sub_elem_cache.template status_fill<_Hessian>())
  {
    push_fwd.template
    transform_2<RefBasis::range,RefBasis::rank,sdim>(
      s_id_,ref_space_elem,phys_domain_elem,sub_elem_cache);
  }
  if (sub_elem_cache.template status_fill<_Divergence>())
  {
    auto &divergences = sub_elem_cache.template get_data<_Divergence>();
    eval_divergences_from_gradients(
      sub_elem_cache.template get_data<_Gradient>(),
      divergences);

    divergences.set_status_filled(true);
  }

  sub_elem_cache.set_filled(true);

}


template<int dim_,int range_,int rank_,int codim_>
void
PhysicalBasisElementHandler<dim_,range_,rank_,codim_>::
set_flags_impl(const topology_variant &topology,
               const typename space_element::Flags &flag)
{
  auto set_flag_dispatcher = SetFlagsDispatcher(
                               flag,
                               phys_space_->get_transformation_type(),
                               *ref_space_handler_,
                               *phys_domain_handler_,
                               this->flags_);
  boost::apply_visitor(set_flag_dispatcher,topology);
}


template<int dim_,int range_,int rank_,int codim_>
void
PhysicalBasisElementHandler<dim_,range_,rank_,codim_>::
init_cache_impl(BaseElem &elem,
                const eval_pts_variant &quad) const
{
  auto init_cache_dispatcher =
    InitCacheDispatcher(*ref_space_handler_,*phys_domain_handler_,this->flags_,elem);
  boost::apply_visitor(init_cache_dispatcher,quad);
}


template<int dim_,int range_,int rank_,int codim_>
void
PhysicalBasisElementHandler<dim_,range_,rank_,codim_>::
fill_cache_impl(const topology_variant &topology,
                BaseElem &elem,
                const int s_id) const
{
  auto fill_cache_dispatcher =
    FillCacheDispatcher(s_id,*ref_space_handler_,*phys_domain_handler_,*this,elem);
  boost::apply_visitor(fill_cache_dispatcher,topology);
}

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/physical_basis_element_handler.inst>
