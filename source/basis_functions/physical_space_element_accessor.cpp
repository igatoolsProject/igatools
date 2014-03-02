//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2014  by the igatools authors (see authors.txt).
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

#include <igatools/basis_functions/physical_space_element_accessor.h>
#include <igatools/base/exceptions.h>

using std::array;
using std::vector;
using std::shared_ptr;

IGA_NAMESPACE_OPEN

template< class PhysSpace >
PhysicalSpaceElementAccessor<PhysSpace>::
PhysicalSpaceElementAccessor(const PhysSpace &phys_space, const Index index)
    :
    RefElemAccessor(*phys_space.get_reference_space(), index),
    PfElemAccessor(*phys_space.push_forward_, index),
    phys_space_(&phys_space)
{
    Assert(phys_space_ != nullptr, ExcNullPtr());
}



template< class PhysSpace >
int
PhysicalSpaceElementAccessor<PhysSpace>::
get_flat_index() const
{
    return RefElemAccessor::get_flat_index();
}



template< class PhysSpace >
void
PhysicalSpaceElementAccessor<PhysSpace>::
ElementValuesCache::
reset(const int n_basis_per_element,
      const QuadratureType &quad,
      const ValueFlags fill_flag)
{
    n_points_ = quad.get_num_points();

    if (contains(fill_flag , ValueFlags::value))
    {
        if (D0phi_.get_num_points() != n_points_ ||
            D0phi_.get_num_functions() != n_basis_per_element)
            D0phi_.resize(n_basis_per_element,n_points_);

        D0phi_.zero();
        fill_values_ = true;
    }
    else
    {
        D0phi_.clear();
        fill_values_ = false;
    }

    if (contains(fill_flag , ValueFlags::gradient))
    {
        if (D1phi_.get_num_points() != n_points_ ||
            D1phi_.get_num_functions() != n_basis_per_element)
            D1phi_.resize(n_basis_per_element,n_points_);

        D1phi_.zero();
        fill_gradients_ = true;
    }
    else
    {
        D1phi_.clear();
        fill_gradients_ = false;
    }

    if (contains(fill_flag , ValueFlags::hessian))
    {
        if (D2phi_.get_num_points() != n_points_ ||
            D2phi_.get_num_functions() != n_basis_per_element)
            D2phi_.resize(n_basis_per_element,n_points_);

        D2phi_.zero();
        fill_hessians_ = true;
    }
    else
    {
        D2phi_.clear();
        fill_hessians_ = false;
    }

    this->set_initialized(true);
}



template< class PhysSpace >
ValueFlags
PhysicalSpaceElementAccessor<PhysSpace>::
get_reference_space_accessor_fill_flags(const ValueFlags fill_flag) const
{
    bool fill_values = false;
    bool fill_gradients = false;
    bool fill_hessians = false;

    if (contains(fill_flag , ValueFlags::value))
        fill_values = true;

    if (contains(fill_flag , ValueFlags::gradient))
        fill_gradients = true;

    if (contains(fill_flag , ValueFlags::hessian))
        fill_hessians = true;


    bool fill_D0_phi_hat = false;
    bool fill_D1_phi_hat = false;
    bool fill_D2_phi_hat = false;
    if (transformation_type == Transformation::h_grad)
    {
        fill_D0_phi_hat = fill_values;
        fill_D1_phi_hat = fill_gradients || fill_hessians;
        fill_D2_phi_hat = fill_hessians;
    }
    else if (transformation_type == Transformation::h_div  ||
             transformation_type == Transformation::h_curl ||
             transformation_type == Transformation::l_2)
    {
        fill_D0_phi_hat = fill_values || fill_gradients || fill_hessians;
        fill_D1_phi_hat = fill_gradients || fill_hessians;
        fill_D2_phi_hat = fill_hessians;
    }


    ValueFlags reference_space_accessor_fill_flags = ValueFlags::none;
    if (fill_D0_phi_hat)
        reference_space_accessor_fill_flags |= ValueFlags::value;

    if (fill_D1_phi_hat)
        reference_space_accessor_fill_flags |= ValueFlags::gradient;

    if (fill_D2_phi_hat)
        reference_space_accessor_fill_flags |= ValueFlags::hessian;

    return reference_space_accessor_fill_flags;
}



template< class PhysSpace >
ValueFlags
PhysicalSpaceElementAccessor<PhysSpace>::
get_push_forward_accessor_fill_flags(const ValueFlags fill_flag) const
{
    const ValueFlags common_flag =
        ValueFlags::point|ValueFlags::map_value|ValueFlags::map_gradient|ValueFlags::map_hessian|
        ValueFlags::w_measure;

    ValueFlags pf_flags = fill_flag & common_flag;

    if (contains(fill_flag , ValueFlags::value))
        pf_flags |= ValueFlags::tran_value;

    if (contains(fill_flag , ValueFlags::gradient))
        pf_flags |= ValueFlags::tran_gradient;

    if (contains(fill_flag , ValueFlags::hessian))
        pf_flags |= ValueFlags::tran_hessian;

    return pf_flags;
}



template< class PhysSpace >
void
PhysicalSpaceElementAccessor<PhysSpace>::
init_values(const ValueFlags fill_flag,
            const QuadratureType &quad)
{
    const ValueFlags ref_sp_flag =
        get_reference_space_accessor_fill_flags(fill_flag);
    RefElemAccessor::init_values(ref_sp_flag, quad);

    const ValueFlags pf_flag = get_push_forward_accessor_fill_flags(fill_flag);
    PfElemAccessor::init_values(pf_flag, quad);


    const Size n_basis = phys_space_->get_reference_space()->get_num_basis_per_element();
    elem_values_.reset(n_basis, quad, fill_flag);
}



template< class PhysSpace >
void PhysicalSpaceElementAccessor<PhysSpace>::
fill_values()
{
    Assert(elem_values_.is_initialized(), ExcNotInitialized());

    PfElemAccessor::fill_values();
    RefElemAccessor::fill_values();

    if (elem_values_.fill_values_)
    {
        PfElemAccessor::template transform_values<RefSpace::dim_range,RefSpace::rank>
        (RefElemAccessor::get_basis_values(), elem_values_.D0phi_);
    }

    if (elem_values_.fill_gradients_)
    {
        if (transformation_type == Transformation::h_grad)
        {
            ValueTable<typename RefElemAccessor::Value> dummy;
            PfElemAccessor::
            template transform_gradients<PhysSpace::dim_range,PhysSpace::rank>(
                dummy,
                RefElemAccessor::get_basis_gradients(),
                elem_values_.D1phi_);
        }
        else
        {
            PfElemAccessor::
            template transform_gradients<PhysSpace::dim_range,PhysSpace::rank>(
                RefElemAccessor::get_basis_values(),
                RefElemAccessor::get_basis_gradients(),
                elem_values_.D1phi_);
        }
    }

    elem_values_.set_filled(true);
}



template< class PhysSpace >
const ValueVector<Real> &
PhysicalSpaceElementAccessor<PhysSpace>::
get_w_measures() const
{
    return PfElemAccessor::get_w_measures();
}



template< class PhysSpace >
auto
PhysicalSpaceElementAccessor<PhysSpace>::
get_point(const Index qp) const -> const Point<space_dim> &
{
    Assert(elem_values_.is_filled(), ExcCacheNotFilled());
    return (PfElemAccessor::get_values_map())[qp];
}



template< class PhysSpace >
auto
PhysicalSpaceElementAccessor<PhysSpace>::
evaluate_field(const Vector &coefs) const -> ValueVector< Value >
{
    Assert(phys_space_->get_reference_space()->get_num_basis() == coefs.size(),
    ExcDimensionMismatch(phys_space_->get_reference_space()->get_num_basis(), coefs.size()));

    auto field_hat = RefElemAccessor::evaluate_field(coefs);

    ValueVector< Value > field(field_hat.size());

    PfElemAccessor::template
    transform_values<PhysSpace::dim_range,PhysSpace::rank>(field_hat, field);

    return field;
}



template< class PhysSpace >
auto
PhysicalSpaceElementAccessor<PhysSpace>::
evaluate_field_gradients(const Vector &coefs) const -> ValueVector< Derivative<1> >
{
    Assert(phys_space_->get_reference_space()->get_num_basis() == coefs.size(),
    ExcDimensionMismatch(phys_space_->get_reference_space()->get_num_basis(), coefs.size()));

    auto D1field_hat = RefElemAccessor::evaluate_field_gradients(coefs);

    const auto n_quad_points = D1field_hat.size();

    ValueVector< typename RefElemAccessor::Value > D0field_hat(n_quad_points);
    if (transformation_type != Transformation::h_grad)
        D0field_hat = RefElemAccessor::evaluate_field(coefs);

    ValueVector< Derivative<1> > D1field(n_quad_points);
    PfElemAccessor::
    template transform_gradients<PhysSpace::dim_range,PhysSpace::rank>(
        D0field_hat, D1field_hat, D1field);

    return D1field;
}



template< class PhysSpace >
auto
PhysicalSpaceElementAccessor<PhysSpace>::
evaluate_field_hessians(const Vector &coefs) const -> ValueVector< Derivative<2> >
{
    AssertThrow(false,ExcNotImplemented());

    Assert(elem_values_.is_filled(), ExcCacheNotFilled());
    Assert(elem_values_.fill_hessians_ == true, ExcInvalidState());


    ValueVector< Derivative<2> > D2field;

    return D2field;
}



template< class PhysSpace >
auto
PhysicalSpaceElementAccessor<PhysSpace>::
get_basis_values() const -> ValueTable<Value> const &
{
    Assert(elem_values_.is_filled(), ExcCacheNotFilled());
    Assert(elem_values_.fill_values_,ExcInvalidState());
    return elem_values_.D0phi_;
}



template< class PhysSpace >
auto
PhysicalSpaceElementAccessor<PhysSpace>::
get_basis_values(int iFunc) const -> typename ValueTable<Value>::const_function_view
{
    return this->get_basis_values().get_function(iFunc);
}



template< class PhysSpace >
auto
PhysicalSpaceElementAccessor<PhysSpace>::
get_basis_value(const Index func, const Index qp) const -> const Value &
{
    Assert(qp >= 0 && qp < elem_values_.n_points_,
           ExcIndexRange(qp,0,elem_values_.n_points_));

    return this->get_basis_values(func)[qp];
}



template< class PhysSpace >
Real
PhysicalSpaceElementAccessor<PhysSpace>::
get_basis_divergence(const Index func, const Index qp) const
{
    return (trace(this->get_basis_gradient(func,qp)));
}



template< class PhysSpace >
auto
PhysicalSpaceElementAccessor<PhysSpace>::
get_basis_gradients() const -> ValueTable< Derivative<1> > const &
{
    Assert(elem_values_.is_filled(), ExcCacheNotFilled());
    Assert(elem_values_.fill_gradients_,ExcInvalidState());
    return elem_values_.D1phi_;
}



template< class PhysSpace >
auto
PhysicalSpaceElementAccessor<PhysSpace>::
get_basis_gradients(int iFunc) const -> typename ValueTable< Derivative<1> >::const_function_view
{
    return this->get_basis_gradients().get_function(iFunc);
}



template< class PhysSpace >
auto
PhysicalSpaceElementAccessor<PhysSpace>::
get_basis_gradient(const Index func, const Index qp) const -> const Derivative<1> &
{
    Assert(qp >= 0 && qp < elem_values_.n_points_,
    ExcIndexRange(qp,0,elem_values_.n_points_));

    return this->get_basis_gradients(func)[qp];
}



template< class PhysSpace >
auto
PhysicalSpaceElementAccessor<PhysSpace>::
get_basis_hessians() const -> ValueTable< Derivative<2> > const &
{
    Assert(elem_values_.is_filled(), ExcCacheNotFilled());
    Assert(elem_values_.fill_hessians_,ExcInvalidState());
    return elem_values_.D2phi_;
}



template< class PhysSpace >
auto
PhysicalSpaceElementAccessor<PhysSpace>::
get_basis_hessians(int iFunc) const -> typename ValueTable< Derivative<2> >::const_function_view
{
    return this->get_basis_hessians().get_function(iFunc);
}



template< class PhysSpace >
auto
PhysicalSpaceElementAccessor<PhysSpace>::
get_basis_hessian(const Index func, const Index qp) const -> const Derivative<2> &
{
    Assert(qp >= 0 && qp < elem_values_.n_points_,
    ExcIndexRange(qp,0,elem_values_.n_points_));

    return this->get_basis_hessians(func)[qp];
}



template< class PhysSpace >
auto
PhysicalSpaceElementAccessor<PhysSpace>::
get_points() const ->
const ValueVector< Point<space_dim> > &
{
    Assert(elem_values_.is_filled(), ExcCacheNotFilled());
    return PfElemAccessor::get_values_map();
}



template< class PhysSpace >
auto
PhysicalSpaceElementAccessor<PhysSpace>::
get_map_gradient_at_points() const ->
const ValueVector< typename Mapping<dim, codim>::GradientType > &
{
    Assert(elem_values_.is_filled(), ExcCacheNotFilled());
    return PfElemAccessor::get_gradients_map();
}



template< class PhysSpace >
bool
PhysicalSpaceElementAccessor<PhysSpace>::
is_boundary() const
{
    return PfElemAccessor::is_boundary();
}


template< class PhysSpace >
bool
PhysicalSpaceElementAccessor<PhysSpace>::
is_boundary(int face) const
{
    return PfElemAccessor::is_boundary(face);
}



template< class PhysSpace >
const PhysSpace *
PhysicalSpaceElementAccessor<PhysSpace>::
get_physical_space() const
{
    return phys_space_;
}



template< class PhysSpace >
void
PhysicalSpaceElementAccessor<PhysSpace>::
operator++()
{
    RefElemAccessor::operator++();
    PfElemAccessor::operator++();
}



template< class PhysSpace >
bool
PhysicalSpaceElementAccessor<PhysSpace>::
operator==(const PhysicalSpaceElementAccessor <PhysSpace> &a) const
{
    return RefElemAccessor::get_flat_index() ==
           static_cast<RefElemAccessor>(a).get_flat_index();
}



template< class PhysSpace >
bool
PhysicalSpaceElementAccessor<PhysSpace>::
operator!=(const PhysicalSpaceElementAccessor <PhysSpace> &a) const
{
    return RefElemAccessor::get_flat_index() !=
           static_cast<RefElemAccessor>(a).get_flat_index();
}

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/physical_space_element_accessor.inst>
