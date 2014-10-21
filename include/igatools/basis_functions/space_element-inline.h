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


#ifndef SPACE_ELEMENT_INLINE_H_
#define SPACE_ELEMENT_INLINE_H_

// TODO (pauletti, Oct 20, 2014): the inline name is missleading,
// there is no need for inlining, it is only ignorance on how to
// provide proper instantiation, needs to be improved

#include <igatools/basis_functions/space_element_accessor.h>

IGA_NAMESPACE_OPEN

template<class Space>
inline
SpaceElement<Space>::
SpaceElement(const std::shared_ptr<const Space> space,
                     const Index elem_index)
    :
    base_t(space->get_grid(), elem_index),
    space_(space),
    n_basis_direction_(space->get_num_all_element_basis())
{
    Assert(space_ != nullptr, ExcNullPtr());

    using Indexer = CartesianProductIndexer<dim>;
    for (int comp_id : basis_functions_indexer_.get_active_components_id())
    {
        // creating the objects for fast conversion from flat-to-tensor indexing
        // (in practice it is an hash-table from flat to tensor indices)
        basis_functions_indexer_[comp_id] =
            std::shared_ptr<Indexer>(new Indexer(n_basis_direction_[comp_id]));
    }

    comp_offset_[0] = 0;
    for (int comp_id = 1; comp_id < Space::n_components; ++comp_id)
        comp_offset_[comp_id] = comp_offset_[comp_id-1] +
        n_basis_direction_.comp_dimension[comp_id-1];
}



template<class Space>
inline
SpaceElement<Space>::
SpaceElement(const std::shared_ptr<const Space> space,
             const TensorIndex<dim> &elem_index)
    :
    SpaceElement(space->get_grid(), space->get_grid()->tensor_to_flat(elem_index))
{}



template<class Space>
inline
SpaceElement<Space>::
SpaceElement(const SpaceElement<Space> &elem,
                     const CopyPolicy &copy_policy)
    :
    CartesianGridElement<Space::dim>(elem,copy_policy),
    space_(elem.space_),
    n_basis_direction_(elem.n_basis_direction_),
    basis_functions_indexer_(elem.basis_functions_indexer_),
    comp_offset_(elem.comp_offset_)
{
    if (elem.local_cache_ != nullptr)
    {
        if (copy_policy == CopyPolicy::shallow)
        {
            local_cache_ = elem.local_cache_;
        }
        else
        {
            local_cache_ = std::shared_ptr<LocalCache>(new LocalCache(*elem.local_cache_));
        }
    }
}



template<class Space>
void
SpaceElement<Space>::
copy_from(const SpaceElement<Space> &elem,
          const CopyPolicy &copy_policy)
{
    CartesianGridElement<Space::dim>::copy_from(elem,copy_policy);
    if (this != &elem)
    {
        space_ = elem.space_;
        n_basis_direction_ = elem.n_basis_direction_;
        basis_functions_indexer_ = elem.basis_functions_indexer_;
        comp_offset_ = elem.comp_offset_;

        if (copy_policy == CopyPolicy::deep)
        {
            Assert(elem.local_cache_ != nullptr, ExcNullPtr());
            local_cache_ = std::shared_ptr<LocalCache>(new LocalCache(*elem.local_cache_));
        }
        else if (copy_policy == CopyPolicy::shallow)
        {
            local_cache_ = elem.local_cache_;
        }
        else
        {
            Assert(false,ExcNotImplemented());
            AssertThrow(false,ExcNotImplemented());
        }
    }
}



template<class Space>
void
SpaceElement<Space>::
deep_copy_from(const SpaceElement<Space> &elem)
{
    this->copy_from(elem,CopyPolicy::deep);
}



template<class Space>
void
SpaceElement<Space>::
shallow_copy_from(const SpaceElement<Space> &elem)
{
    this->copy_from(elem,CopyPolicy::shallow);
}



template<class Space>
SpaceElement<Space> &
SpaceElement<Space>::
operator=(const SpaceElement<Space> &element)
{
    this->shallow_copy_from(element);
    return (*this);
}



template<class Space>
inline
auto
SpaceElement<Space>::
as_cartesian_grid_element_accessor() -> CartesianGridElement<dim> &
{
    return static_cast<CartesianGridElement<dim> &>(*this);
}



template<class Space>
inline
auto
SpaceElement<Space>::
as_cartesian_grid_element_accessor() const -> const CartesianGridElement<dim> &
{
    return static_cast<const CartesianGridElement<dim> &>(*this);
}



template<class Space>
inline
auto
SpaceElement<Space>::
as_derived_element_accessor() -> DerivedElementAccessor &
{
    return static_cast<DerivedElementAccessor &>(*this);
}



template<class Space>
inline
auto
SpaceElement<Space>::
as_derived_element_accessor() const -> const DerivedElementAccessor &
{
    return static_cast<const DerivedElementAccessor &>(*this);
}



//template<class Space>
//template<int skel_dim, int der_order>
//auto
//SpaceElement<Space>::
//get_basis_ders(const int j) const
//{
//    const auto &cache = local_cache_->template get_value_cache<skel_dim>(j);
//    Assert(cache.is_filled() == true, ExcCacheNotFilled());
//   // Assert(cache.flags_handler_.values_filled(), ExcCacheNotFilled());
//
//    return cache.template get_der<der_order>();
//}


#if 0
template<class Space>
inline
auto
SpaceElement<Space>::
get_face_basis_values(const Index face_id) const
-> ValueTable<Value> const &
{
    return this->get_basis_values(FaceTopology<dim>(face_id));
}



template<class Space>
inline
auto
SpaceElement<Space>::
get_basis_values(const Index i,) const
-> typename ValueTable<Value>::const_view
{
    return this->get_basis_values(topology_id).get_function_view(i);
}


template<class Space>
inline
auto
SpaceElement<Space>::
get_basis_value(const Index basis, const Index qp,) const -> Value const &
{
    Assert(qp >= 0 && qp < this->get_values_cache(topology_id).quad_.get_num_points_direction().flat_size(),
           ExcIndexRange(qp,0,this->get_values_cache(topology_id).quad_.get_num_points_direction().flat_size()));
    return this->get_basis_values(basis,topology_id)[qp];
}




template<class Space>
inline
auto
SpaceElement<Space>::
get_basis_gradients() const -> ValueTable<Derivative<1>> const &
{
    const auto &cache = this->get_values_cache(topology_id);
    Assert(cache.is_filled() == true, ExcCacheNotFilled());
    Assert(cache.flags_handler_.gradients_filled(), ExcCacheNotFilled());

    return cache.D1phi_;
}

template<class Space>
inline
auto
SpaceElement<Space>::
get_face_basis_gradients(const Index face_id) const -> ValueTable<Derivative<1>> const &
{
    return this->get_basis_gradients(FaceTopology<dim>(face_id));
}



template<class Space>
inline
auto
SpaceElement<Space>::
get_basis_gradients(const Index i,) const -> typename ValueTable<Derivative<1>>::const_view
{
    return this->get_basis_gradients(topology_id).get_function_view(i);
}

template<class Space>
inline
auto
SpaceElement<Space>::
get_basis_gradient(const Index basis, const Index qp,) const -> Derivative<1> const &
{
    Assert(qp >= 0 && qp < this->get_values_cache(topology_id).quad_.get_num_points_direction().flat_size(),
           ExcIndexRange(qp,0,this->get_values_cache(topology_id).quad_.get_num_points_direction().flat_size()));
    return this->get_basis_gradients(basis,topology_id)[qp];
}





template<class Space>
inline
auto
SpaceElement<Space>::
get_basis_hessians() const -> ValueTable<Derivative<2>> const &
{
    const auto &cache = this->get_values_cache(topology_id);
    Assert(cache.is_filled() == true, ExcCacheNotFilled());
    Assert(cache.flags_handler_.hessians_filled(), ExcCacheNotFilled());

    return cache.D2phi_;
}

template<class Space>
inline
auto
SpaceElement<Space>::
get_face_basis_hessians(const Index face_id) const -> ValueTable<Derivative<2>> const &
{
    return this->get_basis_hessians(FaceTopology<dim>(face_id));
}



template<class Space>
inline
auto
SpaceElement<Space>::
get_basis_hessians(const Index i,) const -> typename ValueTable<Derivative<2>>::const_view
{
    return this->get_basis_hessians(topology_id).get_function_view(i);
}

template<class Space>
inline
auto
SpaceElement<Space>::
get_basis_hessian(const Index basis, const Index qp,) const -> Derivative<2> const &
{
    Assert(qp >= 0 && qp < this->get_values_cache(topology_id).quad_.get_num_points_direction().flat_size(),
           ExcIndexRange(qp,0,this->get_values_cache(topology_id).quad_.get_num_points_direction().flat_size()));
    return this->get_basis_hessians(basis,topology_id)[qp];
}




template<class Space>
inline
auto
SpaceElement<Space>::
get_basis_divergences() const -> ValueTable<Div> const &
{
    const auto &cache = this->get_values_cache(topology_id);
    Assert(cache.is_filled() == true, ExcCacheNotFilled());
    Assert(cache.flags_handler_.divergences_filled(), ExcCacheNotFilled());

    return cache.div_phi_;
}

template<class Space>
inline
auto
SpaceElement<Space>::
get_face_basis_divergences(const Index face_id) const -> ValueTable<Div> const &
{
    return this->get_basis_divergences(FaceTopology<dim>(face_id));
}



template<class Space>
inline
auto
SpaceElement<Space>::
get_basis_divergences(const Index i,) const -> typename ValueTable<Div>::const_view
{
    return this->get_basis_divergences(topology_id).get_function_view(i);
}

template<class Space>
inline
auto
SpaceElement<Space>::
get_basis_divergence(const Index basis, const Index qp,) const -> Div const &
{
    Assert(qp >= 0 && qp < this->get_values_cache(topology_id).quad_.get_num_points_direction().flat_size(),
           ExcIndexRange(qp,0,this->get_values_cache(topology_id).quad_.get_num_points_direction().flat_size()));
    return this->get_basis_divergences(basis,topology_id)[qp];
}






template<class Space>
inline
auto
SpaceElement<Space>::
evaluate_field(const vector<Real> &local_coefs,) const
-> ValueVector<Value>
{
    Assert(this->get_values_cache(topology_id).is_filled() == true, ExcCacheNotFilled());
    Assert(this->get_values_cache(topology_id).flags_handler_.fill_values() == true, ExcCacheNotFilled());
    Assert(this->get_num_basis() == local_coefs.size(),
    ExcDimensionMismatch(this->get_num_basis(),local_coefs.size()));

    const auto &D0phi_hat = this->get_basis_values(topology_id) ;
    Assert(D0phi_hat.get_num_functions() == this->get_num_basis(),
    ExcDimensionMismatch(D0phi_hat.get_num_functions(), this->get_num_basis())) ;

    return D0phi_hat.evaluate_linear_combination(local_coefs) ;
}

template<class Space>
inline
auto
SpaceElement<Space>::
evaluate_face_field(const Index face_id, const vector<Real> &local_coefs) const
-> ValueVector<Value>
{
    return this->evaluate_field(local_coefs,FaceTopology<dim>(face_id));
}


template<class Space>
inline
auto
SpaceElement<Space>::
evaluate_field_gradients(const vector<Real> &local_coefs,) const
-> ValueVector< Derivative<1> >
{
    Assert(this->get_values_cache(topology_id).is_filled() == true, ExcCacheNotFilled());
    Assert(this->get_values_cache(topology_id).flags_handler_.fill_gradients() == true, ExcCacheNotFilled());
    Assert(this->get_num_basis() == local_coefs.size(),
    ExcDimensionMismatch(this->get_num_basis(),local_coefs.size()));

    const auto &D1phi_hat = this->get_basis_gradients(topology_id) ;
    Assert(D1phi_hat.get_num_functions() == this->get_num_basis(),
    ExcDimensionMismatch(D1phi_hat.get_num_functions(), this->get_num_basis())) ;

    return D1phi_hat.evaluate_linear_combination(local_coefs) ;
}

template<class Space>
inline
auto
SpaceElement<Space>::
evaluate_face_field_gradients(const Index face_id, const vector<Real> &local_coefs) const
-> ValueVector< Derivative<1> >
{
    return this->evaluate_field_gradients(local_coefs,FaceTopology<dim>(face_id));
}

template<class Space>
inline
auto
SpaceElement<Space>::
evaluate_field_divergences(
    const vector<Real> &local_coefs,
    ) const -> ValueVector<Div>
{
    Assert(this->get_values_cache(topology_id).is_filled() == true, ExcCacheNotFilled());
    Assert(this->get_values_cache(topology_id).flags_handler_.fill_divergences() == true, ExcCacheNotFilled());
    Assert(this->get_num_basis() == local_coefs.size(),
    ExcDimensionMismatch(this->get_num_basis(),local_coefs.size()));

    const auto &div_phi_hat = this->get_basis_divergences(topology_id) ;
    Assert(div_phi_hat.get_num_functions() == this->get_num_basis(),
    ExcDimensionMismatch(div_phi_hat.get_num_functions(), this->get_num_basis())) ;

    return div_phi_hat.evaluate_linear_combination(local_coefs) ;
}

template<class Space>
inline
auto
SpaceElement<Space>::
evaluate_face_field_divergences(const Index face_id,
                                const vector<Real> &local_coefs) const
-> ValueVector<Div>
{
    return this->evaluate_field_divergences(local_coefs,FaceTopology<dim>(face_id));
}

template<class Space>
inline
auto
SpaceElement<Space>::
evaluate_field_hessians(const vector<Real> &local_coefs,) const -> ValueVector< Derivative<2> >
{
    Assert(this->get_values_cache(topology_id).is_filled() == true, ExcCacheNotFilled());
    Assert(this->get_values_cache(topology_id).flags_handler_.fill_hessians() == true, ExcCacheNotFilled());
    Assert(this->get_num_basis() == local_coefs.size(),
    ExcDimensionMismatch(this->get_num_basis(),local_coefs.size()));

    const auto &D2phi_hat = this->get_basis_hessians(topology_id) ;
    Assert(D2phi_hat.get_num_functions() == this->get_num_basis(),
    ExcDimensionMismatch(D2phi_hat.get_num_functions(), this->get_num_basis())) ;

    return D2phi_hat.evaluate_linear_combination(local_coefs) ;
}

template<class Space>
inline
auto
SpaceElement<Space>::
evaluate_face_field_hessians(const Index face_id,
                             const vector<Real> &local_coefs) const
-> ValueVector< Derivative<2> >
{
    return this->evaluate_field_hessians(local_coefs,FaceTopology<dim>(face_id));
}


#endif


template<class Space>
inline
Size
SpaceElement<Space>::
get_num_basis() const
{
    Index total_num_basis = 0;
    for (const auto n_basis_direction_component : n_basis_direction_)
        total_num_basis += n_basis_direction_component.flat_size();

    return total_num_basis;
}

template<class Space>
inline
auto
SpaceElement<Space>::
get_basis_offset() const -> ComponentContainer<int>
{
    return comp_offset_;
}


template<class Space>
inline
int
SpaceElement<Space>::
get_num_basis(const int i) const
{
    return n_basis_direction_[i].flat_size();
}


template<class Space>
inline
auto
SpaceElement<Space>::
get_local_to_global() const -> vector<Index>
{
    return space_->get_loc_to_global(*this);
}


template<class Space>
inline
auto
SpaceElement<Space>::
get_local_to_patch() const -> vector<Index>
{
    return space_->get_loc_to_patch(*this);
}


template<class Space>
inline
auto
SpaceElement<Space>::
get_space() const -> std::shared_ptr<const Space>
{
    return space_;
}


template<class Space>
inline
auto
SpaceElement<Space>::
get_quad_points() const -> const Quadrature<dim> &
{
    const auto &cache = this->get_values_cache(topology_id);
    Assert(cache.is_initialized(), ExcNotInitialized());

    return cache.quad_;
}

template<class Space>
inline
void
SpaceElement<Space>::
ValuesCache::
resize(const FunctionFlags &flags_handler,
       const Quadrature<dim> &quad,
       const SpaceDimensionTable &n_basis_)
{
    flags_handler_ = flags_handler;

    const auto n_points_direction = quad.get_num_points_direction();
    const auto total_n_points = n_points_direction.flat_size();
    const auto total_n_basis = n_basis_direction.total_dimension;

    Assert(total_n_points > 0, ExcLowerRange(total_n_points,1));
    Assert(total_n_basis > 0, ExcLowerRange(total_n_basis,1));

    if (flags_handler_.fill_values())
        resize_der<0>(total_n_basis,total_n_points);
    if (flags_handler_.fill_gradients())
        resize_der<1>(total_n_basis,total_n_points);
    if (flags_handler_.fill_hessians())
        resize_der<2>(total_n_basis,total_n_points);

#if 0
    if (flags_handler_.fill_divergences())
    {
        Assert(flags_handler_.fill_gradients(),
               ExcMessage("Divergence requires gradient to be filled."));

        if (div_phi_.get_num_points() != total_n_points ||
            div_phi_.get_num_functions() != total_n_basis)
        {
            div_phi_.resize(total_n_basis,total_n_points);
            div_phi_.zero();
        }
    }
    else
    {
        div_phi_.clear();
    }
#endif

    this->set_initialized(true);
}



template<class Space>
inline
auto
SpaceElement<Space>::
ValuesCache::print_info(LogStream &out) const -> void
{
    out.begin_item("Fill flags:");
    flags_handler_.print_info(out);
    out.end_item();

    out.begin_item("Values:");
    get_der<0>().print_info(out);
    out.end_item();

    out.begin_item("Gradients:");
    get_der<1>().print_info(out);
    out.end_item();

    out.begin_item("Hessians:");
    get_der<2>().print_info(out);
    out.end_item();

//    out.begin_item("Divergeces:");
//    div_phi_.print_info(out);
//    out.end_item();
}

template<class Space>
void
SpaceElement<Space>::print_info(LogStream &out) const
{
    base_t::print_info(out);
    out.begin_item("Number of element basis: ");
    n_basis_direction_.print_info(out);
    out.end_item();
}

template<class Space>
void
SpaceElement<Space>::
print_cache_info(LogStream &out) const
{
    base_t::print_cache_info(out);

    Assert(local_cache_ != nullptr,ExcNullPtr());
    local_cache_->print_info(out);
}


template<class Space>
void
SpaceElement<Space>::
LocalCache::
print_info(LogStream &out) const
{
    out.begin_item("Element Cache:");
    get_value_cache<0>(0).print_info(out);
    out.end_item();

    for (int i = 0 ; i < n_faces ; ++i)
    {
        out.begin_item("Face "+ std::to_string(i) + " Cache:");
        get_value_cache<1>(i).print_info(out);
        out.end_item();
    }

}


IGA_NAMESPACE_CLOSE

#endif
