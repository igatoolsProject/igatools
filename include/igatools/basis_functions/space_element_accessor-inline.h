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


#ifndef SPACE_ELEMENT_ACCESSOR_INLINE_H_
#define SPACE_ELEMENT_ACCESSOR_INLINE_H_

#include <igatools/basis_functions/space_element_accessor.h>


IGA_NAMESPACE_OPEN



template<class Space>
inline
SpaceElementAccessor<Space>::
SpaceElementAccessor(const std::shared_ptr<const Space> space,
                     const Index elem_index)
    :
    CartesianGridElementAccessor<dim>(space->get_grid(), elem_index),
    space_(space)
{
    Assert(space_ != nullptr, ExcNullPtr());

    using Indexer = CartesianProductIndexer<dim>;
    auto n_basis = space->get_num_basis_per_element_table();
    for (int comp_id : basis_functions_indexer_.get_active_components_id())
    {
        // creating the objects for fast conversion from flat-to-tensor indexing
        // (in practice it is an hash-table from flat to tensor indices)
        basis_functions_indexer_(comp_id) =
            std::shared_ptr<Indexer>(new Indexer(n_basis(comp_id)));
    }

    comp_offset_(0) = 0;
    for (int comp_id = 1; comp_id < Space::n_components; ++comp_id)
        comp_offset_(comp_id)= comp_offset_(comp_id-1) + n_basis.comp_dimension(comp_id);

}



template<class Space>
inline
SpaceElementAccessor<Space>::
SpaceElementAccessor(const std::shared_ptr<const Space> space,
                     const TensorIndex<dim> &elem_index)
    :
    CartesianGridElementAccessor<dim>(space->get_grid(), elem_index),
    space_(space)
{
    Assert(space_ != nullptr, ExcNullPtr());

    using Indexer = CartesianProductIndexer<dim>;
    auto n_basis = space->get_num_basis_per_element_table();
    for (int comp_id : basis_functions_indexer_.get_active_components_id())
    {
        // creating the objects for fast conversion from flat-to-tensor indexing
        // (in practice it is an hash-table from flat to tensor indices)
        basis_functions_indexer_(comp_id) =
            std::shared_ptr<Indexer>(new Indexer(n_basis(comp_id)));
    }

    comp_offset_(0) = 0;
    for (int comp_id = 1; comp_id < Space::n_components; ++comp_id)
        comp_offset_(comp_id)= comp_offset_(comp_id-1) + n_basis.comp_dimension(comp_id);

}



template<class Space>
inline
auto
SpaceElementAccessor<Space>::
as_cartesian_grid_element_accessor() -> CartesianGridElementAccessor<dim> &
{
    return static_cast<CartesianGridElementAccessor<dim> &>(*this);
}

template<class Space>
inline
auto
SpaceElementAccessor<Space>::
as_cartesian_grid_element_accessor() const -> const CartesianGridElementAccessor<dim> &
{
    return static_cast<const CartesianGridElementAccessor<dim> &>(*this);
}


template<class Space>
inline
auto
SpaceElementAccessor<Space>::
as_derived_element_accessor() -> DerivedElementAccessor &
{
    return static_cast<DerivedElementAccessor &>(*this);
}

template<class Space>
inline
auto
SpaceElementAccessor<Space>::
as_derived_element_accessor() const -> const DerivedElementAccessor &
{
    return static_cast<const DerivedElementAccessor &>(*this);
}

template<class Space>
inline
auto
SpaceElementAccessor<Space>::
evaluate_basis_values_at_points(const vector<Point> &points) const -> ValueTable<Value>
{
    return this->as_derived_element_accessor().template evaluate_basis_derivatives_at_points<0>(points);
}

template<class Space>
inline
auto
SpaceElementAccessor<Space>::
evaluate_basis_gradients_at_points(const vector<Point> &points) const -> ValueTable<Derivative<1> >
{
    return this->as_derived_element_accessor().template evaluate_basis_derivatives_at_points<1>(points);
}

template<class Space>
inline
auto
SpaceElementAccessor<Space>::
evaluate_basis_hessians_at_points(const vector<Point> &points) const -> ValueTable<Derivative<2> >
{
    return this->as_derived_element_accessor().template evaluate_basis_derivatives_at_points<2>(points);
}

template<class Space>
inline
auto
SpaceElementAccessor<Space>::
get_basis_values(const Quadrature<dim> &quad) const -> ValueTable< Value >
{
    this->evaluate_basis_values_at_points(quad.get_points());
}

template<class Space>
inline
auto
SpaceElementAccessor<Space>::
get_basis_gradients(const Quadrature<dim> &quad) const -> ValueTable< Derivative<1> >
{
    this->evaluate_basis_gradients_at_points(quad.get_points());
}

template<class Space>
inline
auto
SpaceElementAccessor<Space>::
get_basis_hessians(const Quadrature<dim> &quad) const -> ValueTable< Derivative<2> >
{
    this->evaluate_basis_hessians_at_points(quad.get_points());
}




template<class Space>
template <int deriv_order>
inline
auto
SpaceElementAccessor<Space>::
evaluate_field_derivatives_at_points(
    const vector<Real> &local_coefs,
    const vector<Point> &points) const ->
ValueVector< Conditional< deriv_order==0,Value,Derivative<deriv_order> > >
{
    const auto &derived_element_accessor = this->as_derived_element_accessor();
    Assert(derived_element_accessor.get_num_basis() == local_coefs.size(),
    ExcDimensionMismatch(derived_element_accessor.get_num_basis(),local_coefs.size()));

    const auto derivatives_phi_hat =
    derived_element_accessor.template evaluate_basis_derivatives_at_points<deriv_order>(points);
    Assert(derivatives_phi_hat.get_num_functions() == derived_element_accessor.get_num_basis(),
    ExcDimensionMismatch(derivatives_phi_hat.get_num_functions(), derived_element_accessor.get_num_basis())) ;

    return derivatives_phi_hat.evaluate_linear_combination(local_coefs) ;
}


template<class Space>
inline
auto
SpaceElementAccessor<Space>::
evaluate_field_values_at_points(
    const vector<Real> &local_coefs,
    const vector<Point> &points) const -> ValueVector<Value>
{
    return this->evaluate_field_derivatives_at_points<0>(local_coefs,points);
}

template<class Space>
inline
auto
SpaceElementAccessor<Space>::
evaluate_field_gradients_at_points(
    const vector<Real> &local_coefs,
    const vector<Point> &points) const -> ValueVector<Derivative<1> >
{
    return this->evaluate_field_derivatives_at_points<1>(local_coefs,points);
}

template<class Space>
inline
auto
SpaceElementAccessor<Space>::
evaluate_field_hessians_at_points(
    const vector<Real> &local_coefs,
    const vector<Point> &points) const -> ValueVector<Derivative<2> >
{
    return this->evaluate_field_derivatives_at_points<2>(local_coefs,points);
}


template<class Space>
inline
auto
SpaceElementAccessor<Space>::
get_basis_values(const TopologyId<dim> &topology_id) const -> ValueTable<Value> const &
{
    const auto &cache = this->get_values_cache(topology_id);
    Assert(cache.is_filled() == true, ExcCacheNotFilled());
    Assert(cache.flags_handler_.values_filled(), ExcCacheNotFilled());

    return cache.phi_;
}

template<class Space>
inline
auto
SpaceElementAccessor<Space>::
get_face_basis_values(const Index face_id) const -> ValueTable<Value> const &
{
    return this->get_basis_values(FaceTopology<dim>(face_id));
}



template<class Space>
inline
auto
SpaceElementAccessor<Space>::
get_basis_values(const Index i,const TopologyId<dim> &topology_id) const -> typename ValueTable<Value>::const_view
{
    return this->get_basis_values(topology_id).get_function_view(i);
}

template<class Space>
inline
auto
SpaceElementAccessor<Space>::
get_basis_value(const Index basis, const Index qp,const TopologyId<dim> &topology_id) const -> Value const &
{
    Assert(qp >= 0 && qp < this->get_values_cache(topology_id).quad_.get_num_points_direction().flat_size(),
           ExcIndexRange(qp,0,this->get_values_cache(topology_id).quad_.get_num_points_direction().flat_size()));
    return this->get_basis_values(basis,topology_id)[qp];
}




template<class Space>
inline
auto
SpaceElementAccessor<Space>::
get_basis_gradients(const TopologyId<dim> &topology_id) const -> ValueTable<Derivative<1>> const &
{
    const auto &cache = this->get_values_cache(topology_id);
    Assert(cache.is_filled() == true, ExcCacheNotFilled());
    Assert(cache.flags_handler_.gradients_filled(), ExcCacheNotFilled());

    return cache.D1phi_;
}

template<class Space>
inline
auto
SpaceElementAccessor<Space>::
get_face_basis_gradients(const Index face_id) const -> ValueTable<Derivative<1>> const &
{
    return this->get_basis_gradients(FaceTopology<dim>(face_id));
}



template<class Space>
inline
auto
SpaceElementAccessor<Space>::
get_basis_gradients(const Index i,const TopologyId<dim> &topology_id) const -> typename ValueTable<Derivative<1>>::const_view
{
    return this->get_basis_gradients(topology_id).get_function_view(i);
}

template<class Space>
inline
auto
SpaceElementAccessor<Space>::
get_basis_gradient(const Index basis, const Index qp,const TopologyId<dim> &topology_id) const -> Derivative<1> const &
{
    Assert(qp >= 0 && qp < this->get_values_cache(topology_id).quad_.get_num_points_direction().flat_size(),
           ExcIndexRange(qp,0,this->get_values_cache(topology_id).quad_.get_num_points_direction().flat_size()));
    return this->get_basis_gradients(basis,topology_id)[qp];
}





template<class Space>
inline
auto
SpaceElementAccessor<Space>::
get_basis_hessians(const TopologyId<dim> &topology_id) const -> ValueTable<Derivative<2>> const &
{
    const auto &cache = this->get_values_cache(topology_id);
    Assert(cache.is_filled() == true, ExcCacheNotFilled());
    Assert(cache.flags_handler_.hessians_filled(), ExcCacheNotFilled());

    return cache.D2phi_;
}

template<class Space>
inline
auto
SpaceElementAccessor<Space>::
get_face_basis_hessians(const Index face_id) const -> ValueTable<Derivative<2>> const &
{
    return this->get_basis_hessians(FaceTopology<dim>(face_id));
}



template<class Space>
inline
auto
SpaceElementAccessor<Space>::
get_basis_hessians(const Index i,const TopologyId<dim> &topology_id) const -> typename ValueTable<Derivative<2>>::const_view
{
    return this->get_basis_hessians(topology_id).get_function_view(i);
}

template<class Space>
inline
auto
SpaceElementAccessor<Space>::
get_basis_hessian(const Index basis, const Index qp,const TopologyId<dim> &topology_id) const -> Derivative<2> const &
{
    Assert(qp >= 0 && qp < this->get_values_cache(topology_id).quad_.get_num_points_direction().flat_size(),
           ExcIndexRange(qp,0,this->get_values_cache(topology_id).quad_.get_num_points_direction().flat_size()));
    return this->get_basis_hessians(basis,topology_id)[qp];
}




template<class Space>
inline
auto
SpaceElementAccessor<Space>::
get_basis_divergences(const TopologyId<dim> &topology_id) const -> ValueTable<Div> const &
{
    const auto &cache = this->get_values_cache(topology_id);
    Assert(cache.is_filled() == true, ExcCacheNotFilled());
    Assert(cache.flags_handler_.divergences_filled(), ExcCacheNotFilled());

    return cache.div_phi_;
}

template<class Space>
inline
auto
SpaceElementAccessor<Space>::
get_face_basis_divergences(const Index face_id) const -> ValueTable<Div> const &
{
    return this->get_basis_divergences(FaceTopology<dim>(face_id));
}



template<class Space>
inline
auto
SpaceElementAccessor<Space>::
get_basis_divergences(const Index i,const TopologyId<dim> &topology_id) const -> typename ValueTable<Div>::const_view
{
    return this->get_basis_divergences(topology_id).get_function_view(i);
}

template<class Space>
inline
auto
SpaceElementAccessor<Space>::
get_basis_divergence(const Index basis, const Index qp,const TopologyId<dim> &topology_id) const -> Div const &
{
    Assert(qp >= 0 && qp < this->get_values_cache(topology_id).quad_.get_num_points_direction().flat_size(),
           ExcIndexRange(qp,0,this->get_values_cache(topology_id).quad_.get_num_points_direction().flat_size()));
    return this->get_basis_divergences(basis,topology_id)[qp];
}


template<class Space>
inline
void
SpaceElementAccessor<Space>::
fill_face_cache(const Index face_id)
{
    this->as_derived_element_accessor().fill_cache(FaceTopology<dim>(face_id));
}


template<class Space>
inline
auto
SpaceElementAccessor<Space>::
get_values_cache(const TopologyId<dim> &topology_id) const -> const ValuesCache &
{
    Assert(topology_id.is_element() || topology_id.is_face(),
           ExcMessage("Only element or face topology is allowed."));
    if (topology_id.is_element())
    {
        return elem_values_;
    }
    else
    {
        Assert(topology_id.get_id()>=0 && topology_id.get_id() < n_faces,
               ExcIndexRange(topology_id.get_id(),0,n_faces));

        Assert(this->is_boundary(topology_id.get_id()),
               ExcMessage("The requested face_id=" +
                          std::to_string(topology_id.get_id()) +
                          " is not a boundary for the element"));
        return face_values_[topology_id.get_id()];
    }
}


template<class Space>
inline
auto
SpaceElementAccessor<Space>::
get_values_cache(const TopologyId<dim> &topology_id) -> ValuesCache &
{
    Assert(topology_id.is_element() || topology_id.is_face(),
    ExcMessage("Only element or face topology is allowed."));
    if (topology_id.is_element())
    {
        return elem_values_;
    }
    else
    {
        Assert(topology_id.get_id()>=0 && topology_id.get_id() < n_faces,
        ExcIndexRange(topology_id.get_id(),0,n_faces));

        Assert(this->is_boundary(topology_id.get_id()),
        ExcMessage("The requested face_id=" +
        std::to_string(topology_id.get_id()) +
        " is not a boundary for the element"));
        return face_values_[topology_id.get_id()];
    }
}


template<class Space>
inline
void
SpaceElementAccessor<Space>::
ValuesCache::
reset(const BasisElemValueFlagsHandler &flags_handler,
      const ComponentContainer<TensorSize<dim> > &n_basis_direction,
      const Quadrature<dim> &quad)
{
    quad_ = quad;
    const TensorSize<dim> n_points_direction = quad_.get_num_points_direction();

    flags_handler_ = flags_handler;

    //--------------------------------------------------------------------------
    // computing the total number of basis functions and the number of evaluation points
    const int total_n_points = n_points_direction.flat_size();

    int total_n_basis = 0;
    for (int i = 0; i < Space::n_components; ++i)
        total_n_basis += n_basis_direction(i).flat_size();

    Assert(total_n_points > 0, ExcLowerRange(total_n_points,1));
    Assert(total_n_basis > 0, ExcLowerRange(total_n_basis,1));
    //--------------------------------------------------------------------------


    //--------------------------------------------------------------------------
    // resizing the containers for the basis functions
    if (flags_handler_.fill_values())
    {
        if (phi_.get_num_points() != total_n_points ||
            phi_.get_num_functions() != total_n_basis)
        {
            phi_.resize(total_n_basis,total_n_points);
            phi_.zero();
        }
    }
    else
    {
        phi_.clear();
    }

    if (flags_handler_.fill_gradients())
    {
        if (D1phi_.get_num_points() != total_n_points ||
            D1phi_.get_num_functions() != total_n_basis)
        {
            D1phi_.resize(total_n_basis,total_n_points);
            D1phi_.zero();
        }
    }
    else
    {
        D1phi_.clear();
    }


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



    if (flags_handler_.fill_hessians())
    {
        if (D2phi_.get_num_points() != total_n_points ||
            D2phi_.get_num_functions() != total_n_basis)
        {
            D2phi_.resize(total_n_basis,total_n_points);
            D2phi_.zero();
        }
    }
    else
    {
        D2phi_.clear();
    }
    //--------------------------------------------------------------------------


    this->set_initialized(true);
}


template<class Space>
inline
void
SpaceElementAccessor<Space>::
ElementValuesCache::
reset(const BasisElemValueFlagsHandler &flags_handler,
      const ComponentContainer<TensorSize<dim> > &n_basis_direction,
      const Quadrature<dim> &quad)
{
    ValuesCache::reset(flags_handler, n_basis_direction,quad);
}

template<class Space>
inline
void
SpaceElementAccessor<Space>::
FaceValuesCache::
reset(const Index face_id,
      const BasisFaceValueFlagsHandler &flags_handler,
      const ComponentContainer<TensorSize<dim> > &n_basis_direction,
      const Quadrature<dim> &quad1)
{
    Assert(face_id < n_faces && face_id >= 0, ExcIndexRange(face_id,0,n_faces));

    const auto quad = quad1.collapse_to_face(face_id);

    ValuesCache::reset(flags_handler, n_basis_direction,quad);
}



template<class Space>
inline
void
SpaceElementAccessor<Space>::
FaceValuesCache::
reset(const Index face_id,
      const BasisFaceValueFlagsHandler &flags_handler,
      const ComponentContainer<TensorSize<dim> > &n_basis_direction,
      const Quadrature<dim-1> &quad)
{
    Assert(false,ExcNotImplemented()) ;
    AssertThrow(false,ExcNotImplemented()) ;
}



template<class Space>
inline
void
SpaceElementAccessor<Space>::
reset_element_and_faces_cache(const ValueFlags fill_flag,
                              const Quadrature<dim> &quad)
{
    //--------------------------------------------------------------------------
    BasisElemValueFlagsHandler elem_flags_handler(fill_flag);
    BasisFaceValueFlagsHandler face_flags_handler(fill_flag);


    Assert(!elem_flags_handler.fill_none() ||
           !face_flags_handler.fill_none(),
           ExcMessage("Nothing to reset"));

    auto n_basis = space_->get_num_basis_per_element_table();
    if (!elem_flags_handler.fill_none())
        this->elem_values_.reset(elem_flags_handler, n_basis, quad);


    if (!face_flags_handler.fill_none())
    {
        Index face_id = 0 ;
        for (auto &face_value : this->face_values_)
            face_value.reset(face_id++, face_flags_handler, n_basis, quad);
    }
    //--------------------------------------------------------------------------
}



template<class Space>
inline
auto
SpaceElementAccessor<Space>::
ValuesCache::
get_values() const -> const ValueTable<Value> &
{
    return phi_;
}

template<class Space>
inline
auto
SpaceElementAccessor<Space>::
ValuesCache::
get_gradients() const -> const ValueTable<Derivative<1>> &
{
    return D1phi_;
}


template<class Space>
inline
auto
SpaceElementAccessor<Space>::
ValuesCache::
get_hessians() const -> const ValueTable<Derivative<2>> &
{
    return D2phi_;
}

template<class Space>
inline
auto
SpaceElementAccessor<Space>::
ValuesCache::
get_divergences() const -> const ValueTable<Div> &
{
    return div_phi_;
}


template<class Space>
inline
auto
SpaceElementAccessor<Space>::
evaluate_field(const vector<Real> &local_coefs,const TopologyId<dim> &topology_id) const
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
SpaceElementAccessor<Space>::
evaluate_face_field(const Index face_id, const vector<Real> &local_coefs) const
-> ValueVector<Value>
{
    return this->evaluate_field(local_coefs,FaceTopology<dim>(face_id));
}


template<class Space>
inline
auto
SpaceElementAccessor<Space>::
evaluate_field_gradients(const vector<Real> &local_coefs,const TopologyId<dim> &topology_id) const
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
SpaceElementAccessor<Space>::
evaluate_face_field_gradients(const Index face_id, const vector<Real> &local_coefs) const
-> ValueVector< Derivative<1> >
{
    return this->evaluate_field_gradients(local_coefs,FaceTopology<dim>(face_id));
}

template<class Space>
inline
auto
SpaceElementAccessor<Space>::
evaluate_field_divergences(
    const vector<Real> &local_coefs,
    const TopologyId<dim> &topology_id) const -> ValueVector<Div>
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
SpaceElementAccessor<Space>::
evaluate_face_field_divergences(const Index face_id, const vector<Real> &local_coefs) const
-> ValueVector<Div>
{
    return this->evaluate_field_divergences(local_coefs,FaceTopology<dim>(face_id));
}

template<class Space>
inline
auto
SpaceElementAccessor<Space>::
evaluate_field_hessians(const vector<Real> &local_coefs,const TopologyId<dim> &topology_id) const -> ValueVector< Derivative<2> >
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
SpaceElementAccessor<Space>::
evaluate_face_field_hessians(const Index face_id, const vector<Real> &local_coefs) const
-> ValueVector< Derivative<2> >
{
    return this->evaluate_field_hessians(local_coefs,FaceTopology<dim>(face_id));
}





template<class Space>
inline
Size
SpaceElementAccessor<Space>::
get_num_basis() const
{
    return this->space_->get_num_basis_per_element();
}


template<class Space>
inline
int
SpaceElementAccessor<Space>::
get_num_basis(const int i) const
{
    return space_->get_num_basis_per_element(i);
}


template<class Space>
inline
auto
SpaceElementAccessor<Space>::
get_local_to_global() const -> vector<Index>
{
    return space_->get_loc_to_global(this->get_tensor_index());
}


template<class Space>
inline
auto
SpaceElementAccessor<Space>::
get_space() const -> std::shared_ptr<const Space>
{
    return space_;
}


template<class Space>
inline
auto
SpaceElementAccessor<Space>::
get_quad_points(const TopologyId<dim> &topology_id) const -> const Quadrature<dim> &
{
    const auto &cache = this->get_values_cache(topology_id);
    Assert(cache.is_initialized(), ExcNotInitialized());

    return cache.quad_;
}





IGA_NAMESPACE_CLOSE

#endif // #ifndef SPACE_ELEMENT_ACCESSOR_INLINE_H_
