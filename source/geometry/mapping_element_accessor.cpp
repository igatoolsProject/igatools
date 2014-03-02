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


#include <igatools/geometry/mapping_element_accessor.h>

#include <igatools/base/exceptions.h>
#include <igatools/geometry/unit_element.h>


using std::array ;
using std::vector ;
using std::shared_ptr ;

IGA_NAMESPACE_OPEN

template<int dim_ref_, int codim_ >
MappingElementAccessor<dim_ref_,codim_>::
MappingElementAccessor(Mapping<dim,codim> &mapping, const int index)
    :
    CartesianGridElementAccessor<dim>(*mapping.get_grid(), index),
    mapping_(&mapping)
{
    Assert(mapping_->get_grid() != nullptr, ExcNullPtr());
}



template< int dim_ref_, int codim_ >
void
MappingElementAccessor<dim_ref_,codim_>::
ElementValuesCache::
reset(const ValueFlags fill_flag,
      const Quadrature<dim> &quad)
{
    quad_ = quad;
    num_points_ = quad.get_num_points();

    Assert(contains(fill_flag, ValueFlags::none),
           ExcMessage("nothing to reset"));


    if (contains(fill_flag, ValueFlags::point) ||
        contains(fill_flag, ValueFlags::map_value))
    {
        fill_values_ = true ;
        if (values_.size() != num_points_)
            values_.resize(num_points_);

        values_.zero();
    }
    else
    {
        fill_values_ = false ;
        values_.clear();
    }

    if (contains(fill_flag, ValueFlags::map_gradient))
    {
        fill_gradients_ = true ;
        if (gradients_.size() != num_points_)
            gradients_.resize(num_points_);

        gradients_.zero();
    }
    else
    {
        fill_gradients_ = false ;
        gradients_.clear();
    }

    if (contains(fill_flag, ValueFlags::map_hessian))
    {
        fill_hessians_ = true ;
        if (hessians_.size() != num_points_)
            hessians_.resize(num_points_);

        hessians_.zero();
    }
    else
    {
        fill_hessians_ = false ;
        hessians_.clear();
    }

    if (contains(fill_flag, ValueFlags::map_inv_gradient))
    {
        fill_inv_gradients_ = true ;
        if (inv_gradients_.size() != num_points_)
            inv_gradients_.resize(num_points_);

        inv_gradients_.zero();
    }
    else
    {
        fill_inv_gradients_ = false ;
        inv_gradients_.clear();
    }

    if (contains(fill_flag, ValueFlags::map_inv_hessian))
    {
        fill_inv_hessians_ = true ;
        if (inv_hessians_.size() != num_points_)
            inv_hessians_.resize(num_points_);

        inv_hessians_.zero();
    }
    else
    {
        fill_inv_hessians_ = false ;
        inv_hessians_.clear();
    }

    if (contains(fill_flag, ValueFlags::measure))
    {
        Assert(fill_gradients_, ExcNotInitialized());
        fill_dets_ = true ;
        if (dets_.size() != num_points_)
            dets_.resize(num_points_);

        dets_.zero();
    }
    else
    {
        fill_dets_ = false ;
        dets_.clear();
    }


    if (contains(fill_flag, ValueFlags::w_measure))
    {
        Assert(fill_dets_, ExcNotInitialized());
        fill_w_measures_ = true ;
        if (w_measures_.size() != num_points_)
            w_measures_.resize(num_points_);

        w_measures_.zero();
    }
    else
    {
        fill_w_measures_ = false ;
        w_measures_.clear();
    }

    this->set_initialized(true);
}



template< int dim_ref_, int codim_ >
void
MappingElementAccessor<dim_ref_,codim_>::
FaceValuesCache::
reset(const Index face_id,
      const ValueFlags fill_flag,
      const Quadrature<dim> &quad)
{
    Assert(face_id < n_faces && face_id >= 0, ExcIndexRange(face_id,0,n_faces));

    quad_ = quad.collapse_to_face(face_id);
    num_points_ = quad_.get_num_points();

    Assert(contains(fill_flag, ValueFlags::none),
           ExcMessage("nothing to reset"));


    if (contains(fill_flag, ValueFlags::face_point) ||
        contains(fill_flag, ValueFlags::map_face_value))
    {
        fill_values_ = true ;
        if (values_.size() != num_points_)
            values_.resize(num_points_);

        values_.zero();
    }
    else
    {
        fill_values_ = false ;
        values_.clear();
    }

    if (contains(fill_flag, ValueFlags::map_face_gradient))
    {
        fill_gradients_ = true ;
        if (gradients_.size() != num_points_)
            gradients_.resize(num_points_);

        gradients_.zero();
    }
    else
    {
        fill_gradients_ = false ;
        gradients_.clear();
    }

    if (contains(fill_flag, ValueFlags::map_face_hessian))
    {
        fill_hessians_ = true ;
        if (hessians_.size() != num_points_)
            hessians_.resize(num_points_);

        hessians_.zero();
    }
    else
    {
        fill_hessians_ = false ;
        hessians_.clear();
    }

    if (contains(fill_flag, ValueFlags::map_face_inv_gradient))
    {
        fill_inv_gradients_ = true ;
        if (inv_gradients_.size() != num_points_)
            inv_gradients_.resize(num_points_);

        inv_gradients_.zero();
    }
    else
    {
        fill_inv_gradients_ = false ;
        inv_gradients_.clear();
    }

    if (contains(fill_flag, ValueFlags::map_face_inv_hessian))
    {
        fill_inv_hessians_ = true ;
        if (inv_hessians_.size() != num_points_)
            inv_hessians_.resize(num_points_);

        inv_hessians_.zero();
    }
    else
    {
        fill_inv_hessians_ = false ;
        inv_hessians_.clear();
    }

    if (contains(fill_flag, ValueFlags::face_measure))
    {
        Assert(fill_gradients_, ExcNotInitialized());
        fill_dets_ = true ;
        if (dets_.size() != num_points_)
            dets_.resize(num_points_);

        dets_.zero();
    }
    else
    {
        fill_dets_ = false ;
        dets_.clear();
    }


    if (contains(fill_flag, ValueFlags::face_w_measure))
    {
        Assert(fill_dets_, ExcNotInitialized());
        fill_w_measures_ = true ;
        if (w_measures_.size() != num_points_)
            w_measures_.resize(num_points_);

        w_measures_.zero();
    }
    else
    {
        fill_w_measures_ = false ;
        w_measures_.clear();
    }

    if (contains(fill_flag, ValueFlags::face_normal))
    {
        fill_normals_ = true ;
        if (normals_.size() != num_points_)
            normals_.resize(num_points_);

        normals_.zero();
    }
    else
    {
        fill_normals_ = false ;
        normals_.clear();
    }

    this->set_initialized(true);
}



template< int dim_ref_, int codim_ >
void
MappingElementAccessor<dim_ref_,codim_>::
FaceValuesCache::
reset(const Index face_id,
      const ValueFlags fill_flag,
      const Quadrature<dim-1> &quad1)
{
    Assert(false, ExcNotImplemented());
}



template< int dim_ref_, int codim_ >
void
MappingElementAccessor<dim_ref_,codim_>::
init_values(const ValueFlags fill_flag,
            const Quadrature<dim> &quad)
{
//  Assert(elem_values_.initialized_ == false, ExcMessage("Already initialized")) ;
    Assert((fill_flag|admisible_flag) == admisible_flag,
           typename CartesianGridElementAccessor<dim_ref_>::ExcFillFlagNotSupported(admisible_flag, fill_flag));

    // initalizing the cache of the CartesianGridElementAccessor
    {
        ValueFlags grid_flag = mapping_->required_flags();;
        if (contains(fill_flag , ValueFlags::point))
            grid_flag |= ValueFlags::point;
        if (contains(fill_flag , ValueFlags::w_measure))
            grid_flag |= ValueFlags::w_measure;
        if (contains(fill_flag , ValueFlags::face_point))
            grid_flag |= ValueFlags::face_point;
        if (contains(fill_flag , ValueFlags::face_w_measure))
            grid_flag |= ValueFlags::face_w_measure;
        CartesianGridElementAccessor<dim_ref_>::init_values(grid_flag,quad);
    }

    auto f_flag = fill_flag;
    if (contains(f_flag , ValueFlags::w_measure))
        f_flag |= ValueFlags::measure;
    if (contains(f_flag , ValueFlags::measure))
        f_flag |= ValueFlags::map_gradient;
    if (contains(f_flag , ValueFlags::face_w_measure))
        f_flag |= ValueFlags::face_measure;
    if (contains(f_flag , ValueFlags::face_measure))
        f_flag |= ValueFlags::map_face_gradient;

    elem_values_.reset(f_flag, quad);

    Index face_id = 0 ;
    for (auto& face_value : face_values_)
        face_value.reset(face_id++, f_flag, quad);

    mapping_->init_element(fill_flag, quad);
}


template< int dim_ref_, int codim_ >
void
MappingElementAccessor<dim_ref_,codim_>::
init_face_values(const Index face_id,
                 const ValueFlags fill_flag,
                 const Quadrature<dim-1> &quad)
{
    Assert(false, ExcNotImplemented());
}



template< int dim_ref_, int codim_ >
void
MappingElementAccessor<dim_ref_,codim_>::
fill_values()
{
    CartesianGridElementAccessor<dim_ref_>::fill_values();
    mapping_->set_element(*this);

    Assert(elem_values_.is_initialized(), ExcNotInitialized()) ;

    Assert(elem_values_.num_points_ ==
           elem_values_.quad_.get_points().flat_size(),
           ExcDimensionMismatch(elem_values_.num_points_,
                                elem_values_.quad_.get_points().flat_size()));

    if (elem_values_.fill_values_)
        mapping_->evaluate(elem_values_.values_);

    if (elem_values_.fill_gradients_)
        mapping_->evaluate_gradients(elem_values_.gradients_);

    if (elem_values_.fill_hessians_)
        mapping_->evaluate_hessians(elem_values_.hessians_);

    if (elem_values_.fill_inv_gradients_)
        for (Index i = 0; i < elem_values_.num_points_; i++)
        {
            elem_values_.dets_[i] =
                inverse<dim,space_dim>(elem_values_.gradients_[i],
                                       elem_values_.inv_gradients_[i]);
        }
    /*
     * To fill the hessian of F{^-1}, we use the formula
     * D2F{^-1} [u] = DF{^-1} * D2F[u] * DF{^-1},
     * This formula can be obtained by differentiating the identity
     * DF * DF{^-1} = I
     */
    if (elem_values_.fill_inv_hessians_)
        for (Index i = 0; i < elem_values_.num_points_; i++)
        {
            const auto &DF_inv = elem_values_.inv_gradients_[i];
            const auto &D2F = elem_values_.hessians_[i];
            for (int u=0; u<dim; ++u) //TODO: should we define a compose in tensor for this?
            {
                const auto temp = compose(DF_inv, D2F[u]);
                elem_values_.inv_hessians_[i][u] = compose(temp, DF_inv);
            }
        }


    if (elem_values_.fill_dets_ ||
        elem_values_.fill_w_measures_)
    {
        for (Index i = 0; i < elem_values_.num_points_; i++)
            elem_values_.dets_[i] =
                determinant<dim,space_dim>(elem_values_.gradients_[i]);


        if (elem_values_.fill_w_measures_)
        {
            const ValueVector<Real> &dets_map = elem_values_.dets_ ;
            auto weights = CartesianGridElementAccessor<dim_ref_>::get_w_measures();

            for (Index i = 0; i < elem_values_.num_points_; i++)
                elem_values_.w_measures_[i] = dets_map[i] * weights[i] ;
        }
    }

    elem_values_.set_filled(true);
}



template< int dim_ref_, int codim_ >
void
MappingElementAccessor<dim_ref_,codim_>::
fill_face_values(const Index face_id)
{
    Assert(face_id < n_faces && face_id >= 0, ExcIndexRange(face_id,0,n_faces));

    CartesianGridElementAccessor<dim_ref_>::fill_face_values(face_id);
    mapping_->set_face_element(face_id, *this);

    auto &face_value = face_values_[face_id] ;

    const auto &num_points = face_value.num_points_ ;

    Assert(face_value.is_initialized(), ExcNotInitialized()) ;

    Assert(num_points ==
           face_value.quad_.get_points().flat_size(),
           ExcDimensionMismatch(num_points,
                                face_value.quad_.get_points().flat_size()));

    if (face_value.fill_values_)
        mapping_->evaluate_face(face_id, face_value.values_);

    if (face_value.fill_gradients_)
        mapping_->evaluate_face_gradients(face_id, face_value.gradients_);

    if (face_value.fill_hessians_)
        mapping_->evaluate_face_hessians(face_id, face_value.hessians_);

//    TODO: to be solved
//    if (face_value.fill_inv_gradients_)
//        for (Index i = 0; i < num_points; i++)
//        {
//            face_value.dets_[i] =
//                inverse<dim,space_dim>(face_value.gradients_[i],
//                                       face_value.inv_gradients_[i]);
//        }
    /*
     * To fill the hessian of F{^-1}, we use the formula
     * D2F{^-1} [u] = DF{^-1} * D2F[u] * DF{^-1},
     * This formula can be obtained by differentiating the identity
     * DF * DF{^-1} = I
     */
    if (face_value.fill_inv_hessians_)
        for (Index i = 0; i < num_points; i++)
        {
            const auto &DF_inv = face_value.inv_gradients_[i];
            const auto &D2F = face_value.hessians_[i];
            for (int u=0; u<dim; ++u) //TODO: should we define a compose in tensor for this?
            {
                const auto temp = compose(DF_inv, D2F[u]);
                face_value.inv_hessians_[i][u] = compose(temp, DF_inv);
            }
        }


    if (face_value.fill_dets_ ||
        face_value.fill_w_measures_)
    {
//    TODO: to be solved
//        for (Index i = 0; i < num_points; i++)
//            face_value.dets_[i] =
//                determinant<dim,space_dim>(face_value.gradients_[i]);


        if (face_value.fill_w_measures_)
        {
            const ValueVector<Real> &dets_map = face_value.dets_ ;
            auto weights = CartesianGridElementAccessor<dim_ref_>::get_face_w_measures(face_id);

            for (Index i = 0; i < num_points; i++)
                face_value.w_measures_[i] = dets_map[i] * weights[i] ;
        }
    }

    if (face_value.fill_normals_)
    {
        // Obtain n_hat from UnitElement
        Point<dim_ref_> n_hat = UnitElement<dim_ref_>::face_normal[face_id] ;
        for (Index i = 0; i < num_points; i++)
        {
            const auto DF_inv_t = co_tensor(transpose(face_value.inv_gradients_[i])) ;
            face_value.normals_[i] = action(DF_inv_t, n_hat);
            face_value.normals_[i] /= face_value.normals_[i].norm();
        }
    }

    face_value.set_filled(true);
}



template< int dim_ref_, int codim_ >
auto
MappingElementAccessor<dim_ref_,codim_>::
get_values_map() const -> const ValueVector<ValueMap> &
{
    Assert(elem_values_.is_filled(), ExcMessage("The cache is not filled."));
    Assert(elem_values_.fill_values_, ExcNotInitialized());
    return elem_values_.values_;
}



template< int dim_ref_, int codim_ >
auto
MappingElementAccessor<dim_ref_,codim_>::
get_gradients_map() const -> const ValueVector<GradientMap> &
{
    Assert(elem_values_.is_filled(), ExcMessage("The cache is not filled."));
    Assert(elem_values_.fill_gradients_, ExcNotInitialized());
    return elem_values_.gradients_;
}



template< int dim_ref_, int codim_ >
auto
MappingElementAccessor<dim_ref_,codim_>::
get_hessians_map() const -> const ValueVector<HessianMap> &
{
    Assert(elem_values_.is_filled(), ExcMessage("The cache is not filled."));
    Assert(elem_values_.fill_hessians_, ExcNotInitialized());
    return elem_values_.hessians_;
}



template< int dim_ref_, int codim_ >
auto
MappingElementAccessor<dim_ref_,codim_>::
get_inv_gradients_map() const -> const ValueVector<Derivatives<space_dim,dim,1,1>> &
{
    Assert(elem_values_.is_filled(), ExcMessage("The cache is not filled."));
    Assert(elem_values_.fill_inv_gradients_, ExcNotInitialized());
    return elem_values_.inv_gradients_;
}



template< int dim_ref_, int codim_ >
auto
MappingElementAccessor<dim_ref_,codim_>::
get_inv_hessians_map() const -> const ValueVector<Derivatives<space_dim,dim,1,2>> &
{
    Assert(elem_values_.is_filled(), ExcMessage("The cache is not filled."));
    Assert(elem_values_.fill_inv_hessians_, ExcNotInitialized());
    return elem_values_.inv_hessians_;
}



template< int dim_ref_, int codim_ >
auto
MappingElementAccessor<dim_ref_,codim_>::
get_dets_map() const -> const ValueVector<Real> &
{
    Assert(elem_values_.is_filled(), ExcMessage("The cache is not filled."));
    Assert(elem_values_.fill_dets_, ExcNotInitialized());
    return elem_values_.dets_;
}



template< int dim_ref_, int codim_ >
auto
MappingElementAccessor<dim_ref_,codim_>::
get_w_measures() const -> const ValueVector<Real> &
{
    Assert(elem_values_.is_filled(), ExcMessage("The cache is not filled."));
    Assert(elem_values_.fill_w_measures_, ExcNotInitialized());
    return elem_values_.w_measures_;
}



template< int dim_ref_, int codim_ >
auto
MappingElementAccessor<dim_ref_,codim_>::
get_face_values_map(const Index face_id) const -> const ValueVector<ValueMap> &
{
    Assert(face_id < n_faces && face_id >= 0, ExcIndexRange(face_id,0,n_faces));
    Assert(face_values_[face_id].is_filled(), ExcMessage("The cache is not filled."));
    Assert(face_values_[face_id].fill_values_, ExcNotInitialized());
    return face_values_[face_id].values_;
}



template< int dim_ref_, int codim_ >
auto
MappingElementAccessor<dim_ref_,codim_>::
get_face_gradients_map(const Index face_id) const ->
const ValueVector<GradientFaceMap> &
{
    Assert(face_id < n_faces && face_id >= 0, ExcIndexRange(face_id,0,n_faces));
    Assert(face_values_[face_id].is_filled(), ExcMessage("The cache is not filled."));
    Assert(face_values_[face_id].fill_gradients_, ExcNotInitialized());
    return face_values_[face_id].gradients_;
}



template< int dim_ref_, int codim_ >
auto
MappingElementAccessor<dim_ref_,codim_>::
get_face_hessians_map(const Index face_id) const ->
const ValueVector<HessianFaceMap> &
{
    Assert(face_id < n_faces && face_id >= 0, ExcIndexRange(face_id,0,n_faces));
    Assert(face_values_[face_id].is_filled(), ExcMessage("The cache is not filled."));
    Assert(face_values_[face_id].fill_hessians_, ExcNotInitialized());
    return face_values_[face_id].hessians_;
}



template< int dim_ref_, int codim_ >
auto
MappingElementAccessor<dim_ref_,codim_>::
get_face_inv_gradients_map(const Index face_id) const ->
const ValueVector<Derivatives<space_dim,face_dim,1,1>> &
{
    Assert(face_id < n_faces && face_id >= 0, ExcIndexRange(face_id,0,n_faces));
    Assert(face_values_[face_id].is_filled(), ExcMessage("The cache is not filled."));
    Assert(face_values_[face_id].fill_inv_gradients_, ExcNotInitialized());
    return face_values_[face_id].inv_gradients_;
}



template< int dim_ref_, int codim_ >
auto
MappingElementAccessor<dim_ref_,codim_>::
get_face_inv_hessians_map(const Index face_id) const ->
const ValueVector<Derivatives<space_dim,face_dim,1,2>> &
{
    Assert(face_id < n_faces && face_id >= 0, ExcIndexRange(face_id,0,n_faces));
    Assert(face_values_[face_id].is_filled(), ExcMessage("The cache is not filled."));
    Assert(face_values_[face_id].fill_inv_hessians_, ExcNotInitialized());
    return face_values_[face_id].inv_hessians_;
}



template< int dim_ref_, int codim_ >
auto
MappingElementAccessor<dim_ref_,codim_>::
get_face_dets_map(const Index face_id) const -> const ValueVector<Real> &
{
    Assert(face_id < n_faces && face_id >= 0, ExcIndexRange(face_id,0,n_faces));
    Assert(face_values_[face_id].is_filled(), ExcMessage("The cache is not filled."));
    Assert(face_values_[face_id].fill_dets_, ExcNotInitialized());
    return face_values_[face_id].dets_;
}



template< int dim_ref_, int codim_ >
auto
MappingElementAccessor<dim_ref_,codim_>::
get_face_w_measures(const Index face_id) const -> const ValueVector<Real> &
{
    Assert(face_id < n_faces && face_id >= 0, ExcIndexRange(face_id,0,n_faces));
    Assert(face_values_[face_id].is_filled(), ExcMessage("The cache is not filled."));
    Assert(face_values_[face_id].fill_w_measures_, ExcNotInitialized());
    return face_values_[face_id].w_measures_;
}



template< int dim_ref_, int codim_ >
auto
MappingElementAccessor<dim_ref_,codim_>::
get_face_normals(const Index face_id) const -> const ValueVector<ValueMap> &
{
    Assert(face_id < n_faces && face_id >= 0, ExcIndexRange(face_id,0,n_faces));
    Assert(face_values_[face_id].is_filled(), ExcMessage("The cache is not filled."));
    Assert(face_values_[face_id].fill_normals_, ExcNotInitialized());
    return face_values_[face_id].normals_;
}



template< int dim_ref_, int codim_ >
auto
MappingElementAccessor<dim_ref_,codim_>::
transform_external_normals() const -> array< ValueVector<ValueMap>, codim >
{
    Assert(elem_values_.is_filled(), ExcMessage("The cache is not filled."));
    Assert(elem_values_.fill_gradients_, ExcNotInitialized());

    array<ValueVector<ValueMap>, codim> normals ;
    normals.fill(ValueVector<Point<space_dim>>(elem_values_.num_points_));


    for (Index i = 0; i < elem_values_.num_points_; i++)
    {
        Assert(false, ExcNotImplemented());
    }

    Assert(false, ExcNotImplemented());

    return normals;
}


template< int dim_ref_, int codim_ >
int
MappingElementAccessor<dim_ref_,codim_>::
get_num_points() const
{
    return elem_values_.num_points_ ;
}


template< int dim_ref_, int codim_ >
int
MappingElementAccessor<dim_ref_,codim_>::
get_num_face_points(const Index face_id) const
{
    return face_values_[face_id].num_points_ ;
}

template< int dim_ref_, int codim_ >
void
MappingElementAccessor<dim_ref_,codim_>::
print_info(LogStream &out) const
{
    using std::endl ;
    out << "MappingElementAccessor info" << endl ;

    out.push("\t") ;
    out << "num. points = " << elem_values_.num_points_ << endl ;
    out.pop() ;
}


template< int dim_ref_, int codim_ >
void
MappingElementAccessor<dim_ref_,codim_>::
print_memory_info(LogStream &out) const
{
    using std::endl ;
    out << "MappingElementAccessor memory info" << endl ;
    out << "this address = " << this << endl;

    out.push("\t") ;
    out << "mapping_ address = " << &mapping_ << endl ;
    out << "data_ memory address = " << &elem_values_ << endl ;
    out.pop() ;
}

IGA_NAMESPACE_CLOSE


#include <igatools/geometry/mapping_element_accessor.inst>


