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

#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/nurbs_space.h>
#include <igatools/geometry/ig_mapping.h>

using std::array;
using std::vector;
using std::shared_ptr;
using std::static_pointer_cast;
using std::dynamic_pointer_cast;

IGA_NAMESPACE_OPEN

template<int dim_ref_, int codim_ >
MappingElementAccessor<dim_ref_,codim_>::
MappingElementAccessor(const shared_ptr<ContainerType> mapping,
                       const Index index)
    :
    CartesianGridElementAccessor<dim>(mapping->get_grid(), index),
    mapping_(mapping)
{
    using BSplineSp = BSplineSpace<dim,dim+codim,1>;
    using BSplineMapping = IgMapping<BSplineSp>;
    if (dynamic_pointer_cast<const BSplineMapping>(mapping_))
    {
        auto ig_mapping = dynamic_pointer_cast<const BSplineMapping>(mapping_);
        mapping_.reset(new BSplineMapping(ig_mapping->get_data()));
    }



    using NURBSSp = NURBSSpace<dim,dim+codim,1>;
    using NURBSMapping = IgMapping<NURBSSp>;
    if (dynamic_pointer_cast<const NURBSMapping>(mapping_))
    {
        auto ig_mapping = dynamic_pointer_cast<const NURBSMapping>(mapping_);
        mapping_ .reset(new NURBSMapping(ig_mapping->get_data()));
    }

    Assert(mapping_->get_grid() != nullptr, ExcNullPtr());
}

template<>
MappingElementAccessor<0, 0>::
MappingElementAccessor(const shared_ptr<ContainerType> mapping,
                       const int index)
    :
    CartesianGridElementAccessor<0>(mapping->get_grid(), index),
    mapping_(mapping)
{}



template<int dim_ref_, int codim_ >
MappingElementAccessor<dim_ref_,codim_>::
MappingElementAccessor(const shared_ptr<ContainerType> mapping,
                       const TensorIndex<dim> &index)
    :
    CartesianGridElementAccessor<dim>(mapping->get_grid(), index),
    mapping_(mapping)
{
    using BSplineSp = BSplineSpace<dim,dim+codim,1>;
    using BSplineMapping = IgMapping<BSplineSp>;
    if (dynamic_pointer_cast<const BSplineMapping>(mapping_))
    {
        auto ig_mapping = dynamic_pointer_cast<const BSplineMapping>(mapping_);
        mapping_.reset(new BSplineMapping(ig_mapping->get_data()));
    }



    using NURBSSp = NURBSSpace<dim,dim+codim,1>;
    using NURBSMapping = IgMapping<NURBSSp>;
    if (dynamic_pointer_cast<const NURBSMapping>(mapping_))
    {
        auto ig_mapping = dynamic_pointer_cast<const NURBSMapping>(mapping_);
        mapping_ .reset(new NURBSMapping(ig_mapping->get_data()));
    }

    Assert(mapping_->get_grid() != nullptr, ExcNullPtr());
}



template<>
MappingElementAccessor<0, 0>::
MappingElementAccessor(const shared_ptr<ContainerType> mapping,
                       const TensorIndex<0> &index)
    :
    CartesianGridElementAccessor<0>(mapping->get_grid(), index),
    mapping_(mapping)
{}



template<int dim_ref_, int codim_ >
auto
MappingElementAccessor<dim_ref_,codim_>::
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

template< int dim_ref_, int codim_ >
void
MappingElementAccessor<dim_ref_,codim_>::
ValuesCache::
reset(const MappingElemValueFlagsHandler &flags_handler,
      const Quadrature<dim> &quad)
{
    flags_handler_ = flags_handler;

    this->quad_ = quad;
    this->num_points_ = this->quad_.get_num_points();

    if (flags_handler_.fill_values())
    {
        if (this->values_.size() != this->num_points_)
            this->values_.resize(this->num_points_);

        this->values_.zero();
    }
    else
    {
        this->values_.clear();
    }
    flags_handler_.set_values_filled(false);
    flags_handler_.set_points_filled(false);

    if (flags_handler_.fill_gradients())
    {
        if (this->gradients_.size() != this->num_points_)
            this->gradients_.resize(this->num_points_);

        this->gradients_.zero();
    }
    else
    {
        this->gradients_.clear();
    }
    flags_handler_.set_gradients_filled(false);

    if (flags_handler_.fill_hessians())
    {
        if (this->hessians_.size() != this->num_points_)
            this->hessians_.resize(this->num_points_);

        this->hessians_.zero();
    }
    else
    {
        this->hessians_.clear();
    }
    flags_handler_.set_hessians_filled(false);

    if (flags_handler_.fill_inv_gradients())
    {
        if (this->inv_gradients_.size() != this->num_points_)
            this->inv_gradients_.resize(this->num_points_);

        this->inv_gradients_.zero();
    }
    else
    {
        this->inv_gradients_.clear();
    }
    flags_handler_.set_inv_gradients_filled(false);

    if (flags_handler_.fill_inv_hessians())
    {
        if (this->inv_hessians_.size() != this->num_points_)
            this->inv_hessians_.resize(this->num_points_);

        this->inv_hessians_.zero();
    }
    else
    {
        this->inv_hessians_.clear();
    }
    flags_handler_.set_inv_hessians_filled(false);

    if (flags_handler_.fill_measures())
    {
        Assert(flags_handler_.fill_gradients(), ExcNotInitialized());
        if (this->measures_.size() != this->num_points_)
            this->measures_.resize(this->num_points_);

        this->measures_.zero();
    }
    else
    {
        this->measures_.clear();
    }
    flags_handler_.set_measures_filled(false);


    if (flags_handler_.fill_w_measures())
    {
        Assert(flags_handler_.fill_measures(), ExcNotInitialized());
        if (this->w_measures_.size() != this->num_points_)
            this->w_measures_.resize(this->num_points_);

        this->w_measures_.zero();
    }
    else
    {
        this->w_measures_.clear();
    }
    flags_handler_.set_w_measures_filled(false);
}






template< int dim_ref_, int codim_ >
void
MappingElementAccessor<dim_ref_,codim_>::
ElementValuesCache::
reset(const MappingElemValueFlagsHandler &flags_handler,
      const Quadrature<dim> &quad)
{
    ValuesCache::reset(flags_handler,quad);

    this->set_initialized(true);
}



template< int dim_ref_, int codim_ >
void
MappingElementAccessor<dim_ref_,codim_>::
FaceValuesCache::
reset(const Index face_id,
      const MappingFaceValueFlagsHandler &flags_handler,
      const Quadrature<dim> &quad)
{
    Assert(face_id < n_faces && face_id >= 0, ExcIndexRange(face_id,0,n_faces));

    ValuesCache::reset(flags_handler,quad.collapse_to_face(face_id));

    if (flags_handler.fill_normals())
    {
        if (normals_.size() != this->num_points_)
            normals_.resize(this->num_points_);

        normals_.zero();
    }
    else
    {
        normals_.clear();
    }

    this->set_initialized(true);
}



template< int dim_ref_, int codim_ >
void
MappingElementAccessor<dim_ref_,codim_>::
FaceValuesCache::
reset(const Index face_id,
      const MappingFaceValueFlagsHandler &flags_handler,
      const Quadrature<dim-1> &quad1)
{
    Assert(false, ExcNotImplemented());
}



template< int dim_ref_, int codim_ >
void
MappingElementAccessor<dim_ref_,codim_>::
init_cache(const ValueFlags fill_flag,
           const Quadrature<dim> &quad)
{
    Assert((fill_flag|admisible_flag) == admisible_flag,
           typename CartesianGridElementAccessor<dim_ref_>::ExcFillFlagNotSupported(admisible_flag, fill_flag));

    // initializing the cache of the CartesianGridElementAccessor
    {
        ValueFlags grid_flag = mapping_->required_flags();
        if (contains(fill_flag , ValueFlags::point))
            grid_flag |= ValueFlags::point;
        if (contains(fill_flag , ValueFlags::w_measure))
            grid_flag |= ValueFlags::w_measure;
        if (contains(fill_flag , ValueFlags::face_point))
            grid_flag |= ValueFlags::face_point;
        if (contains(fill_flag , ValueFlags::face_w_measure))
            grid_flag |= ValueFlags::face_w_measure;
        if (contains(fill_flag , ValueFlags::face_normal))
            grid_flag |= ValueFlags::face_point;
        CartesianGridElementAccessor<dim_ref_>::init_cache(grid_flag,quad);
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
    if (contains(f_flag , ValueFlags::face_normal))
    {
        f_flag |= ValueFlags::map_face_gradient;
        f_flag |= ValueFlags::map_face_inv_gradient;
    }

    const MappingElemValueFlagsHandler elem_flags_handler(f_flag);
    const MappingFaceValueFlagsHandler face_flags_handler(f_flag);
    Assert(!elem_flags_handler.fill_none() ||
           !face_flags_handler.fill_none(),
           ExcMessage("Nothing to initialize/reset."))


    if (!elem_flags_handler.fill_none())
        elem_values_.reset(elem_flags_handler, quad);

    if (!face_flags_handler.fill_none())
    {
        Index face_id = 0;
        for (auto &face_value : face_values_)
        {
            // TODO: this is temporary and must be removed.
            if (contains(f_flag , ValueFlags::face_normal))
                face_value.fill_normals_ = true;
            face_value.reset(face_id++, face_flags_handler, quad);
        }
    }

    mapping_->init_element(f_flag, quad);
}


template< int dim_ref_, int codim_ >
void
MappingElementAccessor<dim_ref_,codim_>::
init_face_cache(const Index face_id,
                const ValueFlags fill_flag,
                const Quadrature<dim-1> &quad)
{
    Assert(false, ExcNotImplemented());
    AssertThrow(false, ExcNotImplemented());
}


template< int dim_ref_, int codim_ >
void
MappingElementAccessor<dim_ref_,codim_>::
fill_cache()
{
    CartesianGridElementAccessor<dim_ref_>::fill_cache();
    mapping_->set_element(GridIterator(*this));

    Assert(elem_values_.is_initialized(), ExcNotInitialized());

    Assert(elem_values_.num_points_ ==
           elem_values_.quad_.get_points().flat_size(),
           ExcDimensionMismatch(elem_values_.num_points_,
                                elem_values_.quad_.get_points().flat_size()));

    if (elem_values_.flags_handler_.fill_values())
    {
        mapping_->evaluate(elem_values_.values_);
        elem_values_.flags_handler_.set_values_filled(true);
        elem_values_.flags_handler_.set_points_filled(true);
    }

    if (elem_values_.flags_handler_.fill_gradients())
    {
        mapping_->evaluate_gradients(elem_values_.gradients_);
        elem_values_.flags_handler_.set_gradients_filled(true);
    }

    if (elem_values_.flags_handler_.fill_hessians())
    {
        mapping_->evaluate_hessians(elem_values_.hessians_);
        elem_values_.flags_handler_.set_hessians_filled(true);
    }

    elem_values_.fill_composite_values();

    if (elem_values_.flags_handler_.fill_measures() || elem_values_.flags_handler_.fill_w_measures())
    {
        Assert(elem_values_.flags_handler_.gradients_filled(),ExcMessage("Gradients not filled."));
        for (Index i = 0; i < elem_values_.num_points_; i++)
            elem_values_.measures_[i] = determinant<dim,space_dim>(elem_values_.gradients_[i]);
        elem_values_.flags_handler_.set_measures_filled(true);

        if (elem_values_.flags_handler_.fill_w_measures())
        {
            Assert(elem_values_.flags_handler_.measures_filled(),ExcMessage("Measures not filled."));
            const ValueVector<Real> &dets_map = elem_values_.measures_;
            const auto weights = CartesianGridElementAccessor<dim_ref_>::get_w_measures();

            for (Index i = 0; i < elem_values_.num_points_; i++)
                elem_values_.w_measures_[i] = dets_map[i] * weights[i];
            elem_values_.flags_handler_.set_w_measures_filled(true);
        }
    }

    elem_values_.set_filled(true);
}



template< int dim_ref_, int codim_ >
void
MappingElementAccessor<dim_ref_,codim_>::
fill_face_cache(const Index face_id)
{
    Assert(face_id < n_faces && face_id >= 0, ExcIndexRange(face_id,0,n_faces));

    CartesianGridElementAccessor<dim_ref_>::fill_face_cache(face_id);
    mapping_->set_face_element(face_id, *this);

    auto &face_values = face_values_[face_id];

    const auto &num_points = face_values.num_points_;

    Assert(face_values.is_initialized(), ExcNotInitialized());

    Assert(num_points ==
           face_values.quad_.get_points().flat_size(),
           ExcDimensionMismatch(num_points,
                                face_values.quad_.get_points().flat_size()));

    if (face_values.flags_handler_.fill_values())
    {
        mapping_->evaluate_face(face_id, face_values.values_);
        face_values.flags_handler_.set_values_filled(true);
        face_values.flags_handler_.set_points_filled(true);
    }

    if (face_values.flags_handler_.fill_gradients())
    {
        mapping_->evaluate_face_gradients(face_id, face_values.gradients_);
        face_values.flags_handler_.set_gradients_filled(true);
    }

    if (face_values.flags_handler_.fill_hessians())
    {
        mapping_->evaluate_face_hessians(face_id, face_values.hessians_);
        face_values.flags_handler_.set_hessians_filled(true);
    }

    face_values.fill_composite_values();

    if (face_values.flags_handler_.fill_measures() || face_values.flags_handler_.fill_w_measures())
    {
        Assert(face_values.flags_handler_.gradients_filled(),ExcMessage("Gradients not filled."));

        const auto active_directions = UnitElement<dim>::face_active_directions[face_id];
        const auto face_dim = UnitElement<dim>::face_dim;
        Derivatives<face_dim, space_dim, 1, 1> face_gradient;
        for (Index i = 0; i < face_values.num_points_; i++)
        {
            auto &gradient = face_values.gradients_[i];
            for (int dir = 0; dir < face_dim; ++dir)
                face_gradient[dir] = gradient[active_directions[dir]];

            face_values.measures_[i] = determinant<face_dim,space_dim>(face_gradient);
        }
        face_values.flags_handler_.set_measures_filled(true);

        if (face_values.flags_handler_.fill_w_measures())
        {
            Assert(face_values.flags_handler_.measures_filled(),ExcMessage("Measures not filled."));
            const ValueVector<Real> &dets_map = face_values.measures_;
            const auto weights =
                CartesianGridElementAccessor<dim_ref_>::get_w_measures(FaceTopology<dim_ref_>(face_id));

            for (Index i = 0; i < face_values.num_points_; i++)
                face_values.w_measures_[i] = dets_map[i] * weights[i];
            face_values.flags_handler_.set_w_measures_filled(true);
        }
    }

    if (face_values.fill_normals_)
    {
        Assert(face_values.flags_handler_.inv_gradients_filled(),ExcMessage("Inverse gradients not filled."));
//        Assert(false, ExcMessage("The computation of face normals must be tested before used."));
//        AssertThrow(false, ExcMessage("The computation of face normals must be tested before used."));
        // Obtain n_hat from UnitElement
        Points<dim_ref_> n_hat = UnitElement<dim_ref_>::face_normal[face_id];
        for (Index i = 0; i < num_points; i++)
        {
            const auto DF_inv_t = co_tensor(transpose(face_values.inv_gradients_[i]));
            face_values.normals_[i] = action(DF_inv_t, n_hat);
            face_values.normals_[i] /= face_values.normals_[i].norm();
        }
        face_values.normals_filled_ = true;
    }

    face_values.set_filled(true);
}

template< int dim_ref_, int codim_ >
void
MappingElementAccessor<dim_ref_,codim_>::
ValuesCache::
fill_composite_values()
{
    if (flags_handler_.fill_inv_gradients())
    {
        Assert(flags_handler_.gradients_filled(),ExcMessage("Gradients not filled."));
        for (Index i = 0; i < num_points_; i++)
            MappingElementAccessor<dim_ref_,codim_>::evaluate_inverse_gradient(
                gradients_[i],
                inv_gradients_[i]);

        flags_handler_.set_inv_gradients_filled(true);
    }

    if (flags_handler_.fill_inv_hessians())
    {
        Assert(flags_handler_.hessians_filled(),ExcMessage("Hessians not filled."));
        Assert(flags_handler_.inv_gradients_filled(),ExcMessage("Hessians not filled."));

        for (Index i = 0; i < num_points_; i++)
        {
            /*
             * To fill the hessian of F{^-1}, we use the formula
             * D2F{^-1} [u][v] = - DF{^-1}[ D2F[ DF{^-1}[u] ][ DF{^-1}[v] ] ],
             * This formula can be obtained by differentiating the identity
             * DF * DF{^-1} = I
             */
            MappingElementAccessor<dim_ref_,codim_>::evaluate_inverse_hessian(
                hessians_[i],
                inv_gradients_[i],
                inv_hessians_[i]);
        }
        flags_handler_.set_inv_hessians_filled(true);

    }
}


template< int dim_ref_, int codim_ >
void
MappingElementAccessor<dim_ref_,codim_>::
evaluate_inverse_hessian(
    const HessianMap &D2F,
    const Derivatives<space_dim,dim,1,1> &DF_inv,
    Derivatives<space_dim,dim,1,2> &D2F_inv)
{
    for (int u=0; u<dim; ++u)
    {
        const auto tmp_u = action(D2F,DF_inv[u]);
        for (int v=0; v<dim; ++v)
        {
            const auto tmp_u_v = action(tmp_u,DF_inv[v]);

            D2F_inv[u][v] = - action(DF_inv,tmp_u_v);
        }
    }
}


template< int dim_ref_, int codim_ >
void
MappingElementAccessor<dim_ref_,codim_>::
evaluate_inverse_gradient(const GradientMap &DF, Derivatives<space_dim,dim,1,1> &DF_inv)
{
    inverse<dim,space_dim>(DF,DF_inv);
}

template< int dim_ref_, int codim_ >
auto
MappingElementAccessor<dim_ref_,codim_>::
get_map_values(const TopologyId<dim> &topology_id) const -> const ValueVector<ValueMap> &
{
    const auto &cache = this->get_values_cache(topology_id);
    Assert(cache.is_filled(), ExcCacheNotFilled());
    Assert(cache.flags_handler_.values_filled(), ExcMessage("Values not filled."));
    return cache.values_;
}

template< int dim_ref_, int codim_ >
auto
MappingElementAccessor<dim_ref_,codim_>::
get_face_values(const Index face_id) const -> const ValueVector<ValueMap> &
{
    return this->get_map_values(FaceTopology<dim>(face_id));
}


template< int dim_ref_, int codim_ >
auto
MappingElementAccessor<dim_ref_,codim_>::
get_map_gradients(const TopologyId<dim> &topology_id) const -> const ValueVector<GradientMap> &
{
    const auto &cache = this->get_values_cache(topology_id);
    Assert(cache.is_filled(), ExcCacheNotFilled());
    Assert(cache.flags_handler_.gradients_filled(), ExcMessage("Gradients not filled."));
    return cache.gradients_;
}



template< int dim_ref_, int codim_ >
auto
MappingElementAccessor<dim_ref_,codim_>::
get_map_hessians(const TopologyId<dim> &topology_id) const -> const ValueVector<HessianMap> &
{
    const auto &cache = this->get_values_cache(topology_id);
    Assert(cache.is_filled(), ExcCacheNotFilled());
    Assert(cache.flags_handler_.hessians_filled(), ExcMessage("Hessians not filled."));
    return cache.hessians_;
}



template< int dim_ref_, int codim_ >
auto
MappingElementAccessor<dim_ref_,codim_>::
get_inv_gradients(const TopologyId<dim> &topology_id) const -> const ValueVector<Derivatives<space_dim,dim,1,1>> &
{
    const auto &cache = this->get_values_cache(topology_id);
    Assert(cache.is_filled(), ExcCacheNotFilled());
    Assert(cache.flags_handler_.inv_gradients_filled(), ExcMessage("Inverse gradients not filled."));
    return cache.inv_gradients_;
}



template< int dim_ref_, int codim_ >
auto
MappingElementAccessor<dim_ref_,codim_>::
get_inv_hessians(const TopologyId<dim> &topology_id) const -> const ValueVector<Derivatives<space_dim,dim,1,2>> &
{
    const auto &cache = this->get_values_cache(topology_id);
    Assert(cache.is_filled(), ExcCacheNotFilled());
    Assert(cache.flags_handler_.inv_hessians_filled(), ExcMessage("Inverse hessians not filled."));
    return cache.inv_hessians_;
}



template< int dim_ref_, int codim_ >
auto
MappingElementAccessor<dim_ref_,codim_>::
get_measures(const TopologyId<dim> &topology_id) const -> const ValueVector<Real> &
{
    const auto &cache = this->get_values_cache(topology_id);
    Assert(cache.is_filled(), ExcCacheNotFilled());
    Assert(cache.flags_handler_.measures_filled(), ExcMessage("Measures not filled."));
    return cache.measures_;
}

template< int dim_ref_, int codim_ >
auto
MappingElementAccessor<dim_ref_,codim_>::
get_face_measures(const Index face_id) const -> const ValueVector<Real> &
{
    return this->get_measures(FaceTopology<dim>(face_id));
}


template< int dim_ref_, int codim_ >
auto
MappingElementAccessor<dim_ref_,codim_>::
get_w_measures(const TopologyId<dim> &topology_id) const -> const ValueVector<Real> &
{
    const auto &cache =this->get_values_cache(topology_id);
    Assert(cache.is_filled(), ExcCacheNotFilled());
    Assert(cache.flags_handler_.w_measures_filled(), ExcMessage("W*Measures not filled."));
    return cache.w_measures_;
}

template< int dim_ref_, int codim_ >
auto
MappingElementAccessor<dim_ref_,codim_>::
get_face_w_measures(const Index face_id) const -> const ValueVector<Real> &
{
    return this->get_w_measures(FaceTopology<dim>(face_id));
}


template< int dim_ref_, int codim_ >
auto
MappingElementAccessor<dim_ref_,codim_>::
get_face_normals(const Index face_id) const -> const ValueVector<ValueMap> &
{
    Assert(face_id < n_faces && face_id >= 0, ExcIndexRange(face_id,0,n_faces));
    Assert(face_values_[face_id].is_filled(), ExcCacheNotFilled());
    Assert(face_values_[face_id].normals_filled_, ExcMessage("Normals not filled."));
    return face_values_[face_id].normals_;
}



template< int dim_ref_, int codim_ >
auto
MappingElementAccessor<dim_ref_,codim_>::
transform_external_normals() const -> array< ValueVector<ValueMap>, codim >
{
    Assert(elem_values_.is_filled(), ExcMessage("The cache is not filled."));
    Assert(elem_values_.flags_handler_.fill_gradients(), ExcNotInitialized());

    array<ValueVector<ValueMap>, codim> normals;
    normals.fill(ValueVector<Points<space_dim>>(elem_values_.num_points_));


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
get_num_points(const TopologyId<dim> &topology_id) const
{
    const auto &cache =this->get_values_cache(topology_id);
    return cache.num_points_;
}





template< int dim_ref_, int codim_ >
auto
MappingElementAccessor<dim_ref_,codim_>::
evaluate_values_at_points(const ValueVector<Point> &points) const ->
ValueVector< ValueMap >
{
    const int n_points = points.size();
    Assert(n_points >= 0, ExcEmptyObject());

    const auto points_ref_domain = this->transform_points_unit_to_reference(points);

    ValueVector<ValueMap> map_value(n_points);

    mapping_->evaluate_at_points(points_ref_domain,map_value);

    return map_value;
}

template< int dim_ref_, int codim_ >
auto
MappingElementAccessor<dim_ref_,codim_>::
evaluate_gradients_at_points(const ValueVector<Point> &points) const ->
ValueVector< GradientMap >
{
    const int n_points = points.size();
    Assert(n_points >= 0, ExcEmptyObject());

    const auto points_ref_domain = this->transform_points_unit_to_reference(points);

    ValueVector<GradientMap> map_gradient(n_points);

    mapping_->evaluate_gradients_at_points(points_ref_domain,map_gradient);

    return map_gradient;
}

template< int dim_ref_, int codim_ >
auto
MappingElementAccessor<dim_ref_,codim_>::
evaluate_hessians_at_points(const ValueVector<Point> &points) const ->
ValueVector< HessianMap >
{
    const int n_points = points.size();
    Assert(n_points >= 0, ExcEmptyObject());

    const auto points_ref_domain = this->transform_points_unit_to_reference(points);

    ValueVector<HessianMap> map_hessian(n_points);

    mapping_->evaluate_hessians_at_points(points_ref_domain,map_hessian);

    return map_hessian;
}

//*/

template< int dim_ref_, int codim_ >
void
MappingElementAccessor<dim_ref_,codim_>::
print_info(LogStream &out,const VerbosityLevel verbosity_level) const
{
    using std::endl;

    std::string tab = "   ";

    out << "MappingElementAccessor info" << endl;

    out.push(tab);


    CartesianGridElementAccessor<dim_ref_>::print_info(out);
    out << "num. points = " << elem_values_.num_points_ << endl;


    if (contains(verbosity_level,VerbosityLevel::debug))
    {
        elem_values_.flags_handler_.print_info(out);

        for (int face_id = 0; face_id < n_faces; ++face_id)
        {
            face_values_[face_id].flags_handler_.print_info(out);
        }
    }


    out.pop();
}


template< int dim_ref_, int codim_ >
void
MappingElementAccessor<dim_ref_,codim_>::
print_memory_info(LogStream &out) const
{
    using std::endl;
    out << "MappingElementAccessor memory info" << endl;
    out << "this address = " << this << endl;

    out.push("\t");
    out << "mapping_ address = " << &mapping_ << endl;
    out << "data_ memory address = " << &elem_values_ << endl;
    out.pop();
}

IGA_NAMESPACE_CLOSE


#include <igatools/geometry/mapping_element_accessor.inst>


