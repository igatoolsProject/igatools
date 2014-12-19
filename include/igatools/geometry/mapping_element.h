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

#ifndef MAPPING_ELEMENT_H_
#define MAPPING_ELEMENT_H_

#include <igatools/utils/array.h>
#include <igatools/geometry/mapping.h>
#include <igatools/base/function_element.h>

IGA_NAMESPACE_OPEN


template<int dim_, int codim_ = 0>
class MappingElement
    : public FunctionElement<dim_, 0, dim_+codim_>
{
private:
    using self_t  = MappingElement<dim_, codim_>;
    using parent_t = FunctionElement<dim_, 0, dim_+codim_>;
    using Map = Mapping<dim_, codim_>;

public:
    using ContainerType = Map;
    static const int dim = dim_;
    static const int codim = codim_;
    static const int space_dim = dim_+codim_;


    /** @name Constructors */
    ///@{
    /**
     * Default constructor.
     */
    MappingElement()
    {
        Assert(false,ExcNotImplemented());
    }

    /**
     * Construct an accessor pointing to the element with
     * flat index @p elem_index of the CartesianGrid @p grid.
     */
    MappingElement(const std::shared_ptr<const CartesianGrid<dim_>> grid,
                   const Index elem_index)
        :
        parent_t(grid,elem_index)
    {}

    /**
     * Copy constructor.
     * It can be used with different copy policies
     * (i.e. deep copy or shallow copy).
     * The default behaviour (i.e. using the proper interface of a
     * classic copy constructor)
     * uses the deep copy.
     */
    MappingElement(const self_t &elem,
                   const CopyPolicy &copy_policy = CopyPolicy::deep)
    {
        Assert(false,ExcNotImplemented());
    }

    /**
     * Move constructor.
     */
    MappingElement(self_t &&elem) = default;

    /**
     * Destructor.
     */
    ~MappingElement() = default;
    ///@}

    template<int order>
    using InvDerivative = typename Map::template InvDerivative<order>;

    template <int order>
    using Derivative = typename Map::template Derivative<order>;

    template<int k>
    ValueVector<Real> const &get_measures(const int j) const
    {
        const auto &cache = local_cache_->template get_value_cache<k>(j);
        return cache.measures_;
    }

    template<int k>
    ValueVector<Real> const &get_w_measures(const int j) const
    {
        const auto &cache = local_cache_->template get_value_cache<k>(j);
        return cache.w_measures_;
    }

    template<int order = 0, int k = dim>
    auto
    get_inverse_values(const int j = 0) const
    {
        const auto &cache = local_cache_->template get_value_cache<k>(j);
        Assert(cache.is_filled() == true, ExcCacheNotFilled());
        return cache.template get_inv_values<order>();
    }


//    ValueVector<special_array<Points<space_dim>, codim>>
    ValueVector<Points<space_dim> > get_external_normals() const;

    using MetricTensor =
        Tensor<dim, 1, tensor::covariant, Tensor<dim, 1, tensor::contravariant, Tdouble> >;

    ValueVector<MetricTensor> compute_inv_first_fundamental_form() const;

    ValueVector<MetricTensor> compute_second_fundamental_form() const;

    ValueVector< Derivative<1> > get_D_external_normals() const;

    ValueVector<vector<Real> > get_principal_curvatures() const;


    template<int sub_dim>
    ValueVector<Points<space_dim> >
    get_boundary_normals(const int s_id) const
    {
        Assert(dim==sub_dim+1, ExcNotImplemented());
        ValueVector<Points<space_dim>> res;
        const auto &DF_inv = get_inverse_values<1, sub_dim>(s_id);
        const auto n_hat  = this->get_grid()->template get_boundary_normals<sub_dim>(s_id)[0];

        const auto n_points = DF_inv.get_num_points();
        res.resize(n_points);
        for (int i = 0; i< n_points; ++i)
        {
            const auto DF_inv_t = co_tensor(transpose(DF_inv[i]));
            res[i] = action(DF_inv_t, n_hat);
            res[i] /= res[i].norm();
        }
        return res;
    }

private:
    class ValuesCache : public CacheStatus
    {
    public:

        void resize(const MappingFlags &flags_handler,
                    const int n_points)
        {
            //TODO(pauletti, Oct 11, 2014): missing all necesary clears
            flags_handler_ = flags_handler;

            if (flags_handler_.fill_measures())
                measures_.resize(n_points);
            if (flags_handler_.fill_w_measures())
                w_measures_.resize(n_points);
            if (flags_handler_.fill_inv_gradients())
                std::get<1>(inv_values_).resize(n_points);
            if (flags_handler_.fill_inv_hessians())
                std::get<2>(inv_values_).resize(n_points);

            set_initialized(true);
        }

        void print_info(LogStream &out) const
        {
            flags_handler_.print_info(out);
            measures_.print_info(out);
        }


        template<int k>
        auto &get_inv_values()
        {
            return std::get<k>(inv_values_);
        }

        template<int k>
        const auto &get_inv_values() const
        {
            return std::get<k>(inv_values_);
        }

        MappingFlags flags_handler_;

        ValueVector<Real> measures_;
        ValueVector<Real> w_measures_;

        std::tuple<ValueVector<InvDerivative<0>>,
            ValueVector<InvDerivative<1>>,
            ValueVector<InvDerivative<2>>> inv_values_;


    };

    class LocalCache
    {
    public:
        LocalCache() = default;

        LocalCache(const LocalCache &in) = default;

        LocalCache(LocalCache &&in) = default;

        ~LocalCache() = default;


        LocalCache &operator=(const LocalCache &in) = delete;

        LocalCache &operator=(LocalCache &&in) = delete;

        void print_info(LogStream &out) const;

        template <int k>
        ValuesCache &
        get_value_cache(const int j)
        {
            return std::get<k>(values_)[j];
        }

        template <int k>
        const ValuesCache &
        get_value_cache(const int j) const
        {
            return std::get<k>(values_)[j];
        }

        CacheList<ValuesCache, dim> values_;
    };

    std::shared_ptr<LocalCache> local_cache_;

public:
    using CacheType = LocalCache;

private:
    template <class Accessor> friend class CartesianGridIteratorBase;
    friend class Mapping<dim, codim>;

    /**
     * Creates a new object performing a deep copy of the current object using the MappingElement
     * copy constructor.
     */
    std::shared_ptr<MappingElement<dim_,codim_> > clone() const
    {
        auto elem = std::shared_ptr<MappingElement<dim_,codim_> >(
                        new MappingElement(*this,CopyPolicy::deep));
        Assert(elem != nullptr, ExcNullPtr());
        return elem;
    }

};


IGA_NAMESPACE_CLOSE

#endif // MAPPING_ELEMENT_ACCESSOR_H_


#if 0
/** Type required by the CartesianGridIterator templated iterator */
using ContainerType = const Mapping<dim_,codim_>;

using GridIterator = typename ContainerType::GridIterator;

/** Dimension of the reference domain */
using CartesianGridElement<dim_>::dim;

/** Codimension of the deformed domain */
static const auto codim = ContainerType::codim;

/** Dimension of the deformed domain embedding space. */
static const auto space_dim = ContainerType::space_dim;

private:
static const Size n_faces = UnitElement<dim>::n_faces;

public:
/**
 * see Mapping<dim, codim>::Value
 */
using Point       = typename ContainerType::Point;

//TODO(pauletti, Jun 21, 2014): we should use Value instead of ValueMap
using ValueMap    = typename ContainerType::Value;
using GradientMap = typename ContainerType::Gradient;
using HessianMap  = typename ContainerType::Hessian;

public:


private:
/**
 * Copy constructor.
 * It can be used with different copy policies (i.e. deep copy or shallow copy).
 * The default behaviour (i.e. using the proper interface of a classic copy constructor)
 * uses the deep copy.
 */
MappingElement(const self_t &element, const CopyPolicy &copy_policy = CopyPolicy::deep);

public:

/**
 * Move constructor.
 */
MappingElement(self_t &&element) = default;

/**
 * Destructor.
 */
~MappingElement() = default;
///@}

/** @name Assignment operators */
///@{
/**
 * Copy assignment operator. Performs a <b>shallow copy</b> of the input @p element.
 *
 * @note Internally it uses the function shallow_copy_from().
 */
self_t &operator=(const self_t &element);

/**
 * Move assignment operator.
 */
self_t &operator=(self_t &&element) = default;
///@}


/**
 * @name Functions for performing different kind of copy.
 */
///@{
/**
 * Performs a deep copy of the input @p element,
 * i.e. a new local cache is built using the copy constructor on the local cache of @p element.
 *
 * @note In DEBUG mode, an assertion will be raised if the input local cache is not allocated.
 */
void deep_copy_from(const self_t &element);


/**
 * Performs a shallow copy of the input @p element. The current object will contain a pointer to the
 * local cache used by the input @p element.
 */
void shallow_copy_from(const self_t &element);
///@}


/**
 * @name Getting the mapping-related values stored in the cache.
 * @note In order to use these functions, the cache of the
 * MappingElementAccessor must be properly filled.
 */
///@{
/** Returns the value of the map at the dilated quadrature points.*/
const ValueVector<ValueMap> &
get_map_values(const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

/**
 * Returns the value of the map at the dilated quadrature points
 * on the face specified by @p face_id.
 */
const ValueVector<ValueMap> &
get_face_values(const Index face_id) const;

/** Returns the gradient of the map at the dilated quadrature points.*/
const ValueVector<GradientMap> &
get_map_gradients(const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

/** Returns the hessian of the map at the dilated quadrature points. */
const ValueVector<HessianMap> &
get_map_hessians(const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

/** Returns the inverse of the gradient of the map at the dilated quadrature points. */
const ValueVector< Derivatives< space_dim, dim,1,1 > > &
get_inv_gradients(const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

/** Returns the inverse of the hessian of the map at the dilated quadrature points. */
const ValueVector< Derivatives< space_dim, dim,1,2 > > &
get_inv_hessians(const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

/** Returns the gradient determinant of the map at the dilated quadrature points. */
const ValueVector< Real > &
get_measures(const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

/**
 * Returns the gradient determinant of the map at the dilated quadrature points
 * on the face specified by @p face_id.
 */
const ValueVector<Real> &get_face_measures(const Index face_id) const;

/**
 * Returns the quadrature weights multiplied by the
 * gradient determinant of the map at the dilated quadrature points.
 */
const ValueVector<Real> &get_w_measures(const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

/**
 * Returns the quadrature weights multiplied by the
 * gradient determinant of the map at the dilated quadrature points
 * on the face specified by @p face_id.
 */
const ValueVector<Real> &get_face_w_measures(const Index face_id) const;

/**
 * Returns the face normals for every quadrature point for the
 * specified face.
 */
const ValueVector< ValueMap > &get_face_normals(const Index face_id) const;
///@}

/**
 * Returns the number of evaluation points currently used
 * in the element cache.
 */
Size get_num_points(const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

/**
 * Prints some internal information.
 * @note Mostly used for testing and debugging.
 */
void print_info(LogStream &out,const VerbosityLevel verbosity_level = VerbosityLevel::normal) const;

/**
 * Prints some internal memory information.
 * @note Mostly used for testing and debugging.
 */
void print_cache_info(LogStream &out) const;


private:
// TODO (pauletti, Mar 21, 2014): Document this class
class ValuesCache : public CacheStatus
{
public:
    void reset(const MappingFlags &flags_handler,
               const Quadrature<dim> &quad);

    //TODO: the next member variables should be protected
public:
    /**
     * Fills the following cache values in accordance with the
     * flag specifications used in the reset() function.
     *
     * @pre Before invoking this function, values_, gradients_ and hessians_
     * must be properly filled.
     */
    void fill_composite_values();

    MappingFlags flags_handler_;

    ValueVector< ValueMap > values_;
    ValueVector< GradientMap > gradients_;
    ValueVector< HessianMap > hessians_;
    ValueVector< Derivatives< space_dim,dim,1,1 > > inv_gradients_;
    ValueVector< Derivatives< space_dim,dim,1,2 > > inv_hessians_;
    ValueVector< Real > measures_;
    ValueVector< Real > w_measures_;

    Size num_points_ = 0;
    Quadrature<dim> quad_;

    void print_info(LogStream &out) const
    {
        out.begin_item("Fill flags:");
        flags_handler_.print_info(out);
        out.end_item();

        out.begin_item("Values:");
        values_.print_info(out);
        out.end_item();

        out.begin_item("Gradients:");
        gradients_.print_info(out);
        out.end_item();

        out.begin_item("Hessians:");
        hessians_.print_info(out);
        out.end_item();

        out.begin_item("Inverse Gradients:");
        inv_gradients_.print_info(out);
        out.end_item();

        out.begin_item("Inverse Hessians:");
        inv_hessians_.print_info(out);
        out.end_item();

        out.begin_item("Measures:");
        measures_.print_info(out);
        out.end_item();

        out.begin_item("Weights * Measures:");
        w_measures_.print_info(out);
        out.end_item();
    }

};

/**
 * Structure to cache mapping information used by the transform_* member
 * functions.
 *
 */
struct ElementValuesCache : ValuesCache
{
    void reset(const MappingFlags &flags_handler,
               const Quadrature<dim> &quad);

};

// TODO (pauletti, Mar 21, 2014): Document this class
struct FaceValuesCache : ValuesCache
{
    void reset(const Index face_id,
               const MappingFaceValueFlagsHandler &flags_handler,
               const Quadrature<dim> &quad);

    void reset(const Index face_id,
               const MappingFaceValueFlagsHandler &flags_handler,
               const Quadrature<dim-1> &quad);


    ValueVector< ValueMap > normals_;
    bool fill_normals_ = false;
    bool normals_filled_ = false;

    void print_info(LogStream &out) const
    {
        ValuesCache::print_info(out);

        out.begin_item("Normals:");
        normals_.print_info(out);
        out.end_item();
    }
};

const ValuesCache &get_values_cache(const TopologyId<dim> &topology_id) const;


class LocalCache
{
public:
    LocalCache() = default;

    LocalCache(const LocalCache &in) = default;
    LocalCache(LocalCache &&in) = default;

    ~LocalCache() = default;


    LocalCache &operator=(const LocalCache &in) = delete;
    LocalCache &operator=(LocalCache &&in) = delete;

    void print_info(LogStream &out) const;

    /** Element values cache */
    ElementValuesCache elem_values_;

    /** Face values cache */
    std::array<FaceValuesCache, n_faces> face_values_;

};

std::shared_ptr<LocalCache> local_cache_;



std::shared_ptr<ContainerType> mapping_;

/**
 * @todo implement this function
 *TODO
 */
std::array<ValueVector<ValueMap>, codim> transform_external_normals() const;

protected:
/**
 * Performs a copy of the input @p element.
 * The type of copy (deep or shallow) is specified by the input parameter @p copy_policy.
 */
void copy_from(const MappingElementAccessor<dim_,codim_> &element,
               const CopyPolicy &copy_policy);


/**
 * Evaluates the gradient of F^{-1} (also when the mapping has codim > 0) using the formula
 * D(F^{-1}) = (DF^t * DF)^{-1} * DF^t
 */
static void evaluate_inverse_gradient(const GradientMap &DF, Derivatives<space_dim,dim,1,1> &DF_inv);

/**
 * Evaluates the inverse hessian of F^{-1} (also when the mapping has codim > 0) using the formula
 * D2F{^-1} [u][v] = - DF{^-1}[ D2F[ DF{^-1}[u] ][ DF{^-1}[v] ] ],
 * This formula can be obtained by differentiating the identity
 * DF * DF{^-1} = I
 */
static void evaluate_inverse_hessian(const HessianMap &D2F,
                                     const Derivatives<space_dim,dim,1,1> &DF_inv,
                                     Derivatives<space_dim,dim,1,2> &D2F_inv);
};

#endif

