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

#ifndef MAPPING_ELEMENT_H_
#define MAPPING_ELEMENT_H_
#if 0
#include <igatools/utils/safe_stl_array.h>
#include <igatools/geometry/mapping.h>
#include <igatools/functions/function_element.h>

IGA_NAMESPACE_OPEN

/**
 *
 * @ingroup elements
 */
template<int dim_, int codim_ = 0>
class MappingElement
{
private:
    using self_t  = MappingElement<dim_,codim_>;
    using Map = Mapping<dim_,codim_>;
    using Func = MapFunction_new<dim_,codim_>;

    using Point = typename Func::Point;
    using Value = typename Func::Value;
    using Gradient = typename Func::Gradient;
    using Hessian  = typename Func::Hessian;
    using Div      = typename Func::Div;

public:
    using ContainerType = Map;
    static const int dim = dim_;
    static const int codim = codim_;
    static const int space_dim = dim_+codim_;


    /** @name Constructors */
    ///@{
    /**
     * Default constructor. It does nothing but it is needed for the
     * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     * mechanism.
     */
    MappingElement() = default;

    /**
     * Construct an accessor pointing to the element with
     * flat index @p elem_index of the Function @p func.
     */
    MappingElement(const std::shared_ptr<const Func> func,
                   const Index elem_index);

    /**
     * Copy constructor.
     * It can be used with different copy policies
     * (i.e. deep copy or shallow copy).
     * The default behaviour (i.e. using the proper interface of a
     * classic copy constructor)
     * uses the deep copy.
     */
    MappingElement(const self_t &elem,
                   const CopyPolicy &copy_policy = CopyPolicy::deep);

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

    template <class ValueType, int topology_dim = dim>
    auto &get_values_from_cache(const int topology_id = 0) const
    {
        Assert(local_cache_ != nullptr,ExcNullPtr());
        const auto &cache = local_cache_->template get_sub_elem_cache<topology_dim>(topology_id);
        return cache.template get_data<ValueType>();
    }

    template<int k>
    ValueVector<Real> const &get_measures(const int j) const
    {
        return get_values_from_cache<_Measure,k>(j);
    }

    template<int k>
    ValueVector<Real> const &get_w_measures(const int j) const
    {
        return get_values_from_cache<_W_Measure,k>(j);
    }

    const ValueVector<Points<space_dim> > &get_external_normals() const;

    using MetricTensor =
        Tensor<dim, 1, tensor::covariant, Tensor<dim, 1, tensor::contravariant, Tdouble> >;

    ValueVector<MetricTensor> compute_inv_first_fundamental_form() const;

    ValueVector<MetricTensor> compute_second_fundamental_form() const;

    ValueVector< Derivative<1> > get_D_external_normals() const;

    const ValueVector<SafeSTLVector<Real> > &get_principal_curvatures() const;


    template<int sub_dim>
    const ValueVector<Points<space_dim> > &
    get_boundary_normals(const int s_id) const
    {
#if 0
        Assert(dim==sub_dim+1, ExcNotImplemented());
        ValueVector<Points<space_dim>> res;
        const auto &DF_inv = get_values_from_cache<_InvGradient, sub_dim>(s_id);
        const auto n_hat  = this->get_grid()->template get_boundary_normals<sub_dim>(s_id)[0];

        const auto n_points = DF_inv.get_num_points();
        res.resize(n_points);
        for (int pt = 0; pt < n_points; ++pt)
        {
            const auto DF_inv_t = co_tensor(transpose(DF_inv[pt]));
            res[pt] = action(DF_inv_t, n_hat);
            res[pt] /= res[pt].norm();
        }
        return res;
#endif
        return get_values_from_cache<_BoundaryNormal,sub_dim>(s_id);
    }

    template<class ValueType, int k>
    auto
    get_values(const int j) const
    {
        Assert(local_cache_ != nullptr,ExcNullPtr());
        const auto &cache = local_cache_->template get_sub_elem_cache<k>(j);
        return cache.template get_data<ValueType>();
    }


    /**
     * @name Comparison operators.
     *
     * @brief The comparison operators compares the <em>position</em> of the element in the grid.
     *
     * @warning To be comparable, two Mapping objects must be defined using the same Function
     * (and therefore on the same grid),
     * otherwise an assertion will be raised (in Debug mode).
     */
    ///@{
    /** Returns TRUE if the two elements have the same index on the grid. */
    bool operator==(const self_t &a) const;


    /** Returns TRUE if the two elements have different indices on the grid. */
    bool operator!=(const self_t &a) const;

    /**
     * Returns TRUE if the the index of the element on the left of the operator <tt> < </tt>
     * is smaller than the the index of the element on the right.
     * */
    bool operator<(const self_t &a) const;

    /**
     * Returns TRUE if the the index of the element on the left of the operator <tt> < </tt>
     * is bigger than the the index of the element on the right.
     * */
    bool operator>(const self_t &a) const;
    ///@}

    /**
     * Sets the index of the element using the flatten representation.
     * @note This function also updates the index for the tensor representation.
     * @warning This may be a dangerous function, be careful when using it
     * as it is easy to use incorrectly. Only use it if you know what you
     * are doing.
     */
    void move_to(const Index flat_index) ;


    /** @name Functions related to the indices of the element in the cartesian grid. */
    ///@{
    /** Returns the index of the element in its flatten representation. */
    Index get_flat_index() const;

    /** Returns the index of the element in its tensor representation. */
    TensorIndex<dim> get_tensor_index() const;
    ///@}

    /** Return the cartesian grid from which the element belongs.*/
    std::shared_ptr<const CartesianGrid<dim> > get_grid() const;

private:

    using CType = boost::fusion::map<
                  boost::fusion::pair<         _Point,DataWithFlagStatus<ValueVector<Points<space_dim>>>>,
                  boost::fusion::pair<      _Gradient,DataWithFlagStatus<ValueVector<Derivative<1>>>>,
                  boost::fusion::pair<       _Hessian,DataWithFlagStatus<ValueVector<Derivative<2>>>>,
                  boost::fusion::pair<       _Measure,DataWithFlagStatus<ValueVector<Real>>>,
                  boost::fusion::pair<     _W_Measure,DataWithFlagStatus<ValueVector<Real>>>,
                  boost::fusion::pair<   _InvGradient,DataWithFlagStatus<ValueVector<InvDerivative<1>>>>,
                  boost::fusion::pair<    _InvHessian,DataWithFlagStatus<ValueVector<InvDerivative<2>>>>,
                  boost::fusion::pair<_BoundaryNormal,DataWithFlagStatus<ValueVector<Points<space_dim>>>>,
                  boost::fusion::pair<   _OuterNormal,DataWithFlagStatus<ValueVector<Points<space_dim>>>>,
                  boost::fusion::pair<     _Curvature,DataWithFlagStatus<ValueVector<SafeSTLVector<Real>>>>
                  >;


    /**
     * Returns the flags that are valid to be used with this class.
     *
     * @note The valid flags are defined to be the ones that can be inferred from the ValueType(s)
     * used as key of the boost::fusion::map in CType.
     */
    static ValueFlags get_valid_flags()
    {
        return cacheutils::get_valid_flags_from_cache_type(CType());
    }

    using Cache = FuncValuesCache<dim,CType>;


public:
    using CacheType = AllSubElementsCache<Cache>;



private:

    using FuncElem = FunctionElement<dim_, 0, dim_+codim_>;

    std::shared_ptr<FuncElem> func_elem_;

    std::shared_ptr<CacheType> local_cache_;


    template <class Accessor> friend class CartesianGridIteratorBase;
    friend class Mapping<dim, codim>;

    /**
     * Creates a new object performing a deep copy of the current object using the MappingElement
     * copy constructor.
     */
    std::shared_ptr<MappingElement<dim_,codim_> > clone() const;

public:
    FuncElem &get_func_element();

    const FuncElem &get_func_element() const;

    void print_info(LogStream &out) const;

    void print_cache_info(LogStream &out) const;

private:

#ifdef SERIALIZATION
    /**
     * @name Functions needed for boost::serialization
     * @see <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     */
    ///@{
    friend class boost::serialization::access;

    template<class Archive>
    void
    serialize(Archive &ar, const unsigned int version)
    {
        Assert(false,ExcNotImplemented());
        AssertThrow(false,ExcNotImplemented());
    }
    ///@}
#endif

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


class AllSubElementsCache
{
public:
    AllSubElementsCache() = default;

    AllSubElementsCache(const AllSubElementsCache &in) = default;
    AllSubElementsCache(AllSubElementsCache &&in) = default;

    ~AllSubElementsCache() = default;


    AllSubElementsCache &operator=(const AllSubElementsCache &in) = delete;
    AllSubElementsCache &operator=(AllSubElementsCache &&in) = delete;

    void print_info(LogStream &out) const;

    /** Element values cache */
    ElementValuesCache elem_values_;

    /** Face values cache */
    SafeSTLArray<FaceValuesCache, n_faces> face_values_;

};

std::shared_ptr<AllSubElementsCache> local_cache_;



std::shared_ptr<ContainerType> mapping_;

/**
 * @todo implement this function
 *TODO
 */
SafeSTLArray<ValueVector<ValueMap>, codim> transform_external_normals() const;

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
#endif
