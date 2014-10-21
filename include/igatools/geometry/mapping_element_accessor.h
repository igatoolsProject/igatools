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

#ifndef MAPPING_ELEMENT_ACCESSOR_H_
#define MAPPING_ELEMENT_ACCESSOR_H_

#if 0
#include <igatools/base/config.h>
#include <igatools/base/cache_status.h>
#include <igatools/base/value_flags_handler.h>
#include <igatools/utils/value_vector.h>
#include <igatools/utils/value_table.h>
#include <igatools/base/quadrature.h>
#include <igatools/geometry/topology.h>
#include <igatools/geometry/cartesian_grid_element_accessor.h>
#include <memory>

IGA_NAMESPACE_OPEN

template <int,int> class Mapping;

/**
 * See module on @ref accessors_iterators for a general overview.
 * @ingroup accessors
 *
 * @todo document me
 */
//TODO(pauletti, Sep 12, 2014): change dim_ref_ by dim_
template<int dim_ref_, int codim_ = 0>
class MappingElementAccessor
    : public CartesianGridElement<dim_ref_>
{
private:
    using self_t =  MappingElementAccessor<dim_ref_,codim_>;

public:
    /** Type required by the GridForwardIterator templated iterator */
    using ContainerType = const Mapping<dim_ref_,codim_>;

    using GridIterator = typename ContainerType::GridIterator;

    /** Dimension of the reference domain */
    using CartesianGridElement<dim_ref_>::dim;

    /** Codimension of the deformed domain */
    static const auto codim = ContainerType::codim;

    /** Dimension of the deformed domain embedding space. */
    static const auto space_dim = ContainerType::space_dim;

private:
    static const Size n_faces = UnitElement<dim>::faces_per_element;

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
    /** Fill flags supported by this iterator */
    static const ValueFlags admisible_flag =
        ValueFlags::measure |
        ValueFlags::w_measure |
        ValueFlags::face_point |
        ValueFlags::face_measure |
        ValueFlags::face_w_measure |
        ValueFlags::face_normal |
        ValueFlags::map_value |
        ValueFlags::map_gradient |
        ValueFlags::map_hessian |
        ValueFlags::map_inv_gradient |
        ValueFlags::map_inv_hessian |
        ValueFlags::map_face_value |
        ValueFlags::map_face_gradient |
        ValueFlags::map_face_hessian |
        ValueFlags::map_face_inv_gradient |
        ValueFlags::map_face_inv_hessian;

    /** @name Constructors */
    ///@{
    /**
     * Default constructor. Not allowed to be used.
     */
    MappingElementAccessor() = delete;

    /**
     * Constructs an accessor to element number @p index of a
     * Mapping.
     */
    MappingElementAccessor(const std::shared_ptr<ContainerType> mapping,
                           const Index index);

    MappingElementAccessor(const std::shared_ptr<ContainerType> mapping,
                           const TensorIndex<dim> &index);
    /**
     * Copy constructor.
     * It can be used with different copy policies (i.e. deep copy or shallow copy).
     * The default behaviour (i.e. using the proper interface of a classic copy constructor)
     * uses the deep copy.
     */
    MappingElementAccessor(const self_t &element, const CopyPolicy &copy_policy = CopyPolicy::deep);


    /**
     * Move constructor.
     */
    MappingElementAccessor(self_t &&element) = default;

    /**
     * Destructor.
     */
    ~MappingElementAccessor() = default;
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
     * @name Query information that requires the use of the cache
     */
    ///@{
    /**
     * Initializes the internal cache for the efficient
     * computation of the values requested in
     * the fill_flag on the given @p quadrature points.
     * This implies a uniform quadrature scheme
     * (i.e. the same for all elements).
     * @note This function should be called before fill_values()
     */
    void init_cache(const ValueFlags fill_flag,
                    const Quadrature<dim> &quadrature);

    /**
     * Initializes the internal cache for the efficient
     * computation of the values requested in
     * the fill_flag on the given @p quadrature points,
     * on the face specified by @p face_id.
     * This implies a uniform quadrature scheme
     * (i.e. the same for all elements).
     * @note This function should be called before fill_face_values()
     */
    void init_face_cache(const Index face_id,
                         const ValueFlags fill_flag,
                         const Quadrature<dim-1> &quadrature);

    /**
     * Fills the cache in accordance with the flag specifications used in the
     * init_values() function.
     *
     * Precondition Before invoking this function, you must call init_values().
     */
    void fill_cache();

    /**
     * Fills the cache in accordance with the flag specifications used in the
     * init_face_values() function.
     * @param face_id Face identifier.
     *
     * Precondition Before invoking this function, you must call
     * init_face_values().
     */
    void fill_face_cache(const Index face_id);
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
    void print_memory_info(LogStream &out) const;


    /** @name Functions for the mapping evaluations without the use of the cache */
    ///@{
    /**
     * Returns the value of the map
     * at each point (in the unit domain) specified by the input argument <tt>points</tt>.
     * @note This function does not use the cache and therefore can be called any time without
     * needing to pre-call init_values() / fill_values().
     * @warning The evaluation <tt>points</tt> must belong to the unit hypercube
     * \f$ [0,1]^{\text{dim}} \f$ otherwise, in Debug mode, an assertion will be raised.
     */
    ValueVector< ValueMap >
    evaluate_values_at_points(const ValueVector<Point> &points) const;

    /**
     * Returns the gradient of the map (i.e. the Jacobian)
     * at each point (in the unit domain) specified by the input argument <tt>points</tt>.
     * @note This function does not use the cache and therefore can be called any time without
     * needing to pre-call init_values() / fill_values().
     * @warning The evaluation <tt>points</tt> must belong to the unit hypercube
     * \f$ [0,1]^{\text{dim}} \f$ otherwise, in Debug mode, an assertion will be raised.
     */
    ValueVector< GradientMap >
    evaluate_gradients_at_points(const ValueVector<Point> &points) const;


    /**
     * Returns the hessian of the map
     * at each point (in the unit domain) specified by the input argument <tt>points</tt>.
     * @note This function does not use the cache and therefore can be called any time without
     * needing to pre-call init_values() / fill_values().
     * @warning The evaluation <tt>points</tt> must belong to the unit hypercube
     * \f$ [0,1]^{\text{dim}} \f$ otherwise, in Debug mode, an assertion will be raised.
     */
    ValueVector< HessianMap >
    evaluate_hessians_at_points(const ValueVector<Point> &points) const;

    ///@}


    void print_cache_info(LogStream &out) const;

private:
    // TODO (pauletti, Mar 21, 2014): Document this class
    class ValuesCache : public CacheStatus
    {
    public:
        void reset(const MappingElemValueFlagsHandler &flags_handler,
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

        MappingElemValueFlagsHandler flags_handler_;

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
        void reset(const MappingElemValueFlagsHandler &flags_handler,
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
    void copy_from(const MappingElementAccessor<dim_ref_,codim_> &element,
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

IGA_NAMESPACE_CLOSE
#endif
#endif // MAPPING_ELEMENT_ACCESSOR_H_
