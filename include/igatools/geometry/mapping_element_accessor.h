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
 * @ingroup accessors_iterators
 *
 * @todo document me
 */
template<int dim_ref_, int codim_>
class MappingElementAccessor
    : public CartesianGridElementAccessor<dim_ref_>
{
private:
    using self_t =  MappingElementAccessor<dim_ref_,codim_>;

public:
    /** Type required by the GridForwardIterator templated iterator */
    using ContainerType = Mapping<dim_ref_,codim_>;

    /** Dimension of the reference domain */
    using CartesianGridElementAccessor<dim_ref_>::dim;

    /** Codimension of the deformed domain */
    static const auto codim = ContainerType::codim;

    /** Dimension of the deformed domain embedding space. */
    static const auto space_dim = ContainerType::space_dim;

    // TODO (pauletti, Mar 21, 2014): why do we need this? should it be private?
    /** Dimension of the face.*/
    static const auto face_dim = ContainerType::face_dim ;

    // TODO (pauletti, Mar 21, 2014): should this be private?
    /**
     * see UnitElement<dim_>::faces_per_element
     */
    static const Size n_faces = UnitElement<dim>::faces_per_element;

    /**
     * see Mapping<dim, codim>::Value
     */
    using ValueMap        = typename ContainerType::ValueType;
    using GradientMap     = typename ContainerType::GradientType;
    using HessianMap      = typename ContainerType::HessianType;

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
        ValueFlags::map_face_inv_hessian ;

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
    MappingElementAccessor(Mapping<dim,codim> &mapping,
                           const int index);

    /**
     * Copy constructor.
     * Performs a deep copy of the MappingElementAccessor object.
     * Its cache is also deeply copied.
     */
    MappingElementAccessor(const self_t &element) = default;

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
     * Copy assignment operator.
     * @note Performs a deep copy of the MappingElementAccessor object.
     * Its cache is also deeply copied.
     */
    self_t &operator=(const self_t &element) = default;

    /**
     * Move assignment operator.
     */
    self_t &operator=(self_t &&element) = default;
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
    void init_values(const ValueFlags fill_flag,
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
    void init_face_values(const Index face_id,
                          const ValueFlags fill_flag,
                          const Quadrature<dim-1> &quadrature);

    /**
     * Fills the cache in accordance with the flag specifications used in the
     * init_values() function.
     *
     * Precondition Before invoking this function, you must call init_values().
     */
    void fill_values();

    /**
     * Fills the cache in accordance with the flag specifications used in the
     * init_face_values() function.
     * @param face_id Face identifier.
     *
     * Precondition Before invoking this function, you must call
     * init_face_values().
     */
    void fill_face_values(const Index face_id);
    ///@}

    /**
     * @name Getting the mapping-related values stored in the cache.
     * @note In order to use these functions, the cache of the MappingElementAccessor must be properly filled.
     */
    ///@{
    /** Returns the value of the map at the dilated quadrature points.*/
    const ValueVector<ValueMap> &
    get_values(const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

    /**
     * Returns the value of the map at the dilated quadrature points
     * on the face specified by @p face_id.
     */
    const ValueVector<ValueMap> &
    get_face_values(const Index face_id) const;

    /** Returns the gradient of the map at the dilated quadrature points.*/
    const ValueVector<GradientMap> &
    get_gradients(const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

    /** Returns the hessian of the map at the dilated quadrature points. */
    const ValueVector<HessianMap> &
    get_hessians(const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

    /** Returns the inverse of the gradient of the map at the dilated quadrature points. */
    const ValueVector< Derivatives< space_dim, dim,1,1 > > &
    get_inv_gradients(const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

    /** Returns the inverse of the hessian of the map at the dilated quadrature points. */
    const ValueVector< Derivatives< space_dim, dim,1,2 > > &
    get_inv_hessians(const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

    /** Returns the gradient determinant of the map at the dilated quadrature points. */
    const ValueVector< Real > &
    get_dets(const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

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

//        std::shared_ptr<MappingFaceValueFlagsHandler> get_flags_handler() const;
    };

    const ValuesCache &get_values_cache(const TopologyId<dim> &topology_id) const;

    ElementValuesCache elem_values_;

    std::array<FaceValuesCache, n_faces> face_values_;

    Mapping<dim, codim> *mapping_;

    /**
     * @todo implement this function
     *TODO
     */
    std::array<ValueVector<ValueMap>, codim> transform_external_normals() const;
};

IGA_NAMESPACE_CLOSE

#endif // MAPPING_ELEMENT_ACCESSOR_H_
