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
#include <igatools/utils/value_vector.h>
#include <igatools/utils/value_table.h>
#include <igatools/base/quadrature.h>
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
public:

    /** Type required by the GridForwardIterator templated iterator */
    using ContainerType = Mapping<dim_ref_,codim_>;


    /** Dimension of the reference domain */
    using CartesianGridElementAccessor<dim_ref_>::dim;


    /** Codimension of the deformed domain */
    static const auto codim = ContainerType::codim;


    /** Dimension of the deformed domain embedding space. */
    static const auto space_dim = ContainerType::space_dim;


    /** Dimension of the face.*/
    static const auto face_dim = ContainerType::face_dim ;

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
    using GradientFaceMap = typename ContainerType::GradientFaceType;
    using HessianFaceMap  = typename ContainerType::HessianFaceType;


public:
    /** Fill flags supported by this iterator */
    static const ValueFlags admisible_flag =
        ValueFlags::point|
        ValueFlags::measure |
        ValueFlags::w_measure |
        ValueFlags::face_point |
        ValueFlags::face_measure |
        ValueFlags::face_w_measure |
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
    MappingElementAccessor(const MappingElementAccessor<dim_ref_,codim_> &element) = default;

    /**
     * Move constructor.
     */
    MappingElementAccessor(MappingElementAccessor<dim_ref_,codim_> &&element) = default;

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
    MappingElementAccessor<dim_ref_,codim_> &
    operator=(const MappingElementAccessor<dim_ref_,codim_> &element) = default;

    /**
     * Move assignment operator.
     */
    MappingElementAccessor<dim_ref_,codim_> &
    operator=(MappingElementAccessor<dim_ref_,codim_> &&element) = default;
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
     * @precondition Before invoking this function, you must call init_values().
     */
    void fill_values();


    /**
     * Fills the cache in accordance with the flag specifications used in the
     * init_face_values() function.
     * @param face_id Face identifier.
     * @precondition Before invoking this function, you must call init_face_values().
     */
    void fill_face_values(const Index face_id);
    ///@}

    /**
     * @name Getting the mapping-related values stored in the cache.
     * @note In order to use these functions, the cache of the MappingElementAccessor must be properly filled.
     */
    ///@{
    /** Returns the value of the map at the dilated quadrature points.*/
    const ValueVector<ValueMap> &get_values_map() const;


    /** Returns the gradient of the map at the dilated quadrature points.*/
    const ValueVector<GradientMap> &get_gradients_map() const;


    /** Returns the hessian of the map at the dilated quadrature points. */
    const ValueVector<HessianMap> &get_hessians_map() const;


    /** Returns the inverse of the gradient of the map at the dilated quadrature points. */
    const ValueVector< Derivatives< space_dim, dim,1,1 > > &get_inv_gradients_map() const;


    /** Returns the inverse of the hessian of the map at the dilated quadrature points. */
    const ValueVector< Derivatives< space_dim, dim,1,2 > > &get_inv_hessians_map() const;


    /** Returns the gradient determinant of the map at the dilated quadrature points. */
    const ValueVector< Real > &get_dets_map() const;


    /**
     * Returns the quadrature weights multiplied by the
     * gradient determinant of the map at the dilated quadrature points.
     */
    const ValueVector< Real > &get_w_measures() const;


    /**
     * Returns the face value of the map at the dilated quadrature points
     * at the specified face.
     */
    const ValueVector<ValueMap> &get_face_values_map(const Index face_id) const;


    /**
     * Returns the face gradient of the map at the dilated quadrature points
     * at the specified face.
     */
    const ValueVector<GradientFaceMap> &get_face_gradients_map(const Index face_id) const;


    /**
     * Returns the face hessian of the map at the dilated quadrature points
     * at the specified face.
     */
    const ValueVector<HessianFaceMap> &get_face_hessians_map(const Index face_id) const;


    /**
     * Returns the inverse of the face gradient of the map at the dilated quadrature points
     * at the specified face.
     */
    const ValueVector< Derivatives<space_dim,face_dim,1,1> > &get_face_inv_gradients_map(const Index face_id) const;


    /**
     * Returns the inverse of the face hessian of the map at the dilated quadrature points
     * at the specified face.
     */
    const ValueVector< Derivatives<space_dim,face_dim,1,2> > &get_face_inv_hessians_map(const Index face_id) const;

    /**
     * Returns the face gradient determinant of the map at the dilated quadrature points.
     */
    const ValueVector<Real> &get_face_dets_map(const Index face_id) const;

    /**
     * Returns the quadrature weights multiplied by the
     * gradient determinant of the face map at the dilated quadrature points.
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
    Size get_num_points() const;

    /**
     * Returns the number of evaluation points currently used
     * in the face cache.
     */
    Size get_num_face_points(const Index face_id) const;

    /**
     * Prints some internal information.
     * @note Mostly used for testing and debugging.
     */
    void print_info(LogStream &out) const;

    /**
     * Prints some internal memory information.
     * @note Mostly used for testing and debugging.
     */
    void print_memory_info(LogStream &out) const;

private:

    /**
     * Structure to cache mapping information used by the transform_* member
     * functions.
     *
     */
    struct ElementValuesCache : CacheStatus
    {
        void reset(const ValueFlags fill_flag,
                   const Quadrature<dim> &quad);

        Size num_points_ = 0;

        bool fill_values_ = false;
        bool fill_gradients_ = false;
        bool fill_hessians_ = false;
        bool fill_inv_gradients_ = false;
        bool fill_inv_hessians_ = false;
        bool fill_dets_ = false;
        bool fill_w_measures_ = false;

        ValueVector< ValueMap > values_;
        ValueVector< GradientMap > gradients_;
        ValueVector< HessianMap > hessians_;
        ValueVector< Derivatives< space_dim, dim,1,1 > > inv_gradients_;
        ValueVector< Derivatives< space_dim, dim,1,2 > > inv_hessians_;
        ValueVector< Real > dets_;
        ValueVector< Real > w_measures_;

        Quadrature<dim> quad_;
    };

    struct FaceValuesCache : CacheStatus
    {
        void reset(const Index face_id,
                   const ValueFlags fill_flag,
                   const Quadrature<dim> &quad);

        void reset(const Index face_id,
                   const ValueFlags fill_flag,
                   const Quadrature<dim-1> &quad);

        Size num_points_ = 0;

        bool fill_values_ = false;
        bool fill_gradients_ = false;
        bool fill_hessians_ = false;
        bool fill_inv_gradients_ = false;
        bool fill_inv_hessians_ = false;
        bool fill_dets_ = false;
        bool fill_w_measures_ = false;
        bool fill_normals_ = false;

        ValueVector< ValueMap > values_;
        ValueVector< GradientFaceMap > gradients_;
        ValueVector< HessianFaceMap > hessians_;
        ValueVector< Derivatives< space_dim, face_dim,1,1 > > inv_gradients_;
        ValueVector< Derivatives< space_dim, face_dim,1,2 > > inv_hessians_;
        ValueVector< Real > dets_;
        ValueVector< Real > w_measures_;
        ValueVector< ValueMap > normals_;

        Quadrature<dim> quad_;
    };

    ElementValuesCache elem_values_;

    std::array<FaceValuesCache, n_faces> face_values_;

    Mapping<dim, codim> *mapping_;

    /**
     * @todo implement this function
     *TODO
     */
    std::array< ValueVector<ValueMap>, codim> transform_external_normals() const;
};



IGA_NAMESPACE_CLOSE

#endif // MAPPING_ELEMENT_ACCESSOR_H_
