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

#ifndef CARTESIAN_GRID_ELEMENT_ACCESSORS_H_
#define CARTESIAN_GRID_ELEMENT_ACCESSORS_H_

#include <igatools/base/config.h>
#include <igatools/base/cache_status.h>
#include <igatools/base/value_flags_handler.h>
#include <igatools/base/quadrature.h>
#include <igatools/geometry/topology.h>
#include <igatools/geometry/cartesian_grid_element.h>
#include <igatools/geometry/grid_forward_iterator.h>
#include <igatools/utils/value_vector.h>


IGA_NAMESPACE_OPEN



/**
 * @brief Element accessor for the CartesianGrid.
 *
 * The element can be queried for informations
 * that can be generated on-the-fly
 * (i.e. without the use of a cache) and for informations
 * that are obtained throuch a cache mechanism
 *
 *
 * See module on \ref accessors_iterators for a general overview.
 * @ingroup accessors_iterators
 *
 * @author S.Pauletti, 2012, 2013, 2014
 * @author M.Martinelli, 2013, 2014
 */
template <int dim_>
class CartesianGridElementAccessor : public CartesianGridElement<dim_>
{
public:
    /** Type required by the GridForwardIterator templated iterator */
    using ContainerType = CartesianGrid<dim_>;

    /** Dimension of the grid like container */
    static const auto dim = ContainerType::dim;


    /** Fill flags supported by this iterator */
    static const ValueFlags admisible_flag =
        ValueFlags::point|
        ValueFlags::measure |
        ValueFlags::w_measure |
        ValueFlags::face_point |
        ValueFlags::face_measure |
        ValueFlags::face_w_measure |
        ValueFlags::face_normal;

public:
    /** @name Constructors */
    ///@{
    /**
     * Default constructor. Not allowed to be used.
     */
    CartesianGridElementAccessor() = delete;

    /**
     * Construct an accessor pointing to the element with
     * flat index @p elem_index of the CartesianGrid @p grid.
     */
    CartesianGridElementAccessor(const CartesianGrid<dim_> &grid,
                                 const Index elem_index);

    /**
     * Copy constructor.
     * @note For the constructed object it
     * creates a new element cache, but it shares
     * the one dimensional cache with the copied element.
     */
    CartesianGridElementAccessor(const CartesianGridElementAccessor<dim_> &elem)
        = default;

    /**
     * Move constructor.
     */
    CartesianGridElementAccessor(CartesianGridElementAccessor<dim_> &&elem)
        = default;

    /**
     * Destructor.
     */
    ~CartesianGridElementAccessor() = default;
    ///@}

    /** @name Assignment operators */
    ///@{
    /**
     * Copy assignment operator.
     * Creates a new element cache, but it shares
     * the one dimensional length cache with the copied element.
     */
    CartesianGridElementAccessor<dim_>
    &operator=(const CartesianGridElementAccessor<dim_> &elem) = default;

    /**
     * Move assignment operator.
     */
    CartesianGridElementAccessor<dim_>
    &operator=(CartesianGridElementAccessor<dim_> &&elem) = default;
    ///@}

    /** @name Functions for the cache initialization and filling. */
    ///@{
    /**
     * Initializes the internal cache for the efficient
     * computation of the values requested in
     * the @p fill_flag at the given quadrature points.
     *
     * For the face values, it allows the reuse of the element
     * cache, i.e. it is like using a projected quadrature on
     * the faces.
     */
    void init_values(const ValueFlags flag,
                     const Quadrature<dim_> &quad);

    /**
     * Initializes the internal cache for the efficient
     * computation of the values requested in
     * the @p fill_flag when no quadrature point is necessary
     */
    void init_values(const ValueFlags flag);

    /**
     * To use a different quadrature on the face instead of
     * the projected element quadrature
     */
    void init_face_values(const Index face_id,
                          const ValueFlags flag,
                          const Quadrature<dim_-1> &quad);

    /**
     * Fills the element values cache according to the evaluation points
     * and fill flags specifies in init_values.
     */
    void fill_values();

    /**
     * Fills the i-th face values cache according to the evaluation points
     * and fill flags specified in init_values.
     */
    void fill_face_values(const Index face_id);
    ///@}


    ///@name Query information that requires the use of the cache
    ///@{

    /**
     * Returns the measure of the element in the CartesianGrid<dim_>
     * object referred by this accessor.
     * This is the length for dim_==1,
     * the area for dim_==2 and the volume for dim_==3.
     */
    Real get_measure(const TopologyId &topology_id = ElemTopology()) const;


    /**
     * Returns the element measure multiplied by the weights of the quadrature scheme
     * used to initialize the accessor's cache.
     */
    ValueVector<Real> const &get_w_measures(const TopologyId &topology_id = ElemTopology()) const;


    /**
     * Return a const reference to the one-dimensional container with the
     * values of the map at the evaluation points.
     */
    std::vector<Point<dim>> const get_points(const TopologyId &topology_id = ElemTopology()) const;


    /**
     * Returns the lengths of the coordinate sides of the cartesian element.
     * For example in 2 dimensions
     * \code
     * length = coordinate_lenths();
     * length[0] is the length of the x-side of the element and
     * length[1] the length of the y-side of the element.
     * \endcode
     */
    std::array<Real, dim_> get_coordinate_lengths() const;

    ///@}





    static const Size n_faces = UnitElement<dim_>::faces_per_element;

private:

    /**
     * helper function to compute the measure of element
     */
    Real measure() const;

    /**
     * helper function to compute the measure of j-th face
     */
    Real face_measure(const int j) const;

protected:
    bool operator==(const CartesianGridElementAccessor<dim_> &a) const;

    bool operator!=(const CartesianGridElementAccessor<dim_> &a) const;

    void operator++();

private:


    /**
     * @brief Global CartesianGrid cache, storing the interval length in each direction.
     *
     * For now only a uniform quad is taken care of.
     */
    class LengthCache : public CacheStatus
    {
    public:
        /**
         * Allocates space for the cache
         */
        void reset(const CartesianGrid<dim_> &grid);

        /** pointer to the current entry of of length,
         *  it could be used for optimization of uniform grid
         */
        CartesianProductArray<Real *, dim_> length_;

        /** stores the interval length */
        CartesianProductArray<Real , dim_> length_data_;

    };


    /**
     * @brief Base class for cache of CartesianGridElementAccessor.
     */
    class ValuesCache : public CacheStatus
    {
    public:
        /**
         * Allocate space for the values at quadrature points
         */
        void reset(const GridElemValueFlagsHandler &flags_handler,const Quadrature<dim_> &quad);

        /**
         * Fill the cache member.
         * @note The @p measure is an input argument because of the different function calls
         * between element-measure and face-measure.
         */
        void fill(const Real measure);


        GridElemValueFlagsHandler flags_handler_;

        ///@name The "cache" properly speaking
        ///@{
        /** Measure of the element in the grid (equal to the product of the element lenghts). */
        Real measure_ = 0.0;

        /** Element measure multiplied by the quadrature weights. */
        ValueVector<Real> w_measure_;

        TensorProductArray<dim_> unit_points_;

        ValueVector<Real> unit_weights_;
        ///@}

    };

    /**
     * @brief Cache for the element values at quadrature points
     */
    class ElementValuesCache : public ValuesCache
    {
    public:
        /**
         * Allocate space for the values at quadrature points
         */
        void reset(const GridElemValueFlagsHandler &flags_handler,const Quadrature<dim_> &quad);
    };


    /**
     * @brief Cache for the face values at quadrature points
     */
    class FaceValuesCache : public ValuesCache
    {
    public:
        void reset(const GridFaceValueFlagsHandler &flags_handler,
                   const Quadrature<dim_> &quad,
                   const Index face_id);

        void reset(const GridFaceValueFlagsHandler &flags_handler,
                   const Quadrature<dim_-1> &quad,
                   const Index face_id);
    };

    /**
     * @todo Document this function
     */
    const ValuesCache &get_values_cache(const TopologyId &topology_id) const;


    /** Grid (global) lengths cache */
    LengthCache length_cache_;

    /** Element values cache */
    ElementValuesCache elem_values_;

    /** Face values cache */
    std::array<FaceValuesCache, n_faces> face_values_;

private:
    template <typename Accessor> friend class GridForwardIterator;

protected:
    /**
     * ExceptionUnsupported Value Flag.
     */
    DeclException2(ExcFillFlagNotSupported, ValueFlags, ValueFlags,
                   << "The passed ValueFlag " << arg2
                   << " contains a non admissible flag " << (arg1 ^arg2));
};

IGA_NAMESPACE_CLOSE

#endif /* __CARTESIAN_GRID_ELEMENT_ACCESSORS_H_ */
