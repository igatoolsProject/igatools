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

#ifndef CARTESIAN_GRID_ELEMENT_H_
#define CARTESIAN_GRID_ELEMENT_H_

#include <igatools/base/config.h>
#include <igatools/base/tensor.h>
#include <igatools/geometry/cartesian_grid.h>
#include <igatools/geometry/topology.h>
#include <igatools/utils/value_vector.h>

IGA_NAMESPACE_OPEN




/**
 * @brief This class represents an element within a CartesianGrid.
 *
 * The element can be queried for information
 * that is generated on-the-fly
 * (i.e. without the use of a cache).
 *
 * It is used as base class for CartesianGridElementAccessor which
 * generate the element information using caches techniques.
 *
 * @tparam dim Dimensionality of the grid.
 *
 * @author M.Martinelli, 2014
 */
template <int dim_>
class CartesianGridElement
{
public:
    /** Type required by the GridForwardIterator templated iterator */
    using ContainerType = const CartesianGrid<dim_>;

    /** Dimension of the grid like container */
    static const auto dim = ContainerType::dim;


    /** @name Constructors */
    ///@{
    /**
     * Default constructor.
     */
    CartesianGridElement() = default;

    /**
     * Construct an object pointing to the element with
     * flat index @p elem_index of the CartesianGrid @p grid.
     */
    CartesianGridElement(const std::shared_ptr<ContainerType> grid,
                         const Index elem_index);

    CartesianGridElement(const std::shared_ptr<ContainerType> grid,
                         const TensorIndex<dim> &elem_index);

    /**
     * Copy constructor.
     */
    CartesianGridElement(const CartesianGridElement<dim_> &elem)
        = default;

    /**
     * Move constructor.
     */
    CartesianGridElement(CartesianGridElement<dim_> &&elem)
        = default;

    /**
     * Destructor.
     */
    ~CartesianGridElement() = default;
    ///@}

    /** @name Assignment operators */
    ///@{
    /**
     * Copy assignment operator.
     */
    CartesianGridElement<dim>
    &operator=(const CartesianGridElement<dim_> &elem)
    {
        if ((*this) != elem)
        {
            Assert(grid_ == elem.grid_, ExcMessage("should be same mesh"));
            flat_index_ = elem.flat_index_;
            tensor_index_ = elem.tensor_index_;
        }
        return *this;
    }


    /**
     * Move assignment operator. Not allowed to be used.
     */
    CartesianGridElement<dim>
    &operator=(CartesianGridElement<dim_> &&elem) = default;
    ///@}

    /** Return the cartesian grid from which the element belongs.*/
    const std::shared_ptr<ContainerType> get_grid() const;


    /** @name Functions related to the indices of the element in the cartesian grid. */
    ///@{
    /** Returns the index of the element in its flatten representation. */
    Index get_flat_index() const;

    /** Returns the index of the element in its tensor representation. */
    TensorIndex<dim> get_tensor_index() const;
    ///@}

public:
    /** @name Query geometrical/topological information without use of cache */
    ///@{
    /**
     * Return the @p i-th vertex
     */
    Points<dim> vertex(const int i) const;

    /**
     * Return the center of the element.
     */
    Points<dim> center() const;

    /**
     * Returns the lengths of the coordinate sides of the cartesian element.
     * For example in 2 dimensions
     * \code{.cpp}
       auto length = elem.coordinate_lenths();
       // length[0] is the length of the x-side of the element and
       // length[1] the length of the y-side of the element.
       \endcode
     */
    std::array<Real,dim> get_coordinate_lengths() const;

    /**
     * Returns measure of the element or of the element-face in the
     * CartesianGrid.
     * @note The topology for which the measure is computed is specified by
     * the input argument @p topology_id.
     */
    Real get_measure(const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

    /**
     * Returns measure of j-th face.
     */
    Real get_face_measure(const int j) const;

    /**
     * Test if the point is inside the element.
     */
    bool is_point_inside(const Points<dim> &point) const;

    /**
     * Test if the point is on the element boundary.
     */
    bool is_point_on_boundary(const Points<dim> &point) const;

    /**
     * Test if the element has a boundary face.
     */
    bool is_boundary() const;

    /**
     * Test if the face identified by @p face_id on the current element is on the
     * boundary of the cartesian grid.
     */
    bool is_boundary(const Index face_id) const;
    ///@}


    /**
     * This function take as input a vector of points in the unitary hypercube [0,1]^{dim}
     * and returns the points mapped over the domain (in the parametric coordinate system)
     * represented by this GridElementAccessor.
     */
    ValueVector<Points<dim> >
    transform_points_unit_to_reference(const ValueVector<Points<dim>> &point_unit_domain) const;

    /**
     * This function takes as input argument a vector of points over the element
     * represented by this GridElementAccessor
     * (i.e. the point coordinates are in the parametric coordinate system)
     * and returns the points mapped over the
     * points unitary hypercube [0,1]^{dim}.
     */
    ValueVector<Points<dim> >
    transform_points_reference_to_unit(const ValueVector<Points<dim>> &point_reference_domain) const;


    /**
     * Prints internal information about the CartesianGridElement.
     * Its main use is for testing and debugging.
     */
    void print_info(LogStream &out,
                    const VerbosityLevel verbosity = VerbosityLevel::normal) const;


    bool is_influence() const;
    void set_influence(const bool influence_flag);

    bool is_active() const;
    void set_active(const bool active_flag);


    /**
     * True if the element index valid (i.e. inside the grid from which the element belongs from).
     *
     * @note A typical case in which this function returns false is when the element is moved outside
     * the grid from which belongs from.
     */
    bool is_valid() const;

    /** @name Functions/operators for moving the element in the CartesianGrid.*/
    ///@{
    /**
     * Moves the element to the position that differs from the current one
     * for the quantity given by @p increment.
     *
     * If the resulting position after the movement is valid (i.e. within the grid), then the function
     * returns true, otherwise it returns false.
     */
    bool jump(const TensorIndex<dim> &increment);

    /**
     * Sets the index of the element using the flatten representation.
     * @note This function also updates the index for the tensor representation.
     * @warning This may be a dangerous function, be careful when using it
     * as it is easy to use incorrectly. Only use it if you know what you
     * are doing.
     */
    void move_to(const Index flat_index);


    /**
     * Sets the index of the element using the tensor representation.
     * @note This function also updates the index for the flatten representation.
     * @warning this may be a dangerous function, be careful when using it
     * as it is easy to use incorrectly. Only use it if you know what you
     * are doing.
     */
    void move_to(const TensorIndex<dim> &tensor_index);

    // TODO (pauletti, Aug 21, 2014): the next operators should be protected
    // someone made them public due to hackish code in NURBSelementaccessor
    // we must rethink that code

    /** Moves the element to the next valid element in the CartesianGrid. */
    void operator++();
    ///@}

    /** @name Comparison operators*/
    ///@{
    bool operator==(const CartesianGridElement<dim_> &elem) const;

    /**
     * Returns true if the the CartesianGridElement @p elem has a different flat index w.r.t.
     * the calling object.
     * @note The calling object and the CartesianGridElement @p elem must refers to the same
     * CartesianGrid, otherwise an exception will be raised (in Debug mode).
     */
    bool operator!=(const CartesianGridElement<dim_> &elem) const;
    ///@}

protected:
    /** Cartesian grid from which the element belongs.*/
    std::shared_ptr<ContainerType> grid_;

    /** Flat (linear) index assigned to the current (sub)-element. */
    Index flat_index_;

    /** Tensor product indices of the current struct index @p flat_index_. */
    TensorIndex<dim> tensor_index_;
};

IGA_NAMESPACE_CLOSE

#endif /* CARTESIAN_GRID_ELEMENT_H_ */
