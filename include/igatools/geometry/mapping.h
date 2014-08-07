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

#ifndef MAPPING_H_
#define MAPPING_H_

#include <igatools/base/config.h>
#include <igatools/base/quadrature.h>
#include <igatools/base/function.h>
#include <igatools/geometry/grid_wrapper.h>
#include <igatools/geometry/cartesian_grid.h>

IGA_NAMESPACE_OPEN

//Forward declaration to avoid including header file.
template <int> class CartesianGridElementAccessor;

//Forward declaration to avoid including header file.
template <int, int> class MappingElementAccessor;

/**
 * @brief The mapping is a deformation \f$ F : \hat\Omega \to \Omega\f$
 * which maps the reference domain \f$\hat\Omega \in \mathbb{R}^{dim}\f$ to the
 * physical domain \f$\Omega \in \mathbb{R}^{dim+codim}\f$.
 *
 * In igatools the mapping is a grid-like container, so that on top
 * of being a deformation is aware of an underlying grid structure.
 *
 * Most of the methods in this class are virtual,
 * because we want to use Mapping as general interface for user-defined Mapping
 * specializations (unknown at compile time).
 *
 * @ingroup containers
 *
 * @author M.S. Pauletti 2012, 2013, 2014
 * @author M. Martinelli 2013, 2014
 * @author P. Antolin 2014
 * @author N. Cavallini  2012
 *
 */
template<int dim_, int codim_ = 0>
class Mapping
    :  public std::enable_shared_from_this<Mapping<dim_, codim_> >,
       public GridWrapper<CartesianGrid<dim_>>
{
public:
    /** Type of the Grid */
	//TODO(pauletti, Aug 4, 2014): for some reason the current compiler doesn't
	// understand something like GridIterator if the
	// following commneted using is used
    //using typename GridWrapper<CartesianGrid<dim_>>::GridType;
	using GridType = CartesianGrid<dim_>;

	using GridIterator = typename GridType::ElementAccessor;

	/** Dimension of the reference domain */
    static const int dim = dim_;

    /** Codimension of the deformed domain. */
    static const int codim = codim_;

    /** Dimension of the deformed domain embedding space. */
    static const int space_dim = dim + codim;

private:
    using self_t = Mapping<dim, codim>;

    /** Function type of the mapping as a Function */
    using Func = Function<dim, space_dim>;

private:
    /** Type for the given order derivatives of the
     *  the mapping. */
    template<int order>
    using Derivative = typename Func::template Derivative<order>;

    /** Type for the diferent order derivatives of the inverse of
     * the mapping
     */
    template<int order>
    using InvDerivative = Derivatives<space_dim, dim, 1, order>;

public:
    /** Type of the mapping evaluation point. */
    using Point = typename Func::Point;

    /** Type of the mapping return value. */
    using Value = typename Func::Value;

    /** Type of the mapping gradient. */
    using Gradient = typename Func::Gradient;

    /** Typedef for the mapping hessian. */
    using Hessian = typename Func::Hessian;

public:
    using FaceMapping = Conditional<(dim>0), Mapping<dim-1, codim+1>, self_t >;

    /** Type of the element accessor */
    using ElementAccessor = MappingElementAccessor<dim, codim>;

    /** Type of the element iterator */
    using ElementIterator = GridForwardIterator<ElementAccessor>;

public:
    /** @name Constructors and destructor */
    ///@{
    /** Default constructor.*/
    Mapping() = delete;

    /** Constructs map over grid. */
    Mapping(const std::shared_ptr<GridType> grid);

    /** Destructor */
    virtual ~Mapping();

    /**
     * Copy constructor. The new object has a deep copy (i.e. a new instance)
     * of the grid held by the copied object @p map.
     */
    Mapping(const Mapping<dim_,codim_> &map);
    ///@}

    /** @name Assignment operators. */
    ///@{

    /** Copy assignment operator. Not allowed to be used. */
    Mapping<dim_,codim_> &operator=(const Mapping<dim_,codim_> &map) = delete;
    ///@}
//TODO(pauletti, Aug 6, 2014): use template<order> evaluate_derivative()

    /** @name Mapping as a standard function (using the cache).*/
    ///@{
    virtual void evaluate(std::vector<Value> &values) const;

    virtual void evaluate_gradients(std::vector<Gradient> &gradients) const;

    virtual void evaluate_hessians(std::vector<Hessian> &hessians) const;

    virtual void evaluate_face(const Index face_id,
                               std::vector<Value> &values) const;

    virtual void evaluate_face_gradients(const Index face_id,
                                         std::vector<Gradient> &gradients) const;

    virtual void evaluate_face_hessians(const Index face_id,
                                        std::vector<Hessian> &hessians) const;
    ///@}


    /** @name Mapping as a standard function (without the use of the cache).*/
    ///@{
    virtual void evaluate_at_points(const std::vector<Point> &points,
                                    std::vector<Value> &values) const ;

    virtual void evaluate_gradients_at_points(const std::vector<Point> &points,
                                              std::vector<Gradient> &gradients) const;

    virtual void evaluate_hessians_at_points(const std::vector<Point> &points,
                                             std::vector<Hessian> &hessians) const;
    ///@}


    /** @name Virtual user functions to define the map */
    ///@{
    /**
     * An element based mapping may require some initialization.
     *
     * @warning This function must be reimplemented by in every concrete child
     * class of Mapping.
     */
    virtual void init_element(const ValueFlags flag,
                              const Quadrature<dim> &quad) const = 0;

    virtual void set_element(const GridIterator &elem) const = 0;

    /**
     *
     * @todo evaluate if index should be flat index, tensor index, and
     * or GridElement iterator
     */
    virtual void set_face_element(const Index face_id,
                                  const GridIterator &elem) const = 0;


    virtual std::vector<Value> values() const;

    virtual std::vector<Gradient> gradients() const;

    virtual std::vector<Hessian> hessians() const;
    ///@}

    /** @name Dealing with the element-based iterator. */
    ///@{
    /**
     * Returns a element iterator to the first element of the patch.
     */
    virtual ElementIterator begin() const;

    /**
     * Returns a element iterator to the last element of the patch.
     */
    virtual ElementIterator last() const;

    /**
     * Returns a element iterator to one-pass the end of patch.
     */
    virtual ElementIterator end() const;
    ///@}

    /**
     * Prints internal information about the mapping.
     * @note Mostly used for debugging and testing.
     * @warning Calling this function will throw an assertion.
     * Try to call the same function on a derived class.
     */
    virtual void print_info(LogStream &out) const;

private:

    //TODO(pauletti, Jun 21, 2014): the reason of this functions is unclear
    /**
     * Returns the flag required to evaluate this mapping.
     * This is used from the mapping accessor.
     */
    virtual ValueFlags required_flags() const;

    friend ElementAccessor;
};

IGA_NAMESPACE_CLOSE

#endif /* MAPPING_H_ */
