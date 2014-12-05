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


#ifndef BSPLINE_ELEMENT_H_
#define BSPLINE_ELEMENT_H_

#include <igatools/base/config.h>
#include <igatools/basis_functions/reference_element.h>
#include <igatools/basis_functions/bspline_element_handler.h>

#include <igatools/linear_algebra/dense_matrix.h>
#include <igatools/basis_functions/bernstein_basis.h>
#include <igatools/basis_functions/bspline_element_scalar_evaluator.h>

IGA_NAMESPACE_OPEN

template <int dim, int range, int rank> class BSplineSpace;
template <typename Accessor> class GridForwardIterator;





/**
 * See module on \ref accessors_iterators for a general overview.
 * @ingroup accessors
 */
template <int dim, int range, int rank>
class BSplineElement :
    public ReferenceElement<dim,range,rank>
{
private:
    using self_t = BSplineElement<dim,range,rank>;
    using parent_t = ReferenceElement<dim,range,rank>;

public:
    /** Type for the grid accessor. */
    using GridAccessor = CartesianGridElement<dim>;

    /** Type required by the GridForwardIterator templated iterator */
    using ContainerType = const BSplineSpace<dim, range, rank> ;

    /** Type required for the generic algorithm on the spaces (plots??) */
    using Space = BSplineSpace<dim, range, rank> ;

    using ValuesCache = typename parent_t::ValuesCache;

public:
    template <int order>
    using Derivative = typename parent_t::template Derivative<order>;
    using typename parent_t::Point;
    using typename parent_t::Value;
    //using typename parent_t::Div;


public:
    /** @name Constructors */
    ///@{
    /**
     * Default constructor
     */
    BSplineElement() = default;

    /**
     * Constructs an accessor to element number index of a
     * BsplineSpace space.
     */
    BSplineElement(const std::shared_ptr<ContainerType> space,
                   const Index elem_index);

    BSplineElement(const std::shared_ptr<ContainerType> space,
                   const TensorIndex<dim> &elem_index);
    /**
     * Copy constructor.
     * It can be used with different copy policies (i.e. deep copy or shallow copy).
     * The default behaviour (i.e. using the proper interface of a classic copy constructor)
     * uses the deep copy.
     */
    BSplineElement(const self_t &elem,
                   const CopyPolicy &copy_policy = CopyPolicy::deep);

    /**
     * Move constructor.
     */
    BSplineElement(self_t &&elem) = default;

    /**
     * Destructor.
     */
    ~BSplineElement() = default;
    ///@}

    /** @name Assignment operators */
    ///@{
    /**
     * Copy assignment operator.
     * @note Creates a new element cache, but it shares
     * the one dimensional cache with the copied element.
     */
    self_t &operator=(const self_t &elem) = default;

    /**
     * Move assignment operator.
     */
    self_t &operator=(self_t &&elem) = default;
    ///@}

public:
    /** @name Functions for the basis and field evaluations without the use of
     * the cache */
    ///@{
    /**
     * Returns a ValueTable with the <tt>deriv_order</tt>-th derivatives of all local basis function
     * at each point (in the unit domain) specified by the input argument <tt>points</tt>.
     * @note This function does not use the cache and therefore can be called any time without
     * needing to pre-call init_cache()/fill_cache().
     * @warning The evaluation <tt>points</tt> must belong to the unit hypercube
     * \f$ [0,1]^{\text{dim}} \f$ otherwise, in Debug mode, an assertion will be raised.
     */
    template <int deriv_order>
    ValueTable< Conditional< deriv_order==0,Value,Derivative<deriv_order> > >
    evaluate_basis_derivatives_at_points(const ValueVector<Point> &points) const;

    ///@}



private:
    /**
     * @name Containers for the cache of the element values and for the
     * cache of the face values
     */
    ///@{

    /**
     * This type store the values, first second derivatives
     * of a 1D Bspline functions, i.e BasisValues1d[k]
     * stores the values of the k-th derivative of the (p+1) basis function on a given interval
     * at the quadrature points.
     * BasisValues1d[k] is a (p+1) x n_qp matrix
     */
    using BasisValues1d = vector<DenseMatrix>;


protected:

    /**
     * For each component gives a product array of the dimension
     */
    template<class T>
    using ComponentContainer = typename Space::BaseSpace::template ComponentContainer<T>;

    /**
     * For each component gives a product array of the dimension
     */
    template<class T>
    using ComponentDirectionTable = ComponentContainer<CartesianProductArray<T,dim>>;


private:


    template <typename Accessor> friend class GridForwardIterator;
    friend class BSplineElementHandler<dim, range, rank>;

    std::shared_ptr<const Space> space_;

public:
    /*
        const ComponentContainer<DynamicMultiArray<std::shared_ptr<BSplineElementScalarEvaluator<dim>>,dim> >
                &get_scalar_evaluators() const;
    //*/

#if 0
    ComponentContainer<std::array<ValueTable<Real>,dim> >
    get_univariate_derivatives(const int deriv_order) const;
#endif

    /*
     * Returns a component table with the derivatives (of order @p deriv_order)
     * of the 1D basis function in each direction.
     * @warning The evaluation <tt>points</tt> must belong to the unit hypercube
     * \f$ [0,1]^{\text{dim}} \f$ otherwise, in Debug mode, an assertion will be raised.
     */
    ComponentContainer<std::array<ValueTable<Real>,dim> >
    evaluate_univariate_derivatives_at_points(const int deriv_order, const Quadrature<dim> &quad) const;

    /*
     * Returns a component table with the derivatives (of order @p deriv_order)
     * of the 1D basis function in each direction.
     * @warning The evaluation <tt>points</tt> must belong to the unit hypercube
     * \f$ [0,1]^{\text{dim}} \f$ otherwise, in Debug mode, an assertion will be raised.
     */
    ComponentContainer<std::array<ValueTable<Real>,dim> >
    evaluate_univariate_derivatives_at_points(const int deriv_order, const ValueVector<Point> &points) const;


private:
    ComponentContainer<std::array<ValueTable<Real>,dim> >
    evaluate_univariate_derivatives_at_points(
        const int deriv_order,
        const std::array<vector<Real>,dim> &points) const;
};

IGA_NAMESPACE_CLOSE

#endif
