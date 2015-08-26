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

#ifndef BSPLINE_ELEMENT_H_
#define BSPLINE_ELEMENT_H_

#include <igatools/base/config.h>
#include <igatools/basis_functions/reference_element.h>

#include <igatools/linear_algebra/dense_matrix.h>
#include <igatools/basis_functions/bernstein_basis.h>
#include <igatools/basis_functions/bspline_element_scalar_evaluator.h>

IGA_NAMESPACE_OPEN

template <int dim, int range, int rank> class BSplineSpace;
template <int dim, int range, int rank> class BSplineElementHandler;
template <class Accessor> class GridIterator;

/**
 * See module on \ref accessors_iterators for a general overview.
 * @ingroup elements
 * @ingroup serializable
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
    using GridAccessor = GridElement<dim>;

    /** Type required by the CartesianGridIterator templated iterator */
    using ContainerType = const BSplineSpace<dim, range, rank> ;

    /** Type required for the generic algorithm on the spaces (plots??) */
    using Space = BSplineSpace<dim, range, rank> ;

//    using ValuesCache = typename parent_t::ValuesCache;

    using Grid = CartesianGrid<dim>;
    using IndexType = typename Grid::IndexType;
    using List = typename Grid::List;
    using ListIt = typename Grid::ListIt;

public:
    template <int order>
    using Derivative = typename parent_t::template Derivative<order>;
    using typename parent_t::Point;
    using typename parent_t::Value;

    /** @name Constructors */
    ///@{
public:
    /**
     * Default constructor. It does nothing but it is needed for the
     * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     * mechanism.
     */
    BSplineElement() = default;



    /**
     * Constructs an accessor to element number index of a
     * BsplineSpace space.
     */
    BSplineElement(const std::shared_ptr<ContainerType> space,
                   const ListIt &index,
                   const PropId &prop = ElementProperties::active);

    /**
     * Copy constructor.
     * It can be used with different copy policies (i.e. deep copy or shallow
     *  copy).
     * The default behaviour (i.e. using the proper interface of a classic copy
     * constructor) uses the deep copy.
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
    virtual ~BSplineElement() = default;
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



private:
    /**
     * @name Containers for the cache of the element values and for the
     * cache of the face values
     */
    ///@{

    /**
     * This type store the values, first second derivatives
     * of a 1D Bspline functions, i.e BasisValues1d[k]
     * stores the values of the k-th derivative of the (p+1) basis function
     * on a given interval
     * at the quadrature points.
     * BasisValues1d[k] is a (p+1) x n_qp matrix
     */
    using BasisValues1d = SafeSTLVector<DenseMatrix>;

protected:
    /**
     * For each component gives a product array of the dimension
     */
    template<class T>
    using ComponentContainer =
        typename Space::BaseSpace::template ComponentContainer<T>;

    /**
     * For each component gives a product array of the dimension
     */
    template<class T>
    using ComponentDirectionTable =
        ComponentContainer<CartesianProductArray<T,dim>>;


private:
    /**
     * Returns the BSplineSpace upon which the element is defined.
     */
    std::shared_ptr<const Space> get_bspline_space() const;


    template <class Accessor> friend class GridIterator;
    friend class BSplineElementHandler<dim, range, rank>;

public:
    ComponentContainer<SafeSTLArray<ValueTable<Real>,dim> >
    evaluate_univariate_derivatives_at_points(const int deriv_order,
                                              const Quadrature<dim> &pts) const
    {
        ComponentContainer<SafeSTLArray<ValueTable<Real>,dim> > values;

        Assert(false,ExcNotImplemented());

        return values;
    }

    virtual std::shared_ptr<SpaceElement<dim,0,range,rank,Transformation::h_grad> >
    clone() const override final;


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
    serialize(Archive &ar, const unsigned int version);
    ///@}
#endif // SERIALIZATION
};

IGA_NAMESPACE_CLOSE

#endif
