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


#ifndef BSPLINE_ELEMENT_ACCESSOR_H_
#define BSPLINE_ELEMENT_ACCESSOR_H_

#include <igatools/base/config.h>
#include <igatools/basis_functions/space_element_accessor.h>
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
class BSplineElementAccessor :
    public SpaceElementAccessor<BSplineSpace<dim,range,rank>>
{
public:
    using parent_t = SpaceElementAccessor<BSplineSpace<dim,range,rank>>;

    /** Type for the grid accessor. */
    using GridAccessor = CartesianGridElementAccessor<dim>;

    /** Type required by the GridForwardIterator templated iterator */
    using ContainerType = const BSplineSpace<dim, range, rank> ;

    /** Type required for the generic algorithm on the spaces (plots??) */
    using Space = BSplineSpace<dim, range, rank> ;

    /** Number of faces of the element. */
    using parent_t::n_faces;

    using ValuesCache = typename parent_t::ValuesCache;

    using parent_t::admisible_flag;


public:
    template <int order>
    using Derivative = typename parent_t::template Derivative<order>;
    using typename parent_t::Point;
    using typename parent_t::Value;
    using typename parent_t::RefPoint;
    //using typename parent_t::Div;


public:
    /** @name Constructors */
    ///@{
    /**
     * Default constructor. Not allowed to be used.
     */
    BSplineElementAccessor() = delete;

    /**
     * Constructs an accessor to element number index of a
     * BsplineSpace space.
     */
    BSplineElementAccessor(const std::shared_ptr<ContainerType> space,
                           const int elem_index);

    /**
     * Copy constructor.
     * @note For the constructed object it
     * creates a new element cache, but it shares
     * the one dimensional cache with the copied element.
     */
    BSplineElementAccessor(const BSplineElementAccessor<dim, range, rank> &elem)
        = default;

    /**
     * Move constructor.
     */
    BSplineElementAccessor(BSplineElementAccessor<dim, range, rank> &&elem)
        = default;

    /**
     * Destructor.
     */
    ~BSplineElementAccessor() = default;
    ///@}


    /** @name Assignment operators */
    ///@{

    /**
     * Copy assignment operator.
     * @note Creates a new element cache, but it shares
     * the one dimensional cache with the copied element.
     */
    BSplineElementAccessor<dim, range, rank> &
    operator=(const BSplineElementAccessor<dim, range, rank> &elem)
        = default;

    /**
     * Move assignment operator.
     */
    BSplineElementAccessor<dim, range, rank> &
    operator=(BSplineElementAccessor<dim, range, rank> &&elem)
        = default;
    ///@}



    /** @name Cache initialization and filling */
    ///@{

    /**
     * Initializes the internal cache for the efficient
     * computation of the values requested in
     * the fill_flag on the given quadrature points.
     * This implies a uniform quadrature scheme
     * (i.e. the same for all elements).
     * @note This function should be called before fill_cache()
     */
    void init_cache(const ValueFlags fill_flag,
                    const Quadrature<dim> &quad);

    /**
     * For a given face quadrature.
     */
    void init_face_cache(const Index face_id,
                         const ValueFlags fill_flag,
                         const Quadrature<dim-1> &quad);

    /**
     * Fills the element values cache according to the evaluation points
     * and fill flags specifies in init_cache.
     *
     * @note The topology for which the measure is computed is specified by
     * the input argument @p topology_id.
     */
    void fill_cache(const TopologyId<dim> &topology_id = ElemTopology<dim>());
    ///@}


public:

    /** @name Functions for the basis and field evaluations without the use of the cache */
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
    evaluate_basis_derivatives_at_points(const std::vector<RefPoint> &points) const;

    ///@}

    /**
     * Prints internal information about the BSplineElementAccessor.
     * Its main use is for testing and debugging.
     */
    void print_info(LogStream &out, const VerbosityLevel verbosity_level = VerbosityLevel::normal) const;

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
    using BasisValues1d = std::vector<DenseMatrix>;


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

    ComponentContainer<DynamicMultiArray<std::shared_ptr<BSplineElementScalarEvaluator<dim>>,dim>> scalar_evaluators_;


    using univariate_values_t = ComponentContainer <std::array<const BasisValues1d *,dim>>;

    /**
     * Fills the cache (accordingly with the flags_handler status)
     * from the univariate values (and derivatives up to the order
     * specified by @p max_deriv_order).
     *
     *
     * @note The BSplineElementAccessor @p elem is needed in order to call the function
     * elem.evaluate_bspline_derivatives<p>()
     * that computes the @p p-th order derivatives of a BSpline from the univariate values.
     */
    void fill_values_cache_from_univariate(const int max_deriv_order,
                                           const univariate_values_t &values_1D,
                                           ValuesCache &cache);




    ///@}



    /**
     * Computes the k-th order derivative of the non-zero B-spline basis
     * functions over the current element,
     *   at the evaluation points pre-allocated in the cache.
     *
     * \warning If the output result @p derivatives_phi_hat is not correctly pre-allocated,
     * an exception will be raised.
     */
    template <int deriv_order>
    void evaluate_bspline_derivatives(const ComponentContainer<std::array<const BasisValues1d *, dim> > &elem_values,
                                      const ValuesCache &cache,
                                      ValueTable<
                                      Conditional<(deriv_order==0),Value,Derivative<deriv_order> >
                                      > &derivatives_phi_hat) const;





    class GlobalCache : public CacheStatus
    {
    public:
        /**
         * univariate B-splines values and derivatives at
         * quadrature points
         * splines1d_cache_data_[comp][dir][interval][order][function][point]
         */
        ComponentDirectionTable<BasisValues1d> splines1d_cache_data_;

        // ComponentDirectionTable<const BasisValues1d *> splines1d_cache_;

        int max_deriv_order_ = 0;

        // TODO (pauletti, May 30, 2014): the type for intervals_id should be CartesianProductArray
    protected:
        void reset(const Space &space,
                   const Quadrature<dim> &quad,
                   const int max_der,
                   const std::array<std::vector<int>,dim> &intervals_id);

        TensorSize<dim> n_intervals_;
    };

    /**
     * Cache for the efficient use of Bspline basis on a uniform
     * quadrature scheme on all elements.
     */
    class GlobalElemCache : public GlobalCache
    {
    public:
        /**
         * Allocates space for and compute the values and
         * derivatives at quadrature points for one
         * dimensional B-splines for each component
         * and each direction of a component and each interval
         * of a direction.
         */
        void reset(const Space &space,
                   const Quadrature<dim> &quad,
                   const int max_der);
    };

    class GlobalFaceCache : public GlobalCache
    {
    public:
        /**
         * Allocates space for and compute the values and
         * derivatives at quadrature points for one
         * dimensional B-splines for each component
         * and each direction of a component and each interval
         * of a direction.
         */
        void reset(const Space &space,
                   const Quadrature<dim> &quad1,
                   const Index face_id,
                   const int max_der);
    };



    /** Reset the global cache */
    void reset_global_cache();

    /**
     * Initilizes (reserves memory) the
     * univariate basis values
     * and derivatives at quadrature points cache
     * for efficient use of computations.
     * This function implies a uniform quadrature schemes
     * (the same for each element).
     * The fill_flag provides what information to compute.
     */
    void reset_univariate_cache(const Quadrature<dim> &quad,
                                const int max_der);


    /**
     * Tensor product style space sized cache for
     * storing the values and derivatives of the
     * one dimensional values of the basis
     * function at the quadrature points
     */
    std::shared_ptr< GlobalElemCache > values_1d_elem_ = nullptr;

    std::array<std::shared_ptr<GlobalFaceCache>, n_faces> values_1d_faces_;



protected:

#if 0
    /** Returns the Bezier extraction operator relative to the current element. */
    ComponentTable< std::array< const DenseMatrix *,dim> > get_bezier_extraction_operator() const;
#endif

private:


    template <typename Accessor> friend class GridForwardIterator;




protected:
    const ComponentContainer<DynamicMultiArray<std::shared_ptr<BSplineElementScalarEvaluator<dim>>,dim> >
            &get_scalar_evaluators() const;

};

IGA_NAMESPACE_CLOSE

#endif
