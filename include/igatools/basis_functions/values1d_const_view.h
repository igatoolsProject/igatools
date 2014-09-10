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


#ifndef VALUES_1D_CONST_VIEW_H_
#define VALUES_1D_CONST_VIEW_H_


#include <igatools/base/config.h>
#include <igatools/linear_algebra/dense_matrix.h>
#include <igatools/utils/vector.h>
#include <igatools/utils/tensor_index.h>

IGA_NAMESPACE_OPEN

/**
 * Set of functions evaluation (values and derivatives) which
 * are defined as tensor product of scalar functions over
 * points also defined as tensor products.
 */
class BasisValues1d
{
public:
    BasisValues1d()
{}
    BasisValues1d(const int max_der_order, const int n_func, const int n_points)
    :
        values_(max_der_order, DenseMatrix(n_func, n_points))
    {}

    Size get_num_points() const
    {
    	return values_[0].size2();
    }
    Size get_num_functions() const
    {
    	return values_[0].size1();
    }
    void resize(const int max_der_order, const int n_func, const int n_points)
    {
        values_.resize(max_der_order);
        for (auto matrix: values_)
            matrix.resize(n_func, n_points);
    }

    void print_info(LogStream &out) const
    {
        values_.print_info(out);
    }

    auto &get_derivative(const int order)
    {
        return values_[order];
    }

    auto const &get_derivative(const int order) const
    {
        return values_[order];
    }

private:
    vector<DenseMatrix> values_;
};



/**
 * Reference (View) of a BasisValues1d
 */
class BasisValues1dConstView
{
public:
    /** @name Constructors */
    ///@{
    /** Default constructor. It does nothing. */
    BasisValues1dConstView() = default;

    /**
     * Constructor. Builds the const view on the <tt>func_id</tt>-th row
     * of the DenseMatrix @p funcs.
     */
    BasisValues1dConstView(const BasisValues1d &val)
    :funcs_(&val)
    {}

    auto const &get_derivative(const int order) const
    {
        return funcs_->get_derivative(order);
    }

    /** Copy constructor. */
    BasisValues1dConstView(const BasisValues1dConstView &view) = default ;

    /** Move constructor. */
    BasisValues1dConstView(BasisValues1dConstView &&view) = default ;

    /** Destructor. */
    ~BasisValues1dConstView() = default;
    ///@}

    /** Assignment operators */
    ///@{
    /** Copy assignment operator. */
    BasisValues1dConstView &operator=(const BasisValues1dConstView &view) = default;

    /** Move assignment operator. */
    BasisValues1dConstView &operator=(BasisValues1dConstView &&view) = default;
    ///@}
//
//    /** Returns the value of the fucntion at the <tt>point_id</tt>-th point. */
//    Real operator()(const Index point_id) const;
//
//    /** Return the number of points for which the function is evaluated. */
//    Size get_num_points() const;

private:
    BasisValues1d const* funcs_;
};


template <int dim>
using ElemFuncValues = std::array<BasisValues1dConstView, dim>;

template <int dim>
class TensorProductFunctionEvaluator :
		public ElemFuncValues<dim>
{
public:
	TensorProductFunctionEvaluator(TensorSize<dim> n_func, TensorSize<dim> n_pts)
:f_size(n_func),
 p_size(n_pts)
{}
	/**
     * Evaluate and returns one partial derivative in one point.
     * The order of the partial derivative is specified by the tensor-index
     * @p order_tensor_id,
     * while the point is specified by its tensor-index @p point_tensor_id.
     */
    Real evaluate( const TensorIndex<dim> &order,
                   const TensorIndex<dim> &func,
                   const TensorIndex<dim> &pt) const
    {
        Real res = dim>0 ? (*this)[0].get_derivatives(order[0])[func[0]][pt[0]] : 1.;
        for (int i = 1; i < dim; ++i)
            res *= (*this)[i].get_derivatives(order[i])[func[i]][pt[i]];
        return res;
    }
private:
    TensorSizedContainer<dim> f_size;
    TensorSizedContainer<dim> p_size;
};



/**
 * @brief Const view to one-dimensional BSpline function over an interval.
 *
 * In the igatools library, the one-dimensional BSpline values and derivatives
 * are computed at the evaluation points over each grid interval and
 * then stored using a DenseMatrix object,
 * where:
 * - the rows refers to the non-zero BSpline functions over the interval;
 * - the columns refers to the evaluation points over the interval.
 * This means that if we want to refer to the value of the
 * function @p fn at the point @p pt
 * using a DenseMatrix we must use two indices, e.g.
 * @code{.cpp}
   DenseMatrix values; // this represents all BSpline functions over an interval
   Real value_fn_pt = values(fn,pt); // value of the function fn at point pt
   @endcode
 *
 * This can be inefficient if we need to access the values on the point of a
 * specific function more then one time
 * (e.g. in a loop). Then, the purpose of this class is to represent a single
 * BSpline over an interval providing
 * the access operator operator(const Index point) that can be used to get the
 * values of the function at the different points.
 * The unique constructor
 * Values1DConstView(const DenseMatrix &funcs,const Index func_id) takes the
 * DenseMatrix holding the values and the index of the function we want to
 * represent and than the Values1DConstView object is built using just the
 * memory address of the DenseMatrix plus the fucntion index, resulting
 * in a lightweight object with minimal memory footprint and no expensive
 * copy of function values.
 *
 * @todo Document more
 * @author M. Martinelli
 * @date 29 Mar 2014
 */
class Values1DConstView
{
public:
    /** Type for the container of one dimensional values on a single interval for a single scalar function.*/
    using Values1D = typename DenseMatrix::MatrixRowType ;

    using const_iterator = typename Values1D::const_iterator;

    /** @name Constructors */
    ///@{
    /** Default constructor. It does nothing. */
    Values1DConstView() = default;

    /**
     * Constructor. Builds the const view on the <tt>func_id</tt>-th row
     * of the DenseMatrix @p funcs.
     */
    Values1DConstView(const DenseMatrix &funcs,const Index func_id);

    /** Copy constructor. */
    Values1DConstView(const Values1DConstView &view) = default ;

    /** Move constructor. */
    Values1DConstView(Values1DConstView &&view) = default ;

    /** Destructor. */
    ~Values1DConstView() = default;
    ///@}

    /** Assignment operators */
    ///@{
    /** Copy assignment operator. */
    Values1DConstView &operator=(const Values1DConstView &view) = default;

    /** Move assignment operator. */
    Values1DConstView &operator=(Values1DConstView &&view) = default;
    ///@}

    /** Returns the value of the fucntion at the <tt>point_id</tt>-th point. */
    Real operator()(const Index point_id) const;

    /** Return the number of points for which the function is evaluated. */
    Size get_num_points() const;

private:
    const DenseMatrix *funcs_ = nullptr;
    Index func_id_ = 0;
};

IGA_NAMESPACE_CLOSE


#endif // #ifndef VALUES_1D_CONST_VIEW_H_

