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


#ifndef VALUES_1D_CONST_VIEW_H_
#define VALUES_1D_CONST_VIEW_H_


#include <igatools/base/config.h>
#include <igatools/linear_algebra/dense_matrix.h>
#include <igatools/utils/vector.h>
#include <igatools/utils/tensor_index.h>
#include <igatools/utils/tensor_sized_container.h>

IGA_NAMESPACE_OPEN

/**
 * Container for scalar function values and derivatives
 * computed  over points in an interval.
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

    auto &get_dataivative(const int order)
    {
        return values_[order];
    }

    auto const &get_dataivative(const int order) const
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

    /** Copy constructor. */
    BasisValues1dConstView(const BasisValues1dConstView &view) = default ;

    /** Move constructor. */
    BasisValues1dConstView(BasisValues1dConstView &&view) = default ;

    /** Destructor. */
    ~BasisValues1dConstView() = default;
    ///@}


    const BasisValues1d *operator->() const
    {
        return funcs_;
    }

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
    BasisValues1d const *funcs_;
};


template <int dim>
using ElemFuncValues = std::array<BasisValues1dConstView, dim>;

template <int dim>
class TensorProductFunctionEvaluator :
    public ElemFuncValues<dim>
{
public:
    void update_size(bool points_have_tensor_product_struct = true)
    {
        points_have_tensor_product_struct_ = points_have_tensor_product_struct;

        TensorSize<dim> n_func;
        TensorSize<dim> n_pts;
        for (int i = 0; i < dim; ++i)
        {
            n_func[i] = (*this)[i]->get_num_functions();
            n_pts[i] = (*this)[i]->get_num_points();
        }
        f_size_ = TensorSizedContainer<dim>(n_func);
        p_size_ = TensorSizedContainer<dim>(n_pts);

#ifndef NDEBUG
        if (!points_have_tensor_product_struct_)
        {
            // if the points have not a tensor product structure,
            // they must be the same number in all directions
            for (int i = 1; i < dim; ++i)
                Assert(n_pts[i] = n_pts[0],ExcDimensionMismatch(n_pts[i],n_pts[0]));
        }
#endif
    }
    /**
     * Evaluate and returns one partial derivative in one point.
     * The order of the partial derivative is specified by the tensor-index
     * @p order_tensor_id,
     * while the point is specified by its tensor-index @p point_tensor_id.
     */
    Real evaluate(const TensorIndex<dim> &order,
                  const TensorIndex<dim> &func,
                  const TensorIndex<dim> &pt) const
    {
        Real res = (dim>0) ? (*this)[0]->get_dataivative(order[0])(func[0],pt[0]) : 1.0;
        for (int i = 1; i < dim; ++i)
            res *= (*this)[i]->get_dataivative(order[i])(func[i], pt[i]);
        return res;
    }

    auto func_flat_to_tensor(const Index func_id) const
    {
        return f_size_.flat_to_tensor(func_id);
    }

    auto points_flat_to_tensor(const Index p_flat_id) const
    {
        if (points_have_tensor_product_struct_)
        {
            return p_size_.flat_to_tensor(p_flat_id);
        }
        else
        {
            return TensorIndex<dim>(p_flat_id);
            Assert(false,ExcNotImplemented());
        }
    }

private:
    TensorSizedContainer<dim> f_size_;
    TensorSizedContainer<dim> p_size_;

    /**
     * TRUE if the points are arranged in tensor product way.
     * In this case the total number of points is
     * <t>p_size_[0] * p_size_[1] * ... * p_size_[dim-1]</t>
     *
     * FALSE if the points are not arranged in tensor product way.
     * In this case it must hold
     * <t>p_size_[0] == p_size_[1] == ... == p_size_[dim-1]</t>
     * and each value along a specific direction refers to a single point.
     */
    bool points_have_tensor_product_struct_ = true;
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

