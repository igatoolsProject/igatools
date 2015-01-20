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

#ifndef QUADRATURE_H_
#define QUADRATURE_H_

#include <igatools/base/config.h>
#include <igatools/utils/array.h>
#include <igatools/utils/tensor_product_array.h>

IGA_NAMESPACE_OPEN


template <int dim_>
class EvaluationPoints
{
public:

    using Point = Points<dim_>;

    static const int dim = dim_;

    /**
     * Returns the total number of evaluation points.
     */
    int get_num_points() const;

    /**
     * Returns TRUE if the evaluation points have a tensor-product structure.
     */
    virtual bool is_tensor_product_struct() const;

    /**
     * Returns the coordinates indices relative to the point with (flat) index <p>point_id</p>.
     */
    TensorIndex<dim_> get_coords_id_from_point_id(const int point_id) const;


    /**
     * Returns the number of point coordinates along each direction.
     */
    TensorIndex<dim_> get_num_coords_direction() const noexcept;


    /**
     * @name Constructors.
     */
    ///@{
    /**
     * Default constructor. It sets a unit bounding-box (i.e. the hypercube \f$ [0,1]^{dim}\f$) with no points inside.
     */
    EvaluationPoints();

    /**
     * Construct the object given a vector of points in the <t>dim_</t>-dimensional space.
     */
    EvaluationPoints(const ValueVector<Point> &pts);

    /**
     * Copy constructor.
     */
    EvaluationPoints(const EvaluationPoints<dim_> &pts) = default;

    /**
     * Move constructor.
     */
    EvaluationPoints(EvaluationPoints<dim_> &&pts) = default;

    /**
     * Destructor.
     */
    virtual ~EvaluationPoints() = default;
    ///@}


    /**
     * @name Assignment operators.
     */
    ///@{
    /**
     * Copy assignment operator.
     */
    EvaluationPoints<dim_> &operator=(const EvaluationPoints<dim_> &pts) = default;

    /**
     * Move assignment operator.
     */
    EvaluationPoints<dim_> &operator=(EvaluationPoints<dim_> &&pts) = default;
    ///@}



    /**
     * Returns coordinates of the points along the <p>i</p>-th direction.
     */
    const vector<Real> &get_coords_direction(const int i) const;


    /**
     * Returns the coordinates of the the point with (flat) index <p>pt_id</p>
     */
    Point get_point(const int pt_id) const;

    /**
     * Returns the weight of the the point with (flat) index <p>pt_id</p>
     */
    Real get_weight(const int pt_id) const;


    /**
     * Returns the bounding box in which the points are located.
     */
    const BBox<dim_> &get_bounding_box() const;

    /**
     * @name Functions for performing dilation and translation of the points (and weights).
     */
    //@{
    /**
     * Dilation of the points (and of the corresponding bounding box)
     */
    void dilate(const Point &dilate);

    /**
     * Translation of the points (and of the corresponding bounding box)
     */
    void translate(const Point &translate);

    /**
     * Dilation followed by a translation of the points (and of the corresponding bounding box).
     */
    void dilate_translate(const Point &dilate, const Point &translate);
    ///@}


    /**
     * Reset the bounding box in which the points must be located.
     */
    void reset_bounding_box(const BBox<dim_> &bounding_box);

    /**
     * This function performs the following task:
     *   - reset the value of the points coordinates and the map point_id_to_coords_is,
     * given a vector of points in the <t>dim_</t>-dimensional space
     *   - reset the weight value associated to each point.
     *
     * @note In DEBUG mode the points are tested if they are within the bounding box.
     * Points on the side of the bounding box are still valid points.
     */
    void reset_points_coordinates_and_weights(const ValueVector<Point> &pts,const ValueVector<Real> &weights);


private:

    /**
     * Coordinates of the points.
     *
     * It does not contain multiple values.
     */
    special_array< vector<Real>, dim_> coordinates_;


    /**
     * Weights of the points. By default are set to be equal to one.
     */
    ValueVector<Real> weights_;

    /**
     * Map between the point (flat) ids and its coordinates ids.
     */
    vector<TensorIndex<dim_>> map_point_id_to_coords_id_;

    BBox<dim_> bounding_box_;

    bool is_tensor_product_struct_;



};


/**
 * @brief Base class for tensor product quadrature formulas with tensor-product structure,
 *  in @p dim dimensions.
 *
 * This class stores quadrature points and weights on a dim-dimensional unit
 * hypercube
 * \f$ [0,1]^\text{dim} \f$.
 *
 * It also provides methods for projecting the points on the face and for building
 * quadrature extension in the hypercube
 * \f$ [0,1]^\text{dim+1} \f$.
 *
 */
template<int dim_>
class QuadratureTensorProduct
    : public EvaluationPoints<dim_>
{
private:
    using self_t = QuadratureTensorProduct<dim_>;
public:

    using typename EvaluationPoints<dim_>::Point;
    using EvaluationPoints<dim_>::dim;

    using WeigthArray = TensorProductArray<dim>;
    using PointArray  = CartesianProductArray<Real, dim>;
public:
    ///@name Constructors
    ///@{
    /**
     * Default constructor. It does nothing.
     */
    QuadratureTensorProduct() = default;

    /**
     * Creates a tensor product quadrature rule
     * on the unit d-dimensional hypercube \f$ [0,1]^d \f$,
     * with @p num_points[i] number of points in the i-th dimension.
     *
     * The purpose of the function pointer argument <p>compute_coords_and_weight_1d</p>
     * is to let a derived class to pass the method to compute the 1D coordinates and weights.
     *
     * The <p>eps</p> argument allows to perform a local scaling of the quadrature points.
     */
    explicit QuadratureTensorProduct(
        const TensorSize<dim> num_points,
        void (*compute_coords_and_weight_1d)(const int n_pts_1d, vector<Real> &coords,vector<Real> &weights),
        const Real eps_scaling = 0.0);

    /**
     * Creates a tensor product quadrature rule with the points, the
     * weights and the domain coordinates of the d-dimensional hypercube
     * upon which the quadrature is referred to.
     */
    explicit QuadratureTensorProduct(const PointArray &points,
                                     const WeigthArray &weights);

    /**
     * Destructor.
     */
    virtual ~QuadratureTensorProduct() = default;

    /**
     * Copy constructor.
     * It performs a deep copy of the Quadrature object.
     */
    QuadratureTensorProduct(const self_t &quad_scheme) = default;

    /**
    * Move constructor.
    */
    QuadratureTensorProduct(self_t &&quad_scheme) = default;
    ///@}

    ///@name Assignment operators
    ///@{
    /**
     * Copy assignment operator.
     * It performs a deep copy of the Quadrature object.
     */
    QuadratureTensorProduct<dim> &operator=(const self_t &quad_scheme) = default;

    /**
     * Move assignment operator.
     */
    QuadratureTensorProduct<dim> &operator=(self_t  &&quad_scheme) = default;
    ///@}

    ///@name Getting informations about the points and the domain.
    ///@{
    /**
     * Returns TRUE.
     */
    virtual bool is_tensor_product_struct() const override final;

public:
    ///@}

    ///@name Getting a copy of the internal variables.
    ///@{
    /**
     * Return all quadrature weights in
     * a tensor-product structure.
     */
    WeigthArray get_weights() const noexcept ;

    /**
     * Return all quadrature points in
     * a tensor-product structure.
     */
    PointArray get_points() const noexcept;
    ///@}

    /**
     * Returns a dim dimensional quadrature obtained by using
     * a single point on the active face direction.
     * @todo write example
     * Usually use for face values
     */
    template<int k>
    self_t collapse_to_sub_element(const int id) const;

    /**
     * Prints internal information about the quadrature scheme.
     * \note Mostly used for debugging and testing.
     */
    void print_info(LogStream &out) const;

protected:
    /**
     * Quadrature points.
     */
    PointArray points_;

    /**
     * Quadrature weights.
     */
    WeigthArray weights_;

};


/**
 * Given evaluation points on a dim dimensional face, of a dim+1
 * domain, this functions creates an extended dimension
 * version of this points applicable to the volume.
 *
 * @relates Quadrature
 */
template<int k, int dim>
EvaluationPoints<dim>
extend_sub_elem_quad(const EvaluationPoints<k> &quad, const int sub_elem_id);



/**
 * Given a quadrature rule on a dim dimensional face, of a dim+1
 * domain, this functions creates an extended dimension
 * version of this quadrature applicable to the volume.
 *
 * @relates Quadrature
 */
template<int k, int dim>
QuadratureTensorProduct<dim>
extend_sub_elem_quad(const QuadratureTensorProduct<k> &quad, const int sub_elem_id);


IGA_NAMESPACE_CLOSE

#endif // #ifndef QUADRATURE_H_
