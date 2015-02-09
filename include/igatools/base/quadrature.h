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
#include <igatools/base/tensor.h>
#include <igatools/utils/value_vector.h>
//#include <igatools/utils/tensor_product_array.h>

IGA_NAMESPACE_OPEN

/**
 * @brief This class represents a container for an user-definable set of evaluation points.
 *
 *
 * Its main constructor takes as input arguments the following parameters:
 * - a vector of points in the space;
 * - the weights associated to the coordinates of the points;
 * - the bounding box enclosing the points
 *   (this is useful if some affine transformation are to be performed).
 *
 * The way in which this class is implemented allows to be used both for points having
 * a tensor-product structure (as for the multi-dimensional Gaussian quadrature scheme)
 * or for general points over a <tt>dim</tt>-dimensional domain.
 *
 * In order to do so, the class does not store the points itself but,
 * after an analysis/preprocessing stage (done in the constructor),
 * it stores the points coordinates and a multi-index for each point, that represents the
 * coordinates id for a given point.
 * For example, in 2D case, let suppose that the points used in the constructor argument are
 *
 * \f[
 * \mathbf{P}_0 = (0.1,0.2) \quad , \quad
 * \mathbf{P}_1 = (0.9,0.6) \quad , \quad
 * \mathbf{P}_2 = (0.5,0.8) \quad , \quad
 * \mathbf{P}_3 = (0.1,0.8) \quad , \quad
 * \mathbf{P}_4 = (0.9,0.2)
 * \f]
 * then this class stores internally the two vectors of (sorted) coordinates
 * \f{align*}{
 * \mathtt{coordinates\_[0]} &= \{0.1,0.5,0.9\} \\
 * \mathtt{coordinates\_[1]} &= \{0.2,0.6,0.8\}
 * \f}
 * and the coordinates id of each point
 * \f{align*}{
 * \mathtt{map\_point\_id\_to\_coords\_id\_[0]} &= \{ 1, 1\} \\
 * \mathtt{map\_point\_id\_to\_coords\_id\_[1]} &= \{ 3, 2\} \\
 * \mathtt{map\_point\_id\_to\_coords\_id\_[2]} &= \{ 2, 3\} \\
 * \mathtt{map\_point\_id\_to\_coords\_id\_[3]} &= \{ 1, 3\} \\
 * \mathtt{map\_point\_id\_to\_coords\_id\_[4]} &= \{ 3, 1\}
 * \f}
 *
 * @ingroup eval_pts_scheme
 *
 * @author M. Martinelli
 * @date 2014, 2015
 */
template <int dim_>
class EvaluationPoints
{
public:

    /**
     * @brief Alias for the point-type that is returned by the function EvaluationPoints::get_point()
     */
    using Point = Points<dim_>;


    static const int dim = dim_;

    /**
     * @name Functions returning information about the number of points and coordinates.
     */
    ///@{
    /**
     * Returns the total number of evaluation points.
     */
    int get_num_points() const;

    /**
     * Returns the number of point coordinates along each direction.
     */
    TensorIndex<dim_> get_num_coords_direction() const noexcept;
    ///@}



    /**
     * @name Functions returning informations about the arrangement structure of the points and weights.
     */
    ///@{
    /**
     * Returns TRUE if the evaluation points have a tensor-product structure.
     */
    bool have_points_tensor_product_struct() const;

    /**
     * Returns TRUE if the weights have a tensor-product structure.
     */
    bool have_weights_tensor_product_struct() const;
    ///@}




    /**
     * Prints internal information about the evaluation points.
     * \note Mostly used for debugging and testing.
     */
    void print_info(LogStream &out) const;


    /**
     * @name Constructors.
     */
    ///@{
    /**
     * Default constructor. It sets the bounding-box to be the hypercube \f$ [0,1]^{dim}\f$ with no points inside.
     */
    EvaluationPoints();

protected:
    /**
     * Construct the object with a user-defined bounding-box, with no points inside.
     */
    EvaluationPoints(const BBox<dim> &bounding_box);

public:
    /**
     * Construct the object given a vector of <tt>points</tt>
     * in the <tt>dim_</tt>-dimensional space,
     * and assign the weights value to be equal to 1.
     *
     * @note It sets the bounding-box to be the hypercube \f$ [0,1]^{dim}\f$.
     */
    EvaluationPoints(const ValueVector<Point> &points);

    /**
     * Construct the object given:
     * - a vector of <tt>points</tt> in the <tt>dim_</tt>-dimensional space;
     * - the <tt>weights_1d</tt> are the weights associated to the points coordinates;
     * - the <tt>bounding_box</tt> in which the points are defined.
     */
    EvaluationPoints(
        const ValueVector<Point> &points,
        const special_array<vector<Real>,dim_> &weights_1d,
        const BBox<dim> &bounding_box);

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
    ~EvaluationPoints() = default;
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
     * @name Functions returning internal data (points, weights, coordinates, bounding box, etc.)
     */
    ///@{
    /**
     * Returns all the points.
     */
    ValueVector<Point> get_points() const;

    /**
     * Returns the coordinates of the the point with (flat) index <tt>pt_id</tt>
     */
    Point get_point(const int pt_id) const;

    /**
     * Returns coordinates of the points along the <tt>i</tt>-th direction.
     */
    const vector<Real> &get_coords_direction(const int i) const;

    /**
     * Returns all the weights.
     */
    ValueVector<Real> get_weights() const;

    /**
     * Returns the weight of the the point with (flat) index <tt>pt_id</tt>
     */
    Real get_weight(const int pt_id) const;


    const special_array<vector<Real>,dim_> &get_weights_1d() const;

    /**
     * Returns the bounding box in which the points are located.
     */
    const BBox<dim_> &get_bounding_box() const;

    /**
     * Returns the coordinates indices relative to the point with (flat) index <tt>point_id</tt>.
     */
    TensorIndex<dim_> get_coords_id_from_point_id(const int point_id) const;
    ///@}


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
     * Returns a dim dimensional quadrature obtained by using
     * a single point on the active face direction.
     * @todo write example
     * Usually use for face values
     */
    template<int k>
    EvaluationPoints<dim_> collapse_to_sub_element(const int id) const;


private:

    /**
     * Reset the <tt>bounding_box</tt> in which the points must be located.
     */
    void reset_bounding_box(const BBox<dim_> &bounding_box);

protected:

    /**
     * This function performs the following task:
     *   - reset the value of the points coordinates and the map point_id_to_coords_is,
     * given a vector of points in the <t>dim_</t>-dimensional space
     *   - reset the weight value associated to each point.
     *
     * @note In DEBUG mode the points are tested if they are within the bounding box.
     * Points on the side of the bounding box are still valid points.
     */
    void reset_points_coordinates_and_weights(
        const ValueVector<Point> &pts,
        const special_array<vector<Real>,dim_> &weights_1d);



    /**
     * Coordinates of the points.
     *
     * It does not contain multiple values.
     */
    special_array< vector<Real>, dim_> coordinates_;


    /**
     * Weights associated to the points coordinates. By default are set to be equal to one.
     */
    special_array<vector<Real>,dim_> weights_1d_;


    /**
     * Map between the point (flat) ids and its coordinates ids.
     */
    vector<TensorIndex<dim_>> map_point_id_to_coords_id_;


    BBox<dim_> bounding_box_;

    bool  points_have_tensor_product_struct_ = false;

    bool weights_have_tensor_product_struct_ = false;
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
 * @note This class in not intended to be used/instantiated directly .
 * It's main purpose is to be the base class for QGauss, QGaussLobatto, QUniform.
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

    ///@name Constructors
    ///@{
    /**
     * Default constructor. It sets the bounding-box to be the hypercube \f$ [0,1]^{dim}\f$ with no points inside.
     */
    QuadratureTensorProduct();

protected:

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

public:
    /**
     * Creates a tensor product quadrature rule with the points, the
     * weights and the domain coordinates of the d-dimensional hypercube
     * upon which the quadrature is referred to.
     */
    explicit QuadratureTensorProduct(
        const ValueVector<Point> &points,
        const special_array<vector<Real>,dim_> &weights_1d,
        const BBox<dim> &bounding_box);


    /**
     * Destructor.
     */
    ~QuadratureTensorProduct() = default;


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

protected:

    /**
     * This function pointer is pointing to the function that performs
     * the 1D points and weights computation.
     * @param[in] n_pts_1d Number of points along one direction.
     * @param[out] coords Point coordinates.
     * @param[out] weights Weights associated to the 1D point
     *
     * @note It is assumed that the 1D points and weights along each coordinate component
     * are computed using the same algorithm (but possibly with different number of points).
     */
    void (*compute_coords_and_weight_1d)(const int n_pts_1d, vector<Real> &coords,vector<Real> &weights);
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



IGA_NAMESPACE_CLOSE

#endif // #ifndef QUADRATURE_H_
