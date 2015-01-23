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
    bool have_points_tensor_product_struct() const;

    /**
     * Returns TRUE if the weights have a tensor-product structure.
     */
    bool have_weights_tensor_product_struct() const;

    /**
     * Returns the coordinates indices relative to the point with (flat) index <p>point_id</p>.
     */
    TensorIndex<dim_> get_coords_id_from_point_id(const int point_id) const;


    /**
     * Returns the number of point coordinates along each direction.
     */
    TensorIndex<dim_> get_num_coords_direction() const noexcept;


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
     * Construct the object given a vector of points in the <t>dim_</t>-dimensional space,
     * and assign the weights value to be equal to 1.
     *
     * @note It sets the bounding-box to be the hypercube \f$ [0,1]^{dim}\f$.
     */
    EvaluationPoints(const ValueVector<Point> &pts);

    /**
     * Construct the object given:
     *   - a vector of <p>points</p> in the <t>dim_</t>-dimensional space;
     *   - the <p>weights_1d</p> are the weights associated to the points coordinates;
     *   - the bounding box in which the points are defined.
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
     * Returns all the points.
     */
    ValueVector<Point> get_points() const;

    /**
     * Returns all the weights.
     */
    ValueVector<Real> get_weights() const;

    const special_array<vector<Real>,dim_> &get_weights_1d() const;

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
     * Returns a dim dimensional quadrature obtained by using
     * a single point on the active face direction.
     * @todo write example
     * Usually use for face values
     */
    template<int k>
    EvaluationPoints<dim_> collapse_to_sub_element(const int id) const;


protected:

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

    using WeightArray = TensorProductArray<dim>;
    using PointArray  = CartesianProductArray<Real, dim>;

    ///@name Constructors
    ///@{
    /**
     * Default constructor. It sets the bounding-box to be the hypercube \f$ [0,1]^{dim}\f$ with no points inside.
     */
    QuadratureTensorProduct();

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
                                     const WeightArray &weights,
                                     const BBox<dim> &bounding_box);


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

#if 0
    ///@name Getting informations if the points and weights have the tensor-product structure.
    ///@{
    /**
     * Returns TRUE.
     */
    virtual bool have_points_tensor_product_struct() const override final;

    /**
     * Returns TRUE.
     */
    virtual bool have_weights_tensor_product_struct() const override final;
    ///@}
#endif

#if 0
    /**
     * Returns a dim dimensional quadrature obtained by using
     * a single point on the active face direction.
     * @todo write example
     * Usually use for face values
     */
    template<int k>
    self_t collapse_to_sub_element(const int id) const;
#endif


    /**
     * Return the weights with their underlying tensor-product structure.
     * @return
     */
    const WeightArray &get_weights_tensor_product() const;

protected:

    /**
     * Quadrature weights.
     *
     * @note We store this information keeping its tensor-product structure
     * because we need it when the points must be extended/collapsed to a face.
     */
    WeightArray weights_;

    /**
     * This function pointer is pointing to the function that performs
     * the 1D points and weights computation.
     * @param[in] n_pts_1d Number of points along one direction.
     * @param[out] coords Point coordinates.
     * @param[out] weights Weights associated to the 1D point
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


#if 0
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
#endif

IGA_NAMESPACE_CLOSE

#endif // #ifndef QUADRATURE_H_
