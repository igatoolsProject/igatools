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

#ifndef QUADRATURE_TENSOR_PRODUCT_H_
#define QUADRATURE_TENSOR_PRODUCT_H_

#include <igatools/base/evaluation_points.h>

IGA_NAMESPACE_OPEN
#if 0

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
 * @note This class in not intended to be used/instantiated directly
 * (in fact, it has no public constructors).
 * It's main purpose is to be the base class for QGauss, QGaussLobatto, QUniform.
 */
template<int dim_>
class Quadrature
    : public Quadrature<dim_>
{
private:
    using parent_t = Quadrature<dim_>;
    using self_t = Quadrature<dim_>;
public:

    using typename Quadrature<dim_>::Point;
    using typename Quadrature<dim_>::PointVector;
    using typename Quadrature<dim_>::PointArray;
    using typename Quadrature<dim_>::WeightArray;

protected:

    ///@name Constructors
    ///@{
    /**
     * Default constructor. It sets the bounding-box to be the hypercube \f$ [0,1]^{dim}\f$ with no points inside.
     */
    Quadrature();

#if 0
    /**
     * Creates a tensor product quadrature rule
     * on the unit d-dimensional hypercube \f$ [0,1]^d \f$,
     * with @p num_points[i] number of points in the i-th dimension.
     *
     * The purpose of the function pointer argument
     * <p>compute_coords_and_weight_1d</p>
     * is to let a derived class to pass the method to compute
     * the 1D coordinates and weights.
     *
     * The <p>eps</p> argument allows to perform a local scaling
     * of the quadrature points.
     */
    explicit Quadrature(
        const TensorSize<dim_> num_points,
        void (*compute_coords_and_weight_1d)
        (const int n_pts_1d, vector<Real> &coords,vector<Real> &weights),
        const Real eps_scaling = 0.0);

public:
    /**
     * Creates a tensor product quadrature rule with the points, the
     * weights and the domain coordinates of the d-dimensional hypercube
     * upon which the quadrature is referred to.
     */
    explicit Quadrature(
        const PointVector &points,
        const WeightArray &weights_1d,
        const BBox<dim_> &bounding_box);
#endif

public:
    explicit Quadrature(
    /**
     * Destructor.
     */
    ~Quadrature() = default;


    /**
     * Copy constructor.
     * It performs a deep copy of the Quadrature object.
     */
    Quadrature(const self_t &quad_scheme) = default;

    /**
    * Move constructor.
    */
    Quadrature(self_t &&quad_scheme) = default;
    ///@}

    ///@name Assignment operators
    ///@{
    /**
     * Copy assignment operator.
     * It performs a deep copy of the Quadrature object.
     */
    Quadrature<dim_> &operator=(const self_t &quad_scheme) = default;

    /**
     * Move assignment operator.
     */
    Quadrature<dim_> &operator=(self_t  &&quad_scheme) = default;
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
    void (*compute_coords_and_weight_1d)
       (const int n_pts_1d, vector<Real> &coords,vector<Real> &weights);
};



#endif

IGA_NAMESPACE_CLOSE

#endif // #ifndef QUADRATURE_TENSOR_PRODUCT_H_
