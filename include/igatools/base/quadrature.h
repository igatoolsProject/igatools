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
#ifndef QUADRATURE_H_
#define QUADRATURE_H_

#include <igatools/base/config.h>
#include <igatools/utils/tensor_product_array.h>

IGA_NAMESPACE_OPEN

/**
 * @brief Base class for tensor product quadrature formulas in @p dim dimensions.
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
template<int dim>
class Quadrature
{
public:
    using QuadArray = TensorProductArray<dim>;

    ///@name Constructors
    ///@{
    /**
     * Default constructor. It does nothing.
     */
    Quadrature() = default;

    /**
     * Creates a tensor product quadrature rule
     * on the unit d-dimensional hypercube \f$ [0,1]^d \f$,
     * with @p num_points[i] number of points in the i-th dimension.
     */
    explicit Quadrature(const TensorSize<dim > num_points);

    /**
     * Creates a tensor product quadrature rule
     * on the unit d-dimensional hypercube \f$ [0,1]^d \f$,
     * with @p num_points points in each direction.
     */
    explicit Quadrature(const Index num_points);

    /**
     * Creates a tensor product quadrature rule with the points, the
     * weights and the domain coordinates of the d-dimensional hypercube
     * upon which the quadrature is referred to.
     */
    explicit Quadrature(const CartesianProductArray<Real,dim> &points,
                        const TensorProductArray<dim> &weights);

    /**
     * Destructor.
     */
    ~Quadrature() = default;

    /**
     * Copy constructor.
     * It performs a deep copy of the Quadrature object.
     */
    Quadrature(const Quadrature<dim> &quad_scheme) = default;

    /**
    * Move constructor.
    */
    Quadrature(Quadrature<dim> &&quad_scheme) = default;
    ///@}

    ///@name Assignment operators
    ///@{
    /**
     * Copy assignment operator.
     * It performs a deep copy of the Quadrature object.
     */
    Quadrature<dim> &operator=(const Quadrature< dim > &quad_scheme) = default;

    /**
     * Move assignment operator.
     */
    Quadrature<dim> &operator=(Quadrature< dim > &&quad_scheme) = default;
    ///@}

    ///@name Getting informations about the points and the domain.
    ///@{
    /**
     * Returns the total number of quadrature points.
     */
    Size get_num_points() const noexcept;

    /**
     * Returns the number of quadrature points along each coordinate direction.
     */
    TensorSize<dim> get_num_points_direction() const noexcept;
    ///@}

    ///@name Getting a copy of the internal variables.
    ///@{
    /**
     * Return all quadrature weights in
     * a tensor-product structure.
     */
    TensorProductArray<dim> get_weights() const noexcept ;

    /**
     * Return all quadrature points in
     * a tensor-product structure.
     */
    CartesianProductArray<Real,dim> get_points() const noexcept;

    ///@}


    /**
     * @name Function returning restriction and extension quadrature (useful for face evaluation)
     */
    ///@{

    /**
     * Returns a quadrature obtained by projecting the quadrature to the face
     * identified by the input argument @p face_id.
     */
    Quadrature<dim> get_restriction(const int face_id) const;


    /**
     * Returns an extended dimension version of this quadrature.
     */
    Quadrature<dim+1> get_extension(const int face_id) const;
    ///@}

    /**
     * Prints internal information about the quadrature scheme.
     * \note Mostly used for debugging and testing.
     */
    void print_info(LogStream &out) const;

protected:
    /**
     * Quadrature points.
     */
    CartesianProductArray<Real,dim> points_ ;

    /**
     * Quadrature weights.
     */
    TensorProductArray<dim> weights_ ;


//    /**
//     * Coordinates of the domain (i.e. d-dimensional hypercube) in which the quadrature is referred to.
//     */
//    std::array< std::array<Real,2>, dim> domain_coordinates_ ;
//
//    /**
//     * This flags specifies if (and in which component) the point coordinates have constant value.
//     *
//     * * If all the entries of this variable are false, then the all points do not belong to any lower dimension element.
//     * * If the entry @b i is true and the others are false, then the points coordinates have the i-th coordinate fixed, i.e. the points
//     * belongs to an element with dimension dim-1, and so on.
//     *
//     * @note This variable is useful for knowing if the points belong to an element face.
//     */
//    std::array<bool,dim> info_constant_coordinates_ ;

} ;





IGA_NAMESPACE_CLOSE

#endif // #ifndef QUADRATURE_H_
