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
#ifndef QUADRATURE_LIB_H_
#define QUADRATURE_LIB_H_

#include <igatools/base/config.h>
#include <igatools/base/quadrature.h>

//#include <memory>

IGA_NAMESPACE_OPEN

/**
 * @brief Gauss-Legendre quadrature of arbitrary order.
 *
 * The coefficients of these quadrature rules are computed by the
 * function found in <tt>Numerical Recipies</tt>.
 *
 * @note this code was adapted from the dealii library for use
 * in igatools
 *
 * @ingroup eval_pts_scheme
 */
template< int dim >
class QGauss :
  public Quadrature< dim >
{
public:
  using typename Quadrature<dim>::Point;

  /**
   * Default constructor. Not allowed to be used.
   */
  QGauss() = delete ;

  /**
   * Constructor.
   * Builds a Gauss-Legendre quadrature scheme on the \f$ d \f$-dimensional hypercube \f$ [0,1]^d \f$
   * with \p num_points points in each coordinate direction.
   */
  explicit QGauss(const Size num_points);

  /**
   * Constructor.
   * Builds a Gauss-Legendre quadrature scheme on the \f$ d \f$-dimensional cube \f$ [0,1]^d \f$
   * with a (possibly) different number of points in each coordinate direction.
   * The number of points along the \p i-th coordinate direction is specified by \p num_points[i].
   */
  explicit QGauss(const TensorSize<dim> num_points);

  /**
   * Copy constructor. Performs a deep copy of the QGauss<dim> object.
   */
  QGauss(const QGauss< dim > &quad_scheme) = default ;

  /**
   * Copy assignment operator. Performs a deep copy of the QGauss<dim> object.
   */
  QGauss<dim> &operator=(const QGauss< dim > &quad_scheme) = default ;

  /**
   * Move constructor.
   */
  QGauss(QGauss< dim > &&quad_scheme) = default ;

  /**
   * Move assignment operator.
   */
  QGauss<dim> &operator=(QGauss< dim > &&quad_scheme) = default ;

  /**
   * Destructor.
   */
  ~QGauss() = default ;

  /**
   * Returns a Gauss-Legendre quadrature scheme (wrapped by a std::shared_ptr)
   * on the \f$ d \f$-dimensional hypercube \f$ [0,1]^d \f$
   * with \p num_points points in each coordinate direction.
   */
  static std::shared_ptr< QGauss< dim > >
  create(const Size num_points);

  /**
   * Returns a Gauss-Legendre quadrature scheme (wrapped by a std::shared_ptr)
   * on the \f$ d \f$-dimensional cube \f$ [0,1]^d \f$
   * with a (possibly) different number of points in each coordinate direction.
   * The number of points along the \p i-th coordinate direction is specified by \p num_points[i].
   */
  static std::shared_ptr< QGauss< dim > >
  create(const TensorSize<dim> num_points);


  /**
   * Returns a Gauss-Legendre quadrature scheme (wrapped by a std::shared_ptr)
   * on the \f$ d \f$-dimensional hypercube \f$ [0,1]^d \f$
   * with \p num_points points in each coordinate direction.
   */
  static std::shared_ptr< const QGauss< dim > >
  const_create(const Size num_points);

  /**
   * Returns a Gauss-Legendre quadrature scheme (wrapped by a std::shared_ptr)
   * on the \f$ d \f$-dimensional cube \f$ [0,1]^d \f$
   * with a (possibly) different number of points in each coordinate direction.
   * The number of points along the \p i-th coordinate direction is specified by \p num_points[i].
   */
  static std::shared_ptr< const QGauss< dim > >
  const_create(const TensorSize<dim> num_points);

} ;





/**
 * @brief Gauss-Lobatto quadrature of arbitrary order.
 *
 * @note this code was adapted from the dealii library for use
 * in igatools
 *
 * @ingroup eval_pts_scheme
 *
 */
template< int dim >
class QGaussLobatto :
  public Quadrature< dim >
{
public:
  using typename Quadrature<dim>::Point;

  /**
   * Default constructor. Not allowed to be used.
   */
  QGaussLobatto() = delete ;

  /**
   * Constructor.
   * Builds a Gauss-Lobatto quadrature scheme on the \f$ d \f$-dimensional cube \f$ [0,1]^d \f$
   * with \p num_points points in each coordinate direction.
   * The <p>eps</p> argument allows to perform a local scaling of the quadrature points.
   */
  explicit QGaussLobatto(const Size num_points);

  /**
   * Constructor.
   * Builds a Gauss-Lobatto quadrature scheme on the \f$ d \f$-dimensional cube \f$ [0,1]^d \f$
   * with a (possibly) different number of points in each coordinate direction.
   * The number of points along the \p i-th coordinate direction is specified by \p num_points[i].
   * The <p>eps</p> argument allows to perform a local scaling of the quadrature points.
   */
  explicit QGaussLobatto(const TensorSize<dim> num_points);

  /**
   * Copy constructor. Performs a deep copy of the QGaussLobatto<dim> object.
   */
  QGaussLobatto(const QGaussLobatto< dim > &quad_scheme) = default ;

  /**
   * Copy assignment operator. Performs a deep copy of the QGaussLobatto<dim> object.
   */
  QGaussLobatto< dim > &operator=(const QGaussLobatto< dim > &quad_scheme) = default ;

  /**
   * Move constructor.
   */
  QGaussLobatto(QGaussLobatto< dim > &&quad_scheme) = default ;

  /**
   * Move assignment operator.
   */
  QGaussLobatto< dim > &operator=(QGaussLobatto< dim > &&quad_scheme) = default ;

  /**
   * Destructor.
   */
  ~QGaussLobatto() = default;

  /**
   * Returns a Gauss-Lobatto quadrature scheme (wrapped by a std::shared_ptr)
   * on the \f$ d \f$-dimensional hypercube \f$ [0,1]^d \f$
   * with \p num_points points in each coordinate direction.
   */
  static std::shared_ptr< QGaussLobatto< dim > >
  create(const Size num_points);

  /**
   * Returns a Gauss-Lobatto quadrature scheme (wrapped by a std::shared_ptr)
   * on the \f$ d \f$-dimensional cube \f$ [0,1]^d \f$
   * with a (possibly) different number of points in each coordinate direction.
   * The number of points along the \p i-th coordinate direction is specified by \p num_points[i].
   */
  static std::shared_ptr< QGaussLobatto< dim > >
  create(const TensorSize<dim> num_points);
} ;



/**
 * @brief Uniform quadrature rule.
 *
 * @ingroup eval_pts_scheme
 */
template< int dim >
class QUniform :
  public Quadrature<dim>
{
public:
  using typename Quadrature<dim>::Point;

  /**
   * Default constructor. Not allowed to be used.
   */
  QUniform() = delete ;

  /**
   * Constructor.
   * Builds a uniform quadrature scheme on the \f$ d \f$-dimensional cube \f$ [0,1]^d \f$
   * with \p num_points points in each coordinate direction.
   * The eps argument allows to do a local scaling of the quadrature points.
   */
  explicit QUniform(const Size num_points) ;

  /**
   * Constructor.
   * Builds a uniform quadrature scheme on the \f$ d \f$-dimensional cube \f$ [0,1]^d \f$
   * with a (possibly) different number of points in each coordinate direction.
   * The number of points along the \p i-th coordinate direction is specified by \p num_points[i].
   * The eps argument allows to do a local scaling of the quadrature points.
   */
  explicit QUniform(const TensorSize<dim> num_points) ;

  /**
   * Copy constructor. Performs a deep copy of the QUniform<dim> object.
   */
  QUniform(const QUniform< dim > &quad_scheme) = default ;

  /**
   * Copy assignment operator. Performs a deep copy of the QUniform<dim> object.
   */
  QUniform< dim > &operator=(const QUniform< dim > &quad_scheme) = default ;

  /**
   * Move constructor.
   */
  QUniform(QUniform< dim > &&quad_scheme) = default ;

  /**
   * Move assignment operator.
   */
  QUniform< dim > &operator=(QUniform< dim > &&quad_scheme) = default ;

  /**
   * Destructor.
   */
  ~QUniform() = default ;

  /**
   * Returns a uniform quadrature scheme (wrapped by a std::shared_ptr) on
   * the \f$ d \f$-dimensional cube \f$ [0,1]^d \f$ with \p num_points
   * points in each coordinate direction.
   * The eps argument allows to do a local scaling of the quadrature points.
   */
  static std::shared_ptr< QUniform< dim > >
  create(const Size num_points) ;

  /**
   * Returns a uniform quadrature scheme (wrapped by a std::shared_ptr) on
   * the \f$ d \f$-dimensional cube \f$ [0,1]^d \f$ with a (possibly)
   * different number of points in each coordinate direction.
   * The number of points along the \p i-th coordinate direction is specified by \p num_points[i].
   * The eps argument allows to do a local scaling of the quadrature points.
   */
  static std::shared_ptr< QUniform< dim > >
  create(const TensorSize<dim> num_points) ;
} ;


/**
 * @brief Trapezoidal quadrature rule, exact for linear polynomials.
 *
 * @ingroup eval_pts_scheme
 */
template <int dim>
class QTrapez :
  public QUniform<dim>
{
public:
  /**
   * Constructor.
   * Builds a trapezoidal quadrature scheme on the \f$ d \f$-dimensional cube \f$ [0,1]^d \f$,
   * exact for linear polynomials.
   * The eps argument allows to do a local scaling of the quadrature points.
   */
  explicit QTrapez();

  /**
   * Copy constructor. Performs a deep copy of the QTrapez<dim> object.
   */
  QTrapez(const QTrapez< dim > &quad_scheme) = default ;

  /**
   * Copy assignment operator. Performs a deep copy of the QTrapez<dim> object.
   */
  QTrapez< dim > &operator=(const QTrapez< dim > &quad_scheme) = default ;

  /**
   * Move constructor.
   */
  QTrapez(QTrapez< dim > &&quad_scheme) = default ;

  /**
   * Move assignment operator.
   */
  QTrapez< dim > &operator=(QTrapez< dim > &&quad_scheme) = default ;

  /**
   * Destructor.
   */
  ~QTrapez() = default ;

  /**
   * Returns a uniform quadrature scheme (wrapped by a std::shared_ptr) on
   * the \f$ d \f$-dimensional cube \f$ [0,1]^d \f$, exact for linear
   * polynomials.
   * The eps argument allows to do a local scaling of the quadrature points.
   */
  static std::shared_ptr< QTrapez< dim > >
  create() ;

} ;


IGA_NAMESPACE_CLOSE

#endif // #ifndef QUADRATURE_LIB_H_
