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

#ifndef __GRID_FUNCTION_LIB_H_
#define __GRID_FUNCTION_LIB_H_

#include <igatools/geometry/formula_grid_function.h>

IGA_NAMESPACE_OPEN

/**
 * Collection of useful functions derived from the GridFunction class.
 */
namespace grid_functions
{

/**
 * y = b
 */
template<int dim, int space_dim>
class ConstantGridFunction :
  public FormulaGridFunction<dim,space_dim>
{
  using base_t = GridFunction<dim,space_dim>;
  using parent_t = FormulaGridFunction<dim,space_dim>;
  using self_t = ConstantGridFunction<dim,space_dim>;
  using typename base_t::GridType;
public:
  using typename parent_t::Value;
  using typename parent_t::GridPoint;
  template <int order>
  using Derivative = typename parent_t::template Derivative<order>;

public:
  static std::shared_ptr<self_t>
  create(const std::shared_ptr<GridType> &domain,
         const Value &b);

  static std::shared_ptr<const self_t>
  const_create(const std::shared_ptr<const GridType> &domain,
               const Value &b);

  ConstantGridFunction(const self_t &) = default;

  virtual ~ConstantGridFunction() = default;

  virtual void print_info(LogStream &out) const override final;

  const Value &get_constant_value() const;

protected:
  ConstantGridFunction(
    const SharedPtrConstnessHandler<GridType> &domain,
    const Value &b);

private:
  void evaluate_0(const ValueVector<GridPoint> &points,
                  ValueVector<Value> &values) const override;

  void evaluate_1(const ValueVector<GridPoint> &points,
                  ValueVector<Derivative<1>> &values) const override;

  void evaluate_2(const ValueVector<GridPoint> &points,
                  ValueVector<Derivative<2>> &values) const override;

  const Value b_;

#ifdef MESH_REFINEMENT
  /**
   * Rebuild the internal state of the object after an Grid::insert_knots() function is invoked.
   *
   * @pre Before invoking this function, must be invoked the function grid_->insert_knots().
   * @note This function is connected to the Grid's signal for the refinement, and
   * it is necessary in order to avoid infinite loops in the insert_knots() function calls.
   *
   * @ingroup h_refinement
   */
  virtual void rebuild_after_insert_knots(
    const SafeSTLArray<SafeSTLVector<Real>,dim> &knots_to_insert,
    const Grid<dim> &old_grid) override final;
#endif // MESH_REFINEMENT

};


//------------------------------------------------------------------------------
/**
 * F(x) = A * x + b
 */
template<int dim, int space_dim>
class LinearGridFunction :
  public FormulaGridFunction<dim,space_dim>
{
  using base_t = GridFunction<dim,space_dim>;
  using parent_t = FormulaGridFunction<dim,space_dim>;
  using self_t = LinearGridFunction<dim,space_dim>;
  using typename base_t::GridType;
public:
  using typename parent_t::Value;
  using typename parent_t::GridPoint;
  template <int order>
  using Derivative = typename parent_t::template Derivative<order>;

public:
  static std::shared_ptr<self_t>
  create(const std::shared_ptr<GridType> &domain,
         const Derivative<1> &A,
         const Value &b);

  static std::shared_ptr<const self_t>
  const_create(const std::shared_ptr<const GridType> &domain,
               const Derivative<1> &A,
               const Value &b);

  LinearGridFunction(const self_t &) = default;

  virtual ~LinearGridFunction() = default;

  virtual void print_info(LogStream &out) const override final;

protected:
  LinearGridFunction(
    const SharedPtrConstnessHandler<GridType> &domain,
    const Derivative<1> &A,
    const Value &b);

private:
  void evaluate_0(const ValueVector<GridPoint> &points,
                  ValueVector<Value> &values) const override;

  void evaluate_1(const ValueVector<GridPoint> &points,
                  ValueVector<Derivative<1>> &values) const override;

  void evaluate_2(const ValueVector<GridPoint> &points,
                  ValueVector<Derivative<2>> &values) const override;

  const Derivative<1> A_;
  const Value    b_;

#ifdef MESH_REFINEMENT
  /**
   * Rebuild the internal state of the object after an Grid::insert_knots() function is invoked.
   *
   * @pre Before invoking this function, must be invoked the function grid_->insert_knots().
   * @note This function is connected to the Grid's signal for the refinement, and
   * it is necessary in order to avoid infinite loops in the insert_knots() function calls.
   *
   * @ingroup h_refinement
   */
  virtual void rebuild_after_insert_knots(
    const SafeSTLArray<SafeSTLVector<Real>,dim> &knots_to_insert,
    const Grid<dim> &old_grid) override final;
#endif // MESH_REFINEMENT

};


//------------------------------------------------------------------------------
/**
 * F(x) = x
 */
template<int dim>
class IdentityGridFunction :
  public FormulaGridFunction<dim,dim>
{
  using base_t = GridFunction<dim,dim>;
  using parent_t = FormulaGridFunction<dim,dim>;
  using self_t = IdentityGridFunction<dim>;
  using typename base_t::GridType;
public:
  using typename parent_t::Value;
  using typename parent_t::GridPoint;
  template <int order>
  using Derivative = typename parent_t::template Derivative<order>;

public:
  static std::shared_ptr<self_t>
  create(const std::shared_ptr<GridType> &domain);

  static std::shared_ptr<const self_t>
  const_create(const std::shared_ptr<const GridType> &domain);


  IdentityGridFunction(const self_t &) = default;

  virtual ~IdentityGridFunction() = default;

  virtual void print_info(LogStream &out) const override final;

protected:
  IdentityGridFunction(
    const SharedPtrConstnessHandler<GridType> &domain);

private:
  void evaluate_0(const ValueVector<GridPoint> &points,
                  ValueVector<Value> &values) const override;

  void evaluate_1(const ValueVector<GridPoint> &points,
                  ValueVector<Derivative<1>> &values) const override;

  void evaluate_2(const ValueVector<GridPoint> &points,
                  ValueVector<Derivative<2>> &values) const override;


#ifdef MESH_REFINEMENT
  /**
   * Rebuild the internal state of the object after an Grid::insert_knots() function is invoked.
   *
   * @pre Before invoking this function, must be invoked the function grid_->insert_knots().
   * @note This function is connected to the Grid's signal for the refinement, and
   * it is necessary in order to avoid infinite loops in the insert_knots() function calls.
   *
   * @ingroup h_refinement
   */
  virtual void rebuild_after_insert_knots(
    const SafeSTLArray<SafeSTLVector<Real>,dim> &knots_to_insert,
    const Grid<dim> &old_grid) override final;
#endif // MESH_REFINEMENT

};



//------------------------------------------------------------------------------
/**
 * Maps a hyper rectangle into a spherical ball sector using the
 * dim-dimensional spherical coordinates, maps a hyper-rectangle
 * r in [0,R], phi_1 in [0, 2 pi], and phi_2, phi_dim-1 in [0,pi]
 * such that
 * x1 = r cos (phi_1)
 * x2 = r sin (phi_1) cos (phi_2)
 * etc
 *
 */
template<int dim>
class BallGridFunction : public FormulaGridFunction<dim, dim>
{
private:
  using base_t = GridFunction<dim, dim>;
  using parent_t =  FormulaGridFunction<dim, dim>;
  using self_t = BallGridFunction<dim>;
  using typename base_t::GridType;
public:
  using typename parent_t::Value;
  using typename parent_t::GridPoint;
  template <int order>
  using Derivative = typename parent_t::template Derivative<order>;

  static std::shared_ptr<self_t>
  create(std::shared_ptr<GridType> grid);

  static std::shared_ptr<const self_t>
  const_create(std::shared_ptr<const GridType> grid);

  BallGridFunction(const self_t &) = default;
  virtual ~BallGridFunction() = default;

  virtual void print_info(LogStream &out) const override final;


protected:
  BallGridFunction(const SharedPtrConstnessHandler<GridType> &grid);

private:
  template<int order>
  auto
  get_aux_vals(const ValueVector<GridPoint> &points) const;

private:
  void evaluate_0(const ValueVector<GridPoint> &points,
                  ValueVector<Value> &values) const override final;

  void evaluate_1(const ValueVector<GridPoint> &points,
                  ValueVector<Derivative<1>> &values) const override final;

  void evaluate_2(const ValueVector<GridPoint> &points,
                  ValueVector<Derivative<2>> &values) const override final;


#ifdef MESH_REFINEMENT
  /**
   * Rebuild the internal state of the object after an Grid::insert_knots() function is invoked.
   *
   * @pre Before invoking this function, must be invoked the function grid_->insert_knots().
   * @note This function is connected to the Grid's signal for the refinement, and
   * it is necessary in order to avoid infinite loops in the insert_knots() function calls.
   *
   * @ingroup h_refinement
   */
  virtual void rebuild_after_insert_knots(
    const SafeSTLArray<SafeSTLVector<Real>,dim> &knots_to_insert,
    const Grid<dim> &old_grid) override final;
#endif // MESH_REFINEMENT
};
//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
/**
 * Maps a hyper rectangle into a spherical sector using the
 * dim-dimensional spherical coordinates.
 * It maps a hyper-rectangle
 * r in [0,R], phi_1 in [0, 2 pi], and phi_2, phi_dim-1 in [0,pi]
 * such that
 * x1 = r cos (phi_1)
 * x2 = r sin (phi_1) cos (phi_2)
 * etc
 *
 */
template<int dim>
class SphereGridFunction : public FormulaGridFunction<dim, dim+1>
{
private:
  static const int space_dim = dim + 1;
  using base_t = GridFunction<dim, dim+1>;
  using parent_t = FormulaGridFunction<dim, dim+1>;
  using self_t = SphereGridFunction<dim>;
  using typename base_t::GridType;
public:
  using typename parent_t::GridPoint;
  using typename parent_t::Value;
  template <int order>
  using Derivative = typename parent_t::template Derivative<order>;

  static std::shared_ptr<self_t>
  create(std::shared_ptr<GridType> grid);

  static std::shared_ptr<const self_t>
  const_create(std::shared_ptr<const GridType> grid);

  SphereGridFunction(const self_t &) = default;
  virtual ~SphereGridFunction() = default;

  virtual void print_info(LogStream &out) const override final;

protected:
  SphereGridFunction(const SharedPtrConstnessHandler<GridType> &grid);

private:
  template<int order>
  auto
  get_aux_vals(const ValueVector<GridPoint> &points) const;

private:
  void evaluate_0(const ValueVector<GridPoint> &points,
                  ValueVector<Value> &values) const override final;

  void evaluate_1(const ValueVector<GridPoint> &points,
                  ValueVector<Derivative<1>> &values) const override final;

  void evaluate_2(const ValueVector<GridPoint> &points,
                  ValueVector<Derivative<2>> &values) const override final;

  static const int R = 1.;

#ifdef MESH_REFINEMENT
  /**
   * Rebuild the internal state of the object after an Grid::insert_knots() function is invoked.
   *
   * @pre Before invoking this function, must be invoked the function grid_->insert_knots().
   * @note This function is connected to the Grid's signal for the refinement, and
   * it is necessary in order to avoid infinite loops in the insert_knots() function calls.
   *
   * @ingroup h_refinement
   */
  virtual void rebuild_after_insert_knots(
    const SafeSTLArray<SafeSTLVector<Real>,dim> &knots_to_insert,
    const Grid<dim> &old_grid) override final;
#endif // MESH_REFINEMENT
};
//------------------------------------------------------------------------------




//------------------------------------------------------------------------------
/**
 * \brief This class represent a cylindrical annulus mapping.
 *
 * The mapping is
 * \f{equation*}{
 *    \begin{aligned}
 *    F(\hat{\theta},\hat{r},\hat{z}) \colon [0,1] \times[0,1]\times[0,1] & \to \Omega \\
 *    (\hat{\theta},\hat{r},\hat{z}) & \mapsto
 *    F(\hat{\theta},\hat{r},\hat{z}) =
 *    \begin{pmatrix}
 *      \bigl[ (r_1-r_0) \hat{r} + r_0 \bigr] \cos\bigl[ (\theta_1-\theta_0) \hat{\theta} \bigr] \\
 *      \bigl[ (r_1-r_0) \hat{r} + r_0 \bigr] \sin\bigl[ (\theta_1-\theta_0) \hat{\theta} \bigr] \\
 *      h_0 + (h_1-h_0) \hat{z}
 *    \end{pmatrix}
 *    \end{aligned}
 * \f}
 * where \f$ \Omega \f$ is a section of cylindrical annulus with the following characteristics:
 *   section angle \f$ \theta \in [\theta_0,\theta_1] \f$,
 *   radius \f$ r \in [r_0,r_1] \f$,
 *   height \f$ z \in[h_0,h_1] \f$.
 *
 * \author martinelli
 * \date 31 Jan 2013
 */
class CylindricalAnnulusGridFunction : public FormulaGridFunction<3,3>
{
private:
  using base_t = GridFunction<3,3>;
  using parent_t = FormulaGridFunction<3,3>;
  using self_t = CylindricalAnnulusGridFunction;
  using typename base_t::GridType;
public:
  using typename parent_t::GridPoint;
  using typename parent_t::Value;

  template <int order>
  using Derivative = typename parent_t::template Derivative<order>;

  static std::shared_ptr<self_t>
  create(std::shared_ptr<GridType> grid,
         const Real r0,
         const Real r1,
         const Real h0,
         const Real h1,
         const Real theta0,
         const Real theta1);

  static std::shared_ptr<const self_t>
  const_create(std::shared_ptr<const GridType> grid,
               const Real r0,
               const Real r1,
               const Real h0,
               const Real h1,
               const Real theta0,
               const Real theta1);


  CylindricalAnnulusGridFunction(const self_t &) = default;
  virtual ~CylindricalAnnulusGridFunction() = default;


  virtual void print_info(LogStream &out) const override final;

protected:
  CylindricalAnnulusGridFunction(const SharedPtrConstnessHandler<GridType> &grid,
                                 const Real r0,
                                 const Real r1,
                                 const Real h0,
                                 const Real h1,
                                 const Real theta0,
                                 const Real theta1);

private:
  void evaluate_0(const ValueVector<GridPoint> &points,
                  ValueVector<Value> &values) const override final;

  void evaluate_1(const ValueVector<GridPoint> &points,
                  ValueVector<Derivative<1>> &values) const override final;

  void evaluate_2(const ValueVector<GridPoint> &points,
                  ValueVector<Derivative<2>> &values) const override final;


  const Real r0_;
  const Real r1_;
  const Real h0_;
  const Real h1_;
  const Real theta0_;
  const Real theta1_;

  const Real dR_;
  const Real dT_;
  const Real dH_;

#ifdef MESH_REFINEMENT
  /**
   * Rebuild the internal state of the object after an Grid::insert_knots() function is invoked.
   *
   * @pre Before invoking this function, must be invoked the function grid_->insert_knots().
   * @note This function is connected to the Grid's signal for the refinement, and
   * it is necessary in order to avoid infinite loops in the insert_knots() function calls.
   *
   * @ingroup h_refinement
   */
  virtual void rebuild_after_insert_knots(
    const SafeSTLArray<SafeSTLVector<Real>,3> &knots_to_insert,
    const Grid<3> &old_grid) override final;
#endif // MESH_REFINEMENT
};

#endif
} // of namespace functions.

IGA_NAMESPACE_CLOSE


