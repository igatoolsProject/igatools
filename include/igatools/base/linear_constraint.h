//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2016  by the igatools authors (see authors.txt).
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

#ifndef LINEAR_CONSTRAINT_H
#define LINEAR_CONSTRAINT_H

#include <igatools/base/config.h>
#include <igatools/utils/safe_stl_vector.h>

IGA_NAMESPACE_OPEN


/**
 * @brief This class represent a linear constraint, i.e. a relation of the type:
 * \f$ c_{i_1} \mathbf{x}_{i_1} + \dots + c_{i_n} \mathbf{x}_{i_n} = f \f$
 * where the indices \f$ (i_1,\dots,i_n) \f$ are the dofs indices,
 * \f$ (c_{i_1},\dots,c_{i_n}) \in \mathbb{R}^n \f$ are the coefficients and \f$ f \in \mathbb{R} \f$
 * is the right hand side that defines the linear constraint.
 *
 * To the LinearConstraint is associated a LinearConstraintType member variable
 * that indicates the type of the constraint.
 *
 * @author M. Martinelli
 * @date 2014
 */
class LinearConstraint
{
public:

  /**
   * It builds and returns a std::shared_ptr wrapping a
   * LinearConstraint (associated to the @p global_dof_id) of the specified @type,
   * from a SafeSTLVector of @p dofs defining the constraint,
   * their relative @p coefficients and the tight hand side value @p rhs.
   */
  static std::shared_ptr<LinearConstraint>
  create(const Index global_dof_id,
         const LinearConstraintType &type,
         const SafeSTLVector<Index> &dofs,const SafeSTLVector<Real> &coeffs,const Real rhs);

  /** Destructor. */
  ~LinearConstraint() = default;


  /** Returns the value of the right hand side of the LinearConstraint. */
  Real get_rhs() const;


  /** Sets the value of the right hand side of the LinearConstraint. */
  void set_rhs(const Real rhs);

  /**
   * Returns the number of terms involved in the left hand side of the linear constraint.
   * This number is also the number of dofs (and therefore coefficients) involved in the linear constraint.
   */
  Index get_num_lhs_terms() const;

  /**
   * Returns the number of dofs involved in the left hand side of the linear constraint.
   */
  Index get_num_dofs() const;

  /**
   * Returns the number of coefficients involved in the left hand side of the linear constraint.
   */
  Index get_num_coeffs() const;

  /**
   * Return the index of the <tt>p</tt>-th dof involved to define the linear constraint.
   */
  Index get_dof_index(const int i) const;

  /**
   * Return the coefficient associated to the the <tt>p</tt>-th dof involved to define the linear constraint.
   */
  Real get_coeff(const int i) const;


  /**
   * Returns the dofs involved in the definition of the LinearConstraint.
   */
  SafeSTLVector<Index> get_dofs_id() const;

  /**
   * Returns the coefficients of the dofs involved in the definition of the LinearConstraint.
   */
  SafeSTLVector<Real> get_coefficients() const;

  /**
   * Returns the type of the LinearConstraint.
   */
  LinearConstraintType get_type() const;

  /**
   * Returns the global dof id associated to the linear constraint.
   */
  Index get_global_dof_id() const;

  /**
   * @name Printing info
   */
  ///@{

  /**
   * Prints the content of the LinearConstraint on the LogStream @p out.
   */
  void print_info(LogStream &out) const ;
  ///@}


  /**
   * Returns TRUE if the @p dof is present in the definition of the LinearConstraint.
   */
  bool is_dof_present(const Index dof) const;


  /**
   * Returns the absolute error of the constraint if applied to the global dof coefficients @p
   * dof_coeffs.
   */
  Real eval_absolute_error(const SafeSTLVector<Real> &dof_coeffs) const;



private:
  /** @name Constructors ad destructor */
  ///@{
  /** Default constructor. Not allowed to be used. */
  LinearConstraint() = delete;

  /**
   * Constructor. It builds the LinearConstraint (associated to the @p global_dof_id)
   * of the specified @type, from a SafeSTLVector of @p dofs defining the constraint,
   * their relative @p coefficients and the tight hand side value @p rhs.
   *
   * @note Not allowed to be used. Use the associated create() function.
   */
  LinearConstraint(const Index global_dof_id,
                   const LinearConstraintType &type,
                   const SafeSTLVector<Index> &dofs,const SafeSTLVector<Real> &coeffs,const Real rhs);

  /**
   * Copy constructor. Not allowed to be used.
   */
  LinearConstraint(const LinearConstraint &lc) = delete;


  /**
   * Move constructor. Not allowed to be used.
   */
  LinearConstraint(LinearConstraint &&lc) = delete;


  /** Copy assignment operator. Not allowed to be used. */
  LinearConstraint &operator=(const LinearConstraint &lc) = delete;


  /** Move assignment operator. Not allowed to be used. */
  LinearConstraint &operator=(LinearConstraint &&lc) = delete;
  ///@}



  /** Returns a const reference of the right hand side of the linear constraint. */
  const std::pair<Index,Real> &get_lhs_term(const int i) const;

  /**
   * Left hand side of the linear equation.
   *
   * Each entry of this SafeSTLVector is a pair (global dof id/coefficient).
   */
  SafeSTLVector<std::pair<Index,Real>> lhs_;

  /**
   * Right hand side of the linear equation.
   */
  Real rhs_;

  /** Type of the linear constraint. */
  LinearConstraintType type_;


  /** Global dof id associated to the linear constraint. */
  Index global_dof_id_;
};


IGA_NAMESPACE_CLOSE



#endif //#ifndef LINEAR_CONSTRAINT_H

