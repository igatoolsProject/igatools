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

#ifndef EQUALITY_CONSTRAINT_H_
#define EQUALITY_CONSTRAINT_H_

#include <igatools/base/config.h>
#include <igatools/base/logstream.h>

IGA_NAMESPACE_OPEN

/**
 * @brief Equality constraint between two degrees of freedom.
 *
 * @author M. Martinelli, 2014
 */
class EqualityConstraint
{
public:

  /** @name Constructor and destructor. */
  ///@{
  /**
   * Default constructor. Not allowed to be used.
   */
  EqualityConstraint() = delete;

  /**
   * Builds the equality constraint between the two degrees of freedom specified by
   * @p dof_id_master and @p dof_id_slave.
   *
   * @warning The two degrees of freedom must be different and <tt> >= 0</tt>,
   * otherwise in DEBUG mode an assertion will be raised.
   */
  EqualityConstraint(const Index dof_id_master,const Index dof_id_slave);

  /**
   * Copy constructor.
   */
  EqualityConstraint(const EqualityConstraint &constraint) = default;

  /**
   * Move constructor.
   */
  EqualityConstraint(EqualityConstraint &&constraint) = default;

  /**
   * Destructor.
   */
  ~EqualityConstraint() = default;
  ///@}


  /**
   * Returns the id of the master dof.
   */
  Index get_dof_id_master() const;

  /**
   * Returns the id of the slave dof.
   */
  Index get_dof_id_slave() const;

  /**
   * Prints the dofs id defining the equality constraint.
   */
  void print_info(LogStream &out) const;

private:

  /**
   * Id of the master dof.
   */
  Index dof_id_master_;

  /**
   * Id of the slave dof.
   */
  Index dof_id_slave_;
};


IGA_NAMESPACE_CLOSE

#endif //#ifndef EQUALITY_CONSTRAINT_H_
