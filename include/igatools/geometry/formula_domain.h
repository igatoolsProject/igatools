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

#ifndef __FORMULA_DOMAIN_H_
#define __FORMULA_DOMAIN_H_

#include <igatools/geometry/physical_domain.h>
#include <igatools/geometry/physical_domain_element.h>
#include <igatools/base/value_types.h>

IGA_NAMESPACE_OPEN

/**
 *
 */
template<int dim, int codim>
class FormulaDomain :
  public PhysicalDomain<dim, codim>
{
private:
  using parent_t =  PhysicalDomain<dim, codim>;
  using self_t = FormulaDomain<dim, codim>;
protected:
  using typename parent_t::GridType;
public:
  using typename parent_t::Point;
  using typename parent_t::GridPoint;

  template <int order>
  using Derivative = typename parent_t::template Derivative<order>;

  FormulaDomain(std::shared_ptr<GridType> grid);

  virtual ~FormulaDomain() = default;

public:

  virtual void evaluate_0(const ValueVector<GridPoint> &points,
                          ValueVector<Point> &values) const = 0;

  virtual void evaluate_1(const ValueVector<GridPoint> &points,
                          ValueVector<Derivative<1>> &values) const = 0;

  virtual void evaluate_2(const ValueVector<GridPoint> &points,
                          ValueVector<Derivative<2>> &values) const = 0;

};

IGA_NAMESPACE_CLOSE

#endif
