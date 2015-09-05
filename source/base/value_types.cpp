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

#include <igatools/base/value_types.h>

using std::string;

IGA_NAMESPACE_OPEN



//---------------------------------------------------------------------
const grid_element::Flags grid_element::_Point::flag;
const string grid_element::_Point::name = "Element Quadrature Points";

const grid_element::Flags grid_element::_Weight::flag;
const string grid_element::_Weight::name = "Element Quadrature Weights";
//---------------------------------------------------------------------



//---------------------------------------------------------------------
const domain_element::CacheFlags domain_element::_Point::flag;
const string domain_element::_Point::name = "Element Quadrature Points";

const domain_element::CacheFlags domain_element::_Measure::flag;
const string domain_element::_Measure::name = "Element measure";

const domain_element::CacheFlags domain_element::_Gradient::flag;
const string domain_element::_Gradient::name = "domain gradients";

domain_element::activate::FlagsToCache  domain_element::activate::domain =
{
  {domain_element::Flags::point, domain_element::CacheFlags::point},
  {
    domain_element::Flags::w_measure, domain_element::CacheFlags::gradient|
    domain_element::CacheFlags::measure
  },
  {
    domain_element::Flags::measure, domain_element::CacheFlags::gradient|
    domain_element::CacheFlags::measure
  }
};
domain_element::activate::FlagsToGrid domain_element::activate::grid =
{
  {domain_element::Flags::point, grid_element::Flags::point},
  {domain_element::Flags::w_measure, grid_element::Flags::weight},
  {domain_element::Flags::measure, grid_element::Flags::none}
};
//---------------------------------------------------------------------



//---------------------------------------------------------------------
const function_element::Flags function_element::_Value::flag;
const string function_element::_Value::name = "Function values";

const function_element::Flags function_element::_Gradient::flag;
const string function_element::_Gradient::name = "Function gradients";
//---------------------------------------------------------------------



//---------------------------------------------------------------------
const space_element::Flags space_element::_Value::flag;
const string space_element::_Value::name = "Basis function values";

const space_element::Flags space_element::_Gradient::flag;
const string space_element::_Gradient::name = "Basis function gradients";

const space_element::Flags space_element::_Hessian::flag;
const string space_element::_Hessian::name = "Basis function hessians";

const space_element::Flags space_element::_Divergence::flag;
const string space_element::_Divergence::name = "Basis function divergences";
//---------------------------------------------------------------------


IGA_NAMESPACE_CLOSE
