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


namespace grid_element
{
//---------------------------------------------------------------------
const CacheFlags _Point::flag;
const string _Point::name = "Element Quadrature Points";

const CacheFlags _Weight::flag;
const string _Weight::name = "Element Quadrature Weights";

activate::FlagsToCache  activate::grid =
{
  {Flags::point, CacheFlags::point},
  {Flags::weight, CacheFlags::weight}
};
//---------------------------------------------------------------------
};


//---------------------------------------------------------------------
namespace domain_element
{

const CacheFlags _Measure::flag;
const string _Measure::name = "Element measure";

const CacheFlags _InvJacobian::flag;
const string _InvJacobian::name = "Element inverse jacobian";

const CacheFlags _InvHessian::flag;
const string _InvHessian::name = "Element inverse hessian";


const CacheFlags _BoundaryNormal::flag;
const string _BoundaryNormal::name = "Element boundary normal";


const CacheFlags _ExtNormal::flag;
const string _ExtNormal::name = "Element exterior normal";


activate::FlagsToCache  activate::domain =
{
  {Flags::point, CacheFlags::none},
  {Flags::w_measure, CacheFlags::measure},
  {Flags::measure,   CacheFlags::measure},
  {Flags::ext_normal, CacheFlags::ext_normal},
  {Flags::jacobian, CacheFlags::none},
  {Flags::inv_jacobian, CacheFlags::inv_jacobian},
  {Flags::hessian, CacheFlags::none},
  {Flags::inv_hessian, CacheFlags::inv_hessian},
  {Flags::boundary_normal, CacheFlags::boundary_normal | CacheFlags::inv_jacobian},
};

activate::FlagsToGridFunc activate::grid_func =
{
  {Flags::point, grid_function_element::Flags::D0},
  {Flags::w_measure, grid_function_element::Flags::D1},
  {Flags::measure, grid_function_element::Flags::D1},
  {Flags::ext_normal, grid_function_element::Flags::D1},
  {Flags::jacobian, grid_function_element::Flags::D1},
  {Flags::inv_jacobian, grid_function_element::Flags::D1},
  {Flags::hessian, grid_function_element::Flags::D2},
  {Flags::inv_hessian, grid_function_element::Flags::D2},
  {Flags::boundary_normal, grid_function_element::Flags::none}
};

activate::FlagsToGrid activate::grid =
{
  {Flags::point, grid_element::Flags::none},
  {Flags::w_measure, grid_element::Flags::weight},
  {Flags::measure, grid_element::Flags::none},
  {Flags::ext_normal, grid_element::Flags::none},
  {Flags::jacobian, grid_element::Flags::none},
  {Flags::inv_jacobian, grid_element::Flags::none},
  {Flags::hessian, grid_element::Flags::none},
  {Flags::inv_hessian, grid_element::Flags::none},
  {Flags::boundary_normal, grid_element::Flags::none}
};

};
//---------------------------------------------------------------------



//---------------------------------------------------------------------
namespace grid_function_element
{
//TODO: to be generated with mpl? or cog?
template<> const CacheFlags _D<0>::flag = CacheFlags::D0;
template<> const string _D<0>::name = "Grid Function D0";

template<> const CacheFlags _D<1>::flag = CacheFlags::D1;
template<> const string _D<1>::name = "Grid Function D1";

template<> const CacheFlags _D<2>::flag = CacheFlags::D2;
template<> const string _D<2>::name = "Grid Function D2";

template<> const CacheFlags _D<3>::flag = CacheFlags::D3;
template<> const string _D<3>::name = "Grid Function D3";


activate::FlagsToCache  activate::grid_function =
{
  {Flags::D0, CacheFlags::D0},
  {Flags::D1, CacheFlags::D1},
  {Flags::D2, CacheFlags::D2},
  {Flags::D3, CacheFlags::D3}
};

activate::FlagsToGrid activate::grid =
{
  {Flags::D0, grid_element::Flags::none},
  {Flags::D1, grid_element::Flags::none},
  {Flags::D2, grid_element::Flags::none},
  {Flags::D3, grid_element::Flags::none}
};

};
//---------------------------------------------------------------------



//---------------------------------------------------------------------
namespace function_element
{
const CacheFlags _Value::flag;
const string _Value::name = "Function Values";

const CacheFlags _Gradient::flag;
const string _Gradient::name = "Function gradients";

const CacheFlags _D2::flag;
const string _D2::name = "Function D2";

activate::FlagsToCache  activate::function =
{
  {Flags::value, CacheFlags::value},
  {Flags::gradient, CacheFlags::gradient},
  {Flags::D2, CacheFlags::D2}
};

function_element::activate::FlagsToDomain activate::domain =
{
  {Flags::value, domain_element::Flags::none},
  {Flags::gradient, domain_element::Flags::none},
  {Flags::D2, domain_element::Flags::none}
};

};
//---------------------------------------------------------------------



//---------------------------------------------------------------------
namespace space_element
{
const space_element::Flags space_element::_Value::flag;
const string space_element::_Value::name = "Basis function values";

const space_element::Flags space_element::_Gradient::flag;
const string space_element::_Gradient::name = "Basis function gradients";

const space_element::Flags space_element::_Hessian::flag;
const string space_element::_Hessian::name = "Basis function hessians";

const space_element::Flags space_element::_Divergence::flag;
const string space_element::_Divergence::name = "Basis function divergences";
};
//---------------------------------------------------------------------

IGA_NAMESPACE_CLOSE
