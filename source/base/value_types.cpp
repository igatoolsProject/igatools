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


constexpr int _Value::id;
constexpr int _Value::order;
const string  _Value::name = "Value";


constexpr int _Gradient::id;
constexpr int _Gradient::order;
const string  _Gradient::name = "Gradient";


constexpr int _Hessian::id;
constexpr int _Hessian::order;
const string  _Hessian::name = "Hessian";


constexpr int _Divergence::id;
constexpr int _Divergence::order;
const string  _Divergence::name = "Divergence";


constexpr int _Point::id;
constexpr int _Point::order;
const string  _Point::name = "Point";

IGA_NAMESPACE_CLOSE
