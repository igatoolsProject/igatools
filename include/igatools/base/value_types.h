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

#ifndef VALUE_TYPES_H_
#define VALUE_TYPES_H_

#include <igatools/base/config.h>

#include <string>

IGA_NAMESPACE_OPEN


template <int value_type_id>
class ValueType;

template<>
class ValueType<0>
{
public:
    static constexpr int type_id = 0;
    static constexpr int order = 0;
    static const std::string name;
};
using _Value = ValueType<0>;


template<>
class ValueType<1>
{
public:
    static constexpr int type_id = 1;
    static constexpr int order = 1;
    static const std::string name;
};
using _Gradient = ValueType<1>;

template<>
class ValueType<2>
{
public:
    static constexpr int type_id = 2;
    static constexpr int order = 2;
    static const std::string name;
};
using _Hessian = ValueType<2>;


template<>
class ValueType<-1>
{
public:
    static constexpr int type_id = -1;
    static constexpr int order = 1;
    static const std::string name;
};
using _Divergence = ValueType<-1>;



IGA_NAMESPACE_CLOSE

#endif //#ifndef  VALUE_TYPES_H_
