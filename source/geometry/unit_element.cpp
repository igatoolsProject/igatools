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

#include <igatools/geometry/unit_element.h>
#include <igatools/base/array_utils.h>

IGA_NAMESPACE_OPEN

template <>
const int
UnitElement<0>::vertex_to_component[vertices_per_element][0]
=
{{}};

template <>
const int
UnitElement<1>::vertex_to_component[vertices_per_element][1]
=
{
    { 0 },
    { 1 }
};

template <>
const int
UnitElement<2>::vertex_to_component[vertices_per_element][2]
=
{
    { 0, 0 },
    { 1, 0 },
    { 0, 1 },
    { 1, 1 }
};


template <>
const int
UnitElement<3>::vertex_to_component[vertices_per_element][3]
=
{
    { 0, 0, 0 },
    { 1, 0, 0 },
    { 0, 1, 0 },
    { 1, 1, 0 },
    { 0, 0, 1 },
    { 1, 0, 1 },
    { 0, 1, 1 },
    { 1, 1, 1 }
};

template <>
const int
UnitElement<4>::vertex_to_component[vertices_per_element][4]
=
{
    { 0, 0, 0, 0},
    { 1, 0, 0, 0},
    { 0, 1, 0, 0},
    { 1, 1, 0, 0},
    { 0, 0, 1, 0},
    { 1, 0, 1, 0},
    { 0, 1, 1, 0},
    { 1, 1, 1, 0},
    { 0, 0, 0, 1},
    { 1, 0, 0, 1},
    { 0, 1, 0, 1},
    { 1, 1, 0, 1},
    { 0, 0, 1, 1},
    { 1, 0, 1, 1},
    { 0, 1, 1, 1},
    { 1, 1, 1, 1}
};



template <int dim>
const std::array<int, UnitElement<dim>::faces_per_element>
UnitElement<dim>::faces = sequence<UnitElement<dim>::faces_per_element>();


template <>
const int
UnitElement<0>::face_to_component[faces_per_element][2]
    = {};


template <>
const int
UnitElement<1>::face_to_component[faces_per_element][2]
=
{
    {0, -1} ,
    {0,  1}
};


template <>
const int
UnitElement<2>::face_to_component[faces_per_element][2]
=
{
    {0, -1} ,
    {0, 1} ,
    {1, -1} ,
    {1, 1}
};


template <>
const int
UnitElement<3>::face_to_component[faces_per_element][2]
=
{
    {0, -1} ,
    {0, 1} ,
    {1, -1} ,
    {1, 1} ,
    {2, -1} ,
    {2, 1}
};

template <>
const int
UnitElement<4>::face_to_component[faces_per_element][2]
=
{
    {0, -1} ,
    {0, 1} ,
    {1, -1} ,
    {1, 1} ,
    {2, -1} ,
    {2, 1},
    {3, -1} ,
    {3, 1}
};

template <>
const TensorIndex<0>
UnitElement<0>::face_active_directions[faces_per_element]
    = {};

template <>
const TensorIndex<0>
UnitElement<1>::face_active_directions[faces_per_element]
= {{},{}};

template <>
const TensorIndex<1>
UnitElement<2>::face_active_directions[faces_per_element]=
{{{1}},{{1}},{{0}},{{0}}};

template <>
const TensorIndex<2>
UnitElement<3>::face_active_directions[faces_per_element]=
{
    {{1, 2}} ,
    {{1, 2}} ,
    {{0, 2}} ,
    {{0, 2}} ,
    {{0, 1}} ,
    {{0, 1}}
};

template <>
const TensorIndex<3>
UnitElement<4>::face_active_directions[faces_per_element]=
{
    {{1, 2, 3}} ,
    {{1, 2, 3}} ,
    {{0, 2, 3}} ,
    {{0, 2, 3}} ,
    {{0, 1, 3}} ,
    {{0, 1, 3}} ,
    {{0, 1, 2}} ,
    {{0, 1, 2}}
};


template <>
const int_array<0>
UnitElement<1>::active_directions[1]
= {{}};


template <>
const int_array<1>
UnitElement<2>::active_directions[2]=
{{{1}},{{0}}};


template <>
const int_array<2>
UnitElement<3>::active_directions[3]=
{
    {{1, 2}} ,
    {{0, 2}} ,
    {{0, 1}}
};

template <>
const int_array<3>
UnitElement<4>::active_directions[4]=
{
    {{1, 2, 3}} ,
    {{0, 2, 3}} ,
    {{0, 1, 3}} ,
    {{0, 1, 2}}
};


template <>
const int
UnitElement<0>::face_side[faces_per_element]=
    {};

template <>
const int
UnitElement<1>::face_side[faces_per_element]=
{0,1};


template <>
const int
UnitElement<2>::face_side[faces_per_element]=
{0,1,0,1};

template <>
const int
UnitElement<3>::face_side[faces_per_element]=
{0,1,0,1,0,1};

template <>
const int
UnitElement<4>::face_side[faces_per_element]=
{0,1,0,1,0,1,0,1};

template <>
const int
UnitElement<0>::face_constant_direction[faces_per_element]=
    {};

template <>
const int
UnitElement<1>::face_constant_direction[faces_per_element]=
{0,0};

template <>
const int
UnitElement<2>::face_constant_direction[faces_per_element]=
{0,0,1,1};

template <>
const int
UnitElement<3>::face_constant_direction[faces_per_element]=
{0,0,1,1,2,2};

template <>
const int
UnitElement<4>::face_constant_direction[faces_per_element]=
{0,0,1,1,2,2,3,3};

template <>
const Real
UnitElement<0>::face_constant_coordinate[faces_per_element]=
    {};

template <>
const Real
UnitElement<1>::face_constant_coordinate[faces_per_element]=
{0.,1.};

template <>
const Real
UnitElement<2>::face_constant_coordinate[faces_per_element]=
{0.,1.,0.,1.};


template <>
const Real
UnitElement<3>::face_constant_coordinate[faces_per_element]=
{0.,1.,0.,1.,0.,1.};

template <>
const Real
UnitElement<4>::face_constant_coordinate[faces_per_element]=
{0.,1.,0.,1.,0.,1.,0.,1.};

template <>
const int
UnitElement<0>::face_normal_direction[faces_per_element]
    = {};

template <>
const int
UnitElement<1>::face_normal_direction[faces_per_element]
    = {-1, 1};

template <>
const int
UnitElement<2>::face_normal_direction[faces_per_element]
    = {-1, 1, -1, 1};

template <>
const int
UnitElement<3>::face_normal_direction[faces_per_element]
    = {-1, 1, -1, 1, -1, 1};

template <>
const int
UnitElement<4>::face_normal_direction[faces_per_element]
    = {-1, 1, -1, 1, -1, 1, -1, 1};

template <>
const int
UnitElement<0>::opposite_vertex[vertices_per_element]
    = {0};

template <>
const int
UnitElement<1>::opposite_vertex[vertices_per_element]
    = {1,0};


template <>
const int
UnitElement<2>::opposite_vertex[vertices_per_element]
    = {3,2,1,0};


template <>
const int
UnitElement<3>::opposite_vertex[vertices_per_element]
    = {7,6,5,4,3,2,1,0};

template <>
const int
UnitElement<4>::opposite_vertex[vertices_per_element]
    = {15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0};

template <>
const Points<0>
UnitElement<0>::face_normal[faces_per_element]=
    {};


template <>
const Points<1>
UnitElement<1>::face_normal[faces_per_element]=
{
    Points<1>({Real(-1.0)}),
    Points<1>({Real(1.0)})
};

template <>
const Points<2>
UnitElement<2>::face_normal[faces_per_element]=
{
    Points<2>({Real(-1.0),Real(0.0)}),
    Points<2>({Real(1.0),Real(0.0)}),
    Points<2>({Real(0.0),Real(-1.0)}),
    Points<2>({Real(0.0),Real(1.0)})
};

template <>
const Points<3>
UnitElement<3>::face_normal[faces_per_element]=
{
    Points<3>({Real(-1.0),Real(0.0),Real(0.0)}),
    Points<3>({Real(1.0),Real(0.0),Real(0.0)}),
    Points<3>({Real(0.0),Real(-1.0),Real(0.0)}),
    Points<3>({Real(0.0),Real(1.0),Real(0.0)}),
    Points<3>({Real(0.0),Real(0.0),Real(-1.0)}),
    Points<3>({Real(0.0),Real(0.0),Real(1.0)})
};

template <>
const Points<4>
UnitElement<4>::face_normal[faces_per_element]=
{
    Points<4>({Real(-1.0),Real(0.0),Real(0.0),Real(0.0)}),
    Points<4>({Real(1.0),Real(0.0),Real(0.0),Real(0.0)}),
    Points<4>({Real(0.0),Real(-1.0),Real(0.0),Real(0.0)}),
    Points<4>({Real(0.0),Real(1.0),Real(0.0),Real(0.0)}),
    Points<4>({Real(0.0),Real(0.0),Real(-1.0),Real(0.0)}),
    Points<4>({Real(0.0),Real(0.0),Real(1.0),Real(0.0)}),
    Points<4>({Real(0.0),Real(0.0),Real(0.0),Real(1.0)}),
    Points<4>({Real(0.0),Real(0.0),Real(0.0),Real(-1.0)})
};


IGA_NAMESPACE_CLOSE


#include <igatools/geometry/unit_element.inst>

