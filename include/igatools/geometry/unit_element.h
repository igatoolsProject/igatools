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


#ifndef UNIT_ELEMENT_H_
#define UNIT_ELEMENT_H_

#include <igatools/base/config.h>
#include <igatools/base/tensor.h>


IGA_NAMESPACE_OPEN

/**
 * @brief This class provides dimension independent information of all topological
 * structures that make up the elements in the reference patch or knotspans.
 *
 */
template <int dim>
struct UnitElement
{
    /** Number of vertices of a element. */
    static const int vertices_per_element = 1 << dim;


    /**
     * Converts the local vertex index to the unit hypercube coordinates.
     * For example in dim==2 the hypercube is the square
     * and:
     * - vertex 0 has coordinates (0,0)
     * - vertex 1 has coordinates (1,0)
     * - vertex 2 has coordinates (0,1)
     * - vertex 3 has coordinates (1,1)
     */
    static const int vertex_to_component[vertices_per_element][dim];

    /** Dimension of the face. */
    static const int face_dim = (dim >= 1)?dim-1:0;

    /** Number of faces per element.*/
    static constexpr int faces_per_element = 2 * dim;

    static const std::array<int,faces_per_element> faces;
    /**
     * Converts the local face index of the unit element
     * to the hyperplane it belongs to.
     * More specifically it gives the constant coordinate
     * and its value.
     * For example in dim==2 the element is the unit square
     * and:
     * - face 0  is given by x=0 represented by {0,-1}
     * - face 1  is given by x=1 represented by {0, 1}, etc.
     *
     *todo: call it  face_to_plane or/and replace by face_constant_direction
     */
    static const int face_to_component[faces_per_element][2];


    /**
     * Given a constant direction, in dimension <em>dim</em>, there are <em>dim-1</em>
     * active directions.
     */
    static const
    Conditional< dim != 0, std::array<int,dim-1>, std::array<int,0> > active_directions[dim];


    /**
     * For each face id, in dimension <em>dim</em>, there are <em>dim-1</em>
     * active directions.
     */
    static const
    Conditional<dim != 0, TensorIndex<dim-1>, TensorIndex<0> >
    face_active_directions[faces_per_element];

    /** Direction along which the face coordinates are constant. */
    static const int face_constant_direction[faces_per_element];

    /** For each face gives the side (0 or 1). */
    static const int face_side[faces_per_element];


    /** Value of the constant coordinate locating the face. */
    static const Real face_constant_coordinate[faces_per_element];


    /**
     * For each face x_i=constant, the unit normal is (0,..0,+-1,0,...0)
     * here we store the sign (-1 or +1) corresponding to the outer
     * direction.
     */
    static const int face_normal_direction[faces_per_element];


    /** For each vertex, gives the opposite vertex index. */
    static const int opposite_vertex[vertices_per_element];


    /** Gives the outer boundary normal for every face. */
    static const Points<dim> face_normal[faces_per_element];
};

IGA_NAMESPACE_CLOSE

#endif /* UNIT_ELEMENT_H_ */
