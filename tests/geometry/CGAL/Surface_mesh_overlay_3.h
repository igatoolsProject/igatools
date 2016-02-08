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
// Copyright (c) 2014  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Sebastien Loriot <sebastien.loriot@cgal.org>
//

#ifndef CGAL_SURFACE_MESH_OVERLAY_3_H
#define CGAL_SURFACE_MESH_OVERLAY_3_H

/// \todo Maybe we should have an overlay of 2D patches to prevent points to be projected
///       in another region if we have a segmentation of the surface. I have in mind what happen
///       to the elephant near a saddle point (at the junction of the trump and the tusks).
/// \todo See how to remove the vertex index property map by changing the point ppmap to directly
///       using the vertex instead of its index, and by having a normal vertex ppmap

//#include <CGAL/Surface_mesh_overlay_3/Surface_mesh_overlay_3_impl.h>
#include "Surface_mesh_overlay_3/Surface_mesh_overlay_3_impl.h"

namespace CGAL
{

template <class Input_kernel, class Exact_kernel, class Polyhedron>
void surface_mesh_overlay_3(
  Polyhedron &P1,
  Polyhedron &P2,
  double min_half_diameter,
  std::vector< std::pair<std::size_t, std::size_t> > &input_facet_ids,
  bool allow_partial_overlay = false
)
{
  Surface_mesh_overlay_3::Surface_mesh_overlay_3_impl<
  Input_kernel,Exact_kernel,Polyhedron>                   algorithm;
  algorithm.run(P1, P2, min_half_diameter, allow_partial_overlay, input_facet_ids);
}

template <class Input_kernel, class Exact_kernel, class Polyhedron>
void surface_mesh_overlay_3(
  Polyhedron &P1,
  Polyhedron &P2,
  double min_half_diameter,
  const std::list<
  std::pair< std::vector<typename Polyhedron::Halfedge_handle>,
  std::vector<typename Polyhedron::Halfedge_handle> >
  >& matching_features,
  std::vector< std::pair<std::size_t, std::size_t> > &input_facet_ids,
  bool allow_partial_overlay = false
)
{
  Surface_mesh_overlay_3::Surface_mesh_overlay_3_impl<
  Input_kernel,Exact_kernel,Polyhedron>                   algorithm;
  algorithm.run(P1, P2, min_half_diameter, matching_features, allow_partial_overlay, input_facet_ids);
}

template <class EdgeSkipper, class Input_kernel, class Exact_kernel, class Polyhedron>
void surface_mesh_overlay_with_edge_skipper_3(
  Polyhedron &P1,
  Polyhedron &P2,
  double min_half_diameter,
  std::vector< std::pair<std::size_t, std::size_t> > &input_facet_ids,
  bool allow_partial_overlay = false
)
{
  Surface_mesh_overlay_3::Surface_mesh_overlay_3_impl<
  Input_kernel,Exact_kernel,Polyhedron, EdgeSkipper>   algorithm;
  algorithm.run(P1, P2, min_half_diameter, allow_partial_overlay, input_facet_ids);
}

template <class EdgeSkipper, class Input_kernel, class Exact_kernel, class Polyhedron>
void surface_mesh_overlay_with_edge_skipper_3(
  Polyhedron &P1,
  Polyhedron &P2,
  double min_half_diameter,
  const std::list<
  std::pair< std::vector<typename Polyhedron::Halfedge_handle>,
  std::vector<typename Polyhedron::Halfedge_handle> >
  >& matching_features,
  std::vector< std::pair<std::size_t, std::size_t> > &input_facet_ids,
  bool allow_partial_overlay = false
)
{
  Surface_mesh_overlay_3::Surface_mesh_overlay_3_impl<
  Input_kernel,Exact_kernel,Polyhedron, EdgeSkipper>   algorithm;
  algorithm.run(P1, P2, min_half_diameter, matching_features, allow_partial_overlay, input_facet_ids);
}

} // end of namespace CGAL

#endif // CGAL_SURFACE_MESH_OVERLAY_3_H
