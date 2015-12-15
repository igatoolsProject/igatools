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

#ifndef CGAL_SURFACE_MESH_OVERLAY_3_POLYHEDRON_UTILS_H
#define CGAL_SURFACE_MESH_OVERLAY_3_POLYHEDRON_UTILS_H

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/HalfedgeDS_face_max_base_with_id.h>
#include <CGAL/HalfedgeDS_vertex_max_base_with_id.h>
#include <CGAL/HalfedgeDS_halfedge_max_base_with_id.h>
#include "../Polyhedron_simplex_type.h"
#include <CGAL/array.h>

#include<set>

namespace CGAL
{

namespace Surface_mesh_overlay_3
{

template<class Refs>
struct Facet_with_input_id :
  CGAL::HalfedgeDS_face_max_base_with_id< Refs, CGAL::Tag_false, std::size_t>
{
  typedef CGAL::HalfedgeDS_face_max_base_with_id
  < Refs, CGAL::Tag_false, std::size_t> Base;

  Facet_with_input_id()
    : Base(),
      m_input_id((std::numeric_limits<std::size_t>::max)()) {}

  std::size_t m_input_id;
  std::size_t &input_id()
  {
    return m_input_id;
  }
  std::size_t  input_id() const
  {
    return m_input_id;
  }
};

class Poly_items_with_id_3
{
public:
  template < class Refs, class Traits>
  struct Vertex_wrapper
  {
    typedef typename Traits::Point_3 Point;
    typedef CGAL::HalfedgeDS_vertex_max_base_with_id
    < Refs, Point, std::size_t> Vertex;
  };
  template < class Refs, class Traits>
  struct Halfedge_wrapper
  {
    typedef CGAL::HalfedgeDS_halfedge_max_base_with_id
    <Refs, std::size_t> Halfedge;
  };
  template < class Refs, class Traits>
  struct Face_wrapper
  {
    typedef Facet_with_input_id< Refs >  Face;
  };
};

/// utility functions
// checks if the line (p1,p2) intersects `f` in its interior, on an edge or a vertex
template <class Polyhedron, class Exact_kernel>
std::pair<typename Polyhedron::Halfedge_handle, Polyhedron_simplex_type>
intersection_type(typename Polyhedron::Face_handle f,
                  const typename Exact_kernel::Point_3 &p1,
                  const typename Exact_kernel::Point_3 &p2)
{
  typedef std::pair<typename Polyhedron::Halfedge_handle, Polyhedron_simplex_type> result_type;
  CGAL::Cartesian_converter<typename Polyhedron::Traits, Exact_kernel> to_exact;
  cpp11::array<typename Exact_kernel::Point_3, 3> pts;
  typename Polyhedron::Halfedge_handle h=f->halfedge();
  for (int i=0; i<3; ++i)
  {
    pts[i]=to_exact(h->opposite()->vertex()->point());
    h=h->next();
  }
  // indicates if h, h->next() , h->prev and p1,p2 are coplanar
  cpp11::array<bool, 3> is_coplanar;
  int nb_coplanar=0;
  for (int i=0; i<3; ++i)
    if ((is_coplanar[i]=coplanar(pts[i], pts[(i+1)%3], p1, p2)))
      ++nb_coplanar;
  CGAL_assertion(nb_coplanar<3);

  switch (nb_coplanar)
  {
    case 0:
      return result_type(h, POLYHEDRON_FACET);
    case 1:
      for (int i=0; i<3; ++i)
      {
        if (is_coplanar[i])
          return result_type(h, POLYHEDRON_EDGE);
        h=h->next();
      }
    case 2:
      for (int i=0; i<3; ++i)
      {
        if (!is_coplanar[i])
          return result_type(h->next(), POLYHEDRON_VERTEX);
        h=h->next();
      }
    default:
      CGAL_assertion(!"We should never be here");
      return result_type(h, POLYHEDRON_NONE);
  }
}

template <class Halfedge_handle>
Halfedge_handle
canonical_edge(Halfedge_handle h)
{
  return h<h->opposite()?h:h->opposite();
}

template <class Polyhedron, class VertexIndexPropertyMap>
void set_halfedgeds_items_id(
  Polyhedron &P,
  std::vector<typename Polyhedron::Facet_handle> &facets,
  std::vector<typename Polyhedron::Halfedge_handle> &edges,
  std::vector<typename Polyhedron::Vertex_handle> &vertices,
  VertexIndexPropertyMap vipmap)
{
  std::size_t vertex_id   = 0 ;
  std::size_t edge_id = 0 ;
  std::size_t face_id     = 0 ;

  for (typename Polyhedron::Vertex_iterator vit = P.vertices_begin(),
       evit = P.vertices_end()
              ; vit != evit ; ++  vit
      )
    put(vipmap, vit, vertex_id ++);

  for (typename Polyhedron::Edge_iterator hit = P.edges_begin(),
       ehit = P.edges_end()
              ; hit != ehit ; ++  hit
      )
  {
    hit->id() = edge_id ++ ;
    hit->opposite()->id() = edge_id ;
  }

  for (typename Polyhedron::Facet_iterator fit = P.facets_begin(),
       efit = P.facets_end()
              ; fit != efit ; ++  fit
      )
    fit->id() = face_id ++ ;

  facets.reserve(face_id);
  vertices.reserve(vertex_id);
  edges.reserve(edge_id);

  for (typename Polyhedron::Vertex_iterator vit = P.vertices_begin(),
       evit = P.vertices_end()
              ; vit != evit ; ++  vit
      )
  {
    vertices.push_back(vit);
  }

  for (typename Polyhedron::Edge_iterator hit = P.edges_begin(),
       ehit = P.edges_end()
              ; hit != ehit ; ++  hit
      )
    edges.push_back(canonical_edge(hit));

  for (typename Polyhedron::Facet_iterator fit = P.facets_begin(),
       efit = P.facets_end()
              ; fit != efit ; ++  fit
      )
    facets.push_back(fit);
}

template <class Halfedge_handle>
bool nested_edge_same_orientation(Halfedge_handle includer, Halfedge_handle included)
{
  CGAL_precondition(
    CGAL::collinear_are_ordered_along_line(includer->vertex()->point(),
                                           included->vertex()->point(),
                                           includer->opposite()->vertex()->point()) &&
    CGAL::collinear_are_ordered_along_line(includer->vertex()->point(),
                                           included->opposite()->vertex()->point(),
                                           includer->opposite()->vertex()->point())
  );

  return CGAL::collinear_are_ordered_along_line(
           included->opposite()->vertex()->point(),
           included->vertex()->point(),
           includer->vertex()->point());
}


} //end of namespace Surface_mesh_overlay_3

} //end of namespace CGAL

#endif //CGAL_SURFACE_MESH_OVERLAY_3_POLYHEDRON_UTILS_H
