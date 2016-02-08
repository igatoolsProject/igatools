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
// Copyright (c) 2014 GeometryFactory (France).
// All rights reserved.
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
//
// Author(s)     : Sébastien Loriot

#ifndef CGAL_WALK_ON_POLYHEDRON_VISITOR_H
#define CGAL_WALK_ON_POLYHEDRON_VISITOR_H

namespace CGAL
{
namespace Surface_mesh_overlay_3
{

template <class Polyhedron, class Modifier, class Root_of_2, class VertexIndexPropertyMap>
class Walk_on_polyhedron_visitor
{
  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
  typedef typename Polyhedron::Vertex_handle   Vertex_handle;
  typedef typename Polyhedron::Facet_handle    Facet_handle;

  /// `includer` and `included` are two halfedges from two different polyhedron.
  /// the edge of `includer` fully contains in its interior `included`. This
  /// predicates indicate if the halfedges point in the same direction
  bool nested_edge_same_orientation(Halfedge_handle includer,
                                    Halfedge_handle included) const
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

//data members
  std::size_t m_source_index;
  std::size_t m_target_index;
  Halfedge_handle m_source; // halfedge of source (used only when source and target are on an edge)
  Halfedge_handle m_target;
  Polyhedron_simplex_type m_tgt_inter_type;
  Halfedge_handle m_hp2;
  Modifier &m_modifier;
  VertexIndexPropertyMap m_vertex_ipmap;
  Root_of_2 m_last_beta;   // barycentric coordinate of the last intersection
  // point for the edge in P2, m_hp2

public:
  Walk_on_polyhedron_visitor(std::size_t source_index,
                             std::size_t target_index,
                             Halfedge_handle h_src,
                             Halfedge_handle h_tgt,
                             Polyhedron_simplex_type tgt_inter_type,
                             Halfedge_handle hp2,
                             Modifier &modifier,
                             VertexIndexPropertyMap vipmap) :
    m_source_index(source_index), m_target_index(target_index),
    m_source(h_src), m_target(h_tgt), m_tgt_inter_type(tgt_inter_type),
    m_hp2(hp2), m_modifier(modifier), m_vertex_ipmap(vipmap)
  {}

  // intersected_edge points in the facet into which the edge should be added
  // Warning inside_edge_info must have been computed using
  // `canonical_edge(intersected_edge)`
  void found_edge_intersection(Halfedge_handle intersected_edge,
                               const Root_of_2 &inside_edge_info)
  {
    std::size_t newpt_index =
      m_modifier.add_point_on_edge(canonical_edge(intersected_edge), inside_edge_info);
    //insert the current edge
    m_modifier.add_edge(m_source_index,newpt_index, intersected_edge->facet(), m_hp2);
    m_source_index=newpt_index;
  }

  void found_vertex_intersection(Halfedge_handle h)
  {
    m_modifier.add_edge(m_source_index, get(m_vertex_ipmap, h->vertex()), h->facet(), m_hp2);
    m_source_index=get(m_vertex_ipmap, h->vertex());
  }

  void on_walk_end(Facet_handle facet)
  {
    m_modifier.add_edge(
      m_source_index,
      m_target_index,
      facet,
      m_hp2);
  }

  /// indicates that hp1 is on m_hp2, with the same orientation.
  /// if hp1 is not initialized, it indicates that m_source and m_target
  /// are both in the interior of the edge and a geometric predicates
  /// need to be used to get a correct hp1
  void on_input_edge(Halfedge_handle hp1)
  {
    if (hp1==Halfedge_handle())
    {
      hp1 = nested_edge_same_orientation(m_source, m_hp2)?
            m_source:m_source->opposite();
      CGAL_precondition(nested_edge_same_orientation(hp1, m_hp2));
    }

    std::size_t h_src_index = get(m_vertex_ipmap, hp1->opposite()->vertex());
    std::size_t h_tgt_index = get(m_vertex_ipmap, hp1->vertex());
    bool walk_ends_in_hp1_interior =
      (m_tgt_inter_type==POLYHEDRON_EDGE &&
       (m_target==hp1 || m_target==hp1->opposite()));

    /// update the edge incidence only if hp1 is fully contained in the projected version of m_hp2
    if (m_source_index==h_src_index   // entry point of the walk is the source of the edge
        &&
        (m_target_index==h_tgt_index ||  // walks ends at hp1 target
         !walk_ends_in_hp1_interior)
       )
      m_modifier.update_edge_incidence(hp1, m_hp2);
    else
    {
      /// the projection of m_hp2 overlaps hp1 (not all of it)
      m_modifier.record_edge_overlap(hp1,
                                     m_source_index,
                                     walk_ends_in_hp1_interior ? m_target_index :h_tgt_index,
                                     m_hp2);
    }
    m_source_index=get(m_vertex_ipmap, hp1->vertex());
  }

  ///called when the walk is over but the last segment is on an existing edge
  void on_walk_end_on_edge(Halfedge_handle hp1)
  {
    on_input_edge(hp1);
  }

  std::size_t last_seen_vertex_index() const
  {
    return m_source_index;
  }

  void set_beta(const Root_of_2 &beta)
  {
    m_last_beta=beta;
  }

  const Root_of_2 &last_beta() const
  {
    return m_last_beta;
  }
};


}
} //end of namespace CGAL::Surface_mesh_overlay_3

#endif // CGAL_WALK_ON_POLYHEDRON_VISITOR_H
