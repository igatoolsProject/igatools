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

#ifndef CGAL_WALK_ON_POLYHEDRON_H
#define CGAL_WALK_ON_POLYHEDRON_H

#include <CGAL/assertions.h>
#include <CGAL/Polyhedron_simplex_type.h>
#include <boost/optional.hpp>
#include <boost/variant.hpp>
#include <set>
#include <vector>
#include <stdexcept>

/// \todo make sure the code does not suppose we have a quad mesh

#ifdef CGAL_DEBUG_WALK_ON_POLYHEDRON
#include <iostream>
#define CGAL_WALKER_TRACE(X) std::cout << X;
#else
#define CGAL_WALKER_TRACE(X)
#endif

namespace CGAL
{

namespace internal
{
struct Walk_on_polyhedron_border_crossing_exception : std::logic_error
{
  Walk_on_polyhedron_border_crossing_exception()
    : std::logic_error("Trying to walking outside of the polyhedron")
  {}
};
}

template< class Polyhedron,
          class EdgeIntersectionPredicate,
          class EdgeIntersectionVisitor >
class Walk_on_polyhedron
{
  typedef typename  Polyhedron::Vertex_handle Vertex_handle;
  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
  typedef typename Polyhedron::Facet_handle Facet_handle;

  Halfedge_handle
  canonical_edge(Halfedge_handle h) const
  {
    return h<h->opposite()?h:h->opposite();
  }

  /// get the set of facets incident to a given simplex
  void get_incident_facets(Polyhedron_simplex_type simplex_type,
                           Halfedge_handle h,
                           std::set<Facet_handle> &facets) const
  {
    switch (simplex_type)
    {
      case POLYHEDRON_FACET:
        facets.insert(h->facet());
        break;
      case POLYHEDRON_EDGE:
        if (!h->is_border()) facets.insert(h->facet());
        if (!h->opposite()->is_border()) facets.insert(h->opposite()->facet());
        break;
      case POLYHEDRON_VERTEX:
      {
        Halfedge_handle start=h;
        do
        {
          if (!h->is_border())
            facets.insert(h->facet());
          h=h->next()->opposite();
        }
        while (start!=h);
      }
      break;
      case POLYHEDRON_NONE:
#if 0
        CGAL_assertion(!"Should never be here");
#endif
        break;
    }
  }

  Halfedge_handle
  common_edge(Facet_handle f1,
              Facet_handle f2) const
  {
    CGAL_assertion(f1!=Facet_handle());
    CGAL_assertion(f2!=Facet_handle());
    Halfedge_handle curr = f1->halfedge();
    do
    {
      if (curr->opposite()->face() == f2) return curr;
      curr=curr->next();
    }
    while (curr!=f1->halfedge());

    CGAL_assertion(!"We should not get here");

    return typename Polyhedron::Halfedge_handle();
  }

  /// info returned by the predicate when an edge of the input polyhedra
  /// is intersected in its interior
  typedef typename EdgeIntersectionPredicate::Barycentric_NT Barycentric_NT;

public:

  /// We walk from (source, source_type) to (target, target_type) intersecting
  /// all edges encountered. The predicate is responsible to indicate the correct
  /// direction to move in case it's a cycle.
  bool operator()(
    const Polyhedron &,
    Halfedge_handle source,
    Polyhedron_simplex_type source_type,
    Halfedge_handle target,
    Polyhedron_simplex_type target_type,
    const EdgeIntersectionPredicate &inside_edge_pred,
    EdgeIntersectionVisitor &visitor
  ) const
  {
    /// First check that the projected vertices do not fall in a common simplex
    Halfedge_handle common_facet;
    bool new_edge_to_add = true;
    switch (source_type)
    {
      case POLYHEDRON_FACET:
        CGAL_WALKER_TRACE("\n1-POLYHEDRON_FACET\n")
        switch (target_type)
        {
          case POLYHEDRON_FACET:
            CGAL_WALKER_TRACE("2-POLYHEDRON_FACET\n")
            if (source->facet() == target->facet())
              common_facet=source;
            break;
          case POLYHEDRON_EDGE:
            CGAL_WALKER_TRACE("2-POLYHEDRON_EDGE\n")
            if (source->facet()==target->facet() ||
                source->facet()==target->opposite()->facet())
              common_facet=source;
            break;
          case POLYHEDRON_VERTEX:
          {
            CGAL_WALKER_TRACE("2-POLYHEDRON_VERTEX\n")
            Halfedge_handle start=target;
            do
            {
              if (target->facet() == source->facet())
              {
                common_facet=source;
                break;
              }
              target=target->next()->opposite();
            }
            while (target!=start);
            target=start;
          }
          break;
          case POLYHEDRON_NONE:
            CGAL_WALKER_TRACE("2-POLYHEDRON_NONE\n")
#if 0
            CGAL_assertion("Should never be here\n");
#endif
            break;
        }
        break;
      case POLYHEDRON_EDGE:
        CGAL_WALKER_TRACE("\n1-POLYHEDRON_EDGE\n")
        switch (target_type)
        {
          case POLYHEDRON_FACET:
            CGAL_WALKER_TRACE("2-POLYHEDRON_FACET\n")
            if (target->facet()==source->facet() ||
                target->facet()==source->opposite()->facet())
              common_facet=target;
            break;
          case POLYHEDRON_EDGE:
            CGAL_WALKER_TRACE("2-POLYHEDRON_EDGE\n")
            if (canonical_edge(source)==canonical_edge(target))
            {
              /// we explicitly do not update common_facet
              /// so that the visitor knows source and target are completly
              /// inside an edge
              new_edge_to_add=false;
            }
            else
            {
              if (!source->is_border() &&
                  (source->facet() == target->facet() ||
                   source->facet() == target->opposite()->facet()))
                common_facet=source;
              else
              {
                if (!source->opposite()->is_border())
                {
                  if (target->facet() == source->opposite()->facet())
                    common_facet=target;
                  else if (target->opposite()->facet() == source->opposite()->facet())
                    common_facet=target->opposite();
                }
              }
            }
            break;
          case POLYHEDRON_VERTEX:
            CGAL_WALKER_TRACE("2-POLYHEDRON_VERTEX\n")
            if (source->vertex()==target->vertex() ||
                source->opposite()->vertex()==target->vertex())
            {
              common_facet = source->vertex()==target->vertex()?source:source->opposite();
              new_edge_to_add=false;
            }
            else
            {
              if (!source->is_border() &&
                  source->next()->vertex() == target->vertex())
                common_facet=source;
              else if (!source->opposite()->is_border() &&
                       source->opposite()->next()->vertex() == target->vertex())
                common_facet=source->opposite();
            }
            break;
          case POLYHEDRON_NONE:
            CGAL_WALKER_TRACE("2-POLYHEDRON_NONE\n")
#if 0
            CGAL_assertion("Should never be here\n");
#endif
            break;
        }
        break;
      case POLYHEDRON_VERTEX:
        CGAL_WALKER_TRACE("\n1-POLYHEDRON_VERTEX\n")
        switch (target_type)
        {
          case POLYHEDRON_FACET:
          {
            CGAL_WALKER_TRACE("2-POLYHEDRON_FACET\n")
            Halfedge_handle start=source;
            do
            {
              if (target->facet() == source->facet())
              {
                common_facet=target;
                break;
              }
              source=source->next()->opposite();
            }
            while (source!=start);
            source=start;
          }
          break;
          case POLYHEDRON_EDGE:
            CGAL_WALKER_TRACE("2-POLYHEDRON_EDGE\n")
            if (target->vertex()==source->vertex() ||
                target->opposite()->vertex()==source->vertex())
            {
              common_facet=target->vertex()==source->vertex()?
                           target->opposite():target;
              new_edge_to_add=false;
            }
            else
            {
              if (!target->is_border() &&
                  target->next()->vertex() == source->vertex())
                common_facet=target;
              else if (!target->opposite()->is_border() &&
                       target->opposite()->next()->vertex() == source->vertex())
                common_facet=target->opposite();
            }
            break;
          case POLYHEDRON_VERTEX:
          {
            /// \todo Handle this case and the similar one on an edge and one a face
            ///       (two points with the same projection)
            CGAL_assertion(source->vertex()!=target->vertex() || !"NOT ALREADY HANDLED");
            CGAL_WALKER_TRACE("2-POLYHEDRON_VERTEX\n")
            Halfedge_handle start=source;
            do
            {
              if (source->opposite()->vertex()==target->vertex())
              {
                common_facet=source->opposite();
                new_edge_to_add=false;
                break;
              }
              source=source->next()->opposite();
            }
            while (source!=start);
          }
          break;
          case POLYHEDRON_NONE:
            CGAL_WALKER_TRACE("2-POLYHEDRON_NONE\n")
#if 0
            CGAL_assertion("Should never be here\n");
#endif
            break;
        }
        break;
      case POLYHEDRON_NONE:
        CGAL_assertion("Should never be here\n");
    }

    if (!new_edge_to_add)
    {
      // In this case, the edge to add is already on an edge of the input mesh
      visitor.on_walk_end_on_edge(common_facet);
      return true;
    }
    if (common_facet != Halfedge_handle())
    {
      CGAL_assertion(!common_facet->is_border());
      visitor.on_walk_end(common_facet->facet());
      return true;
    }


    /// We now do the walk in the triangulation and find edge-edge intersection points
    std::set <Facet_handle> target_facets; /// \todo try using a boost::flat_set instead
    get_incident_facets(target_type, target, target_facets);
    std::set <Facet_handle> facets_already_visited;

    while (true)
    {
      std::vector<Halfedge_handle> hedges_to_test;
      switch (source_type)
      {
        case POLYHEDRON_FACET:
          CGAL_WALKER_TRACE("   i-POLYHEDRON_FACET\n")
          hedges_to_test.push_back(source);
          hedges_to_test.push_back(source->next());
          hedges_to_test.push_back(source->next()->next());
          break;
        case POLYHEDRON_EDGE:
          CGAL_WALKER_TRACE("   i-POLYHEDRON_EDGE\n")
          // if the v_src is on an edge, we don't know in which direction to go
          // This is only true the first time of the loop
          if (facets_already_visited.empty() &&
              !source->opposite()->is_border())
          {
            hedges_to_test.push_back(source->opposite()->next());
            hedges_to_test.push_back(source->opposite()->next()->next());
          }
          hedges_to_test.push_back(source->next());
          hedges_to_test.push_back(source->next()->next());
          break;
        case POLYHEDRON_VERTEX:
        {
          CGAL_WALKER_TRACE("   i-POLYHEDRON_VERTEX\n")
          Halfedge_handle start=source;
          do
          {
            //we don't want to go back, so we filter facets
            if (!facets_already_visited.count(source->facet()))
              if (!source->is_border())
              {
                hedges_to_test.push_back(source->next()->next());
              }
            source=source->next()->opposite();
          }
          while (source!=start);
        }
        break;
        case POLYHEDRON_NONE:
          CGAL_assertion("Should never be here\n");
      }

      Halfedge_handle intersected_edge;
      boost::optional< std::pair<Barycentric_NT, boost::variant<bool, Barycentric_NT> > > inter_res;

      CGAL_WALKER_TRACE("   Looking for intersections hedges_to_test.size()="<< hedges_to_test.size() << "\n")
      std::size_t k=0;
      for (; k<hedges_to_test.size(); ++k)
      {
        CGAL_WALKER_TRACE("     k="<< k << ", testing "
                          << hedges_to_test[k]->opposite()->vertex()->point()
                          << " "
                          << hedges_to_test[k]->vertex()->point()
                          << "\n")
        intersected_edge=hedges_to_test[k];
        inter_res =  inside_edge_pred(canonical_edge(intersected_edge));
        if (inter_res)
        {
          visitor.set_beta(inter_res->first);
          CGAL_WALKER_TRACE("     intersection found\n")
#if 0
          if (intersected_edge->is_border_edge())
            throw internal::Walk_on_polyhedron_border_crossing_exception();
          // \todo use the boolean in the API about incomplete overlay
#endif
          break;
        }
      }

      if (k==hedges_to_test.size())
      {
        // this handles the case when the source is on the border
        // and the walk directly goes out
        if (target_type==POLYHEDRON_NONE) break;
        // this handle the case when the source is on the border
        // and the target is inside or on the boundary of the supporting mesh
        // but the path is going out
        CGAL_WALKER_TRACE("WARNING: Halting but target not reached!\n")
        return false;


        CGAL_WALKER_TRACE("ERROR DID NOT FIND THE INTERSECTION!!!!!!!!!!!!!!!\n")
        CGAL_assertion(!"NO INTERSECTION FOUND");
        break;
      }

      std::set <Facet_handle> current_incident_facets;
      get_incident_facets(source_type, source, current_incident_facets);
      facets_already_visited.clear();

      const Barycentric_NT *barycentric_coord =
        boost::get<Barycentric_NT>(&(inter_res->second));
      if (!barycentric_coord)
      {
        // we reached a vertex
        Halfedge_handle hedge=canonical_edge(intersected_edge);
        bool is_target = boost::get<bool>(inter_res->second);
        if (!is_target) hedge=hedge->opposite();
        Vertex_handle vh=hedge->vertex();

        Polyhedron_simplex_type prev_type=source_type;
        Halfedge_handle prev_source=source;
        source_type=POLYHEDRON_VERTEX;
        source=hedge;

        std::vector<Halfedge_handle> common_facets;

        Halfedge_handle start=hedge;
        bool did_break=false;
        do
        {
          if (!hedge->is_border())    //note that target_facets and current_incident_facets do not contain any NULL face
          {
            if (!did_break && target_facets.count(hedge->facet()))
              did_break=true;
            if (current_incident_facets.count(hedge->facet()))
            {
              common_facets.push_back(hedge);
              facets_already_visited.insert(hedge->facet());
            }
          }
          hedge=hedge->next()->opposite();
        }
        while (hedge!=start);

        //insert current edge
        if (common_facets.size()==1)
        {
          // vertex-vertex case (on the mesh border)
          if (prev_type==POLYHEDRON_VERTEX)
          {
            CGAL_assertion(common_facets[0]->vertex()==vh);
            visitor.on_input_edge(
              common_facets[0]->opposite()->vertex()==prev_source->vertex()
              ? common_facets[0]
              : common_facets[0]->next()->opposite());
          }
          else
          {
            // edge->vertex case (on the mesh border)
            if (prev_type==POLYHEDRON_EDGE && (prev_source->vertex()==vh || prev_source->opposite()->vertex()==vh))
              visitor.on_input_edge(prev_source->vertex()==vh
                                    ? prev_source:prev_source->opposite());
            else
              // we reached a vertex opposite to the edge or in a face
              visitor.found_vertex_intersection(common_facets[0]);
          }
        }
        else
        {
          CGAL_assertion(common_facets.size()==2);
          //this is a vertex-vertex case or a edge-vertex
          Halfedge_handle edge =
            common_edge(common_facets[0]->facet(), common_facets[1]->facet());
          visitor.on_input_edge(edge->vertex()==vh?edge:edge->opposite());
        }

        if (did_break)
        {
          CGAL_assertion(hedge==start);
          common_facets.clear();
          do
          {
            if (!hedge->is_border() && target_facets.count(hedge->facet()))
              common_facets.push_back(hedge);
            hedge=hedge->next()->opposite();
          }
          while (hedge!=start);
          //insert the last edge
          if (common_facets.size()==1)
          {
            // the last visited edge is a complete border edge
            if (target_type==POLYHEDRON_VERTEX)
            {
              CGAL_assertion(common_facets[0]->vertex()==vh);
              visitor.on_walk_end_on_edge(
                common_facets[0]->opposite()->vertex()==target->vertex()
                ? common_facets[0]->opposite()
                : common_facets[0]->next());
            }
            else
            {
              // the walk ends on a border edge
              if (target_type==POLYHEDRON_EDGE &&
                  (target->vertex()==vh||target->opposite()->vertex()==vh))
              {
                visitor.on_walk_end_on_edge(target->vertex()==vh
                                            ? target->opposite():target);
              }
              else
                visitor.on_walk_end(common_facets[0]->facet());
            }
          }
          else
          {
            CGAL_assertion(common_facets.size()==2);
            //this is a vertex-vertex case or a edge-vertex
            Halfedge_handle edge =
              common_edge(common_facets[0]->facet(), common_facets[1]->facet());
            visitor.on_walk_end_on_edge(edge->vertex()==vh?edge->opposite():edge);
          }
          break;
        }
      }
      else
      {
        // we crossed the interior of an edge
        facets_already_visited.insert(source->facet());
        source_type=POLYHEDRON_EDGE;
        source=intersected_edge->opposite();
        // barycentric_coord must have been computed using `canonical_edge(intersected_edge)`
        visitor.found_edge_intersection(intersected_edge, *barycentric_coord);

        if (target_facets.find(source->facet())
            != target_facets.end())
        {
          CGAL_assertion(target_type != POLYHEDRON_VERTEX ||
                         (source->vertex()!=target->vertex() &&
                          source->opposite()->vertex() != target->vertex()));
          //insert the last edge
          if (!intersected_edge->opposite()->is_border())
            visitor.on_walk_end(intersected_edge->opposite()->facet());
          break;
        }
      }
      if (intersected_edge->is_border_edge())
      {
        if (target_type!=POLYHEDRON_NONE)
        {
          CGAL_WALKER_TRACE("WARNING: Halting but target not reached!\n")
          return false;
        }
        break; // the walk reached the boundary of the domain
      }
    }
    return true;
  }
};

}

#endif // CGAL_WALK_ON_POLYHEDRON_H
