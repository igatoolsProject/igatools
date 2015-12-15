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

#ifndef CGAL_SURFACE_MESH_OVERLAY_3_POLYHEDRON_REFINEMENT_MODIFIER_3_H
#define CGAL_SURFACE_MESH_OVERLAY_3_POLYHEDRON_REFINEMENT_MODIFIER_3_H

#include "Polyhedron_utils.h"

#include <CGAL/HalfedgeDS_decorator.h>
#include <CGAL/Modifier_base.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Root_of_traits.h>
#include <CGAL/tuple.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/linear_least_squares_fitting_3.h>

#include <boost/foreach.hpp>
#include <boost/math/special_functions/next.hpp>

//For faces triangulation
struct CDT_TMP_DEBUG {};
#define CGAL_CT2_WANTS_TO_HAVE_EXTRA_ACTION_FOR_INTERSECTING_CONSTRAINTS
#define CGAL_CDT2_EXTRA_ACTION_FOR_INTERSECTING_CONSTRAINTS throw CDT_TMP_DEBUG();

#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_2_projection_traits_3.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>

#include <vector>
#include <utility>
#if !defined(CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT) || \
    !defined(CGAL_NO_DEBUG_OVERLAY_3_TRIANGULATION)
#include <fstream>
#endif

namespace CGAL
{

namespace internal
{
struct First_of_pair_less
{
  template<class P1, class P2>
  bool operator()(const P1 &p1, const P2 &p2) const
  {
    return p1.first < p2.first;
  }
};
} //end of internal namespace


/// Important note to user of this class:
/// There is two different indexing of the edges and halfedges:
/// The polyhedron that is going to be walked on (and refined) require its
/// halfedge to be numbered from 0 to nb_edges with opposite edges having
/// the same id (like P1 in the mesh overlay).
/// The halfedges that are used by to define where we walk (P2 in the mesh overlay)
/// must be indexed from 0 to nb_halfedges_considered (which is the number which is passed
/// to the constructor of this class)
template < class Polyhedron,
           class Exact_kernel,
           class Input_kernel,
           class ExactPointPropertyMap,
           class VertexIndexPropertyMap,
           bool AllowPathIntersections = false>
class Polyhedron_refinement_modifier_3 :
  public Modifier_base<typename Polyhedron::HDS>
{
//typedef's
  typedef typename Polyhedron::Halfedge_handle                  Halfedge_handle;
  typedef typename Polyhedron::Facet_handle                        Facet_handle;
  typedef typename Polyhedron::Vertex_handle                      Vertex_handle;
  typedef typename Polyhedron::HDS                                          HDS;
  typedef typename HDS::Halfedge                                       Halfedge;
  typedef typename HDS::Face                                              Facet;
  typedef          HalfedgeDS_decorator<HDS>                          Decorator;
  typedef typename HDS::Halfedge::Base                                    HBase;
  typedef          std::size_t                                               ID;
  typedef          Exact_kernel                                              EK;
  typedef typename Root_of_traits<typename EK::FT>::Root_of_2         Root_of_2;
  typedef typename Input_kernel::Point_3                                Point_3;

//private members
  const std::size_t max_id;
  //associate to each halfedge of P2, the halfedges corresponding in P1
  std::vector< std::vector<Halfedge_handle> > m_P2_edge_to_P1_edges;
  ExactPointPropertyMap m_exact_ppmap;
  VertexIndexPropertyMap m_vertex_ipmap;
//converters
  CGAL::Cartesian_converter<EK,Input_kernel> to_input;

public:
  std::vector< Halfedge_handle > edges;
  std::vector< Facet_handle > facets;
  std::vector< Vertex_handle > vertices;


  std::vector< std::vector<std::pair<Root_of_2,std::size_t> > > vertices_on_edges;
  std::vector< std::vector<std::pair<Point_3,std::size_t> > > vertices_in_facets;
  std::vector< std::vector<std::pair<ID,ID> > > edges_in_facets;
  std::vector< std::vector<Halfedge_handle> > edges_in_facets_P2_halfedges;
  // in case a halfedge of P1 is partially but not entirely covered
  // by the projection of P2, we register the correspond edge portion
  // using a pair of vertex ids
  typedef cpp11::tuple<ID, ID, Halfedge_handle>         Vertex_ids_and_halfedge;
  std::vector< Vertex_ids_and_halfedge >    hedges_of_P1_covered_by_hedge_of_P2;

  // used to associated to a halfedge of P2 a halfedge of P1
  // that will be split to embed the vertices which indices
  // are in the tuple. The hedge in the tuple is from P2.
  // The key is the id of the original halfedge of P1
// containers to handle partial overlays
  std::vector< std::pair<Vertex_handle, std::size_t> >
  vertices_with_no_projection; // store vertices of P2 without projection
  // together with their indices in `vertices`

  typedef cpp11::tuple<Halfedge_handle, std::size_t, std::size_t>
  Hedge_and_P1_vid; // <Halfedge from P2, src id in P1, tgt id in P1>
  // edges from P2 without projection, or with partial projection
  std::vector< Hedge_and_P1_vid > halfedges_with_parts_outside_P1;
  std::size_t nb_faces_with_no_full_projections;
  // each halfedge of P2 with parts with no projection in P1 is associated
  // a ccb id
  std::map<Halfedge_handle, std::size_t> ccb_ids_in_P2;

#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
  std::ofstream debug_output; //only for debug purpose
  std::vector<Point_3> debug_points; //only for debug purpose
#endif

  std::size_t new_vertex_index;
  std::size_t nb_input_vertices;

  const std::vector< std::vector< Halfedge_handle > > &
  halfedge_map() const
  {
    return m_P2_edge_to_P1_edges;
  }

  Polyhedron_refinement_modifier_3(
    std::size_t max_nb_hedges,
    ExactPointPropertyMap ppmap,
    VertexIndexPropertyMap ipmap=VertexIndexPropertyMap()
  ) : max_id((std::numeric_limits<std::size_t>::max)())
    , m_P2_edge_to_P1_edges(max_nb_hedges)
    , m_exact_ppmap(ppmap)
    , m_vertex_ipmap(ipmap)
    , nb_faces_with_no_full_projections(0)
  {}

  void init()
  {
    vertices_on_edges.resize(edges.size());
    vertices_in_facets.resize(facets.size());
    edges_in_facets.resize(facets.size());
    edges_in_facets_P2_halfedges.resize(facets.size());
    new_vertex_index=vertices.size();
    nb_input_vertices=new_vertex_index;
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
    debug_output.open("debug/output.cgal");
#endif
  }

  std::size_t
  add_point_on_edge(Halfedge_handle h, const Root_of_2 &alpha)
  {
    /// \todo in case AllowPathIntersections==true, we need to handle the case when a point on an edge is added more than once
    CGAL_assertion(Surface_mesh_overlay_3::canonical_edge(h)==h);
    vertices_on_edges[ h->id() ].push_back(std::make_pair(alpha, new_vertex_index++));
    put(m_exact_ppmap,new_vertex_index-1, typename EK::Point_3());
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
    debug_points.push_back(to_input(approx_rational_point_on_edge(h,alpha)));
#endif
    return new_vertex_index-1;
  }

  std::size_t
  add_point_in_facet(Facet_handle f, const Point_3 &pt)
  {
    /// \todo in case AllowPathIntersections==true, we need to handle the case when a point in a facet is added more than once
    vertices_in_facets[ f->id() ].push_back(std::make_pair(pt, new_vertex_index++));
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
    debug_points.push_back(pt);
#endif
    return new_vertex_index-1;
  }

  std::size_t
  add_vertex_with_no_projection(Vertex_handle vh)
  {
    vertices_with_no_projection.push_back(std::make_pair(vh,new_vertex_index));
    ++new_vertex_index;
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
    debug_points.push_back(vh->point());
#endif
    return new_vertex_index-1;
  }

  template <class EdgeSkipper>
  void set_ccb_id(Halfedge_handle h, const EdgeSkipper &skip)
  {
    if (ccb_ids_in_P2.insert(
          std::make_pair(h,nb_faces_with_no_full_projections)).second)
    {
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
      int nb=1;
#endif
      //set the ccb id of all halfedges
      Halfedge_handle next=h;
      do
      {
        next=next->next();
        while (skip(next))
          next=next->opposite()->next();
        if (next==h) break;
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
        ++nb;
#endif
        CGAL_assertion_code(bool insert_ok =)
        ccb_ids_in_P2.insert(
          std::make_pair(next,nb_faces_with_no_full_projections)).second;
        CGAL_assertion(insert_ok);
      }
      while (true);
      ++nb_faces_with_no_full_projections;
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
      if (h->is_border())
        std::cout << "New CCB (border) ";
      else
        std::cout << "New CCB ";
      std::cout << nb_faces_with_no_full_projections-1 << " with " << nb << " hedges\n";
#endif
    }
  }

  template <class EdgeSkipper>
  void add_edges_with_partial_projection(
    Halfedge_handle h, std::size_t src_id, std::size_t tgt_id,
    const EdgeSkipper &skip)
  {
    halfedges_with_parts_outside_P1.push_back(
      Hedge_and_P1_vid(h,src_id,tgt_id));
    halfedges_with_parts_outside_P1.push_back(
      Hedge_and_P1_vid(h->opposite(),tgt_id,src_id));
    set_ccb_id(h, skip);
    set_ccb_id(h->opposite(), skip);
  }

  template <class EdgeSkipper>
  void add_edges_with_no_projection(Halfedge_handle h, const EdgeSkipper &skip)
  {
    add_edges_with_partial_projection(h, max_id, max_id, skip);
  }

  /// add a new edge in facet `f` with vertex `i` as source
  /// and vertex `j` as target
  /// This new edge correspond to (a part) of the edge `supporting_halfedge` from P2
  void add_edge(std::size_t i, std::size_t j,
                Facet_handle f,
                Halfedge_handle supporting_halfedge)
  {
    CGAL_assertion(f!= Facet_handle());
    edges_in_facets[f->id()].push_back(std::make_pair(i,j));
    edges_in_facets_P2_halfedges[f->id()].push_back(supporting_halfedge);
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
    debug_output << "2 " << get_debug_point(i) << " "
                 << get_debug_point(j) << "\n";
#endif
  }

  /// The projection of `supporting_edge` entirely contains `hedge`.
  /// `hedge` is always an input edge.
  void update_edge_incidence(Halfedge_handle hedge,
                             Halfedge_handle supporting_hedge)
  {
    m_P2_edge_to_P1_edges[supporting_hedge->id()].push_back(hedge);
    m_P2_edge_to_P1_edges[supporting_hedge->opposite()->id()].push_back(hedge->opposite());
  }

  /// the projection of h_P2 is overlapping with h_P1 but the intersection is not h_P1
  void record_edge_overlap(Halfedge_handle /*h_P1*/,
                           ID v1_id, ID v2_id, Halfedge_handle h_P2)
  {
    hedges_of_P1_covered_by_hedge_of_P2.push_back(
      Vertex_ids_and_halfedge(v1_id, v2_id, h_P2));
  }

#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
  const Point_3 &
  get_debug_point(std::size_t i) const
  {
    if (i< vertices.size())
      return vertices[i]->point();
    CGAL_assertion(i-vertices.size() < debug_points.size());
    return debug_points[ i-vertices.size() ];
  }
#endif

  typename EK::Point_3
  approx_rational_point_on_edge(Halfedge_handle edge, const Root_of_2 &alpha)
  {
    typename EK::Point_3 b0 = m_exact_ppmap[ get(m_vertex_ipmap, edge->opposite()->vertex()) ];
    typename EK::Point_3 b1 = m_exact_ppmap[ get(m_vertex_ipmap, edge->vertex()) ];

    return b0 + typename EK::FT(to_double(alpha)) * (b1-b0);
  }

  void create_interior_hedges_and_triangulate_facet(
    HDS &hds,
    Decorator &decorator,
    std::size_t facet_id,
    bool only_border_vertices_and_non_convex_face=false)
  {
    typedef Triangulation_2_projection_traits_3<EK>      Traits;
    typedef Triangulation_vertex_base_with_info_2<std::size_t,Traits> Vb;
    typedef Constrained_triangulation_face_base_2<Traits>            Fbb;
    typedef Triangulation_face_base_with_info_2<bool,Traits,Fbb>      Fb;
    typedef Triangulation_data_structure_2<Vb,Fb>                    TDS;
    typedef No_intersection_tag                                     Itag;
    typedef Constrained_Delaunay_triangulation_2<Traits,
            TDS,
            Itag>               CDT;

#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
    debug_output.close();
#endif

    std::vector< std::pair<typename EK::Point_3,std::size_t> > points;
    std::vector< typename EK::Point_3 > pfn;
    //map a pair of vertices in the CDT to the halfedge in the output
    // and the original supporting edge in P2
    typedef std::map< std::pair<std::size_t, std::size_t>,
            std::pair<Halfedge_handle,Halfedge_handle> > Halfedge_map;
    Halfedge_map halfedges;

    Halfedge_handle h=facets[facet_id]->halfedge(), start=h;

    std::size_t nb_border_vertices=0;
    do
    {
      ++nb_border_vertices;
      h=h->next();
    }
    while (h!=start);

    if (nb_border_vertices==3 &&
        (only_border_vertices_and_non_convex_face ||
         vertices_in_facets[facet_id].empty()))
    {
      CGAL_assertion(facets[facet_id]->id()==facet_id);
      return;
    }

#ifndef CGAL_NO_DEBUG_OVERLAY_3_TRIANGULATION
    std::cout << "Triangulating facet " << facet_id << "\n";
#endif

    for (std::size_t i=0; i<nb_border_vertices; ++i)
    {
#ifndef CGAL_NO_DEBUG_OVERLAY_3_TRIANGULATION
      std::cout << "BORDER PT  " << h->vertex()->point() << " - " << get(m_vertex_ipmap, h->vertex()) << "\n";
#endif
      const typename EK::Point_3 &pt=m_exact_ppmap[ get(m_vertex_ipmap, h->vertex()) ];
      /// \todo this might be incorrect if one uses reusable indices and input-created vertices are mixed up...
      /// we look for input vertices because they define the original triangle and aren't collinear
      if (get(m_vertex_ipmap, h->vertex()) < nb_input_vertices)
        pfn.push_back(pt);
      points.push_back(std::make_pair(pt,get(m_vertex_ipmap, h->vertex())));
      halfedges[
        std::make_pair(get(m_vertex_ipmap, h->opposite()->vertex()),get(m_vertex_ipmap, h->vertex()))
      ]=std::make_pair(h, Halfedge_handle());
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
      std::cout << "register border edge "
                << get(m_vertex_ipmap, h->opposite()->vertex()) << " "
                << get(m_vertex_ipmap, h->vertex()) << "\n";
#endif
      h=h->next();
    }

    if (!only_border_vertices_and_non_convex_face)
    {
      std::size_t nb_pt_inside=vertices_in_facets[facet_id].size();
      for (std::size_t i=0; i<nb_pt_inside; ++i)
      {
#ifndef CGAL_NO_DEBUG_OVERLAY_3_TRIANGULATION
        std::cout << "INTERIOR PT  " << m_exact_ppmap[ vertices_in_facets[facet_id][i].second ] << " - " << vertices_in_facets[facet_id][i].second << "\n";
#endif

        points.push_back(
          std::make_pair(
            m_exact_ppmap[ vertices_in_facets[facet_id][i].second ],
            vertices_in_facets[facet_id][i].second));
      }
    }

    typename EK::Vector_3 face_normal;
    if (only_border_vertices_and_non_convex_face)
    {
      std::vector<typename Input_kernel::Point_3> pts;
      for (std::size_t i=0; i<nb_border_vertices; ++i, h=h->next())
        pts.push_back(h->vertex()->point());
      typename Input_kernel::Plane_3 plane;
      linear_least_squares_fitting_3(pts.begin(), pts.end(),plane,Dimension_tag<0>());
      Cartesian_converter<Input_kernel,EK> to_exact;
      face_normal=to_exact(plane.orthogonal_vector());
    }
    else
    {
      CGAL_assertion(pfn.size()>=3);
      face_normal=normal(pfn[0], pfn[1], pfn[2]);
    }
    Traits traits(face_normal);
    CDT cdt(traits);
    cdt.insert(points.begin(), points.end());

    CGAL_assertion(points.size()==cdt.number_of_vertices());   //check no vertex is duplicated

    /// \todo you can do better than that!!!
    std::map<std::size_t, typename CDT::Vertex_handle> cdt_vertices;
    for (typename CDT::Finite_vertices_iterator
         vit=cdt.finite_vertices_begin(), vit_end=cdt.finite_vertices_end();
         vit!=vit_end; ++vit)
    {
      cdt_vertices[vit->info()]=vit;
    }
/// \todo Question: what should we do in the case of failure?
    try
    {
      std::size_t nb_interior_edges = only_border_vertices_and_non_convex_face?
                                      0:edges_in_facets[facet_id].size();
#ifndef CGAL_NO_DEBUG_OVERLAY_3_TRIANGULATION
      std::ofstream output("debug/interior_constraints.cgal");
      for (std::size_t i=0; i<nb_interior_edges; ++i)
      {
        std::pair<std::size_t, std::size_t> indices=edges_in_facets[facet_id][i];
        CGAL_assertion(cdt_vertices.count(indices.first) && cdt_vertices.count(indices.second));
        output << "2 " << cdt_vertices[indices.first]->point() << " " <<
               cdt_vertices[indices.second]->point() << "\n";
      }
      output.close();
#endif

      for (std::size_t i=0; i<nb_interior_edges; ++i)
      {
        std::pair<std::size_t, std::size_t> indices=edges_in_facets[facet_id][i];
        CGAL_assertion(cdt_vertices[indices.first]!=typename CDT::Vertex_handle());
        CGAL_assertion(cdt_vertices[indices.second]!=typename CDT::Vertex_handle());

        std::pair<typename Halfedge_map::iterator, bool> insert_res =
          halfedges.insert(std::make_pair(indices,
                                          std::make_pair(Halfedge_handle(),
                                                         Halfedge_handle()))),
          insert_res_opp=
            halfedges.insert(std::make_pair(std::make_pair(indices.second, indices.first),
                                            std::make_pair(Halfedge_handle(),
                                                           Halfedge_handle())));
        // check the halfedge is not a facet border halfedge
        if (insert_res.second && insert_res_opp.second)
        {
          Halfedge_handle nh
            =hds.edges_push_back(Halfedge(), Halfedge());

          decorator.set_vertex(nh, vertices[indices.second]);
          decorator.set_vertex(nh->opposite(), vertices[indices.first]);
          decorator.set_vertex_halfedge(vertices[indices.first],nh->opposite());
          decorator.set_vertex_halfedge(vertices[indices.second],nh);
          Halfedge_handle support_P2 = edges_in_facets_P2_halfedges[facet_id][i];

          insert_res.first->second=std::make_pair(nh, support_P2);
          insert_res_opp.first->second=std::make_pair(nh->opposite(),
                                                      support_P2->opposite());

          /// insert the contraints in the CDT
          cdt.insert_constraint(cdt_vertices[indices.first],
                                cdt_vertices[indices.second]);
        }


        //check the edge is there
        CGAL_assertion(cdt.is_edge(cdt_vertices[indices.first],
                                   cdt_vertices[indices.second]));
      }
    }
    catch (CDT_TMP_DEBUG)
    {
      std::cout << "Intersection of constraints!!!!!\n";
      CGAL_assertion(!"Intersection of constraints!!!!!\n");
      return;
    }

    // constraint the border, it is important if the face triangulated is
    //  not convex (only if only_border_vertices_and_non_convex_face==true)
    for (std::size_t i=0; i<nb_border_vertices; ++i)
    {
      CGAL_assertion(cdt_vertices[points[i].second] != typename CDT::Vertex_handle());
      CGAL_assertion(cdt_vertices[points[(i+1)%nb_border_vertices].second] != typename CDT::Vertex_handle());

      cdt.insert_constraint(cdt_vertices[points[i].second],
                            cdt_vertices[points[(i+1)%nb_border_vertices].second]);
    }

    //mark faces outside the face domain if not convex
    if (only_border_vertices_and_non_convex_face)
    {
      for (typename CDT::All_faces_iterator fit = cdt.all_faces_begin(),
           end = cdt.all_faces_end();
           fit != end; ++fit)
      {
        fit->info()=true;
      }
      // starting from an unbounded face, flood the faces outside the domain
      // defined by the polyhedron facet
      std::vector<typename CDT::Face_handle> queue;
      queue.push_back(cdt.infinite_face());
      while (! queue.empty())
      {
        typename CDT::Face_handle fh = queue.back();
        queue.pop_back();
        if (fh->info() == true)
        {
          fh->info()=false;
          for (int i = 0; i < 3; ++i)
            if (!fh->is_constrained(i)) queue.push_back(fh->neighbor(i));
        }
      }
    }

    // add additionnal edges that are here to have a triangulated polyhedron
    // if one need to mark them, it's here that it should be done
    for (typename CDT::Finite_edges_iterator
         eit = cdt.finite_edges_begin(),
         end = cdt.finite_edges_end();
         eit != end; ++eit)
    {
      //need to add edges not already added
      if (!cdt.is_constrained(*eit) &&
          (!only_border_vertices_and_non_convex_face || eit->first->info()))
      {
        Halfedge_handle h=hds.edges_push_back(Halfedge(), Halfedge());

        std::size_t i1=eit->first->vertex(CDT::ccw(eit->second))->info();
        std::size_t i2=eit->first->vertex(CDT::cw(eit->second))->info();
        halfedges[ std::make_pair(i1,i2) ] = std::make_pair(h, Halfedge_handle());
        halfedges[ std::make_pair(i2,i1) ] = std::make_pair(h->opposite(), Halfedge_handle());
        decorator.set_vertex(h, vertices[i2]);
        decorator.set_vertex(h->opposite(), vertices[i1]);
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
        std::cout << "registering edge " << i1 << " " << i2 <<"\n";
#endif
      }
    }

    // In case the face is non-convex we might have used points
    // that define a normal that is opposite to the real face
    // normal. In that case, triangle orientations should be reversed
    bool reverse_triangulation=false;
    if (only_border_vertices_and_non_convex_face)
    {
      std::size_t i1 = get(m_vertex_ipmap, h->opposite()->vertex()),
                  i2 = get(m_vertex_ipmap, h->vertex());
      typename CDT::Face_handle fh;
      int i;
      CGAL_assertion_code(bool is_edge =)
      cdt.is_edge(cdt_vertices[i1], cdt_vertices[i2],fh,i);
      CGAL_assertion(is_edge);
      if (fh->info()) reverse_triangulation=true;
    }

    bool first_run=true;
    Facet_handle the_face=facets[facet_id];
    for (typename CDT::Finite_faces_iterator fit=cdt.finite_faces_begin(),
         fit_end=cdt.finite_faces_end();
         fit!=fit_end; ++fit)
    {
      if (only_border_vertices_and_non_convex_face && !fit->info()) continue;
      if (!first_run)
      {
        the_face=hds.faces_push_back(*the_face);
        the_face->id()=facets.size();
        facets.push_back(the_face);
      }

      first_run=false;
      std::size_t indices[]=
      {
        fit->vertex(0)->info(),
        fit->vertex(1)->info(),
        fit->vertex(2)->info()
      };
      if (reverse_triangulation) std::swap(indices[0], indices[1]);
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
      std::cout << "requesting edge " << indices[0] << " " << indices[1] <<"\n";
      std::cout << "requesting edge " << indices[1] << " " << indices[2] <<"\n";
      std::cout << "requesting edge " << indices[2] << " " << indices[0] <<"\n";
#endif

      Halfedge_handle h1, h2, h3; //output halfedges
      Halfedge_handle sh1, sh2, sh3; //supporting halfedges from P2
      cpp11::tie(h1,sh1)=halfedges[std::make_pair(indices[0],indices[1])];
      cpp11::tie(h2,sh2)=halfedges[std::make_pair(indices[1],indices[2])];
      cpp11::tie(h3,sh3)=halfedges[std::make_pair(indices[2],indices[0])];

      if (sh1!=Halfedge_handle())
        m_P2_edge_to_P1_edges[sh1->id()].push_back(h1);
      if (sh2!=Halfedge_handle())
        m_P2_edge_to_P1_edges[sh2->id()].push_back(h2);
      if (sh3!=Halfedge_handle())
        m_P2_edge_to_P1_edges[sh3->id()].push_back(h3);


      CGAL_assertion(h1 != Halfedge_handle());
      CGAL_assertion(h2 != Halfedge_handle());
      CGAL_assertion(h3 != Halfedge_handle());

      //set next-prev relations
      h1->HBase::set_next(h2);
      h2->HBase::set_next(h3);
      h3->HBase::set_next(h1);
      decorator.set_prev(h2, h1);
      decorator.set_prev(h3, h2);
      decorator.set_prev(h1, h3);
      //set face-halfedge relation
      decorator.set_face(h1,the_face);
      decorator.set_face(h2,the_face);
      decorator.set_face(h3,the_face);
      decorator.set_face_halfedge(the_face,h1);
    }
  }

  struct Less_ptid_along_edge
  {
    typedef typename boost::property_traits< ExactPointPropertyMap >::value_type Exact_point_3;
    std::size_t m_src_id;
    ExactPointPropertyMap m_exact_ppmap;

    Less_ptid_along_edge(std::size_t src_id, ExactPointPropertyMap exact_ppmap)
      : m_src_id(src_id)
      , m_exact_ppmap(exact_ppmap)
    {}

    bool operator()(std::size_t i, std::size_t j) const
    {
      return compare_distance_to_point(get(m_exact_ppmap, m_src_id),
                                       get(m_exact_ppmap, i),
                                       get(m_exact_ppmap, j));
    }
  };

  /// In case we are in a scenario where edges to insert into facets
  /// might intersect, we split the edges intersecting at the intersection
  /// points to guarantee that no constrained edges will intersect when
  /// retriangulating the facets using the CDT_2
  void handle_intersection_of_constraints(Decorator &decorator)
  {
    std::size_t nb_facets = edges_in_facets.size();
    for (std::size_t fid=0; fid<nb_facets; ++fid)
    {
      std::vector<std::pair<ID,ID> > &constraints = edges_in_facets[fid];
      if (constraints.empty()) continue;

      std::size_t nb_edges = constraints.size();
      typedef typename EK::Segment_3 Exact_segment_3;
      typedef typename EK::Point_3 Exact_point_3;
      // container maintaining indices sorted along the constrained edge from its source to its target
      typedef std::set<std::size_t, Less_ptid_along_edge> Index_set;

      std::vector< Index_set > intersections_per_edge;
      intersections_per_edge.reserve(nb_edges);
      std::map<Exact_point_3, std::size_t> intersection_points;

      //init the index containers
      for (std::size_t i=0; i< nb_edges; ++i)
      {
        Less_ptid_along_edge less(constraints[i].first, m_exact_ppmap);
        intersections_per_edge.push_back(Index_set(less));
      }

      for (std::size_t i=0; i< nb_edges-1; ++i)
      {
        // create s1 representing the edge inside the facet
        Exact_segment_3 s1(m_exact_ppmap[constraints[i].first], m_exact_ppmap[constraints[i].second]);

        // look for points inside s1
        typedef std::pair<Point_3, std::size_t> Pair;
        BOOST_FOREACH(const Pair& p, vertices_in_facets[fid])
        {
          if (p.second == constraints[i].first ||
              p.second==constraints[i].second) continue;
          const Exact_point_3 &pt = m_exact_ppmap[p.second];
          if (s1.has_on(pt))
          {
            std::pair<
            typename std::map<Exact_point_3, std::size_t>::iterator,
                     bool
                     > insert_res =
                       intersection_points.insert(std::make_pair(pt, p.second));
            CGAL_assertion(insert_res.first->second==p.second);
            intersections_per_edge[i].insert(p.second);
          }
        }

        //look for intersection points with another segment
        for (std::size_t j=i+1; j< nb_edges; ++j)
        {
          Exact_segment_3 s2(m_exact_ppmap[constraints[j].first], m_exact_ppmap[constraints[j].second]);

          /// skip pair of segments sharing a common end-point
          /// \todo check this is really what we want
          if (constraints[i].first == constraints[j].first ||
              constraints[i].first == constraints[j].second ||
              constraints[i].second == constraints[j].first ||
              constraints[i].second == constraints[j].second) continue;

          if (CGAL::do_intersect(s1, s2))
          {
            boost::optional< boost::variant<Exact_point_3, Exact_segment_3> > res =
              CGAL::intersection(s1,s2);
#if 0 // commented by martinelli, Dec 15, 2015
            CGAL_assertion(res);
#endif
            if (const Exact_point_3 *pt_ptr = get<Exact_point_3>(&(*res)))
            {
              std::pair<
              typename std::map<Exact_point_3, std::size_t>::iterator,
                       bool
                       > insert_res =
                         intersection_points.insert(std::make_pair(*pt_ptr, new_vertex_index));
              if (insert_res.second)
              {
                std::size_t vid = new_vertex_index++;
                put(m_exact_ppmap, vid, *pt_ptr);
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
                debug_points.push_back(to_input(*pt_ptr));
#endif
                CGAL_assertion(vertices.size()==vid);
                vertices.push_back(decorator.vertices_push_back(
                                     typename HDS::Vertex(to_input(*pt_ptr))));
                put(m_vertex_ipmap, vertices.back(), vid);
                vertices_in_facets[ fid ].push_back(std::make_pair(vertices.back()->point(), vid));
              }

              intersections_per_edge[i].insert(insert_res.first->second);
              intersections_per_edge[j].insert(insert_res.first->second);
            }
            else
            {
              CGAL_assertion(get<Exact_segment_3>(&(*res)));
              CGAL_assertion(!"HANDLE ME");
              /// \todo Handle this case
              /// all will be fine for the CDT2 if we do nothing but not for porting
              /// the corresponding projected edges
              /// This becomes arbitrarily compilated as a pair of halfedges overlapping might
              /// be intersecting by another edge (in their intersection), thus leading to
              /// a duplicated intersection point (shall we have a skip list?)
            }
          }
        }
      }

      // update the set of constrained edges to insert in the facet
      // split each constrained at intersection points and inside points
      for (std::size_t i=0; i< nb_edges; ++i)
      {
        Index_set &interpt_ids = intersections_per_edge[i];
        if (interpt_ids.empty()) continue;

        std::pair<ID,ID> id_pair = constraints[i];
        Halfedge_handle hprojected = edges_in_facets_P2_halfedges[fid][i];

        typename Index_set::iterator iti=interpt_ids.begin(), iti_end=interpt_ids.end();
        constraints[i] = std::pair<ID, ID>(id_pair.first, *iti);   // reuse the space in the vector
        std::size_t prev=*iti++;
        for (; iti!=iti_end; ++iti)
        {
          constraints.push_back(std::pair<ID, ID>(prev, *iti));
          edges_in_facets_P2_halfedges[fid].push_back(hprojected);
          prev=*iti;
        }
        constraints.push_back(std::pair<ID, ID>(prev, id_pair.second));
        edges_in_facets_P2_halfedges[fid].push_back(hprojected);
      }
    }
  }

  template <class EdgeOutputIterator>
  EdgeOutputIterator
  dump_refinement_edges(EdgeOutputIterator out)
  {
    std::vector<Point_3> points(new_vertex_index);
    BOOST_FOREACH(Vertex_handle vh, vertices)
    points[get(m_vertex_ipmap,vh)] = vh->point();

    typedef std::vector<std::pair<Point_3,std::size_t> > Facet_point_vector;
    BOOST_FOREACH(const Facet_point_vector& facet_points, vertices_in_facets)
    {
      typedef std::pair<Point_3,std::size_t> Pair;
      BOOST_FOREACH(const Pair& pair, facet_points)
      points[pair.second]=pair.first;
    }

    std::size_t i=0;
    typedef std::vector<std::pair<Root_of_2,std::size_t> > Edge_point_vector;
    BOOST_FOREACH(const Edge_point_vector& edge_points, vertices_on_edges)
    {
      ++i;
      if (edge_points.empty()) continue;
      const Point_3 &source=edges[i-1]->opposite()->vertex()->point();
      const Point_3 &target=edges[i-1]->vertex()->point();
      typedef std::pair<Root_of_2,std::size_t> Pair;
      BOOST_FOREACH(const Pair& pair, edge_points)
      points[pair.second]=source + CGAL::to_double(pair.first) * (target-source);
    }

    typedef std::vector<std::pair<ID,ID> > Facet_edge_vector;
    BOOST_FOREACH(const Facet_edge_vector& facet_edges, edges_in_facets)
    {
      typedef std::pair<ID,ID> Pair;
      BOOST_FOREACH(const Pair& edge, facet_edges)
      *out++=std::pair<Point_3, Point_3>(points[edge.first], points[edge.second]);
    }

    return out;
  }

  void operator()(HDS &hds)
  {
    Decorator decorator(hds);

    //redimension the vertex vector
    vertices.resize(new_vertex_index);

    // insert the vertices in the facets
    std::size_t nbf=vertices_in_facets.size();
    for (std::size_t i=0; i<nbf; ++i)
    {
      std::size_t nbv=vertices_in_facets[i].size();
      for (std::size_t j=0; j<nbv; ++j)
      {
        std::pair<Point_3,std::size_t> pair=vertices_in_facets[i][j];
        vertices[ pair.second ] = decorator.vertices_push_back(
                                    typename HDS::Vertex(pair.first));
        put(m_vertex_ipmap, vertices[ pair.second ], pair.second);
      }
    }

    // insert the vertices on the edges
    std::size_t nbe=vertices_on_edges.size();
    for (std::size_t i=0; i<nbe; ++i)
    {
      if (vertices_on_edges[i].empty()) continue;
      std::size_t nbv=vertices_on_edges[i].size();
      Halfedge_handle hedge=edges[i];
      //sort the points along the halfedge
      std::sort(vertices_on_edges[i].begin(),
                vertices_on_edges[i].end(),
                internal::First_of_pair_less());

      std::vector<double> approx_alphas(nbv);
      for (std::size_t iv=0; iv!=nbv; ++iv)
      {
        approx_alphas[iv]=to_double(vertices_on_edges[i][iv].first);
        CGAL_assertion(iv==nbv-1 ||
                       vertices_on_edges[i][iv].first != vertices_on_edges[i][iv+1].first);
      }

      //we do another sort to ensure the alpha approximations are sorted
      //like they should be. Note that the values might no longer lie in their
      //original intervals.
      std::sort(approx_alphas.begin(), approx_alphas.end());

      // Ensure that alpha's are in ]0,1[
      //  --> First handle values smaller or equal to 0
      std::size_t k=0;
      while (k < nbv && approx_alphas[k]<=0) ++k;
      if (k!=0)
      {
        double previous = 0;
        for (std::size_t j=0; j<k; ++j)
        {
          approx_alphas[j] = boost::math::nextafter(previous,10);
          previous=approx_alphas[j];
        }
        while (k < nbv && approx_alphas[k]<=approx_alphas[k-1])
        {
          approx_alphas[k] = boost::math::nextafter(approx_alphas[k-1],10);
          ++k;
        }
      }
      //  --> then handle values greater or equal to 1
      k=nbv;
      while (k > 0 && approx_alphas[k-1]>=1) --k;

      if (k!=nbv)
      {
        double previous = 1;
        for (std::size_t j=nbv; j>k; --j)
        {
          approx_alphas[j-1] = boost::math::nextafter(previous,-10);
          previous=approx_alphas[j-1];
        }
        while (k > 0 && approx_alphas[k-1]>=approx_alphas[k])
        {
          approx_alphas[k-1] = boost::math::nextafter(approx_alphas[k],-10);
          --k;
        }
      }

      //check for equal values and update them in another desperate try to avoid pbs
      for (std::size_t iv=1; iv!=nbv; ++iv)
      {
        if (approx_alphas[iv]<=approx_alphas[iv-1])
          approx_alphas[iv] = boost::math::nextafter(approx_alphas[iv-1],10);
      }

      /// if we passed 1 then the update in the previous loop has gone wrong.
      if (approx_alphas.back()>=1)
      {
        /// \todo Question: Do you want me to implement snap-rounding?
        /// another solution is to sample the double to remove equality but this will work
        /// provided enough double exists
        /// I could always find a set of rational alpha that matches the order
        /// but would it be meaningful for the double embedding for your application?
        std::cerr<< "Two intersection points are too close on an edge, you need snap-rounding!\n";
        CGAL_error();
      }

      CGAL_assertion(Surface_mesh_overlay_3::canonical_edge(hedge)==hedge);
      typename EK::Point_3 b0 = m_exact_ppmap[get(m_vertex_ipmap, hedge->opposite()->vertex())];
      typename EK::Point_3 b1 = m_exact_ppmap[get(m_vertex_ipmap, hedge->vertex())];
      typename EK::Vector_3 b0b1 = b1 - b0;

      //set the corresponding rational approximate points
      for (std::size_t iv=0; iv!=nbv; ++iv)
      {
        CGAL_assertion(approx_alphas[iv]>0);
        CGAL_assertion(approx_alphas[iv]<1);

        m_exact_ppmap[ vertices_on_edges[i][iv].second ] =
          b0 + typename EK::FT(approx_alphas[iv]) * b0b1;
      }

      for (std::size_t j=0; j<nbv; ++j)
      {
        std::size_t vertex_index=vertices_on_edges[i][j].second;
        /// \todo replace by split_edge!!
        Vertex_handle v = decorator.vertices_push_back(
                            typename HDS::Vertex(to_input(m_exact_ppmap[vertex_index])));
        vertices[ vertex_index ] = v;
        put(m_vertex_ipmap, v, vertex_index);

        //   new_hedge    hedge
        //  ----------->   ----------->
        //               v
        //  <-----------   <-----------
        //   new_opposite     opposite
        //

        Halfedge_handle opposite=hedge->opposite();

        Halfedge_handle new_hedge=hds.edges_push_back(*hedge);
        Halfedge_handle new_opposite=new_hedge->opposite();

        //update next-prev relations
        new_hedge->HBase::set_next(hedge);
        new_hedge->prev()->HBase::set_next(new_hedge);
        decorator.set_prev(hedge, new_hedge);

        opposite->HBase::set_next(new_opposite);
        decorator.set_prev(new_opposite, opposite);
        decorator.set_prev(new_opposite->next(),new_opposite);

        //update vertex-hedge relation
        decorator.set_vertex(opposite,v);
        decorator.set_vertex(new_hedge,v);
        decorator.set_vertex_halfedge(v,new_hedge);
        decorator.set_vertex_halfedge(new_opposite->vertex(),new_opposite);
      }
    }

    // update the P2 halfedge corresponding when a edge of P2 is totally or
    // partially included into an edge of P1.
    // At this point they are subdivided, thus we can update the correspondance
    BOOST_FOREACH(const Vertex_ids_and_halfedge& iih,
                  hedges_of_P1_covered_by_hedge_of_P2)
    {
      Vertex_handle v1 = vertices[ cpp11::get<0>(iih)];
      Vertex_handle v2 = vertices[ cpp11::get<1>(iih)];
      Halfedge_handle h_P2 = cpp11::get<2>(iih);

      CGAL_assertion(cpp11::get<0>(iih) < vertices.size());
      CGAL_assertion(cpp11::get<1>(iih) < vertices.size());
      CGAL_assertion(v1!=Vertex_handle());
      CGAL_assertion(v2!=Vertex_handle());

      // recover the halfedge v1->v2 in P1 to associate it h from P2
      Halfedge_handle h_P1 = v2->halfedge();
      while (h_P1->opposite()->vertex() != v1)
      {
        h_P1=h_P1->next()->opposite();
        CGAL_assertion(h_P1!=v2->halfedge());
      }
      // do the association
      m_P2_edge_to_P1_edges[h_P2->id()].push_back(h_P1);
      m_P2_edge_to_P1_edges[h_P2->opposite()->id()].push_back(h_P1->opposite());
    }

/// \todo comment l'overlap partiel d'arete est influence par handle_intersection_of_constraints?
///       il faudrait probablement mettre  jour les paires en function des intersections...
///       un truc simple et de regarder si on trouve pas une arete, alors a veut dire qu'un
///       point a  ajout et donc il faut ajouter les points au milieu (mais il faut appeler
///       cette fonction la ligne d'avant)
    if (AllowPathIntersections) handle_intersection_of_constraints(decorator);

    for (std::size_t i=0; i<nbf; ++i)
      /// \todo instead of using exact points to do the cdt, we could use the
      ///       original mesh (P2) connectivity to get how the new edges
      ///       should be. However, we would still need to triangulate the
      ///       facets after
      create_interior_hedges_and_triangulate_facet(hds, decorator, i);


    // now handle facets that has no projection
    // create the vertices
    typedef std::pair<Vertex_handle, std::size_t> Vertex_and_id;
    std::map<Vertex_handle, std::size_t> vid_P2;
    BOOST_FOREACH(const Vertex_and_id& v_i, vertices_with_no_projection)
    {
      vertices[v_i.second]=decorator.vertices_push_back(
                             typename HDS::Vertex(v_i.first->point()));
      vid_P2.insert(std::make_pair(v_i.first, v_i.second));
      put(m_vertex_ipmap, vertices[v_i.second], v_i.second);
    }

    typedef std::pair<std::size_t,std::size_t> Indices;
    typedef std::map<Indices, Halfedge_handle> Edge_map;
    Edge_map edge_map; // edges from P1 vertex indices
    typedef std::map<std::size_t, Halfedge_handle> Src_id_to_incident_hedge;
    std::vector< Src_id_to_incident_hedge > halfedges_per_ccb(nb_faces_with_no_full_projections); //map src id -> halfedge per ccb

    std::vector<bool> is_border_ccb(nb_faces_with_no_full_projections, false);

    // create the halfedges in P1, and set the association with the parent
    // halfedge from P2. Vertices are also associated to the halfedge
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
    std::ofstream edges_outside_domain("debug/edges_outside_domain.cgal");
#endif
    BOOST_FOREACH(const Hedge_and_P1_vid& tpl, halfedges_with_parts_outside_P1)
    {
      // get P2 halfedge and endpoint vertex indices in P1
      Halfedge_handle h_P2=get<0>(tpl);
      Indices indices(get<1>(tpl), get<2>(tpl));
      if (indices.first==max_id) indices.first=vid_P2[h_P2->opposite()->vertex()];
      if (indices.second==max_id) indices.second=vid_P2[h_P2->vertex()];
      //create the halfedge if it not already done
      std::pair<typename Edge_map::iterator, bool> insert_res =
        edge_map.insert(std::make_pair(indices,Halfedge_handle()));
      if (insert_res.second)
      {
        Halfedge_handle h=hds.edges_push_back(Halfedge(), Halfedge());
        insert_res.first->second=h;
        edge_map.insert(std::make_pair(Indices(indices.second, indices.first),
                                       h->opposite()));
        //update vertex<->halfedge pointers
        decorator.set_vertex(h,vertices[indices.second]);
        decorator.set_vertex(h->opposite(),vertices[indices.first]);
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
        edges_outside_domain << "2 " << h->opposite()->vertex()->point()
                             << " " << h->vertex()->point() << "\n";
#endif
        if (vertices[indices.second]->halfedge()==Halfedge_handle())
          decorator.set_vertex_halfedge(vertices[indices.second],h);
        if (vertices[indices.first]->halfedge()==Halfedge_handle())
          decorator.set_vertex_halfedge(vertices[indices.first],h->opposite());
        // associate the newly created halfedges in P1 to its parent in P2
        m_P2_edge_to_P1_edges[h_P2->id()].push_back(h);
        m_P2_edge_to_P1_edges[h_P2->opposite()->id()].push_back(h->opposite());
      }
      // insert the halfedges in the per ccb map
      std::size_t ccb_id=ccb_ids_in_P2[h_P2];
      if (h_P2->is_border()) is_border_ccb[ccb_id]=true;
      halfedges_per_ccb[ ccb_id ].insert(
        std::make_pair(indices.first, insert_res.first->second));
      ccb_id=ccb_ids_in_P2[h_P2->opposite()];
      if (h_P2->opposite()->is_border()) is_border_ccb[ccb_id]=true;
      halfedges_per_ccb[ ccb_id ].insert(
        std::make_pair(indices.second, insert_res.first->second->opposite()));
    }
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
    edges_outside_domain.close();
#endif

#ifdef CGAL_DEBUG_OVERLAY_OUTSIDE_CCBS
    std::size_t nb_ccb=halfedges_per_ccb.size();
    for (std::size_t i=0; i<nb_ccb; ++i)
    {
      std::stringstream sstr;
      sstr << "debug/ccb-" << i << ".cgal";
      std::ofstream ouput(sstr.str().c_str());
      typedef std::pair<std::size_t,Halfedge_handle> Pair_type;
      BOOST_FOREACH(Pair_type p,halfedges_per_ccb[ i ])
      ouput << "2 "<< p.second->opposite()->vertex()->point() << " " << p.second->vertex()->point() << "\n";
      ouput.close();
    }
#endif

    // set face pointer for each halfedge of each ccb and collect how
    // next/prev pointers should be set
    typedef std::pair<Halfedge_handle, Halfedge_handle> Hedge_pair;
    std::vector< Hedge_pair > prev_next;
    std::vector< std::size_t > faces_to_triangulate;
    prev_next.reserve(halfedges_with_parts_outside_P1.size());
    std::size_t ccb_id=0;
    BOOST_FOREACH(Src_id_to_incident_hedge& map, halfedges_per_ccb)
    {
      while (!map.empty())
      {
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
        std::cout << "Handling ccb " << ccb_id << "\n";
        std::cout << "   map.size() " << map.size() << "\n";
        std::cout << "   is_border_ccb[ccb_id] " << is_border_ccb[ccb_id] << "\n";
#endif
        Facet_handle f;
        // create a new face only for non-border ccb
        if (!is_border_ccb[ccb_id])
        {
          f=decorator.faces_push_back(Facet());
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
          std::cout << "   new face " << &(*f) << "\n";
#endif
          f->id()=facets.size();
          facets.push_back(f);
        }
        std::size_t size_of_border=0;
        // we use the map associating a vertex id to the halfedge in the ccb
        // having that vertex as source. When the vertex is not one of the
        // source vertices of the ccb registered, it indicates that a border
        // halfedge from P1 should be used
        std::pair<std::size_t, Halfedge_handle> first=*map.begin();
        Halfedge_handle next=first.second; // starting point
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
        std::cout << "   starting point 2 " << next->opposite()->vertex()->point() << " " << next->vertex()->point() << "\n";
#endif
        if (f!=Facet_handle()) decorator.set_face_halfedge(f, next);
        do
        {
          if (f!=Facet_handle()) decorator.set_face(next, f);
          std::size_t tgt_index=get(m_vertex_ipmap,next->vertex());
          typename Src_id_to_incident_hedge::iterator find_res=map.find(tgt_index);
          if (find_res==map.end())
          {
            CGAL_assertion(next->vertex()->halfedge()!=Halfedge_handle());
            // look for a halfedge on the border of P1 with the correct src vertex
            Halfedge_handle candidate=next->vertex()->halfedge()->opposite();
            CGAL_assertion_code(Halfedge_handle first=candidate;)
            while (!candidate->is_border())
            {
              candidate=candidate->opposite()->next();
              CGAL_assertion(candidate!=first);
            }
            CGAL_assertion(get(m_vertex_ipmap,candidate->opposite()->vertex())==tgt_index);
            prev_next.push_back(Hedge_pair(next, candidate));
            next=candidate;
            ++size_of_border;
          }
          else
          {
            prev_next.push_back(Hedge_pair(next, find_res->second));
            map.erase(find_res); // remove the element handled
            next=find_res->second;
            ++size_of_border;
          }
        }
        while (next!=first.second);
        if (size_of_border!=3 && f!=Facet_handle())
          faces_to_triangulate.push_back(f->id());
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
        std::cout << "size_of_border " << size_of_border << "\n";
#endif
      }
      ++ccb_id;
    }
    //now update the pointers
    BOOST_FOREACH(const Hedge_pair& p, prev_next)
    {
      p.first->HBase::set_next(p.second);
      decorator.set_prev(p.second,p.first);
    }

    // triangulate faces outside P1 domain
    BOOST_FOREACH(std::size_t i, faces_to_triangulate)
    {
      try
      {
        create_interior_hedges_and_triangulate_facet(hds, decorator,i,true);
      }
      catch (CDT_TMP_DEBUG)
      {
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
        std::cout << "Falling back on the hole filling to triangulate a face\n";
#endif
        typedef CGAL::Triple<std::size_t, std::size_t, std::size_t> Triangle_int;
        std::vector<Triangle_int> triangles;
        std::vector<Point_3> polyline;
        std::vector<Point_3> third_points;
        std::vector<Halfedge_handle> hedges;
        Halfedge_handle first=facets[i]->halfedge(), h=first;
        do
        {
          polyline.push_back(h->vertex()->point());
          hedges.push_back(h);
          h=h->next();
          if (h->opposite()->is_border())
            third_points.push_back(
              midpoint(h->opposite()->vertex()->point(),h->vertex()->point())
            );
          else
            third_points.push_back(h->opposite()->next()->vertex()->point());
        }
        while (h!=first);
        CGAL::Polygon_mesh_processing::triangulate_hole_polyline(
          polyline, third_points,
          std::back_inserter(triangles));
        CGAL_assertion(!triangles.empty());
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
        std::cout << "   face with " << hedges.size() << " vertices decomposed in "
                  << triangles.size() << " triangles\n";
#endif
        typedef std::pair<std::size_t, std::size_t> Indices;
        std::set< Indices > diagonals;
        BOOST_FOREACH(Triangle_int& t, triangles)
        {
          std::size_t indices[3]= {t.first, t.second, t.third};
          // look for a diagonal
          std::size_t k=0;
          for (k=0; k<3; ++k)
            if ((indices[k]+1)%hedges.size()!=indices[(k+1)%3])
              diagonals.insert(
                make_sorted_pair(indices[k], indices[(k+1)%3]));
        }

        // indices are kept sorted so that the first halfedge of the pair
        // is always on the correct side. Only the second halfedge might
        // need an update
        BOOST_FOREACH(Indices p, diagonals)
        {
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
          std::cout << "Adding diagonal " << p.first << " " << p.second << "\n";
#endif
          Halfedge_handle h=hedges[p.first], g=hedges[p.second];
          while (h->face()!=g->face()) g=g->next()->opposite();
          CGAL_assertion(h->face()==g->face());
          CGAL_assertion_code(Facet_handle old_face = h->face();)
          decorator.split_face(h, g);
          CGAL_assertion(old_face==h->face());
          g->face()->id()=facets.size();
          facets.push_back(g->face());
          hedges[p.second]=g->next()->opposite();
        }
      }
    }
  }
};

} // end of namespace CGAL

#endif // CGAL_SURFACE_MESH_OVERLAY_3_POLYHEDRON_REFINEMENT_MODIFIER_3_H
