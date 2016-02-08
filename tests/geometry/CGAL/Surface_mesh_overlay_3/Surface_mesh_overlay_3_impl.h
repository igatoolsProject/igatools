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

#ifndef CGAL_SURFACE_MESH_OVERLAY_3_SURFACE_MESH_OVERLAY_3_H
#define CGAL_SURFACE_MESH_OVERLAY_3_SURFACE_MESH_OVERLAY_3_H

#include "Polyhedron_utils.h"
#include "Polyhedron_refinement_modifier_3.h"
#include "Constructions.h"
#include "../Walk_on_polyhedron.h"
#include "Walk_on_polyhedron_visitor.h"

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include "AABB_traits.h"
#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_halfedge_graph_segment_primitive.h>
#include "../Polyhedron_simplex_type.h"
#include <CGAL/property_map.h>
#include <CGAL/Default.h>
#include <CGAL/tuple.h>

#include <boost/foreach.hpp>

#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
#define CGAL_DEBUG_OVERLAY_3_REFINEMENT(M) std::cout << M;
#else
#include <iostream>
#define CGAL_DEBUG_OVERLAY_3_REFINEMENT(M)
#endif

namespace CGAL
{

namespace Surface_mesh_overlay_3
{

template <class Vertex_handle>
struct Default_vertex_index_pmap
{
  //classical typedefs
  typedef Vertex_handle key_type;
  typedef std::size_t value_type;
  typedef value_type reference;
  typedef boost::read_write_property_map_tag category;

  friend
  reference get(Default_vertex_index_pmap, key_type v)
  {
    return v->id();
  }

  friend
  void put(Default_vertex_index_pmap, key_type v, value_type k)
  {
    v->id()=k;
  }
};

template <class Point_3>
class Vector_ppmap
{
  std::vector<Point_3> *points_ptr;
public:
  //classical typedefs
  typedef typename std::vector<Point_3>::size_type key_type;
  typedef Point_3 value_type;
  typedef value_type &reference;
  typedef boost::lvalue_property_map_tag category;
  Vector_ppmap(std::vector<Point_3> &points_):points_ptr(&points_) {}

  reference operator[](std::size_t i)
  {
    CGAL_assertion(i<points_ptr->size());
    return (*points_ptr)[i];
  }

  friend
  const reference get(const Vector_ppmap &vmap, key_type i)
  {
    CGAL_assertion(i<vmap.points_ptr->size());
    return (*vmap.points_ptr)[i];
  }

  friend
  void put(const Vector_ppmap &vmap, key_type i, const Point_3 &pt)
  {
    CGAL_assertion(i==vmap.points_ptr->size());
    vmap.points_ptr->push_back(pt);
  }
};

template<class Polyhedron>
struct No_edge_skipt
{
  typedef typename Polyhedron::Edge_iterator Iterator;
  /// returns true if the edge must be considered
  bool operator()(Iterator) const
  {
    return false;
  }
};

template <class Input_kernel,
          class Exact_kernel,
          class Polyhedron,
          class ES=Default,
          class VIPM=Default
          >
struct Surface_mesh_overlay_3_impl
{
  typedef Input_kernel K;
  typedef Exact_kernel EK;

  typedef typename Default::Get<ES,
          No_edge_skipt<Polyhedron> >::type   EdgeSkipper;
  typedef typename Default::Get<
  VIPM,
  Default_vertex_index_pmap<
  typename Polyhedron::Vertex_handle>
  >::type                               VertexIndexPropertyMap;

  typedef typename Input_kernel::Point_3                            Input_point;
  typedef typename Input_kernel::Vector_3                          Input_vector;
  typedef typename Exact_kernel::Segment_3                            Segment_3;
  typedef typename Exact_kernel::Point_3                                Point_3;
  typedef typename Exact_kernel::Vector_3                              Vector_3;
  typedef typename Root_of_traits<typename EK::FT>::Root_of_2         Root_of_2;

  typedef typename Polyhedron::Halfedge_handle                  Halfedge_handle;
  typedef typename Polyhedron::Halfedge_handle                               HH;
  typedef typename Polyhedron::Facet_handle                        Facet_handle;
  typedef typename Polyhedron::Vertex_handle                      Vertex_handle;

  typedef AABB_face_graph_triangle_primitive<Polyhedron>              Primitive;
  typedef Surface_mesh_overlay_3::AABB_traits<Input_kernel, Primitive>   Traits;
  typedef AABB_tree<Traits>                                                Tree;

  // typedef for the edge intersection
  typedef Edge_intersection<Polyhedron,
          Exact_kernel,
          VertexIndexPropertyMap>      Intersection_predicate;
  typedef typename Intersection_predicate::Barycentric_NT        Barycentric_NT;
  typedef typename Intersection_predicate::result_type                Inter_res;

  //
  typedef Vector_ppmap<typename EK::Point_3>                        Exact_ppmap;
  typedef Polyhedron_refinement_modifier_3< Polyhedron,
          EK, K,
          Exact_ppmap,
          VertexIndexPropertyMap>    Modifier;
  typedef Walk_on_polyhedron_visitor<Polyhedron, Modifier, Barycentric_NT,
          VertexIndexPropertyMap>       Walk_visitor;
  typedef Walk_on_polyhedron<Polyhedron,Intersection_predicate,
          Walk_visitor>                               Walker;
  /// contains the intersection type, the vertex id and a halfedge
  /// corresponding to the simplex the point is on.
  typedef CGAL::cpp11::tuple<Polyhedron_simplex_type,
          std::size_t,
          Halfedge_handle >                      Inter_info;

  CGAL::Cartesian_converter<K,EK> to_exact;
  CGAL::Cartesian_converter<EK,K> to_input;

  VertexIndexPropertyMap m_vertex_ipmap;

  Surface_mesh_overlay_3_impl(VertexIndexPropertyMap vipmap =
                                VertexIndexPropertyMap())
    :m_vertex_ipmap(vipmap) {}

  void get_incident_facets(Polyhedron_simplex_type inter_type,
                           Halfedge_handle h,
                           std::set<Facet_handle> &facets)
  {
    switch (inter_type)
    {
      case POLYHEDRON_FACET:
        facets.insert(h->facet());
        break;
      case POLYHEDRON_EDGE:
        facets.insert(h->facet());
        facets.insert(h->opposite()->facet());
        break;
      case POLYHEDRON_VERTEX:
      {
        Halfedge_handle start=h;
        do
        {
          facets.insert(h->facet());
          h=h->next()->opposite();
        }
        while (start!=h);
      }
      break;
      case POLYHEDRON_NONE:
        CGAL_assertion(!"Should never be here");
        break;
    }
  }

  //mostly for debug: returns (1-alpha) * source + alpha * target
  Point_3
  compute_point_on_edge(Halfedge_handle h, const typename EK::FT &alpha)
  {
    Vector_3 p0 = to_exact(h->opposite()->vertex()->point()) - CGAL::ORIGIN;
    Vector_3 p1 = to_exact(h->vertex()->point()) - CGAL::ORIGIN;
    return CGAL::ORIGIN + (1-alpha) * p0 + alpha * p1;
  }

  // compute the projection of a point onto P1 using a pair of segments
  // attached to `pt`.
  // the projection vector is also updated the match the normal used
  // to compute the intersection point
  cpp11::tuple<Halfedge_handle, Polyhedron_simplex_type, Point_3>
  compute_projection(const Input_point &input_pt,
                     Vector_3 &projection_normal,
                     double min_half_diameter,
                     const Tree &tree)
  {
    /// the query type is not a segment but a pair of segments to avoid
    /// missing the projection point while defining the inexact segment
    Surface_mesh_overlay_3::Query_type<K> query;
    const Point_3 exact_pt = to_exact(input_pt);
    const Point_3 target1 =  exact_pt - min_half_diameter*projection_normal;
    const Point_3 target2 =  exact_pt + min_half_diameter*projection_normal;

    query.s1 = typename Input_kernel::Segment_3(input_pt, to_input(target1));
    query.s2 = typename Input_kernel::Segment_3(input_pt, to_input(target2));

    boost::optional<typename Primitive::Id> primitive =
      tree.any_intersected_primitive(query);

    // nothing to do if no intersection was found
    if (primitive==boost::none)
      return cpp11::make_tuple(Halfedge_handle(), POLYHEDRON_NONE, Point_3());

    CGAL_assertion(tree.traits().last_segment_intersected!=0);
    bool is_on_s1 = tree.traits().last_segment_intersected==1;
    // update the target point to match the one used in the AABB-tree
    const Point_3 &target = to_exact(is_on_s1 ? query.s1.target()
                                     : query.s2.target());
    // look at the intersection type
    std::pair<Halfedge_handle, Polyhedron_simplex_type>  res =
      intersection_type<Polyhedron, Exact_kernel>(*primitive, exact_pt, target);
    if (res.second==POLYHEDRON_VERTEX)
      return cpp11::make_tuple(res.first, res.second, Point_3());

    CGAL_assertion(res.second==POLYHEDRON_EDGE || res.second==POLYHEDRON_FACET);

    ///  \todo do snap-rounding in the triangle and snap interpt on edge or vertex
    ///        also don't forget to update the normal vertor. Also do this in edge_case
    ///        but be aware that you cannot snap if the point before projection
    ///        is already on the edge or the face

    //exact segment used to compute the intersection point
    typename EK::Segment_3 s(exact_pt, target);
    const Point_3 &interpt  = compute_inter_pt<EK,Polyhedron>(res.first,s);
    if (interpt!=s.source())
    {
      //define the new normal that is the one used to compute the intersection
      //the orientation of the normal must be preserved (thus the `one` trick)
      typename EK::FT one = is_on_s1 ? -1:1;
      typename EK::Vector_3 new_normal(s.source(), s.target());
      projection_normal = (one / min_half_diameter) * new_normal;
    }
    return cpp11::make_tuple(res.first, res.second, interpt);
  }

  // compute the projection of a point onto P1 using a line passing though
  // `pt`, `projection_normal` as direction
  // `projection_normal` is not updated as the original normal vector was
  // computed using the input kernel and the fact that we use a line
  // makes that the intersection predicates are exact
  cpp11::tuple<Halfedge_handle, Polyhedron_simplex_type, Point_3>
  compute_projection(const Input_point &input_pt,
                     const Vector_3 &projection_normal,
                     const Tree &tree)
  {
    typename Input_kernel::Line_3 query(input_pt,to_input(projection_normal));

    std::set<typename Primitive::Id> primitives;
    tree.all_intersected_primitives(query, std::inserter(primitives, primitives.begin()));

    // look for the closest intersection point
    cpp11::tuple<Halfedge_handle, Polyhedron_simplex_type, Point_3>  result =
      cpp11::make_tuple(Halfedge_handle(), POLYHEDRON_NONE, Point_3());
    const Point_3 exact_pt = to_exact(input_pt);
    const typename EK::Line_3 exact_line(exact_pt, projection_normal);
    const Point_3 p0 = exact_line.point(0);
    const Point_3 p1 = exact_line.point(1);

    while (!primitives.empty())
    {
      typename Primitive::Id fh = *primitives.begin();
      primitives.erase(*primitives.begin());
      // look at the intersection type
      std::pair<Halfedge_handle, Polyhedron_simplex_type>  itype =
        intersection_type<Polyhedron, Exact_kernel>(fh, p0, p1);
      // remove faces which intersection point will be found several times
      if (itype.second==POLYHEDRON_EDGE)
        primitives.erase(itype.first->opposite()->face());
      else if (itype.second==POLYHEDRON_VERTEX)
      {
        Halfedge_handle h=itype.first;
        do
        {
          primitives.erase(h->face());
          h=h->next()->opposite();
        }
        while (h!=itype.first);
      }
      Point_3 interpt;

      if (itype.second==POLYHEDRON_EDGE || itype.second==POLYHEDRON_FACET)
        interpt=compute_inter_pt<EK,Polyhedron>(itype.first,to_exact(query));
      else
      {
        CGAL_assertion(itype.second==POLYHEDRON_VERTEX);
        interpt=to_exact(itype.first->vertex()->point());
      }
      if (cpp11::get<1>(result)==POLYHEDRON_NONE ||
          CGAL::compare_distance_to_point(exact_pt, interpt, cpp11::get<2>(result))==CGAL::SMALLER)
        result = cpp11::make_tuple(itype.first, itype.second, interpt);
    }

    return result;
  }

  struct Less_edge_intersection_alphas
  {
    bool operator()(const std::pair<Halfedge_handle, Inter_res> &p1,
                    const std::pair<Halfedge_handle, Inter_res> &p2)
    {
      return p1.second->first < p2.second->first;
    }
  };

  template<class OutputIterator, class HalfedgeRange>
  OutputIterator
  compute_intersections(Halfedge_handle h_p2,
                        const HalfedgeRange &hedges_p1,
                        double min_half_diameter,
                        const std::vector<Vector_3> normals,
                        OutputIterator out)
  {

    Intersection_predicate intersection(normals, h_p2, m_vertex_ipmap);

/// \todo use the AABB-tree!
    std::set<Halfedge_handle> hedge_to_skip;
    BOOST_FOREACH(Halfedge_handle h_p1, hedges_p1)
    {
      if (hedge_to_skip.count(h_p1)) continue; // intersection at a vertex
      // already detected
      Inter_res inter_res = intersection(h_p1);


      if (inter_res!=boost::none)
      {
        // filter intersection that are too far
        // \todo this is a delicate part that is not good yet (threshold choice)
        double beta=to_double(inter_res->first), alpha=0;

        // in case we hit a vertex we need to filter out the other
        // adjacent border edge
        const bool *is_target_ptr=boost::get<bool>(&(inter_res->second));
        if (is_target_ptr)
        {
          hedge_to_skip.insert(
            canonical_edge(h_p1->is_border()
                           ? (*is_target_ptr?h_p1->next():h_p1->prev())
                           : (*is_target_ptr ? h_p1->opposite()->prev()
                              : h_p1->opposite()->next())
                          ));
          if (*is_target_ptr) alpha=1;
        }
        else
          alpha=to_double(*boost::get<Barycentric_NT>(&(inter_res->second)));

        const Point_3 &p1=compute_point_on_edge(h_p1,alpha);
        const Point_3 &p2=compute_point_on_edge(h_p2,beta);

        double threshold = min_half_diameter>0
                           ? min_half_diameter*min_half_diameter
                           : squared_distance(h_p1->opposite()->vertex()->point(),
                                              h_p1->vertex()->point());
        if (squared_distance(p1,p2) <= threshold)
          *out++=std::make_pair(h_p1, inter_res);
      }
    }
    return out;
  }

  void partial_edge_overlay(
    Halfedge_handle h_p2,
    const std::vector< std::pair<Halfedge_handle, Inter_res> > &intersections,
    std::size_t prev_id,
    std::size_t final_id,
    Polyhedron &P1,
    const std::vector<Vector_3> &normals,
    Modifier &modifier)
  {
    const std::size_t max_id = (std::numeric_limits<std::size_t>::max)();
    EdgeSkipper skip;

    std::size_t i=0;
    std::size_t src_id=max_id;
    Polyhedron_simplex_type src_inter_type=POLYHEDRON_EDGE;
    do
    {
      //handle tangency at the end
      if (intersections.size()==i+1)
      {
        Halfedge_handle h_p1=intersections.back().first;
        const Inter_res &inter_res = intersections.back().second;
        const Barycentric_NT *barycentric_coord =
          boost::get<Barycentric_NT>(&(inter_res->second));

        if (!barycentric_coord)
        {
          // a vertex was reached
          bool is_target = boost::get<bool>(inter_res->second);
          std::size_t id = get(m_vertex_ipmap, is_target?h_p1->vertex():h_p1->opposite()->vertex());
          modifier.add_edges_with_partial_projection(h_p2, prev_id, id, skip);
          prev_id=id;
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
          std::cout << "handle a vertex tangency\n";
#endif
        }
        else
          CGAL_assertion(!"Should not get here, pb filtering intersections?");
        break;
      }

      const std::pair<Halfedge_handle, Inter_res> &src=intersections[i];
      const std::pair<Halfedge_handle, Inter_res> &tgt=intersections[i+1];

      std::size_t tgt_id=max_id;
      Halfedge_handle h_src=src.first;
      Halfedge_handle h_tgt=tgt.first;
      Polyhedron_simplex_type tgt_inter_type=POLYHEDRON_EDGE;
      // source (smallest intersection point)
      if (src_id==max_id)
      {
        const Inter_res &inter_src = src.second;
        const Barycentric_NT *src_barycentric_coord =
          boost::get<Barycentric_NT>(&(inter_src->second));

        if (!src_barycentric_coord)
        {
          // a vertex was reached
          bool is_target = boost::get<bool>(inter_src->second);
          src_id = get(m_vertex_ipmap, is_target?h_src->vertex()
                       :h_src->opposite()->vertex());
          src_inter_type=POLYHEDRON_VERTEX;
          if (!is_target) h_src=h_src->opposite();
        }
        else
        {
          src_id = modifier.add_point_on_edge(h_src, *src_barycentric_coord);
          src_inter_type=POLYHEDRON_EDGE;
        }
      }
      // target (largest intersection point)
      const Inter_res &inter_tgt = tgt.second;
      const Barycentric_NT *tgt_barycentric_coord =
        boost::get<Barycentric_NT>(&(inter_tgt->second));

      if (!tgt_barycentric_coord)
      {
        // a vertex was reached
        bool is_target = boost::get<bool>(inter_tgt->second);
        tgt_id = get(m_vertex_ipmap, is_target?h_tgt->vertex()
                     :h_tgt->opposite()->vertex());
        tgt_inter_type=POLYHEDRON_VERTEX;
        if (!is_target) h_tgt=h_tgt->opposite();
      }
      else
        tgt_id = modifier.add_point_on_edge(h_tgt, *tgt_barycentric_coord);

      // register a part of `h_p2` as outside the domain
      modifier.add_edges_with_partial_projection(h_p2, prev_id, src_id, skip);

      Intersection_predicate predicate(normals, h_p2, m_vertex_ipmap);
      Walker walker;
      Walk_visitor walk_visitor(src_id, tgt_id,
                                h_src, h_tgt, tgt_inter_type,
                                h_p2, modifier, m_vertex_ipmap);

      bool target_reached = walker(P1, h_src, src_inter_type,
                                   h_tgt, tgt_inter_type,
                                   predicate, walk_visitor);

      if (!target_reached)
      {
        //tangency case
        CGAL_assertion(src_inter_type==POLYHEDRON_VERTEX);
        ++i;
        prev_id=src_id;
        // avoid recomputing everything for the next src
        src_id=tgt_id;
        src_inter_type=tgt_inter_type;
      }
      else
      {
        i+=2;
        prev_id=tgt_id;
        src_id=max_id;
      }
    }
    while (i<intersections.size());
    modifier.add_edges_with_partial_projection(h_p2, prev_id, final_id, skip);
  }

  void run(Polyhedron &P1,
           Polyhedron &P2,
           double min_half_diameter,
           const std::vector< Vertex_handle> &vertices_of_P2_already_matching,
           const std::vector< std::pair<Halfedge_handle, double> > &simplices_of_P1_already_matching,
           bool allow_partial_overlay,
           std::vector< std::pair<std::size_t, std::size_t> > &input_facet_ids)
  {
    /// set a unique id for all halfedges of P2
    std::size_t hedge_index_P2=0;
    for (typename Polyhedron::Halfedge_iterator hedge_it=P2.halfedges_begin(),
         hedge_end=P2.halfedges_end();
         hedge_it!=hedge_end; ++hedge_it
        ) hedge_it->id()=hedge_index_P2++;

    /// set a unique id for all vertices of P2
    std::size_t v_index_P2=0;
    for (typename Polyhedron::Vertex_iterator vit=P2.vertices_begin(),
         vit_end=P2.vertices_end();
         vit!=vit_end; ++vit, ++v_index_P2
        ) put(m_vertex_ipmap, vit, v_index_P2);

    std::vector<typename EK::Point_3> exact_points;
    Modifier modifier(hedge_index_P2, Exact_ppmap(exact_points), m_vertex_ipmap);
    set_halfedgeds_items_id(P1,
                            modifier.facets,
                            modifier.edges,
                            modifier.vertices,
                            m_vertex_ipmap
                           );

    /// init the exact version of the points
    exact_points.reserve(modifier.vertices.size());
    for (typename Polyhedron::Vertex_iterator vit = P1.vertices_begin(),
         evit = P1.vertices_end();
         vit != evit ; ++  vit)
    {
      exact_points.push_back(to_exact(vit->point()));
    }

    modifier.init();

    Tree tree(P1.facets_begin(),P1.facets_end(),P1);

    //store the projection normal
    std::vector<Vector_3> normals(v_index_P2);

    // store the projected points and the simplex of P1 on which they are
    std::vector< Inter_info > inter_infos(v_index_P2);

    /// simple map to catch identical points
    std::map<Input_point, Vertex_handle> p1_point_to_vertex;
    for (typename Polyhedron::Vertex_iterator pit=P1.vertices_begin(),
         pit_end=P1.vertices_end();
         pit!=pit_end; ++pit)
    {
      p1_point_to_vertex[pit->point()]=pit;
    }

#ifndef CGAL_NO_DEBUG_OVERLAY_3_PROJECTED_POINTS
    std::ofstream output_proj_segments("debug/projection_vectors.cgal");
    std::ofstream output_no_proj_segments("debug/projection_vectors_none.cgal");
    std::ofstream output_edge_proj_segments("debug/projection_vectors_edge.cgal");
    std::ofstream output_face_proj_segments("debug/projection_vectors_face.cgal");
    std::ofstream output_vertex_proj_segments("debug/projection_vectors_vertex.cgal");
#endif
    namespace PMP = CGAL::Polygon_mesh_processing;

    /// Init the projections of P2 vertices with those given as matching.
    /// The normal vector of a vertex is the mean of the normal vectors of facets incident to
    /// the simplex the vertex falls on (to avoid having orthogonal vectors while walking
    /// along an edge)
    std::vector<bool> vertices_handled(v_index_P2, false);
    CGAL_assertion(vertices_of_P2_already_matching.size()==simplices_of_P1_already_matching.size());
    std::size_t nb_v_matching=vertices_of_P2_already_matching.size();
    for (std::size_t i=0; i< nb_v_matching; ++i)
    {
      Vertex_handle vh=vertices_of_P2_already_matching[i];
      Halfedge_handle h;
      double alpha;
      boost::tie(h, alpha)=simplices_of_P1_already_matching[i];
      CGAL_assertion(h==canonical_edge(h));
      v_index_P2=get(m_vertex_ipmap, vh);

      if (alpha==0)
      {
        std::size_t v_index_P1 = get(m_vertex_ipmap, h->opposite()->vertex());
        inter_infos[v_index_P2] = Inter_info(POLYHEDRON_VERTEX, v_index_P1, h->opposite());

        normals[v_index_P2] = to_exact(PMP::compute_vertex_normal(h->opposite()->vertex(),P2));
      }
      else if (alpha==1)
      {
        std::size_t v_index_P1 = get(m_vertex_ipmap, h->vertex());
        inter_infos[v_index_P2] = Inter_info(POLYHEDRON_VERTEX, v_index_P1, h);
        normals[v_index_P2] = to_exact(PMP::compute_vertex_normal(h->vertex(), P2));
      }
      else
      {
        CGAL_assertion(alpha > 0 && alpha < 1);
        typename EK::FT exact_alpha = alpha;
        std::size_t vertex_index = modifier.add_point_on_edge(h, exact_alpha);
        inter_infos[v_index_P2] = Inter_info(POLYHEDRON_EDGE, vertex_index, h);
        /// \todo TAG_SL_GOOD_PROJECTION: we have a check to ensure the projection using the AABB-tree
        ///       is correct. It it is not the case, points on feature edges might not be correctly sorted
        ///       and some alpha permutation might need to be done. Note that the issue might go behind
        ///       edges (order swapped in adjacent edges for example)
        /// \todo exact point on edge computed several times. Store it ?
        Input_vector normal = CGAL::NULL_VECTOR;
        if (!h->is_border()) normal = normal + PMP::compute_face_normal(h->facet(),P2);
        if (!h->opposite()->is_border()) normal = normal + PMP::compute_face_normal(h->opposite()->facet(),P2);
        if (h->is_border()==h->opposite()->is_border()) normal = normal / 2;
        normals[v_index_P2] = to_exact(normal);
      }

      // in case the vertex was exactly on an edge or on a vertex
      if (normals[v_index_P2] == CGAL::NULL_VECTOR)
        normals[v_index_P2] = to_exact(PMP::compute_vertex_normal(vh,P2));
      else
      {
        // make the normal vector larger to avoid having a norm too small
        typename EK::FT x = normals[v_index_P2][0],
                        y = normals[v_index_P2][1],
                        z = normals[v_index_P2][2];
        if (x != 0)
          normals[v_index_P2] = Vector_3(1, y/x, z/x);
        else if (y != 0)
          normals[v_index_P2] = Vector_3(x/y, 1, z/y);
        else
          normals[v_index_P2] = Vector_3(x/z, y/z, 1);

        // make the normal almost unitary
        double norm = sqrt(CGAL::to_double(normals[v_index_P2].squared_length()));
        normals[v_index_P2] = normals[v_index_P2] / norm;
        // make sure the normal points in the correct direction
        Vector_3 surface_normal = to_exact(PMP::compute_vertex_normal(vh,P2));
        if (normals[v_index_P2] * surface_normal < 0)
          normals[v_index_P2] = - normals[v_index_P2];
        CGAL_assertion(normals[v_index_P2] * surface_normal > 0);
      }

      //mark the vertex as handled
      vertices_handled[v_index_P2]=true;
    }

    /// compute the projection of vertices
    for (typename Polyhedron::Vertex_iterator pit=P2.vertices_begin(),
         pit_end=P2.vertices_end();
         pit!=pit_end; ++pit)
    {
      v_index_P2 = get(m_vertex_ipmap,pit);
      if (vertices_handled[v_index_P2]) continue;   // already handled

      //compute_vertex_normal returns a normalized vector using doubles
      normals[v_index_P2] =
        to_exact(PMP::compute_vertex_normal(pit,P2));

#ifndef CGAL_NO_DEBUG_OVERLAY_3_PROJECTED_POINTS
      const Point_3 exact_pit_point = to_exact(pit->point());
      const Point_3 target1 =  exact_pit_point - min_half_diameter*normals[v_index_P2];
      const Point_3 target2 =  exact_pit_point + min_half_diameter*normals[v_index_P2];
      output_proj_segments << "2 " << target1 << " " << target2 << "\n";
#endif

      //check for exact common point first
      typename std::map<Input_point, Vertex_handle>::iterator it_map=
        p1_point_to_vertex.find(pit->point());
      if (it_map!=p1_point_to_vertex.end())
      {
        inter_infos[v_index_P2] =
          Inter_info(POLYHEDRON_VERTEX,
                     get(m_vertex_ipmap, it_map->second),
                     it_map->second->halfedge());
#ifndef CGAL_NO_DEBUG_OVERLAY_3_PROJECTED_POINTS
        output_vertex_proj_segments << "2 " << target1 << " " << target2 << "\n";
#endif
        continue;
      }

      // compute the projection point onto P1, and only in the segment case
      // update the projection normal to match the one used to compute
      // the projection (in the other case, no update is needed)
      cpp11::tuple<Halfedge_handle, Polyhedron_simplex_type, Point_3> inter_res
        = min_half_diameter > 0
          ? compute_projection(pit->point(),
                               normals[v_index_P2],
                               min_half_diameter,
                               tree) // use a pair of segments
          : compute_projection(pit->point(),
                               normals[v_index_P2],
                               tree); // use a line

      switch (cpp11::get<1>(inter_res))
      {
        case POLYHEDRON_FACET:
        {
          Halfedge_handle h = cpp11::get<0>(inter_res);
          exact_points.push_back(cpp11::get<2>(inter_res));

          std::size_t vertex_index =
            modifier.add_point_in_facet(h->facet(), to_input(cpp11::get<2>(inter_res)));
          inter_infos[v_index_P2] = Inter_info(POLYHEDRON_FACET, vertex_index, h);
#ifndef CGAL_NO_DEBUG_OVERLAY_3_PROJECTED_POINTS
          output_face_proj_segments << "2 " << target1 << " " << target2 << "\n";
#endif
        }
        break;
        case POLYHEDRON_EDGE: //this is an edge
        {
#ifndef CGAL_NO_DEBUG_OVERLAY_3_PROJECTED_POINTS
          std::cout << "Warning: Projection into an edge\n";
#endif
          Halfedge_handle h = canonical_edge(cpp11::get<0>(inter_res));
          /// \todo this intersection point is going to be recomputed, see
          ///       if it's worth not redoing it

          Point_3 exact_source = to_exact(h->opposite()->vertex()->point());
          Vector_3 source_pt = cpp11::get<2>(inter_res) - exact_source;
          Vector_3 source_target =
            to_exact(h->vertex()->point()) - exact_source;
          //get alpha such that pt=source + alpha * (target - source)
          typename
          EK::FT alpha =
            source_target.x() != 0 ?
            source_pt.x() / source_target.x() :
            source_target.y() != 0 ?
            source_pt.y() / source_target.y():
            source_pt.z() / source_target.z();

          std::size_t vertex_index = modifier.add_point_on_edge(h, alpha);

          inter_infos[v_index_P2] = Inter_info(POLYHEDRON_EDGE, vertex_index, h);
#ifndef CGAL_NO_DEBUG_OVERLAY_3_PROJECTED_POINTS
          output_edge_proj_segments << "2 " << target1 << " " << target2 << "\n";
#endif
        }
        break;
        case POLYHEDRON_VERTEX: //this is a vertex
        {
#ifndef CGAL_NO_DEBUG_OVERLAY_3_PROJECTED_POINTS
          std::cout << "Warning: Projection on a vertex\n";
#endif

          Halfedge_handle h = cpp11::get<0>(inter_res);
          inter_infos[v_index_P2] = Inter_info(POLYHEDRON_VERTEX, get(m_vertex_ipmap, h->vertex()), h);
#ifndef CGAL_NO_DEBUG_OVERLAY_3_PROJECTED_POINTS
          output_vertex_proj_segments << "2 " << target1 << " " << target2 << "\n";
#endif
        }
        break;
        case POLYHEDRON_NONE:
          if (!allow_partial_overlay)
          {
            inter_infos[v_index_P2] = Inter_info(POLYHEDRON_NONE, -1, Halfedge_handle());
            std::cout << "WARNING: the point " << pit->point() << " has no projection,";
            std::cout << "try setting a larger value for the distance between the meshes.\n";
          }
          else
          {
            std::size_t index=modifier.add_vertex_with_no_projection(pit);
            exact_points.push_back(to_exact(pit->point()));
            inter_infos[v_index_P2] =
              Inter_info(POLYHEDRON_NONE, index, Halfedge_handle());
          }

#ifndef CGAL_NO_DEBUG_OVERLAY_3_PROJECTED_POINTS
          output_no_proj_segments << "2 " << target1 << " " << target2 << "\n";
#endif
      }
    }

#ifndef CGAL_NO_DEBUG_OVERLAY_3_PROJECTED_POINTS
    output_no_proj_segments.close();
    output_face_proj_segments.close();
    output_edge_proj_segments.close();
    output_vertex_proj_segments.close();
    output_proj_segments.close();
#endif

    /// \todo add a check to ensure 2 points are not projected on the same points
    ///       (I'm not talking about rounding issue necessarily)

#ifndef CGAL_NO_DEBUG_OVERLAY_3_PROJECTED_POINTS
    {
      std::ofstream output("debug/pts.xyz");
      for (std::size_t k=0; k<inter_infos.size(); ++k)
        if (cpp11::get<0>(inter_infos[k]) != POLYHEDRON_NONE)
          output << modifier.get_debug_point(cpp11::get<1>(inter_infos[k])) << "\n";
    }
#endif

#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
    std::ofstream debug_full_edge_output("debug/full_edges_output.cgal");
    std::ofstream debug_missing_output("debug/missing_output.cgal");
#endif

    std::vector<Halfedge_handle> P1_border;
    if (allow_partial_overlay && !P1.is_closed())
    {
      BOOST_FOREACH(Halfedge_handle h, halfedges(P1))
      if (h->is_border())
        P1_border.push_back(canonical_edge(h));
    }

    const std::size_t max_id = (std::numeric_limits<std::size_t>::max)();

    /// compute the intersection of edges
    EdgeSkipper skip;
    for (typename Polyhedron::Edge_iterator edge_it=P2.edges_begin(),
         edge_end=P2.edges_end();
         edge_it!=edge_end; ++edge_it)
    {
      /// skip edges filtered
      if (skip(edge_it)) continue;

      Halfedge_handle hedge = edge_it;

      /// We walk from v_src to v_tgt and intersect all edges encountered
      Vertex_handle v_src = hedge->opposite()->vertex();
      Vertex_handle v_tgt = hedge->vertex();

      Polyhedron_simplex_type inter_type_source = cpp11::get<0>(inter_infos[get(m_vertex_ipmap, v_src)]);
      Polyhedron_simplex_type inter_type_target = cpp11::get<0>(inter_infos[get(m_vertex_ipmap, v_tgt)]);

      if (inter_type_source == POLYHEDRON_NONE)
      {
        if (inter_type_target != POLYHEDRON_NONE)
        {
          std::swap(inter_type_source, inter_type_target);
          std::swap(v_src, v_tgt);
          hedge=hedge->opposite();
        }
        else
        {
          // compute all the intersections with P1 border edges
          std::vector< std::pair<Halfedge_handle, Inter_res> > intersections;
          compute_intersections(hedge, P1_border, min_half_diameter, normals,
                                std::back_inserter(intersections));

          if (intersections.empty())
            modifier.add_edges_with_no_projection(hedge, skip);
          else
          {
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
            std::cout << "NONE-NONE - Handling edge " << hedge->opposite()->vertex()->point()
                      << " " << hedge->vertex()->point() << "\n";
            std::cout << "{\n";
            std::cout << "intersections.size() " << intersections.size() << "\n";
#endif

            // here we need to walk from an intersection point to another
            // we extract the left most and right most intersections and
            // trigger the walk using those.
            std::sort(intersections.begin(),
                      intersections.end(),
                      Less_edge_intersection_alphas());
            partial_edge_overlay(hedge, intersections, max_id, max_id,
                                 P1, normals, modifier);
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
            std::cout << "}\n";
#endif
          }
          continue;
        }
      }

      Halfedge_handle h_src = cpp11::get<2>(inter_infos[get(m_vertex_ipmap, v_src)]);
      Halfedge_handle h_tgt = cpp11::get<2>(inter_infos[get(m_vertex_ipmap, v_tgt)]);

      std::size_t source_index = cpp11::get<1>(inter_infos[get(m_vertex_ipmap, v_src)]);
      std::size_t target_index = cpp11::get<1>(inter_infos[get(m_vertex_ipmap, v_tgt)]);

      // now walk in P1 along hedge
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
      {
        Vertex_handle v_src = hedge->opposite()->vertex();
        Vertex_handle v_tgt = hedge->vertex();
        std::ofstream debug_full_edge_output("debug/normal_interpolation.cgal");
        Input_point source = v_src->point();
        Input_point target = v_tgt->point();
        std::cout << "==> Walking " << source << " " << target <<"\n";
        const int nb_steps=100;
        typename Input_kernel::Vector_3 source_to_target = (target-source)/(nb_steps-1);
        typename Input_kernel::Vector_3 n_source = to_input(normals[get(m_vertex_ipmap, v_src)]);
        typename Input_kernel::Vector_3 n_target = to_input(normals[get(m_vertex_ipmap, v_tgt)]);
        typename Input_kernel::Vector_3 n_intrpl = (n_target-n_source)/(nb_steps-1);

        for (int i=0; i<nb_steps; ++i)
        {
          Input_point intrpl_pt = source + i * source_to_target;
          typename Input_kernel::Vector_3 intrpl_vect = n_source + i * n_intrpl;
          debug_full_edge_output << "2 " << intrpl_pt + (- intrpl_vect)
                                 << "  " << intrpl_pt + (intrpl_vect) << "\n";
        }
      }
#endif

      Intersection_predicate predicate(normals, hedge, m_vertex_ipmap);
      Walker walker;
      Walk_visitor walk_visitor_1(source_index, target_index,
                                  h_src, h_tgt,inter_type_target,
                                  hedge, modifier, m_vertex_ipmap);

      bool target_reached = walker(P1, h_src, inter_type_source,
                                   h_tgt, inter_type_target,
                                   predicate, walk_visitor_1);

      // indices of endpoints of a edge with partial projection
      std::size_t vid1=walk_visitor_1.last_seen_vertex_index(),
                  vid2=(std::numeric_limits<std::size_t>::max)();
      Barycentric_NT beta_1=0;

      CGAL_assertion(inter_type_source!=POLYHEDRON_NONE || target_reached);

      if (!target_reached)
      {
        CGAL_assertion(inter_type_target!=POLYHEDRON_NONE);
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
        std::cout << "Target not reached!\n";
#endif
        std::swap(inter_type_source, inter_type_target);
        std::swap(h_src, h_tgt);
        std::swap(source_index, target_index);
        hedge=hedge->opposite();

        Intersection_predicate predicate(normals, hedge, m_vertex_ipmap);
        Walker walker;
        Walk_visitor walk_visitor_2(source_index, target_index,
                                    h_src, h_tgt, inter_type_target,
                                    hedge, modifier, m_vertex_ipmap);

        bool target_reached = walker(P1, h_src, inter_type_source,
                                     h_tgt, inter_type_target,
                                     predicate, walk_visitor_2);

        CGAL_assertion(!target_reached);
        // id must be swapt too!
        vid2=vid1;
        vid1=walk_visitor_2.last_seen_vertex_index();
        beta_1=walk_visitor_2.last_beta();
      }

      CGAL_assertion(inter_type_source!=POLYHEDRON_NONE);

      if (inter_type_target==POLYHEDRON_NONE || !target_reached)
      {
        std::vector< std::pair<Halfedge_handle, Inter_res> > intersections;
        compute_intersections(hedge, P1_border, min_half_diameter, normals,
                              std::back_inserter(intersections));

        bool remove_first=false, remove_last=false;
        std::size_t final_id=max_id;
        if (inter_type_target==POLYHEDRON_NONE)
        {
          if (intersections.size()>1)
          {
            remove_first=true;
            beta_1=walk_visitor_1.last_beta();
          }
        }
        else
        {
          if (intersections.size()>2)
          {
            remove_first=true;
            remove_last=true;
          }
          final_id=vid2;
        }

        if (remove_first)
        {
          std::sort(intersections.begin(),
                    intersections.end(),
                    Less_edge_intersection_alphas());
          //remove already discovered intersections
          intersections.erase(intersections.begin());
          if (remove_last) intersections.pop_back();

          // filter more intersections in case border vertices were encountered
          while (!intersections.empty() &&
                 beta_1 >= intersections[0].second->first)
            intersections.erase(intersections.begin());
          if (!target_reached)
          {
            const Barycentric_NT &beta_2=1-walk_visitor_1.last_beta(); // hedge swapt
            while (!intersections.empty() &&
                   beta_2 <= intersections.back().second->first)
              intersections.pop_back();
          }
          if (!intersections.empty())
          {
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
            std::cout << "Handling extra walks...\n{\n";
#endif
            partial_edge_overlay(hedge, intersections, vid1, final_id,
                                 P1, normals, modifier);
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
            std::cout << "}\n";
#endif
            continue;
          }
        }
        modifier.add_edges_with_partial_projection(
          hedge, vid1, final_id, EdgeSkipper());
      }
    }

    P1.delegate(modifier);

    /// init input_facet_ids
    std::size_t final_nb_facets=modifier.facets.size();
    input_facet_ids.reserve(final_nb_facets);
    for (std::size_t i=0; i<final_nb_facets; ++i)
    {
      CGAL_assertion(modifier.facets[i]->id()==input_facet_ids.size());
      input_facet_ids.push_back(
        std::make_pair(modifier.facets[i]->input_id(), max_id)
      );
    }

#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
    {
      CGAL_assertion(P1.is_valid());
      std::ofstream output("debug/out.off");
      output << P1;
    }
#endif

    // Only used when allow_partial_overlay is true. Contains
    // halfedges which are on the boundary of the region which is covered
    // by the projection of P2
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
    std::ofstream output_not_to_cross("debug/not_to_cross.cgal");
#endif
    /// \todo add a property into edges
    /// \todo set this for all edges that are projection of P2, that way we can remove interior_halfedges
    std::set<Halfedge_handle> hedge_not_to_cross;

    /// now transfert the facet id from P2 facets to P1 facets
    for (typename Polyhedron::Halfedge_iterator
         hedge_it=P2.halfedges_begin(),hedge_end=P2.halfedges_end();
         hedge_it!=hedge_end; ++hedge_it)
    {
      std::size_t hindex=hedge_it->id();
      const std::vector<Halfedge_handle> &hedges=modifier.halfedge_map()[hindex];

      if (hedge_it->is_border())
      {
        if (allow_partial_overlay)
          for (std::size_t k=0, end=hedges.size(); k<end; ++k)
          {
            hedge_not_to_cross.insert(hedges[k]->opposite());
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
            output_not_to_cross << "2 " << hedges[k]->opposite()->vertex()->point()
                                << "  " << hedges[k]->vertex()->point() << "\n";
#endif
          }
        continue;
      }

      for (std::size_t k=0, end=hedges.size(); k<end; ++k)
      {
        if (hedges[k]->is_border()) continue;

        CGAL_assertion(
          input_facet_ids[ hedges[k]->facet()->id() ].second==max_id ||
          input_facet_ids[ hedges[k]->facet()->id() ].second==
          hedge_it->facet()->input_id());

        input_facet_ids[ hedges[k]->facet()->id() ].second =
          hedge_it->facet()->input_id();
      }
    }
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
    output_not_to_cross.close();
#endif

    // If a facet of P1 is totally included into a facet of P2, the P2 index is not set.
    // It is safe to consider all neighboring facets since edges are
    // not of the boundary of the projection of a facet of P2
    //   We start by collecting facet with a non-initialized id which have a neighbor with an initialized id
    std::vector<Halfedge_handle> stack;
    for (std::size_t i=0; i< final_nb_facets; ++i)
    {
      if (input_facet_ids[i].second  != max_id)
      {
        Halfedge_handle hedge = modifier.facets[i]->halfedge();
        for (int i=0; i<3; ++i)
        {
          if (!hedge->opposite()->is_border() &&
              !hedge_not_to_cross.count(hedge) &&
              input_facet_ids[hedge->opposite()->facet()->id()].second==max_id)
            stack.push_back(hedge->opposite());
          hedge=hedge->next();
        }
      }
    }

    // we now go through the stack and update the id of facets, putting
    // in the stack incident facets with non-initialized ids
    while (!stack.empty())
    {
      Halfedge_handle hedge=stack.back();
      stack.pop_back();
      CGAL_assertion(!hedge->is_border_edge());
      if (input_facet_ids[hedge->facet()->id()].second==max_id)
      {
        CGAL_assertion(input_facet_ids[hedge->opposite()->facet()->id()].second!=max_id);

        input_facet_ids[hedge->facet()->id()].second =
          input_facet_ids[hedge->opposite()->facet()->id()].second;

        for (int i=0; i<2; ++i)
        {
          hedge=hedge->next();
          if (!hedge->opposite()->is_border() &&
              !hedge_not_to_cross.count(hedge) &&
              input_facet_ids[hedge->opposite()->facet()->id()].second==max_id)
            stack.push_back(hedge->opposite());
        }
      }
      else
        CGAL_assertion(input_facet_ids[hedge->facet()->id()].second ==
                       input_facet_ids[hedge->opposite()->facet()->id()].second);
    }
    // from here facets with max_id for P2_id are not covered by the projection of P2
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
    if (!allow_partial_overlay)
      for (std::size_t i=0; i< final_nb_facets; ++i)
      {
        if (input_facet_ids[i].second == max_id)
        {
          std::cout << "Looks like we have a bug here!!";
          std::cerr << "4 "
                    << modifier.facets[i]->halfedge()->vertex()->point() << " "
                    << modifier.facets[i]->halfedge()->next()->vertex()->point() << " "
                    << modifier.facets[i]->halfedge()->next()->next()->vertex()->point() << " "
                    << modifier.facets[i]->halfedge()->vertex()->point() << "\n";
        }
      }
    else
    {
      // debug: print face with no cell id from P1
      std::vector<Facet_handle> faces;
      for (std::size_t i=0; i< final_nb_facets; ++i)
        if (modifier.facets[i]->input_id() == max_id)
          faces.push_back(modifier.facets[i]);
      std::ofstream out("debug/no_cellid1.off");
      out << "OFF\n" << faces.size()*3 << " " << faces.size() << " 0\n";
      BOOST_FOREACH(Facet_handle fh, faces)
      {
        out << fh->halfedge()->vertex()->point() << "\n";
        out << fh->halfedge()->next()->vertex()->point() << "\n";
        out << fh->halfedge()->prev()->vertex()->point() << "\n";
      }
      for (std::size_t i=0; i<faces.size(); ++i)
        out << "3 " << 3*i << " " << 3*i+1 << " " << 3*i+2 << "\n";
      out.close();
      // debug: print face with no cell id from P2
      faces.clear();
      for (std::size_t i=0; i< final_nb_facets; ++i)
        if (input_facet_ids[i].second == max_id)
          faces.push_back(modifier.facets[i]);
      out.open("debug/no_cellid2.off");
      out << "OFF\n" << faces.size()*3 << " " << faces.size() << " 0\n";
      BOOST_FOREACH(Facet_handle fh, faces)
      {
        out << fh->halfedge()->vertex()->point() << "\n";
        out << fh->halfedge()->next()->vertex()->point() << "\n";
        out << fh->halfedge()->prev()->vertex()->point() << "\n";
      }
      for (std::size_t i=0; i<faces.size(); ++i)
        out << "3 " << 3*i << " " << 3*i+1 << " " << 3*i+2 << "\n";
    }
#endif
  }

  void run(Polyhedron &P1,
           Polyhedron &P2,
           double min_half_diameter,
           const std::list<
           std::pair< std::vector<Halfedge_handle>,
           std::vector<Halfedge_handle> >
           >& matching_features,
           bool allow_partial_overlay,
           std::vector< std::pair<std::size_t, std::size_t> > &input_facet_ids
          )
  {
    typedef CGAL::AABB_halfedge_graph_segment_primitive<Polyhedron> Primitive;
    typedef CGAL::AABB_traits<Input_kernel, Primitive> Traits;
    typedef CGAL::AABB_tree<Traits> Tree;
    typedef typename boost::graph_traits<Polyhedron>::edge_descriptor edge_descriptor;

    std::vector< Vertex_handle> vertices_of_P2_already_matching;
    std::vector< std::pair<Halfedge_handle, double> > simplices_of_P1_already_matching;

    /// Compute the projection of each feature of P2 onto the corresponding feature of P1
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
    std::ofstream output("debug/feature_matching.cgal");
#endif
    typedef
    std::pair< std::vector<Halfedge_handle>,
        std::vector<Halfedge_handle> > Vector_pair;
    BOOST_FOREACH(const Vector_pair& vector_pair, matching_features)
    {
      Tree tree;
      BOOST_FOREACH(Halfedge_handle h, vector_pair.first)
      {
        edge_descriptor ed(h);
        tree.insert(Primitive(ed, P1));
      }
      tree.build();

      tree.accelerate_distance_queries();

      BOOST_FOREACH(Halfedge_handle h, vector_pair.second)
      {
        Vertex_handle vh = h->vertex();
        Input_point interpt;
        edge_descriptor e;
        boost::tie(interpt,e) = tree.closest_point_and_primitive(vh->point());
        h = canonical_edge(halfedge(e,P1));
        const Input_point &source = h->opposite()->vertex()->point();
        const Input_point &target = h->vertex()->point();

        Input_vector source_pt = interpt - source;
        Input_vector source_target = target - source;
        //get alpha such that pt=source + alpha * (target - source)
        double alpha =
          source_target.x() != 0 ?
          source_pt.x() / source_target.x() :
          source_target.y() != 0 ?
          source_pt.y() / source_target.y():
          source_pt.z() / source_target.z();
        if (alpha < 0) alpha=0;
        else if (alpha>1) alpha=1;

#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
        output << "2 " << vh->point() <<  " " << source+alpha * (target-source) << "\n";
#endif
        vertices_of_P2_already_matching.push_back(vh);
        simplices_of_P1_already_matching.push_back(std::make_pair(h,alpha));
      }

#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
      /// TAG_SL_GOOD_PROJECTION
      /// test to ensure the projection was done correctly
      std::map<Vertex_handle, std::size_t> vid;
      std::size_t index=0;
      BOOST_FOREACH(Vertex_handle vh, vertices_of_P2_already_matching)
      vid[vh]=index++;
      std::map<Vertex_handle, std::vector<std::size_t> > incident_vertex_ids;
      BOOST_FOREACH(Halfedge_handle h, vector_pair.second)
      {
        incident_vertex_ids[h->vertex()].push_back(vid[h->opposite()->vertex()]);
        incident_vertex_ids[h->opposite()->vertex()].push_back(vid[h->vertex()]);
      }
      typedef std::pair<Vertex_handle, std::vector<std::size_t> > Pair;
      BOOST_FOREACH(const Pair& pair, incident_vertex_ids)
      {
        if (pair.second.size()!=2) continue;
        // points on P2
        Point_3 p0 = to_exact(vertices_of_P2_already_matching[ pair.second[0] ]->point());
        Point_3 p1 = to_exact(pair.first->point());
        Point_3 p2 = to_exact(vertices_of_P2_already_matching[ pair.second[1] ]->point());
        // points projected on P1
        std::pair<Halfedge_handle, double> h_and_alpha = simplices_of_P1_already_matching[ pair.second[0] ];
        Point_3 pp0 = compute_point_on_edge(h_and_alpha.first, h_and_alpha.second);
        h_and_alpha = simplices_of_P1_already_matching[ vid[pair.first] ];
        Point_3 pp1 = compute_point_on_edge(h_and_alpha.first, h_and_alpha.second);
        h_and_alpha = simplices_of_P1_already_matching[ pair.second[1] ];
        Point_3 pp2 = compute_point_on_edge(h_and_alpha.first, h_and_alpha.second);

        CGAL_assertion(CGAL::sign((p1- p0)*(p2- p1)) ==
                       CGAL::sign((pp1-pp0)*(pp2-pp1)));
      }
#endif

    }
#ifndef CGAL_NO_DEBUG_OVERLAY_3_REFINEMENT
    output.close();
#endif

    run(P1, P2, min_half_diameter,
        vertices_of_P2_already_matching,
        simplices_of_P1_already_matching,
        allow_partial_overlay,
        input_facet_ids);
  }

  /// \todo make P2 const
  /// \todo turns input_facet_ids into a property_map
  void run(Polyhedron &P1,
           Polyhedron &P2,
           double min_half_diameter,
           bool allow_partial_overlay,
           std::vector< std::pair<std::size_t, std::size_t> > &input_facet_ids
          )
  {
    std::vector< Vertex_handle> vertices_of_P2_already_matching;
    std::vector< std::pair<Halfedge_handle, double> > simplices_of_P1_already_matching;

    run(P1, P2, min_half_diameter,
        vertices_of_P2_already_matching,
        simplices_of_P1_already_matching,
        allow_partial_overlay,
        input_facet_ids);
  }
};

} //end of namespace Surface_mesh_overlay_3

} //end of namespace CGAL

#endif // CGAL_SURFACE_MESH_OVERLAY_3_SURFACE_MESH_OVERLAY_3_H
