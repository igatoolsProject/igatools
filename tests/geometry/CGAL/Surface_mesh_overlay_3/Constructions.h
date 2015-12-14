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

#ifndef CGAL_SURFACE_MESH_OVERLAY_3_CONSTRUCTIONS_H
#define CGAL_SURFACE_MESH_OVERLAY_3_CONSTRUCTIONS_H

#include <CGAL/Kernel_traits.h>
#include <CGAL/Cartesian_converter.h>

#ifdef CGAL_DEBUG_SURFACE_MESH_OVERLAY_CONSTRUCTIONS
#include <iostream>
#define CGAL_OVERLAY_CSTRCT_TRACE(X) std::cout << X;
#else
#define CGAL_OVERLAY_CSTRCT_TRACE(X)
#endif

namespace CGAL
{

namespace Surface_mesh_overlay_3
{

/// compute the intersection point of the face of h with s
template <class Kernel, class Polyhedron>
typename Kernel::Point_3
compute_inter_pt(
  typename Polyhedron::Halfedge_handle h,
  const typename Kernel::Segment_3 &s)
{
  typedef typename Kernel_traits<
  typename Polyhedron::Point_3>::Kernel Input_kernel;
  typedef typename Kernel::Point_3 Point_3;

  CGAL::Cartesian_converter<Input_kernel, Kernel> to_exact;

  typename Kernel::Plane_3 p(
    to_exact(h->vertex()->point()),
    to_exact(h->next()->vertex()->point()),
    to_exact(h->opposite()->vertex()->point()));

  boost::optional< boost::variant<Point_3, typename Kernel::Segment_3> >
  inter = CGAL::intersection(p,s);
  CGAL_assertion(inter!=boost::none);

  const Point_3 *pt = boost::get<Point_3>(&(*inter));

  CGAL_assertion(pt !=NULL);
  return  *pt;
}

/// compute the intersection point of the face of h with l
template <class Kernel, class Polyhedron>
typename Kernel::Point_3
compute_inter_pt(
  typename Polyhedron::Halfedge_handle h,
  const typename Kernel::Line_3 &l)
{
  typedef typename Kernel_traits<
  typename Polyhedron::Point_3>::Kernel Input_kernel;
  typedef typename Kernel::Point_3 Point_3;

  CGAL::Cartesian_converter<Input_kernel, Kernel> to_exact;

  typename Kernel::Plane_3 p(
    to_exact(h->vertex()->point()),
    to_exact(h->next()->vertex()->point()),
    to_exact(h->opposite()->vertex()->point()));

  boost::optional< boost::variant<Point_3, typename Kernel::Line_3> >
  inter = CGAL::intersection(p,l);
  CGAL_assertion(inter);

  const Point_3 *pt = boost::get<Point_3>(&(*inter));

  CGAL_assertion(pt !=NULL);
  return *pt;
}

/// Compute the intersection of the edge b0b1 from P1 and g0g1 from P2
/// The point lies in b0b1.
/// The optional is empty if no intersection exists
/// A bool is returned if it is either b0 or b1 (false or true respectively)
/// otherwise, alpha \in ]0,1[ is returned, with
/// interpt = b0 + alpha * (b1-b0)
/// For the parameter names, we use the same notations as in Section 3.2 of
///   Overlaying Surface Meshes, Part I;
///   International Journal on Computational Geometry and Applications.
///   Vol 14(6), pages 379-402. December 2004
template <class Kernel, class Polyhedron>
boost::optional<
std::pair<
typename Root_of_traits<typename Kernel::FT>::Root_of_2,
         boost::variant<bool, typename Root_of_traits<typename Kernel::FT>::Root_of_2>
         >
         >
         compute_edge_intersection(
           const typename Kernel::Point_3 &b0,
           const typename Kernel::Point_3 &b1,
           const typename Kernel::Point_3 &g0,
           const typename Kernel::Point_3 &g1,
           const typename Kernel::Vector_3 &d0,
           const typename Kernel::Vector_3 &d1)
{
  typedef typename Kernel::FT                                          FT;
  typedef typename Root_of_traits<FT>::Root_of_2                Root_of_2;
  typedef boost::variant<bool,Root_of_2>                          Variant;
  typedef typename Kernel::Point_3                                Point_3;
  typedef typename Kernel::Vector_3                              Vector_3;

  //std::cout << "2 " << b0 << " " << b1 <<"\n";
  //std::cout << "2 " << g0 << " " << g1 <<"\n";
  //std::cout << "2 " << d1-d0 <<"\n";

  //precompute some vectors
  Vector_3 b0b1(b0,b1);
  Vector_3 g0g1(g0,g1);
  Vector_3 b0g0(b0,g0);
  Vector_3 b0b1_c_d0 = cross_product(b0b1,d0);

  // linear case when d1 and d0 are collinear
  if (CGAL::collinear(Point_3(0,0,0),
                      ORIGIN+d1,
                      ORIGIN+d0))
  {
    CGAL_OVERLAY_CSTRCT_TRACE("                      collinear case\n")
    // We assume d1==d0
    // FT c2 = 0
    FT c1 = b0b1_c_d0 * g0g1;
    FT c0 = b0b1_c_d0 * b0g0;

    if (c1 == 0) return boost::none;

    FT beta = -c0/c1;

    if (beta < 0 || beta >1) return boost::none;
    Vector_3 l = CGAL::cross_product(g0g1 , d0);
    FT l_b0b1 = l * b0b1;
    if (l_b0b1 == 0) return boost::none;

    FT alpha = b0g0 * l / l_b0b1;

    if (alpha < 0 || alpha >1) return boost::none;

    if (alpha == 0)  return std::make_pair(Root_of_2(beta),Variant(false)); //b0
    if (alpha == 1)  return std::make_pair(Root_of_2(beta),Variant(true));  //b1

    return std::make_pair(Root_of_2(beta), Variant(Root_of_2(alpha)));
  }

  //precompute some other vectors
  Vector_3 d1_m_d0 = d1 - d0;
  Vector_3 b0b1_c_d1_m_d0 = CGAL::cross_product(b0b1,d1_m_d0);

  FT c2 = b0b1_c_d1_m_d0 * g0g1;
  FT c1 = b0b1_c_d0 * g0g1 + b0b1_c_d1_m_d0 * b0g0;
  FT c0 = b0b1_c_d0 * b0g0;

  // d1_m_d0 is the projection vector for the intersection
  if (c2 == 0)
  {
    CGAL_OVERLAY_CSTRCT_TRACE("                      c2=0\n")
    if (c1 == 0) return boost::none;
    FT beta = -c0/c1;

    if (beta < 0 || beta >1) return boost::none;
    Vector_3 l = CGAL::cross_product(g0g1 , d0 + beta * d1_m_d0);
    FT l_b0b1 = l * b0b1;

    if (l_b0b1 == 0) return boost::none;

    FT alpha = b0g0 * l / l_b0b1;

    if (alpha < 0 || alpha >1) return boost::none;

    if (alpha == 0)  return std::make_pair(Root_of_2(beta),Variant(false)); //b0
    if (alpha == 1)  return std::make_pair(Root_of_2(beta),Variant(true));  //b1

    return std::make_pair(Root_of_2(beta),Variant(Root_of_2(alpha)));
  }

  FT delta = CGAL::square(c1) - 4 * c2 * c0;

  switch (CGAL::sign(delta))
  {
    case CGAL::NEGATIVE:
      CGAL_OVERLAY_CSTRCT_TRACE("                      delta<0\n")
      return boost::none;
    case CGAL::ZERO:
      CGAL_OVERLAY_CSTRCT_TRACE("                      delta=0\n")
      return boost::none;// if the curve is tangent to the edge we do not go out of the face

    case CGAL::POSITIVE:
    {
      CGAL_OVERLAY_CSTRCT_TRACE("                      delta>0\n")
      FT inv_2c2 = FT(1)/FT(2)/c2;
      CGAL::cpp11::array<Root_of_2,2> betas= {{
          CGAL::make_root_of_2(-c1*inv_2c2,inv_2c2,delta),
          CGAL::make_root_of_2(-c1*inv_2c2,-inv_2c2,delta)
        }
      };
      std::vector<boost::optional< std::pair<Root_of_2,Variant> > > possible_res;

      for (int i=0; i<2; ++i)
      {
        CGAL_OVERLAY_CSTRCT_TRACE("                      beta-" << i <<  " = " << to_double(betas[i]) << "\n")
        const Root_of_2 &beta=betas[i];

        if (beta < 0 || beta > 1) continue;

        Root_of_2 d0_p_beta_x_d1_m_d0__x = d0.x() + beta * d1_m_d0.x();
        Root_of_2 d0_p_beta_x_d1_m_d0__y = d0.y() + beta * d1_m_d0.y();
        Root_of_2 d0_p_beta_x_d1_m_d0__z = d0.z() + beta * d1_m_d0.z();

        Root_of_2 lx = g0g1.y()*d0_p_beta_x_d1_m_d0__z -
                       g0g1.z()*d0_p_beta_x_d1_m_d0__y;
        Root_of_2 ly = g0g1.z()*d0_p_beta_x_d1_m_d0__x -
                       g0g1.x()*d0_p_beta_x_d1_m_d0__z;
        Root_of_2 lz = g0g1.x()*d0_p_beta_x_d1_m_d0__y -
                       g0g1.y()*d0_p_beta_x_d1_m_d0__x;

        Root_of_2 l_b0b1 = lx * b0b1.x() + ly * b0b1.y() + lz * b0b1.z();
        Root_of_2 alpha;

        if (l_b0b1 == 0)
        {
          CGAL_OVERLAY_CSTRCT_TRACE("                      l_b0b1=0\n")
          // In that case g0g1^di and b0b1 are orthogonal
          // We cannot compute the intersection as the intersection of the
          // plane P containing g0g1 with l=g0g1^di as normal since b0b1 is parallel
          // to that plane.
          // there is an intersection only if b0 lies inside P
          // this can happen only if b0,b1,g0 and g1 are coplanar since
          // g0 and g1 \in P and b0b1^l=0.
          if (CGAL::orientation(b0, b1, g0, g1)!=CGAL::COPLANAR)
            return boost::none;

          // we use the equation g + gamma*d = b
          // where g=g0+beta*(g1-g0), d= d0+beta(d1-d0) and b=b0+alpha*(b1-b0)
          // we try replacing gamma by an expression computed from one line
          // into another to get a value for alpha
          Root_of_2 gx = g0.x() + beta * g0g1.x();
          Root_of_2 gy = g0.y() + beta * g0g1.y();
          const Root_of_2 &dx = d0_p_beta_x_d1_m_d0__x;
          const Root_of_2 &dy = d0_p_beta_x_d1_m_d0__y;
          const Root_of_2 &dz = d0_p_beta_x_d1_m_d0__z;

          Root_of_2 denum = b0b1.x()*dy-b0b1.y()*dx;
          if (denum != 0)
            alpha = ((gx-b0.x())*dy+(b0.y()-gy)*dx) / denum;
          else
          {
            Root_of_2 gz = g0.z() + beta * (g1.z() - g0.z());
            denum = b0b1.x()*dz-b0b1.z()*dx;
            if (denum!=0)
              alpha = ((gx-b0.x())*dz+(b0.z()-gz)*dx) / denum;
            else
            {
              denum = b0b1.y()*dz-b0b1.z()*dy;
              // the intersection is b0b1 if di ^ b0b1==NULL_VECTOR and P contains b0.
              if (denum == 0) return boost::none; // \todo shall we really return no intersection?
              alpha = ((gy-b0.y())*dz+(b0.z()-gz)*dy) / denum;
            }
          }
        }
        else
        {
          CGAL_OVERLAY_CSTRCT_TRACE("                      l_b0b1!=0\n")
          Root_of_2 l_b0g0 = lx * b0g0.x() + ly * b0g0.y() + lz * b0g0.z();
          alpha = l_b0g0 / l_b0b1;
        }
        CGAL_OVERLAY_CSTRCT_TRACE("                      alpha" <<  " = " << to_double(alpha) << "\n")
        if (alpha < 0 || alpha >1) continue;
        if (alpha == 0)
        {
          possible_res.push_back(std::make_pair(beta, Variant(false)));  //b0
          continue;
        }
        if (alpha == 1)
        {
          possible_res.push_back(std::make_pair(beta, Variant(true)));  //b1
          continue;
        }

        possible_res.push_back(std::make_pair(beta, Variant(alpha)));
      }

      if (possible_res.size()==1)
        return possible_res[0];

      return boost::none;
    }
    default:
      CGAL_assertion(!"Should not get here\n");
      return boost::none;
  }
}

/// Compute the intersection of the edge hp1 from P1 and hp2 from P2
/// The point lies in hp1.
/// The optional is empty if no intersection exists
/// A bool is returned if it is an end point of hp1 (source=false, target=true)
/// otherwise, alpha \in ]0,1[ is returned, with
/// interpt = source(hp1) + alpha * (target(hp1)-source(hp1)
template <class Kernel, class Polyhedron, class VertexIndexPropertyMap>
boost::optional<
std::pair<
typename Root_of_traits<typename Kernel::FT>::Root_of_2,
         boost::variant<bool, typename Root_of_traits<typename Kernel::FT>::Root_of_2>
         >
         >
         compute_edge_intersection(
           typename Polyhedron::Halfedge_handle hp1,
           typename Polyhedron::Halfedge_handle hp2,
           const std::vector<typename Kernel::Vector_3> &normals,
           VertexIndexPropertyMap vipmap)
{
  typedef typename Kernel_traits<
  typename Polyhedron::Point_3>::Kernel                    Input_kernel;
  typedef typename Kernel::Point_3                                Point_3;
  typedef typename Kernel::Vector_3                              Vector_3;


  CGAL::Cartesian_converter<Input_kernel, Kernel> to_exact;

  Point_3 b0 = to_exact(hp1->opposite()->vertex()->point());
  Point_3 b1 = to_exact(hp1->vertex()->point());

  Point_3 g0 = to_exact(hp2->opposite()->vertex()->point());
  Point_3 g1 = to_exact(hp2->vertex()->point());

  Vector_3 d0 = normals[get(vipmap, hp2->opposite()->vertex())];
  Vector_3 d1 = normals[get(vipmap, hp2->vertex())];

  return compute_edge_intersection<Kernel, Polyhedron>(b0, b1, g0, g1, d0, d1);
}

template <class Polyhedron, class Kernel, class VertexIndexPropertyMap>
struct Edge_intersection
{
  const std::vector<typename Kernel::Vector_3> &m_normals;
  typename Polyhedron::Halfedge_handle m_hp2;
  VertexIndexPropertyMap m_vertex_ipmap;
public:

  typedef typename
  Root_of_traits<typename Kernel::FT>::Root_of_2 Barycentric_NT;

  typedef boost::optional< std::pair<Barycentric_NT, boost::variant<bool, Barycentric_NT > > >
  result_type;

  Edge_intersection(const std::vector<typename Kernel::Vector_3> &normals,
                    typename Polyhedron::Halfedge_handle hp2,
                    VertexIndexPropertyMap vipmap)
    :m_normals(normals), m_hp2(hp2), m_vertex_ipmap(vipmap) {}

  result_type
  operator()(typename Polyhedron::Halfedge_handle hp1) const
  {
    return
      compute_edge_intersection<Kernel, Polyhedron>(hp1, m_hp2,
                                                    m_normals, m_vertex_ipmap);
  }
};

} //end of namespace Surface_mesh_overlay_3

} //end of namespace CGAL

#endif //CGAL_SURFACE_MESH_OVERLAY_3_CONSTRUCTIONS_H
