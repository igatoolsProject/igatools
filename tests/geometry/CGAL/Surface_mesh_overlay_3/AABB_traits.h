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


#ifndef CGAL_SURFACE_MESH_OVERLAY_AABB_TRAITS_H
#define CGAL_SURFACE_MESH_OVERLAY_AABB_TRAITS_H

#include <CGAL/AABB_traits.h>
#include <CGAL/Bbox_3.h>

namespace CGAL
{
namespace Surface_mesh_overlay_3
{

template <class GeomTraits>
struct Query_type
{
  typename GeomTraits::Segment_3 s1, s2;
};

template <class  GeomTraits, class AABBPrimitive>
struct AABB_traits :
  public CGAL::AABB_traits<GeomTraits, AABBPrimitive>
{
  // if Do_intersect returned true `last_segment_intersected`
  // is set to 1 if s1 triggered the intersection,
  // 2 if s2 did it, and 0 otherwise.
  mutable int last_segment_intersected;

  typedef CGAL::AABB_traits<GeomTraits, AABBPrimitive> Base;
  typedef AABB_traits<GeomTraits, AABBPrimitive> Self;
  typedef typename Base::Do_intersect Do_intersect_base;

  AABB_traits() : Base(), last_segment_intersected(0)
  {}

  struct Do_intersect : public Do_intersect_base
  {
    const Self &m_traits; /// \todo remove me and make the one of the base
    ///       class protected
    Do_intersect(const Self &traits)
      : Do_intersect_base(traits)
      , m_traits(traits)
    {}

    using Do_intersect_base::operator();

    bool
    operator()(const Query_type<GeomTraits> &query, const CGAL::Bbox_3 &bbox)
    {
      return Do_intersect_base::operator()(query.s1, bbox) ||
             Do_intersect_base::operator()(query.s2, bbox);
    }

    bool
    operator()(const Query_type<GeomTraits> &query, const AABBPrimitive &primitive)
    {
      bool do_inter = Do_intersect_base::operator()(query.s1, primitive);
      if (do_inter)
      {
        m_traits.last_segment_intersected=1;
        return true;
      }
      do_inter =  Do_intersect_base::operator()(query.s2, primitive);
      if (do_inter)
      {
        m_traits.last_segment_intersected=2;
        return true;
      }
      m_traits.last_segment_intersected=0;
      return false;
    }
  };

  Do_intersect do_intersect_object() const
  {
    return Do_intersect(*this);
  }

};

}
};


#endif // CGAL_SURFACE_MESH_OVERLAY_AABB_TRAITS_H
