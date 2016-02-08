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
#ifndef VTKUNSTRUCTURED_GRID_OVERLAY_3
#define VTKUNSTRUCTURED_GRID_OVERLAY_3


#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include "CGAL/Surface_mesh_overlay_3.h"

#ifdef OVERLAY_WITH_DEBUG_INFO
#include <fstream>
#include <CGAL/IO/Polyhedron_iostream.h>
#endif

#include "CGAL/vtkUnstructuredGrid_to_Polyhedron_3.h"
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkIntArray.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>

namespace internal
{

template<class Polyhedron>
vtkSmartPointer<vtkUnstructuredGrid>
convert_to_vtkUnstructuredGrid(
  Polyhedron &P,
  const std::vector< std::pair<std::size_t, std::size_t> > &input_facet_ids)
{
  vtkSmartPointer<vtkPoints> vtk_points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> vtk_cells = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkIntArray> vtk_cell_id1 = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> vtk_cell_id2 = vtkSmartPointer<vtkIntArray>::New();

  vtk_points->Allocate(P.size_of_vertices());
  vtk_cells->Allocate(input_facet_ids.size());

  vtk_cell_id1->SetNumberOfComponents(1);
  vtk_cell_id1->SetName("cellIds1");

  vtk_cell_id2->SetNumberOfComponents(1);
  vtk_cell_id2->SetName("cellIds2");

  std::size_t inum=0;
  for (typename Polyhedron::Vertex_iterator
       vit=P.vertices_begin(), vit_end=P.vertices_end(); vit!=vit_end; ++vit)
  {
    vtk_points->InsertNextPoint(vit->point().x(),
                                vit->point().y(),
                                vit->point().z());
    vit->id() = inum++;
  }

  const std::size_t max_id = (std::numeric_limits<std::size_t>::max)();
  for (typename Polyhedron::Facet_iterator
       fit=P.facets_begin(), fit_end=P.facets_end(); fit!=fit_end; ++fit)
  {
    vtkIdType cell[3];
    typename Polyhedron::Halfedge_handle h=fit->halfedge();
    for (int i = 0; i < 3; ++i, h=h->next())
      cell[i] = h->vertex()->id();
    vtk_cells->InsertNextCell(3, cell);
    if (input_facet_ids[fit->id()].first!=max_id)
      vtk_cell_id1->InsertNextValue(input_facet_ids[fit->id()].first);
    else
      vtk_cell_id1->InsertNextValue(-1);
    if (input_facet_ids[fit->id()].second != max_id)
      vtk_cell_id2->InsertNextValue(input_facet_ids[fit->id()].second);
    else
      vtk_cell_id2->InsertNextValue(-1);
  }

  vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =
    vtkSmartPointer<vtkUnstructuredGrid>::New();

  unstructuredGrid->SetPoints(vtk_points);
  unstructuredGrid->SetCells(VTK_TRIANGLE,vtk_cells);
  unstructuredGrid->GetCellData()->AddArray(vtk_cell_id1);
  unstructuredGrid->GetCellData()->AddArray(vtk_cell_id2);

  return unstructuredGrid;
}

} //end of internal namespace

template <class Polyhedron>
struct Skip_diagonals
{
  typedef typename Polyhedron::Edge_iterator Iterator;
  // skip boundary edges or diagonals
  bool operator()(Iterator it) const
  {
    return  !it->is_border_edge() && it->facet()->input_id() == it->opposite()->facet()->input_id();
  }
};

template <class Polyhedron>
struct Skip_diagonals_and_boundaries
{
  typedef typename Polyhedron::Edge_iterator Iterator;
  // skip boundary edges or diagonals
  bool operator()(Iterator it) const
  {
    return  it->is_border_edge() ||
            it->facet()->input_id() == it->opposite()->facet()->input_id();
  }
};

template <class Polyhedron>
void ensure_consistant_orientation(Polyhedron &P1, Polyhedron &P2,
                                   bool allow_partial_overlay, double min_half_diameter)
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel                 K;
  typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron>        Primitive;
  typedef CGAL::AABB_traits<K, Primitive>                                Traits;
  typedef CGAL::AABB_tree<Traits>                                          Tree;

  Tree tree(P2.facets_begin(),P2.facets_end(),P2);
  tree.accelerate_distance_queries();
  double min_sqd = std::numeric_limits<double>::max();
  bool do_inside_out = false;

  for (typename Polyhedron::Vertex_iterator  vit=P1.vertices_begin(),
       vit_end=P1.vertices_end();
       vit!=vit_end; ++vit)
  {
    std::pair<K::Point_3, typename Tree::Primitive_id> res =
      tree.closest_point_and_primitive(vit->point());

    if (allow_partial_overlay)
    {
      double sqd = CGAL::squared_distance(vit->point(),res.first);
      if (min_half_diameter > 0)
      {
        if (sqd > min_half_diameter*min_half_diameter) continue;
      }
      else
      {
        if (sqd > min_sqd) continue;
        min_sqd=sqd;
      }
    }
    namespace PMP = CGAL::Polygon_mesh_processing;
    K::Vector_3 v1 =  PMP::compute_vertex_normal(vit, P1);
    K::Vector_3 v2 =  PMP::compute_face_normal(res.second, P2);

    do_inside_out = (v1*v2 < 0);

    if (!allow_partial_overlay || min_half_diameter > 0) break;
  }
  if (do_inside_out)
    // reverse orientation of P2
    P2.inside_out();
}

/// returns true if after a call to
/// Polyhedron_3::split_facet(h,h->next()->next()) the triangles
/// are correctly oriented
template <class Kernel, class Halfedge_handle>
bool is_valid_edge_split(Halfedge_handle h)
{
  const typename Kernel::Point_3 &p0 = h->vertex()->point();
  const typename Kernel::Point_3 &p1 = h->next()->vertex()->point();
  const typename Kernel::Point_3 &p2 = h->next()->next()->vertex()->point();
  const typename Kernel::Point_3 &p3 = h->next()->next()->next()->vertex()->point();
  typename Kernel::Vector_3 v1=CGAL::cross_product(p2-p1,p0-p1);
  typename Kernel::Vector_3 v2=CGAL::cross_product(p0-p3,p2-p3);
  return v1 * v2 > 0;
}

vtkSmartPointer<vtkUnstructuredGrid>
vtkUnstructuredGridOverlay_3(
  vtkUnstructuredGrid *UG1,
  vtkUnstructuredGrid *UG2,
  double min_half_diameter,
  bool allow_partial_overlay,
  bool skip_orientation_consistancy_test)
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel                  K;
  typedef CGAL::Exact_predicates_exact_constructions_kernel                   EK;
  typedef CGAL::Surface_mesh_overlay_3::Poly_items_with_id_3          Poly_items;
  typedef CGAL::Polyhedron_3<K,Poly_items>                            Polyhedron;

  Polyhedron P1, P2;
  CGAL::convert_to_polyhedron(UG1, P1);
  CGAL::convert_to_polyhedron(UG2, P2);

  //We will refine P1 mesh
  // set index of P1 facets that need to be preserved when doing the refinement
  std::size_t findex=0;
  for (Polyhedron::Facet_iterator fit=P1.facets_begin(),
       fit_end=P1.facets_end();
       fit!=fit_end; ++fit)
  {
    fit->input_id()=findex++;
  }

  // set index of P2 facets that need to be transfered when doing the refinement
  findex=0;
  for (Polyhedron::Facet_iterator fit=P2.facets_begin(),
       fit_end=P2.facets_end();
       fit!=fit_end; ++fit)
  {
    fit->input_id()=findex++;
  }

  //triangulate the facets
  Polyhedron::Face_handle first_face = P1.facets_begin();
  if (first_face->halfedge()->next()->next()->next()!=first_face->halfedge())
  {
    std::vector< Polyhedron::Face_handle> facets;
    for (Polyhedron::Face_iterator fit=P1.facets_begin(),
         fit_end=P1.facets_end(); fit!=fit_end; ++fit)
    {
      facets.push_back(fit);
    }
    for (std::size_t i=0; i<facets.size(); ++i)
    {
#if 0 // commented by martinelli, Dec 15, 2015
      CGAL_assertion(is_valid_edge_split<K>(facets[i]->halfedge()));
#endif
      P1.split_facet(
        facets[i]->halfedge(),
        facets[i]->halfedge()->next()->next());
    }
  }

  first_face = P2.facets_begin();
  if (first_face->halfedge()->next()->next()->next()!=first_face->halfedge())
  {
    std::vector< Polyhedron::Face_handle> facets;
    for (Polyhedron::Face_iterator fit=P2.facets_begin(),
         fit_end=P2.facets_end(); fit!=fit_end; ++fit)
    {
      facets.push_back(fit);
    }
    for (std::size_t i=0; i<facets.size(); ++i)
    {
#if 0 // commented by martinelli, Dec 15, 2015
      CGAL_assertion(is_valid_edge_split<K>(facets[i]->halfedge()));
#endif
      P2.split_facet(
        facets[i]->halfedge(),
        facets[i]->halfedge()->next()->next());
    }
  }

  //make sure P1 and P2 are oriented consistantly
  if (!skip_orientation_consistancy_test)
    ensure_consistant_orientation(P1, P2,
                                  allow_partial_overlay, min_half_diameter);

#ifdef OVERLAY_WITH_DEBUG_INFO
  std::ofstream output("debug/P1.off");
  output.precision(20);
  output << P1;
  output.close();
  output.open("debug/P2.off");
  output.precision(20);
  output << P2;
  output.close();
#endif

  typedef Polyhedron::Halfedge_handle Halfedge_handle;
  std::list<
  std::pair< std::vector<Halfedge_handle>,
      std::vector<Halfedge_handle> >
      > matching_features;

  // no boundary matching if partial overlay are allowed
  if (!allow_partial_overlay && !P1.is_closed() && !P2.is_closed())
  {
    P1.normalize_border();
    P2.normalize_border();

    matching_features.push_back(
      std::pair< std::vector<Halfedge_handle>,
      std::vector<Halfedge_handle> > ());

    for (Polyhedron::Edge_iterator eit=P1.border_edges_begin(),
         eit_end=P1.edges_end();
         eit!=eit_end; ++eit)
    {
      matching_features.back().first.push_back(eit);
    }

    for (Polyhedron::Edge_iterator eit=P2.border_edges_begin(),
         eit_end=P2.edges_end();
         eit!=eit_end; ++eit)
    {
      matching_features.back().second.push_back(eit);
    }
  }

  // vector giving the id of the input facets corresponding to an output facet
  // input_facet_ids[ facet->id() ] gives the input-id of
  // the facet from P1 and P2 respectively
  std::vector< std::pair<std::size_t, std::size_t> > input_facet_ids;

  /// check if the boundaries are identical, in which case it is useless
  /// to do the matching with the AABB-tree
  bool do_matching_with_aabb_tree=false;
  if (!matching_features.empty())
  {
    if (matching_features.back().first.size() ==
        matching_features.back().second.size())
    {
      std::set<K::Point_3> points;
      BOOST_FOREACH(Halfedge_handle h, matching_features.back().first)
      points.insert(h->vertex()->point());
      BOOST_FOREACH(Halfedge_handle h, matching_features.back().second)
      if (points.find(h->vertex()->point())==points.end())
      {
        do_matching_with_aabb_tree=true;
        break;
      }
    }
    else
      do_matching_with_aabb_tree=true;
  }

  //compute the overlay
  if (do_matching_with_aabb_tree)
    CGAL::surface_mesh_overlay_with_edge_skipper_3
    <Skip_diagonals_and_boundaries<Polyhedron>, K, EK>(
      P1, P2, min_half_diameter, matching_features, input_facet_ids, allow_partial_overlay);
  else
    CGAL::surface_mesh_overlay_with_edge_skipper_3
    <Skip_diagonals<Polyhedron>, K, EK>(
      P1, P2, min_half_diameter, input_facet_ids, allow_partial_overlay);

  //convert P1 to a vtkUnstructuredGrid
  vtkSmartPointer<vtkUnstructuredGrid> result =
    vtkSmartPointer<vtkUnstructuredGrid>::New();

  return internal::convert_to_vtkUnstructuredGrid(P1, input_facet_ids);
}

#endif //VTKUNSTRUCTURED_GRID_OVERLAY_3
