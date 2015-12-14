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

#ifndef CGAL_VTKUNSTRUCTUREDGRID_TO_POLYHEDRON_3
#define CGAL_VTKUNSTRUCTUREDGRID_TO_POLYHEDRON_3

#include <CGAL/Modifier_base.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/array.h>

#include <vtkUnstructuredGrid.h>
#include <vtkCellType.h>
#include <vtkCell.h>
#include <vtkIdTypeArray.h>
#include <vtkSmartPointer.h>

namespace CGAL
{

namespace internal
{
template <int d>
struct vtk_nb_facet_to_word;

template<>
struct vtk_nb_facet_to_word<3>
{
  static const VTKCellType value = VTK_TRIANGLE;
};

template<>
struct vtk_nb_facet_to_word<4>
{
  static const VTKCellType value = VTK_QUAD;
};
}

template < class Polyhedron, int facet_size>
class vtkUnstructuredGrid_to_Polyhedron_3 :
  public Modifier_base<typename Polyhedron::Halfedge_data_structure>
{
  typedef typename Polyhedron::Vertex::Point Point;
  typedef typename Polyhedron::Halfedge_data_structure HDS;
  std::vector< cpp11::array<std::size_t,facet_size> > m_facet_indices;
  std::vector< Point > m_points;
public:

  vtkUnstructuredGrid_to_Polyhedron_3(vtkUnstructuredGrid *UG)
  {
    vtkSmartPointer<vtkIdTypeArray> id_array =
      vtkSmartPointer<vtkIdTypeArray>::New();

    UG->GetIdsOfCellsOfType(
      internal::vtk_nb_facet_to_word<facet_size>::value, id_array);

    std::map<vtkIdType, std::size_t> index_map;
    cpp11::array<std::size_t,facet_size> quad_indices;
    std::size_t index=0;

    vtkIdType nb_quads = id_array->GetNumberOfTuples();
    for (vtkIdType f=0; f<nb_quads; ++f)
    {
      vtkSmartPointer<vtkCell> cell = UG->GetCell(id_array->GetTuple(f)[0]);

      bool insertion_ok;
      std::map<vtkIdType, std::size_t>::iterator it;
      for (int i=0; i < facet_size; ++i)
      {
        CGAL::cpp11::tie(it,insertion_ok) =
          index_map.insert(std::make_pair(cell->GetPointId(i),index));
        if (insertion_ok)
        {
          ++index;
          double *coords = UG->GetPoint(it->first);
          m_points.push_back(Point(coords[0],coords[1],coords[2]));
        }
        quad_indices[i]=it->second;
      }
      m_facet_indices.push_back(quad_indices);
    }
  }

  void operator()(HDS &target)
  {
    std::size_t nbpts = m_points.size();
    std::size_t nbfacets = m_facet_indices.size();

    Polyhedron_incremental_builder_3<HDS> B(target);
    B.begin_surface(nbpts, nbfacets);
    for (std::size_t i=0; i<nbpts; ++i)
    {
      B.add_vertex(m_points[i]);
    }
    for (std::size_t i=0; i<nbfacets; ++i)
    {
      B.begin_facet();
      for (int k=0; k < facet_size; ++k)
      {
        B.add_vertex_to_facet(m_facet_indices[i][k]);
      }
      B.end_facet();
    }
    B.end_surface();
  }
};

template <class Polyhedron>
void convert_quads_to_polyhedron(
  vtkUnstructuredGrid *UG,
  Polyhedron &P)
{
  CGAL::vtkUnstructuredGrid_to_Polyhedron_3<Polyhedron, 4> modifier(UG);
  P.delegate(modifier);
}

template <class Polyhedron>
void convert_triangles_to_polyhedron(
  vtkUnstructuredGrid *UG,
  Polyhedron &P)
{
  CGAL::vtkUnstructuredGrid_to_Polyhedron_3<Polyhedron, 3> modifier(UG);
  P.delegate(modifier);
}

template <class Polyhedron>
void convert_to_polyhedron(
  vtkUnstructuredGrid *UG,
  Polyhedron &P)
{
  vtkSmartPointer<vtkIdTypeArray> id_array_tri =
    vtkSmartPointer<vtkIdTypeArray>::New();
  UG->GetIdsOfCellsOfType(VTK_TRIANGLE, id_array_tri);
  vtkSmartPointer<vtkIdTypeArray> id_array_quad =
    vtkSmartPointer<vtkIdTypeArray>::New();
  UG->GetIdsOfCellsOfType(VTK_QUAD, id_array_quad);

  if (id_array_tri->GetNumberOfTuples()!=0 &&
      id_array_quad->GetNumberOfTuples()==0)
    return convert_triangles_to_polyhedron(UG, P);

  if (id_array_quad->GetNumberOfTuples()!=0 &&
      id_array_tri->GetNumberOfTuples()==0)
    return convert_quads_to_polyhedron(UG, P);

  std::cerr << "CELL TYPE NOT HANDLED OR MIXED\n";
  exit(1);
}

} //end of namespace CGAL

#endif // CGAL_VTKUNSTRUCTUREDGRID_TO_POLYHEDRON_3
