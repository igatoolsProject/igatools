//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2015  by the igatools authors (see authors.txt).
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

#ifdef PARAVIEW_PLUGIN

#include <paraview_plugin/element_vtk_info.h>

// #include <vtkMultiBlockDataSet.h>
// #include <vtkStructuredGrid.h>
// #include <vtkUnstructuredGrid.h>
// #include <vtkSmartPointer.h>
// #include <vtkInformation.h>
// #include <vtkCellArray.h>
// #include <vtkDoubleArray.h>
// #include <vtkPointData.h>
// #include <vtkHexahedron.h>
// #include <vtkQuad.h>
// #include <vtkLine.h>
// #include <vtkPolyLine.h>
// #include <vtkPolyVertex.h>
// 
// #include <igatools/functions/functions_container.h>
// #include <igatools/base/quadrature_lib.h>
// 
// using std::string;
// using std::shared_ptr;
// using std::pair;
// 

IGA_NAMESPACE_OPEN


template <int dim>
ElementVTKInfo<dim>::
ElementVTKInfo(const TensorSize<dim> &num_vis_elements,
               const vtkGridType &grid_type)
  :
  quadrature_(create_quadrature(num_vis_elements, grid_type)),
  points_map_(create_points_map(num_vis_elements, grid_type)),
  points_mask_(create_points_mask(num_vis_elements, grid_type)),
{
#ifndef NDEBUG
  for (int dir = 0; dir < dim; ++dir)
    Assert(num_vis_elements[dir] > 0,
           ExcMessage("The number of visualization elements in each direction "
                      "must be > 0."));
#endif
};



template <int dim>
auto
ElementVTKInfo<dim>::
create(const TensorSize &num_vis_elements,
       const vtkGridType &grid_type) ->
SelfPtr_
{
  return SelfPtr_(new Self_(num_vis_elements, grid_type));
};



template <int dim>
Size
ElementVTKInfo<dim>::
get_num_points() const
{
  return points_mask_.size();
};



template <int dim>
auto
ElementVTKInfo<dim>::
get_quadrature() const ->
QuadConstPtr_
{
  return quadrature;
};



template <int dim>
auto
ElementVTKInfo<dim>::
get_points_map() const ->
const PointsMap_ &
{
  return points_map;
};



template <int dim>
auto
ElementVTKInfo<dim>::
get_points_mask() const ->
const PointsMask_ &
{
  return points_mask;
};

template <int dim>
QuadConstPtr_
create_quadrature(const TensorSize<dim> &num_vis_elements,
                  const vtkGridType &grid_type)
{
    const bool is_quadratic = grid_type == vtkGridType::UnstructuredQuadratic;

    const Size n_pts_per_cell_dir = is_quadratic ? 3 : 2;

    TensorSize<dim> n_points;
    for (int dir = 0; dir < dim; ++dir)
    {
        Assert(num_vis_elements[dir] > 0, ExcMessage("Wrong number of visualization elements."));
        n_points[dir] = num_vis_elements[dir] * (n_pts_per_cell_dir - 1) + 1;
    }

    QUniform<dim> uniform_quad(n_points);

    return std::make_shared<const Quadrature<dim>> (uniform_quad);

#if 0
    // TODO: Right now is not possible to use non tensor product quadrature
    // for cache elements in igatools. It will be in the future.
    // The code below defines a non tensor product quadrature for VTK quadratic
    // cells.

    const auto &uniform_array_points = uniform_quad.get_points_1d();
    const Size n_total_points = uniform_quad.get_num_points();

    // The points that are not in an edge are removed from the quadrature points.
    // The condition to be in an edge is that just one or zero directions are
    // active. (This also works for the 1D case, where no point is removed.)
    SafeSTLVector<Index> point_indices;
    for (int i_pt = 0; i_pt < n_total_points; ++i_pt)
    {
        const auto t_id = uniform_array_points.flat_to_tensor(i_pt);

        int active_directions = 0;
        for (int dir = 0; dir < dim; ++dir)
            if ((t_id[dir] > 0) && (t_id[dir] < (n_points[dir] - 1)))
                ++active_directions;

        if (active_directions < 2)
            point_indices.push_back(i_pt);
    }

    ValueVector<Points<dim>> quad_points(point_indices.size());
    auto qp = quad_points.begin();
    for (const auto &id : point_indices)
        *qp++ = uniform_quad.get_point(id);

    BBox<dim> box;

    return shared_ptr<Quadrature<dim>> (new Quadrature<dim> (quad_points,
                                                             uniform_quad.get_weights_1d(), box));
#endif

};

IGA_NAMESPACE_CLOSE

#endif // PARAVIEW_PLUGIN
