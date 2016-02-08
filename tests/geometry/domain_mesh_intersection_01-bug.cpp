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


#include <igatools/geometry/domain.h>
#include <igatools/functions/grid_function_lib.h>
#include <igatools/io/writer.h>
#include "../tests.h"

#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkQuad.h>
//#include <vtkPolyVertex.h>
//#include <vtkCellArray.h>
//#include <vtkPolyLine.h>

#include "vtkUnstructuredGridOverlay_3.h"


using std::string;
using std::to_string;




void write_grid(const string &filename,const Grid<2> &grid)
{
  auto vtk_grid = vtkSmartPointer <vtkUnstructuredGrid>::New();

  const int dim = 2;
  const auto &n_pts_dir = grid.get_num_knots_dim();
  TensorSizedContainer<dim> pts_t_size(n_pts_dir);


  //---------------------------------------------------------------------------
  // Filling the points -- begin
  const int n_pts = n_pts_dir.flat_size();
  auto points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(n_pts);

  double pt_coords[3] = { 0.0, 0.0, 0.0 };

  TensorIndex<dim> pt_t_id;

  for (int pt_f_id = 0 ; pt_f_id < n_pts ; ++pt_f_id)
  {
    out << "pt_f_id: " << pt_f_id << endl;
    pt_t_id = pts_t_size.flat_to_tensor(pt_f_id);

    for (int dir = 0 ; dir < dim ; ++dir)
    {
      const int knt_id = pt_t_id[dir];
      pt_coords[dir] = grid.get_knot_coordinates(dir)[knt_id];
    }
    points->SetPoint(pt_f_id,pt_coords);
  }
  vtk_grid->SetPoints(points);
  // Filling the points -- end
  //---------------------------------------------------------------------------



  //---------------------------------------------------------------------------
  // Filling the connectivity -- begin
  const int n_vertices = 4;
  using T_ = SafeSTLArray<SafeSTLArray<int,dim>,n_vertices>;
  const T_ delta_idx
  {
    { 0, 0 },
    { 1, 0 },
    { 1, 1 },
    { 0, 1 }};


  for (const auto &elem : grid)
  {
    const auto pt_offset_t_id = elem.get_index().get_tensor_index();

    auto vtk_quad = vtkSmartPointer<vtkQuad>::New();

    for (int vertex = 0 ; vertex < n_vertices ; ++vertex)
    {
      const auto &delta_idx_vertex = delta_idx[vertex];
      for (int dir = 0 ; dir < dim ; ++dir)
        pt_t_id[dir] = pt_offset_t_id[dir] + delta_idx_vertex[dir];

      const auto pt_f_id = pts_t_size.tensor_to_flat(pt_t_id);
      out << "pt_f_id: " << pt_f_id << endl;

      vtk_quad->GetPointIds()->SetId(vertex,pt_f_id);
    }
    vtk_grid->InsertNextCell(vtk_quad->GetCellType(),vtk_quad->GetPointIds());
  }
  // Filling the connectivity -- end
  //---------------------------------------------------------------------------


  auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetDataModeToAscii();
  string outputname = filename + ".vtu";
  writer->SetFileName(outputname.c_str());
  writer->SetInputData(vtk_grid);
  writer->Write();

}


template<int codim>
void write_domain(const string &filename, const Domain<2,codim> &domain)
{
  auto domain_handler = domain.create_cache_handler();

  domain_handler->set_element_flags(domain_element::Flags::point);

  auto elem = domain->begin();
  const auto end = domain->end();

//  const auto quad = QTrapez<2>::create();
//  domain_handler->init_element_cache(elem,quad);

  Assert(false,ExcNotImplemented());

}


void
create_file_grid_to_project(const string &filename, const int n_knots)
{
  TensorSize<2> n_pts_dir(2);

  //---------------------------------------------------------------------------
  // Filling the points -- begin
  const int n_pts = n_pts_dir.flat_size();
  auto points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(n_pts);

  SafeSTLArray<SafeSTLVector<Real>,2> knots
  {
    {0.1, 0.5, 0.9},
    {0.1, 0.5, 0.9}};

  auto grid = Grid<2>::create(knots);

  write_grid(filename,*grid);
#if 0
  double pt_coords[3] = { 0.0, 0.0, 0.0 };

  pt_coords[0] = 0.1;
  pt_coords[1] = 0.1;
  points->SetPoint(0,pt_coords);

  pt_coords[0] = 0.9;
  pt_coords[1] = 0.1;
  points->SetPoint(1,pt_coords);

  pt_coords[0] = 0.1;
  pt_coords[1] = 0.9;
  points->SetPoint(2,pt_coords);

  pt_coords[0] = 0.9;
  pt_coords[1] = 0.9;
  points->SetPoint(3,pt_coords);
  // Filling the points -- end
  //---------------------------------------------------------------------------



  auto grid = vtkSmartPointer <vtkUnstructuredGrid>::New();
  grid->SetPoints(points);


  const int n_elems = 1;
  auto cell_ids = vtkSmartPointer<vtkIdTypeArray>::New();

  const int n_pts_per_elem = 4;
  const int tuple_size = n_pts_per_elem+1;
  cell_ids->SetNumberOfComponents(tuple_size);
  cell_ids->SetNumberOfTuples(n_elems);


  std::vector<vtkIdType> elem_conn(tuple_size);
  elem_conn[0] = n_pts_per_elem;


  auto vtk_quad = vtkSmartPointer<vtkQuad>::New();
  vtk_quad->GetPointIds()->SetId(0,0);
  vtk_quad->GetPointIds()->SetId(1,1);
  vtk_quad->GetPointIds()->SetId(2,3);
  vtk_quad->GetPointIds()->SetId(3,2);

  grid->InsertNextCell(vtk_quad->GetCellType(),vtk_quad->GetPointIds());






//  std::string outputname = "kkk.vtu";
  auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetDataModeToAscii();
  string outputname = filename + ".vtu";
  writer->SetFileName(outputname.c_str());
  writer->SetInputData(grid);
  writer->Write();
#endif

//  AssertThrow(false,ExcNotImplemented());
}

template <int dim>
void create_grid_files(const std::array<string,2> &names,const std::array<int,dim> &n_knots)
{
  using IdentityFunc = grid_functions::IdentityGridFunction<dim>;
  using LinearFunc = grid_functions::LinearGridFunction<dim,dim>;

  const int n_pts_dir = 2;

  {
    create_file_grid_to_project(names[0],n_knots[0]);
    /*
        auto grid_a = Grid<dim>::create(n_knots[0]);

        typename LinearFunc::template Derivative<1> A;
        typename LinearFunc::Value b;
        const Real eps = 0.1;
        for (int i = 0 ; i < dim ; ++i)
        {
          A[i][i] = 1.0 - 2.0 * eps;
          b[i] = eps;
        }

        auto func_a = LinearFunc::create(grid_a,A,b);
        auto domain_a = Domain<dim,0>::create(func_a);

        Writer<dim,0> writer(domain_a,n_pts_dir);
        writer.save(names[0]);
        //*/
  }


  {
    auto grid_b = Grid<dim>::create(n_knots[1]);
    /*
    auto func_b = IdentityFunc::create(grid_b);
    auto domain_b = Domain<dim,0>::create(func_b);

    Writer<dim,0> writer(domain_b,n_pts_dir);
    writer.save(names[1]);
    //*/

    write_grid(names[1],*grid_b);
  }
}

template<int dim>
void test_overlay()
{
  std::array<int,dim> n_knots;
  n_knots[0] = 3;
  n_knots[1] = 4;

  std::array<string,2> names;
  names[0] = "domain_to_proj_" + to_string(n_knots[0]-1) + "_elems";
  names[1] = "domain_support_" + to_string(n_knots[1]-1) + "_elems";

  create_grid_files<dim>(names,n_knots);


  using std::cout;
  using std::endl;

  //read all the data from the files
  using ReaderVTU = vtkXMLUnstructuredGridReader;

  using VTKUGridPtr = vtkSmartPointer<vtkUnstructuredGrid>;
  std::array<VTKUGridPtr,2> unstruct_grids;

  for (int i = 0 ; i <= 1 ; ++i)
  {
    const string filename = names[i] + ".vtu";
    cout << "Reading the file \"" << filename << "\"...";
    auto reader = vtkSmartPointer<ReaderVTU>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();
    unstruct_grids[i] = reader->GetOutput();
    int n_elems = unstruct_grids[i]->GetNumberOfCells();
    cout << "done (n_elements = " << n_elems << ")" << endl;
  }


  auto grid_to_project = unstruct_grids[0];
  auto grid_support = unstruct_grids[1];

  const double mesh_distance = -1.0;
  const bool partial_overlay = true;
  cout << "Computing overlay...";
  VTKUGridPtr overlay =
    vtkUnstructuredGridOverlay_3(grid_support, grid_to_project, mesh_distance, partial_overlay, false);
  int n_elems = overlay->GetNumberOfCells();
  cout << "done (n_elements = " << n_elems << ")" << endl;

  /*
  //This compute the "parent" vertices. Useful in Numea to avoid multiple interpolation
  cout << "BuildParentVertices...";
  BuildParentVertices(res, ugA.GetPointer(), ugB.GetPointer(), 0.0000001);
  cout << "done" << endl;
  //*/

  // Write the unstructured grid (different treatment for .vtk and .vtk file)
  std::string outputname = "overlay_" + names[0] + "_" + names[1] + ".vtu";
  auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetDataModeToAscii();
  writer->SetFileName(outputname.c_str());
  writer->SetInputData(overlay);
  writer->Write();

}


int main()
{
//  test_overlay<1>();
  test_overlay<2>();
//  test_overlay<3>();

  return 0;
}

