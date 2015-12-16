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


#define USE_IGATOOLS

#ifdef USE_IGATOOLS
#include <igatools/geometry/domain.h>
#include <igatools/geometry/grid_function_lib.h>
#include <igatools/io/writer.h>
#include "../tests.h"
#endif

#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include "vtkUnstructuredGridOverlay_3.h"


using std::string;
using std::to_string;

#ifdef USE_IGATOOLS
template <int dim>
void create_grid_files(const std::array<string,2> &names,const std::array<int,dim> &n_knots)
{
  using IdentityFunc = grid_functions::IdentityGridFunction<dim>;
  using LinearFunc = grid_functions::LinearGridFunction<dim,dim>;

  const int n_pts_dir = 2;

  {
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
  }


  {
    auto grid_b = Grid<dim>::create(n_knots[1]);
    auto func_b = IdentityFunc::create(grid_b);
    auto domain_b = Domain<dim,0>::create(func_b);

    Writer<dim,0> writer(domain_b,n_pts_dir);
    writer.save(names[1]);
  }
}
#endif

template<int dim>
void test_overlay()
{
  std::array<int,dim> n_knots;
  n_knots[0] = 2;
  n_knots[1] = 3;

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

