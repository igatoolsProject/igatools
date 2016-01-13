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

/*
 *  Test for BSpline constructors
 *
 *  author: pauletti
 *  date: 2014-10-23
 *
 */

#include "../tests.h"

#include <igatools/basis_functions/bspline.h>


template <int dim>
void serialize_deserialize(const std::shared_ptr<BSpline<dim>> space_in)
{
  std::shared_ptr<ReferenceSpaceBasis<dim>> basis = space_in;
//  auto basis = std::make_shared<SpaceBase<dim>>(space_in->get_grid());
//  std::shared_ptr<Space<dim,0,1,1,Transformation::h_grad>> basis = space_in;
//  std::shared_ptr<Space<dim,0,1,1,Transformation::h_grad>> basis = std::make_shared<TTT<dim>>();
  out.begin_item("Original BSpline:");
  basis->print_info(out);
  out.end_item();


  std::string filename = "bspline_space_dim" + std::to_string(dim) + ".xml";
  std::string tag_name = "BSpline_dim" + std::to_string(dim);
  {
    // serialize the BSpline object to an xml file
    std::ofstream xml_ostream(filename);
    OArchive xml_out(xml_ostream);

    xml_out << basis;
  }

  basis.reset();
  {
    // de-serialize the BSpline object from an xml file
    std::ifstream xml_istream(filename);
    IArchive xml_in(xml_istream);

    xml_in >> basis;
  }
  out.begin_item("BSpline after serialize-deserialize:");
  basis->print_info(out);
  out.end_item();

}



namespace grid
{
template<int dim>
shared_ptr<Grid<dim> >
uniform(const int n_knots)
{
  return Grid<dim>::create(n_knots);
}


};


template<int dim>
void uniform_degree(const int deg, shared_ptr<Grid<dim>> grid)
{
  OUTSTART
  auto space = SplineSpace<dim>::create(deg, grid);
  auto basis = BSpline<dim>::create(space);
//  std::shared_ptr<BSpline<dim>> basis = BSpline<dim>::create(deg, grid);

  serialize_deserialize(basis);

  OUTEND
}


template<int dim>
void direction_degree(const TensorIndex<dim> &deg,
                      shared_ptr<Grid<dim>> grid)
{
  OUTSTART
  auto space = SplineSpace<dim>::create(deg, grid);
  auto basis = BSpline<dim>::create(space);
//  std::shared_ptr<BSpline<dim>> basis = BSpline<dim>::create(deg, grid);

  serialize_deserialize(basis);

  OUTEND
}


int main()
{
  const int deg = 1;
  const int n_knots = 2;
//  uniform_degree<0>(deg, grid::uniform<0>(n_knots));
  uniform_degree<1>(deg, grid::uniform<1>(n_knots));
  uniform_degree<2>(deg, grid::uniform<2>(n_knots));
  uniform_degree<3>(deg, grid::uniform<3>(n_knots));

//  TensorIndex<0> deg0;
//  direction_degree<0>(deg0, grid::uniform<0>(n_knots));

  TensorIndex<1> deg1 = {1};
  direction_degree<1>(deg1, grid::uniform<1>(n_knots));

  TensorIndex<2> deg2 = {2,3};
  direction_degree<2>(deg2, grid::uniform<2>(n_knots));

  TensorIndex<3> deg3 = {3,4,5};
  direction_degree<3>(deg3, grid::uniform<3>(n_knots));

  return 0;
}
