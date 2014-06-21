//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2014  by the igatools authors (see authors.txt).
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

#include <igatools/base/function.h>
#include <igatools/base/exceptions.h>

using std::vector;

IGA_NAMESPACE_OPEN

template< int dim, int range, int rank >
Function< dim, range, rank >::~Function()
{}



template< int dim, int range, int rank >
void Function< dim, range, rank >::
evaluate_gradients(const vector<PointType> & ,
                   vector<Gradient> &) const
{
    Assert(false, ExcNotImplemented());
}



template< int dim, int range, int rank >
void Function< dim, range, rank >::
evaluate_hessians(const vector<PointType> &,
                  vector<Hessian> &) const
{
    Assert(false, ExcNotImplemented());
}



template< int dim, int range, int rank >
void Function< dim, range, rank >::
evaluate_values_and_gradients(
    const vector<PointType> &points,
    vector<Value> &values,
    vector<Gradient> &gradients) const
{
    this->evaluate(points, values) ;
    this->evaluate_gradients(points, gradients) ;
}


IGA_NAMESPACE_CLOSE

#include <igatools/base/function.inst>

