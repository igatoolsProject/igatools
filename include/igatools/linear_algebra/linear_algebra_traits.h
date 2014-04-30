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

#ifndef LINEAR_ALGEBRA_TRAITS_H_
#define LINEAR_ALGEBRA_TRAITS_H_

#include <igatools/base/config.h>

IGA_NAMESPACE_OPEN

template <LinearAlgebraPackage linear_algebra_package>
struct LinearAlgebraTraits
{
    using MatrixType = Matrix<linear_algebra_package>;
    using VectorType = Vector<linear_algebra_package>;
    using LinearSolverType = LinearSolver<linear_algebra_package>;
};



IGA_NAMESPACE_CLOSE

#endif // #ifndef LINEAR_ALGEBRA_TRAITS_H_
