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

#include <igatools/basis_functions/values1d_const_view.h>.h>
#include <igatools/base/exceptions.h>


IGA_NAMESPACE_OPEN

Values1DConstView::
Values1DConstView(const DenseMatrix &funcs,const Index func_id)
    :
    funcs_(&funcs),
    func_id_(func_id)
{
    Assert(func_id >= 0 && func_id < Size(funcs_->size1()),
           ExcIndexRange(func_id,0,Size(funcs_->size1())))
}


Size
Values1DConstView::
get_num_points() const
{
	Assert(funcs_ != nullptr,ExcNullPtr());
	return funcs_->size2();
}

Real
Values1DConstView::
operator()(const Index point_id) const
{
	Assert(funcs_ != nullptr,ExcNullPtr());
    Assert(point_id >= 0 && point_id < this->get_num_points(),
           ExcIndexRange(point_id,0,this->get_num_points()));

    return (*funcs_)(func_id_,point_id);
}

IGA_NAMESPACE_CLOSE
