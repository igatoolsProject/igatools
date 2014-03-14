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


#include <igatools/base/value_flags_handler.h>
#include <igatools/base/exceptions.h>


IGA_NAMESPACE_OPEN

//====================================================
ValueFlagsHandler::
ValueFlagsHandler()
    :
    fill_values_(false),
    fill_gradients_(false),
    fill_hessians_(false)
{}


bool
ValueFlagsHandler::
fill_values() const
{
    return fill_values_;
}

bool
ValueFlagsHandler::
fill_gradients() const
{
    return fill_gradients_;
}

bool
ValueFlagsHandler::
fill_hessians() const
{
    return fill_hessians_;
}
//====================================================



//====================================================
GridElemValueFlagsHandler::
GridElemValueFlagsHandler()
    :
    fill_points_(false),
    fill_measures_(false),
    fill_w_measures_(false)
{}

bool
GridElemValueFlagsHandler::
fill_points() const
{
    return fill_points_;
}


bool
GridElemValueFlagsHandler::
fill_measures() const
{
    return fill_measures_;
}


bool
GridElemValueFlagsHandler::
fill_w_measures() const
{
    return fill_w_measures_;
}
//====================================================



//====================================================
MappingValueFlagsHandler::
MappingValueFlagsHandler()
    :
    ValueFlagsHandler(),
    fill_inv_gradients_(false),
    fill_inv_hessians_(false)
{}


MappingValueFlagsHandler::
MappingValueFlagsHandler(const ValueFlags &flags)
{
    Assert(false,ExcNotImplemented());
    AssertThrow(false,ExcNotImplemented());
}


bool
MappingValueFlagsHandler::
fill_inv_gradients() const
{
    return fill_inv_gradients_;
}


bool
MappingValueFlagsHandler::
fill_inv_hessians() const
{
    return fill_inv_hessians_;
}
//====================================================











IGA_NAMESPACE_CLOSE




