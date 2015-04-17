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


#include <igatools/base/flags_handler.h>
#include <igatools/base/exceptions.h>



IGA_NAMESPACE_OPEN

//====================================================
#if 0
GridFlags::
GridFlags(const ValueFlags &flags)
{
    if (contains(flags, ValueFlags::point))
        points_flags_.fill_  = true;

    if (contains(flags, ValueFlags::w_measure))
        w_measures_flags_.fill_ = true;
}



bool
GridFlags::
fill_none() const
{
    bool fill_none = true;
    if (points_flags_.fill_ || w_measures_flags_.fill_)
        fill_none = false;
    return fill_none;
}



bool
GridFlags::
fill_points() const
{
    return points_flags_.fill_;
}



bool
GridFlags::
points_filled() const
{
    return points_flags_.filled_;
}



void
GridFlags::
set_points_filled(const bool status)
{
    points_flags_.filled_ = status;
}





bool
GridFlags::
fill_w_measures() const
{
    return w_measures_flags_.fill_;
}



bool
GridFlags::
w_measures_filled() const
{
    return w_measures_flags_.filled_;
}



void
GridFlags::
set_w_measures_filled(const bool status)
{
    w_measures_flags_.filled_ = status;
}




void
GridFlags::
print_info(LogStream &out) const
{
    out.begin_item("points");
    points_flags_.print_info(out);
    out.end_item();

    out.begin_item("w_measures");
    w_measures_flags_.print_info(out);
    out.end_item();
}
#endif

#if 0
//====================================================


/**
 * Exception used when a ValueFlag is not admissibile for the caller object.
 */
DeclException2(ExcFillFlagNotSupported, ValueFlags, ValueFlags,
               << "The passed ValueFlag " << arg2
               << " contains a non admissible flag " << (arg1 ^arg2));

#endif


GridFlags::
GridFlags(const ValueFlags &flags)
    :
    GridFlags()
{
    this->set_fill_status_from_value_flags(flags);
}


FunctionFlags::
FunctionFlags(const ValueFlags &flags)
    :
    FunctionFlags()
{
    this->set_fill_status_from_value_flags(flags);
}


ValueFlags
FunctionFlags::to_grid_flags(const ValueFlags &flags)
{
    ValueFlags transfer_flag = ValueFlags::measure |
                               ValueFlags::w_measure |
                               ValueFlags::boundary_normal;
    ValueFlags g_flag = flags & transfer_flag;
    if (contains(flags, ValueFlags::point) || contains(flags, ValueFlags::value))
    {
        g_flag |= ValueFlags::point;
    }
    return g_flag;
}


//====================================================







bool
MappingFlags::
fill_none() const
{
    bool fill_none = true;

    if (inv_gradients_flags_.fill_ ||
        inv_hessians_flags_.fill_ ||
        !FunctionFlags::fill_none())
        fill_none = false;

    return fill_none;
}


MappingFlags::
MappingFlags(const ValueFlags &flags)
    :
    FunctionFlags::FunctionFlags(to_function_flags(flags))
{
    if (contains(flags, ValueFlags::inv_gradient)    ||
        contains(flags, ValueFlags::boundary_normal) ||
        contains(flags, ValueFlags::curvature))
        inv_gradients_flags_.fill_ = true;

    if (contains(flags, ValueFlags::inv_hessian))
        inv_hessians_flags_.fill_ = true;

    if (contains(flags, ValueFlags::measure))
        measures_flags_.fill_ = true;

    if (contains(flags, ValueFlags::w_measure))
    {
        measures_flags_.fill_ = true;
        w_measures_flags_.fill_ = true;
    }
}



ValueFlags
MappingFlags::to_function_flags(const ValueFlags &flags)
{
    ValueFlags transfer_flag = ValueFlags::measure |
                               ValueFlags::w_measure |
                               ValueFlags::boundary_normal |
                               FunctionFlags::valid_flags;


    ValueFlags f_flag = flags & transfer_flag;

    if (contains(flags, ValueFlags::measure) ||
        contains(flags, ValueFlags::w_measure) ||
        contains(flags, ValueFlags::inv_gradient) ||
        contains(flags, ValueFlags::outer_normal))
        f_flag |=  ValueFlags::gradient;

    if (contains(flags, ValueFlags::inv_hessian) ||
        contains(flags, ValueFlags::curvature))
        f_flag |=  ValueFlags::gradient | ValueFlags::hessian;

    return f_flag;
}



bool
MappingFlags::
fill_inv_gradients() const
{
    return inv_gradients_flags_.fill_;
}

bool
MappingFlags::
inv_gradients_filled() const
{
    return inv_gradients_flags_.filled_;
}

void
MappingFlags::
set_inv_gradients_filled(const bool status)
{
    inv_gradients_flags_.filled_ = status;
}

bool
MappingFlags::
fill_inv_hessians() const
{
    return inv_hessians_flags_.fill_;
}

bool
MappingFlags::
inv_hessians_filled() const
{
    return inv_hessians_flags_.filled_;
}

void
MappingFlags::
set_inv_hessians_filled(const bool status)
{
    inv_hessians_flags_.filled_ = status;
}
bool
MappingFlags::
fill_measures() const
{
    return measures_flags_.fill_;
}



bool
MappingFlags::
measures_filled() const
{
    return measures_flags_.filled_;
}



void
MappingFlags::
set_measures_filled(const bool status)
{
    measures_flags_.filled_ = status;
}



bool
MappingFlags::
fill_w_measures() const
{
    return w_measures_flags_.fill_;
}



bool
MappingFlags::
w_measures_filled() const
{
    return w_measures_flags_.filled_;
}



void
MappingFlags::
set_w_measures_filled(const bool status)
{
    w_measures_flags_.filled_ = status;
}


void
MappingFlags::
print_info(LogStream &out) const
{
    FunctionFlags::print_info(out);

    out.begin_item("inv gradients");
    inv_gradients_flags_.print_info(out);
    out.end_item();

    out.begin_item("inv hessians");
    inv_hessians_flags_.print_info(out);
    out.end_item();

    out.begin_item("measures");
    measures_flags_.print_info(out);
    out.end_item();

    out.begin_item("w * measures");
    w_measures_flags_.print_info(out);
    out.end_item();
}
//====================================================




IGA_NAMESPACE_CLOSE




