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
GridFlags::
GridFlags(const ValueFlags &flags)
{
    /*
     * meas -> lengths
     * points -> lengths
     * general_points -> meas
     * w_meas -> meas
     */
    if (contains(flags, ValueFlags::length))
    {
        fill_lengths_ = true;
    }
    if (contains(flags, ValueFlags::point))
    {
        fill_points_  = true;
        fill_lengths_ = true;
    }
    if (contains(flags, ValueFlags::measure))
    {
        fill_measures_ = true;
        fill_lengths_  = true;
    }
    if (contains(flags, ValueFlags::w_measure))
    {
        fill_measures_ = true;
        fill_w_measures_ = true;
    }
}



bool
GridFlags::
fill_none() const
{
    bool fill_none = true;
    if (fill_lengths_ || fill_points_ || fill_measures_ || fill_w_measures_)
        fill_none = false;
    return fill_none;
}



bool
GridFlags::
fill_points() const
{
    return fill_points_;
}



bool
GridFlags::
points_filled() const
{
    return points_filled_;
}



void
GridFlags::
set_points_filled(const bool status)
{
    points_filled_ = status;
}



bool
GridFlags::
fill_measures() const
{
    return fill_measures_;
}



bool
GridFlags::
measures_filled() const
{
    return measures_filled_;
}



void
GridFlags::
set_measures_filled(const bool status)
{
    measures_filled_ = status;
}



bool
GridFlags::
fill_w_measures() const
{
    return fill_w_measures_;
}



bool
GridFlags::
w_measures_filled() const
{
    return w_measures_filled_;
}



void
GridFlags::
set_w_measures_filled(const bool status)
{
    w_measures_filled_ = status;
}



bool
GridFlags::
fill_lengths() const
{
    return fill_lengths_;
}



bool
GridFlags::
lengths_filled() const
{
    return lengths_filled_;
}



void
GridFlags::
set_lengths_filled(const bool status)
{
    lengths_filled_ = status;
}



void
GridFlags::
print_info(LogStream &out) const
{
    out.begin_item("lengths");
    out << "   fill = " << fill_lengths_ << "    filled = " << lengths_filled_;
    out.end_item();

    out.begin_item("points");
    out << "   fill = " << fill_points_  << "    filled = " << points_filled_;
    out.end_item();

    out.begin_item("measures");
    out << "   fill = " << fill_measures_<< "    filled = " << measures_filled_;
    out.end_item();

    out.begin_item("w_measures");
    out << "   fill = " << fill_w_measures_ << "    filled = " << w_measures_filled_;
    out.end_item();
}

#if 0
//====================================================


/**
 * Exception used when a ValueFlag is not admissibile for the caller object.
 */
DeclException2(ExcFillFlagNotSupported, ValueFlags, ValueFlags,
               << "The passed ValueFlag " << arg2
               << " contains a non admissible flag " << (arg1 ^arg2));

#endif


FunctionFlags::
FunctionFlags(const ValueFlags &flags)
{
    if (contains(flags, ValueFlags::point))
        fill_points_ = true;

    if (contains(flags, ValueFlags::value))
        fill_values_ = true;

    if (contains(flags, ValueFlags::gradient))
        fill_gradients_ = true;

    if (contains(flags, ValueFlags::hessian))
        fill_hessians_ = true;
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


bool
FunctionFlags::
fill_none() const
{
    bool fill_none = true;

    if (fill_values_ || fill_gradients_ || fill_hessians_)
        fill_none = false;

    return fill_none;
}


bool
FunctionFlags::
fill_points() const
{
    return fill_points_;
}



bool
FunctionFlags::
points_filled() const
{
    return points_filled_;
}



void
FunctionFlags::
set_points_filled(const bool status)
{
    points_filled_ = status;
}



bool
FunctionFlags::
fill_values() const
{
    return fill_values_;
}


bool
FunctionFlags::
values_filled() const
{
    return values_filled_;
}

void
FunctionFlags::
set_values_filled(const bool status)
{
    values_filled_ = status;
}

bool
FunctionFlags::
fill_gradients() const
{
    return fill_gradients_;
}

bool
FunctionFlags::
gradients_filled() const
{
    return gradients_filled_;
}

void
FunctionFlags::
set_gradients_filled(const bool status)
{
    gradients_filled_ = status;
}

bool
FunctionFlags::
fill_hessians() const
{
    return fill_hessians_;
}

bool
FunctionFlags::
hessians_filled() const
{
    return hessians_filled_;
}

void
FunctionFlags::
set_hessians_filled(const bool status)
{
    hessians_filled_ = status;
}


void
FunctionFlags::
print_info(LogStream &out) const
{
    out.begin_item("values");
    out << "   fill = " << fill_values_ << "    filled = " << values_filled_;
    out.end_item();

    out.begin_item("gradients");
    out << "   fill = " << fill_gradients_ << "    filled = " << gradients_filled_;
    out.end_item();

    out.begin_item("hessians");
    out << "   fill = " << fill_hessians_ << "    filled = " << hessians_filled_;
    out.end_item();
}

//====================================================







bool
MappingFlags::
fill_none() const
{
    bool fill_none = true;

    if (fill_inv_gradients_ ||
        fill_inv_hessians_ ||
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
        fill_inv_gradients_ = true;

    if (contains(flags, ValueFlags::inv_hessian))
        fill_inv_hessians_ = true;

    if (contains(flags, ValueFlags::measure))
        fill_measures_ = true;

    if (contains(flags, ValueFlags::w_measure))
    {
        fill_measures_ = true;
        fill_w_measures_ = true;
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
    return fill_inv_gradients_;
}

bool
MappingFlags::
inv_gradients_filled() const
{
    return inv_gradients_filled_;
}

void
MappingFlags::
set_inv_gradients_filled(const bool status)
{
    inv_gradients_filled_ = status;
}

bool
MappingFlags::
fill_inv_hessians() const
{
    return fill_inv_hessians_;
}

bool
MappingFlags::
inv_hessians_filled() const
{
    return inv_hessians_filled_;
}

void
MappingFlags::
set_inv_hessians_filled(const bool status)
{
    inv_hessians_filled_ = status;
}
bool
MappingFlags::
fill_measures() const
{
    return fill_measures_;
}



bool
MappingFlags::
measures_filled() const
{
    return measures_filled_;
}



void
MappingFlags::
set_measures_filled(const bool status)
{
    measures_filled_ = status;
}



bool
MappingFlags::
fill_w_measures() const
{
    return fill_w_measures_;
}



bool
MappingFlags::
w_measures_filled() const
{
    return w_measures_filled_;
}



void
MappingFlags::
set_w_measures_filled(const bool status)
{
    w_measures_filled_ = status;
}


void
MappingFlags::
print_info(LogStream &out) const
{
    FunctionFlags::print_info(out);

    out.begin_item("inv gradients");
    out << "   fill = " << fill_inv_gradients_ << "    filled = " << inv_gradients_filled_;
    out.end_item();

    out.begin_item("inv hessians");
    out << "   fill = " << fill_inv_hessians_ << "    filled = " << inv_hessians_filled_;
    out.end_item();



}
//====================================================



#if 0

//====================================================
MappingFaceValueFlagsHandler::
MappingFaceValueFlagsHandler()
    :
    MappingFlags(),
    fill_normals_(false),
    normals_filled_(false)
{}


bool
MappingFaceValueFlagsHandler::
fill_none() const
{
    bool fill_none = true;

    if (fill_normals_ || !MappingFlags::fill_none())
        fill_none = false;

    return fill_none;
}

MappingFaceValueFlagsHandler::
MappingFaceValueFlagsHandler(const ValueFlags &flags)
{
    if (contains(flags, ValueFlags::face_point) ||
        contains(flags, ValueFlags::map_face_value))
    {
        GridFlags::fill_points_ = true;
        fill_values_ = true;
    }

    if (contains(flags, ValueFlags::map_face_gradient))
    {
        fill_gradients_ = true;
    }

    if (contains(flags, ValueFlags::map_face_hessian))
    {
        fill_hessians_ = true;
    }

    if (contains(flags, ValueFlags::map_face_inv_gradient))
    {
        fill_gradients_ = true;
        fill_measures_ = true;
        fill_inv_gradients_ = true;
    }

    if (contains(flags, ValueFlags::map_face_inv_hessian))
    {
        fill_hessians_ = true;
        fill_inv_hessians_ = true;
    }

    if (contains(flags, ValueFlags::face_measure))
    {
        fill_gradients_ = true;
        fill_measures_ = true;
    }

    if (contains(flags, ValueFlags::face_w_measure))
    {
        fill_gradients_ = true;
        fill_measures_ = true;
        fill_w_measures_ = true;
    }

    if (contains(flags, ValueFlags::face_normal))
    {
        fill_normals_ = true;
    }
}



bool
MappingFaceValueFlagsHandler::
fill_normals() const
{
    return fill_normals_;
}

bool
MappingFaceValueFlagsHandler::
normals_filled() const
{
    return normals_filled_;
}

void
MappingFaceValueFlagsHandler::
set_normals_filled(const bool status)
{
    normals_filled_ = status;
}
//====================================================


//====================================================
BasisElemValueFlagsHandler::
BasisElemValueFlagsHandler()
    :
    ValueFlagsHandler(),
    fill_divergences_(false),
    divergences_filled_(false)
{}


bool
BasisElemValueFlagsHandler::
fill_none() const
{
    bool fill_none = true;

    if (fill_divergences_ || !ValueFlagsHandler::fill_none())
        fill_none = false;

    return fill_none;
}



BasisElemValueFlagsHandler::
BasisElemValueFlagsHandler(const ValueFlags &flags)
{
    if (contains(flags, ValueFlags::value))
    {
        fill_values_ = true;
    }

    if (contains(flags, ValueFlags::gradient))
    {
        fill_gradients_ = true;
    }

    if (contains(flags, ValueFlags::hessian))
    {
        fill_hessians_ = true;
    }

    if (contains(flags, ValueFlags::divergence))
    {
        fill_gradients_ = true;
        fill_divergences_ = true;
    }
}



bool
BasisElemValueFlagsHandler::
fill_divergences() const
{
    return fill_divergences_;
}

bool
BasisElemValueFlagsHandler::
divergences_filled() const
{
    return divergences_filled_;
}

void
BasisElemValueFlagsHandler::
set_divergences_filled(const bool status)
{
    divergences_filled_ = status;
}


void
BasisElemValueFlagsHandler::
print_info(LogStream &out) const
{
    ValueFlagsHandler::print_info(out);
    out << "divergences -->    fill = "
        << fill_divergences_ << "    filled = " << divergences_filled_ << std::endl;
}
//====================================================




//====================================================
BasisFaceValueFlagsHandler::
BasisFaceValueFlagsHandler(const ValueFlags &flags)
{
    if (contains(flags, ValueFlags::face_value))
    {
        fill_values_ = true;
    }

    if (contains(flags, ValueFlags::face_gradient))
    {
        fill_gradients_ = true;
    }

    if (contains(flags, ValueFlags::face_hessian))
    {
        fill_hessians_ = true;
    }

    if (contains(flags, ValueFlags::face_divergence))
    {
        fill_gradients_ = true;
        fill_divergences_ = true;
    }
}
//====================================================

#endif

IGA_NAMESPACE_CLOSE




