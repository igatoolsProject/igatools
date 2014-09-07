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


/**
 * Exception used when a ValueFlag is not admissibile for the caller object.
 */
DeclException2(ExcFillFlagNotSupported, ValueFlags, ValueFlags,
               << "The passed ValueFlag " << arg2
               << " contains a non admissible flag " << (arg1 ^arg2));




//====================================================
ValueFlagsHandler::
ValueFlagsHandler()
    :
    fill_values_(false),
    values_filled_(false),
    fill_gradients_(false),
    gradients_filled_(false),
    fill_hessians_(false),
    hessians_filled_(false)
{}

bool
ValueFlagsHandler::
fill_none() const
{
    bool fill_none = true;

    if (fill_values_ || fill_gradients_ || fill_hessians_)
        fill_none = false;

    return fill_none;
}


bool
ValueFlagsHandler::
fill_values() const
{
    return fill_values_;
}


bool
ValueFlagsHandler::
values_filled() const
{
    return values_filled_;
}

void
ValueFlagsHandler::
set_values_filled(const bool status)
{
    values_filled_ = status;
}

bool
ValueFlagsHandler::
fill_gradients() const
{
    return fill_gradients_;
}

bool
ValueFlagsHandler::
gradients_filled() const
{
    return gradients_filled_;
}

void
ValueFlagsHandler::
set_gradients_filled(const bool status)
{
    gradients_filled_ = status;
}

bool
ValueFlagsHandler::
fill_hessians() const
{
    return fill_hessians_;
}

bool
ValueFlagsHandler::
hessians_filled() const
{
    return hessians_filled_;
}

void
ValueFlagsHandler::
set_hessians_filled(const bool status)
{
    hessians_filled_ = status;
}


void
ValueFlagsHandler::
print_info(LogStream &out) const
{
    using std::endl;

    out.begin_item("Flags:");
    out << "   values -->    fill = "
        << fill_values_ << "    filled = " << values_filled_ << endl;
    out << "gradients -->    fill = "
        << fill_gradients_ << "    filled = " << gradients_filled_ << endl;
    out << " hessians -->    fill = "
        << fill_hessians_ << "    filled = " << hessians_filled_ << endl;
    out.end_item();
}

//====================================================



//====================================================
GridElemValueFlagsHandler::
GridElemValueFlagsHandler()
    :
    fill_points_(false),
    points_filled_(false),
    fill_measures_(false),
    measures_filled_(false),
    fill_w_measures_(false),
    w_measures_filled_(false)
{}

bool
GridElemValueFlagsHandler::
fill_none() const
{
    bool fill_none = true;

    if (fill_points_ || fill_measures_ || fill_w_measures_)
        fill_none = false;

    return fill_none;
}

GridElemValueFlagsHandler::
GridElemValueFlagsHandler(const ValueFlags &flags)
{
    if (contains(flags, ValueFlags::point))
    {
        fill_points_ = true;
    }

    if (contains(flags, ValueFlags::measure))
    {
        fill_measures_ = true ;
    }

    if (contains(flags, ValueFlags::w_measure))
    {
        fill_measures_ = true ;
        fill_w_measures_ = true ;
    }
}

bool
GridElemValueFlagsHandler::
fill_points() const
{
    return fill_points_;
}

bool
GridElemValueFlagsHandler::
points_filled() const
{
    return points_filled_;
}

void
GridElemValueFlagsHandler::
set_points_filled(const bool status)
{
    points_filled_ = status;
}

bool
GridElemValueFlagsHandler::
fill_measures() const
{
    return fill_measures_;
}

bool
GridElemValueFlagsHandler::
measures_filled() const
{
    return measures_filled_;
}

void
GridElemValueFlagsHandler::
set_measures_filled(const bool status)
{
    measures_filled_ = status;
}

bool
GridElemValueFlagsHandler::
fill_w_measures() const
{
    return fill_w_measures_;
}

bool
GridElemValueFlagsHandler::
w_measures_filled() const
{
    return w_measures_filled_;
}

void
GridElemValueFlagsHandler::
set_w_measures_filled(const bool status)
{
    w_measures_filled_ = status;
}

void
GridElemValueFlagsHandler::
print_info(LogStream &out) const
{
    using std::endl;

    out.begin_item("Fill Flags:");
    out << "    points -->    fill = "
        << fill_points_ << "    filled = " << points_filled_ << endl;
    out << "  measures -->    fill = "
        << fill_measures_ << "    filled = " << measures_filled_ << endl;
    out << "w_measures -->    fill = "
        << fill_w_measures_ << "    filled = " << w_measures_filled_ << endl;
    out.end_item();
}

//====================================================


//====================================================
GridFaceValueFlagsHandler::
GridFaceValueFlagsHandler()
    :
    GridElemValueFlagsHandler(),
    fill_normals_(false),
    normals_filled_(false)
{}

bool
GridFaceValueFlagsHandler::
fill_none() const
{
    bool fill_none = true;

    if (fill_normals_ || !GridElemValueFlagsHandler::fill_none())
        fill_none = false;

    return fill_none;
}


GridFaceValueFlagsHandler::
GridFaceValueFlagsHandler(const ValueFlags &flags)
{
    if (contains(flags, ValueFlags::face_point))
        fill_points_ = true;

    if (contains(flags, ValueFlags::face_measure))
        fill_measures_ = true ;

    if (contains(flags, ValueFlags::face_w_measure))
    {
        fill_measures_ = true ;
        fill_w_measures_ = true ;
    }

    if (contains(flags, ValueFlags::face_normal))
    {
        fill_normals_ = true ;
    }
}


bool
GridFaceValueFlagsHandler::
fill_normals() const
{
    return fill_normals_;
}

bool
GridFaceValueFlagsHandler::
normals_filled() const
{
    return normals_filled_;
}

void
GridFaceValueFlagsHandler::
set_normals_filled(const bool status)
{
    normals_filled_ = status;
}
//====================================================


//====================================================
MappingElemValueFlagsHandler::
MappingElemValueFlagsHandler()
    :
    ValueFlagsHandler(),
    GridElemValueFlagsHandler(),
    fill_inv_gradients_(false),
    inv_gradients_filled_(false),
    fill_inv_hessians_(false),
    inv_hessians_filled_(false)
{}


bool
MappingElemValueFlagsHandler::
fill_none() const
{
    bool fill_none = true;

    if (fill_inv_gradients_ ||
        fill_inv_hessians_ ||
        !ValueFlagsHandler::fill_none() ||
        !GridElemValueFlagsHandler::fill_none())
        fill_none = false;

    return fill_none;
}

MappingElemValueFlagsHandler::
MappingElemValueFlagsHandler(const ValueFlags &flags)
{
    if (contains(flags, ValueFlags::point) ||
        contains(flags, ValueFlags::map_value))
    {
        fill_points_ = true;
        fill_values_ = true;
    }

    if (contains(flags, ValueFlags::map_gradient))
    {
        fill_gradients_ = true ;
    }

    if (contains(flags, ValueFlags::map_hessian))
    {
        fill_hessians_ = true ;
    }

    if (contains(flags, ValueFlags::map_inv_gradient))
    {
        fill_gradients_ = true ;
        fill_measures_ = true ;
        fill_inv_gradients_ = true ;
    }

    if (contains(flags, ValueFlags::map_inv_hessian))
    {
        fill_hessians_ = true ;
        fill_inv_hessians_ = true ;
    }

    if (contains(flags, ValueFlags::measure))
    {
        fill_gradients_ = true ;
        fill_measures_ = true ;
    }

    if (contains(flags, ValueFlags::w_measure))
    {
        fill_gradients_ = true ;
        fill_measures_ = true ;
        fill_w_measures_ = true ;
    }
}


bool
MappingElemValueFlagsHandler::
fill_inv_gradients() const
{
    return fill_inv_gradients_;
}

bool
MappingElemValueFlagsHandler::
inv_gradients_filled() const
{
    return inv_gradients_filled_;
}

void
MappingElemValueFlagsHandler::
set_inv_gradients_filled(const bool status)
{
    inv_gradients_filled_ = status;
}

bool
MappingElemValueFlagsHandler::
fill_inv_hessians() const
{
    return fill_inv_hessians_;
}

bool
MappingElemValueFlagsHandler::
inv_hessians_filled() const
{
    return inv_hessians_filled_;
}

void
MappingElemValueFlagsHandler::
set_inv_hessians_filled(const bool status)
{
    inv_hessians_filled_ = status;
}


void
MappingElemValueFlagsHandler::
print_info(LogStream &out) const
{
	ValueFlagsHandler::print_info(out);
    GridElemValueFlagsHandler::print_info(out);
}
//====================================================





//====================================================
MappingFaceValueFlagsHandler::
MappingFaceValueFlagsHandler()
    :
    MappingElemValueFlagsHandler(),
    fill_normals_(false),
    normals_filled_(false)
{}


bool
MappingFaceValueFlagsHandler::
fill_none() const
{
    bool fill_none = true;

    if (fill_normals_ || !MappingElemValueFlagsHandler::fill_none())
        fill_none = false;

    return fill_none;
}

MappingFaceValueFlagsHandler::
MappingFaceValueFlagsHandler(const ValueFlags &flags)
{
    if (contains(flags, ValueFlags::face_point) ||
        contains(flags, ValueFlags::map_face_value))
    {
        fill_points_ = true;
        fill_values_ = true;
    }

    if (contains(flags, ValueFlags::map_face_gradient))
    {
        fill_gradients_ = true ;
    }

    if (contains(flags, ValueFlags::map_face_hessian))
    {
        fill_hessians_ = true ;
    }

    if (contains(flags, ValueFlags::map_face_inv_gradient))
    {
        fill_gradients_ = true ;
        fill_measures_ = true ;
        fill_inv_gradients_ = true ;
    }

    if (contains(flags, ValueFlags::map_face_inv_hessian))
    {
        fill_hessians_ = true ;
        fill_inv_hessians_ = true ;
    }

    if (contains(flags, ValueFlags::face_measure))
    {
        fill_gradients_ = true ;
        fill_measures_ = true ;
    }

    if (contains(flags, ValueFlags::face_w_measure))
    {
        fill_gradients_ = true ;
        fill_measures_ = true ;
        fill_w_measures_ = true ;
    }

    if (contains(flags, ValueFlags::face_normal))
    {
        fill_normals_ = true ;
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
        fill_gradients_ = true ;
    }

    if (contains(flags, ValueFlags::hessian))
    {
        fill_hessians_ = true ;
    }

    if (contains(flags, ValueFlags::divergence))
    {
        fill_gradients_ = true ;
        fill_divergences_ = true ;
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
        fill_gradients_ = true ;
    }

    if (contains(flags, ValueFlags::face_hessian))
    {
        fill_hessians_ = true ;
    }

    if (contains(flags, ValueFlags::face_divergence))
    {
        fill_gradients_ = true ;
        fill_divergences_ = true ;
    }
}
//====================================================


IGA_NAMESPACE_CLOSE




