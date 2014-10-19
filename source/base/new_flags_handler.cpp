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


#include <igatools/base/new_flags_handler.h>
#include <igatools/base/exceptions.h>

IGA_NAMESPACE_OPEN

//====================================================
GridFlags::
GridFlags(const NewValueFlags &flags)
{
    /*
     * meas -> lengths
     * points -> lengths
     * general_points -> meas
     * w_meas -> meas
     */
    if (contains(flags, NewValueFlags::point))
    {
        fill_points_  = true;
        fill_lengths_ = true;
    }
    if (contains(flags, NewValueFlags::measure))
    {
        fill_measures_ = true;
        fill_lengths_  = true;
    }
    if (contains(flags, NewValueFlags::w_measure))
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
    if (fill_points_ || fill_measures_ || fill_w_measures_)
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
    using std::endl;
    out.push("  ");
    out << "points     -->  fill = "
        << fill_points_ << "    filled = " << points_filled_ << endl;
    out << "measures   -->  fill = "
        << fill_measures_ << "    filled = " << measures_filled_ << endl;
    out << "w_measures -->  fill = "
        << fill_w_measures_ << "    filled = " << w_measures_filled_;
    out.pop();
}

#if 0
//====================================================


/**
 * Exception used when a ValueFlag is not admissibile for the caller object.
 */
DeclException2(ExcFillFlagNotSupported, NewValueFlags, NewValueFlags,
               << "The passed ValueFlag " << arg2
               << " contains a non admissible flag " << (arg1 ^arg2));



ValueFlagsHandler::
ValueFlagsHandler(const NewValueFlags &flags)
{
    if (contains(flags, NewValueFlags::point))
        fill_points_ = true;

    if (contains(flags, NewValueFlags::value))
        fill_values_ = true;

    if (contains(flags, NewValueFlags::gradient))
        fill_gradients_ = true;

    if (contains(flags, NewValueFlags::hessian))
        fill_hessians_ = true;
}

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
fill_points() const
{
    return fill_points_;
}



bool
ValueFlagsHandler::
points_filled() const
{
    return points_filled_;
}



void
ValueFlagsHandler::
set_points_filled(const bool status)
{
    points_filled_ = status;
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
    out << "   values -->    fill = "
        << fill_values_ << "    filled = " << values_filled_ << endl;
    out << "gradients -->    fill = "
        << fill_gradients_ << "    filled = " << gradients_filled_ << endl;
    out << " hessians -->    fill = "
        << fill_hessians_ << "    filled = " << hessians_filled_ << endl;

}

//====================================================





//====================================================
GridFaceValueFlagsHandler::
GridFaceValueFlagsHandler()
    :
    GridFlags(),
    fill_normals_(false),
    normals_filled_(false)
{}

bool
GridFaceValueFlagsHandler::
fill_none() const
{
    bool fill_none = true;

    if (fill_normals_ || !GridFlags::fill_none())
        fill_none = false;

    return fill_none;
}


GridFaceValueFlagsHandler::
GridFaceValueFlagsHandler(const NewValueFlags &flags)
{
    if (contains(flags, NewValueFlags::face_point))
    {
        fill_points_  = true;
        fill_lengths_ = true;
    }
    if (contains(flags, NewValueFlags::face_measure))
    {
        fill_measures_ = true;
        fill_lengths_  = true;
    }
    if (contains(flags, NewValueFlags::face_w_measure))
    {
        fill_measures_ = true;
        fill_lengths_  = true;
        fill_w_measures_ = true;
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
    GridFlags(),
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
        !GridFlags::fill_none())
        fill_none = false;

    return fill_none;
}

MappingElemValueFlagsHandler::
MappingElemValueFlagsHandler(const NewValueFlags &flags)
{
    if (contains(flags, NewValueFlags::point) ||
        contains(flags, NewValueFlags::map_value))
    {
        GridFlags::fill_points_ = true;
        fill_values_ = true;
    }

    if (contains(flags, NewValueFlags::map_gradient))
    {
        fill_gradients_ = true;
    }

    if (contains(flags, NewValueFlags::map_hessian))
    {
        fill_hessians_ = true;
    }

    if (contains(flags, NewValueFlags::map_inv_gradient))
    {
        fill_gradients_ = true;
        fill_measures_ = true;
        fill_inv_gradients_ = true;
    }

    if (contains(flags, NewValueFlags::map_inv_hessian))
    {
        fill_hessians_ = true;
        fill_inv_hessians_ = true;
    }

    if (contains(flags, NewValueFlags::measure))
    {
        fill_gradients_ = true;
        fill_measures_ = true;
    }

    if (contains(flags, NewValueFlags::w_measure))
    {
        fill_gradients_ = true;
        fill_measures_ = true;
        fill_w_measures_ = true;
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
    GridFlags::print_info(out);
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
MappingFaceValueFlagsHandler(const NewValueFlags &flags)
{
    if (contains(flags, NewValueFlags::face_point) ||
        contains(flags, NewValueFlags::map_face_value))
    {
        GridFlags::fill_points_ = true;
        fill_values_ = true;
    }

    if (contains(flags, NewValueFlags::map_face_gradient))
    {
        fill_gradients_ = true;
    }

    if (contains(flags, NewValueFlags::map_face_hessian))
    {
        fill_hessians_ = true;
    }

    if (contains(flags, NewValueFlags::map_face_inv_gradient))
    {
        fill_gradients_ = true;
        fill_measures_ = true;
        fill_inv_gradients_ = true;
    }

    if (contains(flags, NewValueFlags::map_face_inv_hessian))
    {
        fill_hessians_ = true;
        fill_inv_hessians_ = true;
    }

    if (contains(flags, NewValueFlags::face_measure))
    {
        fill_gradients_ = true;
        fill_measures_ = true;
    }

    if (contains(flags, NewValueFlags::face_w_measure))
    {
        fill_gradients_ = true;
        fill_measures_ = true;
        fill_w_measures_ = true;
    }

    if (contains(flags, NewValueFlags::face_normal))
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
BasisElemValueFlagsHandler(const NewValueFlags &flags)
{
    if (contains(flags, NewValueFlags::value))
    {
        fill_values_ = true;
    }

    if (contains(flags, NewValueFlags::gradient))
    {
        fill_gradients_ = true;
    }

    if (contains(flags, NewValueFlags::hessian))
    {
        fill_hessians_ = true;
    }

    if (contains(flags, NewValueFlags::divergence))
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
BasisFaceValueFlagsHandler(const NewValueFlags &flags)
{
    if (contains(flags, NewValueFlags::face_value))
    {
        fill_values_ = true;
    }

    if (contains(flags, NewValueFlags::face_gradient))
    {
        fill_gradients_ = true;
    }

    if (contains(flags, NewValueFlags::face_hessian))
    {
        fill_hessians_ = true;
    }

    if (contains(flags, NewValueFlags::face_divergence))
    {
        fill_gradients_ = true;
        fill_divergences_ = true;
    }
}
//====================================================

#endif

IGA_NAMESPACE_CLOSE




