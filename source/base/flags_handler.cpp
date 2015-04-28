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


void
FlagStatus::
print_info(LogStream &out) const
{
    out << "   fill = " << fill_ << "    filled = " << filled_;
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

#if 0
GridFlags::
GridFlags(const ValueFlags &flags)
    :
    GridFlags()
{
    this->set_fill_status_true_from_value_flags(flags);
}
#endif

#if 0
FunctionFlags::
FunctionFlags(const ValueFlags &flags)
    :
    FunctionFlags()
{
    const auto valid_flags = this->get_valid_flags();
    auto f_flags = flags & valid_flags;
    if (contains(f_flags, ValueFlags::divergence))
        f_flags |= ValueFlags::gradient;

    this->set_fill_status_true_from_value_flags(f_flags);
}
#endif

#if 0
ValueFlags
FunctionFlags::to_grid_flags(const ValueFlags &flags)
{
    ValueFlags transfer_flag = ValueFlags::w_measure |
                               ValueFlags::boundary_normal;
    ValueFlags g_flag = flags & transfer_flag;
    if (contains(flags, ValueFlags::point) || contains(flags, ValueFlags::value))
    {
        g_flag |= ValueFlags::point;
    }
    return g_flag;
}
#endif


#if 0
MappingFlags::
MappingFlags(const ValueFlags &flags)
    :
    MappingFlags()
{
    const auto valid_flags = this->get_valid_flags();
    auto m_flags = flags & valid_flags;

    if (contains(flags, ValueFlags::boundary_normal) ||
        contains(flags, ValueFlags::curvature))
        m_flags |= ValueFlags::inv_gradient;

    if (contains(flags, ValueFlags::w_measure))
        m_flags |= ValueFlags::measure;

    this->set_fill_status_true_from_value_flags(m_flags);
}
#endif

#if 0
ValueFlags
MappingFlags::to_function_flags(const ValueFlags &flags)
{
    FunctionFlags func_flags;

    ValueFlags transfer_flag = ValueFlags::measure |
                               ValueFlags::w_measure |
                               ValueFlags::boundary_normal |
                               func_flags.get_valid_flags();


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
#endif



IGA_NAMESPACE_CLOSE




