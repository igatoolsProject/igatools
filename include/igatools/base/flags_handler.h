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

#ifndef __FLAGS_HANDLER_H_
#define __FLAGS_HANDLER_H_

#include <igatools/base/config.h>
#include <igatools/base/logstream.h>
#include <igatools/base/value_types.h>


#include <boost/fusion/include/make_map.hpp>
#include <boost/fusion/include/at_key.hpp>
#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/include/any.hpp>

IGA_NAMESPACE_OPEN

/**
 * @brief This structure describes the possible two states (<tt>fill</tt> and <tt>filled</tt>)
 * of a cache associated to a given ValueType.
 */
struct FlagStatus
{
    bool fill_ = false;
    bool filled_ = false;

    void print_info(LogStream &out) const
    {
        out << "   fill = " << fill_ << "    filled = " << filled_;
    }
    /*
        bool fill() const
        {
            return fill_;
        };

        void set_fill(const bool fill_status)
        {
            fill_ = fill_status;
        };

        bool filled() const
        {
            return filled_;
        };

        void set_filled(const bool filled_status)
        {
            filled_ = filled_status;
        };
    //*/
};


inline
ValueFlags
mapping_to_function_flags(const ValueFlags &flags)
{
    ValueFlags valid_func_flags = ValueFlags::value |
                                  ValueFlags::gradient |
                                  ValueFlags::hessian |
                                  ValueFlags::divergence |
                                  ValueFlags::point;

    ValueFlags transfer_flags = ValueFlags::measure |
                                ValueFlags::w_measure |
                                ValueFlags::boundary_normal |
                                valid_func_flags;


    ValueFlags f_flags = flags & transfer_flags;

    if (contains(flags, ValueFlags::measure) ||
        contains(flags, ValueFlags::w_measure) ||
        contains(flags, ValueFlags::inv_gradient) ||
        contains(flags, ValueFlags::outer_normal))
        f_flags |=  ValueFlags::gradient;

    if (contains(flags, ValueFlags::inv_hessian) ||
        contains(flags, ValueFlags::curvature))
        f_flags |=  ValueFlags::gradient | ValueFlags::hessian;

    return f_flags;
}


IGA_NAMESPACE_CLOSE

#endif
