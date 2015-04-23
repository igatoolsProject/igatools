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
 * @brief This structure describes the possible two states of a cache associated to a given ValueType.
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
auto
create_grid_flags_data()
{
    return boost::fusion::make_map<_Point,_W_Measure>(
               FlagStatus(), FlagStatus());
}

inline
auto
create_function_flags_data()
{
    return boost::fusion::make_map<_Point,_Value,_Gradient,_Hessian,_Divergence>(
               FlagStatus(), FlagStatus(), FlagStatus(), FlagStatus(), FlagStatus());
}

inline
auto
create_map_flags_data()
{
    return boost::fusion::make_map<
           _Point,
           _Value,_Gradient,_Hessian,
           _Measure, _W_Measure,
           _InvGradient, _InvHessian>(
               FlagStatus(), FlagStatus(), FlagStatus(), FlagStatus(),
               FlagStatus(), FlagStatus(), FlagStatus(), FlagStatus());
}

inline
ValueFlags
function_to_grid_flags(const ValueFlags &function_flags)
{
    ValueFlags transfer_flags = ValueFlags::w_measure |
                                ValueFlags::boundary_normal;
    ValueFlags g_flags = function_flags & transfer_flags;
    if (contains(function_flags, ValueFlags::point) || contains(function_flags, ValueFlags::value))
    {
        g_flags |= ValueFlags::point;
    }
    return g_flags;

}


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


/**
 * @brief Base class for defining an association between ValueType(s) and FlagStatus.
 *
 * Basically it is an associative container (a boost::fusion::map) in which the <tt>key</tt> is a ValueType
 * and the associated <tt>value</tt> is a FlagStatus.
 *
 * The type defining the map is passed as template argument <tt>FusionMap_ValueType_FlagStatus</tt>.
 */
template<class FusionMap_ValueType_FlagStatus>
class Flags
{
protected:
    /**
     * Constructor. Sets all boolean flags to false.
     * @note This constructor is intended to be called from a derived class only.
     */
    Flags(const FusionMap_ValueType_FlagStatus &map_value_types_and_flag_status)
        :
        map_value_types_and_flag_status_ {map_value_types_and_flag_status}
    {};


public:
    /**
     * @name Functions used to query or modify the Flag status for a given ValueType
     */
    ///@{
    /** Returns true if the quantity associated to @p ValueType must be filled. */
    template<class ValueType>
    bool fill() const
    {
        return boost::fusion::at_key<ValueType>(map_value_types_and_flag_status_).fill_;
    }

    /** Returns true if the quantity associated to @p ValueType is filled. */
    template<class ValueType>
    bool filled() const
    {
        return boost::fusion::at_key<ValueType>(map_value_types_and_flag_status_).filled_;
    }

    /** Sets the filled @p status the quantity associated to @p ValueType. */
    template<class ValueType>
    void set_filled(const bool status)
    {
        boost::fusion::at_key<ValueType>(map_value_types_and_flag_status_).filled_ = status;
    }
    ///@}

#if 0
    /**
     * Prints internal information about the FlagStatus for all valid ValueType(s).
     * Its main use is for testing and debugging.
     */
    void print_info(LogStream &out) const
    {
        boost::fusion::for_each(map_value_types_and_flag_status_,
                                [&out](const auto & type_and_status) -> void
        {
            using ValueType_Status = typename std::remove_reference<decltype(type_and_status)>::type;
            using ValueType = typename ValueType_Status::first_type;

            out.begin_item(ValueType::name);
            type_and_status.second.print_info(out);
            out.end_item();
        } // end lambda function
                               );
    }

    /** Returns true if the nothing must be filled. */
    bool fill_none() const
    {
        const bool fill_someone = boost::fusion::any(map_value_types_and_flag_status_,
                                                     [](const auto & type_and_status) -> bool
        {
            return type_and_status.second.fill_ == true;
        } // end lambda function
                                                    );

        return !fill_someone;
    }
#endif

protected:
    /**
     * Map used to realize the association between the ValueType and the relative FlagStatus.
     */
    FusionMap_ValueType_FlagStatus map_value_types_and_flag_status_;

    /**
     * Sets the fill status to TRUE for the types corresponing to the input @p flags.
     */
    void set_fill_status_true_from_value_flags(const ValueFlags &flags)
    {
        boost::fusion::for_each(map_value_types_and_flag_status_,
                                [&](auto & type_and_status) -> void
        {
            using ValueType_Status = typename std::remove_reference<decltype(type_and_status)>::type;
            using ValueType = typename ValueType_Status::first_type;

            if (contains(flags, ValueType::flag))
                type_and_status.second.fill_ = true;
        } // end lambda function
                               );
    }

public:

    /**
     * Returns the flags that are valid to be used with this class.
     */
    ValueFlags get_valid_flags() const
    {
        ValueFlags valid_flags = ValueFlags::none;

        boost::fusion::for_each(map_value_types_and_flag_status_,
                                [&](const auto & type_and_status) -> void
        {
            using ValueType_Status = typename std::remove_reference<decltype(type_and_status)>::type;
            using ValueType = typename ValueType_Status::first_type;

            valid_flags |= ValueType::flag;
        } // end lambda function
                               );
        return valid_flags;
    }
};


class FunctionFlags : public Flags<decltype(create_function_flags_data())>
{
    using parent_t = Flags<decltype(create_function_flags_data())>;
public:

//    static ValueFlags to_grid_flags(const ValueFlags &flag);

    /** @name Constructors */
    ///@{
    /**
     * Default constructor. Sets all boolean flags to false.
     */
    FunctionFlags()
        :
        parent_t(create_function_flags_data())
    {};

    /**
     * Constructor. It sets the fill status of the valid ValueType(s) given the input @p flags.
     */
    FunctionFlags(const ValueFlags &flag);

    /** Copy constructor. */
    FunctionFlags(const FunctionFlags &in) = default;

    /** Move constructor. */
    FunctionFlags(FunctionFlags &&in) = default;


    /** Destructor. */
    ~FunctionFlags() = default;
    ///@}


    /** @name Assignment operators */
    ///@{
    /** Copy assignment operator. */
    FunctionFlags &operator=(const FunctionFlags &in) = default;


    /** Move assignment operator. */
    FunctionFlags &operator=(FunctionFlags &&in) = default;
    ///@}
};


#if 0
class GridFlags : public Flags<decltype(create_grid_flags_data())>
{
    using parent_t = Flags<decltype(create_grid_flags_data())>;
public:

    /** @name Constructors */
    ///@{
    /**
     * Default constructor. Sets all boolean flags to false.
     */
    GridFlags()
        :
        parent_t(create_grid_flags_data())
    {};

    /**
     * Constructor. It sets the fill status of the valid ValueType(s) given the input @p flags.
     */
    GridFlags(const ValueFlags &flag);

    /** Copy constructor. */
    GridFlags(const GridFlags &in) = default;

    /** Move constructor. */
    GridFlags(GridFlags &&in) = default;


    /** Destructor. */
    ~GridFlags() = default;
    ///@}


    /** @name Assignment operators */
    ///@{
    /** Copy assignment operator. */
    GridFlags &operator=(const GridFlags &in) = default;


    /** Move assignment operator. */
    GridFlags &operator=(GridFlags &&in) = default;
    ///@}
};
#endif

class MappingFlags : public Flags<decltype(create_map_flags_data())>
{
    using parent_t = Flags<decltype(create_map_flags_data())>;
public:
#if 0
    static const ValueFlags valid_flags =
        FunctionFlags::valid_flags |
        ValueFlags::inv_gradient|
        ValueFlags::inv_hessian |
        ValueFlags::measure|
        ValueFlags::w_measure|
        ValueFlags::boundary_normal|
        ValueFlags::outer_normal|
        ValueFlags::curvature;
#endif

//    static ValueFlags to_function_flags(const ValueFlags &flag);

    /** @name Constructors */
    ///@{
    /** Default constructor. Sets all boolean flags to false. */
    MappingFlags()
        :
        parent_t(create_map_flags_data())
    {};

    /**
     * Constructor. Transforms the value flags for the mapping in the correspondent booleans
     * that specify the quantities that must be computed/filled.
     */
    MappingFlags(const ValueFlags &flags);

    /** Copy constructor. */
    MappingFlags(const MappingFlags &in) = default;

    /** Move constructor. */
    MappingFlags(MappingFlags &&in) = default;


    /** Destructor. */
    ~MappingFlags() = default;
    ///@}


    /** @name Assignment operators */
    ///@{
    /** Copy assignment operator. */
    MappingFlags &operator=(const MappingFlags &in) = default;


    /** Move assignment operator. */
    MappingFlags &operator=(MappingFlags &&in) = default;
    ///@}
};


IGA_NAMESPACE_CLOSE

#endif
