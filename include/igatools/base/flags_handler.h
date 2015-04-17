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



#ifndef NEW_FLAGS_HANDLER_H_
#define NEW_FLAGS_HANDLER_H_

#include <igatools/base/config.h>
#include <igatools/base/logstream.h>
#include <igatools/base/value_types.h>


#include <boost/fusion/include/make_map.hpp>
#include <boost/fusion/include/at_key.hpp>
#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/include/any.hpp>

IGA_NAMESPACE_OPEN

struct FlagStatus
{
    bool fill_ = false;
    bool filled_ = false;

    void print_info(LogStream &out) const
    {
        out << "   fill = " << fill_ << "    filled = " << filled_;
    }
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
     * @name Functions used to query or modify the Flag status for a give ValueType
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

protected:
    /**
     * Map used to realize the association between the ValueType and the relative FlagStatus.
     */
    FusionMap_ValueType_FlagStatus map_value_types_and_flag_status_;

    /**
     * Sets the fill status given the input @p flags.
     */
    void set_fill_status_from_value_flags(const ValueFlags &flags)
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
};


class FunctionFlags : public Flags<decltype(create_function_flags_data())>
{
    using parent_t = Flags<decltype(create_function_flags_data())>;
public:

    static const ValueFlags valid_flags =
        ValueFlags::point|
        ValueFlags::value|
        ValueFlags::gradient|
        ValueFlags::hessian|
        ValueFlags::divergence;

    static ValueFlags to_grid_flags(const ValueFlags &flag);

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



class GridFlags : public Flags<decltype(create_grid_flags_data())>
{
    using parent_t = Flags<decltype(create_grid_flags_data())>;
public:

    static const ValueFlags valid_flags =
        ValueFlags::point|
        ValueFlags::w_measure;

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


class MappingFlags :
    public FunctionFlags
{
public:

    static const ValueFlags valid_flags =
        FunctionFlags::valid_flags |
        ValueFlags::inv_gradient|
        ValueFlags::inv_hessian |
        ValueFlags::measure|
        ValueFlags::w_measure|
        ValueFlags::boundary_normal|
        ValueFlags::outer_normal|
        ValueFlags::curvature;

    static ValueFlags to_function_flags(const ValueFlags &flag);

    /** @name Constructors */
    ///@{
    /** Default constructor. Sets all boolean flags to false. */
    MappingFlags() = default;

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

    /** Returns true if the nothing must be filled. */
    bool fill_none() const;

    /** Returns true if the gradients inverse must be filled. */
    bool fill_inv_gradients() const;

    /** Returns true if the gradients are filled. */
    bool inv_gradients_filled() const;

    /** Sets the filled status for gradients. */
    void set_inv_gradients_filled(const bool status);

    /** Returns true if the hessians inverse must be filled. */
    bool fill_inv_hessians() const;

    /** Returns true if the hessians are filled. */
    bool inv_hessians_filled() const;

    /** Sets the filled status for hessians. */
    void set_inv_hessians_filled(const bool status);

    /** Returns true if the element measure must be filled. */
    bool fill_measures() const;

    /** Returns true if the measures are filled. */
    bool measures_filled() const;

    /** Sets the filled status for measures. */
    void set_measures_filled(const bool status);

    /** Returns true if the quadrature weight multiplied by the element measure must be filled. */
    bool fill_w_measures() const;

    /** Returns true if the w_measures are filled. */
    bool w_measures_filled() const;

    /** Sets the filled status for w_measures. */
    void set_w_measures_filled(const bool status);
    /**
     * Prints internal information about the ElementValuesCache.
     * Its main use is for testing and debugging.
     */
    void print_info(LogStream &out) const;

protected:
    FlagStatus inv_gradients_flags_;

    FlagStatus inv_hessians_flags_;

    FlagStatus measures_flags_;

    FlagStatus w_measures_flags_;
};


IGA_NAMESPACE_CLOSE

#endif
