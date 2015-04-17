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


class GridFlags
{
public:
    static const ValueFlags valid_flags =
        ValueFlags::point|
        ValueFlags::w_measure;

    /** @name Constructors */
    ///@{
    /** Default constructor. Sets all boolean flags to false. */
    GridFlags() = default;

    /**
     * Constructor. Transforms the value flags for grid-like element accessor in
     * the correspondent booleans
     * that specify the quantities that must be computed/filled.
     */
    GridFlags(const ValueFlags &flags);

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

    /** Returns true if the nothing must be filled. */
    bool fill_none() const;

    /** Returns true if the quadrature points on the element must be filled. */
    bool fill_points() const;

    /** Returns true if the points are filled. */
    bool points_filled() const;

    /** Sets the filled status for points. */
    void set_points_filled(const bool status);


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
    FlagStatus points_flags_;

    FlagStatus w_measures_flags_;
};



inline
auto
create_function_flags_type_and_status()
{
    return boost::fusion::make_map<_Point,_Value,_Gradient,_Hessian,_Divergence>(
               FlagStatus(), FlagStatus(), FlagStatus(), FlagStatus(), FlagStatus());
}


class FunctionFlags
{
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
    FunctionFlags() = default;

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

    /** Returns true if the quantity associated to @p ValueType must be filled. */
    template<class ValueType>
    bool fill() const
    {
//        return value_type_flags_.at(ValueType::id).fill_;
        return boost::fusion::at_key<ValueType>(flags_type_and_status_).fill_;
    }

    /** Returns true if the quantity associated to @p ValueType is filled. */
    template<class ValueType>
    bool filled() const
    {
//        return value_type_flags_.at(ValueType::id).filled_;
        return boost::fusion::at_key<ValueType>(flags_type_and_status_).filled_;
    }

    /** Sets the filled @p status the quantity associated to @p ValueType. */
    template<class ValueType>
    void set_filled(const bool status)
    {
//        value_type_flags_[ValueType::id].filled_ = status;
        boost::fusion::at_key<ValueType>(flags_type_and_status_).filled_ = status;
    }


    /** Returns true if the nothing must be filled. */
    bool fill_none() const;


    /**
     * Prints internal information about the ElementValuesCache.
     * Its main use is for testing and debugging.
     */
    void print_info(LogStream &out) const;

protected:

    /**
     * Map used to realize the association between the ValueType and the relative FlagStatus.
     */
    decltype(create_function_flags_type_and_status()) flags_type_and_status_ = create_function_flags_type_and_status();
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
