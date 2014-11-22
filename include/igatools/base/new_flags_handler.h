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



#ifndef NEW_FLAGS_HANDLER_H_
#define NEW_FLAGS_HANDLER_H_

#include <igatools/base/config.h>
#include <igatools/base/logstream.h>

IGA_NAMESPACE_OPEN


class GridFlags
{
public:
    static const NewValueFlags valid_flags =
        NewValueFlags::point|
        NewValueFlags::measure |
        NewValueFlags::w_measure |
        NewValueFlags::length;

    /** @name Constructors */
    ///@{
    /** Default constructor. Sets all boolean flags to false. */
    GridFlags() = default;

    /**
     * Constructor. Transforms the value flags for grid-like element accessor in
     * the correspondent booleans
     * that specify the quantities that must be computed/filled.
     */
    GridFlags(const NewValueFlags &flags);

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

    bool fill_lengths() const;
    bool lengths_filled() const;
    void set_lengths_filled(const bool status);

    /**
     * Prints internal information about the ElementValuesCache.
     * Its main use is for testing and debugging.
     */
    void print_info(LogStream &out) const;


protected:
    bool fill_points_ = false;

    bool points_filled_ = false;

    bool fill_measures_ = false;

    bool measures_filled_ = false;

    bool fill_w_measures_ = false;

    bool w_measures_filled_ = false;

    bool fill_lengths_   = false;

    bool lengths_filled_ = false;
};



class FunctionFlags
{
public:
    static const NewValueFlags valid_flags =
        NewValueFlags::point|
        NewValueFlags::value|
        NewValueFlags::gradient|
        NewValueFlags::hessian;

    static NewValueFlags to_grid_flags(const NewValueFlags &flag);

    /** @name Constructors */
    ///@{
    /**
     * Default constructor. Sets all boolean flags to false.
     */
    FunctionFlags() = default;

    FunctionFlags(const NewValueFlags &flag);

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

    /** Returns true if the nothing must be filled. */
    bool fill_none() const;

    /** Returns true if the values must be filled. */
    bool fill_points() const;

    /** Returns true if the values are filled. */
    bool points_filled() const;

    /** Sets the filled status for values. */
    void set_points_filled(const bool status);

    /** Returns true if the values must be filled. */
    bool fill_values() const;

    /** Returns true if the values are filled. */
    bool values_filled() const;

    /** Sets the filled status for values. */
    void set_values_filled(const bool status);

    /** Returns true if the gradients must be filled. */
    bool fill_gradients() const;

    /** Returns true if the gradients are filled. */
    bool gradients_filled() const;

    /** Sets the filled status for gradients. */
    void set_gradients_filled(const bool status);

    /** Returns true if the hessians must be filled. */
    bool fill_hessians() const;

    /** Returns true if the hessians are filled. */
    bool hessians_filled() const;

    /** Sets the filled status for hessians. */
    void set_hessians_filled(const bool status);


    /**
     * Prints internal information about the ElementValuesCache.
     * Its main use is for testing and debugging.
     */
    void print_info(LogStream &out) const;

protected:
    bool fill_points_ = false;
    bool points_filled_ = false;

    bool fill_values_ = false;
    bool values_filled_ = false;

    bool fill_gradients_ = false;
    bool gradients_filled_ = false;

    bool fill_hessians_ = false;
    bool hessians_filled_ = false;
};




class MappingFlags :
    public FunctionFlags
{
public:

    static const NewValueFlags valid_flags =
        FunctionFlags::valid_flags |
        NewValueFlags::inv_gradient|
        NewValueFlags::inv_hessian |
        NewValueFlags::measure|
        NewValueFlags::w_measure|
        NewValueFlags::boundary_normal|
        NewValueFlags::outer_normal|
        NewValueFlags::curvature;

    static NewValueFlags to_function_flags(const NewValueFlags &flag);

    /** @name Constructors */
    ///@{
    /** Default constructor. Sets all boolean flags to false. */
    MappingFlags() = default;

    /**
     * Constructor. Transforms the value flags for the mapping in the correspondent booleans
     * that specify the quantities that must be computed/filled.
     */
    MappingFlags(const NewValueFlags &flags);

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
    bool fill_inv_gradients_ = false;

    bool inv_gradients_filled_ = false;

    bool fill_inv_hessians_ = false;

    bool inv_hessians_filled_ = false;

    bool fill_measures_ = false;

    bool measures_filled_ = false;

    bool fill_w_measures_ = false;

    bool w_measures_filled_ = false;

};

#if 0

class PushFowardFlags
{
public:
    /** @name Constructors */
    ///@{
    /** Default constructor. Sets all boolean flags to false. */
    PushFowardFlags();

    /**
     * Constructor. Transforms the value flags for the mapping in the correspondent booleans
     * that specify the quantities that must be computed/filled.
     */
    PushFowardFlags(const NewValueFlags &flags);

    /** Copy constructor. */
    PushFowardFlags(const PushFowardFlags &in) = default;

    /** Move constructor. */
    PushFowardFlags(PushFowardFlags &&in) = default;


    /** Destructor. */
    ~PushFowardFlags() = default;
    ///@}


    NewValueFlags to_mapping_flags(const NewValueFlags &flag);

    /** @name Assignment operators */
    ///@{
    /** Copy assignment operator. */
    MappingFlags &operator=(const PushFowardFlags &in) = default;


    /** Move assignment operator. */
    PushFowardFlags &operator=(PushFowardFlags &&in) = default;
    ///@}

    /** Returns true if the nothing must be filled. */
    bool fill_none() const;

    /** Returns true if the gradients inverse must be filled. */
    bool fill_trans() const;

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


    /**
     * Prints internal information about the ElementValuesCache.
     * Its main use is for testing and debugging.
     */
    void print_info(LogStream &out) const;

protected:
    bool fill_inv_gradients_ = false;

    bool inv_gradients_filled_ = false;

    bool fill_inv_hessians_ = false;

    bool inv_hessians_filled_ = false;
};


class SpaceFlags
{
public:
    static const NewValueFlags valid_flags =
        NewValueFlags::value|
        NewValueFlags::gradient|
        NewValueFlags::hessian |
        NewValueFlags::w_measure;


    /** @name Constructors */
    ///@{
    /**
     * Default constructor. Sets all boolean flags to false.
     */
    SpaceFlags();

    SpaceFlags(const NewValueFlags &flag);

    /** Copy constructor. */
    SpaceFlags(const SpaceFlags &in) = default;

    /** Move constructor. */
    SpaceFlags(SpaceFlags &&in) = default;


    /** Destructor. */
    ~SpaceFlags() = default;
    ///@}

    NewValueFlags to_grid_flags(const NewValueFlags &flag);
    NewValueFlags to_ref_space_flag(const NewValueFlags &flag);
    NewValueFlags to_push_forward_flag(const NewValueFlags &flag);

    /** @name Assignment operators */
    ///@{
    /** Copy assignment operator. */
    SpaceFlags &operator=(const SpaceFlags &in) = default;


    /** Move assignment operator. */
    SpaceFlags &operator=(SpaceFlags &&in) = default;
    ///@}


    /** Returns true if the nothing must be filled. */
    bool fill_none() const;

    /** Returns true if the values must be filled. */
    bool fill_values() const;

    /** Returns true if the values are filled. */
    bool values_filled() const;

    /** Sets the filled status for values. */
    void set_values_filled(const bool status);

    /** Returns true if the gradients must be filled. */
    bool fill_gradients() const;

    /** Returns true if the gradients are filled. */
    bool gradients_filled() const;

    /** Sets the filled status for gradients. */
    void set_gradients_filled(const bool status);

    /** Returns true if the hessians must be filled. */
    bool fill_hessians() const;

    /** Returns true if the hessians are filled. */
    bool hessians_filled() const;

    /** Sets the filled status for hessians. */
    void set_hessians_filled(const bool status);


    /**
     * Prints internal information about the ElementValuesCache.
     * Its main use is for testing and debugging.
     */
    void print_info(LogStream &out) const;

protected:
    bool fill_points_ = false;
    bool points_filled_ = false;

    bool fill_values_ = false;
    bool values_filled_ = false;

    bool fill_gradients_ = false;
    bool gradients_filled_ = false;

    bool fill_hessians_ = false;
    bool hessians_filled_ = false;
};

#endif

IGA_NAMESPACE_CLOSE

#endif
