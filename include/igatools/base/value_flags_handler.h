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



#ifndef VALUE_FLAGS_HANDLER_H_
#define VALUE_FLAGS_HANDLER_H_

#include <igatools/base/config.h>

IGA_NAMESPACE_OPEN




class ValueFlagsHandler
{
public:
    /** @name Constructors */
    ///@{
    /**
     * Default constructor. Sets all boolean flags to false
     * (except fill_none_ that is set to true).
     */
    ValueFlagsHandler();

    /** Copy constructor. */
    ValueFlagsHandler(const ValueFlagsHandler &in) = default;

    /** Move constructor. */
    ValueFlagsHandler(ValueFlagsHandler &&in) = default;


    /** Destructor. */
    ~ValueFlagsHandler() = default;
    ///@}


    /** @name Assignment operators */
    ///@{
    /** Copy assignment operator. */
    ValueFlagsHandler &operator=(const ValueFlagsHandler &in) = default;


    /** Move assignment operator. */
    ValueFlagsHandler &operator=(ValueFlagsHandler &&in) = default;
    ///@}


    /** Returns true if the values must be filled. */
    bool fill_values() const;


    /** Returns true if the gradients must be filled. */
    bool fill_gradients() const;


    /** Returns true if the hessians must be filled. */
    bool fill_hessians() const;


protected:
    bool fill_none_ = true;

    bool fill_values_ = false;

    bool fill_gradients_ = false;

    bool fill_hessians_ = false;
};



/**
 * This is an helper class that is intended to be used as a filter for the flags that
 * refers to a grid-like element accessor.
 *
 * The enum class ValueFlags is a bitmask that implements a lot of different flags,
 * also referring to different concepts, and is therefore difficult to manage.
 * This is why this class is useful: the unique constructor take as input argument a ValueFlags
 * entry and filters the values that have valid meaning for a grid-like element accessor,
 * setting the corresponding boolean entries properly.
 *
 * @author M. Martinelli
 * @date 14 Mar 2014
 */
class GridElemValueFlagsHandler
{
public:
    /** @name Constructors */
    ///@{
    /** Default constructor. Sets all boolean flags to false. */
    GridElemValueFlagsHandler();

    /**
     * Constructor. Transforms the value flags for grid-like element accessor in
     * the correspondent booleans
     * that specify the quantities that must be computed/filled.
     */
    GridElemValueFlagsHandler(const ValueFlags &flags);

    /** Copy constructor. */
    GridElemValueFlagsHandler(const GridElemValueFlagsHandler &in) = default;

    /** Move constructor. */
    GridElemValueFlagsHandler(GridElemValueFlagsHandler &&in) = default;


    /** Destructor. */
    ~GridElemValueFlagsHandler() = default;
    ///@}


    /** @name Assignment operators */
    ///@{
    /** Copy assignment operator. */
    GridElemValueFlagsHandler &operator=(const GridElemValueFlagsHandler &in) = default;


    /** Move assignment operator. */
    GridElemValueFlagsHandler &operator=(GridElemValueFlagsHandler &&in) = default;
    ///@}

    /** Returns true if the quadrature points on the element must be filled. */
    bool fill_points() const;


    /** Returns true if the element measure must be filled. */
    bool fill_measures() const;


    /** Returns true if the quadrature weight multiplied by the element measure must be filled. */
    bool fill_w_measures() const;


protected:
    bool fill_points_ = false;

    bool fill_measures_ = false;

    bool fill_w_measures_ = false;
};



/**
 * This is an helper class that is intended to be used as a filter for the flags that
 * refers to the mapping on the element.
 *
 * The enum class ValueFlags is a bitmask that implements a lot of different flags,
 * also referring to different concepts, and is therefore difficult to manage.
 * This is why this class is useful: the unique constructor take as input argument a ValueFlags
 * entry and filters the values that have valid meaning for the mapping on the element, setting the
 * corresponding boolean entries properly.
 *
 * @author M. Martinelli
 * @date 14 Mar 2014
 */
class MappingValueFlagsHandler :
    public ValueFlagsHandler,
    public GridElemValueFlagsHandler
{
public:
    /** @name Constructors */
    ///@{
    /** Default constructor. Sets all boolean flags to false. */
    MappingValueFlagsHandler();

    /**
     * Constructor. Transforms the value flags for the mapping in the correspondent booleans
     * that specify the quantities that must be computed/filled.
     */
    MappingValueFlagsHandler(const ValueFlags &flags);

    /** Copy constructor. */
    MappingValueFlagsHandler(const MappingValueFlagsHandler &in) = default;

    /** Move constructor. */
    MappingValueFlagsHandler(MappingValueFlagsHandler &&in) = default;


    /** Destructor. */
    ~MappingValueFlagsHandler() = default;
    ///@}


    /** @name Assignment operators */
    ///@{
    /** Copy assignment operator. */
    MappingValueFlagsHandler &operator=(const MappingValueFlagsHandler &in) = default;


    /** Move assignment operator. */
    MappingValueFlagsHandler &operator=(MappingValueFlagsHandler &&in) = default;
    ///@}
    /** Returns true if the gradients inverse must be filled. */
    bool fill_inv_gradients() const;


    /** Returns true if the hessians inverse must be filled. */
    bool fill_inv_hessians() const;

protected:
    bool fill_inv_gradients_ = false;

    bool fill_inv_hessians_ = false;
};



/**
 * This is an helper class that is intended to be used as a filter for the flags that
 * refers to the face of a mapping.
 *
 * The enum class ValueFlags is a bitmask that implements a lot of different flags,
 * also referring to different concepts, and is therefore difficult to manage.
 * This is why this class is useful: the unique constructor take as input argument a ValueFlags
 * entry and filters the values that have valid meaning for the face of a mapping, setting the
 * corresponding boolean entries properly.
 *
 * @author M. Martinelli
 * @date 14 Mar 2014
 */
class MappingFaceValueFlagsHandler :
    public MappingValueFlagsHandler
{
public:
    /** @name Constructors */
    ///@{
    /** Default constructor. Sets all boolean flags to false. */
    MappingFaceValueFlagsHandler();

    /**
     * Constructor. Transforms the value flags for the mapping face
     * in the correspondent booleans
     * that specify the quantities that must be computed/filled.
     */
    MappingFaceValueFlagsHandler(const ValueFlags &flags);

    /** Copy constructor. */
    MappingFaceValueFlagsHandler(const MappingFaceValueFlagsHandler &in) = default;

    /** Move constructor. */
    MappingFaceValueFlagsHandler(MappingFaceValueFlagsHandler &&in) = default;


    /** Destructor. */
    ~MappingFaceValueFlagsHandler() = default;
    ///@}


    /** @name Assignment operators */
    ///@{
    /** Copy assignment operator. */
    MappingFaceValueFlagsHandler &operator=(const MappingFaceValueFlagsHandler &in) = default;


    /** Move assignment operator. */
    MappingFaceValueFlagsHandler &operator=(MappingFaceValueFlagsHandler &&in) = default;
    ///@}


    /** Returns true if the gradients inverse must be filled. */
    bool fill_normals() const;

protected:
    bool fill_normals_ = false;
};


IGA_NAMESPACE_CLOSE



#endif // #ifndef VALUE_FLAGS_HANDLER_H_
