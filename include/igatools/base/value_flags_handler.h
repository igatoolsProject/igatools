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
     * Default constructor. Sets all boolean flags to false.
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


protected:
    bool fill_values_ = false;
    bool values_filled_ = false;

    bool fill_gradients_ = false;
    bool gradients_filled_ = false;

    bool fill_hessians_ = false;
    bool hessians_filled_ = false;
};



/**
 * @brief This is an helper class that is intended to be used as a filter for the flags that
 * refers to a grid-like element accessor.
 *
 * The enum class ValueFlags is a bitmask that implements a lot of different flags,
 * also referring to different concepts, and is therefore difficult to manage.
 * This is the reason that makes this class useful: the unique constructor
 * GridElemValueFlagsHandler(const ValueFlags &flags) takes as input argument a ValueFlags
 * entry and filters the values that have valid meaning for a grid-like element accessor,
 * setting the corresponding boolean entries properly.
 *
 * The ValueFlags filtered by this class are:
 * - ValueFlags::point
 * - ValueFlags::measure
 * - ValueFlags::w_measure
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


protected:
    bool fill_points_ = false;

    bool points_filled_ = false;

    bool fill_measures_ = false;

    bool measures_filled_ = false;

    bool fill_w_measures_ = false;

    bool w_measures_filled_ = false;
};



/**
 * @brief This is an helper class that is intended to be used as a filter for the flags that
 * refers to a grid-like element-face accessor.
 *
 * The enum class ValueFlags is a bitmask that implements a lot of different flags,
 * also referring to different concepts, and is therefore difficult to manage.
 * This is the reason that makes this class useful: the unique constructor
 * GridFaceValueFlagsHandler(const ValueFlags &flags) takes as input argument a ValueFlags
 * entry and filters the values that have valid meaning for a grid-like element accessor,
 * setting the corresponding boolean entries properly.
 *
 * The ValueFlags filtered by this class are:
 * - ValueFlags::face_point
 * - ValueFlags::face_measure
 * - ValueFlags::face_w_measure
 * - ValueFlags::face_normal
 *
 * @author M. Martinelli
 * @date 14 Mar 2014
 */
class GridFaceValueFlagsHandler : public GridElemValueFlagsHandler
{
public:
    /** @name Constructors */
    ///@{
    /** Default constructor. Sets all boolean flags to false. */
    GridFaceValueFlagsHandler();

    /**
     * Constructor. Transforms the value flags for grid-like element-face accessor in
     * the correspondent booleans
     * that specify the quantities that must be computed/filled.
     */
    GridFaceValueFlagsHandler(const ValueFlags &flags);

    /** Copy constructor. */
    GridFaceValueFlagsHandler(const GridFaceValueFlagsHandler &in) = default;

    /** Move constructor. */
    GridFaceValueFlagsHandler(GridFaceValueFlagsHandler &&in) = default;


    /** Destructor. */
    ~GridFaceValueFlagsHandler() = default;
    ///@}


    /** @name Assignment operators */
    ///@{
    /** Copy assignment operator. */
    GridFaceValueFlagsHandler &operator=(const GridFaceValueFlagsHandler &in) = default;


    /** Move assignment operator. */
    GridFaceValueFlagsHandler &operator=(GridFaceValueFlagsHandler &&in) = default;
    ///@}

    /** Returns true if the nothing must be filled. */
    bool fill_none() const;

    /** Returns true if the gradients inverse must be filled. */
    bool fill_normals() const;

    /** Returns true if the normals are filled. */
    bool normals_filled() const;

    /** Sets the filled status for normals. */
    void set_normals_filled(const bool status);


protected:
    bool fill_normals_ = false;

    bool normals_filled_ = false;
};


/**
 * @brief This is an helper class that is intended to be used as a filter for the flags that
 * refers to the mapping on the element.
 *
 * The enum class ValueFlags is a bitmask that implements a lot of different flags,
 * also referring to different concepts, and is therefore difficult to manage.
 * This is the reason that makes this class useful: the unique constructor
 * MappingElemValueFlagsHandler(const ValueFlags &flags) takes as input argument a ValueFlags
 * entry and filters the values that have valid meaning for the mapping on the element, setting the
 * corresponding boolean entries properly.
 *
 * The ValueFlags filtered by this class are:
 * - ValueFlags::point
 * - ValueFlags::map_value
 * - ValueFlags::map_gradient
 * - ValueFlags::map_hessian
 * - ValueFlags::map_inv_gradient
 * - ValueFlags::map_inv_hessian
 * - ValueFlags::measure
 * - ValueFlags::w_measure
 *
 * @author M. Martinelli
 * @date 14 Mar 2014
 */
class MappingElemValueFlagsHandler :
    public ValueFlagsHandler,
    public GridElemValueFlagsHandler
{
public:
    /** @name Constructors */
    ///@{
    /** Default constructor. Sets all boolean flags to false. */
    MappingElemValueFlagsHandler();

    /**
     * Constructor. Transforms the value flags for the mapping in the correspondent booleans
     * that specify the quantities that must be computed/filled.
     */
    MappingElemValueFlagsHandler(const ValueFlags &flags);

    /** Copy constructor. */
    MappingElemValueFlagsHandler(const MappingElemValueFlagsHandler &in) = default;

    /** Move constructor. */
    MappingElemValueFlagsHandler(MappingElemValueFlagsHandler &&in) = default;


    /** Destructor. */
    ~MappingElemValueFlagsHandler() = default;
    ///@}


    /** @name Assignment operators */
    ///@{
    /** Copy assignment operator. */
    MappingElemValueFlagsHandler &operator=(const MappingElemValueFlagsHandler &in) = default;


    /** Move assignment operator. */
    MappingElemValueFlagsHandler &operator=(MappingElemValueFlagsHandler &&in) = default;
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

protected:
    bool fill_inv_gradients_ = false;

    bool inv_gradients_filled_ = false;

    bool fill_inv_hessians_ = false;

    bool inv_hessians_filled_ = false;
};



/**
 * @brief This is an helper class that is intended to be used as a filter for the flags that
 * refers to the mapping on a face.
 *
 * The enum class ValueFlags is a bitmask that implements a lot of different flags,
 * also referring to different concepts, and is therefore difficult to manage.
 * This is the reason that makes this class useful: the unique constructor
 * MappingFaceValueFlagsHandler(const ValueFlags &flags) takes as input argument a ValueFlags
 * entry and filters the values that have valid meaning for the face of a mapping, setting the
 * corresponding boolean entries properly.
 *
 * The ValueFlags filtered by this class are:
 * - ValueFlags::face_point
 * - ValueFlags::map_face_value
 * - ValueFlags::map_face_gradient
 * - ValueFlags::map_face_hessian
 * - ValueFlags::map_face_inv_gradient
 * - ValueFlags::map_face_inv_hessian
 * - ValueFlags::face_measure
 * - ValueFlags::face_w_measure
 * - ValueFlags::face_normal
 *
 * @author M. Martinelli
 * @date 14 Mar 2014
 */
class MappingFaceValueFlagsHandler :
    public MappingElemValueFlagsHandler
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

    /** Returns true if the nothing must be filled. */
    bool fill_none() const;

    /** Returns true if the gradients inverse must be filled. */
    bool fill_normals() const;

    /** Returns true if the normals are filled. */
    bool normals_filled() const;

    /** Sets the filled status for normals. */
    void set_normals_filled(const bool status);


protected:
    bool fill_normals_ = false;

    bool normals_filled_ = false;
};






/**
 * @brief This is an helper class that is intended to be used as a filter for the flags that
 * refers to the basis function on the element.
 *
 * The enum class ValueFlags is a bitmask that implements a lot of different flags,
 * also referring to different concepts, and is therefore difficult to manage.
 * This is the reason that makes this class useful: the unique constructor
 * BasisElemValueFlagsHandler(const ValueFlags &flags) takes as input argument a ValueFlags
 * entry and filters the values that have valid meaning for the mapping on the element, setting the
 * corresponding boolean entries properly.
 *
 * The ValueFlags filtered by this class are:
 * - ValueFlags::value
 * - ValueFlags::gradient
 * - ValueFlags::hessian
 * - ValueFlags::divergence
 *
 * @author M. Martinelli
 * @date 17 Mar 2014
 */
class BasisElemValueFlagsHandler : public ValueFlagsHandler
{
public:
    /** @name Constructors */
    ///@{
    /**
     * Default constructor. Sets all boolean flags to false
     * (except fill_none_ that is set to true).
     */
    BasisElemValueFlagsHandler();

    /**
     * Constructor. Transforms the value flags for the basis functions on the element
     * in the correspondent booleans
     * that specify the quantities that must be computed/filled.
     */
    BasisElemValueFlagsHandler(const ValueFlags &flags);


    /** Copy constructor. */
    BasisElemValueFlagsHandler(const BasisElemValueFlagsHandler &in) = default;

    /** Move constructor. */
    BasisElemValueFlagsHandler(BasisElemValueFlagsHandler &&in) = default;


    /** Destructor. */
    ~BasisElemValueFlagsHandler() = default;
    ///@}


    /** @name Assignment operators */
    ///@{
    /** Copy assignment operator. */
    BasisElemValueFlagsHandler &operator=(const BasisElemValueFlagsHandler &in) = default;


    /** Move assignment operator. */
    BasisElemValueFlagsHandler &operator=(BasisElemValueFlagsHandler &&in) = default;
    ///@}

    /** Returns true if the nothing must be filled. */
    bool fill_none() const;

    /** Returns true if the divergences must be filled. */
    bool fill_divergences() const;

    /** Returns true if the divergences are filled. */
    bool divergences_filled() const;

    /** Sets the filled status for divergences. */
    void set_divergences_filled(const bool status);

protected:
    bool fill_divergences_ = false;

    bool divergences_filled_ = false;
};






/**
 * @brief This is an helper class that is intended to be used as a filter for the flags that
 * refers to the basis function on the face if an element.
 *
 * The enum class ValueFlags is a bitmask that implements a lot of different flags,
 * also referring to different concepts, and is therefore difficult to manage.
 * This is the reason that makes this class useful: the unique constructor
 * BasisFaceValueFlagsHandler(const ValueFlags &flags) takes as input argument a ValueFlags
 * entry and filters the values that have valid meaning for the mapping on the element, setting the
 * corresponding boolean entries properly.
 *
 * The ValueFlags filtered by this class are:
 * - ValueFlags::face_value
 * - ValueFlags::face_gradient
 * - ValueFlags::face_hessian
 * - ValueFlags::face_divergence
 *
 * @author M. Martinelli
 * @date 17 Mar 2014
 */
class BasisFaceValueFlagsHandler : public BasisElemValueFlagsHandler
{
public:
    /** @name Constructors */
    ///@{
    /**
     * Default constructor. Sets all boolean flags to false
     * (except fill_none_ that is set to true).
     */
    BasisFaceValueFlagsHandler();

    /**
     * Constructor. Transforms the value flags for the basis functions on the face
     * of the element in the correspondent booleans
     * that specify the quantities that must be computed/filled.
     */
    BasisFaceValueFlagsHandler(const ValueFlags &flags);


    /** Copy constructor. */
    BasisFaceValueFlagsHandler(const BasisFaceValueFlagsHandler &in) = default;

    /** Move constructor. */
    BasisFaceValueFlagsHandler(BasisFaceValueFlagsHandler &&in) = default;


    /** Destructor. */
    ~BasisFaceValueFlagsHandler() = default;
    ///@}


    /** @name Assignment operators */
    ///@{
    /** Copy assignment operator. */
    BasisFaceValueFlagsHandler &operator=(const BasisFaceValueFlagsHandler &in) = default;


    /** Move assignment operator. */
    BasisFaceValueFlagsHandler &operator=(BasisFaceValueFlagsHandler &&in) = default;
    ///@}
};



IGA_NAMESPACE_CLOSE



#endif // #ifndef VALUE_FLAGS_HANDLER_H_
