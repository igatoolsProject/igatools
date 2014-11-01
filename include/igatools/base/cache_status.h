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

#ifndef CACHE_STATUS_H_
#define CACHE_STATUS_H_

#include <igatools/base/config.h>

#include<igatools/geometry/unit_element.h>
#include<igatools/base/quadrature.h>
#include<tuple>

IGA_NAMESPACE_OPEN

/**
 * @brief This class is used to manage and query the status of a cache.
 *
 * The possible status (not mutually exclusive) of a cache are:
 * - <b>initialized</b>, true if the cache is initialized at least one time;
 * - <b>initialized_once</b>, true if the cache is initialized exactly one time;
 * - <b>filled</b>, true if the cache is filled at least one time and ready to be used;
 * - <b>filled_once</b>, true if the cache is filled exactly one time and ready to be used;
 * - <b>copied</b>, true if the cache is copied from another object.
 *
 * @author M.Martinelli
 * @date 2013
 */
class CacheStatus
{
public:
    /** @name Constructor and destructor. */
    ///@{

    /** Default constructor. Sets all the status to false. */
    CacheStatus() = default;


    /** Copy constructor. */
    CacheStatus(const CacheStatus &cache);

    /** Move constructor. */
    CacheStatus(CacheStatus &&cache) = default;


    /** Destructor. */
    ~CacheStatus() = default;
    ///@}

    /** @name Assignment operator. */
    ///@{

    /** Copy assignment operator. */
    CacheStatus &operator=(const CacheStatus &cache);

    /** Move assignment operator. */
    CacheStatus &operator=(CacheStatus &&cache) = default;

    ///@}

    /** @name Functions for querying the status of the cache. */
    ///@{
    /** Return true if the cache is initialized at least one time. */
    bool is_initialized() const ;

    /** Return true if the cache is initialized exactly one time. */
    bool is_initialized_once() const ;

    /** Return true if the cache is filled at least one time and ready to be used. */
    bool is_filled() const ;

    /** Return true if the cache is filled exactly one time and ready to be used. */
    bool is_filled_once() const ;

    /** Return true if the cache is copied from another object. */
    bool is_copied() const ;
    ///@}


    /** @name Function for setting the status of the cache. */
    ///@{

    /**
     * This function is used to set the status of the cache to be "initialized" or not.
     * If <tt>status == true</tt> then the cache is considered to be initialized.
     * If <tt>status == false</tt> then the cache is considered to be not initialized.
     * @post After calling this function the status of the cache is set to "not filled".
     */
    void set_initialized(const bool status);

    /**
     * This function is used to set the status of the cache to be "filled" or not.
     * If <tt>status == true</tt> then the cache is considered to be filled.
     * If <tt>status == false</tt> then the cache is considered to be not filled.
     * @pre In order to call this function, the cache must be in the "initialized" status.
     * In debug mode a check on this precondition is performed.
     */
    void set_filled(const bool status);
    ///@}

private:

    /**
     * Counts the number of times that the cache is initialized.
     * Zero means that the cache is not initialized.
     */
    int initialized_counter_ = 0;

    /**
     * Counts the number of times that the cache is filled.
     * Zero means that the cache is not filled.
     */
    int filled_counter_ = 0;


    /** True if the cache is copied from another object. */
    bool copied_ = false;
};




IGA_NAMESPACE_CLOSE


#endif // CACHE_STATUS_H_
