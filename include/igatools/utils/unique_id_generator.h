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

// QualityAssurance: martinelli, 04 Sep 2014

#ifndef UNIQUE_ID_GENERATOR_H_
#define UNIQUE_ID_GENERATOR_H_

#include <igatools/base/config.h>




IGA_NAMESPACE_OPEN

/**
 * @brief The purpose of this class is to give an unique id to object that ask for it.
 *
 * Is basically a class with only static methods and data, in which the unique id is generated
 * by a static integer that increases its value each time the get_unique_id() method is invoked.
 *
 * @author M.Martinelli
 * @date 04 September 2014
 */
class UniqueIdGenerator
{
public:

    /** @name Constructors */
    ///@{
    ///@}

    /** @name Assignment operators */
    ///@{

    /**
     * Copy assignment operator
     */
    UniqueIdGenerator &operator=(const UniqueIdGenerator &obj) = delete;

    /**
     * Move assignment operator
     */
    UniqueIdGenerator &operator=(UniqueIdGenerator &&obj) = delete;
    ///@}


    /**
     * Returns a unique id each time this function is invoked.
     */
    static Index get_unique_id();

private:

    static Index id_;
};

IGA_NAMESPACE_CLOSE


#endif // UNIQUE_ID_GENERATOR_H_



