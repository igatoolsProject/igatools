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


#ifndef DENSE_VECTOR_H_
#define DENSE_VECTOR_H_

#include <igatools/base/config.h>
#include <igatools/base/logstream.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

IGA_NAMESPACE_OPEN

/**
 * @brief Dense vector with real entries.
 *
 * A typical use of it is the local vector.
 * @note This class inherits (also the constructors) from the vector type in the
 * <a href="http://www.boost.org/doc/libs/1_55_0/libs/numeric/ublas/doc/index.htm">Boost uBlas library</a>
 * Therefore it can be used as a boost::numeric::ublas::vector<Real>.
 *
 * @ingroup linear_algebra
 * @author M. Martinelli, 2012, 2013, 2014
 * @author S. Pauletti, 2012, 2013
 */
class DenseVector : public boost::numeric::ublas::vector<Real>
{
public:

    /** Inherith the constructors of the base class. */
    using vector<Real>::vector;

    /**
     * Assignment operator for assigning zeros to all entries of the vector
     * by writing
     * @code
       DenseVector vec;
       ... // working with the vector
       vector = 0.0; // reset the vector entries to zero
       @endcode
     * @note If used in Debug mode with a @p value different from zero,
     * an assertion will be raised.
     */
    DenseVector &operator=(const Real value);

    void print_info(LogStream &out) const
    {
        Assert(false,ExcNotImplemented());
        // out << SafeSTLVector<Real>::SafeSTLVector(*this);
    }

    /** Returns the number of entries in the DenseVector. */
    int size() const;
};


IGA_NAMESPACE_CLOSE

#endif /* DENSE_VECTOR_H_ */
