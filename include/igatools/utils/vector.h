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

#ifndef VECTOR_H_
#define VECTOR_H_

#include <igatools/base/config.h>
#include <igatools/base/print_info_utils.h>

#include <vector>

IGA_NAMESPACE_OPEN

/**
 * @brief iga version on std::vector
 */
template<class T>
class vector : public std::vector<T>
{
public :
    /** Inherith the constructors of the base class. */
    using std::vector<T>::vector;

    Size size() const noexcept
    {
        return std::vector<T>::size();
    }

    typename std::vector<T>::reference operator[](Size n)
    {
        Assert(n<size(), ExcIndexRange(n, 0, size()));
        return std::vector<T>::operator[](n);
    }

    typename std::vector<T>::const_reference operator[](Size n) const
    {
        Assert(n<size(), ExcIndexRange(n, 0, size()));
        return std::vector<T>::operator[](n);
    }
    /**
     * @name Printing info
     */
    ///@{
    /**
     * Prints the content of the vector on the LogStream @p out.
     * Its use is intended mainly for testing and debugging purpose.
     */
private:
    template <class A>
    EnableIf<has_print_info<A>(0), void>
    t_print_info(LogStream &out) const
    {
        for (auto &entry : *this)
        {
            entry.print_info(out);
            out << std::endl;
        }
    }

    template <class A>
    EnableIf<(!has_print_info<A>(0)), void>
    t_print_info(LogStream &out) const
    {
        out << "[ ";
        for (auto &entry : *this)
        {
            out << entry << " ";
        }
        out << "]";
    }

public:
    void print_info(LogStream &out) const
    {
        t_print_info<T>(out);
    }
    ///@}

} ;


IGA_NAMESPACE_CLOSE


#endif /* VALUE_VECTOR_H_ */
