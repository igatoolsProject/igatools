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

#ifndef ARRAY_H_
#define ARRAY_H_

#include <igatools/base/config.h>
#include <igatools/utils/safe_stl_container.h>

#include <array>

IGA_NAMESPACE_OPEN
#if 0
/**
 * @brief iga version on std::array. It can be used as a std::array but in Debug mode
 * it provides bounds checking.
 */
template<class T, Size N>
class SafeSTLArray
{
private:
    using base_t = std::array<T, N>;
    base_t data_;
public:
    using size_type = Size;
    using reference = typename std::array<T, N>::reference;
    using const_reference = typename std::array<T, N>::const_reference;
    using iterator = typename std::array<T, N>::iterator;
    using const_iterator = typename std::array<T, N>::const_iterator;
    using value_type = typename std::array<T, N>::value_type;

    SafeSTLArray() = default;

    SafeSTLArray(const T &val)
    {
        data_.fill(val);
    }

    /**
     * Constructor using a std::array
     */
    SafeSTLArray(const base_t &data)
        :
        data_(data)
    {};

    /**
     * Initializer-list constructor.
     */
    constexpr SafeSTLArray(const std::initializer_list<T> &list)
        :
        data_(list)
    {};


    iterator begin() noexcept
    {
        return data_.begin();
    }
    const_iterator begin() const noexcept
    {
        return data_.begin();
    }
    iterator end() noexcept
    {
        return data_.end();
    }
    const_iterator end() const noexcept
    {
        return data_.end();
    }
    value_type *data() noexcept
    {
        return data_.data();
    }
    const value_type *data() const noexcept
    {
        return data_.data();
    }

    constexpr size_type size() const noexcept
    {
        return data_.size();
    }

    reference operator[](Size n)
    {
        Assert(n >= 0 && n <size(), ExcIndexRange(n, 0, size()));
        return data_[n];
    }

    const_reference operator[](Size n) const
    {
        Assert(n >= 0 && n < size(), ExcIndexRange(n, 0, size()));
        return data_[n];
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
        for (auto &entry : data_)
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
        for (auto &entry : data_)
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
#endif



/**
 * @brief iga version on std::array.
 * It can be used as a std::array but in Debug mode
 * it provides bounds checking.
 */
template<class T,int N>
class SafeSTLArray :
		public SafeSTLContainer<std::array<T,N>>
{
public :
    /** Inherit the constructors of the base class. */
    using SafeSTLContainer<std::array<T,N>>::SafeSTLContainer;

    /**
     * Default Constructor. It fills the <tt>N</tt> entries of the array with <tt>N</tt> copies
     * of the default-construced object of type <tt>T</tt>.
     */
    SafeSTLArray() {};


    /**
     * Constructor. It fills the <tt>N</tt> entries of the array with <tt>N</tt> copies
     * of the input argument <tt>val</tt>.
     */
    SafeSTLArray(const T &val)
    {
    	std::array<T,N>::fill(val);
    }

} ;

IGA_NAMESPACE_CLOSE


#endif
