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

// QualityAssurance: martinelli, 21 Jan 2014


#include <igatools/base/cache_status.h>
#include <igatools/base/exceptions.h>


IGA_NAMESPACE_OPEN

CacheStatus::
CacheStatus(const CacheStatus &cache)
    :
    initialized_counter_ {cache.initialized_counter_},
    filled_counter_ {cache.filled_counter_},
    copied_ {true}
{}



CacheStatus &
CacheStatus::
operator=(const CacheStatus &cache)
{
    if (this != &cache)
    {
        initialized_counter_ = cache.initialized_counter_;
        filled_counter_ = cache.filled_counter_;
        copied_ = true;
    }
    return *this;
}



bool
CacheStatus::
is_initialized() const
{
    return (initialized_counter_ > 0) ? true : false;
}



bool
CacheStatus::
is_initialized_once() const
{
    return (initialized_counter_ == 1) ? true : false;
}



bool
CacheStatus::
is_filled() const
{
    return (filled_counter_ > 0) ? true : false;
}



bool
CacheStatus::
is_filled_once() const
{
    return (filled_counter_ == 1) ? true : false;
}



bool
CacheStatus::
is_copied() const
{
    return copied_;
}



void
CacheStatus::
set_initialized(const bool status)
{
    if (status == true)
    {
        // Setting the status to initialized/
        initialized_counter_++;
    }
    else
    {
        // Setting the status to not initialized.
        initialized_counter_ = 0;
    }

    // In any case the cache is not filled and therefore not ready to be used.
    filled_counter_ = 0;
}



void
CacheStatus::
set_filled(const bool status)
{
    if (status == true)
    {
        Assert(this->is_initialized(),ExcNotInitialized());

        // setting the status to filled
        filled_counter_++;
    }
    else
    {
        // setting the status to not filled
        filled_counter_ = 0;
    }
}



IGA_NAMESPACE_CLOSE
