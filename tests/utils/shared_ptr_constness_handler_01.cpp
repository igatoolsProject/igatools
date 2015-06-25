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

/*
 * Test for SafeSTLVector, SafeSTLArray and their serialization
 * author: pauletti, martinelli
 * date:   2014-08-26
 * date:   2015-05-05 (added the tests for serialization)
 *
 */

#include "../tests.h"

#include <igatools/geometry/cartesian_grid.h>

IGA_NAMESPACE_OPEN

template<class T>
class SharedPtrConstnessHandler
{
public:
    using Ptr = const std::shared_ptr<T>;
    using PtrToConst = const std::shared_ptr<const T>;


    SharedPtrConstnessHandler(const Ptr &data)
        :
        data_(data),
        data_is_const_(false)
    {};

    SharedPtrConstnessHandler(const PtrToConst &data)
        :
        data_(data),
        data_is_const_(true)
    {};

    void print_info(LogStream &out) const
    {
        if (data_is_const_)
        {
            out << "Pointer to const data" << std::endl;
            boost::get<PtrToConst>(data_)->print_info(out);
        }
        else
        {
            out << "Pointer to non-const data" << std::endl;
            boost::get<Ptr>(data_)->print_info(out);
        }
    }

private:
    boost::variant<Ptr,PtrToConst> data_;

    bool data_is_const_;
};


IGA_NAMESPACE_CLOSE



using namespace iga;

using std::shared_ptr;

using Grid = CartesianGrid<1>;


void do_test_nonconst()
{
    OUTSTART

    shared_ptr<Grid> grid_nonconst = Grid::create();

    SharedPtrConstnessHandler<Grid> const_handler(grid_nonconst);
    const_handler.print_info(out);

    OUTEND
}

int main()
{
    OUTSTART

    do_test_nonconst();

    OUTEND

    return 0;

}
