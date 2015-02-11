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
 * Test for ValueTable and associated iterator and views
 * martinelli
 * 25 Feb 2014
 *
 */

#include "../tests.h"

#include <igatools/base/tensor.h>
#include <igatools/base/function.h>
#include <igatools/utils/value_table.h>

void do_test(Size n_funcs, Size n_pts)
{
    out << "========== BEGIN do_test (n_funcs=" << n_funcs
        << ",n_pts=" << n_pts <<") ==========" << endl ;

    using Grad = typename Function<1,1,1>::Gradient;
    using ClassToTest = ValueTable<Grad>;

    out << "Default constructor" << endl;
    ClassToTest table1;
    table1.print_info(out);
    out << endl;


    out << "Standard constructor (n_funcs=" << n_funcs << ",n_pts=" << n_pts <<")" << endl;
    ClassToTest table2(n_funcs,n_pts);
    table2.print_info(out);
    out << endl;

    Real val = 0.0;
    out << "Filling the table with the iterator" << endl;
    for (auto &v : table2)
    {
        v = val;
        val += 1.0;
    }
    table2.print_info(out);

    out << endl;
    out << "Testing the function view using the begin()/end() functions" << endl;
    for (int i_fn = 0 ; i_fn < n_funcs ; ++i_fn)
    {
        auto func_view = table2.get_function_view(i_fn);
        auto func_view_beg = func_view.begin();
        const auto func_view_end = func_view.end();

        out << "\tFunction["<< i_fn << "] = ";
        for (; func_view_beg != func_view_end ; ++func_view_beg)
            out << *func_view_beg << " ";
        out << endl;
    }
    out << endl;

    out << "Testing the function view using the access operator[]" << endl;
    for (int i_fn = 0 ; i_fn < n_funcs ; ++i_fn)
    {
        auto func_view = table2.get_function_view(i_fn);
        out << "\tFunction["<< i_fn << "] = ";
        for (int i_pt = 0 ; i_pt < n_pts ; ++i_pt)
            out << func_view[i_pt] << " ";
        out << endl;
    }
    out << endl;

    out << "Testing the point view using the begin()/end() functions" << endl;
    for (int i_pt = 0 ; i_pt < n_pts ; ++i_pt)
    {
        auto pt_view = table2.get_point_view(i_pt);
        auto pt_view_beg = pt_view.begin();
        const auto pt_view_end = pt_view.end();

        out << "\tPoint["<< i_pt << "] = ";
        for (; pt_view_beg != pt_view_end ; ++pt_view_beg)
            out << *pt_view_beg << " ";
        out << endl;
    }
    out << endl;

    out << "Testing the point view using the access operator[]" << endl;
    for (int i_pt = 0 ; i_pt < n_pts ; ++i_pt)
    {
        auto pt_view = table2.get_point_view(i_pt);

        out << "\tPoint["<< i_pt << "] = ";
        for (int i_fn = 0 ; i_fn < n_funcs ; ++i_fn)
            out << pt_view[i_fn] << " ";
        out << endl;
    }

    out << "========== END do_test (n_funcs=" << n_funcs
        << ",n_pts=" << n_pts <<") ==========" << endl ;
    out << endl;

}



int main()
{
    do_test(3,2);
    do_test(3,4);

    do_test(6,3);

    return 0;
}
