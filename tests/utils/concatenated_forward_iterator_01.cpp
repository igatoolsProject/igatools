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
/*
 *
 * Test for the Concatenated iterator class, using as base iterator the
 * std::vector<int>::iterator.
 *
 * martinelli
 * 30 May 2014
 *
 */

#include "../tests.h"

#include <igatools/utils/concatenated_iterator.h>
#include <vector>

using std::vector;


using std::endl;


void do_test_1_const()
{
    out << "========== do_test_1_const() --- begin ==========" << endl;
    vector<int> v0 = {1,2,3,4};
    vector<int> v1 = {5,6,7,8,9};
    vector<int> v2 = {10,11,12};
    vector<int> v3 = {13};

    using VecIterator = vector<int>::const_iterator;
    using VecView = ConstView<VecIterator>;

    std::vector<VecView> ranges;
    ranges.push_back(VecView(v0.begin(),v0.end()));
    ranges.push_back(VecView(v1.begin(),v1.end()));
    ranges.push_back(VecView(v2.begin(),v2.end()));
    ranges.push_back(VecView(v3.begin(),v3.end()));


    ConcatenatedForwardConstIterator<VecView> begin(ranges,0);
    ConcatenatedForwardConstIterator<VecView> end(ranges,IteratorState::pass_the_end);



    using std::endl;
    int i = 0;
    for (; begin != end ; ++begin, ++i)
        out << "i = " << i << "     value = " << *begin << endl;

    out << "========== do_test_1_const() --- end ==========" << endl;
    out << endl;
}

void do_test_1_nonconst()
{
    out << "========== do_test_1_nonconst() --- begin ==========" << endl;
    vector<int> v0 = {1,2,3,4};
    vector<int> v1 = {5,6,7,8,9};
    vector<int> v2 = {10,11,12};
    vector<int> v3 = {13};

    using VecIterator = vector<int>::iterator;
    using VecConstIterator = vector<int>::const_iterator;
    using VecView = View<VecIterator,VecConstIterator>;

    std::vector<VecView> ranges;
    ranges.push_back(VecView(v0.begin(),v0.end()));
    ranges.push_back(VecView(v1.begin(),v1.end()));
    ranges.push_back(VecView(v2.begin(),v2.end()));
    ranges.push_back(VecView(v3.begin(),v3.end()));


    ConcatenatedForwardIterator<VecView> begin(ranges,0);
    ConcatenatedForwardIterator<VecView> end(ranges,IteratorState::pass_the_end);



    using std::endl;
    int i = 0;
    for (; begin != end ; ++begin, ++i)
        out << "i = " << i << "     value = " << *begin << endl;

    out << "========== do_test_1_nonconst() --- end ==========" << endl;
    out << endl;
}


void do_test_2()
{

    out << "========== do_test_2() --- begin ==========" << endl;

    vector<int> v0a = {1};
    vector<int> v1a = {2,3,4};
    vector<int> v2a = {5,6,7,8,9};

    vector<int> v0b = {10,11,12};
    vector<int> v1b = {13};

    using ItType_0 = vector<int>::iterator;
    using VecView0 = ConstView<ItType_0>;

    std::vector<VecView0> ranges_a;
    ranges_a.push_back(VecView0(v0a.begin(),v0a.end()));
    ranges_a.push_back(VecView0(v1a.begin(),v1a.end()));
    ranges_a.push_back(VecView0(v2a.begin(),v2a.end()));

    using ItType_1 = ConcatenatedForwardConstIterator<VecView0>;
    ItType_1 begin_a(ranges_a,0);
    ItType_1 end_a(ranges_a,IteratorState::pass_the_end);


    std::vector<VecView0> ranges_b;
    ranges_b.push_back(VecView0(v0b.begin(),v0b.end()));
    ranges_b.push_back(VecView0(v1b.begin(),v1b.end()));

    ItType_1 begin_b(ranges_b,0);
    ItType_1 end_b(ranges_b,IteratorState::pass_the_end);


    std::vector<VecView0> ranges;
    for (const auto &r : ranges_a)
        ranges.push_back(r);
    for (const auto &r : ranges_b)
        ranges.push_back(r);

    ItType_1 begin(ranges,0);
    ItType_1 end(ranges,IteratorState::pass_the_end);

    int i_a = 0;
    for (; begin_a != end_a ; ++begin_a, ++i_a)
        out << "i_a = " << i_a << "     value = " << *begin_a << std::endl;
    out << endl;

    int i_b = 0;
    for (; begin_b != end_b ; ++begin_b, ++i_b)
        out << "i_b = " << i_b << "     value = " << *begin_b << std::endl;
    out << endl;

    int i = 0;
    for (; begin != end ; ++begin, ++i)
        out << "i = " << i << "     value = " << *begin << std::endl;

    out << "========== do_test_2() --- end ==========" << endl;
    out << endl;
}


int main()
{
    do_test_1_const();
    do_test_1_nonconst();
    do_test_2();

    return 0;
}
