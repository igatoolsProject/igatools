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

#include <igatools/utils/concatenated_forward_iterator.h>
#include <vector>

using std::vector;


int main()
{
    vector<int> v0 = {1,2,3,4};
    vector<int> v1 = {5,6,7,8,9};
    vector<int> v2 = {10,11,12};
    vector<int> v3 = {13};

    using VecIterator = vector<int>::iterator;

    std::vector<std::pair<VecIterator,VecIterator>> ranges;
    ranges.push_back(std::make_pair<VecIterator,VecIterator>(v0.begin(),v0.end()));
    ranges.push_back(std::make_pair<VecIterator,VecIterator>(v1.begin(),v1.end()));
    ranges.push_back(std::make_pair<VecIterator,VecIterator>(v2.begin(),v2.end()));
    ranges.push_back(std::make_pair<VecIterator,VecIterator>(v3.begin(),v3.end()));


    ConcatenatedForwardIterator<VecIterator> begin(ranges,0);
    ConcatenatedForwardIterator<VecIterator> end(ranges,IteratorState::pass_the_end);



    using std::endl;
    int i = 0;
    for (; begin != end ; ++begin, ++i)
        out << "i = " << i << "     value = " << *begin << std::endl;
}
