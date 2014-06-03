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
 * martinelli
 * 30 May 2014
 *
 */

#include "../tests.h"

#include <igatools/base/exceptions.h>
#include <vector>

using std::vector;

template <class Iterator>
class ConcatenatedIterator
{
public:
    void push_back(const Iterator &begin,const Iterator &end)
    {
        Assert(begin < end,ExcInvalidIterator());
        ranges_.push_back(std::make_pair(begin,end));
    }



private:

    std::vector<std::pair<Iterator,Iterator>> ranges_;
};



int main()
{
    vector<int> v0 = {1,2,3,4};
    vector<int> v1 = {2,3,4,5};

    using VecIterator = vector<int>::iterator;

    ConcatenatedIterator<VecIterator> concatenated_iterator;
    concatenated_iterator.push_back(v0.begin(),v0.end());
    concatenated_iterator.push_back(v1.begin(),v1.end());

}
