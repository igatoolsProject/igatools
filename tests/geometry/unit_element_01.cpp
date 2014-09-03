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
 *  Test of unit element
 *  author: pauletti
 *  date: 2014-08-18
 *
 */

#include "../tests.h"

#include <igatools/geometry/unit_element.h>



template <int dim>
void skeleton()
{
    OUTSTART

    for (int k = 0; k<=dim; ++k)
    	out << UnitElement<dim>::skeleton_size[k] << endl;

    OUTEND
}

template <int dim, int k>
EnableIf< (dim==k) || (k<0),
std::array<typename UnitElement<dim>::template Skeleton<k>, skel_size(dim, k)>>
fill_skeleton()
{
    std::array<typename UnitElement<dim>::template Skeleton<k>, skel_size(dim, k)> res;
    return res;
}


template <int dim, int k>
EnableIf< (dim>k) && (k>=0),
std::array<typename UnitElement<dim>::template Skeleton<k>, skel_size(dim, k)>>
fill_skeleton()
{
    std::array<typename UnitElement<dim>::template Skeleton<k>, skel_size(dim, k)> res;

    auto skel_dim_1 = fill_skeleton<dim-1, k>();
    auto skel_dim_1_0 = fill_skeleton<dim-1, k-1>();
    auto polygon = res.begin();

    for (auto &polygon_dim_1 : skel_dim_1_0)
    {
    	auto &dirs_dim_1 = polygon_dim_1.constant_directions;
    	auto &dirs       = polygon->constant_directions;
    	std::copy(dirs_dim_1.begin(), dirs_dim_1.end(), dirs.begin());
    	polygon->constant_values = polygon_dim_1.constant_values;
    	++polygon;
    }

    for (auto &polygon_dim_1 : skel_dim_1)
    {
        auto &dirs_dim_1 = polygon_dim_1.constant_directions;
        auto &values_1 = polygon_dim_1.constant_values;
        for (int j = 0; j<2; ++j)
        {
        	auto &dirs       = polygon->constant_directions;
        	auto &values       = polygon->constant_values;
        	std::copy(dirs_dim_1.begin(), dirs_dim_1.end(), dirs.begin());
        	dirs[dim - k -1] = dim-1;
        	std::copy(values_1.begin(), values_1.end(), values.begin());
        	++polygon;
        }
    }
    return res;
}


template<int dim, int sub_dim>
void describe_skeleton()
{
	auto faces = fill_skeleton<dim, sub_dim>();
	for (auto &face : faces)
	{
		out << "polygon const direction: ";
		for (auto &dir : face.constant_directions)
		{
			out << dir << " ";
		}
		out << endl;
	}
}


int main()
{
    out.depth_console(20);

    describe_skeleton<1,1>();
    describe_skeleton<1,0>();

    describe_skeleton<2,2>();
    describe_skeleton<2,1>();
    describe_skeleton<2,0>();

    describe_skeleton<3,3>();
    describe_skeleton<3,2>();
    describe_skeleton<3,1>();
    describe_skeleton<3,0>();

//    skeleton<1>();
 //   skeleton<2>();
//    skeleton<3>();

    return 0;
}
