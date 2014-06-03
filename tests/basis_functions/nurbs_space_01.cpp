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

#include "../tests.h"

#include <igatools/basis_functions/nurbs_space.h>
#include <igatools/base/exceptions.h>


#include <iostream>
using std::cout ;
using std::endl ;

#include <vector>
using std::vector ;

#include <array>
using std::array ;

#include <memory>
using std::shared_ptr ;







template< int dim_domain, int dim_range, int rank >
void do_test()
{
    Assert(dim_domain == 1 || dim_domain == 2 || dim_domain == 3, ExcIndexRange(dim_domain, 1, 4)) ;


//  LogStream cout( rcp( new std::ostream( std::cout.rdbuf() ) ), "\t" ) ;



    //----------------------------------------------------------------------------------------------
    // begin : testing the constructor
    vector< iga::Real > coord_x ;
    coord_x.push_back(0.0) ;
    coord_x.push_back(1.0) ;
    coord_x.push_back(2.0) ;
    coord_x.push_back(3.0) ;
    coord_x.push_back(4.0) ;

    vector< iga::Real > coord_y ;
    coord_y.push_back(5.0) ;
    coord_y.push_back(6.0) ;
    coord_y.push_back(7.0) ;
    coord_y.push_back(8.0) ;

    vector< iga::Real > coord_z ;
    coord_z.push_back(9.0) ;
    coord_z.push_back(10.0) ;
    coord_z.push_back(11.0) ;

    CartesianProductArray< iga::Real, dim_domain> coord ;
    CartesianProductArray<Index , dim_domain>  mult ;
    TensorIndex<dim_domain> degree ;

    if (dim_domain == 1)
    {
        coord.copy_data_direction(0,coord_x) ;
        degree[0] = 3 ;
    }
    else if (dim_domain == 2)
    {
        coord.copy_data_direction(0,coord_x);
        coord.copy_data_direction(1,coord_y);

        degree[0] = 3 ;
        degree[1] = 2 ;
    }
    else if (dim_domain == 3)
    {
        coord.copy_data_direction(0,coord_x);
        coord.copy_data_direction(1,coord_y);
        coord.copy_data_direction(2,coord_z);

        degree[0] = 3 ;
        degree[1] = 2 ;
        degree[2] = 1 ;
    }




    auto  knots = CartesianGrid<dim_domain>::create(coord) ;
    StaticMultiArray<TensorIndex<dim_domain>,dim_range,rank> deg;
    deg.fill(degree);
    NURBSSpace< dim_domain, dim_range, rank > nurbs_space(deg, knots) ;
    nurbs_space.print_info(out);
    out << endl ;
    // end : testing the constructor
    //----------------------------------------------------------------------------------------------

}


int main(int argc, char *argv[])
{
    do_test< 1, 1, 1 >() ;
    do_test< 1, 2, 1 >() ;
    do_test< 1, 3, 1 >() ;

    do_test< 2, 1, 1 >() ;
    do_test< 2, 2, 1 >() ;
    do_test< 2, 3, 1 >() ;

    do_test< 3, 1, 1 >() ;
    do_test< 3, 3, 1 >() ;

    return 0;
}
