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





#ifndef __IG_READER_H_
#define __IG_READER_H_

#include <igatools/base/exceptions.h>

#include <igatools/utils/cartesian_product_array.h>
#include <igatools/geometry/cartesian_grid.h>
#include <igatools/geometry/ig_mapping.h>
#include <igatools/basis_functions/nurbs_space.h>
#include <igatools/basis_functions/physical_space.h>


#include <boost/algorithm/string.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <string>
#include <set>
#include <exception>
#include <vector>
#include <iostream>

#include <array>
#include <memory>

IGA_NAMESPACE_OPEN

//TODO: This class has to be better coded and documented
//TODO: Document the xml format
//TODO: get_mapping_iga should return a base type
template<int dim_ref_domain, int dim_phys_domain=dim_ref_domain>
class IgReader
{
public:
    IgReader();

    void load_xml(const std::string &filename);

    std::shared_ptr<CartesianGrid<dim_ref_domain> > get_cartesian_grid();

    std::shared_ptr<NURBSSpace< dim_ref_domain,dim_phys_domain,1> >
    get_nurbs_space();

    std::shared_ptr<IgMapping<NURBSSpace< dim_ref_domain,dim_phys_domain,1>> >
            get_mapping_iga();

private:
    bool read = false; // use this flag to say that data has been read
    std::vector<int> deg;
    std::vector<std::vector<Index> > mlt;
    std::vector<std::vector<Real> > breack_point;
    std::vector<std::vector<Real> > control_point;
    DynamicMultiArray<Real,dim_ref_domain> weights_;
    TensorSize<dim_ref_domain> cp_per_ref_dir;

    std::array< int, dim_ref_domain > n_knots;
    std::array< int, dim_ref_domain > degree;

    std::array< std::vector< Real >, dim_ref_domain > coord;
    Multiplicity< dim_ref_domain > mult;

    std::shared_ptr<CartesianGrid<dim_ref_domain> > grid_;

};



template <int dim, int codim = 0>
std::shared_ptr< Mapping<dim,codim> >
ig_mapping_reader(const std::string &filename);

IGA_NAMESPACE_CLOSE


#endif

