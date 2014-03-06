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
//TODO put copyright here

#include <igatools/io/ig_reader.h>
#include <igatools/basis_functions/nurbs_element_accessor.h>

using std::vector;
using std::array;
using std::shared_ptr;
using std::string;

IGA_NAMESPACE_OPEN


template<int dim_ref_domain, int dim_phys_domain>
IgReader< dim_ref_domain, dim_phys_domain>::IgReader()
{}

template<int dim_ref_domain, int dim_phys_domain>
void
IgReader< dim_ref_domain, dim_phys_domain>::
load_xml(const std::string &filename)
{

    using boost::property_tree::ptree;
    ptree pt;

    read_xml(filename, pt);

    BOOST_FOREACH(ptree::value_type const& patch, pt.get_child("XMLFile.Patch"))
    {
        if (boost::iequals(patch.first, "<xmlattr>"))
        {
            const int dim_phys = patch.second.get<int>("DimPhysicalDomain");
            const int dim_ref  = patch.second.get<int>("DimReferenceDomain");

            AssertThrow(dim_phys == dim_phys_domain,
                        ExcDimensionMismatch(dim_phys,dim_phys_domain));
            AssertThrow(dim_ref == dim_ref_domain,
                        ExcDimensionMismatch(dim_ref,dim_ref_domain));

        }

        if (boost::iequals(patch.first, "KnotVector"))
        {

            BOOST_FOREACH(boost::property_tree::ptree::value_type const& knot, patch.second)
            {
                if (boost::iequals(knot.first, "<xmlattr>"))
                {
                    int d = knot.second.get<int>("Degree");
                    deg.push_back(d);
                }
                if (boost::iequals(knot.first, "BreakPoints"))
                {
                    std::string iss = knot.second.get<std::string>("");
                    std::vector<Real> row;
                    Real num;
                    std::stringstream line_stream(iss);
                    while (line_stream >> num) row.push_back(num);
                    breack_point.push_back(row);
                }
                if (boost::iequals(knot.first, "Multiplicities"))
                {
                    std::string iss = knot.second.get<std::string>("");
                    std::vector<Index> row;
                    Index num;
                    std::stringstream line_stream(iss);
                    while (line_stream >> num) row.push_back(num);
                    mlt.push_back(row);
                }
            }
        }

        if (boost::iequals(patch.first, "ControlPoints"))
        {
            vector<Real> weights_vector;
            BOOST_FOREACH(boost::property_tree::ptree::value_type const& cp, patch.second)
            {
                if (boost::iequals(cp.first, "NumDir"))
                {
                    std::string iss = cp.second.get<std::string>("");
                    std::stringstream line_stream(iss);
                    int num;
                    int direction_id = 0;
                    while (line_stream >> num)
                    {
                        AssertThrow(direction_id >= 0 && direction_id < dim_ref_domain,
                                    ExcIndexRange(direction_id,0,dim_ref_domain));
                        cp_per_ref_dir[direction_id] = num;
                        direction_id++;
                    }

                }
                if (boost::iequals(cp.first, "Coordinates"))
                {
                    std::string iss = cp.second.get<std::string>("");
                    std::vector<Real> row;
                    Real num;
                    std::stringstream line_stream(iss);
                    while (line_stream >> num) row.push_back(num);
                    control_point.push_back(row);
                }
                if (boost::iequals(cp.first, "Weights"))
                {
                    std::string iss = cp.second.get<std::string>("");
                    std::stringstream line_stream(iss);
                    Real num;
                    while (line_stream >> num) weights_vector.push_back(num);
                }
            }

            weights_.resize(cp_per_ref_dir);
            const int n_entries = weights_.flat_size();
            AssertThrow(weights_.flat_size() == weights_vector.size(),
                        ExcDimensionMismatch(weights_.flat_size(),weights_vector.size()));
            for (Size i=0 ; i < n_entries ; ++i)
                weights_(i) = weights_vector[i];

        }
    }//PATCH LOOP
    read = true;

    for (int i = 0; i<breack_point.size(); ++i)
    {
        n_knots[i] = breack_point[i].size();
        degree[i] = deg[i];
    }

    for (int i = 0 ; i < dim_ref_domain ; ++i)
        coord[i].resize(n_knots[i]);

    for (int i = 0 ; i < dim_ref_domain ; ++i)
        for (int j = 0 ; j < n_knots[i] ; ++j)
            coord[i][j] = breack_point[i][j];

    grid_ =  CartesianGrid<dim_ref_domain>::create(
                 CartesianProductArray< Real, dim_ref_domain>(coord));
}

template<int dim_ref_domain, int dim_phys_domain>
std::shared_ptr<NURBSSpace<dim_ref_domain,dim_phys_domain,1> >
IgReader< dim_ref_domain, dim_phys_domain>::get_nurbs_space()
{
    Assert(read,ExcMessage("If you want to get something from the IgReader, you should read first, use load functions"));

    for (int i = 0; i<breack_point.size(); ++i)
    {
        degree[i] = deg[i];
    }

    //mult.resize(mlt.size()) ;
    for (int i = 0 ; i < mlt.size() ; ++i)
    {
        mult.copy_data_direction(i, mlt[i]) ;
    }

    return NURBSSpace<dim_ref_domain,dim_phys_domain,1>::create(grid_,mult,degree,weights_);
}

template<int dim_ref_domain, int dim_phys_domain>
std::shared_ptr<IgMapping<NURBSSpace< dim_ref_domain,dim_phys_domain,1> > >
IgReader< dim_ref_domain, dim_phys_domain>::get_mapping_iga()
{
    Assert(read,ExcMessage("If you want to get something from the IgReader, you should read first, use load functions"));

    auto space_ptr = get_nurbs_space();

    Index l = 0;
    for (int i = 0; i < control_point.size(); ++i)
    {
        Index s = control_point[i].size();
        l = l+s;
    }

    std::vector<Index> global_dofs(l) ;
    for (Index dof = 0 ; dof < l ; ++dof)
        global_dofs[dof] = dof ;

    std::vector<Real> cp_ptr(global_dofs.size());

    std::array<Index, dim_phys_domain> off_set;

    off_set[0] = 0;
    for (int i = 1; i<control_point.size(); ++i)
        off_set[i] = off_set[i-1] + control_point[i-1].size();

    for (int i = 0; i<control_point.size(); ++i)
        for (int j = 0; j<control_point[i].size(); ++j)
        {
            Index el = off_set[i] + j;
            cp_ptr[el] = control_point[i][j];
        }

    return std::shared_ptr< IgMapping<NURBSSpace< dim_ref_domain,dim_phys_domain,1>> > (new
            IgMapping<NURBSSpace< dim_ref_domain,dim_phys_domain,1>>(space_ptr,cp_ptr));
}


template <int dim, int codim = 0>
std::shared_ptr< Mapping<dim,codim> >
ig_mapping_reader(const std::string &filename)
{
	const int dim_phys = dim + codim;

    using boost::property_tree::ptree;
    ptree xml_tree;

    read_xml(filename, xml_tree);


    using tree_value_t = typename ptree::value_type;


    LogStream out ;
    TensorIndex<dim> degree;
	CartesianProductArray<Real,dim> knots_unique_values;
	Multiplicity<dim> multiplicities;


    TensorSize<dim> n_ctrl_points_dim;
    array<vector<Real>,dim_phys> control_pts_coords;
    DynamicMultiArray<Real,dim> weights;

    bool is_nurbs_mapping = false;

    for (const tree_value_t & patch : xml_tree.get_child("XMLFile.Patch"))
    {
        if (boost::iequals(patch.first,"<xmlattr>"))
        {
            const int dim_ref_domain  = patch.second.get<int>("DimReferenceDomain");
            const int dim_phys_domain = patch.second.get<int>("DimPhysicalDomain");

            AssertThrow(dim == dim_ref_domain,
                        ExcDimensionMismatch(dim,dim_ref_domain));
            AssertThrow(dim_phys == dim_phys_domain,
                        ExcDimensionMismatch(dim_phys,dim_phys_domain));
        }


        if (boost::iequals(patch.first,"KnotVector"))
        {
            Index deg = -1 ;
            Index direction_id = -1 ;
            Size n_break_pts = 0 ;

            vector<Real> knots_unique_values_dir;
            vector<Index> multiplicities_dir;
            for (const auto & knot : patch.second)
            {
                if (boost::iequals(knot.first,"<xmlattr>"))
                {
                    deg = knot.second.get<Index>("Degree");
                    AssertThrow(deg >= 1, ExcLowerRange(deg,1));

                    direction_id = knot.second.get<Index>("Direction");
                    AssertThrow(direction_id >= 0, ExcLowerRange(direction_id,0));

                    n_break_pts = knot.second.get<Index>("NumBreakPoints");
                    AssertThrow(n_break_pts >= 2, ExcLowerRange(n_break_pts,2));
                }

                if (boost::iequals(knot.first, "BreakPoints"))
                {
                    Real knt;
                    std::stringstream line_stream(knot.second.get<std::string>(""));
                    while (line_stream >> knt)
                        knots_unique_values_dir.push_back(knt);

                	for (const auto & brk : knot.second)
                        if (boost::iequals(brk.first,"<xmlattr>"))
                        {
                        	const Index n_brk = brk.second.get<Index>("Num");
                        	AssertThrow(n_brk == knots_unique_values_dir.size(),
                        			ExcDimensionMismatch(n_brk,knots_unique_values_dir.size()));
                        }
                }

                if (boost::iequals(knot.first, "Multiplicities"))
                {
                    Index m;
                    std::stringstream line_stream(knot.second.get<std::string>(""));
                    while (line_stream >> m)
                        multiplicities_dir.push_back(m);

                	for (const auto & mlt : knot.second)
                        if (boost::iequals(mlt.first,"<xmlattr>"))
                        {
                        	const Index n_mlt = mlt.second.get<Index>("Num");
                        	AssertThrow(n_mlt == multiplicities_dir.size(),
                        			ExcDimensionMismatch(n_mlt,multiplicities_dir.size()));
                        }

                }
                //*/
            }
            AssertThrow(n_break_pts == knots_unique_values_dir.size(),
            		ExcDimensionMismatch(n_break_pts,knots_unique_values_dir.size()));
            AssertThrow(n_break_pts == multiplicities_dir.size(),
            		ExcDimensionMismatch(n_break_pts,multiplicities_dir.size()));

            degree[direction_id] = deg;
            knots_unique_values.copy_data_direction(direction_id,knots_unique_values_dir);
            multiplicities.copy_data_direction(direction_id,multiplicities_dir);
        }

        if (boost::iequals(patch.first, "ControlPoints"))
        {
        	Size n_ctrl_pts = 0;

            for (const auto &cp : patch.second)
            {
                if (boost::iequals(cp.first,"<xmlattr>"))
                {
                	n_ctrl_pts = cp.second.get<Index>("Num");
                    AssertThrow(n_ctrl_pts >= 0, ExcLowerRange(n_ctrl_pts,0));
                }

                vector<Size> tmp;
                if (boost::iequals(cp.first, "NumDir"))
                {
                    std::stringstream line_stream(cp.second.get<std::string>(""));
                    Index direction_id = 0;
                    while (line_stream >> n_ctrl_points_dim[direction_id])
                    {
                        AssertThrow(direction_id >= 0 && direction_id < dim,
                                    ExcIndexRange(direction_id,0,dim));
                        direction_id++;
                    }
                }


                if (boost::iequals(cp.first, "Coordinates"))
                {
                    vector<Real> control_pts_coord_dir;
                    Index direction_id = -1 ;

                    Real v;
                    std::stringstream line_stream(cp.second.get<std::string>(""));
                    while (line_stream >> v)
                    	control_pts_coord_dir.push_back(v);

                	for (const auto & coord : cp.second)
                        if (boost::iequals(coord.first,"<xmlattr>"))
                        {
                        	const Index n_coord = coord.second.get<Index>("Num");
                        	AssertThrow(n_coord == control_pts_coord_dir.size(),
                        			ExcDimensionMismatch(n_coord,control_pts_coord_dir.size()));

                            direction_id = coord.second.get<Index>("Dir");
                            AssertThrow(direction_id >= 0, ExcLowerRange(direction_id,0));

                            control_pts_coords[direction_id] = control_pts_coord_dir;
                        }
                }//*/


                if (boost::iequals(cp.first, "Weights"))
                {
                	is_nurbs_mapping = true;

                	vector<Real> weights_vec;
                    std::stringstream line_stream(cp.second.get<std::string>(""));
                    Real w;
                    while (line_stream >> w)
                    	weights_vec.push_back(w);

                	for (const auto & wght : cp.second)
                        if (boost::iequals(wght.first,"<xmlattr>"))
                        {
                        	const Index n_weights = wght.second.get<Index>("Num");
                        	AssertThrow(n_weights == weights_vec.size(),
                        			ExcDimensionMismatch(n_weights,weights_vec.size()));
                        }

                	weights.resize(n_ctrl_points_dim);

                	AssertThrow(weights_vec.size() == weights.flat_size(),
                	     ExcDimensionMismatch(weights_vec.size(),weights.flat_size()));

                	Index flat_id = 0;
                	for ( auto & w : weights)
                		w =  weights_vec[flat_id++];
                }
                //*/
            }
            AssertThrow(n_ctrl_pts == n_ctrl_points_dim.flat_size(),
            		ExcDimensionMismatch(n_ctrl_pts,n_ctrl_points_dim.flat_size()));

            if (is_nurbs_mapping)
            {
            	AssertThrow(n_ctrl_pts == weights.flat_size(),
            	     ExcDimensionMismatch(n_ctrl_pts,weights.flat_size()));
            }
            /*
                        weights_.resize(cp_per_ref_dir);
                        const int n_entries = weights_.flat_size();
                        AssertThrow(weights_.flat_size() == weights_vector.size(),
                                    ExcDimensionMismatch(weights_.flat_size(),weights_vector.size()));
                        for (Size i=0 ; i < n_entries ; ++i)
                            weights_(i) = weights_vector[i];
            //*/
        }

    }//PATCH LOOP

    vector<Real> control_pts;
    for (int i = 0 ; i < dim_phys ; ++i)
        control_pts.insert(control_pts.end(),control_pts_coords[i].begin(),control_pts_coords[i].end());



    auto grid = CartesianGrid<dim>::create(knots_unique_values);
    shared_ptr< Mapping<dim,codim> > map;

    if (is_nurbs_mapping)
    {
    	using space_t = NURBSSpace<dim,dim_phys,1>;
    	auto space = space_t::create(
    			grid,
    			StaticMultiArray<Multiplicity<dim>,dim_phys,1>(multiplicities),
    			StaticMultiArray<TensorIndex<dim>,dim_phys,1>(degree),
    			StaticMultiArray<DynamicMultiArray<Real,dim>,dim_phys,1>(weights));

        map = IgMapping<space_t>::create(space,control_pts);
    }
    else
    {
    	using space_t = BSplineSpace<dim,dim_phys,1>;
    	auto space = space_t::create(
    			grid,
    			StaticMultiArray<Multiplicity<dim>,dim_phys,1>(multiplicities),
    			StaticMultiArray<TensorIndex<dim>,dim_phys,1>(degree));

        map = IgMapping<space_t>::create(space,control_pts);
    }
//    map->print_info(out);


    return map;
}



IGA_NAMESPACE_CLOSE


#include <igatools/io/ig_reader.inst>

