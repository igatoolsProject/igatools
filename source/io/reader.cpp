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

#include <igatools/io/reader.h>
#include <igatools/base/exceptions.h>
//#include <igatools/utils/array.h>
#include <igatools/base/ig_function.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/basis_functions/nurbs_element.h>


#include <boost/algorithm/string.hpp>
#include <boost/property_tree/xml_parser.hpp>


#include <set>

using std::string;
using std::shared_ptr;
using std::set;

using namespace iga::reader_utils;


IGA_NAMESPACE_OPEN


namespace reader_utils
{
boost::property_tree::ptree
get_xml_tree(const string &filename)
{
    using boost::property_tree::ptree;
    ptree xml_tree;

    read_xml(filename, xml_tree);

    return xml_tree;
}


boost::property_tree::ptree
get_xml_element_attributes(const boost::property_tree::ptree &element)
{
    return get_xml_element(element,"<xmlattr>");
}




Size
count_xml_elements_same_tag(
    const boost::property_tree::ptree &tree,
    const string &tag_name)
{
    Size counter = 0;
    for (const auto &leaf : tree)
        if (boost::iequals(leaf.first,tag_name))
            counter++;

    return counter;
}



bool
xml_element_is_present(
    const boost::property_tree::ptree &tree,
    const string &tag_name)
{
    return (count_xml_elements_same_tag(tree,tag_name) > 0) ? true:false;
}



bool
xml_element_is_unique(
    const boost::property_tree::ptree &tree,
    const string &tag_name)
{
    const Size counter = count_xml_elements_same_tag(tree,tag_name);
    AssertThrow(counter >= 1, ExcLowerRange(counter,1));
    return (counter == 1) ? true:false;
}



vector< boost::property_tree::ptree >
get_xml_element_vector(
    const boost::property_tree::ptree &tree,
    const string &tag_name)
{
    vector<boost::property_tree::ptree> element;
    for (const auto &leaf : tree)
        if (boost::iequals(leaf.first,tag_name))
            element.push_back(leaf.second);

    AssertThrow(element.size() > 0,ExcMessage("XML element <"+tag_name+"> not found."));

    return element;
}



boost::property_tree::ptree
get_xml_element(
    const boost::property_tree::ptree &tree,
    const string &tag_name)
{
    const auto &elements = get_xml_element_vector(tree,tag_name);

    AssertThrow(elements.size() == 1,ExcDimensionMismatch(elements.size(),1));

    return elements[0];
}






template <class ScalarType>
vector<ScalarType>
get_vector_data_from_xml(const boost::property_tree::ptree &tree)
{
    vector<ScalarType> data;

    ScalarType v;
    std::stringstream line_stream(tree.get<std::string>(""));
    while (line_stream >> v)
        data.push_back(v);

    return data;
}



template <int dim>
CartesianProductArray<Size,dim>
get_interior_multiplicity_from_xml(const boost::property_tree::ptree &tree)
{
    //-------------------------------------------------------------------------
    // reading the Multiplicity attributes
    const auto &mult_attributes = get_xml_element_attributes(tree);

    const int dim_from_file = mult_attributes.get<int>("Dim");
    AssertThrow(dim == dim_from_file,
                ExcDimensionMismatch(dim,dim_from_file));
    //-------------------------------------------------------------------------


    //-------------------------------------------------------------------------
    // reading the multiplicity along each direction
    const auto &mult_elements = get_xml_element_vector(tree,"InteriorMultiplicity");
    AssertThrow(mult_elements.size() == dim,ExcDimensionMismatch(mult_elements.size(),dim));

    std::array<vector<int>,dim> mlt_data;
    for (const auto mlt_element : mult_elements)
    {
        const auto &mlt_attributes = get_xml_element_attributes(mlt_element);

        const int direction_id = mlt_attributes.get<int>("Direction");
        mlt_data[direction_id] = get_vector_data_from_xml<Size>(mlt_element);

        const int n_mlt_from_file = mlt_attributes.get<int>("Size");
        AssertThrow(mlt_data[direction_id].size() == n_mlt_from_file,
                    ExcDimensionMismatch(mlt_data[direction_id].size(),n_mlt_from_file));
    }
    //-------------------------------------------------------------------------
    CartesianProductArray<Size,dim> mult(mlt_data);

    return mult;
}

} // end of namespace reader_utils


string
get_xml_input_file_format(const std::string &filename)
{
    string file_format_version;

    const auto &xml_tree = get_xml_tree(filename);
    for (const auto &leaf : xml_tree)
    {
        if (boost::iequals(leaf.first,"XMLFile"))
        {
            file_format_version = "1.0";
        }
        else if (boost::iequals(leaf.first,"Igatools"))
        {
            const auto &igatools_attr = get_xml_element_attributes(leaf.second);
            file_format_version = igatools_attr.get<std::string>("FormatVersion");
        }
        else
        {
            string err_msg("XML input file is not in a valid format for igatools.");
            AssertThrow(false,ExcMessage(err_msg));
        }
    }

    return file_format_version;
}



template <int dim, int codim = 0>
std::shared_ptr< MapFunction<dim,dim+codim> >
get_mapping_from_file(const std::string &filename)
{
    const string file_format_version = get_xml_input_file_format(filename) ;

    shared_ptr< MapFunction<dim,dim+codim> > map;
    if (file_format_version == "1.0")
    {
        Assert(false,ExcNotImplemented());
        AssertThrow(false,ExcNotImplemented());
        // use the reader for format version 1.0
//        map = ig_mapping_reader_version_1_0<dim,codim>(filename);
    }
    else if (file_format_version == "2.0")
    {
        const auto &file_tree = get_xml_tree(filename);

        const auto &igatools_tree = get_xml_element(file_tree,"Igatools");

        // use the format version 2.0
        map = get_mapping_from_xml<dim,codim>(igatools_tree);
    }
    else
    {
        string err_msg("Input file format version "+file_format_version+" not supported.");
        AssertThrow(false,ExcMessage(err_msg));
    }
    AssertThrow(map != nullptr,ExcNullPtr());

    return map;
}







template <int dim, int codim = 0>
std::shared_ptr< MapFunction<dim,dim+codim> >
get_ig_mapping_from_xml(const boost::property_tree::ptree &igatools_tree)
{
    AssertThrow(xml_element_is_unique(igatools_tree,"IgMapping"),
                ExcMessage("The IgMapping tag is not unique."));

    const auto &mapping_tree = get_xml_element(igatools_tree,"IgMapping");


    //-------------------------------------------------------------------------
    // reading the IgMapping attributes
    const auto &mapping_attributes = get_xml_element_attributes(mapping_tree);

    const int dim_from_file = mapping_attributes.get<int>("Dim");
    AssertThrow(dim == dim_from_file,
                ExcDimensionMismatch(dim,dim_from_file));

    const int codim_from_file = mapping_attributes.get<int>("Codim");
    AssertThrow(codim == codim_from_file,
                ExcDimensionMismatch(codim,codim_from_file));

    const string ref_space_type = mapping_attributes.get<string>("RefSpaceType");
    AssertThrow(ref_space_type == "BSplineSpace" || ref_space_type == "NURBSSpace",
                ExcMessage("Unknown reference space type: "+ref_space_type));

    const bool is_nurbs_space = (ref_space_type == "NURBSSpace")?true:false;
    //-------------------------------------------------------------------------


    //-------------------------------------------------------------------------
    // reading the control points
    const auto &ctrl_pts_tree = get_xml_element(mapping_tree,"ControlPoints");

    const auto &ctrl_pts_attributes = get_xml_element_attributes(ctrl_pts_tree);

    const int dim_ctrl_pts = ctrl_pts_attributes.get<int>("Dim");
    AssertThrow(dim_ctrl_pts == 1,ExcDimensionMismatch(dim_ctrl_pts,1));

    const int n_ctrl_pts = ctrl_pts_attributes.get<int>("Size");
    AssertThrow(n_ctrl_pts >= 1,ExcLowerRange(n_ctrl_pts,1));

    vector<Real> cntrl_pts = get_vector_data_from_xml<Real>(ctrl_pts_tree);
    AssertThrow(cntrl_pts.size() == n_ctrl_pts,ExcDimensionMismatch(cntrl_pts.size(),n_ctrl_pts));
    //-------------------------------------------------------------------------


    //-------------------------------------------------------------------------
    // reading the reference space
    shared_ptr< MapFunction<dim,dim+codim> > map;

    const int dim_phys = dim + codim;
    if (is_nurbs_space)
    {
#ifdef NURBS
        using ref_space_t = NURBSSpace<dim,dim_phys,1>;
        auto ref_space = get_nurbs_space_from_xml<dim,dim_phys,1>(mapping_tree);

        map = IgFunction<ref_space_t>::create(ref_space,cntrl_pts);
#else
        Assert(false,ExcMessage("NURBS support disabled from configuration cmake parameters."));
        AssertThrow(false,ExcMessage("NURBS support disabled from configuration cmake parameters."));
#endif
    }
    else
    {
        using ref_space_t = BSplineSpace<dim,dim_phys,1>;
        auto ref_space = get_bspline_space_from_xml<dim,dim_phys,1>(mapping_tree);

        map = IgFunction<ref_space_t>::create(ref_space,cntrl_pts);
    }
    //-------------------------------------------------------------------------
    AssertThrow(map != nullptr,ExcNullPtr());

    return map;
}


template <int dim, int codim = 0>
std::shared_ptr< MapFunction<dim,dim+codim> >
get_mapping_from_xml(const boost::property_tree::ptree &igatools_tree)
{
    shared_ptr< MapFunction<dim,dim+codim> > map;

    if (xml_element_is_unique(igatools_tree,"IgMapping"))
    {
        // use the reader for format version 2.0
        map = get_ig_mapping_from_xml<dim,codim>(igatools_tree);
    }
    else
    {
        AssertThrow(false,ExcNotImplemented());
    }
    AssertThrow(map != nullptr,ExcNullPtr());

    return map;
}



template <int dim>
shared_ptr< CartesianGrid<dim> >
get_cartesian_grid_from_xml(const boost::property_tree::ptree &tree)
{
    AssertThrow(xml_element_is_unique(tree,"CartesianGrid"),
                ExcMessage("The CartesianGrid tag is not unique."));

    const auto &grid_tree = get_xml_element(tree,"CartesianGrid");

    //-------------------------------------------------------------------------
    // reading the CartesianGrid attributes
    const auto &grid_attributes = get_xml_element_attributes(grid_tree);

    const int dim_from_file = grid_attributes.get<int>("Dim");
    AssertThrow(dim == dim_from_file,
                ExcDimensionMismatch(dim,dim_from_file));
    //-------------------------------------------------------------------------


    //-------------------------------------------------------------------------
    // reading the knots
    const auto &knots_elements = get_xml_element_vector(grid_tree,"Knots");
    AssertThrow(knots_elements.size() == dim,ExcDimensionMismatch(knots_elements.size(),dim));

    std::array<vector<Real>,dim> knots;
    for (int i = 0 ; i < dim ; ++i)
    {
        const auto &knots_attributes = get_xml_element_attributes(knots_elements[i]);

        const int direction_id = knots_attributes.get<int>("Direction");
        knots[direction_id] = get_vector_data_from_xml<Real>(knots_elements[i]);

        const int n_knts_from_file = knots_attributes.get<int>("Size");
        AssertThrow(knots[direction_id].size() == n_knts_from_file,
                    ExcDimensionMismatch(knots[direction_id].size(),n_knts_from_file));
    }
    //-------------------------------------------------------------------------

    return CartesianGrid<dim>::create(CartesianProductArray<Real,dim>(knots));
}





template <int dim, int range, int rank>
shared_ptr< BSplineSpace<dim,range,rank> >
get_bspline_space_from_xml(const boost::property_tree::ptree &tree)
{
    AssertThrow(xml_element_is_unique(tree,"BSplineSpace"),
                ExcMessage("The BSplineSpace tag is not unique."));

    const auto &ref_space_tree = get_xml_element(tree,"BSplineSpace");

    using space_t = BSplineSpace<dim,range,rank>;

    //-------------------------------------------------------------------------
    // reading the BSplineSpace attributes
    const auto &ref_space_attributes = get_xml_element_attributes(ref_space_tree);

    const int dim_from_file = ref_space_attributes.get<int>("Dim");
    AssertThrow(dim == dim_from_file,
                ExcDimensionMismatch(dim,dim_from_file));

    const int range_from_file = ref_space_attributes.get<int>("Range");
    AssertThrow(range == range_from_file,
                ExcDimensionMismatch(range,range_from_file));

    const int rank_from_file = ref_space_attributes.get<int>("Rank");
    AssertThrow(rank == rank_from_file,
                ExcDimensionMismatch(range,rank_from_file));
    //-------------------------------------------------------------------------


    //-------------------------------------------------------------------------
    // reading the CartesianGrid
    auto grid = get_cartesian_grid_from_xml<dim>(ref_space_tree);
    //-------------------------------------------------------------------------


    //-------------------------------------------------------------------------
    // reading the ScalarComponents
    const auto &scalar_components_tree = get_xml_element(ref_space_tree,"BSplineSpaceScalarComponents");
    const auto &scalar_components_attributes = get_xml_element_attributes(scalar_components_tree);
    const int n_sc_components_from_file = scalar_components_attributes.get<int>("Size");
    AssertThrow(n_sc_components_from_file >= 1 && n_sc_components_from_file <= space_t::n_components,
                ExcIndexRange(n_sc_components_from_file,1,space_t::n_components+1));


    //-----
    // components_map
    const auto &components_map_tree = get_xml_element(scalar_components_tree,"ComponentsMap");
    const auto &components_map_attributes = get_xml_element_attributes(components_map_tree);
    const int n_components_map_entries = components_map_attributes.get<int>("Size");
    AssertThrow(space_t::n_components == n_components_map_entries,
                ExcDimensionMismatch(space_t::n_components,n_components_map_entries));
    vector<Index> components_map_vec = get_vector_data_from_xml<Index>(components_map_tree);
    std::array<Index,space_t::n_components> components_map;
    for (Index i = 0 ; i < space_t::n_components ; ++i)
        components_map[i] = components_map_vec[i];

    set<Index> active_components_id(components_map_vec.begin(),components_map_vec.end());
    const int n_active_components = active_components_id.size();
    //-----



    const auto &scalar_component_vector = get_xml_element_vector(scalar_components_tree,"BSplineSpaceScalarComponent");
    AssertThrow(scalar_component_vector.size() == n_sc_components_from_file,
                ExcDimensionMismatch(scalar_component_vector.size(),n_sc_components_from_file));
    AssertThrow(scalar_component_vector.size() == n_active_components,
                ExcDimensionMismatch(scalar_component_vector.size(),n_active_components));

    using       DegreeTable = typename space_t::DegreeTable;
    using MultiplicityTable = typename space_t::MultiplicityTable;

    DegreeTable degrees(components_map);
    shared_ptr<MultiplicityTable> multiplicities(new MultiplicityTable(components_map));

    for (const auto &comp_element : scalar_component_vector)
    {
        const auto &component_attributes = get_xml_element_attributes(comp_element);
        const int comp_id = component_attributes.get<int>("Id");
        AssertThrow(comp_id >=0 && comp_id < space_t::n_components,
                    ExcIndexRange(comp_id,0,space_t::n_components));


        //---------------------------------------------
        // getting the dofs tensor size
        const auto &comp_dofs_tens_size_element = get_xml_element(comp_element,"DofsTensorSize");
        const auto &comp_dofs_tens_size_attributes = get_xml_element_attributes(comp_dofs_tens_size_element);
        const int dofs_rank = comp_dofs_tens_size_attributes.get<int>("Dim");
        AssertThrow(dofs_rank == dim,ExcDimensionMismatch(dofs_rank,dim));

        vector<Index> dofs_size_vec = get_vector_data_from_xml<Index>(comp_dofs_tens_size_element);
        AssertThrow(dofs_rank == dofs_size_vec.size(),ExcDimensionMismatch(dofs_rank,dofs_size_vec.size()));


        TensorSize<dim> dofs_size;
        for (int i = 0 ; i < dim ; ++i)
            dofs_size[i] = dofs_size_vec[i];
        //---------------------------------------------


        //---------------------------------------------
        // getting the degrees
        const auto &comp_degree_element = get_xml_element(comp_element,"Degrees");
        const auto &comp_degree_attributes = get_xml_element_attributes(comp_degree_element);
        const int n_degree = comp_degree_attributes.get<int>("Dim");
        AssertThrow(n_degree == dim,ExcDimensionMismatch(n_degree,dim));

        vector<Index> degrees_vec = get_vector_data_from_xml<Index>(comp_degree_element);
        AssertThrow(n_degree == degrees_vec.size(),ExcDimensionMismatch(n_degree,degrees_vec.size()));

        for (int i = 0 ; i < dim ; ++i)
            degrees[comp_id][i] = degrees_vec[i];
        //---------------------------------------------


        //---------------------------------------------
        // getting the multiplicities
        const auto &comp_multiplicities_element = get_xml_element(comp_element,"InteriorMultiplicities");
        (*multiplicities)[comp_id] = get_interior_multiplicity_from_xml<dim>(comp_multiplicities_element);
        //---------------------------------------------

    } // end loop over the scalar components
    //-------------------------------------------------------------------------


    // TODO (pauletti, Dec 26, 2014): read periodic, end_behaviour and boundary knots from file
    typename space_t::EndBehaviourTable
	end_behaviour(components_map, filled_array<BasisEndBehaviour, dim>(BasisEndBehaviour::interpolatory));
    typename space_t::PeriodicTable periodic(components_map, filled_array<bool, dim>(false));

    auto ref_space = space_t::create(degrees, grid, multiplicities, periodic, end_behaviour);

    return ref_space;
}

#ifdef NURBS
template <int dim, int range, int rank>
shared_ptr< NURBSSpace<dim,range,rank> >
get_nurbs_space_from_xml(const boost::property_tree::ptree &tree)
{
//    Assert(false,ExcMessage("This function must be checked (MM:26/11/2014!"));
//    AssertThrow(false,ExcMessage("This function must be checked (MM:26/11/2014!"));

    AssertThrow(xml_element_is_unique(tree,"NURBSSpace"),
                ExcMessage("The NURBSSpace tag is not unique."));

    const auto &ref_space_tree = get_xml_element(tree,"NURBSSpace");

    using space_t = NURBSSpace<dim,range,rank>;

    //-------------------------------------------------------------------------
    // reading the NURBSSpace attributes
    const auto &ref_space_attributes = get_xml_element_attributes(ref_space_tree);

    const int dim_from_file = ref_space_attributes.get<int>("Dim");
    AssertThrow(dim == dim_from_file,
                ExcDimensionMismatch(dim,dim_from_file));

    const int range_from_file = ref_space_attributes.get<int>("Range");
    AssertThrow(range == range_from_file,
                ExcDimensionMismatch(range,range_from_file));

    const int rank_from_file = ref_space_attributes.get<int>("Rank");
    AssertThrow(rank == rank_from_file,
                ExcDimensionMismatch(range,rank_from_file));
    //-------------------------------------------------------------------------


    //-------------------------------------------------------------------------
    // reading the CartesianGrid
    auto grid = get_cartesian_grid_from_xml<dim>(ref_space_tree);
    //-------------------------------------------------------------------------


    //-------------------------------------------------------------------------
    // reading the ScalarComponents
    const auto &scalar_components_tree = get_xml_element(ref_space_tree,"NURBSSpaceScalarComponents");
    const auto &scalar_components_attributes = get_xml_element_attributes(scalar_components_tree);
    const int n_sc_components_from_file = scalar_components_attributes.get<int>("Size");
    AssertThrow(n_sc_components_from_file >= 1 && n_sc_components_from_file <= space_t::n_components,
                ExcIndexRange(n_sc_components_from_file,1,space_t::n_components+1));


    //-----
    // components_map
    const auto &components_map_tree = get_xml_element(scalar_components_tree,"ComponentsMap");
    const auto &components_map_attributes = get_xml_element_attributes(components_map_tree);
    const int n_components_map_entries = components_map_attributes.get<int>("Size");
    AssertThrow(space_t::n_components == n_components_map_entries,
                ExcDimensionMismatch(space_t::n_components,n_components_map_entries));
    vector<Index> components_map_vec = get_vector_data_from_xml<Index>(components_map_tree);
    array<Index,space_t::n_components> components_map;
    for (Index i = 0 ; i < space_t::n_components ; ++i)
        components_map[i] = components_map_vec[i];

    set<Index> active_components_id(components_map_vec.begin(),components_map_vec.end());
    const int n_active_components = active_components_id.size();
    //-----



    const auto &scalar_component_vector = get_xml_element_vector(scalar_components_tree,"NURBSSpaceScalarComponent");
    AssertThrow(scalar_component_vector.size() == n_sc_components_from_file,
                ExcDimensionMismatch(scalar_component_vector.size(),n_sc_components_from_file));
    AssertThrow(scalar_component_vector.size() == n_active_components,
                ExcDimensionMismatch(scalar_component_vector.size(),n_active_components));

    using       DegreeTable = typename space_t::DegreeTable;
    using      WeightsTable = typename space_t::WeightsTable;
    using MultiplicityTable = typename space_t::MultiplicityTable;

    DegreeTable degrees(components_map);
    WeightsTable weights(components_map);
    shared_ptr<MultiplicityTable> multiplicities(new MultiplicityTable(components_map));

    for (const auto &comp_element : scalar_component_vector)
    {
        const auto &component_attributes = get_xml_element_attributes(comp_element);
        const int comp_id = component_attributes.get<int>("Id");
        AssertThrow(comp_id >=0 && comp_id < space_t::n_components,
                    ExcIndexRange(comp_id,0,space_t::n_components));


        //---------------------------------------------
        // getting the dofs tensor size
        const auto &comp_dofs_tens_size_element = get_xml_element(comp_element,"DofsTensorSize");
        const auto &comp_dofs_tens_size_attributes = get_xml_element_attributes(comp_dofs_tens_size_element);
        const int dofs_rank = comp_dofs_tens_size_attributes.get<int>("Dim");
        AssertThrow(dofs_rank == dim,ExcDimensionMismatch(dofs_rank,dim));

        vector<Index> dofs_size_vec = get_vector_data_from_xml<Index>(comp_dofs_tens_size_element);
        AssertThrow(dofs_rank == dofs_size_vec.size(),ExcDimensionMismatch(dofs_rank,dofs_size_vec.size()));


        TensorSize<dim> dofs_size;
        for (int i = 0 ; i < dim ; ++i)
            dofs_size[i] = dofs_size_vec[i];
        //---------------------------------------------


        //---------------------------------------------
        // getting the degrees
        const auto &comp_degree_element = get_xml_element(comp_element,"Degrees");
        const auto &comp_degree_attributes = get_xml_element_attributes(comp_degree_element);
        const int n_degree = comp_degree_attributes.get<int>("Dim");
        AssertThrow(n_degree == dim,ExcDimensionMismatch(n_degree,dim));

        vector<Index> degrees_vec = get_vector_data_from_xml<Index>(comp_degree_element);
        AssertThrow(n_degree == degrees_vec.size(),ExcDimensionMismatch(n_degree,degrees_vec.size()));

        for (int i = 0 ; i < dim ; ++i)
            degrees[comp_id][i] = degrees_vec[i];
        //---------------------------------------------


        //---------------------------------------------
        // getting the multiplicities
        const auto &comp_multiplicities_element = get_xml_element(comp_element,"InteriorMultiplicities");
        (*multiplicities)[comp_id] = get_interior_multiplicity_from_xml<dim>(comp_multiplicities_element);
        //---------------------------------------------



        //-------------------------------------------------------------------------
        // reading the weights
        const auto &weights_tree = get_xml_element(comp_element,"Weights");
        const auto &weights_attributes = get_xml_element_attributes(weights_tree);

        const int n_weights = weights_attributes.get<int>("Size");
        AssertThrow(n_weights == dofs_size.flat_size(),ExcLowerRange(n_weights,dofs_size.flat_size()));

        vector<Real> weights_vec = get_vector_data_from_xml<Real>(weights_tree);
        AssertThrow(!weights_vec.empty(), ExcEmptyObject());
        AssertThrow(weights_vec.size() == n_weights,ExcDimensionMismatch(weights_vec.size(),n_weights));

        weights[comp_id].resize(dofs_size);
        std::copy(weights_vec.begin(), weights_vec.end(), weights[comp_id].begin());
        //-------------------------------------------------------------------------

    } // end loop over the scalar components
    //-------------------------------------------------------------------------


    typename space_t::EndBehaviourTable end_behaviour(components_map);
    for (const auto comp_id : end_behaviour.get_active_components_id())
    {
        end_behaviour[comp_id][0] = BasisEndBehaviour::interpolatory;
        end_behaviour[comp_id][1] = BasisEndBehaviour::interpolatory;
    }

    auto spline_space = space_t::SpSpace::create(degrees,grid,multiplicities,end_behaviour);
    //---------------------------------------------------------------------------------


    //----------------------------------------
    // building the weight function table --- begin
    using ScalarBSplineSpace = BSplineSpace<dim>;
    using WeightFunc = IgFunction<ScalarBSplineSpace>;

    using ScalarDegreeTable = typename ScalarBSplineSpace::DegreeTable;
    const ScalarDegreeTable scalar_degree_table(degrees[0]);

    auto new_grid = CartesianGrid<dim>::create(*grid);

    using ScalarMultiplicityTable = typename ScalarBSplineSpace::MultiplicityTable;
    const auto scalar_mult_table = shared_ptr<const ScalarMultiplicityTable>(new ScalarMultiplicityTable((*multiplicities)[0]));

    auto scalar_spline_space =
    		ScalarBSplineSpace::create(scalar_degree_table, new_grid,
    				scalar_mult_table,
					typename ScalarBSplineSpace::EndBehaviourTable(filled_array<BasisEndBehaviour,dim>(BasisEndBehaviour::interpolatory)));

    using WeightFuncPtr = shared_ptr<WeightFunc>;
    using WeightFuncPtrTable = typename space_t::template ComponentContainer<WeightFuncPtr>;
    WeightFuncPtrTable w_func_table;
    int comp = 0;
    for (const auto &w_coefs : weights)
        w_func_table[comp++] = WeightFuncPtr(
                                   new WeightFunc(scalar_spline_space,vector<Real>(w_coefs.get_data())));


    //*/
    // building the weight function table --- end
    //----------------------------------------


    auto ref_space = space_t::create(spline_space,w_func_table);

    return ref_space;
}
#endif

IGA_NAMESPACE_CLOSE

#include <igatools/io/reader.inst>


#if 0

#ifdef NURBS
#include <igatools/basis_functions/nurbs_element_accessor.h>
#endif
#include <igatools/basis_functions/bspline_element_accessor.h>

#include <igatools/geometry/cartesian_grid.h>
#include <igatools/geometry/ig_mapping.h>
#include <igatools/utils/vector_tools.h>

#include <boost/algorithm/string.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <array>
#include <set>


using std::array;
using std::set;
using std::make_shared;
using namespace iga::reader_utils;

IGA_NAMESPACE_OPEN

namespace reader_utils
{



}; // of namespace reader_utils.







template <int dim, int codim = 0>
std::shared_ptr< Mapping<dim,codim> >
ig_mapping_reader_version_1_0(const std::string &filename)
{
    Assert(false,ExcNotImplemented());
    AssertThrow(false,ExcNotImplemented());
    return nullptr;

#if 0
    const int dim_phys = dim + codim;

    const auto &xml_tree = get_xml_tree(filename);

    LogStream out ;
    TensorIndex<dim> degree;
    CartesianProductArray<Real,dim> knots_unique_values;
    CartesianProductArray<Size, dim> multiplicities;

    TensorSize<dim> n_ctrl_points_dim;
    array<vector<Real>,dim_phys> control_pts_coords;
    DynamicMultiArray<Real,dim> weights;

    bool is_nurbs_mapping = false;

    for (const auto &patch : xml_tree.get_child("XMLFile.Patch"))
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
            for (const auto &knot : patch.second)
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

                    for (const auto &brk : knot.second)
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

                    for (const auto &mlt : knot.second)
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

                    for (const auto &coord : cp.second)
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

                    for (const auto &wght : cp.second)
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
                    for (auto &w : weights)
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


    using SpSpace = SplineSpace<dim,dim_phis,1>;
    using DegreeTable = typename SpSpace::DegreeTable;
    using MultiplicityTable = typename SpSpace::MultiplicityTable;
    using EndBehaviourTable = typename SpSpace::EndBehaviourTable;

    auto grid = CartesianGrid<dim>::create(knots_unique_values);
    shared_ptr< Mapping<dim,codim> > map;

    if (is_nurbs_mapping)
    {
        Assert(false,ExcNotImplemented());
        AssertThrow(false,ExcNotImplemented());
#if 0

        using space_t = NURBSSpace<dim,dim_phys,1>;
        auto space = space_t::create(
                         grid,
                         StaticMultiArray<Multiplicity<dim>,dim_phys,1>(multiplicities),
                         StaticMultiArray<TensorIndex<dim>,dim_phys,1>(degree),
                         StaticMultiArray<DynamicMultiArray<Real,dim>,dim_phys,1>(weights));

        map = IgMapping<space_t>::create(space,control_pts);
#endif
    }
    else
    {
        using space_t = BSplineSpace<dim,dim_phys,1>;
        auto space = space_t::create(
                         grid,
                         StaticMultiArray<Multiplicity<dim>,dim_phys,1>(multiplicities),
                         StaticMultiArray<TensorIndex<dim>,dim_phys,1>(degree));
//      BSplineSpace (const DegreeTable &deg, std::shared_ptr< GridType > knots, std::shared_ptr< const MultiplicityTable > interior_mult, const EndBehaviourTable &ends)

        map = IgMapping<space_t>::create(space,control_pts);
    }
//    map->print_info(out);

    return map;

#endif

}









IGA_NAMESPACE_CLOSE


#endif
