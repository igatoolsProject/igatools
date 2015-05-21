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
 *  Test for SubFunction class
 *  author: pauletti
 *  date: Oct 12, 2014
 */

#include "../tests.h"

#include <igatools/base/identity_function.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/base/function_lib.h>
#include <igatools/base/function_element.h>
#include <igatools/base/sub_function.h>

#include <igatools/base/ig_function.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/nurbs_space.h>
#include <igatools/basis_functions/physical_space.h>


using std::shared_ptr;
using std::static_pointer_cast;


template <int dim,int codim,int range>
class StuffSameDimAndCodimAndRange
{
public:
    template <int rank>
    using FuncPtr = std::shared_ptr<Function<dim,codim,range,rank>>;

    template <int rank>
    using DictionaryFuncsName = std::map<FuncPtr<rank>,std::string>;

    boost::fusion::map<
    boost::fusion::pair< Topology<1>,DictionaryFuncsName<1> > > data_varying_rank_;

    void print_info(LogStream &out) const
    {
        boost::fusion::for_each(data_varying_rank_,
                                [&](const auto & type_and_data_same_rank)
        {
            using Type_Value = typename std::remove_reference<decltype(type_and_data_same_rank)>::type;
            using Type = typename Type_Value::first_type;

            out.begin_item("Rank : " + std::to_string(Type::value));
            const auto &funcs_with_name = type_and_data_same_rank.second;


            out.begin_item("Functions num. : " + std::to_string(funcs_with_name.size()));
            for (const auto &f : funcs_with_name)
                out << "Function name: " << f.second << std::endl;
            out.end_item();

            out.end_item();
        } // end lambda function
                               );
    }; // end print_info()

private:
#ifdef SERIALIZATION
    /**
     * @name Functions needed for boost::serialization
     * @see <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     */
///@{
    friend class boost::serialization::access;

    template<class Archive>
    void
    serialize(Archive &ar, const unsigned int version)
    {
        boost::fusion::for_each(data_varying_rank_,
                                [&](auto & type_and_data_same_rank)
        {
            using Type_Value = typename std::remove_reference<decltype(type_and_data_same_rank)>::type;
            using Type = typename Type_Value::first_type;

            ar.template register_type<IgFunction<dim,codim,range,Type::value>>();

            const std::string tag_name = "funcs_rank_" + std::to_string(Type::value);
            ar &boost::serialization::make_nvp(tag_name.c_str(),type_and_data_same_rank.second);
        } // end lambda function
                               );

    }
///@}
#endif // SERIALIZATION

};


template <int dim, int codim>
using DataVaryingRange =
    boost::fusion::map<
    boost::fusion::pair<Topology<1>,StuffSameDimAndCodimAndRange<dim,codim,1> >,
    boost::fusion::pair<Topology<dim+codim>,StuffSameDimAndCodimAndRange<dim,codim,dim+codim> >
    >;

template <int dim,int codim>
class StuffSameDimAndCodim
{
public:
    using M = std::shared_ptr<MapFunction<dim,dim+codim>>;

    struct DataAssociatedToMap
    {
        std::string map_name_;
        DataVaryingRange<dim,codim> funcs_;

    private:
#ifdef SERIALIZATION
        /**
         * @name Functions needed for boost::serialization
         * @see <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
         */
        ///@{
        friend class boost::serialization::access;

        template<class Archive>
        void
        serialize(Archive &ar, const unsigned int version)
        {
            ar &boost::serialization::make_nvp("map_name_",map_name_);

            boost::fusion::for_each(funcs_,
                                    [&](auto & func)
            {
                using Type_Value = typename std::remove_reference<decltype(func)>::type;
                using Type = typename Type_Value::first_type;

                const std::string tag_name = "funcs_range_" + std::to_string(Type::value);
                ar &boost::serialization::make_nvp(tag_name.c_str(),func.second);
            } // end lambda function
                                   );

        }
        ///@}
#endif // SERIALIZATION

    };

    std::map<M,DataAssociatedToMap> maps_and_data_varying_range_;


    void print_info(LogStream &out) const
    {
        out.begin_item("Mappings num. : " + std::to_string(maps_and_data_varying_range_.size()));
        for (const auto &map_and_data_varying_range : maps_and_data_varying_range_)
        {
            const auto &data_varying_range = map_and_data_varying_range.second;
            out.begin_item("Map name: " + data_varying_range.map_name_);

            boost::fusion::for_each(data_varying_range.funcs_,
                                    [&](const auto & type_and_data_same_range)
            {
                using Type_Value = typename std::remove_reference<decltype(type_and_data_same_range)>::type;
                using Type = typename Type_Value::first_type;

                out.begin_item("Range : " + std::to_string(Type::value));
                type_and_data_same_range.second.print_info(out);
                out.end_item();
            } // end lambda function
                                   );
            out.end_item();
        }
        out.end_item();
    }; // end print_info()


private:

#ifdef SERIALIZATION
    /**
     * @name Functions needed for boost::serialization
     * @see <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     */
    ///@{
    friend class boost::serialization::access;

    template<class Archive>
    void
    serialize(Archive &ar, const unsigned int version)
    {
//          ar.template register_type<IdentityFunction<dim,dim+codim>>();
        ar.template register_type<IgFunction<dim,0,dim+codim,1>>();
//        ar.template register_type<IgFunction<dim,0,dim+codim,1>>();
        ar &boost::serialization::make_nvp("maps_and_data_varying_range_",
                                           maps_and_data_varying_range_);
    }
    ///@}
#endif // SERIALIZATION

}; // end class StuffSameDimAndCodim




template <int dim>
struct DataVaryingCodim;

template <>
struct DataVaryingCodim<0>
{
    using type = boost::fusion::map<
                 boost::fusion::pair< Topology<0>,StuffSameDimAndCodim<0,0> >,
                 boost::fusion::pair< Topology<1>,StuffSameDimAndCodim<0,1> >,
                 boost::fusion::pair< Topology<1>,StuffSameDimAndCodim<0,2> >,
                 boost::fusion::pair< Topology<1>,StuffSameDimAndCodim<0,3> > >;
};

template <>
struct DataVaryingCodim<1>
{
    using type = boost::fusion::map<
                 boost::fusion::pair< Topology<0>,StuffSameDimAndCodim<1,0> >,
                 boost::fusion::pair< Topology<1>,StuffSameDimAndCodim<1,1> >,
                 boost::fusion::pair< Topology<1>,StuffSameDimAndCodim<1,2> > >;
};

template <>
struct DataVaryingCodim<2>
{
    using type = boost::fusion::map<
                 boost::fusion::pair< Topology<0>,StuffSameDimAndCodim<2,0> >,
                 boost::fusion::pair< Topology<1>,StuffSameDimAndCodim<2,1> > >;
};

template <>
struct DataVaryingCodim<3>
{
    using type = boost::fusion::map<
                 boost::fusion::pair< Topology<0>,StuffSameDimAndCodim<3,0> > >;
};


template <int dim>
class StuffSameDim
{
public:
    typename DataVaryingCodim<dim>::type data_varying_codim_;


    void print_info(LogStream &out) const
    {
        boost::fusion::for_each(data_varying_codim_,
                                [&](const auto & type_and_data_same_codim)
        {
            using Type_Value = typename std::remove_reference<decltype(type_and_data_same_codim)>::type;
            using Type = typename Type_Value::first_type;

            out.begin_item("Codim : " + std::to_string(Type::value));
            type_and_data_same_codim.second.print_info(out);
            out.end_item();

        } // end lambda function
                               );
    }

private:

#ifdef SERIALIZATION
    /**
     * @name Functions needed for boost::serialization
     * @see <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     */
    ///@{
    friend class boost::serialization::access;

    template<class Archive>
    void
    serialize(Archive &ar, const unsigned int version)
    {
        boost::fusion::for_each(data_varying_codim_,
                                [&](auto & type_and_data_same_codim)
        {
            using Type_Value = typename std::remove_reference<decltype(type_and_data_same_codim)>::type;
            using Type = typename Type_Value::first_type;

            const string tag_name = "data_codim_" + std::to_string(Type::value);
            ar &boost::serialization::make_nvp(tag_name.c_str(),type_and_data_same_codim.second);
        } // end lambda function
                               );
    }
    ///@}
#endif // SERIALIZATION

}; // end StuffSameDim


class FunctionsContainer
{
public:

    template<int dim, int space_dim>
    void insert_map(std::shared_ptr<MapFunction<dim,space_dim>> map, const std::string &map_name)
    {
        using boost::fusion::at_key;
        auto &data_same_dim = at_key<Topology<dim>>(data_varying_dim_);
        auto &data_same_dim_codim = at_key<Topology<space_dim-dim>>(data_same_dim.data_varying_codim_);

        Assert(data_same_dim_codim.maps_and_data_varying_range_.count(map) == 0,
               ExcMessage("Map already added in the container."));

        data_same_dim_codim.maps_and_data_varying_range_[map].map_name_ = map_name;
    };

    template<int dim, int codim,int range,int rank>
    void insert_function(
        std::shared_ptr<MapFunction<dim,dim+codim>> map,
        std::shared_ptr<Function<dim,codim,range,rank>> func,
        const std::string &func_name)
    {
        using boost::fusion::at_key;
        auto &data_same_dim = at_key<Topology<dim>>(data_varying_dim_);
        auto &data_same_dim_codim = at_key<Topology<codim>>(data_same_dim.data_varying_codim_);

        Assert(data_same_dim_codim.maps_and_data_varying_range_.count(map) == 1,
               ExcMessage("Map not present in the container."));

        auto &data_same_map = data_same_dim_codim.maps_and_data_varying_range_[map];

        auto &data_same_dim_codim_range = at_key<Topology<range>>(data_same_map.funcs_);

        auto &data_same_dim_codim_range_rank = at_key<Topology<rank>>(data_same_dim_codim_range.data_varying_rank_);

        Assert(data_same_dim_codim_range_rank.count(func) == 0,
               ExcMessage("Function already added to the container."));
        data_same_dim_codim_range_rank[func] = func_name;
    }


    void print_info(LogStream &out) const
    {
        boost::fusion::for_each(data_varying_dim_,
                                [&](const auto & type_and_data_same_dim)
        {
            using Type_Value = typename std::remove_reference<decltype(type_and_data_same_dim)>::type;
            using Type = typename Type_Value::first_type;

            out.begin_item("Dim : " + std::to_string(Type::value));
            type_and_data_same_dim.second.print_info(out);
            out.end_item();

        } // end lambda function
                               );
    }

private:








    boost::fusion::map<
    boost::fusion::pair< Topology<1>,StuffSameDim<1> >,
          boost::fusion::pair< Topology<2>,StuffSameDim<2> >,
          boost::fusion::pair< Topology<3>,StuffSameDim<3> >
          > data_varying_dim_;


private:

#ifdef SERIALIZATION
    /**
     * @name Functions needed for boost::serialization
     * @see <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     */
    ///@{
    friend class boost::serialization::access;

    template<class Archive>
    void
    serialize(Archive &ar, const unsigned int version)
    {
        boost::fusion::for_each(data_varying_dim_,
                                [&](auto & type_and_data_same_dim)
        {
            using Type_Value = typename std::remove_reference<decltype(type_and_data_same_dim)>::type;
            using Type = typename Type_Value::first_type;

            const string tag_name = "data_dim_" + std::to_string(Type::value);
            ar &boost::serialization::make_nvp(tag_name.c_str(),type_and_data_same_dim.second);
        } // end lambda function
                               );
    }
    ///@}
#endif // SERIALIZATION

};




void serialize_deserialize(std::shared_ptr<FunctionsContainer> funcs_container)
{
    OUTSTART

    out.begin_item("Original FunctionsContainer:");
    funcs_container->print_info(out);
    out.end_item();

    std::string filename = "functions_container.xml";
    std::string tag_name = "FunctionsContainer";
    {
        // serialize the PhysicalSpace object to an xml file
        std::ofstream xml_ostream(filename);
        OArchive xml_out(xml_ostream);

        xml_out << boost::serialization::make_nvp(tag_name.c_str(),funcs_container);
        xml_ostream.close();
    }

    funcs_container.reset();
    {
        // de-serialize the PhysicalSpace object from an xml file
        std::ifstream xml_istream(filename);
        IArchive xml_in(xml_istream);

        xml_in >> BOOST_SERIALIZATION_NVP(funcs_container);
        xml_istream.close();
    }
    out.begin_item("FunctionsContainer after serialize-deserialize:");
    funcs_container->print_info(out);
    out.end_item();
//*/

    OUTEND
}



template <int dim,int codim,int range>
using Func = Function<dim,codim,range,1>;


void do_test()
{
    int n_elem_per_side = 2;
    auto grid_1 = CartesianGrid<1>::create(n_elem_per_side+1);
    auto grid_2 = CartesianGrid<2>::create(n_elem_per_side+1);
    auto grid_3 = CartesianGrid<3>::create(n_elem_per_side+1);
//    create_fun<2, 0, 2>();

    auto func_identity_1_1 = IdentityFunction<1,1>::create(grid_1);
    auto func_identity_2_2 = IdentityFunction<2,2>::create(grid_2);
    auto func_identity_3_3 = IdentityFunction<3,3>::create(grid_3);



    const int deg = 3;
    auto bsp_space_1_1 = BSplineSpace<1,1,1>::create(deg, grid_1);
    auto bsp_space_2_1 = BSplineSpace<2,1,1>::create(deg, grid_2);
    auto bsp_space_3_1 = BSplineSpace<3,1,1>::create(deg, grid_3);
    auto bsp_space_2_2 = BSplineSpace<2,2,1>::create(deg, grid_2);
    auto bsp_space_3_3 = BSplineSpace<3,3,1>::create(deg, grid_3);
    auto bsp_space_2_3 = BSplineSpace<2,3,1>::create(deg, grid_2);


    Epetra_SerialComm comm;
    auto bsp_coeff_1_1 = EpetraTools::create_vector(
                             EpetraTools::create_map(bsp_space_1_1, "active", comm));
    (*bsp_coeff_1_1)[0] = 1.;


    auto bsp_coeff_2_1 = EpetraTools::create_vector(
                             EpetraTools::create_map(bsp_space_2_1, "active", comm));
    (*bsp_coeff_2_1)[0] = 1.;


    auto bsp_coeff_3_1 = EpetraTools::create_vector(
                             EpetraTools::create_map(bsp_space_3_1, "active", comm));
    (*bsp_coeff_3_1)[0] = 1.;


    auto bsp_coeff_2_2 = EpetraTools::create_vector(
                             EpetraTools::create_map(bsp_space_2_2, "active", comm));
    (*bsp_coeff_2_2)[0] = 1.;


    auto bsp_coeff_3_3 = EpetraTools::create_vector(
                             EpetraTools::create_map(bsp_space_3_3, "active", comm));
    (*bsp_coeff_3_3)[0] = 1.;


    auto bsp_coeff_2_3 = EpetraTools::create_vector(
                             EpetraTools::create_map(bsp_space_2_3, "active", comm));
    (*bsp_coeff_2_3)[0] = 1.;


    auto bsp_func_1_1 = IgFunction<1,0,1,1>::create(bsp_space_1_1, bsp_coeff_1_1);
    auto bsp_func_2_1 = IgFunction<2,0,1,1>::create(bsp_space_2_1, bsp_coeff_2_1);
    auto bsp_func_3_1 = IgFunction<3,0,1,1>::create(bsp_space_3_1, bsp_coeff_3_1);
    auto bsp_func_2_2 = IgFunction<2,0,2,1>::create(bsp_space_2_2, bsp_coeff_2_2);
    auto bsp_func_3_3 = IgFunction<3,0,3,1>::create(bsp_space_3_3, bsp_coeff_3_3);
    auto bsp_func_2_3 = IgFunction<2,0,3,1>::create(bsp_space_2_3, bsp_coeff_2_3);



    auto phys_space_1_1_1_0 =
        PhysicalSpace<1,1,1,0,Transformation::h_grad>::create(
            bsp_space_1_1,
            bsp_func_1_1);

    auto phys_space_2_1_1_0 =
        PhysicalSpace<2,1,1,0,Transformation::h_grad>::create(
            bsp_space_2_1,
            bsp_func_2_2);

    auto phys_space_3_1_1_0 =
        PhysicalSpace<3,1,1,0,Transformation::h_grad>::create(
            bsp_space_3_1,
            bsp_func_3_3);

    auto phys_space_2_2_1_0 =
        PhysicalSpace<2,2,1,0,Transformation::h_grad>::create(
            bsp_space_2_2,
            bsp_func_2_2);

    auto phys_space_3_3_1_0 =
        PhysicalSpace<3,3,1,0,Transformation::h_grad>::create(
            bsp_space_3_3,
            bsp_func_3_3);

    auto phys_space_2_1_1_1 =
        PhysicalSpace<2,1,1,1,Transformation::h_grad>::create(
            bsp_space_2_1,
            bsp_func_2_3);

    auto phys_space_2_3_1_1 =
        PhysicalSpace<2,3,1,1,Transformation::h_grad>::create(
            bsp_space_2_3,
            bsp_func_2_3);

    auto phys_coeff_1_1_1_0 = EpetraTools::create_vector(
                                  EpetraTools::create_map(phys_space_1_1_1_0, "active", comm));
    (*phys_coeff_1_1_1_0)[0] = 1.;


    auto phys_coeff_2_1_1_0 = EpetraTools::create_vector(
                                  EpetraTools::create_map(phys_space_2_1_1_0, "active", comm));
    (*phys_coeff_2_1_1_0)[0] = 1.;


    auto phys_coeff_3_1_1_0 = EpetraTools::create_vector(
                                  EpetraTools::create_map(phys_space_3_1_1_0, "active", comm));
    (*phys_coeff_3_1_1_0)[0] = 1.;


    auto phys_coeff_2_2_1_0 = EpetraTools::create_vector(
                                  EpetraTools::create_map(phys_space_2_2_1_0, "active", comm));
    (*phys_coeff_2_2_1_0)[0] = 1.;


    auto phys_coeff_3_3_1_0 = EpetraTools::create_vector(
                                  EpetraTools::create_map(phys_space_3_3_1_0, "active", comm));
    (*phys_coeff_3_3_1_0)[0] = 1.;


    auto phys_coeff_2_1_1_1 = EpetraTools::create_vector(
                                  EpetraTools::create_map(phys_space_2_1_1_1, "active", comm));
    (*phys_coeff_2_1_1_1)[0] = 1.;


    auto phys_coeff_2_3_1_1 = EpetraTools::create_vector(
                                  EpetraTools::create_map(phys_space_2_3_1_1, "active", comm));
    (*phys_coeff_2_3_1_1)[0] = 1.;

    auto phys_func_1_1_1_0 = IgFunction<1,0,1,1>::create(phys_space_1_1_1_0,phys_coeff_1_1_1_0);
    auto phys_func_2_1_1_0 = IgFunction<2,0,1,1>::create(phys_space_2_1_1_0,phys_coeff_2_1_1_0);
    auto phys_func_3_1_1_0 = IgFunction<3,0,1,1>::create(phys_space_3_1_1_0,phys_coeff_3_1_1_0);
    auto phys_func_2_2_1_0 = IgFunction<2,0,2,1>::create(phys_space_2_2_1_0,phys_coeff_2_2_1_0);
    auto phys_func_3_3_1_0 = IgFunction<3,0,3,1>::create(phys_space_3_3_1_0,phys_coeff_3_3_1_0);
    auto phys_func_2_1_1_1 = IgFunction<2,1,1,1>::create(phys_space_2_1_1_1,phys_coeff_2_1_1_1);
    auto phys_func_2_3_1_1 = IgFunction<2,1,3,1>::create(phys_space_2_3_1_1,phys_coeff_2_3_1_1);


    auto funcs_container = std::make_shared<FunctionsContainer>();

    funcs_container->insert_map(
        phys_func_1_1_1_0->get_ig_space()->get_map_func(),
        "map_1_1_1_0");

    funcs_container->insert_map(
        phys_func_2_1_1_0->get_ig_space()->get_map_func(),
        "map_2_1_1_0");

    funcs_container->insert_map(
        phys_func_3_1_1_0->get_ig_space()->get_map_func(),
        "map_3_1_1_0");

    funcs_container->insert_map(
        phys_func_2_2_1_0->get_ig_space()->get_map_func(),
        "map_2_2_1_0");

    funcs_container->insert_map(
        phys_func_3_3_1_0->get_ig_space()->get_map_func(),
        "map_3_3_1_0");

    funcs_container->insert_map(
        phys_func_2_1_1_1->get_ig_space()->get_map_func(),
        "map_2_1_1_1");

    funcs_container->insert_map(
        phys_func_2_3_1_1->get_ig_space()->get_map_func(),
        "map_2_3_1_1");

//    funcs_container->insert_map(func_identity_1_1,"map_identity_1_1");
//    funcs_container->insert_map(func_identity_2_2,"map_identity_2_2");
//    funcs_container->insert_map(func_identity_3_3,"map_identity_3_3");

    funcs_container->insert_function(
        phys_func_1_1_1_0->get_ig_space()->get_map_func(),
        static_pointer_cast<Func<1,0,1>>(phys_func_1_1_1_0),
        "phys_func_1_1_1_0");

    funcs_container->insert_function(
        phys_func_2_1_1_0->get_ig_space()->get_map_func(),
        static_pointer_cast<Func<2,0,1>>(phys_func_2_1_1_0),
        "phys_func_2_1_1_0");

    funcs_container->insert_function(
        phys_func_3_1_1_0->get_ig_space()->get_map_func(),
        static_pointer_cast<Func<3,0,1>>(phys_func_3_1_1_0),
        "phys_func_3_1_1_0");

    funcs_container->insert_function(
        phys_func_2_2_1_0->get_ig_space()->get_map_func(),
        static_pointer_cast<Func<2,0,2>>(phys_func_2_2_1_0),
        "phys_func_2_2_1_0");

    funcs_container->insert_function(
        phys_func_3_3_1_0->get_ig_space()->get_map_func(),
        static_pointer_cast<Func<3,0,3>>(phys_func_3_3_1_0),
        "phys_func_3_3_1_0");

    funcs_container->insert_function(
        phys_func_2_1_1_1->get_ig_space()->get_map_func(),
        static_pointer_cast<Func<2,1,1>>(phys_func_2_1_1_1),
        "phys_func_2_1_1_1");

    funcs_container->insert_function(
        phys_func_2_3_1_1->get_ig_space()->get_map_func(),
        static_pointer_cast<Func<2,1,3>>(phys_func_2_3_1_1),
        "phys_func_2_3_1_1");

//    funcs_container->print_info(out);

    serialize_deserialize(funcs_container);
}



int main()
{
    do_test();

    return 0;
}
