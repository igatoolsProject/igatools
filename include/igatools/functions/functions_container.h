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

#ifndef __FUNCTIONS_CONTAINER_H
#define __FUNCTIONS_CONTAINER_H

#include <igatools/base/config.h>



#include <igatools/base/identity_function.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/base/function_lib.h>
#include <igatools/base/function_element.h>
#include <igatools/base/sub_function.h>

#include <igatools/base/ig_function.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/nurbs_space.h>
#include <igatools/basis_functions/physical_space.h>

IGA_NAMESPACE_OPEN

/**
 * Type alias for the associative container (a std::map) between a pointer to a
 * Function<dim,codim,range,rank> and a string (e.g. the function's name).
 */
template <int dim,int codim,int range,int rank>
using DictionaryFuncPtrName =
    std::map<std::shared_ptr<Function<dim,codim,range,rank>>,std::string>;


/**
 * Type alias for the heterogeneous associative container (a boost::fusion::map)
 * between a type representing the <tt>rank</tt> index and a collection of functions
 * with the same quadruplet <tt><dim,codim,range,rank</tt>.
 */
template <int dim,int codim,int range>
using DataVaryingRank = boost::fusion::map<
                        boost::fusion::pair< Topology<1>,DictionaryFuncPtrName<dim,codim,range,1> > >;


template <int dim,int codim,int range>
class StuffSameDimAndCodimAndRange
{
public:
    DataVaryingRank<dim,codim,range> data_varying_rank_;

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
            int f_id = 0;
            for (const auto &f : funcs_with_name)
            {
                out.begin_item("Function[" + std::to_string(f_id++) + "] name: " + f.second);
//                f.first->print_info(out);
                out.end_item();
            }
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
        using std::to_string;
        out.begin_item("Mappings num. : " + to_string(maps_and_data_varying_range_.size()));
        int map_id = 0;
        for (const auto &map_and_data_varying_range : maps_and_data_varying_range_)
        {
            const auto &data_varying_range = map_and_data_varying_range.second;
            out.begin_item("Map[" + to_string(map_id++) + "] name: " + data_varying_range.map_name_);

            const auto &map = *map_and_data_varying_range.first;
            map.print_info(out);

            boost::fusion::for_each(data_varying_range.funcs_,
                                    [&](const auto & type_and_data_same_range)
            {
                using Type_Value = typename std::remove_reference<decltype(type_and_data_same_range)>::type;
                using Type = typename Type_Value::first_type;

                out.begin_item("Range : " + to_string(Type::value));
                type_and_data_same_range.second.print_info(out);
                out.end_item();
            } // end lambda function
                                   );
            out.end_item();
        } // end loop maps
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
        ar.template register_type<IdentityFunction<dim,dim+codim>>();
        ar.template register_type<IgFunction<dim,0,dim+codim,1>>();
//        ar.template register_type<IgFunction<dim,0,dim+codim,1>>();
        ar &boost::serialization::make_nvp("maps_and_data_varying_range_",
                                           maps_and_data_varying_range_);
    }
    ///@}
#endif // SERIALIZATION

}; // end class StuffSameDimAndCodim




template<int dim, std::size_t... I>
auto
make_DataVaryingCodim(std::index_sequence<I...>)
{
    return boost::fusion::map<
           boost::fusion::pair<Topology<I>,StuffSameDimAndCodim<dim,I> > ...>(
               boost::fusion::pair<Topology<I>,StuffSameDimAndCodim<dim,I> >() ...);
}

/**
 *
 */
template <int dim>
using DataVaryingCodim = decltype(make_DataVaryingCodim<dim>(std::make_index_sequence<4-dim>()));

//*/
/*
template <int dim>
struct DataVaryingCodim;

template <>
struct DataVaryingCodim<0>
{
    using type = boost::fusion::map<
                 boost::fusion::pair< Topology<0>,StuffSameDimAndCodim<0,0> >,
                 boost::fusion::pair< Topology<1>,StuffSameDimAndCodim<0,1> >,
                 boost::fusion::pair< Topology<2>,StuffSameDimAndCodim<0,2> >,
                 boost::fusion::pair< Topology<3>,StuffSameDimAndCodim<0,3> > >;
};

template <>
struct DataVaryingCodim<1>
{
    using type = boost::fusion::map<
                 boost::fusion::pair< Topology<0>,StuffSameDimAndCodim<1,0> >,
                 boost::fusion::pair< Topology<1>,StuffSameDimAndCodim<1,1> >,
                 boost::fusion::pair< Topology<2>,StuffSameDimAndCodim<1,2> > >;
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
//*/

template <int dim>
class StuffSameDim
{
public:
    DataVaryingCodim<dim> data_varying_codim_;


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

            const std::string tag_name = "data_codim_" + std::to_string(Type::value);
            ar &boost::serialization::make_nvp(tag_name.c_str(),type_and_data_same_codim.second);
        } // end lambda function
                               );
    }
    ///@}
#endif // SERIALIZATION

}; // end StuffSameDim


template<std::size_t... I>
auto
make_DataVaryingDim(std::index_sequence<I...>)
{
    return boost::fusion::map<
           boost::fusion::pair<Topology<I+1>,StuffSameDim<I+1> > ...>(
               boost::fusion::pair<Topology<I+1>,StuffSameDim<I+1> >() ...);
}


/**
 * Alias for the type representing the all the data in kept in the FunctionsContainer class.
 */
template <int dim>
using DataVaryingDim = decltype(make_DataVaryingDim(std::make_index_sequence<dim>()));

/**
 * @brief Heterogeneous associative container between geometry parametrizations
 * (the container's "key") and functions (the "value" associated to the "key").
 *
 * The container is heterogeneous in the sense that can hold data referred to different combinations
 * of quadruplets <tt><dim,codim,range,rank></tt> and it is defined as a "4-index" container,
 * and the data is organized following the quadruplets order <tt><dim,codim,range,rank></tt>.
 * Therefore the first index selects the data identified with <tt>dim</tt> and the second
 * index the data identified with the pair <tt>dim,codim</tt>: at this point we use an associative
 * container (a std::map) having as a key a (shared) pointer to the geometry parametrization
 * and as associated value, an object containing the name given to the geometry and the container
 * holding all functions associated to the geometry.
 *
 * Unfortunately the STL library does not provide associative containers between heterogeneous
 * types, so we used the boost::fusion library (and specifically the boost::fusion::map class)
 * in order to have such kind of container.
 *
 * @ingroup serializable
 *
 * @author M. Martinelli, 2015
 */
class FunctionsContainer
{
public:

    template <int dim>
    const auto &get_data_dim() const
    {
        return boost::fusion::at_key<Topology<dim>>(data_varying_dim_);
    }

    /**
     * Adds a @p map (i.e. a geometry paramtrization) with the given @p map_name to the container.
     *
     * @note In Debug mode, an assertion will be raised if
     * the @p map is already present in the container,
     */
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

    /**
     * Adds the Function @p function to the container and creates an association with the
     * MapFunction @p map. At the end, the @p function is tagged with the string @p func_name.
     *
     * @pre 1) The @p map should be present in the container (i.e. should be inserted with insert_map()).
     * @pre 2) The association <tt>map-function</tt> must not be already established.
     *
     * @note In Debug mode an assertion will be raised if some of the two precondition
     * above are unsatisfied.
     */
    template<int dim, int codim,int range,int rank>
    void insert_function(
        std::shared_ptr<MapFunction<dim,dim+codim>> map,
        std::shared_ptr<Function<dim,codim,range,rank>> function,
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

        Assert(data_same_dim_codim_range_rank.count(function) == 0,
               ExcMessage("Function already added to the container."));
        data_same_dim_codim_range_rank[function] = func_name;
    }

    /**
     * Prints some internal information. Mostly used for testing and debugging purposes.
     */
    void print_info(LogStream &out) const;

private:


    /**
     * All the data in the FunctionsContainer class, organized by the @p dim index.
     */
    DataVaryingDim<3> data_varying_dim_;



#ifdef SERIALIZATION
    /**
     * @name Functions needed for boost::serialization
     * @see <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     */
    ///@{
    friend class boost::serialization::access;

    template<class Archive>
    void
    serialize(Archive &ar, const unsigned int version);
    ///@}
#endif // SERIALIZATION

};





IGA_NAMESPACE_CLOSE


#endif // __FUNCTIONS_CONTAINER_H
