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



#include <igatools/base/tuple_utils.h>
//#include <igatools/functions/identity_function.h>
//#include <igatools/base/quadrature_lib.h>
#include <igatools/functions/function_lib.h>
#include <igatools/functions/function_element.h>
//#include <igatools/functions/sub_function.h>
#include <igatools/functions/ig_function.h>

#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/nurbs.h>
#include <igatools/basis_functions/physical_space.h>

#include <unordered_map>

IGA_NAMESPACE_OPEN

/**
 * Type alias for a shared pointer to an object representing the domain
 * \f$ \mathbf{F} \colon \mathbb{R}^{\text{dim}} \to \mathbb{R}^{\text{dim}+\text{codim}} \f$
 *
 */
template<int dim, int codim>
using MappingPtr = SharedPtrConstnessHandler<Domain<dim,codim>>;

#if 0
template <int dim,int codim,int range,int rank>
struct PairFuncPtrName
{
  std::shared_ptr<Function<dim,codim,range,rank>> func_ptr_;

  std::string name_;

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
    ar &boost::serialization::make_nvp("func_ptr_",func_ptr_);
    ar &boost::serialization::make_nvp("name_",name_);
  }
  ///@}
#endif // SERIALIZATION

};
#endif

template <int dim,int codim,int range,int rank>
using FuncPtr = SharedPtrConstnessHandler<Function<dim,codim,range,rank> >;

/**
 * Type alias for the associative container (a std::map) between an object_id and
 * the a shared pointer pointer to a
 * Function<dim,codim,range,rank>.
 */
template <int dim,int codim,int range,int rank>
using DictionaryFuncPtr =
  std::map<Index,FuncPtr<dim,codim,range,rank> >;












/**
 * @brief Heterogeneous associative container between domain
 * (the container's "key") and functions (the "value" associated to the "key").
 *
 * The association is of the type <em>one-to-many</em> in the sense that
 * for each domain (i.e. for each "key") there can be associated any number of
 * functions (including zero, meaning no-function associated to the geometry).
 *
 * The container is heterogeneous in the sense that can hold data referred to different combinations
 * of quadruplets <tt><dim,codim,range,rank></tt> and it is defined as a "4-index" container,
 * The data is organized following the quadruplets order <tt><dim,codim,range,rank></tt>:
 * - the first index selects the data identified with <tt>dim</tt>
 * - the second index selects the data identified with the pair <tt>dim,codim</tt>: at this point we use an associative
 * container (a std::map) having as <tt>key</tt> a (shared) pointer to a geometry parametrization
 * \f$ \mathbf{F}_i \colon \mathbb{R}^{\text{dim}} \to \mathbb{R}^{\text{dim}+\text{codim}} \f$
 * and as associated <tt>value</tt>, an object containing the name given to the geometry and the container
 * holding all functions associated to the geometry with <tt>range</tt> and <tt>rank</tt> values compatible with
 * the <tt>dim,codim</tt> pair.
 *
 * Unfortunately the STL library does not provide associative containers between heterogeneous
 * types, so we used the <a href="http://www.boost.org/libs/fusion/">boost::fusion library</a>
 * (and specifically the boost::fusion::map class)
 * in order to have such kind of container.
 *
 * An example of the use of this class is given by the test tests/functions/functions_container_01.cpp
 *
 * @ingroup serializable
 *
 * @author M. Martinelli, 2015
 */
class FunctionsContainer
{

  /**
   * @brief Class used by FunctionsContainer in order to store all the data identified by the
   * index <tt>dim</tt>.
   *
   * @ingroup serializable
   * @author M. Martinelli, 2015
   */
  template <int dim>
  class FunctionsContainerDataSameDim
  {

    /**
     * @brief Class used by FunctionsContainerSameDim in order to store all the data identified by the
     * index <tt>codim</tt>.
     *
     * @ingroup serializable
     * @author M. Martinelli, 2015
     */
    template <int codim>
    class FunctionsContainerDataSameDimAndCodim
    {
      /**
       *
       * @serializable
       * @author M. Martinelli, 2015
       */
      template <int range>
      class FunctionsContainerDataSameDimAndCodimAndRange
      {
      public:
        /**
         * Type alias for the heterogeneous associative container (a boost::fusion::map)
         * between a type representing the <tt>rank</tt> index and a collection of functions
         * with the same quadruplet <tt><dim,codim,range,rank></tt>.
         *
         * @note Currently the type is defined to contain only data only for <tt>rank == 1</tt>.
         */
        using DataVaryingRank = boost::fusion::map<
                                boost::fusion::pair< Topology<1>,DictionaryFuncPtr<dim,codim,range,1> > >;

        /**
         * Returns a const-reference to the data identified by the index @p rank.
         */
        template <int rank>
        const DictionaryFuncPtr<dim,codim,range,rank> &get_data_rank() const
        {
          return boost::fusion::at_key<Topology<rank>>(data_varying_rank_);
        }

        /**
         * Returns a reference to the data identified by the index @p rank.
         */
        template <int rank>
        DictionaryFuncPtr<dim,codim,range,rank> &get_data_rank()
        {
          return boost::fusion::at_key<Topology<rank>>(data_varying_rank_);
        }


        /**
         * Prints some internal information. Mostly used for testing and debugging purposes.
         */
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
              out.begin_item("Function[" + std::to_string(f_id++) + "]" +
                             "   ID: " + std::to_string(f.first) +
                             "   name: " + f.second->get_name());
              //                f.first->print_info(out);
              out.end_item();
            }
            out.end_item();

            out.end_item();
          } // end lambda function
                                 );
        }; // end print_info()

      private:

        DataVaryingRank data_varying_rank_;

#ifdef SERIALIZATION
        /**
         * @name Functions needed for the serialization
         */
        ///@{
        friend class cereal::access;

        template<class Archive>
        void
        serialize(Archive &ar)
        {
          boost::fusion::for_each(data_varying_rank_,
                                  [&](auto & type_and_data_same_rank)
          {
            using Type_Value = typename std::remove_reference<decltype(type_and_data_same_rank)>::type;
            using Type = typename Type_Value::first_type;

//            ar.template register_type<IgFunction<dim,codim,range,Type::value>>();
//            using F = Function<dim,codim,range,Type::value>;
            const std::string tag_name = "funcs_rank_" + std::to_string(Type::value);
            ar &make_nvp(tag_name.c_str(),type_and_data_same_rank.second);
          } // end lambda function
                                 );

        }
        ///@}
#endif // SERIALIZATION

      }; // class FunctionsContainerDataSameDimAndCodimAndRange

    public:
      using DataVaryingRange =
        boost::fusion::map<
        boost::fusion::pair<Topology<1>,FunctionsContainerDataSameDimAndCodimAndRange<1> >,
        boost::fusion::pair<Topology<dim+codim>,FunctionsContainerDataSameDimAndCodimAndRange<dim+codim> >
        >;

      /**
       * Returns true if the Domain @p domain is already present in the container.
       */
      bool is_domain_present(const Domain<dim,codim> &domain) const
      {
        return (domains_and_data_varying_range_.count(domain.get_object_id()) == 1)?true:false;
      }

      /**
       * Return a const-reference to the data associated with the Domain @p domain.
       *
       * @pre The Domain @p domain MUST be already present in the container,
       * otherwise in Debug mode an assertion will be raised.
       */
      const auto &get_domain_data(const Domain<dim,codim> &domain) const
      {
        Assert(this->is_domain_present(domain),
               ExcMessage("Domain not present in the container."));
        return domains_and_data_varying_range_.at(domain.get_object_id());
      }

      /**
       * Return a reference to the data associated with the Domain @p domain.
       *
       * @note If the Domain @p domain MUST is already present in the container
       * the reference to the associated data will by returned,
       * otherwise the associated data will be created its reference returned.
       */
      auto &get_domain_data(const Domain<dim,codim> &domain)
      {
        return domains_and_data_varying_range_[domain.get_object_id()];
      }



      /**
       * Returns SafeSTLVector containing the shared pointers to the geometry parametrizations
       * \f$ \mathbf{F}_i \colon \mathbb{R}^{\text{dim}} \to \mathbb{R}^{\text{dim}+\text{codim}} \f$
       */
      SafeSTLVector<MappingPtr<dim,codim> >
      get_all_domains() const
      {
        SafeSTLVector<MappingPtr<dim,codim>> domains;
        for (const auto &m : domains_and_data_varying_range_)
          domains.emplace_back(m.second.get_domain());

        return domains;
      }


      /**
       * @brief Class used to store the data associated to a domain (i.e. a geometry parametrization)
       * \f$\mathbf{F} \colon \mathbb{R}^{\text{dim}} \to \mathbb{R}^{\text{dim}+\text{codim}} \f$
       *
       * The stored data are:
       * - a shared pointer, pointing to the domain function;
       * - the functions associated to the domain. The functions are stored using two
       *   (nested) boost::fusion::map containers, one for the index <tt>range</tt> and the other
       *   for the index <tt>rank</tt>.
       *
       * @ingroup serializable
       */
      class FunctionsSameDomain
      {
      public:
        void set_domain(const MappingPtr<dim,codim> &domain)
        {
          domain_ = domain;
        }

        const MappingPtr<dim,codim> &get_domain() const
        {
          return domain_;
        }
        //*/
        /*
                void set_domain_name(const std::string map_name)
                {
                  domain_->set_name(map_name);
                }
        //*/
        const std::string &get_domain_name() const
        {
          return domain_->get_name();
        }

        template <int range>
        const auto &get_funcs_range() const
        {
          return boost::fusion::at_key<Topology<range>>(funcs_);
        }

        template <int range>
        auto &get_funcs_range()
        {
          return boost::fusion::at_key<Topology<range>>(funcs_);
        }

        /**
         * Prints some internal information. Mostly used for testing and debugging purposes.
         */
        void print_info(LogStream &out) const
        {
          using std::to_string;
          boost::fusion::for_each(funcs_,
                                  [&](const auto & type_and_data_same_range)
          {
            using Type_Value = typename std::remove_reference<decltype(type_and_data_same_range)>::type;
            using Type = typename Type_Value::first_type;

            out.begin_item("Range : " + to_string(Type::value));
            type_and_data_same_range.second.print_info(out);
            out.end_item();
          } // end lambda function
                                 );
        }; // end print_info()

      private:

        MappingPtr<dim,codim> domain_;

        DataVaryingRange funcs_;

#ifdef SERIALIZATION
        /**
         * @name Functions needed for the serialization
         */
        ///@{
        friend class cereal::access;

        template<class Archive>
        void
        serialize(Archive &ar)
        {
          ar &make_nvp("domain_",domain_);

          boost::fusion::for_each(funcs_,
                                  [&](auto & func)
          {
            using Type_Value = typename std::remove_reference<decltype(func)>::type;
            using Type = typename Type_Value::first_type;

            const std::string tag_name = "funcs_range_" + std::to_string(Type::value);
            ar &make_nvp(tag_name.c_str(),func.second);
          } // end lambda function
                                 );
        }
        ///@}
#endif // SERIALIZATION
      };



      /**
       * Prints some internal information. Mostly used for testing and debugging purposes.
       */
      void print_info(LogStream &out) const
      {
        using std::to_string;
        out.begin_item("Domains num. : " + to_string(domains_and_data_varying_range_.size()));
        int map_id = 0;
        for (const auto &domain_and_functions : domains_and_data_varying_range_)
        {
          const auto &data_varying_range = domain_and_functions.second;
          out.begin_item("Domain[" + to_string(map_id++) + "]" +
                         "   ID: " + std::to_string(domain_and_functions.first) +
                         "   name: " + data_varying_range.get_domain_name());

          const auto &domain = *data_varying_range.get_domain();
          domain.print_info(out);

          data_varying_range.print_info(out);

          out.end_item();
        } // end loop maps
        out.end_item();
      }; // end print_info()


    private:

      std::map<Index,FunctionsSameDomain> domains_and_data_varying_range_;

#ifdef SERIALIZATION
      /**
       * @name Functions needed for the serialization
       */
      ///@{
      friend class cereal::access;

      template<class Archive>
      void
      serialize(Archive &ar)
      {
        ar &make_nvp("domains_and_data_varying_range_",
                     domains_and_data_varying_range_);
      }
      ///@}
#endif // SERIALIZATION

    }; // end class FunctionsContainerDataSameDimAndCodim


  public:
    /**
     * Returns a const-reference to the internal data.
     *
     * @note The returned object refers to different values of the <tt>codim</tt> parameter.
     */
    const auto &get_data() const
    {
      return data_varying_codim_;
    }


    /**
     * Returns a const-reference to the data identified by the index pair <tt><dim,codim></tt>.
     *
     * @note The returned container holds the geometry parametrizations
     * \f$ \mathbf{F}_i \colon \mathbb{R}^{\text{dim}} \to \mathbb{R}^{\text{dim}+\text{codim}}\f$.
     */
    template <int codim>
    const auto &get_data_codim() const
    {
      return boost::fusion::at_key<Topology<codim>>(data_varying_codim_);
    }

    /**
     * Returns a const-reference to the data identified by the index pair <tt><dim,codim></tt>.
     *
     * @note The returned container holds the geometry parametrizations
     * \f$ \mathbf{F}_i \colon \mathbb{R}^{\text{dim}} \to \mathbb{R}^{\text{dim}+\text{codim}}\f$.
     */
    template <int codim>
    auto &get_data_codim()
    {
      return boost::fusion::at_key<Topology<codim>>(data_varying_codim_);
    }

    /**
     * Prints some internal information. Mostly used for testing and debugging purposes.
     */
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

    /**
     * Data for each <tt>codim</tt> index.
     * The valid <tt>codim</tt> indices run from <tt>0</tt> to <tt>3-dim</tt>, so
     * we have the following table:
     * dim | codim
     * ----|-------
     *  1  | 0, 1, 2
     *  2  | 0, 1
     *  3  | 0
     */
    DataVaryingId<FunctionsContainerDataSameDimAndCodim,0,4-dim> data_varying_codim_;


#ifdef SERIALIZATION
    /**
     * @name Functions needed for the serialization
     */
    ///@{
    friend class cereal::access;

    template<class Archive>
    void
    serialize(Archive &ar)
    {
      boost::fusion::for_each(data_varying_codim_,
                              [&](auto & type_and_data_same_codim)
      {
        using Type_Value = typename std::remove_reference<decltype(type_and_data_same_codim)>::type;
        using Type = typename Type_Value::first_type;

        const std::string tag_name = "data_codim_" + std::to_string(Type::value);
        ar &make_nvp(tag_name.c_str(),type_and_data_same_codim.second);
      } // end lambda function
                             );
    }
    ///@}
#endif // SERIALIZATION
  }; // end FunctionsContainerDataSameDim



public:

  /**
   * Returns a const-reference to the internal data.
   *
   * @note The returned object refers to different values of the <tt>dim</tt> parameter.
   */
  const auto &get_data() const
  {
    return data_varying_dim_;
  }


  /**
   * Returns a const-reference to the data identified by the index @p dim.
   */
  template <int dim>
  const auto &get_data_dim() const
  {
    return boost::fusion::at_key<Topology<dim>>(data_varying_dim_);
  }

  /**
   * Returns a reference to the data identified by the index @p dim.
   */
  template <int dim>
  auto &get_data_dim()
  {
    return boost::fusion::at_key<Topology<dim>>(data_varying_dim_);
  }

  /**
   * Returns a const-reference to the data identified by the index pair <tt><dim,codim></tt>.
   *
   * @note The returned container holds the geometry parametrizations
   * \f$ \mathbf{F}_i \colon \mathbb{R}^{\text{dim}} \to \mathbb{R}^{\text{dim}+\text{codim}}\f$.
   */
  template <int dim,int codim>
  const auto &get_data_dim_codim() const
  {
    return this->template get_data_dim<dim>().template get_data_codim<codim>();
  }

  /**
   * Returns a reference to the data identified by the index pair <tt><dim,codim></tt>.
   *
   * @note The returned container holds the geometry parametrizations
   * \f$ \mathbf{F}_i \colon \mathbb{R}^{\text{dim}} \to \mathbb{R}^{\text{dim}+\text{codim}}\f$.
   */
  template <int dim,int codim>
  auto &get_data_dim_codim()
  {
    return this->template get_data_dim<dim>().template get_data_codim<codim>();
  }


private:
  /**
   * Adds a Domain @p domain to the container.
   *
   * @note In Debug mode, an assertion will be raised if
   * the @p map is already present in the container,
   */
  template<int dim, int codim>
  void insert_domain_impl(
    const MappingPtr<dim,codim> &domain)
  {
    auto &data_same_dim_codim = this->template get_data_dim_codim<dim,codim>();

    const auto &dom = *(domain.get_ptr_const_data());
    Assert(!data_same_dim_codim.is_domain_present(dom),
           ExcMessage("Domain already present in the container."));

    auto &data_same_domain = data_same_dim_codim.get_domain_data(dom);
    data_same_domain.set_domain(domain);
  };

public:
  /**
   * Adds a Domain @p domain to the container.
   *
   * @note In Debug mode, an assertion will be raised if
   * the @p map is already present in the container,
   */
  template<int dim, int codim>
  void insert_domain(const std::shared_ptr<Domain<dim,codim>> &domain)
  {
    this->insert_domain_impl(SharedPtrConstnessHandler<Domain<dim,codim>>(domain));
  }

  /**
   * Adds a Domain @p domain to the container.
   *
   * @note In Debug mode, an assertion will be raised if
   * the @p map is already present in the container,
   */
  template<int dim, int codim>
  void insert_domain(const std::shared_ptr<const Domain<dim,codim>> &domain)
  {
    this->insert_domain_impl(SharedPtrConstnessHandler<Domain<dim,codim>>(domain));
  }

private:
  template<int dim, int codim,int range,int rank>
  void insert_function_impl(
    const FuncPtr<dim,codim,range,rank> &function)
  {
    using boost::fusion::at_key;

    auto &data_same_dim_codim = this->template get_data_dim_codim<dim,codim>();

    const auto &domain = *(function->get_domain());
    Assert(data_same_dim_codim.is_domain_present(domain),
           ExcMessage("Domain not present in the container."));

    auto &data_same_map = data_same_dim_codim.get_domain_data(domain);

    auto &data_same_dim_codim_range = data_same_map.template get_funcs_range<range>();

    auto &data_same_dim_codim_range_rank =
      data_same_dim_codim_range.template get_data_rank<rank>();

    const int func_id = function->get_object_id();
    Assert(data_same_dim_codim_range_rank.count(func_id) == 0,
           ExcMessage("Function already added to the container."));
    data_same_dim_codim_range_rank[func_id] = function;
  }

public:
  /**
   * @brief Adds the Function @p function to the container and creates an association with its Domain.
   *
   * @pre 1) The Domain used to define the Function @p function must be present in the container
   * (i.e. should be inserted with insert_domain()).
   * @pre 2) The association <tt>domain-function</tt> must not be already established
   * (i.e. the Function used as input parameter is already added).
   *
   * @note In Debug mode an assertion will be raised if some of the two precondition
   * above are unsatisfied.
   */
  template<int dim, int codim,int range,int rank>
  void insert_function(
    const std::shared_ptr<Function<dim,codim,range,rank>> &function)
  {
    this->insert_function_impl(
      SharedPtrConstnessHandler<Function<dim,codim,range,rank>>(function));
  }

  /**
   * @brief Adds the Function @p function to the container and creates an association with its Domain.
   *
   * @pre 1) The Domain used to define the Function @p function must be present in the container
   * (i.e. should be inserted with insert_domain()).
   * @pre 2) The association <tt>domain-function</tt> must not be already established
   * (i.e. the Function used as input parameter is already added).
   *
   * @note In Debug mode an assertion will be raised if some of the two precondition
   * above are unsatisfied.
   */
  template<int dim, int codim,int range,int rank>
  void insert_function(
    const std::shared_ptr<const Function<dim,codim,range,rank>> &function)
  {
    this->insert_function_impl(
      SharedPtrConstnessHandler<Function<dim,codim,range,rank>>(function));
  }

  /**
   * Returns a SafeSTLVector containing the domains defined by
   * \f$ \mathbf{F}_i \colon \mathbb{R}^{\text{dim}} \to \mathbb{R}^{\text{dim}+\text{codim}} \f$
   *
   * @code{.cpp}
     FunctionsContainer funcs_container;
     ... // populating the funcs_container with some domains;

     // here we retrieve all the domains in the funcs_container object, with dimension 2 and co-dimension 1 (i.e. surfaces in 3D space)
     const auto all_domains_dim_2_codim_1 = funcs_container.template get_domains_dim_codim<2,1>();
     @endcode
   */
  template <int dim,int codim>
  SafeSTLVector< MappingPtr<dim,codim> >
  get_domains_dim_codim() const
  {
    return this->template get_data_dim_codim<dim,codim>().get_all_domains();
  }


  /**
   * Returns a const reference to the object associated to the Domain @p domain
   * used as input argument.
   */
  template <int dim,int codim>
  const auto &
  get_domain_data(const Domain<dim,codim> &domain) const
  {
    return this->template get_data_dim_codim<dim,codim>().get_domain_data(domain);
  }


  /**
   * Returns all the functions of the type Function<dim,codim,range,rank> associated
   * to the given Domain @p domain.
   */
  template <int dim,int codim,int range,int rank>
  const DictionaryFuncPtr<dim,codim,range,rank> &
  get_functions_with_same_domain(const Domain<dim,codim> &domain) const
  {
    return this->template get_domain_data<dim,codim>(domain).
    template get_funcs_range<range>().
    template get_data_rank<rank>();
  }



  /**
   * Prints some internal information. Mostly used for testing and debugging purposes.
   */
  void print_info(LogStream &out) const;

private:

  /**
   * All the data in the FunctionsContainer class, organized by the @p dim index
   * (starting from 1 to 3)
   */
  DataVaryingId<FunctionsContainerDataSameDim,1,3> data_varying_dim_;


#ifdef SERIALIZATION
  /**
   * @name Functions needed for the serialization
   */
  ///@{
  friend class cereal::access;

  template<class Archive>
  void
  serialize(Archive &ar)
  {
    boost::fusion::for_each(data_varying_dim_,
                            [&](auto & type_and_data_same_dim)
    {
      using Type_Value = typename std::remove_reference<decltype(type_and_data_same_dim)>::type;
      using Type = typename Type_Value::first_type;

      const std::string tag_name = "data_dim_" + std::to_string(Type::value);
      ar &make_nvp(tag_name.c_str(),type_and_data_same_dim.second);
    } // end lambda function
                           );
  }
  ///@}
#endif

};





IGA_NAMESPACE_CLOSE


#endif // __FUNCTIONS_CONTAINER_H
