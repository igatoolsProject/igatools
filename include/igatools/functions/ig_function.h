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

#ifndef __IG_FUNCTION_H
#define __IG_FUNCTION_H

#include <igatools/base/value_types.h>
#include <igatools/functions/function.h>
#include <igatools/functions/ig_coefficients.h>

#include <igatools/basis_functions/spline_space.h>
#include <igatools/linear_algebra/epetra_vector.h>
//#include <igatools/basis_functions/bspline_space.h>
//#include <igatools/basis_functions/nurbs_space.h>


#include <boost/fusion/include/filter_if.hpp>
//#include <boost/fusion/include/iterator.hpp>
#include <boost/fusion/include/tag_of.hpp>
#include <boost/fusion/include/key_of.hpp>
#include <boost/mpl/not_equal_to.hpp>
#include <boost/fusion/include/begin.hpp>
IGA_NAMESPACE_OPEN


//template <int,int,int>
//class ReferenceSpace;
template <int,int,int,int,Transformation>
class Space;

template <int,int,int>
class BSplineSpace;

template <int,int,int>
class NURBSSpace;

//template <int,int,int,int,class>
//class PhysicalSpace;


template <int,int,int,int,Transformation>
class SpaceElementHandler;

template <int,int,int>
class BSplineElementHandler;

template <int,int,int>
class NURBSElementHandler;

template <int,int,int,int,Transformation>
class PhysSpaceElementHandler;

template <int,int,int,int,Transformation>
class SpaceElement;

template <int,int,int>
class BSplineElement;

template <int,int,int>
class NURBSElement;

template <int,int,int,int,Transformation>
class PhysicalSpaceElement;


/**
 *
 * @ingroup serializable
 */
template<int dim,int codim,int range,int rank>
class IgFunction :
    public Function<dim,codim,range,rank>
{
public:

    using CoeffType = IgCoefficients;


private:
    using base_t = Function<dim,codim,range,rank>;
    using parent_t = Function<dim,codim,range,rank>;
    using self_t = IgFunction<dim,codim,range,rank>;
    using Sp = Space<dim,codim,range,rank,Transformation::h_grad>;

public:
    //TODO (pauletti, Mar 23, 2015): should we make this private?
    IgFunction(std::shared_ptr<const Sp> space,
               std::shared_ptr<const EpetraTools::Vector> coeff,
               const std::string &property = DofProperties::active);

    IgFunction(std::shared_ptr<const Sp> space,
               const IgCoefficients &coeff,
               const std::string &property = DofProperties::active);


    IgFunction(const self_t &);

    virtual ~IgFunction() = default;

    using typename parent_t::topology_variant;
    using typename parent_t::eval_pts_variant;
    using typename parent_t::Point;
    using typename parent_t::Value;
    using typename parent_t::Gradient;
    using typename parent_t::ElementIterator;
    using typename parent_t::ElementAccessor;
    template <int order>
    using Derivative = typename parent_t::template Derivative<order>;



public:
    static std::shared_ptr<self_t>
    create(std::shared_ptr<const Sp> space,
           std::shared_ptr<const EpetraTools::Vector> coeff,
           const std::string &property = DofProperties::active);

    static std::shared_ptr<self_t>
    create(std::shared_ptr<const Sp> space,
           const IgCoefficients &coeff,
           const std::string &property = DofProperties::active);


    std::shared_ptr<base_t> clone() const override final
    {
        return std::make_shared<self_t>(self_t(*this));
    }

#if 0
    void reset(const ValueFlags &flag, const eval_pts_variant &eval_pts) override;

    void reset_selected_elements(const ValueFlags &flag,
                                 const eval_pts_variant &eval_pts,
                                 const SafeSTLVector<Index> &elements_flat_id);

    void init_cache(ElementAccessor &elem, const topology_variant &k) const override;

    void fill_cache(ElementAccessor &elem, const topology_variant &k, const int j) const override;
#endif

    std::shared_ptr<const Sp> get_ig_space() const;

    const CoeffType &get_coefficients() const;

    const std::string &get_property() const;

    self_t &operator +=(const self_t &fun);

    virtual void print_info(LogStream &out) const override final;

protected:
    /**
     * Default constructor. It does nothing but it is needed for the
     * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     * mechanism.
     */
    IgFunction() = default;

private:

    std::shared_ptr<const Sp> space_;

    CoeffType coeff_;

    const std::string property_;

    using SpaceElem = SpaceElement<dim,codim,range,rank,Transformation::h_grad>;
    GridIterator<SpaceElem> space_elem_;

    using SpaceElemHandler = SpaceElementHandler<dim,codim,range,rank,Transformation::h_grad>;
    std::shared_ptr<SpaceElemHandler> space_elem_handler_;

private:
#if 0
    struct ResetDispatcher : boost::static_visitor<void>
    {
        ResetDispatcher(const ValueFlags flag_in,
                        const SafeSTLVector<Index> &elements_flat_id,
                        SpaceElemHandler &space_elem_handler,
                        SafeSTLArray<ValueFlags, dim+1> &flags)
            :
            flag_in_(flag_in),
            elements_flat_id_(elements_flat_id),
            space_elem_handler_(space_elem_handler),
            flags_(flags)
        {}

        template<int sub_elem_dim>
        void operator()(const Quadrature<sub_elem_dim> &quad)
        {
            flags_[sub_elem_dim] = flag_in_;
            space_elem_handler_.reset_selected_elements(flag_in_,quad,elements_flat_id_);
        }

        const ValueFlags flag_in_;

        /**
         * Elements to reset.
         */
        const SafeSTLVector<Index> &elements_flat_id_;

        SpaceElemHandler &space_elem_handler_;

        SafeSTLArray<ValueFlags, dim+1> &flags_;

    };


    struct InitCacheDispatcher : boost::static_visitor<void>
    {
        InitCacheDispatcher(
            SpaceElemHandler &space_elem_handler,
            SpaceElem &space_elem)
            :
            space_elem_handler_(space_elem_handler),
            space_elem_(space_elem)
        {}

        template<int sub_elem_dim>
        void operator()(const Topology<sub_elem_dim> &sub_elem)
        {
            space_elem_handler_.template init_cache<sub_elem_dim>(space_elem_);
        }

        SpaceElemHandler &space_elem_handler_;
        SpaceElem &space_elem_;
    };

    struct FillCacheDispatcher : boost::static_visitor<void>
    {
        FillCacheDispatcher(const int sub_elem_id,
                            self_t &function,
                            SpaceElemHandler &space_elem_handler,
                            ElementAccessor &func_elem,
                            SpaceElem &space_elem,
                            SafeSTLVector<Real> &loc_coeff,
                            const std::string  &property)
            :
            sub_elem_id_(sub_elem_id),
            function_(function),
            space_elem_handler_(space_elem_handler),
            func_elem_(func_elem),
            space_elem_(space_elem),
            loc_coeff_(loc_coeff),
            property_(property)
        {}


        template<int sub_elem_dim>
        void operator()(const Topology<sub_elem_dim> &sub_elem)
        {
            space_elem_handler_.template fill_cache<sub_elem_dim>(space_elem_,sub_elem_id_);

            auto &local_cache = function_.get_cache(func_elem_);
            auto &cache = local_cache->template get_sub_elem_cache<sub_elem_dim>(sub_elem_id_);

#if 0
            boost::fusion::for_each(cache.get_values(),
                                    [&](auto & type_and_value) -> void
            {
                using ValueType_ValueContainer = typename std::remove_reference<decltype(type_and_value)>::type;
                using ValueType = typename ValueType_ValueContainer::first_type;
                auto &value = type_and_value.second;

                if (value.status_fill())
                {
                    value = space_elem_->template linear_combination<ValueType,sub_elem_dim>(*loc_coeff_,j, *property_);
                    value.set_status_filled(true);
                }
            } // end lambda function
                                                            );
#endif
//#if 0
            //TODO (martinelli Mar 27,2015): bad style. Use the ValueType mechanism in order to avoid the if-switch
            if (cache.template status_fill<_Value>())
            {
                cache.template get_data<_Value>() =
                    space_elem_.template linear_combination<_Value,sub_elem_dim>(loc_coeff_,sub_elem_id_,property_);
                cache.template set_status_filled<_Value>(true);
            }
            if (cache.template status_fill<_Gradient>())
            {
                cache.template get_data<_Gradient>() =
                    space_elem_.template linear_combination<_Gradient,sub_elem_dim>(loc_coeff_,sub_elem_id_, property_);
                cache.template set_status_filled<_Gradient>(true);
            }
            if (cache.template status_fill<_Hessian>())
            {
                cache.template get_data<_Hessian>() =
                    space_elem_.template linear_combination<_Hessian,sub_elem_dim>(loc_coeff_,sub_elem_id_, property_);
                cache.template set_status_filled<_Hessian>(true);
            }
            if (cache.template status_fill<_Divergence>())
            {
                cache.template get_data<_Divergence>() =
                    space_elem_.template linear_combination<_Divergence,sub_elem_dim>(loc_coeff_,sub_elem_id_, property_);
                cache.template set_status_filled<_Divergence>(true);
            }
//#endif
            cache.set_filled(true);
        }

        const int sub_elem_id_;
        self_t &function_;
        SpaceElemHandler &space_elem_handler_;
        ElementAccessor &func_elem_;
        SpaceElem &space_elem_;
        SafeSTLVector<Real> &loc_coeff_;
        const std::string  &property_;
    };
#endif

#ifdef MESH_REFINEMENT

    void create_connection_for_insert_knots(std::shared_ptr<self_t> ig_function);

    void rebuild_after_insert_knots(
        const SafeSTLArray<SafeSTLVector<Real>,dim> &knots_to_insert,
        const CartesianGrid<dim> &old_grid);

#endif // MESH_REFINEMENT


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

#endif
