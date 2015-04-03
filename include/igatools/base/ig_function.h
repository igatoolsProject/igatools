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

#ifndef IG_FUNCTIONS_H
#define IG_FUNCTIONS_H

#include <igatools/base/function.h>
#include <igatools/basis_functions/spline_space.h>
#include <igatools/linear_algebra/epetra.h>

IGA_NAMESPACE_OPEN

using IgCoefficients = EpetraTools::Vector;

#if 0
class IgCoefficients : public std::map<Index,Real>
{
public:
    using std::map<Index,Real>::map;

    IgCoefficients(const std::set<Index> &global_dofs, const vector<Real> &coeffs)
    {
        Assert(Index(global_dofs.size()) == coeffs.size(),
               ExcDimensionMismatch(global_dofs.size(),coeffs.size()));

        int i = 0;
        for (const auto dof : global_dofs)
            (*this)[dof] = coeffs[i++];
    }

    template <class Space>
    IgCoefficients(
        const Space &space,
        const std::string &dofs_property,
        const vector<Real> &coeffs)
        :
        IgCoefficients(space.get_dof_distribution()->get_dofs_id_same_property(dofs_property),coeffs)
    {}

    template <class Space>
    IgCoefficients(
        const Space &space,
        const std::string &dofs_property)
        :
        IgCoefficients(
            space.get_dof_distribution()->get_dofs_id_same_property(dofs_property),
            vector<Real>(space.get_dof_distribution()->
                         get_dofs_id_same_property(dofs_property).size(),0.0))
    {}


    Real &operator()(const Index &global_dof)
    {
#ifdef NDEBUG
        return (*this)[global_dof];
#else
        return (*this).at(global_dof);
#endif
    }

    const Real &operator()(const Index &global_dof) const
    {
#ifdef NDEBUG
        return (*this)[global_dof];
#else
        return (*this).at(global_dof);
#endif
    }


    Size size() const
    {
        return std::map<Index,Real>::size();
    }

    IgCoefficients &operator+=(const IgCoefficients &coeffs)
    {
        Assert(this->size() == coeffs.size(),ExcDimensionMismatch(this->size(),coeffs.size()));
#ifdef NDEBUG
        for (const auto &c : coeffs)
            (*this)[c.first] += c.second;
#else
        for (const auto &c : coeffs)
            (*this).at(c.first) += c.second;
#endif

        return *this;
    }

    vector<Real> get_local_coeffs(const vector<Index> &elem_dofs) const
    {

       vector<Real> loc_coeff;
       for (const auto &dof : elem_dofs)
           loc_coeff.emplace_back((*this)(dof));
       return  loc_coeff;
    }

    void print_info(LogStream &out) const
    {
        int i = 0;
        for (const auto &c : (*this))
        {
            out << "coeff[local_id=" << i << ", global_id=" << c.first << "] = " << c.second << std::endl;
            ++i;
        }
    }
};
#endif

template<class Space>
class IgFunction :
    public Function<Space::dim, Space::codim, Space::range, Space::rank>,
    public std::enable_shared_from_this<IgFunction<Space>>
{
public:
    static const int dim = Space::dim;
    static const int codim = Space::codim;
    static const int range = Space::range;
    static const int rank = Space::rank;

    using CoeffType = IgCoefficients;


private:
    using base_t = Function<dim, codim, range, rank>;
    using parent_t = Function<dim, codim, range, rank>;
    using self_t = IgFunction<Space>;

    std::shared_ptr<const base_t> shared_from_derived() const override final
    {
        return this->shared_from_this();
    }

public:
    //TODO (pauletti, Mar 23, 2015): should we make this private?
    IgFunction(std::shared_ptr<Space> space,
               const CoeffType &coeff,
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
    create(std::shared_ptr<Space> space, const CoeffType &coeff,
           const std::string &property = DofProperties::active);


    std::shared_ptr<base_t> clone() const override final
    {
        return std::make_shared<self_t>(self_t(*this));
    }


    void reset(const ValueFlags &flag, const eval_pts_variant &eval_pts) override;

    void reset_selected_elements(
        const ValueFlags &flag,
        const eval_pts_variant &eval_pts,
        const vector<Index> &elements_flat_id);

    void init_cache(ElementAccessor &elem, const topology_variant &k) override;

    void fill_cache(ElementAccessor &elem, const topology_variant &k, const int j) override;

    std::shared_ptr<const Space> get_ig_space() const;

    const CoeffType &get_coefficients() const;

    const std::string &get_property() const
    {
        return property_;
    }

    self_t &operator +=(const self_t &fun);

    void print_info(LogStream &out) const;

private:

    std::shared_ptr<Space> space_;

    CoeffType coeff_;

    std::string property_;

    typename Space::ElementIterator elem_;

    std::shared_ptr<typename Space::ElementHandler> space_filler_;

private:
    struct ResetDispatcher : boost::static_visitor<void>
    {
        template<class T>
        void operator()(const T &quad)
        {
            (*flags_)[T::dim] = flag;
            Assert(space_handler_ != nullptr, ExcNullPtr());
            space_handler_->reset_selected_elements(flag, quad,*elements_flat_id_);
        }

        ValueFlags flag;
        typename Space::ElementHandler *space_handler_;
        std::array<FunctionFlags, dim + 1> *flags_;

        /**
         * Elements to reset.
         */
        const vector<Index> *elements_flat_id_;
    };


    struct InitCacheDispatcher : boost::static_visitor<void>
    {
        template<class T>
        void operator()(const T &quad)
        {
            Assert(space_handler_ != nullptr, ExcNullPtr());
            space_handler_->template init_cache<T::k>(*space_elem);
        }

        typename Space::ElementHandler  *space_handler_;
        typename Space::ElementAccessor *space_elem;
    };


    struct FillCacheDispatcher : boost::static_visitor<void>
    {
        template<class T>
        void operator()(const T &quad)
        {
            Assert(space_handler_ != nullptr, ExcNullPtr());
            space_handler_->template fill_cache<T::k>(*space_elem,j);

            auto &local_cache = function->get_cache(*func_elem);
            auto &cache = local_cache->template get_value_cache<T::k>(j);
            auto &flags = cache.flags_handler_;

            if (flags.fill_values())
                std::get<0>(cache.values_) =
                    space_elem->template linear_combination<0,T::k>(*loc_coeff,j, *property);
            if (flags.fill_gradients())
                std::get<1>(cache.values_) =
                    space_elem->template linear_combination<1,T::k>(*loc_coeff,j, *property);
            if (flags.fill_hessians())
                std::get<2>(cache.values_) =
                    space_elem->template linear_combination<2,T::k>(*loc_coeff,j, *property);

            cache.set_filled(true);
        }

        int j;
        self_t *function;
        typename Space::ElementHandler *space_handler_;
        ElementAccessor *func_elem;
        typename Space::ElementAccessor *space_elem;
        vector<Real> *loc_coeff;
        std::string *property;
    };

    ResetDispatcher reset_impl;
    InitCacheDispatcher init_cache_impl;
    FillCacheDispatcher fill_cache_impl;


    void create_connection_for_insert_knots(std::shared_ptr<self_t> ig_function);

    void rebuild_after_insert_knots(
        const special_array<vector<Real>,dim> &knots_to_insert,
        const CartesianGrid<dim> &old_grid);


};

IGA_NAMESPACE_CLOSE

#endif
