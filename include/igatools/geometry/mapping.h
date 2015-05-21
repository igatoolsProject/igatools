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

#ifndef NEW_MAPPING_H_
#define NEW_MAPPING_H_

#include <igatools/base/config.h>
#include <igatools/base/function.h>

IGA_NAMESPACE_OPEN

//Forward declaration to avoid including header file.
template <int, int> class MappingElement;

template <int,int,int,int> class IgFunction;

/**
 * @brief The mapping is a deformation \f$ F : \hat\Omega \to \Omega\f$
 * which maps the reference domain \f$\hat\Omega \in \mathbb{R}^{dim}\f$ to the
 * physical domain \f$\Omega \in \mathbb{R}^{dim+codim}\f$.
 *
 * In igatools the mapping is a grid-like container, so that on top
 * of being a deformation is aware of an underlying grid structure.
 *
 * Most of the methods in this class are virtual,
 * because we want to use Mapping as general interface for user-defined Mapping
 * specializations (unknown at compile time).
 *
 * @ingroup containers
 * @ingroup serializable
 *
 * @author pauletti 2014
 * @author M. Martinelli, 2015
 */
template<int dim_, int codim_ = 0>
class Mapping :
    public std::enable_shared_from_this<Mapping<dim_,codim_> >
{
private:
    using self_t = Mapping<dim_, codim_>;
    using FuncType = MapFunction<dim_, dim_ + codim_>;
public:
    using ElementAccessor = MappingElement<dim_, codim_>;
    using ElementIterator = CartesianGridIterator<ElementAccessor>;

    static const int dim = dim_;
    static const int space_dim = dim_ + codim_;

public:
    /** Type for the given order derivatives of the
     *  the mapping. */
    template<int order>
    using Derivative = typename FuncType::template Derivative<order>;


    /** Type for the diferent order derivatives of the inverse of
     * the mapping
     */
    template<int order>
    using InvDerivative = Derivatives<space_dim, dim_, 1, order>;


    /** Type of the mapping evaluation point. */
    using Point = typename FuncType::Point;

    /** Type of the mapping return value. */
    using Value = typename FuncType::Value;

    /** Type of the mapping gradient. */
    using Gradient = typename FuncType::Gradient;

    /** Typedef for the mapping hessian. */
    using Hessian = typename FuncType::Hessian;

public:

    /**
     * Default constructor. It does nothing but it is needed for the
     * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     * mechanism.
     */
    Mapping() = default;

    Mapping(std::shared_ptr<FuncType> F);


    ~Mapping();

    static std::shared_ptr<self_t>  create(std::shared_ptr<FuncType> F);

public:
    template<int k>
    void reset(const ValueFlags flag, const Quadrature<k> &eval_pts);

    template <int k>
    void init_cache(ElementAccessor &elem);

    template <int k>
    void init_cache(ElementIterator &elem)
    {
        init_cache<k>(*elem);
    }

    template <int k>
    void fill_cache(ElementAccessor &elem, const int j);

    template <int k>
    void fill_cache(ElementIterator &elem, const int j)
    {
        fill_cache<k>(*elem, j);
    }


    std::shared_ptr<const CartesianGrid<dim_> > get_grid() const;

    std::shared_ptr<FuncType> get_function() const;

    std::shared_ptr<ElementAccessor> create_element(const Index flat_index) const;

    ElementIterator begin() const;

    ElementIterator end();


private:
    std::shared_ptr<FuncType> F_;

    SafeSTLArray<ValueFlags, dim_ + 1> flags_;

    friend ElementAccessor;

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
        ar.template register_type<IgFunction<dim_,0,dim_+codim_,1> >();
        ar &boost::serialization::make_nvp("F_",F_);
        ar &boost::serialization::make_nvp("flags_",flags_);
    }
    ///@}
#endif

};

IGA_NAMESPACE_CLOSE

#endif
