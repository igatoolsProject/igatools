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


#ifndef REFERENCE_ELEMENT_H_
#define REFERENCE_ELEMENT_H_

#include <igatools/base/config.h>
#include <igatools/basis_functions/space_element.h>
#include <igatools/basis_functions/reference_element_handler.h>

IGA_NAMESPACE_OPEN


template <int, int, int> class ReferenceSpace;

/**
 *
 * @ingroup elements
 * @ingroup serializable
 */
template <int dim, int range, int rank>
class ReferenceElement : public SpaceElement<dim,0,range,rank,Transformation::h_grad>
{
public:
    /** Type for the grid accessor. */
    using GridAccessor = CartesianGridElement<dim>;

    /** Type required by the GridForwardIterator templated iterator */
    using ContainerType = const ReferenceSpace<dim,range,rank> ;

    using Space = ReferenceSpace<dim,range,rank>;
    using ConstSpace = const ReferenceSpace<dim,range,rank>;

    using parent_t = SpaceElement<dim,0,range,rank,Transformation::h_grad>;

    using RefPoint = typename Space::RefPoint;
    using Point = typename Space::Point;
    using Value = typename Space::Value;

    template <int order>
    using Derivative = typename Space::template Derivative<order>;

    using Div = typename Space::Div;

protected:
    /**
     * Default constructor. It does nothing but it is needed for the
     * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     * mechanism.
     */
    ReferenceElement() = default;

public:
    ReferenceElement(const ReferenceElement<dim,range,rank> &elem,
                     const iga::CopyPolicy &copy_policy = CopyPolicy::deep);

    /**
     * Constructs an accessor to element number index of a
     * ReferenceSpace space.
     */
    ReferenceElement(const std::shared_ptr<ConstSpace> space,
                     const Index elem_index);


    virtual ~ReferenceElement() = default;






    /**
     * Returns the <tt>k</tt> dimensional j-th sub-element measure
     * multiplied by the weights of the quadrature.
     */
    template <int k>
    ValueVector<Real> get_w_measures(const int j) const
    {
        return this->as_cartesian_grid_element_accessor().template get_w_measures<k>(j);
    }

    /**
     * Returns the gradient determinant of the identity map at the dilated quadrature points.
     */
    ValueVector<Real> get_element_w_measures() const;


    using OffsetTable = typename Space::template ComponentContainer<int>;
    using TensorSizeTable = typename Space::TensorSizeTable;

protected:

    /** Number of scalar basis functions along each direction, for all space components. */
    TensorSizeTable n_basis_direction_;

    /** Basis function ID offset between the different components. */
    OffsetTable comp_offset_;

    using Indexer = CartesianProductIndexer<dim>;
    using IndexerPtr = std::shared_ptr<Indexer>;
    using IndexerPtrTable = typename Space::template ComponentContainer<IndexerPtr>;

    /** Hash table for fast conversion between flat-to-tensor basis function ids. */
    IndexerPtrTable basis_functions_indexer_;

    std::shared_ptr<const Space> space_;

public:
    using parent_t::get_num_basis;


    /**
     * Returns the basis function ID offset between the different components.
     */
    OffsetTable get_basis_offset() const;

    /**
     * Number of non-zero scalar basis functions associated
     * with the i-th space component on the element.
     * This makes sense as a reference B-spline space
     * is only allowed to be of the cartesian product type
     * V = V1 x V2 x ... X Vn.
     */
    int get_num_basis_comp(const int i) const;

    virtual void print_info(LogStream &out) const override final;


    /**
     *  Maximum number of non zero basis functions, over the current element.
     *  @note The "true" number of basis functions may differ from this value because
     *  some basis functions may be defined to be "inactive".
     */
    virtual Size get_max_num_basis() const override final;


public:
    std::shared_ptr<const Space> get_space() const;

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
    serialize(Archive &ar, const unsigned int version);
    ///@}
#endif // SERIALIZATION
};



IGA_NAMESPACE_CLOSE


#endif // #ifndef REFERENCE_ELEMENT_H_

