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


#ifndef NURBS_ELEMENT_H_
#define NURBS_ELEMENT_H_

#include <igatools/base/config.h>

#ifdef NURBS

#include <igatools/basis_functions/reference_element.h>


//#include <igatools/linear_algebra/dense_matrix.h>
//#include <igatools/basis_functions/bernstein_basis.h>
//#include <igatools/basis_functions/bspline_element_scalar_evaluator.h>



IGA_NAMESPACE_OPEN




template <int dim, int range, int rank> class NURBSSpace;
template <int,int,int> class NURBSElementHandler;
template <class Accessor> class GridIterator;

/**
 * @brief NURBS element
 *
 * See module on \ref accessors_iterators for a general overview.
 * @ingroup elements
 * @ingroup serialization
 */
template <int dim, int range, int rank>
class NURBSElement :
    public ReferenceElement<dim,range,rank>
{
private:
    using self_t = NURBSElement<dim,range,rank>;
    using parent_t = ReferenceElement<dim,range,rank>;

public:
    /** Type for the grid accessor. */
    using GridAccessor = CartesianGridElement<dim>;

    /** Type required by the CartesianGridIterator templated iterator */
    using ContainerType = const NURBSSpace<dim, range, rank> ;

    /** Type required for the generic algorithm on the spaces (plots??) */
    using Space = NURBSSpace<dim, range, rank> ;

    /** @name Constructors */
    ///@{
protected:
    /**
     * Default constructor. It does nothing but it is needed for the
     * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     * mechanism.
     */
    NURBSElement() = default;

public:
    /**
     * Constructs an accessor to element number index of a
     * BsplineSpace space.
     */
    NURBSElement(const std::shared_ptr<ContainerType> space,
                 const Index elem_index);


    /**
     * Copy constructor.
     * It can be used with different copy policies (i.e. deep copy or shallow copy).
     * The default behaviour (i.e. using the proper interface of a classic copy constructor)
     * uses the deep copy.
     */
    NURBSElement(const self_t &elem,
                 const CopyPolicy &copy_policy = CopyPolicy::deep);
//*/

    /**
     * Move constructor.
     */
    NURBSElement(self_t &&elem) = default;

    /**
     * Destructor.
     */
    virtual ~NURBSElement() = default;
    ///@}

    /** @name Assignment operators */
    ///@{
    /**
     * Copy assignment operator.
     * @note Creates a new element cache, but it shares
     * the one dimensional cache with the copied element.
     */
    self_t &operator=(const self_t &elem) = default;

    /**
     * Move assignment operator.
     */
    self_t &operator=(self_t &&elem) = default;
    ///@}

    /** @name Functions/operators for moving the element in the NURBSSpace.*/
    ///@{
    /**
     * Sets the index of the element using the flatten representation.
     * @note This function also updates the index for the tensor representation.
     * @warning This may be a dangerous function, be careful when using it
     * as it is easy to use incorrectly. Only use it if you know what you
     * are doing.
     */
    virtual void move_to(const Index flat_index) override final;
    ///@}


    /**
     * Returns the NURBSSpace in which the NURBSElement is defined.
     */
    std::shared_ptr<const Space> get_nurbs_space() const;

private:

    using SpSpace = typename Space::SpSpace;

    typename SpSpace::ElementAccessor bspline_elem_;

    using WeightFunction = typename Space::WeightFunction;
    using WeightElem = typename WeightFunction::ElementAccessor;
    WeightElem weight_elem_;

public:

    friend class NURBSElementHandler<dim, range, rank>;


    virtual std::shared_ptr<SpaceElement<dim,0,range,rank,Transformation::h_grad> >
    clone() const override final;

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

#endif // #ifdef NURBS


#endif // end of #ifndef NURBS_ELEMENT_H_



