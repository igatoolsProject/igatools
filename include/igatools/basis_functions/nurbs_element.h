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


#ifndef NURBS_ELEMENT_H_
#define NURBS_ELEMENT_H_

#include <igatools/base/config.h>

#ifdef NURBS

#include <igatools/basis_functions/space_element.h>
#include <igatools/basis_functions/nurbs_element_handler.h>


//#include <igatools/linear_algebra/dense_matrix.h>
//#include <igatools/basis_functions/bernstein_basis.h>
//#include <igatools/basis_functions/bspline_element_scalar_evaluator.h>

IGA_NAMESPACE_OPEN


template <int dim, int range, int rank> class NURBSSpace;
template <typename Accessor> class GridForwardIterator;

/**
 * See module on \ref accessors_iterators for a general overview.
 * @ingroup accessors
 */
template <int dim, int range, int rank>
class NURBSElement :
    public SpaceElement<NURBSSpace<dim,range,rank>>
{
private:
    using self_t = NURBSElement<dim,range,rank>;
    using parent_t = SpaceElement<NURBSSpace<dim,range,rank>>;

public:
    /** Type for the grid accessor. */
    using GridAccessor = CartesianGridElement<dim>;

    /** Type required by the GridForwardIterator templated iterator */
    using ContainerType = const NURBSSpace<dim, range, rank> ;

    /** Type required for the generic algorithm on the spaces (plots??) */
    using Space = NURBSSpace<dim, range, rank> ;

    /** @name Constructors */
    ///@{
    /**
     * Default constructor
     */
    NURBSElement() = default;

    /**
     * Constructs an accessor to element number index of a
     * BsplineSpace space.
     */
    NURBSElement(const std::shared_ptr<ContainerType> space,
                 const Index elem_index)
    {
        Assert(false,ExcNotImplemented());
    }

    NURBSElement(const std::shared_ptr<ContainerType> space,
                 const TensorIndex<dim> &elem_index)
    {
        Assert(false,ExcNotImplemented());
    }

    /**
     * Copy constructor.
     * It can be used with different copy policies (i.e. deep copy or shallow copy).
     * The default behaviour (i.e. using the proper interface of a classic copy constructor)
     * uses the deep copy.
     */
    NURBSElement(const self_t &elem,
                 const CopyPolicy &copy_policy = CopyPolicy::deep)
    {
        Assert(false,ExcNotImplemented());
    }

    /**
     * Move constructor.
     */
    NURBSElement(self_t &&elem) = default;

    /**
     * Destructor.
     */
    ~NURBSElement() = default;
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

};

IGA_NAMESPACE_CLOSE

#endif // #ifdef NURBS


#endif // end of #ifndef NURBS_ELEMENT_H_



