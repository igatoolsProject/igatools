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

#ifndef PHYS_SPACE_ELEMENT_HANDLER_H_
#define PHYS_SPACE_ELEMENT_HANDLER_H_

#include <igatools/base/config.h>
#include <igatools/basis_functions/bspline_element_handler.h>
#include <igatools/basis_functions/nurbs_element_handler.h>
#include <igatools/geometry/push_forward.h>

IGA_NAMESPACE_OPEN

/**
 * Element handler for an isogeometric space
 */
template<class PhysSpace>
class PhysSpaceElementHandler
    : public ElementHandler<PhysSpace>
//      : public PhysSpace::PushForwardType
{

    using RefSpace =  typename PhysSpace::RefSpace;
    using RefPhysSpaceElementHandler = typename PhysSpace::RefSpace::ElementHandler;
    using PFCache = typename PhysSpace::PushForwardType;

    using ElementIterator = typename PhysSpace::ElementIterator;
    using ElementAccessor = typename PhysSpace::ElementAccessor;
    using PfElemAccessor = typename PhysSpace::PushForwardType::ElementAccessor;

    using self_t = PhysSpaceElementHandler<PhysSpace>;

public:
    static const int dim = PhysSpace::dim;

//    using PhysSpace::PushForwardType::type;

    /**
     * @name Constructors
     */
    ///@{
    /**
     * Default constructor. Not allowed to be used.
     */
    PhysSpaceElementHandler() = delete;

    PhysSpaceElementHandler(std::shared_ptr<const PhysSpace> space);

    /**
     * Copy constructor. Not allowed to be used.
     */
    PhysSpaceElementHandler(const self_t &) = delete;

    /**
     * Move constructor. Not allowed to be used.
     */
    PhysSpaceElementHandler(self_t &&) = delete;

    /**
     * Destructor.
     */
    ~PhysSpaceElementHandler() = default;
    ///@}

    /**
     * Assignment operators.
     */
    ///@{
    /**
     * Copy assignment operator. Not allowed to be used.
     */
    self_t &operator=(const self_t &) = delete;

    /**
     * Move assignment operator. Not allowed to be used.
     */
    self_t &operator=(self_t &&) = delete;
    ///@}

    /**
     * @name Creators.
     */
    ///@{
    static std::shared_ptr<self_t> create(std::shared_ptr<const PhysSpace> space);
    ///@}


    template<int k>
    void reset(const ValueFlags flag, const EvaluationPoints<k> &eval_pts);


    template<int k>
    void reset_selected_elements(
        const ValueFlags &flag,
        const EvaluationPoints<k> &eval_pts,
        const vector<Index> &elements_flat_id)
    {
        Assert(false,ExcNotImplemented());
        AssertThrow(false,ExcNotImplemented());
    }



    //protected:
    template <int k>
    void fill_cache(ElementAccessor &elem, const int j);


    template <int k>
    void init_cache(ElementAccessor &elem);

    void print_info(LogStream &out) const;

private:
    std::shared_ptr<const PhysSpace> space_;

    std::shared_ptr<typename PhysSpace::RefSpace::ElementHandler> ref_space_handler_;


    typename PhysSpace::PushForwardType push_fwd_;

    std::array<FunctionFlags, dim + 1> flags_;

};


IGA_NAMESPACE_CLOSE

#endif
