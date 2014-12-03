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

#ifndef SPACE_ELEMENT_HANDLER_H_
#define SPACE_ELEMENT_HANDLER_H_

#include <igatools/base/config.h>
#include <igatools/basis_functions/bspline_element_handler.h>
#include <igatools/basis_functions/nurbs_element_handler.h>
#include <igatools/geometry/new_push_forward.h>

IGA_NAMESPACE_OPEN

/**
 * Element handler for an isogeometric space
 */
template<class PhysSpace>
class SpaceElementHandler :
    public PhysSpace::RefSpace::ElementHandler,
    public PhysSpace::PushForwardType
{

    using RefSpace =  typename PhysSpace::RefSpace;
    using RefSpaceElementHandler = typename PhysSpace::RefSpace::ElementHandler;
    using PFCache = typename PhysSpace::PushForwardType;

    using ElementIterator = typename PhysSpace::ElementIterator;
    using ElementAccessor = typename PhysSpace::ElementAccessor;
    using PfElemAccessor = typename PhysSpace::PushForwardType::ElementAccessor;

public:
    static const int dim = PhysSpace::dim;

    using PhysSpace::PushForwardType::type;

    //Allocates and fill the (global) cache
    SpaceElementHandler(std::shared_ptr<const PhysSpace> space);

    template<int k>
    void reset(const NewValueFlags flag, const Quadrature<k> &quad);

    //protected:
    template <int k>
    void fill_cache(ElementAccessor &elem, const int j);

    template <int k>
    void init_cache(ElementAccessor &elem);

    //    void init_all_caches(ElementAccessor &elem);
public:
    template <int k>
    void fill_cache(ElementIterator &elem, const int j)
    {
        fill_cache<k>(elem.get_accessor(), j);
    }

    template <int k>
    void init_cache(ElementIterator &elem)
    {
        init_cache<k>(elem.get_accessor());
    }

    //    void init_all_caches(ElementIterator &elem)
    //    {
    //        init_all_caches(elem.get_accessor());
    //    }


    void print_info(LogStream &out) const;

private:
    std::shared_ptr<const PhysSpace> space_;

    std::array<FunctionFlags, dim + 1> flags_;

};


IGA_NAMESPACE_CLOSE

#endif
