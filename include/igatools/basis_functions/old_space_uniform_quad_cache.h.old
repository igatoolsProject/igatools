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

#ifndef SPACE_UNIFORM_QUAD_CACHE_H_
#define SPACE_UNIFORM_QUAD_CACHE_H_
#if 0
#include <igatools/base/config.h>
#include <igatools/basis_functions/bspline_uniform_quad_cache.h>
#include <igatools/geometry/push_forward_uniform_quad_cache.h>

IGA_NAMESPACE_OPEN

template<class PhysSpace>
class SpaceUniformQuadCache :
    public PhysSpace::RefSpace::UniformQuadCache,
    public PushFowardUniformQuadCache<typename PhysSpace::PushForwardType>
{
    using RefSpace =  typename PhysSpace::RefSpace;
    using RefSpaceCache = typename PhysSpace::RefSpace::UniformQuadCache;
    using PFCache = PushFowardUniformQuadCache<typename PhysSpace::PushForwardType>;
    using ElementIterator = typename PhysSpace::ElementIterator;
    using PfElemAccessor = typename PhysSpace::PushForwardType::ElementAccessor;

protected:
    using ElementAccessor = typename PhysSpace::ElementAccessor;

    void init_element_cache(ElementAccessor &elem);

    void fill_element_cache(ElementAccessor &elem);

    void fill_face_cache(ElementAccessor &elem, const int face);

public:
    static const int dim = PhysSpace::dim;

    //Allocates and fill the (global) cache
    SpaceUniformQuadCache(std::shared_ptr<const PhysSpace> space,
                          const ValueFlags flag,
                          const Quadrature<dim> &quad);

    //Allocates the ElementIterator element_cache
    void init_element_cache(ElementIterator &elem);

    //Fill the ElementIterator element_cache
    void fill_element_cache(ElementIterator &elem);

    /**
     * Fills the ElementIterator face_cache
     * element dependent part
     */
    void fill_face_cache(ElementIterator &elem, const int face);

    void print_info(LogStream &out) const;

private:
    std::shared_ptr<const PhysSpace> space_;

    ValueFlags flags_;
    // BasisElemValueFlagsHandler flags_;

    //  BasisFaceValueFlagsHandler face_flags_;

    Quadrature<dim> quad_;
};


IGA_NAMESPACE_CLOSE
#endif
#endif
