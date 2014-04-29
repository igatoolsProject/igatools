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

#ifndef ANALYTICAL_MAPPING_H_
#define ANALYTICAL_MAPPING_H_

#include <igatools/base/config.h>
#include <igatools/geometry/mapping.h>

IGA_NAMESPACE_OPEN


/**
 * @brief Base class for analytical mappings.
 *
 * @author M. Martinelli
 * @author 13 Mar 2014
 */
template<int dim_, int codim_>
class AnalyticalMapping : public Mapping <dim_, codim_>
{
public:
    using base_t = Mapping<dim_,codim_>;

    using base_t::dim;
    using base_t::codim;
    using base_t::space_dim;

    using typename base_t::PointType;
    using typename base_t::ValueType;
    using typename base_t::GradientType;
    using typename base_t::HessianType;
    using typename base_t::GridType;

    /** @name Constructor and destructor */
    ///@{
    /** The constructors are inherited from Mapping. */
    using Mapping<dim_,codim_>::Mapping;
    ///@}


    /**
     * It does nothing rather than providing a concrete implementation of the
     * the pure virtual function Mapping::init_element().
     */
    void init_element(const ValueFlags flag,
                      const Quadrature<dim> &quad) const override final {} ;
};


IGA_NAMESPACE_CLOSE


#endif // #ifndef ANALYTICAL_MAPPING_H_
