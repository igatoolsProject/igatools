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

#ifndef MAPPING_SLICE_H_
#define MAPPING_SLICE_H_

#include <igatools/geometry/mapping.h>
#include <igatools/geometry/mapping_element_accessor.h>

IGA_NAMESPACE_OPEN

/**
 * \brief It is a mapping obtained as the restriction of a given mapping
 * over a coordinate plane.
 * For example g(x,y) = F(c,x,y)
 */
template<int dim_, int codim_>
class MappingSlice : public Mapping<dim_, codim_>
{
public:
    using base_t = Mapping<dim_, codim_>;
    using base_t::dim;
    using base_t::codim;
    using base_t::space_dim;

    using typename base_t::GridType;
    using typename base_t::GridIterator;

    using typename base_t::Point;
    using typename base_t::Value;
    using typename base_t::Gradient;
    using typename base_t::Hessian;

private:
    using SupMap = Mapping<dim + 1, codim - 1>;
    using self_t = MappingSlice<dim, codim>;

public:
    /**
     *  Given a <em>map</em>, we restrict to a slice perpendicular to
     *  reference <em>direction</em> and located at reference <em>value</em>
     *  on <em>direction</em>.
     */
    MappingSlice(const std::shared_ptr<const SupMap> map,
                 const int face_id,
                 const std::shared_ptr<GridType> grid,
                 const std::shared_ptr<typename SupMap::GridType::FaceGridMap> elem_map);

    /**
     * Copy constructor
     */
    MappingSlice(const self_t &map_slice);

    /**
     * Copy assignment operator. Not allowed to be used.
     */
    self_t &operator=(const self_t &map) = delete;

    static std::shared_ptr<base_t>
    create(const std::shared_ptr<const SupMap> map,
           const int face_id,
           const std::shared_ptr<GridType> grid,
           const std::shared_ptr<typename SupMap::GridType::FaceGridMap> elem_map);

    static std::shared_ptr<base_t>
    create(const std::shared_ptr<const base_t> map,
           const int face_id,
           const std::shared_ptr<GridType > grid,
           const std::shared_ptr<typename base_t::GridType::FaceGridMap> elem_map)
    {
        AssertThrow(true, ExcImpossibleInDim(-1));
        return std::shared_ptr<base_t>();//Should never reach this
    }

    void evaluate(vector<Value> &values) const override;

    void evaluate_gradients(vector<Gradient> &gradients) const override;

    void init_element(const ValueFlags flag, const Quadrature<dim> &quad)  const override;

    void set_element(const GridIterator &elem) const override ;

    void set_face_element(const Index face_id,
                          const GridIterator &elem) const override;

    /**
     * Prints internal information about the mapping.
     * @note Mostly used for debugging and testing.
     */
    void print_info(LogStream &out) const override;

    //TODO: should be private
public:
    Quadrature<dim+1>
    build_extended_quadrature(const Quadrature<dim> &quad) const;

private:
    const std::shared_ptr<const SupMap> map_;
    const int direction_;
    const int value_;

    // the cache
    mutable typename SupMap::ElementIterator element;

    const std::shared_ptr<typename SupMap::GridType::FaceGridMap> elem_map_;

};



IGA_NAMESPACE_CLOSE

#endif /* MAPPING_SLICE_H_ */
