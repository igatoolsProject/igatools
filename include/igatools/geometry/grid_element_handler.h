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

#ifndef GRID_ELEMENT_HANDLER_H_
#define GRID_ELEMENT_HANDLER_H_

#include <igatools/base/config.h>
#include <igatools/basis_functions/values_cache.h>
#include <igatools/base/tuple_utils.h>
#include <igatools/base/quadrature.h>
#include <igatools/utils/tensor_product_array.h>
#include <igatools/geometry/cartesian_grid.h>
#include <igatools/basis_functions/element_handler.h>

IGA_NAMESPACE_OPEN

/**
 * Grid element value manager
 *
 * computational optimization cache, storing the interval length
 * in each direction.
 *
 */
template <int dim>
class GridElementHandler
    : public ElementHandler<CartesianGrid<dim>>
{
private:
    using self_t = GridElementHandler<dim>;

public:
    using GridType = const CartesianGrid<dim>;

    /**
     * Alias for the (static) class holding the topological information.
     */
    using Topology = UnitElement<dim>;

protected:
    using ElementIterator = typename GridType::ElementIterator;
    using ElementAccessor = typename GridType::ElementAccessor;

public:
    /**
     * @name Creators.
     */
    ///@{
    static std::shared_ptr<GridElementHandler<dim>> create(std::shared_ptr<GridType> grid);
    ///@}

    /**
     * @name Constructors
     */
    ///@{
    /**
     * Default constructor. Not allowed to be used.
     */
    GridElementHandler() = delete;

    /**
     * Constructor.
     */
    GridElementHandler(std::shared_ptr<GridType> grid);

    /**
     * Copy constructor.
     */
    GridElementHandler(const self_t &) = default;

    /**
     * Move constructor.
     */
    GridElementHandler(self_t &&) = default;

    /**
     * Destructor.
     */
    ~GridElementHandler() = default;
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

public:
    /**
     * @name Functions for the cache's reset/init/fill mechanism.
     */
    ///@{
    template<int k>
    void reset(const ValueFlags flag, const Quadrature<k> &quad);

    template <int k>
    void init_cache(ElementAccessor &elem);

#if 0
    void init_all_caches(ElementAccessor &elem);

    void init_all_caches(ElementIterator &elem)
    {
        init_all_caches(*elem);
    }
#endif

    template <int k>
    void fill_cache(ElementAccessor &elem, const int j);
    ///@}



    template <int k = dim>
    Size get_num_points() const
    {
        return cacheutils::extract_sub_elements_data<k>(quad_all_sub_elems_).get_num_points();
    }

public:
    void print_info(LogStream &out) const
    {
        // TODO (pauletti, Apr 24, 2015): line below should be implemented
        //quad_all_sub_elems_.print_info(out);
    }
    std::shared_ptr<const GridType> get_grid() const;

//    const TensorProductArray<dim> &get_lengths() const;

private:
    std::shared_ptr<GridType> grid_;

    std::array<ValueFlags, dim + 1> flags_;

protected:
    QuadList<dim> quad_all_sub_elems_;
};

IGA_NAMESPACE_CLOSE

#endif /* GRID_ELEMENT_HANDLER_H_ */
