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

#ifndef UNIT_ELEMENT_H_
#define UNIT_ELEMENT_H_

#include <igatools/base/config.h>
#include <igatools/base/tensor.h>
#include <igatools/base/array_utils.h>

IGA_NAMESPACE_OPEN

template <int> struct UnitElement;

constexpr int skel_size(int dim, int k)
{
    return dim==k ? 1 :
           (((k==-1)||(k>dim)) ? 0 : (2*skel_size(dim-1, k) + skel_size(dim-1, k-1)));
}


template <int dim, int k>
EnableIf< (dim==0) || (k<0),
          std::array<typename UnitElement<dim>::template SubElement<k>, skel_size(dim, k)>>
                  fill_cube_elements()
{
    std::array<typename UnitElement<dim>::template SubElement<k>, skel_size(dim, k)> res;
    return res;
}



template <int dim, int k>
EnableIf< (dim==k) && (k>0),
          std::array<typename UnitElement<dim>::template SubElement<k>, skel_size(dim, k)>>
                  fill_cube_elements()
{
    std::array<typename UnitElement<dim>::template SubElement<k>, skel_size(dim, k)> res;
    res[0].active_directions = sequence<k>();
    return res;
}


template <int dim, int k>
EnableIf< (dim>k)  &&(k>=0),
          std::array<typename UnitElement<dim>::template SubElement<k>, skel_size(dim, k)>>
                  fill_cube_elements()
{
    std::array<typename UnitElement<dim>::template SubElement<k>, skel_size(dim, k)> elements;

    auto sub_elems_1 = fill_cube_elements<dim-1, k>();
    auto sub_elems_0 = fill_cube_elements<dim-1, k-1>();

    auto elem = elements.begin();

    for (auto &sub_elem_0 : sub_elems_0)
    {
        auto &sub_dirs_0 = sub_elem_0.constant_directions;
        auto &dirs       = elem->constant_directions;
        std::copy(sub_dirs_0.begin(), sub_dirs_0.end(), dirs.begin());
        elem->constant_values = sub_elem_0.constant_values;
        ++elem;
    }

    for (int j = 0; j<2; ++j)
        for (auto &sub_elem_1 : sub_elems_1)
        {
            auto &sub_dirs_1 = sub_elem_1.constant_directions;
            auto &sub_values_1 = sub_elem_1.constant_values;
            {
                auto &dirs       = elem->constant_directions;
                auto &values       = elem->constant_values;
                std::copy(sub_dirs_1.begin(), sub_dirs_1.end(), dirs.begin());
                dirs[dim - k -1] = dim-1;
                std::copy(sub_values_1.begin(), sub_values_1.end(), values.begin());
                values[dim - k -1] = j;
                ++elem;
            }
        }

    for (auto &elem : elements)
    {
        auto all = sequence<dim>();
        std::set_difference(all.begin(), all.end(),
                            elem.constant_directions.begin(),
                            elem.constant_directions.end(),
                            elem.active_directions.begin());
    }

    return elements;
}



template<int dim, std::size_t... I>
auto tuple_of_elements(std::index_sequence<I...>)
-> decltype(std::make_tuple(fill_cube_elements<dim, I>() ...))
{
    return std::make_tuple(fill_cube_elements<dim, I>() ...);
}



template<std::size_t dim, typename Indices = std::make_index_sequence<dim+1>>
auto construct_cube_elements()
-> decltype(tuple_of_elements<dim>(Indices()))
{
    return tuple_of_elements<dim>(Indices());
}


template<int dim>
using AllCubeElements = decltype(tuple_of_elements<dim>(std::make_index_sequence<dim+1>()));


/**
 * @brief This class provides dimension independent information
 * of all topological structures that make up the hypercube
 * of dimension dim, which is the topological structure of the isogeometric
 * "element"
 *
 */
template <int dim>
struct UnitElement
{
    /**
     * Indices for the coordinate directions active on the UnitElement.
     *
     * Its main use to allow to use the range-based for loop
     * \code{.cpp}
       for (const auto &dir : active_directions)
       {
          // do something
       }
       \endcode
     * instead of the traditional for loop
     * \code{.cpp}
       for (int dir = 0  ; dir < dim_ ; ++dir)
       {
          // do something
       }
       \endcode
     *
     */
    static const std::array<Size,dim> active_directions;


    /**
     * Number of elements of dimension k=0,...,dim in the
     * hyper-cube of dimension dim
     */
    static const std::array<Size, dim + 1> sub_elements_size;

    /**
     * Element of dimension <tt>k</tt> in a cube of dimension <tt>dim</tt>.
     */
    template<int k>
    struct SubElement
    {
        SubElement() = default;

        std::array<Size, dim - k> constant_directions;
        std::array<Size, dim - k> constant_values;
        std::array<Size, k>       active_directions;
    };

    /**
     * This tuple of size dim+1 provides the caracterization of all
     * j dimensional skeleton of the unit cube
     */
    static const AllCubeElements<dim> all_elems;

    template<int k>
    static constexpr Size
    num_elem()
    {
        return skel_size(dim, k);//sub_elements_size[k];
    }

    template<int k>
    static const SubElement<k> &get_elem(const int j)
    {
        //TODO: put assetion on j in proper range
        return (std::get<k>(all_elems)[j]);
    }

    template<int k>
    static constexpr std::array<Index,num_elem<k>()>
    elems_ids()
    {
        return sequence<num_elem<k>()>();
    }

    static constexpr auto n_faces = num_elem<dim-1>();

};

IGA_NAMESPACE_CLOSE

#endif /* UNIT_ELEMENT_H_ */
