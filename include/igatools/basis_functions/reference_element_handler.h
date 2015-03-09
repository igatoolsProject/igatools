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

#ifndef REFERENCE_ELEMENT_HANDLER_H_
#define REFERENCE_ELEMENT_HANDLER_H_

#include <igatools/base/config.h>
#include <igatools/base/flags_handler.h>
#include <igatools/base/evaluation_points.h>

#include <igatools/basis_functions/element_handler.h>
#include <igatools/geometry/grid_element_handler.h>
#include <igatools/basis_functions/reference_space.h>


IGA_NAMESPACE_OPEN

template<int dim, int range = 1, int rank = 1>
class ReferenceElementHandler
    : public ElementHandler<ReferenceSpace<dim,range,rank>>
{
public:
    using Space = ReferenceSpace<dim,range,rank>;
    using ElementIterator = typename Space::ElementIterator;
    using ElementAccessor = typename Space::ElementAccessor;

    static const int l = iga::max(0, dim-num_sub_elem);

    using v2 = typename seq<Int, l, dim>::type;
    using topology_variant = typename boost::make_variant_over<v2>::type;


    using v3 = typename seq<Quadrature, l, dim>::type;
    using eval_pts_variant = typename boost::make_variant_over<v3>::type;


    static std::shared_ptr<ReferenceElementHandler<dim,range,rank> >
    create(std::shared_ptr<const Space> space);

    /** @name Constructors */
    ///@{
    /**
     * Default constructor. Not allowed to be used.
     */
    ReferenceElementHandler() = delete;

    /**
     * Destructor.
     */
    virtual ~ReferenceElementHandler() = default;

    /**
     * Copy constructor. Not allowed to be used.
     */
    ReferenceElementHandler(const ReferenceElementHandler<dim,range,rank> &elem_handler) = delete;

    /**
     * Move constructor. Not allowed to be used.
     */
    ReferenceElementHandler(ReferenceElementHandler<dim,range,rank> &&elem_handler) = delete;
    ///@}


    /**
     * @name Reset functions
     */
    ///@{
    /**
     * Resets all the internal data in order to use the
     * same quadrature scheme for each active element of the space.
     */
    void reset(const ValueFlags &flag, const eval_pts_variant &quad);

    /**
     * Resets all the internal data in order to use the
     * same quadrature scheme for the single element of the space with ID specified by
     * the input parameter <tt>elem_flat_id</tt>.
     */
    void reset_one_element(
        const ValueFlags &flag,
        const eval_pts_variant &eval_points,
        const int elem_flat_id);

    /**
     * Resets all the internal data in order to use the
     * same quadrature scheme for the elements of the space with ID specified by
     * the input parameter <tt>elements_flat_id</tt>.
     *
     * @note This function is pure virtual and must be implemented in the class that are derived
     * from ReferenceElementHandler.
     */
    virtual void reset_selected_elements(
        const ValueFlags &flag,
        const eval_pts_variant &eval_points,
        const vector<int> elements_flat_id) = 0;
    ///@}


protected:
    ReferenceElementHandler(std::shared_ptr<const Space> space);


public:

    /**
     * @name Init functions
     */
    ///@{
    virtual void init_cache(ElementAccessor &elem, const topology_variant &topology) = 0;

    template <int k>
    void init_cache(ElementAccessor &elem)
    {
        this->init_cache(elem,Int<k>());
    }

    ///@}
    /**
     * @name Fill functions
     */
    ///@{
    virtual void fill_cache(ElementAccessor &elem, const topology_variant &topology, const int j) = 0;

    template<int k>
    void fill_cache(ElementAccessor &elem, const int j)
    {
        this->fill_cache(elem,Int<k>(),j);
    }

    ///@}


    virtual void print_info(LogStream &out) const = 0;

    template <int k = dim>
    Size get_num_points() const
    {
        return grid_handler_.template get_num_points<k>();
    }

protected:
    GridElementHandler<dim> grid_handler_;

private:
    std::shared_ptr<const Space> space_;

public:
    /**
     * Returns the const reference of the GridElementHandler used by the current ReferenceElementHandler.
     * @return
     */
    const GridElementHandler<dim> &get_grid_handler() const;

    /**
     * Returns the ReferenceSpace associated to the current ReferenceElementHandler (const version).
     */
    std::shared_ptr<const Space> get_space() const;
};



IGA_NAMESPACE_CLOSE


#endif // REFERENCE_ELEMENT_HANDLER_H_

