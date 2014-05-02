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

#ifndef SPACE_SPEC_H_
#define SPACE_SPEC_H_

#include <igatools/base/config.h>
#include <igatools/utils/cartesian_product_array.h>
#include <igatools/utils/static_multi_array.h>
#include <igatools/geometry/cartesian_grid.h>
#include <algorithm>
IGA_NAMESPACE_OPEN

/**
 * @brief Spline space specification class
 *
 *
 * @author pauletti, 2013-2014
 *
 */
template<int dim, int range = 1, int rank = 1>
class SpaceSpec
{

	template<class T>
	using ComponentContainer = StaticMultiArray<T,range,rank>;

	static const int n_components = ComponentContainer<int>::n_entries;
private:
    using Grid = CartesianGrid<dim>;

public:
    using BoundaryKnots = std::array<CartesianProductArray<Real,2>, dim>;
    using Degrees  = TensorIndex<dim>;
    using Multiplicity = CartesianProductArray<Size, dim>;

    using DegreeTable = ComponentContainer<Degrees>;
    using MultiplicityTable = ComponentContainer<Multiplicity>;
    using BoundaryKnotsTable = ComponentContainer<BoundaryKnots>;

    // For the boundary kntos types
    // end-points interpolatory (open knot)
    // periodic
    enum class EndBehaviour {interpolatory, periodic};

    // For the interior multiplicities
    // maximum regularity
    // minimul regularity discontinous
    enum class InteriorReg {maximum, minimun};

public:
    /**
     * Most general constructor
     */
    explicit SpaceSpec(std::shared_ptr<const Grid> knots,
    		std::shared_ptr<const MultiplicityTable> interior_mult,
    		BoundaryKnotsTable &boundary_knots,
    		const DegreeTable &deg);

//    explicit SpaceSpec(std::shared_ptr<const Grid> knots,
//    		MultiplicityTable &interior_mult,
//    		const EndBehaviour boundary_knots,
//    		const DegreeTable &deg);

    explicit SpaceSpec(std::shared_ptr<const Grid> knots,
        		const InteriorReg interior_mult,
        		BoundaryKnotsTable &boundary_knots,
        		const DegreeTable &deg)
    :SpaceSpec(knots, fill_max_regularity(knots) ,boundary_knots, deg)
    {}

//    explicit SpaceSpec(std::shared_ptr<const Grid> knots,
//            		const InteriorReg interior_mult,
//            		const EndBehaviour boundary_knots,
//            		const DegreeTable &deg);


    const DegreeTable &get_degree() const
    {return deg_;}

private:
    /**
     * Fill the multiplicy for the maximum possible regularity
     *  of the given number of knots
     */
    std::shared_ptr<MultiplicityTable> fill_max_regularity(std::shared_ptr<const Grid> grid);
//    /**
//     * Fill the multiplicy for the maximum possible regularity
//     *  of the given number of knots
//     */
//    void fill_max_regularity(const int_array<dim> degree);
//
//    /**
//     * Computes the cumulative multiplicity
//     */
//    parent_t accumulate() ;
public:

    MultiplicityTable compute_index_space_offset() const;

    void print_info(LogStream &out)
    {
    	out << "Knots without repetition:\n";
    	grid_->print_info(out);
    	out << "Repeated knots:\n";
    	for(const auto &v : rep_knots_)
    		v.print_info(out);
    	out << "Interior multiplicities:\n";
    	for(const auto &v : *interior_mult_)
    		v.print_info(out);
    	out << "Boundary knots:\n";
    	for(const auto &v : boundary_knots_)
    		for(const auto &w : v)
    			w.print_info(out);

    	out << "Degrees:\n";
    	deg_.print_info(out);
    }

private:
    std::shared_ptr<const Grid> grid_;
    ComponentContainer<CartesianProductArray<Real, dim>> rep_knots_;
    std::shared_ptr<const MultiplicityTable> interior_mult_;
    MultiplicityTable accumulated_mult_;
    BoundaryKnotsTable boundary_knots_;
    DegreeTable deg_;
};

IGA_NAMESPACE_CLOSE

#endif
