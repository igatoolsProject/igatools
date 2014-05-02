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

	const int n_components = ComponentContainer<int>::n_entries;
private:
    using Grid = CartesianGrid<dim>;

public:
    using BoundaryKnots = std::array<CartesianProductArray<Real,2>, dim>;
    using Degrees  = TensorIndex<dim>;
    using Multiplicity = CartesianProductArray<Size, dim>;

    using DegreeTable = ComponentContainer<Degrees>;
    using MultiplicityTable = ComponentContainer<Multiplicity>;
    using BoundaryKnotsTable = ComponentContainer<BoundaryKnots>;

public:
//    explicit SpaceSpec(parent_t &mult, const DegreeTable &deg)
//    :
//          (mult),
//          deg_(deg)
//        {}
    /**
     * Most general constructor
     */
    explicit SpaceSpec(std::shared_ptr<const Grid> knots,
    		MultiplicityTable &interior_mult,
    		BoundaryKnotsTable &boundary_knots,
    		const DegreeTable &deg)
    :
	  grid_(knots),
	  interior_mult_(interior_mult),
	  boundary_knots_(boundary_knots),
      deg_(deg)
    {
    	auto const knots_size = grid_->get_num_knots_dim();

    	for (int iComp = 0; iComp < n_components; ++iComp)
    	{
    		for (int j = 0; j < dim; ++j)
    		{
    			const auto deg = deg_(iComp)[j];
    			const auto order = deg + 1;
    			const auto &knots = grid_->get_knot_coordinates(j);
    			for(int side=0; side<2; ++side)
    			{
    				Assert(boundary_knots_(iComp)[j].get_data_direction(side).size()
    						== order,
    						ExcMessage("Wrong number of boundary knots"));
    			}
    			const auto a = knots.front();
    			const auto b = knots.back();
    			Assert(a >= boundary_knots_(iComp)[j].get_data_direction(0).back(),
    					ExcMessage("Boundary knots should be smaller or equal a") );
    			Assert(b <= boundary_knots_(iComp)[j].get_data_direction(1).front(),
    					ExcMessage("Boundary knots should be greater or equal b") );

    			//Interior multiplicity check
    			const auto &mult = interior_mult_(iComp).get_data_direction(j);
    			Assert(mult.size() == knots_size[j]-2,
    					ExcMessage("Interior multiplicity size does not match the grid") );
    			auto result = std::minmax_element(mult.begin(), mult.end());
    			Assert( (*result.first > 0) && (*result.second <= order),
    					ExcMessage("multiplicity values not between 0 and p+1") );
    		}
    	}

    	//--------------------------------------------------------------------------------------
    	// filling the knots with repetitions

    	for (int iComp = 0; iComp < n_components; ++iComp)
    	{
    		for (int j = 0; j < dim; ++j)
    		{
    			const auto deg = deg_(iComp)[j];
    			const auto order = deg + 1;
    			const auto &knots = grid_->get_knot_coordinates(j);
    			const auto &mult  = interior_mult_(iComp).get_data_direction(j);
    			const auto &left_knts = boundary_knots_(iComp)[j].get_data_direction(0);
    			const auto &right_knts = boundary_knots_(iComp)[j].get_data_direction(1);

    			int size = 2 * order;
    			for (auto &n: mult)
    				size += n;

    			std::vector<Real> rep_knots;
    			rep_knots.reserve(size);
    			rep_knots.insert(rep_knots.end(), left_knts.begin(), left_knts.end());
    			auto m_it = mult.begin();
    			auto k_it = ++knots.begin();
    			auto end = mult.end();
    			for (;m_it !=end; ++m_it, ++k_it)
    			{
    				for (int iMult = 0; iMult < *m_it; ++iMult)
    					rep_knots.push_back(*k_it);
    			}
    			rep_knots.insert(rep_knots.end(), right_knts.begin(), right_knts.end());

    			rep_knots_(iComp).copy_data_direction(j,rep_knots);
    		}


    	}
    }

    /**
     * Maximun regularity multiplicity vectors associated with
     * the grid and degrees
     */
//    explicit SpaceSpec(std::shared_ptr<const Grid> knots,
//                          const DegreeTable &deg,
//                          const bool max_reg);

    const DegreeTable &get_degree() const
    {return deg_;}

private:
    /**
     * Fill the multiplicy for the maximum possible regularity
     *  of the given number of knots
     */
    void fill_max_regularity();
//
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
    	for(const auto &v : interior_mult_)
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
    MultiplicityTable interior_mult_;
    MultiplicityTable accumulated_mult_;
    BoundaryKnotsTable boundary_knots_;
    DegreeTable deg_;
};

IGA_NAMESPACE_CLOSE

#endif
