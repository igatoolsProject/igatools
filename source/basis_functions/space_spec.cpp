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


#include <igatools/basis_functions/space_spec.h>

using std::vector;
using std::array;
using std::shared_ptr;

IGA_NAMESPACE_OPEN

template<int dim, int range, int rank>
SpaceSpec<dim, range, rank>::
SpaceSpec(std::shared_ptr<const Grid> knots,
		shared_ptr<const MultiplicityTable> interior_mult,
		const DegreeTable &deg,
		const PeriodicTable periodic)
    :
	  grid_(knots),
	  interior_mult_(interior_mult),
	  deg_(deg),
	  periodic_(periodic)
    {
#ifndef NDEBUG
	auto const knots_size = grid_->get_num_knots_dim();
	for (int iComp = 0; iComp < n_components; ++iComp)
	{
		for (int j = 0; j < dim; ++j)
		{
			const auto deg = deg_(iComp)[j];
			const auto order = deg + 1;
			const auto &mult = (*interior_mult_)(iComp).get_data_direction(j);
			Assert(mult.size() == knots_size[j]-2,
					ExcMessage("Interior multiplicity size does not match the grid") );
			auto result = std::minmax_element(mult.begin(), mult.end());
			Assert( (*result.first > 0) && (*result.second <= order),
					ExcMessage("multiplicity values not between 0 and p+1") );
		}
	}
#endif

	for (int iComp = 0; iComp < n_components; ++iComp)
	{
		for (int j = 0; j < dim; ++j)
		{
			const auto deg = deg_(iComp)[j];
			const auto order = deg + 1;
			const auto &mult = (*interior_mult_)(iComp).get_data_direction(j);
			int size = order;
			for (auto &n: mult)
				size += n;
			space_dim_(iComp)[j] = size;
		}
	}
}


template<int dim, int range, int rank>
auto
SpaceSpec<dim, range, rank>::
compute_knots_with_repetition(const BoundaryKnotsTable &boundary_knots)
-> KnotsTable
{
#ifndef NDEBUG
	for (int iComp = 0; iComp < n_components; ++iComp)
	{
		for (int j = 0; j < dim; ++j)
		{
			const auto deg = deg_(iComp)[j];
			const auto order = deg + 1;
			const auto &knots = grid_->get_knot_coordinates(j);
			const auto &left_knts = boundary_knots(iComp)[j].get_data_direction(0);
			const auto &right_knts = boundary_knots(iComp)[j].get_data_direction(1);

			if (periodic_(iComp))
			{
				Assert((left_knts.size()==0) && (right_knts.size()==0),
						ExcMessage("Periodic component has non zero size"));
			}
			else
			{
				Assert((left_knts.size() == order) && (right_knts.size() == order),
						ExcMessage("Wrong number of boundary knots"));
				Assert(knots.front() >= left_knts.back(),
						ExcMessage("Boundary knots should be smaller or equal a") );
				Assert(knots.back() <= right_knts.front(),
						ExcMessage("Boundary knots should be greater or equal b") );
				Assert(std::is_sorted(left_knts.begin(), left_knts.end()),
						ExcMessage("Boundary knots is not sorted") );
				Assert(std::is_sorted(right_knts.begin(), right_knts.end()),
						ExcMessage("Boundary knots is not sorted") );
			}
		}
	}
#endif

	KnotsTable result;

	for (int iComp = 0; iComp < n_components; ++iComp)
	{
		for (int j = 0; j < dim; ++j)
		{
			const auto deg = deg_(iComp)[j];
			const auto order = deg + 1;
			const auto &knots = grid_->get_knot_coordinates(j);
			const auto &mult  = (*interior_mult_)(iComp).get_data_direction(j);
			const auto &left_knts = boundary_knots(iComp)[j].get_data_direction(0);
			const auto &right_knts = boundary_knots(iComp)[j].get_data_direction(1);

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

			result(iComp).copy_data_direction(j,rep_knots);
		}
	}

	return result;
}

//template<int dim, int range, int rank>
//SpaceSpec<dim, range, rank>::
//SpaceSpec(std::shared_ptr<const Grid> knots,
//             const DegreeTable &deg,
//             const bool max_reg)
//:
//parent_t::StaticMultiArray(T(knots->get_num_knots_dim())),
//grid_(knots),
//deg_(deg)
//{
//    fill_max_regularity();
//}

//template <int dim>
//auto
//Multiplicity<dim>::
//accumulate() -> parent_t
//{
//    const TensorSize<dim> size = this->tensor_size();
//    parent_t result(size);
//
//    for (int i = 0; i < dim; ++i)
//    {
//        result.entry(i, 0) =  this->data_[i][0];
//
//        const Size size_i = size(i);
//        for (int k = 1 ; k < size_i ; ++k)
//            result.entry(i, k) = result.entry(i, k-1) + this->data_[i][k];
//    }
//
//    return result;
//}
//

template<int dim, int range, int rank>
auto SpaceSpec<dim, range, rank>::
compute_elements_index_space_mark() const -> MultiplicityTable
{
	MultiplicityTable result;
	for (int iComp = 0; iComp < n_components; ++iComp)
	{
		Assert(!periodic_(iComp), ExcMessage("periodic needs to be implemented"))
		for (int j = 0; j < dim; ++j)
		{
			const auto deg = deg_(iComp)[j];
			const auto order = deg + 1;
			const auto &mult  = (*interior_mult_)(iComp).get_data_direction(j);
			std::vector<Size> accum_mult;
			const int size = mult.size();
			accum_mult.reserve(size + 1);
			accum_mult.push_back(0);
			accum_mult.push_back(order);
			for (int i = 0; i < size-1; ++i)
				accum_mult.push_back(accum_mult[i+1] + mult[i]);

			result(iComp).copy_data_direction(j, accum_mult);

			//TODO(pauletti, May 3, 2014): write some post assertions
		}
	}
	return result;
}



template<int dim, int range, int rank>
auto
SpaceSpec<dim, range, rank>::
fill_max_regularity(std::shared_ptr<const Grid> grid) -> std::shared_ptr<MultiplicityTable>
{
	auto  res = std::make_shared<MultiplicityTable>();

	auto const knots_size = grid->get_num_knots_dim();
	for (int iComp = 0; iComp < n_components; ++iComp)
		for (int j = 0; j < dim; ++j)
		{
			(*res)(iComp).copy_data_direction(j, vector<Size>(knots_size[j]-2, 1));
		}
	return res;
}



template<int dim, int range, int rank>
auto
SpaceSpec<dim, range, rank>::interpolatory_end_knots() -> BoundaryKnotsTable
{
	BoundaryKnotsTable result;
	for (int iComp = 0; iComp < n_components; ++iComp)
	{
		BoundaryKnots bdry_knots;
		for (int j = 0; j < dim; ++j)
		{
			const auto deg = deg_(iComp)[j];
			const auto order = deg + 1;
			const Real a = 0;
			const Real b = 1;
			std::vector<Real> vec_left(order, a);
			std::vector<Real> vec_right(order, b);
			bdry_knots[j].copy_data_direction(0, vec_left);
			bdry_knots[j].copy_data_direction(1, vec_right);
		}
		result(iComp) = bdry_knots;
	}
	return result;
}



template<int dim, int range, int rank>
void
SpaceSpec<dim, range, rank>::
print_info(LogStream &out)
{
	out << "Knots without repetition:\n";
	grid_->print_info(out);
	out << "Degrees:\n";
	deg_.print_info(out);
	out << std::endl;
	out << "Interior multiplicities:\n";
	for(const auto &v : *interior_mult_)
		v.print_info(out);
	out << "Dimensionality:\n";
	space_dim_.print_info(out);
	out << std::endl;
}

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/space_spec.inst>

