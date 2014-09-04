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


#include <igatools/basis_functions/spline_space.h>
#include <igatools/base/array_utils.h>

using std::array;
using std::shared_ptr;
using std::make_shared;
using std::const_pointer_cast;

IGA_NAMESPACE_OPEN

template<int dim, int range, int rank>
const std::array<int, dim>
SplineSpace<dim, range, rank>::dims = sequence<dim>();

template<int dim, int range, int rank>
const std::array<int, SplineSpace<dim, range, rank>::n_components>
SplineSpace<dim, range, rank>::components = sequence<SplineSpace<dim, range, rank>::n_components>();

template<int dim, int range, int rank>
SplineSpace<dim, range, rank>::
SplineSpace(const DegreeTable &deg,
            std::shared_ptr<GridType> knots,
            shared_ptr<const MultiplicityTable> interior_mult,
            const PeriodicTable &periodic)
    :
    GridSpace(knots),
    interior_mult_(interior_mult),
    deg_(deg),
    periodic_(periodic)
{

    //-------------------
    const auto comp_map = interior_mult_->get_comp_map();

    end_behaviour_ = EndBehaviourTable(comp_map);
    for (const auto &comp : end_behaviour_.get_active_components_id())
        end_behaviour_[comp] = filled_array<EndBehaviour,dim>(EndBehaviour::interpolatory);
    //-------------------


    this->init();

    // create a signal and a connection for the grid refinement
    this->connect_refinement_h_function(
        std::bind(&SplineSpace<dim,range,rank>::refine_h_after_grid_refinement, this,
                  std::placeholders::_1,std::placeholders::_2));

}

template<int dim, int range, int rank>
void
SplineSpace<dim, range, rank>::
init()
{
#ifndef NDEBUG
    auto const knots_size = this->get_grid()->get_num_knots_dim();
    for (int iComp = 0; iComp < n_components; ++iComp)
    {
        for (int j = 0; j < dim; ++j)
        {
            const auto deg = deg_[iComp][j];
            const auto order = deg + 1;
            const auto &mult = (*interior_mult_)[iComp].get_data_direction(j);
            Assert(mult.size() == knots_size[j]-2,
                   ExcMessage("Interior multiplicity size does not match the grid"));
            if (!mult.empty())
            {
                auto result = std::minmax_element(mult.begin(), mult.end());
                Assert((*result.first > 0) && (*result.second <= order),
                       ExcMessage("multiplicity values not between 0 and p+1"));
            }
        }
    }
#endif
    Index total_dim = 0;
    for (int iComp = 0; iComp < n_components; ++iComp)
    {
        for (int j = 0; j < dim; ++j)
        {
            const auto deg = deg_[iComp][j];
            const auto &mult = (*interior_mult_)[iComp].get_data_direction(j);

            Index size = periodic_[iComp][j]? 0 : deg + 1;

            for (auto &n: mult)
                size += n;
            space_dim_[iComp][j] = size;
        }
        space_dim_.comp_dimension[iComp] = space_dim_[iComp].flat_size();
        total_dim += space_dim_.comp_dimension[iComp];
    }
    space_dim_.total_dimension = total_dim;
}



template<int dim, int range, int rank>
void
SplineSpace<dim, range, rank>::
refine_h_after_grid_refinement(
    const std::array<bool,dim> &refinement_directions,
    const GridType &grid_old)
{
    Assert(this->get_grid()->get_grid_pre_refinement()!=nullptr,ExcNullPtr());
    auto grid_pre_refinement = const_pointer_cast<CartesianGrid<dim>>(this->get_grid()->get_grid_pre_refinement());

    shared_ptr<const MultiplicityTable> interior_mult_prev_refinement =
        make_shared<const MultiplicityTable>(MultiplicityTable(*this->get_interior_mult()));

    spline_space_previous_refinement_ =
        make_shared<const SplineSpace<dim,range,rank> >(
            SplineSpace<dim,range,rank>(
                this->get_degree(),
                grid_pre_refinement,
                interior_mult_prev_refinement));

    for (int direction_id = 0; direction_id < dim; ++direction_id)
    {
        if (refinement_directions[direction_id])
        {
            // knots in the refined grid along the selected direction
            vector<Real> knots_new = this->get_grid()->get_knot_coordinates(direction_id);
            const int n_elements_new = knots_new.size() - 1;

            // knots in the original (unrefined) grid along the selected direction
            vector<Real> knots_old = grid_old.get_knot_coordinates(direction_id);
            const int n_elements_old = knots_old.size() - 1;

            Assert(n_elements_new % n_elements_old == 0,
                   ExcMessage("The number of new elements along one direction is not an integer multiple of the number of old elements."))
            const int refine_factor = n_elements_new / n_elements_old;
            const int n_extra_multiplicities = refine_factor - 1;

            vector<Real> knots_added(knots_new.size());

            // find the knots in the refined grid that are not present in the old grid
            auto it = std::set_difference(
                          knots_new.begin(),knots_new.end(),
                          knots_old.begin(),knots_old.end(),
                          knots_added.begin());

            knots_added.resize(it-knots_added.begin());

            auto &interior_mult = const_cast<MultiplicityTable &>(*this->get_interior_mult());
            for (int comp_id : interior_mult.get_active_components_id())
            {
                //--------------------------------------------------------
                // creating the new multiplicity
                const vector<int> &mult_old = interior_mult[comp_id].get_data_direction(direction_id);

                vector<int> mult_new(n_extra_multiplicities,1);
                for (const int &m : mult_old)
                {
                    mult_new.push_back(m); // adding the old multiplicity value

                    mult_new.insert(mult_new.end(),n_extra_multiplicities,1); // adding the new multiplicity values
                }

                interior_mult[comp_id].copy_data_direction(direction_id,mult_new);
                //--------------------------------------------------------
            } // end loop comp_id
        } // end if(refinement_directions[direction_id])
    } // end loop direction_id

    this->init();
}

template<int dim, int range, int rank>
auto
SplineSpace<dim, range, rank>::
compute_knots_with_repetition(const BoundaryKnotsTable &boundary_knots) const
-> KnotsTable
{
#ifndef NDEBUG
    for (int iComp = 0; iComp < n_components; ++iComp)
    {
        for (int j = 0; j < dim; ++j)
        {
            const auto deg = deg_[iComp][j];
            const auto order = deg + 1;
            const auto &knots = this->get_grid()->get_knot_coordinates(j);
            const auto &left_knts = boundary_knots[iComp][j].get_data_direction(0);
            const auto &right_knts = boundary_knots[iComp][j].get_data_direction(1);

            if (periodic_[iComp][j])
            {
                Assert((left_knts.size()==0) && (right_knts.size()==0),
                       ExcMessage("Periodic component has non zero size"));
            }
            else
            {
                Assert((left_knts.size() == order) && (right_knts.size() == order),
                       ExcMessage("Wrong number of boundary knots"));
                Assert(knots.front() >= left_knts.back(),
                       ExcMessage("Boundary knots should be smaller or equal a"));
                Assert(knots.back() <= right_knts.front(),
                       ExcMessage("Boundary knots should be greater or equal b"));
                Assert(std::is_sorted(left_knts.begin(), left_knts.end()),
                       ExcMessage("Boundary knots is not sorted"));
                Assert(std::is_sorted(right_knts.begin(), right_knts.end()),
                       ExcMessage("Boundary knots is not sorted"));
            }
        }
    }
#endif

    KnotsTable result;

    for (int iComp = 0; iComp < n_components; ++iComp)
    {
        for (int j = 0; j < dim; ++j)
        {
            const auto deg = deg_[iComp][j];
            const auto order = deg + 1;
            const auto &knots = this->get_grid()->get_knot_coordinates(j);
            const auto &mult  = (*interior_mult_)[iComp].get_data_direction(j);
            const auto &left_knts = boundary_knots[iComp][j].get_data_direction(0);
            const auto &right_knts = boundary_knots[iComp][j].get_data_direction(1);

            int size = 2 * order;
            for (auto &n: mult)
                size += n;

            vector<Real> rep_knots;
            rep_knots.reserve(size);
            rep_knots.insert(rep_knots.end(), left_knts.begin(), left_knts.end());
            auto m_it = mult.begin();
            auto k_it = ++knots.begin();
            auto end = mult.end();
            for (; m_it !=end; ++m_it, ++k_it)
            {
                for (int iMult = 0; iMult < *m_it; ++iMult)
                    rep_knots.push_back(*k_it);
            }
            rep_knots.insert(rep_knots.end(), right_knts.begin(), right_knts.end());

            result[iComp].copy_data_direction(j,rep_knots);
        }
    }

    return result;
}



template<int dim, int range, int rank>
auto
SplineSpace<dim, range, rank>::
get_face_mult(const Index face_id) const
-> shared_ptr<typename FaceSpace::MultiplicityTable>
{
    const auto &v_mult = *interior_mult_;
    const auto &active_dirs = UnitElement<dim>::face_active_directions[face_id];
    auto f_int_mult = make_shared<typename FaceSpace::MultiplicityTable> (v_mult.get_comp_map());
    auto &f_mult = *f_int_mult;
    for (int comp : f_mult.get_active_components_id())
    {
        for (auto j : FaceSpace::dims)
            f_mult[comp].copy_data_direction(j, v_mult[comp].get_data_direction(active_dirs[j]));
    }
    return f_int_mult;
}



template<int dim, int range, int rank>
auto SplineSpace<dim, range, rank>::
get_face_degree(const Index face_id) const
-> typename FaceSpace::DegreeTable
{
    const auto &active_dirs = UnitElement<dim>::face_active_directions[face_id];
    typename FaceSpace::DegreeTable f_degree(deg_.get_comp_map());
    for (int comp : f_degree.get_active_components_id())
    {
        for (auto j : FaceSpace::dims)
            f_degree[comp][j] = deg_[comp][active_dirs[j]];
    }
    return f_degree;
}



template<int dim, int range, int rank>
auto SplineSpace<dim, range, rank>::
accumulated_interior_multiplicities() const -> MultiplicityTable
{
    MultiplicityTable result;
    for (int iComp = 0; iComp < n_components; ++iComp)
    {
        for (int j = 0; j < dim; ++j)
        {
            Assert(!periodic_[iComp][j], ExcMessage("periodic needs to be implemented"));
            const auto &mult  = (*interior_mult_)[iComp].get_data_direction(j);
            vector<Size> accum_mult;
            const int size = mult.size();
            accum_mult.reserve(size + 1);
            accum_mult.push_back(0);
            for (int i = 0; i < size; ++i)
                accum_mult.push_back(accum_mult[i] + mult[i]);

            result[iComp].copy_data_direction(j, accum_mult);

            //TODO(pauletti, May 3, 2014): write some post assertions
        }
    }
    return result;
}



template<int dim, int range, int rank>
auto
SplineSpace<dim, range, rank>::
fill_max_regularity(const DegreeTable &deg, std::shared_ptr<const GridType> grid) -> std::shared_ptr<MultiplicityTable>
{
    auto  res = std::make_shared<MultiplicityTable>(deg.get_comp_map());

    auto const knots_size = grid->get_num_knots_dim();
    for (int iComp : res->get_active_components_id())
        for (int j = 0; j < dim; ++j)
        {
            const auto size = knots_size[j]-2;
            if (size>0)
                (*res)[iComp].copy_data_direction(j, vector<Size>(size, 1));
        }
    return res;
}


#if 0
template<int dim, int range, int rank>
auto
SplineSpace<dim, range, rank>::
interpolatory_end_knots() const -> BoundaryKnotsTable
{
    BoundaryKnotsTable result(deg_.get_comp_map());

    for (int iComp : result.get_active_components_id())
    {
        BoundaryKnots bdry_knots;
        for (int j = 0; j < dim; ++j)
        {
            const auto &knots = this->get_grid()->get_knot_coordinates(j);
            const auto deg = deg_(iComp)[j];
            const auto order = deg + 1;
            const Real a = knots.front();
            const Real b = knots.back();
            vector<Real> vec_left(order, a);
            vector<Real> vec_right(order, b);
            bdry_knots[j].copy_data_direction(0, vec_left);
            bdry_knots[j].copy_data_direction(1, vec_right);
        }
        result(iComp) = bdry_knots;
    }
    return result;
}
#endif

template<int dim, int range, int rank>
auto
SplineSpace<dim, range, rank>::
interpolatory_end_knots(const int comp_id,const int dir) const -> CartesianProductArray<Real,2>
{
    CartesianProductArray<Real,2> bdry_knots_dir;

    const auto &knots = this->get_grid()->get_knot_coordinates(dir);
    const auto deg = deg_[comp_id][dir];
    const auto order = deg + 1;
    const Real a = knots.front();
    const Real b = knots.back();
    vector<Real> vec_left(order, a);
    vector<Real> vec_right(order, b);
    bdry_knots_dir.copy_data_direction(0, vec_left);
    bdry_knots_dir.copy_data_direction(1, vec_right);

    return bdry_knots_dir;
}

template<int dim, int range, int rank>
auto
SplineSpace<dim, range, rank>::
compute_knots_with_repetition(const EndBehaviourTable &ends) const -> KnotsTable
{
    BoundaryKnotsTable bdry_knots_table(deg_.get_comp_map());
    for (int iComp : bdry_knots_table.get_active_components_id())
    {
        for (int j = 0; j < dim; ++j)
        {
            if (ends[iComp][j] == EndBehaviour::interpolatory)
                bdry_knots_table[iComp][j] = interpolatory_end_knots(iComp,j);
            else
            {
                Assert(false,ExcNotImplemented());
                AssertThrow(false,ExcNotImplemented());
            }
        }
    }
    return compute_knots_with_repetition(bdry_knots_table);
}


template<int dim, int range, int rank>
void
SplineSpace<dim, range, rank>::
print_info(LogStream &out) const
{
    out.begin_item("Knots without repetition:");
    this->get_grid()->print_info(out);
    out.end_item();

    out.begin_item("Degrees:");
    deg_.print_info(out);
    out.end_item();

    out.begin_item("Interior multiplicities:");
    const MultiplicityTable &interior_mult_ref = *interior_mult_;
    for (const auto &v : interior_mult_ref)
        v.print_info(out);
    out.end_item();

    out.begin_item("Dimensionality Table:");
    space_dim_.print_info(out);
    out.end_item();

    out.begin_item("Component Dimension:");
    space_dim_.comp_dimension.print_info(out);
    out.end_item();

    out << "Total Dimension: " << space_dim_.total_dimension << std::endl;
}

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/spline_space.inst>

