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


#include <igatools/basis_functions/spline_space.h>
#include <igatools/base/array_utils.h>

using std::shared_ptr;
using std::make_shared;
using std::const_pointer_cast;

IGA_NAMESPACE_OPEN

template<int dim, int range, int rank>
const Size SplineSpace<dim, range, rank>::n_components;


template<int dim, int range, int rank>
const std::array<Size, SplineSpace<dim, range, rank>::n_components>
SplineSpace<dim, range, rank>::components =
    sequence<SplineSpace<dim, range, rank>::n_components>();

template<int dim, int range, int rank>
SplineSpace<dim, range, rank>::
SplineSpace(const DegreeTable &deg,
            std::shared_ptr<GridType> knots,
            shared_ptr<const MultiplicityTable> interior_mult,
            const PeriodicTable &periodic)
    :
    GridWrapper<CartesianGrid<dim>>(knots),
                                 interior_mult_(interior_mult),
                                 deg_(deg),
                                 periodic_(periodic)
{
    this->init();
#if 0
    // create a signal and a connection for the grid refinement
    this->connect_refinement_h_function(
        std::bind(&SplineSpace<dim,range,rank>::refine_h_after_grid_refinement, this,
                  std::placeholders::_1,std::placeholders::_2));
#endif
}


template<int dim, int range, int rank>
std::shared_ptr<SplineSpace<dim,range,rank> >
SplineSpace<dim, range, rank>::
create(const DegreeTable &deg,
       std::shared_ptr<GridType> knots,
       std::shared_ptr<const MultiplicityTable> interior_mult,
       const PeriodicTable &periodic)
{
    using SpSpace = SplineSpace<dim,range,rank>;
    auto sp = std::shared_ptr<SpSpace>(new SpSpace(deg, knots,interior_mult,periodic));
    Assert(sp != nullptr, ExcNullPtr());

    sp->create_connection_for_h_refinement(sp);

    return sp;
}



template<int dim, int range, int rank>
void
SplineSpace<dim, range, rank>::
init()
{
#ifndef NDEBUG
    auto const knots_size = this->get_grid()->get_num_knots_dim();
    for (const auto comp : components)
        for (const auto j : Topology::active_directions)
        {
            const auto deg = deg_[comp][j];
            const auto order = deg + 1;
            const auto &mult = (*interior_mult_)[comp].get_data_direction(j);
            Assert(mult.size() == knots_size[j]-2,
                   ExcMessage("Interior multiplicity size does not match the grid"));
            if (!mult.empty())
            {
                auto result = std::minmax_element(mult.begin(), mult.end());
                Assert((*result.first > 0) && (*result.second <= order),
                       ExcMessage("multiplicity values not between 0 and p+1"));
            }
        }
#endif

    // Determine the dimensionality of the spline space
    typename SpaceDimensionTable::base_t n_basis;
    for (const auto iComp : components)
        for (const auto dir : Topology::active_directions)
        {
            const auto deg = deg_[iComp][dir];
            const auto &mult = (*interior_mult_)[iComp].get_data_direction(dir);

            Index size = periodic_[iComp][dir] ? 0 : deg + 1;
            for (auto &n: mult)
                size += n;
            n_basis[iComp][dir] = size;
        }
    space_dim_ = n_basis;

#ifndef NDEBUG
    for (const auto comp : components)
        for (const auto dir : Topology::active_directions)
            if (periodic_[comp][dir])
            {
                const auto deg = deg_[comp][dir];
                const auto order = deg + 1;
                Assert(n_basis[comp][dir]>order,
                       ExcMessage("Not enough basis functions"));
            }
#endif
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


    //------------------------------------------------------------------------------------------
    // check if the multiplicity of the new knot lines is compatible with the minimum degree of the space --- begin
    std::array<int,dim> mult_new_knot_lines;
    for (int dir = 0 ; dir < dim ; ++dir)
        mult_new_knot_lines[dir] = 1;

#ifndef NDEBUG
    const auto &interior_mult = const_cast<MultiplicityTable &>(*this->get_interior_mult());
    for (int dir = 0; dir < dim; ++dir)
    {
        int min_degree = std::numeric_limits<int>::max();
        for (int comp_id : interior_mult.get_active_components_id())
        {
            const auto &degree_comp = this->get_degree()[comp_id];
            min_degree = std::min(min_degree,degree_comp[dir]);

            Assert(mult_new_knot_lines[dir] > 0 && mult_new_knot_lines[dir] <= min_degree,
                   ExcIndexRange(mult_new_knot_lines[dir],1,min_degree+1))
        }
    }
#endif
    // check if the multiplicity of the new knot lines is compatible with the minimum degree of the space --- end
    //------------------------------------------------------------------------------------------


    shared_ptr<const MultiplicityTable> interior_mult_prev_refinement =
        make_shared<const MultiplicityTable>(MultiplicityTable(*this->get_interior_mult()));

    spline_space_previous_refinement_ =
        make_shared<const SplineSpace<dim,range,rank> >(
            SplineSpace<dim,range,rank>(
                this->get_degree(),
                grid_pre_refinement,
                interior_mult_prev_refinement));

    for (const auto direction_id : Topology::active_directions)
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

                vector<int> mult_new(n_extra_multiplicities,mult_new_knot_lines[direction_id]);
                for (const int &m : mult_old)
                {
                    mult_new.push_back(m); // adding the old multiplicity value

                    mult_new.insert(mult_new.end(),n_extra_multiplicities,mult_new_knot_lines[direction_id]); // adding the new multiplicity values
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
compute_knots_with_repetition(const EndBehaviourTable &ends,
                              const BoundaryKnotsTable &boundary_knots) const
-> KnotsTable
{

#ifndef NDEBUG
    for (auto iComp : components)
    {
        for (const int j : Topology::active_directions)
        {
            const auto &l_knots = boundary_knots[iComp][j].get_data_direction(0);
            const auto &r_knots = boundary_knots[iComp][j].get_data_direction(1);

            if (periodic_[iComp][j])
            {
                Assert(ends[iComp][j] == BasisEndBehaviour::periodic,
                       ExcMessage("Periodic inconsistency"));
                Assert(l_knots.size()==0,
                       ExcMessage("Periodic inconsistency"));
                Assert(l_knots.size()==0,
                       ExcMessage("Periodic inconsistency"));

            }
            else
            {
                if (ends[iComp][j] == BasisEndBehaviour::interpolatory)
                {
                    Assert(l_knots.size()==0,
                           ExcMessage("Interpolatory inconsistency"));
                    Assert(l_knots.size()==0,
                           ExcMessage("Interpolatory inconsistency"));
                }
                if (ends[iComp][j] == BasisEndBehaviour::end_knots)
                {
                    const auto &knots = this->get_grid()->get_knot_coordinates(j);
                    const int m = deg_[iComp][j] + 1;
                    Assert(l_knots.size() == m,
                           ExcMessage("Wrong number of boundary knots"));
                    Assert(l_knots.size() == m,
                           ExcMessage("Wrong number of boundary knots"));
                    Assert(knots.front() >= l_knots.back(),
                           ExcMessage("Boundary knots should be smaller or equal a"));
                    Assert(knots.back() <= r_knots.front(),
                           ExcMessage("Boundary knots should be greater or equal b"));
                    Assert(std::is_sorted(l_knots.begin(), l_knots.end()),
                           ExcMessage("Boundary knots is not sorted"));
                    Assert(std::is_sorted(r_knots.begin(), r_knots.end()),
                           ExcMessage("Boundary knots is not sorted"));
                }
            }
        }
    }
#endif

    KnotsTable result;

    for (int comp = 0; comp < n_components; ++comp)
    {
        for (const auto dir : Topology::active_directions)
        {
            const auto deg = deg_[comp][dir];
            const auto order = deg + 1;
            const auto &knots = this->get_grid()->get_knot_coordinates(dir);
            const auto &mult  = (*interior_mult_)[comp].get_data_direction(dir);

            int size = 2 * order;
            const int m = order;
            int K = 0;
            for (auto &n: mult)
                K += n;
            size += K;

            vector<Real> rep_knots(size);

            auto rep_it = rep_knots.begin() + m;
            auto m_it = mult.begin();
            auto k_it = ++knots.begin();
            auto end = mult.end();
            for (; m_it !=end; ++m_it, ++k_it)
            {
                for (int iMult = 0; iMult < *m_it; ++iMult, ++rep_it)
                    *rep_it = *k_it;
            }


            if (periodic_[comp][dir])
            {
                const Real a = knots.front();
                const Real b = knots.back();
                const Real L = b - a;
                for (int i=0; i<m; ++i)
                {
                    rep_knots[i] = rep_knots[K+i] - L;
                    rep_knots[K+m+i] = rep_knots[m+i] + L;
                }
            }
            else
            {
                const auto &endb = ends[comp][dir];
                switch (endb)
                {
                    case BasisEndBehaviour::interpolatory:
                    {
                        const Real a = knots.front();
                        const Real b = knots.back();

                        for (int i=0; i<m; ++i)
                        {
                            rep_knots[i] = a;
                            rep_knots[m+K+i] = b;
                        }
                    }
                    break;
                    case BasisEndBehaviour::end_knots:
                    {
                        const auto &left_knts = boundary_knots[comp][dir].get_data_direction(0);
                        const auto &right_knts = boundary_knots[comp][dir].get_data_direction(1);

                        for (int i=0; i<m; ++i)
                        {
                            rep_knots[i]     = left_knts[i];
                            rep_knots[K+m+i] = right_knts[i];
                        }
                    }
                    break;
                    case BasisEndBehaviour::periodic:
                        Assert(false, ExcMessage("Impossible"));
                }
            }

            result[comp].copy_data_direction(dir,rep_knots);
        }
    }

    return result;
}



template<int dim, int range, int rank>
template<int k>
auto
SplineSpace<dim, range, rank>::
get_sub_space_mult(const Index sub_elem_id) const
-> std::shared_ptr<typename SubSpace<k>::MultiplicityTable >
{
    using SubMultT = typename SubSpace<k>::MultiplicityTable;
    const auto &v_mult = *interior_mult_;

    auto &k_elem = UnitElement<dim>::template get_elem<k>(sub_elem_id);
    const auto &active_dirs = k_elem.active_directions;

    auto sub_mult_ptr = make_shared<SubMultT> (v_mult.get_comp_map());
    auto &sub_mult = *sub_mult_ptr;
    for (int comp : sub_mult.get_active_components_id())
    {
        for (int j=0; j<k; ++j)
            sub_mult[comp].copy_data_direction(j, v_mult[comp].get_data_direction(active_dirs[j]));
    }
    return sub_mult_ptr;
}



template<int dim, int range, int rank>
template<int k>
auto
SplineSpace<dim, range, rank>::
get_sub_space_degree(const Index sub_elem_id) const
-> typename SubSpace<k>::DegreeTable
{
    using SubDegreeT = typename SubSpace<k>::DegreeTable;
    auto &k_elem = UnitElement<dim>::template get_elem<k>(sub_elem_id);
    const auto &active_dirs = k_elem.active_directions;

    SubDegreeT sub_degree(deg_.get_comp_map());

    for (int comp : sub_degree.get_active_components_id())
        for (int j=0; j<k; ++j)
            sub_degree[comp][j] = deg_[comp][active_dirs[j]];

    return sub_degree;
}



template<int dim, int range, int rank>
template<int k>
auto
SplineSpace<dim, range, rank>::
get_sub_space_periodicity(const Index sub_elem_id) const
-> typename SubSpace<k>::PeriodicTable
{
    using SubPeriodicT = typename SubSpace<k>::PeriodicTable;
    auto &k_elem = UnitElement<dim>::template get_elem<k>(sub_elem_id);
    const auto &active_dirs = k_elem.active_directions;

    SubPeriodicT sub_periodic(periodic_.get_comp_map());

    for (int comp : periodic_.get_active_components_id())
        for (int j=0; j<k; ++j)
            sub_periodic[comp][j] = periodic_[comp][active_dirs[j]];

    return sub_periodic;
}



template<int dim, int range, int rank>
auto SplineSpace<dim, range, rank>::
accumulated_interior_multiplicities() const -> MultiplicityTable
{
    MultiplicityTable result;
    for (int iComp = 0; iComp < n_components; ++iComp)
    {
        for (const auto j : Topology::active_directions)
        {
            // Assert(!periodic_[iComp][j], ExcMessage("periodic needs to be implemented"));
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
get_multiplicity_from_regularity(const InteriorReg reg,
                                 const DegreeTable &deg,
                                 const TensorSize<dim> &n_elem)
-> std::shared_ptr<MultiplicityTable>
{
    auto  res = std::make_shared<MultiplicityTable>(deg.get_comp_map());


    for (int comp : res->get_active_components_id())
        for (const auto dir : Topology::active_directions)
        {
            int val;
            switch (reg)
            {
                case InteriorReg::maximum:
                    val = 1;
                    break;
                case InteriorReg::minimun:
                    val = deg[comp][dir] + 1;
                    break;
            }

            const auto size = n_elem[dir]-1;

            vector<Size> mult(size);
            if (size>0)
                (*res)[comp].copy_data_direction(dir, vector<Size>(size, val));
        }
    return res;
}

template<int dim, int range, int rank>
void
SplineSpace<dim, range, rank>::
create_connection_for_h_refinement(std::shared_ptr<SplineSpace<dim,range,rank>> space)
{
    Assert(space != nullptr, ExcNullPtr());

    auto refinement_func_spline_space =
        std::bind(&SplineSpace<dim,range,rank>::refine_h_after_grid_refinement,
                  space.get(),
                  std::placeholders::_1,
                  std::placeholders::_2);

    using SlotType = typename CartesianGrid<dim>::SignalRefineSlot;
    this->connect_refinement_h_function(
        SlotType(refinement_func_spline_space).track_foreign(space));

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

}

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/spline_space.inst>
