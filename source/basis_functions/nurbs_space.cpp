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

#include <igatools/basis_functions/nurbs_space.h>
#include <igatools/basis_functions/space_tools.h>

#include <igatools/geometry/mapping_slice.h>
#include <igatools/base/exceptions.h>

using std::array;
using std::vector;
using std::endl;
using std::shared_ptr;
using std::bind;
using std::placeholders::_1;
using std::placeholders::_2;

IGA_NAMESPACE_OPEN

//TODO(pauletti, Jan 19, 2014): remove the is_bspline from this class

template <int dim_, int range_, int rank_>
const std::array<int, NURBSSpace<dim_, range_, rank_>::n_components>
NURBSSpace<dim_, range_, rank_>::
components = sequence<spline_space_t::n_components>();


template <int dim_, int range_, int rank_>
const std::array<int, NURBSSpace<dim_, range_, rank_>::dim>
NURBSSpace<dim_, range_, rank_>::
dims = sequence<dim>();

template <int dim_, int range_, int rank_>
void
NURBSSpace<dim_, range_, rank_>::
create_refinement_connection()
{
    // create a signal and a connection for the grid refinement
    this->connect_refinement_h_function(
        bind(&self_t::refine_h_weights, this, std::placeholders::_1, std::placeholders::_2));
}

template <int dim_, int range_, int rank_>
NURBSSpace<dim_, range_, rank_>::
NURBSSpace(const int degree,
           shared_ptr< GridType > knots,
           const WeightsTable &weights)
    :
    NURBSSpace(DegreeTable(TensorIndex<dim>(degree)), knots, weights)
{}



template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
create(const int degree,
       shared_ptr< GridType > knots,
       const WeightsTable &weights) -> shared_ptr<self_t>
{
    return shared_ptr<self_t>(new self_t(degree, knots,  weights));
}



template <int dim_, int range_, int rank_>
NURBSSpace<dim_, range_, rank_>::
NURBSSpace(const DegreeTable &degree,
           shared_ptr<GridType> knots,
           const WeightsTable &weights)
    :
    BaseSpace(knots),
    sp_space_(spline_space_t::create(degree,knots)),
    weights_(weights)
{

    create_refinement_connection();
    perform_post_construction_checks();
}



template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
create(const DegreeTable &degree, shared_ptr<GridType> knots,
       const WeightsTable &weights)
-> shared_ptr<self_t>
{
    return shared_ptr<self_t>(new self_t(degree, knots, weights));
}



template <int dim_, int range_, int rank_>
NURBSSpace<dim_, range_, rank_>::
NURBSSpace(const DegreeTable &deg,
           std::shared_ptr<GridType> knots,
           std::shared_ptr<const MultiplicityTable> interior_mult,
           const EndBehaviourTable &ends,
           const WeightsTable &weights)
    :
    BaseSpace(knots),
    sp_space_(spline_space_t::create(deg, knots, interior_mult, ends)),
    weights_(weights)
{
    create_refinement_connection();
    perform_post_construction_checks();
}



template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
create(const DegreeTable &deg,
       std::shared_ptr<GridType> knots,
       std::shared_ptr<const MultiplicityTable> interior_mult,
       const EndBehaviourTable &ends,
       const WeightsTable &weights) -> shared_ptr<self_t>
{
    return shared_ptr<self_t>(new self_t(deg, knots, interior_mult, ends, weights));
}



template <int dim_, int range_, int rank_>
NURBSSpace<dim_, range_, rank_>::
NURBSSpace(std::shared_ptr<spline_space_t> bs_space,
           const WeightsTable &weights)
    :
    BaseSpace(bs_space->get_grid()),
    sp_space_(bs_space),
    weights_(weights)
{
    create_refinement_connection();
    perform_post_construction_checks();
}


template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
create(std::shared_ptr<spline_space_t> bs_space,
       const WeightsTable &weights) -> shared_ptr<self_t>
{
    return shared_ptr<self_t>(new self_t(bs_space, weights));
}

template <int dim_, int range_, int rank_>
void
NURBSSpace<dim_, range_, rank_>::
perform_post_construction_checks() const
{
#ifndef NDEBUG
    // check that the number of weights is equal to the number of basis functions in the space
    for (auto comp : components)
    {
        Assert(sp_space_->get_num_basis(comp) == weights_(comp).flat_size(),
               ExcDimensionMismatch(sp_space_->get_num_basis(comp),weights_(comp).flat_size()));
    }
#endif
}



template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
begin() const -> ElementIterator
{
    return ElementIterator(std::enable_shared_from_this<NURBSSpace<dim_,range_,rank_>>::shared_from_this(), 0);
}



template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
last() const -> ElementIterator
{
    return ElementIterator(std::enable_shared_from_this<NURBSSpace<dim_,range_,rank_>>::shared_from_this(),
                           this->get_grid()->get_num_elements() - 1);
}



template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
end() const -> ElementIterator
{
    return ElementIterator(std::enable_shared_from_this<NURBSSpace<dim_,range_,rank_>>::shared_from_this(),
                           IteratorState::pass_the_end);
}



template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
get_weights() const -> const WeightsTable &
{
    return weights_;
}



template <int dim_, int range_, int rank_>
void
NURBSSpace<dim_, range_, rank_>::
reset_weights(const WeightsTable &weights)
{
    weights_ = weights;
    perform_post_construction_checks();
}


template<int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
get_ref_face_space(const Index face_id,
                   std::vector<Index> &face_to_element_dofs,
                   std::map<int, int> &elem_map) const
-> std::shared_ptr<RefFaceSpace>
{
    auto f_space = sp_space_->get_ref_face_space(face_id, face_to_element_dofs, elem_map);

    // TODO (pauletti, Jun 11, 2014): this should be put and completed in
    // get_face_weigjts()
    const auto &v_weights = weights_;
    //const auto &active_dirs = UnitElement<dim>::face_active_directions[face_id];
    typename RefFaceSpace::WeightsTable f_weights(v_weights.get_comp_map());

    const auto n_basis = f_space->get_num_basis_table();
    for (int comp : f_weights.get_active_components_id())
    {
        f_weights(comp).resize(n_basis(comp),1.0);
        //        for (auto j : RefFaceSpace::dims)
        //            f_weights(comp).copy_data_direction(j, v_weights(comp).get_data_direction(active_dirs[j]));
    }


    return RefFaceSpace::create(f_space, f_weights);
}



template<int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
get_face_space(const Index face_id,
               std::vector<Index> &face_to_element_dofs) const
-> std::shared_ptr<FaceSpace>
{
    auto elem_map = std::make_shared<std::map<int,int> >();
    auto face_ref_sp = get_ref_face_space(face_id, face_to_element_dofs, *elem_map);
    auto map  = get_push_forward()->get_mapping();

    auto fmap = MappingSlice<FaceSpace::PushForwardType::dim, FaceSpace::PushForwardType::codim>::
    create(map, face_id, face_ref_sp->get_grid(), elem_map);
    auto fpf = FaceSpace::PushForwardType::create(fmap);
    auto face_space = FaceSpace::create(face_ref_sp,fpf);

    return face_space;
}

template <int dim_, int range_, int rank_>
void
NURBSSpace<dim_, range_, rank_>::
refine_h_weights(
    const std::array<bool,dim> &refinement_directions,
    const GridType &grid_old1)
{
    auto grid = this->get_grid();
    auto grid_old = this->get_grid()->get_grid_pre_refinement();

    auto knots_with_repetitions_pre_refinement = sp_space_->get_spline_space_previous_refinement()
                                                 ->compute_knots_with_repetition(
                                                     spline_space_t::BaseSpace::EndBehaviour::interpolatory,
                                                     sp_space_->get_end_behaviour());
    auto knots_with_repetitions = sp_space_->compute_knots_with_repetition(
                                      spline_space_t::BaseSpace::EndBehaviour::interpolatory,
                                      sp_space_->get_end_behaviour());

    for (int direction_id = 0; direction_id < dim; ++direction_id)
    {
        if (refinement_directions[direction_id])
        {
            // knots in the refined grid along the selected direction
            vector<Real> knots_new = grid->get_knot_coordinates(direction_id);

            // knots in the original (unrefined) grid along the selected direction
            vector<Real> knots_old = grid_old->get_knot_coordinates(direction_id);

            vector<Real> knots_added(knots_new.size());

            // find the knots in the refined grid that are not present in the old grid
            auto it = std::set_difference(
                          knots_new.begin(),knots_new.end(),
                          knots_old.begin(),knots_old.end(),
                          knots_added.begin());

            knots_added.resize(it-knots_added.begin());
            /*
                        Assert(false,ExcNotImplemented());
                        AssertThrow(false,ExcNotImplemented());
            //*/
            for (int comp_id = 0; comp_id < n_components; ++comp_id)
            {
                const int p = sp_space_->get_degree()(comp_id)[direction_id];
                const auto &U = knots_with_repetitions_pre_refinement(comp_id).get_data_direction(direction_id);
                const auto &X = knots_added;
                const auto &Ubar = knots_with_repetitions(comp_id).get_data_direction(direction_id);


                const int m = U.size()-1;
                const int r = X.size()-1;
                const int a = space_tools::find_span(p,X[0],U);
                const int b = space_tools::find_span(p,X[r],U)+1;

                const int n = m-p-1;

                const auto Pw = weights_(comp_id);
                const auto old_sizes = Pw.tensor_size();
                Assert(old_sizes(direction_id) == n+1,
                       ExcDimensionMismatch(old_sizes(direction_id), n+1));

                auto new_sizes = old_sizes;
                new_sizes(direction_id) += r+1; // r+1 new weights in the refinement direction
                Assert(new_sizes(direction_id) ==
                       sp_space_->get_num_basis(comp_id,direction_id),
                       ExcDimensionMismatch(new_sizes(direction_id),
                                            sp_space_->get_num_basis(comp_id,direction_id)));

                DynamicMultiArray<Real,dim> Qw(new_sizes);

                for (Index j = 0; j <= a-p; ++j)
                {
                    Qw.copy_slice(direction_id,j,Pw.get_slice(direction_id,j));
                }

                for (Index j = b-1; j <= n; ++j)
                {
                    Qw.copy_slice(direction_id,j+r+1,Pw.get_slice(direction_id,j));
                }

                Index i = b + p - 1;
                Index k = b + p + r;
                for (Index j = r; j >= 0; --j)
                {
                    while (X[j] <= U[i] && i > a)
                    {
                        Qw.copy_slice(direction_id,k-p-1,Pw.get_slice(direction_id,i-p-1));
                        k = k-1;
                        i = i-1;
                    }
                    Qw.copy_slice(direction_id,k-p-1,
                                  Qw.get_slice(direction_id,k-p));

                    for (Index l = 1; l <= p; ++l)
                    {
                        Index ind = k-p+l;

                        Real alfa = Ubar[k+l] - X[j];
                        if (fabs(alfa) == 0.0)
                        {
                            Qw.copy_slice(direction_id,ind-1,Qw.get_slice(direction_id,ind));
                        }
                        else
                        {
                            alfa = alfa / (Ubar[k+l] - U[i-p+l]);

                            Qw.copy_slice(direction_id,ind-1,
                                          alfa  * Qw.get_slice(direction_id,ind-1) +
                                          (1.0-alfa) * Qw.get_slice(direction_id,ind));
                        }
                    } // end loop l
                    k = k-1;
                } // end loop j

                weights_(comp_id) = Qw;

            } // end loop comp_id

        } // end if (refinement_directions[direction_id])

    } // end loop direction_id

    this->perform_post_construction_checks();
}

template <int dim_, int range_, int rank_>
void
NURBSSpace<dim_, range_, rank_>::
print_info(LogStream &out) const
{
    out << "NURBSSpace<" << dim_ << "," << range_ << ">" << endl;

    sp_space_->print_info(out);

//    for (auto w : weights_)
//        w.print_info(out);
    out.push("\t");
    for (int comp_id = 0; comp_id < n_components; comp_id++)
    {
        const auto weights_component = weights_(comp_id).get_data();
        out << "weights[" << comp_id << "] = { ";
        for (const Real &w : weights_component)
            out << w << " ";
        out << "}" << endl;
    }
    out.pop();
}


IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/nurbs_space.inst>


