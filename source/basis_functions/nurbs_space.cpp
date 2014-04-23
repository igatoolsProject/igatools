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

#include <igatools/base/exceptions.h>

using std::array;
using std::vector;
using std::endl;
using std::shared_ptr;

IGA_NAMESPACE_OPEN

//TODO(pauletti, Jan 19, 2014): remove the is_bspline from this class

template <int dim_, int range_, int rank_>
void
NURBSSpace<dim_, range_, rank_>::
create_refinement_connection()
{
    //----------------------------------
    // create a signal and a connection for the grid refinement
    this->connect_refinement_h_function(
        std::bind(&self_t::refine_h_weights,
                  this,
                  std::placeholders::_1,std::placeholders::_2));
    //----------------------------------
}



template <int dim_, int range_, int rank_>
NURBSSpace<dim_, range_, rank_>::
NURBSSpace(shared_ptr< GridType > knots, const int &degree)
    :
    base_t(knots, degree)
{
    //----------------------------------------------------------------------------------------------
    // initialize all the weights to 1.0 (then this space will have the same mathematical structure
    // of a BSpline space)
    const auto n_dofs = this->get_num_dofs();
    for (int comp_id = 0; comp_id < n_components; ++comp_id)
    {
        const auto n_dofs_component = n_dofs(comp_id);

        weights_(comp_id).resize(n_dofs_component);
        weights_(comp_id).fill(1.0);
    }
    //----------------------------------------------------------------------------------------------
    create_refinement_connection();

    perform_post_construction_checks();
}



template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
create(shared_ptr< GridType > knots,
       const int &degree) -> shared_ptr< self_t >
{
    return shared_ptr< self_t >(new self_t(knots, degree));
}



template <int dim_, int range_, int rank_>
NURBSSpace<dim_, range_, rank_>::
NURBSSpace(
    shared_ptr<GridType> knots,
    const DegreeTable &degree)
    :
    base_t(knots, degree)
{

    const auto n_dofs = this->get_num_dofs();
    for (int comp_id = 0; comp_id < n_components; ++comp_id)
    {
        const auto n_dofs_component = n_dofs(comp_id);

        weights_(comp_id).resize(n_dofs_component);
        weights_(comp_id).fill(1.0);
    }
    create_refinement_connection();

    perform_post_construction_checks();
}



template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
create(shared_ptr<GridType> knots,
       const DegreeTable &degree)
-> shared_ptr< self_t >
{
    return shared_ptr< self_t >(new self_t(knots, degree));
}



template <int dim_, int range_, int rank_>
NURBSSpace<dim_, range_, rank_>::
NURBSSpace(
    shared_ptr< GridType > knots,
    const MultiplicityTable &mult_vector,
    const DegreeTable &degree)
    :
    base_t(knots, mult_vector, degree)
{
    //----------------------------------------------------------------------------------------------
    // initialize all the weights to 1.0 (then this space will have the same mathematical structure
    // of a BSpline space)
    const auto n_dofs = this->get_num_dofs();
    for (int comp_id = 0; comp_id < n_components; ++comp_id)
    {
        const auto n_dofs_component = n_dofs(comp_id);

        weights_(comp_id).resize(n_dofs_component);
        weights_(comp_id).fill(1.0);
    }
    //----------------------------------------------------------------------------------------------

    create_refinement_connection();

    perform_post_construction_checks();
}



template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
create(shared_ptr< GridType > knots,
       const MultiplicityTable &mult_vector,
       const DegreeTable &degree) -> shared_ptr< self_t >
{
    return (shared_ptr< self_t >(new self_t(knots, mult_vector, degree)));
}



template <int dim_, int range_, int rank_>
NURBSSpace<dim_, range_, rank_>::
NURBSSpace(
    std::shared_ptr< GridType > knots,
    const MultiplicityTable &mult_vector,
    const DegreeTable &degree,
    const WeightsTable &weights)
    :
    base_t(knots, mult_vector, degree),
    weights_(weights)
{

    create_refinement_connection();
    perform_post_construction_checks();
}



template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
create(std::shared_ptr< GridType > knots,
       const MultiplicityTable &mult_vector,
       const DegreeTable &degree,
       const WeightsTable &weights) -> shared_ptr< self_t >
{
    return shared_ptr< self_t >(new self_t(knots, mult_vector, degree, weights));
}



template <int dim_, int range_, int rank_>
void
NURBSSpace<dim_, range_, rank_>::
perform_post_construction_checks() const
{
    //-------------------------------------------
    // check that the number of weights is equal to the number of basis functions in the space
    using weights_container_t = StaticMultiArray<DynamicMultiArray<Real,dim>,range,rank>;
    Size n_weights = 0;
    for (int comp_id = 0; comp_id < weights_container_t::n_entries; ++comp_id)
        n_weights += weights_(comp_id).flat_size();

    Assert(this->get_num_basis() == n_weights,
           ExcDimensionMismatch(this->get_num_basis(),n_weights));
    //------------------------------------------
}



template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
begin() const -> ElementIterator
{
   // return ElementIterator(std::enable_shared_from_this<NURBSSpace<dim_,range_,rank_>>::shared_from_this(), 0);
}



template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
last() const -> ElementIterator
{
   // return ElementIterator(std::enable_shared_from_this<NURBSSpace<dim_,range_,rank_>>::shared_from_this(),
   //                        this->get_grid()->get_num_elements() - 1);
}



template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
end() const -> ElementIterator
{
   // return ElementIterator(std::enable_shared_from_this<NURBSSpace<dim_,range_,rank_>>::shared_from_this(),
                      //     IteratorState::pass_the_end);
}



template <int dim_, int range_, int rank_>
auto
NURBSSpace<dim_, range_, rank_>::
get_weights() const -> const StaticMultiArray<DynamicMultiArray<Real,dim>, range, rank >
{
    return (weights_);
}



template <int dim_, int range_, int rank_>
void
NURBSSpace<dim_, range_, rank_>::
reset_weights(const StaticMultiArray<DynamicMultiArray<Real,dim>,range,rank> &weights)
{
    //-------------------------------------------------------------------------
    for (int i = 0; i < StaticMultiArray<DynamicMultiArray<Real,dim>,range,rank>::n_entries; ++i)
    {
        Assert(this->get_component_num_basis(i) == int(weights(i).flat_size()),
               ExcDimensionMismatch(this->get_component_num_basis(i), weights(i).flat_size()));
    }
    //-------------------------------------------------------------------------
    weights_ = weights;
}



template <int dim_, int range_, int rank_>
void
NURBSSpace<dim_, range_, rank_>::
refine_h_weights(
    const std::array<bool,dim> &refinement_directions,
    const GridType &grid_old)
{
    for (int direction_id = 0; direction_id < dim; ++direction_id)
    {
        if (refinement_directions[direction_id])
        {
            // knots in the refined grid along the selected direction
            vector<Real> knots_new = this->get_grid()->get_knot_coordinates(direction_id);

            // knots in the original (unrefined) grid along the selected direction
            vector<Real> knots_old = grid_old.get_knot_coordinates(direction_id);

            vector<Real> knots_added(knots_new.size());

            // find the knots in the refined grid that are not present in the old grid
            auto it = std::set_difference(
                          knots_new.begin(),knots_new.end(),
                          knots_old.begin(),knots_old.end(),
                          knots_added.begin());

            knots_added.resize(it-knots_added.begin());


            for (int comp_id = 0; comp_id < n_components; ++comp_id)
            {
                const int p = this->get_degree()(comp_id)[direction_id];

                const auto &U =
                    this->knots_with_repetitions_pre_refinement_(comp_id).get_data_direction(direction_id);
                const auto &X = knots_added;
                const auto &Ubar = this->knots_with_repetitions_(comp_id).get_data_direction(direction_id);

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
                       this->get_component_dir_num_basis(comp_id,direction_id),
                       ExcDimensionMismatch(new_sizes(direction_id),
                                            this->get_component_dir_num_basis(comp_id,direction_id)));

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

    base_t::print_info(out);

    out.push("\t");
    for (int comp_id = 0; comp_id < n_components; comp_id++)
    {
        const auto weights_component = weights_(comp_id).get_data();
        out << "weights[" << comp_id << "] = { ";
        for (const Real & w : weights_component)
            out << w << " ";
        out << "}" << endl;
    }
    out.pop();
}


IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/nurbs_space.inst>


