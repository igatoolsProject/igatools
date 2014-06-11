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

#include <igatools/basis_functions/bspline_space.h>
#include <igatools/geometry/push_forward.h>
#include <igatools/geometry/identity_mapping.h>

using std::endl;
using std::array;
using std::vector;
using std::shared_ptr;
using std::make_shared;

IGA_NAMESPACE_OPEN

template<int dim_, int range_, int rank_>
BSplineSpace<dim_, range_, rank_>::
BSplineSpace(const int degree, shared_ptr<GridType> grid)
    :
    BSplineSpace(TensorIndex<dim>(degree), grid)
{}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
create(const int degree, shared_ptr< GridType > knots) -> shared_ptr<self_t>
{
    return shared_ptr<self_t>(new self_t(degree, knots));
}



template<int dim_, int range_, int rank_>
BSplineSpace<dim_, range_, rank_>::
BSplineSpace(const TensorIndex<dim> &degree, shared_ptr<GridType> knots)
    :
    BSplineSpace(DegreeTable(degree), knots, true)
{}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
create(const TensorIndex<dim> &degree, shared_ptr<GridType> knots)
-> shared_ptr<self_t>
{
    return shared_ptr<self_t>(new self_t(degree, knots));
}



template<int dim_, int range_, int rank_>
BSplineSpace<dim_, range_, rank_>::
BSplineSpace(const DegreeTable &deg,
             std::shared_ptr<GridType> knots,
             const bool homogeneous_range)
             :
             BaseSpace(deg, knots, BaseSpace::InteriorReg::maximum),
             basis_indices_(knots,BaseSpace::accumulated_interior_multiplicities(),
                            BaseSpace::get_num_basis_table(),BaseSpace::get_num_basis_per_element_table()),
                            operators_(knots, BaseSpace::compute_knots_with_repetition(BaseSpace::EndBehaviour::interpolatory),
                            BaseSpace::accumulated_interior_multiplicities(), deg)

{}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
create(const DegreeTable &deg,
       std::shared_ptr<GridType> knots,
       const bool homogeneous_range) -> shared_ptr<self_t>
{
    return shared_ptr<self_t>(new self_t(deg, knots, homogeneous_range));
}



template<int dim_, int range_, int rank_>
BSplineSpace<dim_, range_, rank_>::
BSplineSpace(const DegreeTable &deg,
             std::shared_ptr<GridType> knots,
             std::shared_ptr<const MultiplicityTable> interior_mult,
             const EndBehaviourTable &ends)
    :
    BaseSpace(deg, knots, interior_mult),
    basis_indices_(knots,BaseSpace::accumulated_interior_multiplicities(),
                   BaseSpace::get_num_basis_table(),BaseSpace::get_num_basis_per_element_table()),
    operators_(knots, BaseSpace::compute_knots_with_repetition(BaseSpace::EndBehaviour::interpolatory),
               BaseSpace::accumulated_interior_multiplicities(), deg)
{
}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
create(const DegreeTable &deg,
       std::shared_ptr<GridType> knots,
       std::shared_ptr<const MultiplicityTable> interior_mult,
       const EndBehaviourTable &ends)
       -> shared_ptr<self_t>
{
    return shared_ptr<self_t>(new self_t(deg, knots, interior_mult, ends));
}




//template<int dim_, int range_, int rank_>
//void
//BSplineSpace<dim_, range_, rank_>::
//init()
//{
//    init_dofs();
//
//    // create a signal and a connection for the grid refinement
//    this->connect_refinement_h_function(
//        std::bind(&self_t::refine_h_after_grid_refinement, this,
//                  std::placeholders::_1,std::placeholders::_2));
//}




template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
get_reference_space() const -> shared_ptr<const self_t>
{
    return this->shared_from_this();
}






//template<int dim_, int range_, int rank_>
//bool
//BSplineSpace<dim_, range_, rank_>::
//is_range_homogeneous() const
//{
//    return homogeneous_range_;
//}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::begin() const -> ElementIterator
{
    return ElementIterator(this->shared_from_this(), 0);
}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::last() const -> ElementIterator
{
    return ElementIterator(this->shared_from_this(),
                           this->get_grid()->get_num_elements() - 1);
}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::end() const -> ElementIterator
{
    return ElementIterator(this->shared_from_this(),
                           IteratorState::pass_the_end);
}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
get_face_space(const Index face_id,
		std::vector<Index> &face_to_element_dofs) const
		-> std::shared_ptr<RefFaceSpace>
{
	std::map<int, int> elem_map;
	auto face_grid   = this->get_grid()->get_face_grid(face_id, elem_map);
	auto face_mult   = this->get_face_mult(face_id);
	auto face_degree = this->get_face_degree(face_id);

	// TODO (pauletti, Jun 4, 2014): make sure the face space is compatible with the space end behaviou
	auto f_space = RefFaceSpace::create(face_degree, face_grid, face_mult);

	const auto &active_dirs = UnitElement<dim>::face_active_directions[face_id];
	const auto const_dir = UnitElement<dim>::face_constant_direction[face_id];
	const auto face_side = UnitElement<dim>::face_side[face_id];

	TensorIndex<dim> tensor_index;
	face_to_element_dofs.resize(f_space->get_num_basis());
	int k=0;
	int offset=0;
	for (auto comp : components)
	{
		const int face_n_basis = f_space->get_num_basis(comp);
		for (Index i = 0; i < face_n_basis; ++i, ++k)
		{
			const auto f_tensor_idx = f_space->basis_flat_to_tensor(i,comp);
			const int fixed_idx =
					face_side * (this->get_num_basis(comp,const_dir) - 1);
			for (int j : RefFaceSpace::dims)
				tensor_index[active_dirs[j]] =  f_tensor_idx[j];
			tensor_index[const_dir] = fixed_idx;

			const Index dof = basis_tensor_to_flat(tensor_index, comp);

			face_to_element_dofs[k] = offset + dof;
		}
		offset += this->get_num_basis(comp);
	}
	return f_space;
}





template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
get_push_forward() -> shared_ptr<PushForwardType>
{

    return
    PushForwardType::create(IdentityMapping<dim>::create(this->get_grid()));
}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
get_push_forward() const -> shared_ptr<const PushForwardType>
{
    using PushForwardType1 = PushForward<Transformation::h_grad,dim,0>;
    auto grid = this->get_grid();
    auto push_fwd =
    PushForwardType1::create(
        IdentityMapping<dim>::create(
            make_shared<GridType>(GridType(*grid))));

    return push_fwd;
}


#if 0
template<int dim_, int range_, int rank_>
void
BSplineSpace<dim_, range_, rank_>::
refine_h_after_grid_refinement(
    const std::array<bool,dim> &refinement_directions,
    const GridType &grid_old)
{
    // keeping the original knots (with repetitions) before the h-refinement
    knots_with_repetitions_pre_refinement_ = knots_with_repetitions_;

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



            for (int comp_id = 0; comp_id < self_t::n_components; ++comp_id)
            {
                //--------------------------------------------------------
                // creating the new multiplicity
                const vector<int> &mult_old = mult_(comp_id).get_data_direction(direction_id);
                const int n_mult_old = mult_old.size();

                const int n_mult_to_add = n_mult_old - 1;
                const int n_mult_new = n_mult_old + n_mult_to_add;

                vector<int> mult_new(n_mult_new);
                for (int i = 0; i < n_mult_to_add; ++i)
                {
                    mult_new[2*i  ] = mult_old[i];
                    mult_new[2*i+1] = 1;
                }
                mult_new[n_mult_new-1] = mult_old[n_mult_old-1];

                mult_(comp_id).copy_data_direction(direction_id,mult_new);
                //--------------------------------------------------------


                //--------------------------------------------------------
                const auto &knots_old =
                    knots_with_repetitions_pre_refinement_(comp_id).get_data_direction(direction_id);

                // knots with repetitions after refinement
                vector<Real> Ubar = knots_old;
                Ubar.insert(Ubar.end(),
                            knots_added.begin(),knots_added.end());

                sort(Ubar.begin(),Ubar.end());

                knots_with_repetitions_(comp_id).copy_data_direction(direction_id,Ubar);
                //--------------------------------------------------------
            } // end loop comp_id

        } // end if(refinement_directions[direction_id])

    } // end loop direction_id

    init_dofs();
}
#endif




template<int dim_, int range_, int rank_>
void
BSplineSpace<dim_, range_, rank_>::
print_info(LogStream &out) const
{
    BaseSpace::print_info(out);

    basis_indices_.print_info(out);

    operators_.print_info(out);

#if 0

    //TODO: Do we need to call external functions from this output operator?
    out << "Dofs: " << dof_tools::get_dofs(this->shared_from_this())  << endl;

    const SparsityPattern &sparsity_pattern =
        dof_tools::get_sparsity_pattern(this->shared_from_this());
    out << "Num overlapping funcs: ";
    out << sparsity_pattern.get_num_overlapping_funcs() << endl;
#endif
}

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/bspline_space.inst>

