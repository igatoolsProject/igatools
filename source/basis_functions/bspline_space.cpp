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
#include <igatools/geometry/mapping_slice.h>

using std::endl;
using std::array;
using std::vector;
using std::shared_ptr;
using std::make_shared;
using std::const_pointer_cast;

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
    operators_(knots,
               BaseSpace::compute_knots_with_repetition(this->get_end_behaviour()),
               BaseSpace::accumulated_interior_multiplicities(), deg)
{
    this->create_dofs_manager();


    // create a signal and a connection for the grid refinement
    this->connect_refinement_h_function(
        std::bind(&self_t::refine_h_after_grid_refinement, this,
                  std::placeholders::_1,std::placeholders::_2));
}



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
    operators_(knots,
               BaseSpace::compute_knots_with_repetition(this->get_end_behaviour()),
               BaseSpace::accumulated_interior_multiplicities(), deg)
{
    this->create_dofs_manager();

    // create a signal and a connection for the grid refinement
    this->connect_refinement_h_function(
        std::bind(&self_t::refine_h_after_grid_refinement, this,
                  std::placeholders::_1,std::placeholders::_2));
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



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
get_reference_space() const -> shared_ptr<const self_t>
{
    return this->shared_from_this();
}




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
get_ref_face_space(const Index face_id,
                   std::vector<Index> &face_to_element_dofs,
                   std::map<int, int> &elem_map) const
-> std::shared_ptr<RefFaceSpace>
{

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


//#if 0
template<int dim_, int range_, int rank_>
void
BSplineSpace<dim_, range_, rank_>::
refine_h_after_grid_refinement(
    const std::array<bool,dim> &refinement_directions,
    const GridType &grid_old)
{
    basis_indices_ = DofDistribution<dim, range, rank>(
                         this->get_grid(),
                         BaseSpace::accumulated_interior_multiplicities(),
                         BaseSpace::get_num_basis_table(),
                         BaseSpace::get_num_basis_per_element_table());

    operators_ = BernsteinExtraction<dim, range, rank>(
                     this->get_grid(),
                     BaseSpace::compute_knots_with_repetition(this->get_end_behaviour()),
                     BaseSpace::accumulated_interior_multiplicities(),
                     this->get_degree());


    dofs_manager_.reset();
    this->create_dofs_manager();
}





template<int dim_, int range_, int rank_>
void
BSplineSpace<dim_, range_, rank_>::
add_dofs_offset(const Index offset)
{
    basis_indices_.add_dofs_offset(offset);
}

template<int dim_, int range_, int rank_>
void
BSplineSpace<dim_, range_, rank_>::
create_dofs_manager()
{
    Assert(dofs_manager_ == nullptr,ExcInvalidState());
    dofs_manager_ = make_shared<DofsManager>(DofsManager());

    using DofsComponentContainer = std::vector<Index>;
    using DofsComponentView = ContainerView<DofsComponentContainer>;
    using DofsComponentConstView = ConstContainerView<DofsComponentContainer>;
    using DofsIterator = ConcatenatedIterator<DofsComponentView>;
    using DofsConstIterator = ConcatenatedConstIterator<DofsComponentView,DofsComponentConstView>;
    using SpaceDofsView = View<DofsIterator,DofsConstIterator>;

    //---------------------------------------------------------------------------------------------
    dofs_manager_->dofs_view_open();
    auto &index_space = this->get_basis_indices().get_index_distribution();

    vector<DofsComponentView> space_components_view;
    for (auto &index_space_comp : index_space)
    {
        vector<Index> &index_space_comp_data = const_cast<vector<Index> &>(index_space_comp.get_data());
        DofsComponentView index_space_comp_view(
            index_space_comp_data.begin(),index_space_comp_data.end());

        space_components_view.push_back(index_space_comp_view);
    }

    DofsIterator space_dofs_begin(space_components_view,0);
    DofsIterator space_dofs_end(space_components_view,IteratorState::pass_the_end);
    SpaceDofsView dofs_space_view(space_dofs_begin,space_dofs_end);

    const Index space_id = 0;
    dofs_manager_->add_dofs_space_view(space_id,this->get_num_basis(),dofs_space_view);
    dofs_manager_->dofs_view_close(false);
    //---------------------------------------------------------------------------------------------


    //---------------------------------------------------------------------------------------------
    // getting the views of the dofs on each element of the space
    dofs_manager_->elements_dofs_view_open();

    const auto &elements_view = this->get_basis_indices().get_elements_view();
    for (const auto &elem_view : elements_view)
        dofs_manager_->add_element_dofs_view(elem_view);

    dofs_manager_->elements_dofs_view_close();
    //---------------------------------------------------------------------------------------------
}


template<int dim_, int range_, int rank_>
std::shared_ptr<DofsManager>
BSplineSpace<dim_, range_, rank_>::
get_dofs_manager() const
{
    return dofs_manager_;
}

template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
get_basis_indices() const -> const DofDistribution<dim, range, rank> &
{
    return basis_indices_;
}

template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
get_basis_indices() -> DofDistribution<dim, range, rank> &
{
    return basis_indices_;
}


template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
basis_flat_to_tensor(const Index index, const Index comp) const -> TensorIndex<dim>
{
    return basis_indices_.basis_flat_to_tensor(index,comp);
}


template<int dim_, int range_, int rank_>
Index
BSplineSpace<dim_, range_, rank_>::
basis_tensor_to_flat(const TensorIndex<dim> &tensor_index,
                     const Index comp) const
{
    return basis_indices_.basis_tensor_to_flat(tensor_index, comp);
}

template<int dim_, int range_, int rank_>
std::vector<Index>
BSplineSpace<dim_, range_, rank_>::
get_loc_to_global(const TensorIndex<dim> &j) const
{
    return basis_indices_.get_loc_to_global_indices(j);
}

template<int dim_, int range_, int rank_>
void
BSplineSpace<dim_, range_, rank_>::
print_info(LogStream &out) const
{
    out.begin_item("Spline Space:");
    BaseSpace::print_info(out);
    out.end_item();

    out.begin_item("Basis Indices:");
    basis_indices_.print_info(out);
    out.end_item();

    out.begin_item("Bernstein Extraction:");
    operators_.print_info(out);
    out.end_item();
}

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/bspline_space.inst>
