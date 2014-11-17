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

#include <igatools/basis_functions/new_physical_space.h>
#include <igatools/geometry/mapping_slice.h>
#include <igatools/basis_functions/space_manager.h>


using std::array;
using std::shared_ptr;
using std::make_shared;
using std::const_pointer_cast;
using std::endl;

IGA_NAMESPACE_OPEN


template <class RefSpace_,int codim_, Transformation type_>
const std::array<int, NewPhysicalSpace<RefSpace_, codim_, type_>::n_components>
NewPhysicalSpace<RefSpace_, codim_, type_>::components = sequence<NewPhysicalSpace<RefSpace_, codim_, type_>::n_components>();


template <class RefSpace_,int codim_, Transformation type_>
NewPhysicalSpace<RefSpace_, codim_, type_>::
NewPhysicalSpace(shared_ptr<RefSpace> ref_space,
                 shared_ptr<MapFunc> map_func)
    :
    BaseSpace(ref_space->get_grid()),
    ref_space_(ref_space),
    map_func_(map_func)
{
//TODO(pauletti, Jan 18, 2014): put static assert on h_div, h_curl range and rank
    Assert(ref_space_ != nullptr, ExcNullPtr());
    Assert(map_func_ != nullptr, ExcNullPtr());

//    Assert(ref_space_->get_grid() == push_forward_->get_mapping()->get_grid(),
//           ExcMessage("Reference space and mapping grids are not the same."))
}



template <class RefSpace_,int codim_, Transformation type_>
auto
NewPhysicalSpace<RefSpace_, codim_, type_>::
create(shared_ptr<RefSpace> ref_space,
       shared_ptr<MapFunc> map_func) -> shared_ptr<self_t>
{
    return shared_ptr<self_t>(new self_t(ref_space, map_func));
}



template <class RefSpace_,int codim_, Transformation type_>
auto
NewPhysicalSpace<RefSpace_, codim_, type_>::
begin() const -> ElementIterator
{
    return ElementIterator(this->shared_from_this(), 0);
}



template <class RefSpace_,int codim_, Transformation type_>
auto
NewPhysicalSpace<RefSpace_, codim_, type_>::
last() const -> ElementIterator
{
    return ElementIterator(this->shared_from_this(),
                           ref_space_->get_grid()->get_num_active_elems() - 1);
}



template <class RefSpace_,int codim_, Transformation type_>
auto
NewPhysicalSpace<RefSpace_, codim_, type_>::
end() const -> ElementIterator
{
    return ElementIterator(this->shared_from_this(),
                           IteratorState::pass_the_end);
}

#if 0
template <class RefSpace_,int codim_, Transformation type_>
auto
NewPhysicalSpace<RefSpace_, codim_, type_>::
get_element(const Index elem_flat_id) const -> ElementAccessor
{
    Assert(elem_flat_id >= 0 && elem_flat_id < ref_space_->get_grid()->get_num_active_elems(),
           ExcIndexRange(elem_flat_id,0,ref_space_->get_grid()->get_num_active_elems()));

    auto elem = this->begin();
    for (int i = 0 ; i < elem_flat_id ; ++i)
        ++elem;

    return *elem;
}

#endif

template <class RefSpace_,int codim_, Transformation type_>
Index
NewPhysicalSpace<RefSpace_, codim_, type_>::
get_num_basis() const
{
    return ref_space_->get_num_basis();
}


#if 0
template <class RefSpace_,int codim_, Transformation type_>
int
NewPhysicalSpace<RefSpace_, codim_, type_>::
get_num_basis_per_element() const
{
    return ref_space_->get_num_basis_per_element();
}



template <class RefSpace_,int codim_, Transformation type_>
auto
NewPhysicalSpace<RefSpace_, codim_, type_>::
get_push_forward() const -> shared_ptr<const PushForwardType>
{
    return shared_ptr<const PushForwardType>(push_forward_);
}


#endif
template <class RefSpace_,int codim_, Transformation type_>
auto
NewPhysicalSpace<RefSpace_, codim_, type_>::
get_reference_space() const -> shared_ptr<const RefSpace>
{
    return shared_ptr<const RefSpace>(ref_space_);
}



template <class RefSpace_,int codim_, Transformation type_>
auto
NewPhysicalSpace<RefSpace_, codim_, type_>::
get_reference_space() -> shared_ptr<RefSpace>
{
    return ref_space_;
}



template <class RefSpace_,int codim_, Transformation type_>
template<int k>
auto
NewPhysicalSpace<RefSpace_, codim_, type_>::
get_sub_space(const int s_id, InterSpaceMap<k> &dof_map,
              std::shared_ptr<CartesianGrid<k>> sub_grid,
              std::shared_ptr<InterGridMap<k>> elem_map) const
-> std::shared_ptr<SubSpace<k> >
{
    using SubMap = SubMapFunction<k, dim, space_dim>;
    auto grid =  this->get_grid();

    auto sub_ref_space = ref_space_->get_ref_sub_space(s_id, dof_map, sub_grid);
    auto sub_map_func = SubMap::create(sub_grid,  map_func_, s_id, *elem_map);
    auto sub_space = SubSpace<k>::create(sub_ref_space, sub_map_func);
    return sub_space;
}



#if 0
template <class RefSpace_,int codim_, Transformation type_>
auto
NewPhysicalSpace<RefSpace_, codim_, type_>::
get_face_space(const Index face_id,
               vector<Index> &face_to_element_dofs) const -> shared_ptr<FaceSpace>
{
    auto elem_map = std::make_shared<typename GridType::FaceGridMap >();
    auto face_ref_sp = ref_space_->get_ref_face_space(face_id, face_to_element_dofs, *elem_map);
    auto map  = push_forward_->get_mapping();

    auto fmap = MappingSlice<FaceSpace::PushForwardType::dim, FaceSpace::PushForwardType::codim>::
    create(map, face_id, face_ref_sp->get_grid(), elem_map);
    auto fpf = FaceSpace::PushForwardType::create(fmap);
    auto face_space = FaceSpace::create(face_ref_sp,fpf);

    return face_space;
}


template <class RefSpace_,int codim_, Transformation type_>
Index
NewPhysicalSpace<RefSpace_, codim_, type_>::
get_id() const
{
    return ref_space_->get_id();
}
#endif

template <class RefSpace_,int codim_, Transformation type_>
vector<Index>
NewPhysicalSpace<RefSpace_, codim_, type_>::
get_loc_to_global(const CartesianGridElement<dim> &element) const
{
    return ref_space_->get_loc_to_global(element);
}


template <class RefSpace_,int codim_, Transformation type_>
vector<Index>
NewPhysicalSpace<RefSpace_, codim_, type_>::
get_loc_to_patch(const CartesianGridElement<dim> &element) const
{
    return ref_space_->get_loc_to_patch(element);
}



template <class RefSpace_,int codim_, Transformation type_>
auto
NewPhysicalSpace<RefSpace_, codim_, type_>::
get_space_manager() -> shared_ptr<SpaceManager>
{
    auto space_manager = make_shared<SpaceManager>(SpaceManager());

    auto this_space = this->shared_from_this();

    space_manager->spaces_insertion_open();
    space_manager->add_space(this_space);
    space_manager->spaces_insertion_close();


    space_manager->spaces_connectivity_open();
    space_manager->add_spaces_connection(this_space);
    space_manager->spaces_connectivity_close();

    return space_manager;
}

template <class RefSpace_,int codim_, Transformation type_>
auto
NewPhysicalSpace<RefSpace_, codim_, type_>::
get_space_manager() const -> std::shared_ptr<const SpaceManager>
{
    return const_cast<self_t &>(*this).get_space_manager();
}



template <class RefSpace_,int codim_, Transformation type_>
auto
NewPhysicalSpace<RefSpace_, codim_, type_>::
get_dof_distribution_global() const -> const DofDistribution<dim, range, rank> &
{
    return ref_space_->get_dof_distribution_global();
}



template <class RefSpace_,int codim_, Transformation type_>
auto
NewPhysicalSpace<RefSpace_, codim_, type_>::
get_dof_distribution_global() -> DofDistribution<dim, range, rank> &
{
    return ref_space_->get_dof_distribution_global();
}



template <class RefSpace_,int codim_, Transformation type_>
auto
NewPhysicalSpace<RefSpace_, codim_, type_>::
get_dof_distribution_patch() const -> const DofDistribution<dim, range, rank> &
{
    return ref_space_->get_dof_distribution_patch();
}



template <class RefSpace_,int codim_, Transformation type_>
auto
NewPhysicalSpace<RefSpace_, codim_, type_>::
get_dof_distribution_patch() -> DofDistribution<dim, range, rank> &
{
    return ref_space_->get_dof_distribution_patch();
}

#if 0
template <class RefSpace_,int codim_, Transformation type_>
auto
NewPhysicalSpace<RefSpace_, codim_, type_>::
get_degree() const -> const DegreeTable &
{
    return ref_space_->get_degree();
}

#endif

template <class RefSpace_,int codim_, Transformation type_>
void
NewPhysicalSpace<RefSpace_, codim_, type_>::
print_info(LogStream &out) const
{
    out.begin_item("Reference space:");
    ref_space_->print_info(out);
    out.end_item();
}


IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/new_physical_space.inst>


