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

#include <igatools/basis_functions/physical_space.h>
#include <igatools/functions/function.h>
#include <igatools/functions/sub_function.h>
#include <igatools/basis_functions/phys_space_element_handler.h>


using std::shared_ptr;
using std::make_shared;
using std::const_pointer_cast;
using std::endl;

IGA_NAMESPACE_OPEN


template <int dim_, int range_, int rank_, int codim_, Transformation type_>
const SafeSTLArray<int, PhysicalSpace<dim_, range_, rank_, codim_, type_>::n_components>
PhysicalSpace<dim_, range_, rank_, codim_, type_>::components =
    sequence<PhysicalSpace<dim_, range_, rank_, codim_, type_>::n_components>();


template <int dim_, int range_, int rank_, int codim_, Transformation type_>
PhysicalSpace<dim_, range_, rank_, codim_, type_>::
PhysicalSpace(shared_ptr<RefSpace> ref_space,
              shared_ptr<MapFunc> map_func)
    :
    base_t(ref_space->get_grid()),
    ref_space_(ref_space),
    map_func_(map_func->clone())
{
//TODO(pauletti, Jan 18, 2014): put static assert on h_div, h_curl range and rank
    Assert(ref_space_ != nullptr, ExcNullPtr());
    Assert(map_func_ != nullptr, ExcNullPtr());

    Assert(ref_space_->get_grid() == map_func_->get_grid(),
           ExcMessage("Reference space and mapping grids are not the same."))
}



template <int dim_, int range_, int rank_, int codim_, Transformation type_>
auto
PhysicalSpace<dim_, range_, rank_, codim_, type_>::
create(shared_ptr<RefSpace> ref_space,
       shared_ptr<MapFunc> map_func) -> shared_ptr<self_t>
{
    auto sp = shared_ptr<self_t>(new self_t(ref_space, map_func));

    sp->create_connection_for_insert_knots(sp);

    return sp;
}


template <int dim_, int range_, int rank_, int codim_, Transformation type_>
auto
PhysicalSpace<dim_, range_, rank_, codim_, type_>::
get_this_space() const -> std::shared_ptr<const self_t >
{
    auto sp = const_cast<self_t *>(this)->shared_from_this();
    auto this_space = std::dynamic_pointer_cast<self_t>(sp);
    Assert(this_space != nullptr,ExcNullPtr());

    return this_space;
}


template <int dim_, int range_, int rank_, int codim_, Transformation type_>
auto
PhysicalSpace<dim_, range_, rank_, codim_, type_>::
create_element(const Index flat_index) const
-> std::shared_ptr<SpaceElement<dim_,codim_,range_,rank_>>
{
    auto elem = make_shared<ElementAccessor>(this->get_this_space(),flat_index);
    Assert(elem != nullptr, ExcNullPtr());

    return elem;
}



#if 0
template <int dim_, int range_, int rank_, int codim_, Transformation type_>
auto
PhysicalSpace<dim_, range_, rank_, codim_, type_>::
get_element(const Index elem_flat_id) const -> ElementAccessor
{
    Assert(elem_flat_id >= 0 && elem_flat_id < ref_space_->get_grid()->get_num_all_elems(),
           ExcIndexRange(elem_flat_id,0,ref_space_->get_grid()->get_num_all_elems()));

    auto elem = this->begin();
    for (int i = 0 ; i < elem_flat_id ; ++i)
        ++elem;

    return *elem;
}

#endif

template <int dim_, int range_, int rank_, int codim_, Transformation type_>
Index
PhysicalSpace<dim_, range_, rank_, codim_, type_>::
get_num_basis() const
{
    return ref_space_->get_num_basis();
}



#if 0
template <int dim_, int range_, int rank_, int codim_, Transformation type_>
int
PhysicalSpace<dim_, range_, rank_, codim_, type_>::
get_num_basis_per_element() const
{
    return ref_space_->get_num_basis_per_element();
}



template <int dim_, int range_, int rank_, int codim_, Transformation type_>
auto
PhysicalSpace<dim_, range_, rank_, codim_, type_>::
get_push_forward() const -> shared_ptr<const PushForwardType>
{
    return shared_ptr<const PushForwardType>(push_forward_);
}


#endif
template <int dim_, int range_, int rank_, int codim_, Transformation type_>
auto
PhysicalSpace<dim_, range_, rank_, codim_, type_>::
get_reference_space() const -> shared_ptr<const RefSpace>
{
    return shared_ptr<const RefSpace>(ref_space_);
}



template <int dim_, int range_, int rank_, int codim_, Transformation type_>
auto
PhysicalSpace<dim_, range_, rank_, codim_, type_>::
get_reference_space() -> shared_ptr<RefSpace>
{
    return ref_space_;
}



template <int dim_, int range_, int rank_, int codim_, Transformation type_>
template<int k>
auto
PhysicalSpace<dim_, range_, rank_, codim_, type_>::
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
template <int dim_, int range_, int rank_, int codim_, Transformation type_>
auto
PhysicalSpace<dim_, range_, rank_, codim_, type_>::
get_face_space(const Index face_id,
               SafeSTLVector<Index> &face_to_element_dofs) const -> shared_ptr<FaceSpace>
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


template <int dim_, int range_, int rank_, int codim_, Transformation type_>
Index
PhysicalSpace<dim_, range_, rank_, codim_, type_>::
get_id() const
{
    return ref_space_->get_id();
}
#endif


template <int dim_, int range_, int rank_, int codim_, Transformation type_>
void
PhysicalSpace<dim_, range_, rank_, codim_, type_>::
get_element_dofs(
    const Index element_id,
    SafeSTLVector<Index> &dofs_global,
    SafeSTLVector<Index> &dofs_local_to_patch,
    SafeSTLVector<Index> &dofs_local_to_elem,
    const std::string &dofs_property) const
{
    return ref_space_->get_element_dofs(
               element_id,
               dofs_global,
               dofs_local_to_patch,
               dofs_local_to_elem,dofs_property);
}




template <int dim_, int range_, int rank_, int codim_, Transformation type_>
auto
PhysicalSpace<dim_, range_, rank_, codim_, type_>::
get_dof_distribution() const -> std::shared_ptr<const DofDistribution<dim, range, rank> >
{
    return ref_space_->get_dof_distribution();
}



template <int dim_, int range_, int rank_, int codim_, Transformation type_>
auto
PhysicalSpace<dim_, range_, rank_, codim_, type_>::
get_dof_distribution() -> std::shared_ptr<DofDistribution<dim, range, rank> >
{
    return ref_space_->get_dof_distribution();
}



#if 0
template <int dim_, int range_, int rank_, int codim_, Transformation type_>
auto
PhysicalSpace<dim_, range_, rank_, codim_, type_>::
get_degree() const -> const DegreeTable &
{
    return ref_space_->get_degree();
}

#endif

template <int dim_, int range_, int rank_, int codim_, Transformation type_>
void
PhysicalSpace<dim_, range_, rank_, codim_, type_>::
print_info(LogStream &out) const
{
    out.begin_item("Reference space:");
    ref_space_->print_info(out);
    out.end_item();
}



template <int dim_, int range_, int rank_, int codim_, Transformation type_>
auto
PhysicalSpace<dim_, range_, rank_, codim_, type_>::
get_elem_handler() const -> std::shared_ptr<SpaceElementHandler<dim_,codim_,range_,rank_>>
{
    auto sp = const_cast<self_t *>(this)->shared_from_this();
    auto this_space = std::dynamic_pointer_cast<self_t>(sp);
    Assert(this_space != nullptr,ExcNullPtr());

    return ElementHandler::create(this_space);
}



template <int dim_, int range_, int rank_, int codim_, Transformation type_>
auto
PhysicalSpace<dim_, range_, rank_, codim_, type_>::
get_map_func() const -> std::shared_ptr<MapFunc>
{
    return map_func_;
}


template <int dim_, int range_, int rank_, int codim_, Transformation type_>
std::set<Index>
PhysicalSpace<dim_, range_, rank_, codim_, type_>::
get_interior_dofs() const
{
    return ref_space_->get_interior_dofs();
}


template <int dim_, int range_, int rank_, int codim_, Transformation type_>
void
PhysicalSpace<dim_, range_, rank_, codim_, type_>::
create_connection_for_insert_knots(std::shared_ptr<self_t> space)
{
    Assert(space != nullptr, ExcNullPtr());
    Assert(&(*space) == &(*this), ExcMessage("Different objects."));

    auto func_to_connect =
        std::bind(&self_t::rebuild_after_insert_knots,
                  space.get(),
                  std::placeholders::_1,
                  std::placeholders::_2);

    using SlotType = typename CartesianGrid<dim>::SignalInsertKnotsSlot;
    this->connect_insert_knots_function(
        SlotType(func_to_connect).track_foreign(space));
}



template <int dim_, int range_, int rank_, int codim_, Transformation type_>
void
PhysicalSpace<dim_, range_, rank_, codim_, type_>::
rebuild_after_insert_knots(
    const SafeSTLArray<SafeSTLVector<Real>,dim> &knots_to_insert,
    const CartesianGrid<dim> &old_grid)
{
    Assert(false,ExcNotImplemented());

    /*
    this->ref_space_previous_refinement_ =
        shared_ptr<BSplineSpace<dim_,range_,rank_>>(new
                                                    BSplineSpace(
                                                        const_pointer_cast<SpaceData>(
                                                            this->space_data_->get_spline_space_previous_refinement()),
                                                        this->end_b_));


    this->dof_distribution_ = shared_ptr<DofDistribution<dim_,range_,rank_>>(
                                  new DofDistribution<dim_,range_,rank_>(
                                      this->space_data_->get_num_basis_table(),
                                      this->space_data_->get_degree(),
                                      this->space_data_->get_periodic_table()));

    operators_ = BernsteinExtraction<dim, range, rank>(
                     this->get_grid(),
                     this->space_data_->compute_knots_with_repetition(end_b_),
                     this->space_data_->accumulated_interior_multiplicities(),
                     this->space_data_->get_degree());
                     //*/
}


#ifdef SERIALIZATION
template <int dim_, int range_, int rank_, int codim_, Transformation type_>
template<class Archive>
void
PhysicalSpace<dim_, range_, rank_, codim_, type_>::
serialize(Archive &ar, const unsigned int version)
{
    ar &boost::serialization::make_nvp("PhysicalSpace_base_t",
                                       boost::serialization::base_object<base_t>(*this));

    ar.template register_type<BSplineSpace<dim_,range_,rank_> >();
    ar.template register_type<NURBSSpace<dim_,range_,rank_> >();
    ar &boost::serialization::make_nvp("ref_space_",ref_space_);
    Assert(ref_space_ != nullptr,ExcNullPtr());

    ar.template register_type<IgFunction<dim,0,dim+codim,1> >();
    ar &boost::serialization::make_nvp("map_func_",map_func_);
    Assert(map_func_ != nullptr,ExcNullPtr());

    ar &boost::serialization::make_nvp("phys_space_previous_refinement_",phys_space_previous_refinement_);
}

///@}
#endif


IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/physical_space.inst>


