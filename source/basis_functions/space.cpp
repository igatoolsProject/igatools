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

#include <igatools/basis_functions/space.h>
#include <igatools/utils/unique_id_generator.h>
#include <igatools/basis_functions/space_element.h>
#include <igatools/functions/ig_function.h>
#include <igatools/functions/identity_function.h>
#include <igatools/basis_functions/dof_distribution.h>
#include <igatools/basis_functions/nurbs_space.h>



using std::shared_ptr;
using std::unique_ptr;

IGA_NAMESPACE_OPEN

template <int dim_>
SpaceBase<dim_>::
SpaceBase(const shared_ptr<const CartesianGrid<dim_>> &grid)
    :
    object_id_(UniqueIdGenerator::get_unique_id()),
    grid_(grid)
{};

template <int dim_>
SpaceBase<dim_>::
SpaceBase(const shared_ptr<CartesianGrid<dim_>> &grid)
    :
    object_id_(UniqueIdGenerator::get_unique_id()),
    grid_(grid)
{};

template <int dim_>
Index
SpaceBase<dim_>::
get_object_id() const
{
    return object_id_;
}

template <int dim_>
std::shared_ptr<CartesianGrid<dim_> >
SpaceBase<dim_>::
get_ptr_grid()
{
    return grid_.get_ptr_data();
}

template <int dim_>
std::shared_ptr<const CartesianGrid<dim_> >
SpaceBase<dim_>::
get_ptr_const_grid() const
{
    return grid_.get_ptr_const_data();
}


template <int dim_>
const std::string &
SpaceBase<dim_>::
get_name() const
{
    return name_;
}

template <int dim_>
void
SpaceBase<dim_>::
set_name(const std::string &name)
{
    name_ = name;
}


#ifdef MESH_REFINEMENT

template <int dim_>
void
SpaceBase<dim_>::
refine_h(const Size n_subdivisions)
{
    this->get_ptr_grid()->refine(n_subdivisions);
}

#endif // MESH_REFINEMENT

#ifdef SERIALIZATION
template <int dim_>
template<class Archive>
void
SpaceBase<dim_>::
serialize(Archive &ar, const unsigned int version)
{
    ar &boost::serialization::make_nvp("grid_",grid_);

    ar &boost::serialization::make_nvp("object_id_",object_id_);

    ar &boost::serialization::make_nvp("name_",name_);

}
///@}
#endif // SERIALIZATION






template <int dim_,int codim_,int range_,int rank_,Transformation type_>
Space<dim_,codim_,range_,rank_,type_>::
Space(const shared_ptr<CartesianGrid<dim_>> &grid,
      const shared_ptr<MapFunc> &map_func)
    :
    base_t(grid),
    phys_domain_(std::make_shared<PhysDomain>(map_func))
{
    Assert(map_func != nullptr, ExcNullPtr());
//    map_func_.get_ref_ptr_data().swap(const_cast<shared_ptr<MapFunc> &>(map_func));
//    Assert(map_func_.unique(), ExcNotUnique());

    Assert(phys_domain_ != nullptr,ExcNullPtr());

    Assert(this->get_ptr_grid() == phys_domain_->get_grid(),
           ExcMessage("Reference space and mapping grids are not the same."));

}

template <int dim_,int codim_,int range_,int rank_,Transformation type_>
Space<dim_,codim_,range_,rank_,type_>::
Space(const shared_ptr<const CartesianGrid<dim_>> &grid,
      const shared_ptr<MapFunc> &map_func)
    :
    base_t(grid),
    phys_domain_(std::make_shared<PhysDomain>(map_func))
{
    Assert(map_func != nullptr, ExcNullPtr());
//    map_func_.get_ref_ptr_data().swap(const_cast<shared_ptr<MapFunc> &>(map_func));
//    Assert(map_func_.unique(), ExcNotUnique());


    Assert(phys_domain_ != nullptr,ExcNullPtr());

    Assert(this->get_ptr_const_grid() == phys_domain_->get_grid(),
           ExcMessage("Reference space and mapping grids are not the same."));
}


template <int dim_,int codim_,int range_,int rank_,Transformation type_>
auto
Space<dim_,codim_,range_,rank_,type_>::
begin(const PropId &prop) -> ElementIterator
{
    return ElementIterator(this->shared_from_this(),
    this->get_ptr_grid()->get_element_property(prop).begin(),
    prop);
}



template <int dim_,int codim_,int range_,int rank_,Transformation type_>
auto
Space<dim_,codim_,range_,rank_,type_>::
end(const PropId &prop) -> ElementIterator
{
    return ElementIterator(this->shared_from_this(),
    this->get_ptr_grid()->get_element_property(prop).end(),
    prop);
}



template <int dim_,int codim_,int range_,int rank_,Transformation type_>
auto
Space<dim_,codim_,range_,rank_,type_>::
get_num_basis() const -> Size
{
    return this->get_ptr_const_dof_distribution()->get_num_dofs_table().total_dimension();
}


template <int dim_,int codim_,int range_,int rank_,Transformation type_>
auto
Space<dim_,codim_,range_,rank_,type_>::
get_num_basis(const int comp) const -> Size
{
    return this->get_ptr_const_dof_distribution()->get_num_dofs_table().get_component_size(comp);
}

template <int dim_,int codim_,int range_,int rank_,Transformation type_>
auto
Space<dim_,codim_,range_,rank_,type_>::
get_num_basis(const int comp, const int dir) const -> Size
{
    return this->get_ptr_const_dof_distribution()->get_num_dofs_table()[comp][dir];
}


template <int dim_,int codim_,int range_,int rank_,Transformation type_>
auto
Space<dim_,codim_,range_,rank_,type_>::
get_elem_num_basis() const -> Size
{
    return this->get_ptr_const_dof_distribution()->get_num_dofs_table().total_dimension();
}

template <int dim_,int codim_,int range_,int rank_,Transformation type_>
auto
Space<dim_,codim_,range_,rank_,type_>::
get_global_dof_id(const TensorIndex<dim> &tensor_index,
                  const Index comp) const -> Index
{
    return this->get_ptr_const_dof_distribution()->get_index_table()[comp](tensor_index);
}

template <int dim_,int codim_,int range_,int rank_,Transformation type_>
auto
Space<dim_,codim_,range_,rank_,type_>::
get_interior_dofs() const -> std::set<Index>
{
    return this->get_ptr_const_dof_distribution()->get_interior_dofs();
}

template <int dim_,int codim_,int range_,int rank_,Transformation type_>
auto
Space<dim_,codim_,range_,rank_,type_>::
get_boundary_dofs(const int s_id, const topology_variant &topology) const -> std::set<Index>
{
    return this->get_ptr_const_dof_distribution()->get_boundary_dofs(s_id,topology);
}


#ifdef SERIALIZATION
template <int dim_,int codim_,int range_,int rank_,Transformation type_>
template<class Archive>
void
Space<dim_,codim_,range_,rank_,type_>::
serialize(Archive &ar, const unsigned int version)
{
    ar.template register_type<BSplineSpace<dim_,range_,rank_>>();

#ifdef NURBS
    ar.template register_type<NURBSSpace<dim_,range_,rank_>>();
#endif

    ar.template register_type<PhysicalSpace<dim_,range_,rank_,codim_,Transformation::h_grad>>();

    ar &boost::serialization::make_nvp("Space_base_t",
                                       boost::serialization::base_object<base_t>(*this));

    ar.template register_type<IgFunction<dim_,0,dim_+codim_,1> >();
    ar.template register_type<IdentityFunction<dim_,dim_> >();
    ar &boost::serialization::make_nvp("map_func_",map_func_);
//    Assert(map_func_ != nullptr,ExcNullPtr());

}
///@}
#endif // SERIALIZATION

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/space.inst>
