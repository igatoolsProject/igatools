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

using std::shared_ptr;
using std::unique_ptr;

IGA_NAMESPACE_OPEN

template <int dim_>
SpaceBase<dim_>::
SpaceBase(shared_ptr<CartesianGrid<dim_>> grid)
    :
    base_t(grid),
    space_id_(UniqueIdGenerator::get_unique_id())
{};


template <int dim_>
Index
SpaceBase<dim_>::
get_space_id() const
{
    return space_id_;
}

#ifdef SERIALIZATION
template <int dim_>
template<class Archive>
void
SpaceBase<dim_>::
serialize(Archive &ar, const unsigned int version)
{
    ar &boost::serialization::make_nvp("SpaceBase_base_t",
                                       boost::serialization::base_object<base_t>(*this));

    ar &boost::serialization::make_nvp("space_id_",space_id_);
}
///@}
#endif // SERIALIZATION






template <int dim_,int codim_,int range_,int rank_>
Space<dim_,codim_,range_,rank_>::
Space(shared_ptr<CartesianGrid<dim_>> grid,
      const shared_ptr<MapFunc> &map_func)
    :
    base_t(grid)
{
    Assert(map_func != nullptr, ExcNullPtr());
    Assert(map_func.unique(), ExcNotUnique());

    map_func_.swap(const_cast<shared_ptr<MapFunc> &>(map_func));
    Assert(this->get_grid() == this->map_func_->get_grid(),
           ExcMessage("Reference space and mapping grids are not the same."))
}


template <int dim_,int codim_,int range_,int rank_>
auto
Space<dim_,codim_,range_,rank_>::
begin(const std::string &element_property) const -> ElementIterator
{
    return ElementIterator(
               this->create_element(
                   this->get_grid()->get_first_element_id_same_property(element_property)),
               element_property);
}



template <int dim_,int codim_,int range_,int rank_>
auto
Space<dim_,codim_,range_,rank_>::
last(const std::string &element_property) const -> ElementIterator
{
    return ElementIterator(
               this->create_element(
                   this->get_grid()->get_first_element_id_same_property(element_property)),
               element_property);
}



template <int dim_,int codim_,int range_,int rank_>
auto
Space<dim_,codim_,range_,rank_>::
end(const std::string &element_property) const -> ElementIterator
{
    return ElementIterator(this->create_element(IteratorState::pass_the_end),element_property);
}



template <int dim_,int codim_,int range_,int rank_>
auto
Space<dim_,codim_,range_,rank_>::
get_num_basis() const -> Size
{
    return this->get_dof_distribution()->get_num_dofs_table().total_dimension();
}


template <int dim_,int codim_,int range_,int rank_>
auto
Space<dim_,codim_,range_,rank_>::
get_num_basis(const int comp) const -> Size
{
    return this->get_dof_distribution()->get_num_dofs_table().get_component_size(comp);
}

template <int dim_,int codim_,int range_,int rank_>
auto
Space<dim_,codim_,range_,rank_>::
get_num_basis(const int comp, const int dir) const -> Size
{
    return this->get_dof_distribution()->get_num_dofs_table()[comp][dir];
}


template <int dim_,int codim_,int range_,int rank_>
auto
Space<dim_,codim_,range_,rank_>::
get_elem_num_basis() const -> Size
{
    return this->get_dof_distribution()->get_num_dofs_table().total_dimension();
}

template <int dim_,int codim_,int range_,int rank_>
auto
Space<dim_,codim_,range_,rank_>::
get_global_dof_id(const TensorIndex<dim> &tensor_index,
                  const Index comp) const -> Index
{
    return this->get_dof_distribution()->get_index_table()[comp](tensor_index);
}

template <int dim_,int codim_,int range_,int rank_>
auto
Space<dim_,codim_,range_,rank_>::
get_interior_dofs() const -> std::set<Index>
{
    return this->get_dof_distribution()->get_interior_dofs();
}

template <int dim_,int codim_,int range_,int rank_>
auto
Space<dim_,codim_,range_,rank_>::
get_boundary_dofs(const int s_id, const topology_variant &topology) const -> std::set<Index>
{
    return this->get_dof_distribution()->get_boundary_dofs(s_id,topology);
}


#ifdef SERIALIZATION
template <int dim_,int codim_,int range_,int rank_>
template<class Archive>
void
Space<dim_,codim_,range_,rank_>::
serialize(Archive &ar, const unsigned int version)
{
//  ar.template register_type<BSplineSpace<dim_,range_,rank_>>();
//  ar.template register_type<NURBSSpace<dim_,range_,rank_>>();
    ar &boost::serialization::make_nvp("Space_base_t",
                                       boost::serialization::base_object<base_t>(*this));

    ar.template register_type<IgFunction<dim_,0,dim_+codim_,1> >();
    ar.template register_type<IdentityFunction<dim_,dim_> >();
    ar &boost::serialization::make_nvp("map_func_",map_func_);
    Assert(map_func_ != nullptr,ExcNullPtr());

}
///@}
#endif // SERIALIZATION

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/space.inst>
