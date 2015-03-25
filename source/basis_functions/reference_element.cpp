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


#include <igatools/basis_functions/reference_element.h>
#include <igatools/basis_functions/reference_space.h>
#include <igatools/basis_functions/space_element.h>
#include <igatools/basis_functions/reference_element_handler.h>


IGA_NAMESPACE_OPEN



template <int dim, int range, int rank>
ReferenceElement<dim, range, rank>::
ReferenceElement(const std::shared_ptr<ConstSpace> space,
                 const Index elem_index)
    :
    parent_t(space,elem_index)
{
    Assert(this->get_space() != nullptr,ExcNullPtr());

    //-------------------------------------------------
    const auto &degree_table = space->get_degree();
    TensorSizeTable n_basis(degree_table.get_comp_map());
    for (auto comp : degree_table.get_active_components_id())
        n_basis[comp] = TensorSize<dim>(degree_table[comp]+1);

    n_basis_direction_ = n_basis;
    //-------------------------------------------------


    //----------------------------------------------------------------
    comp_offset_[0] = 0;
    for (int comp = 1; comp < Space::n_components; ++comp)
        comp_offset_[comp] = comp_offset_[comp-1] +
                             this->n_basis_direction_.get_component_size(comp-1);
    //----------------------------------------------------------------


    //----------------------------------------------------------------
    for (int comp : basis_functions_indexer_.get_active_components_id())
    {
        // creating the objects for fast conversion from flat-to-tensor indexing
        // (in practice it is an hash-table from flat to tensor indices)
        basis_functions_indexer_[comp] =
            std::shared_ptr<Indexer>(new Indexer(this->n_basis_direction_[comp]));
    }
    //----------------------------------------------------------------
};


template <int dim, int range, int rank>
ReferenceElement<dim, range, rank>::
ReferenceElement(const std::shared_ptr<ConstSpace> space,
                 const TensorIndex<dim> &elem_index)
    :
    ReferenceElement(space,space->get_grid()->tensor_to_flat(elem_index))
{}


template <int dim, int range, int rank>
ReferenceElement<dim, range, rank>::
ReferenceElement(const ReferenceElement<dim,range,rank> &elem,
                 const iga::CopyPolicy &copy_policy)
    :
    parent_t(elem,copy_policy),
    n_basis_direction_(elem.n_basis_direction_),
    comp_offset_(elem.comp_offset_),
    basis_functions_indexer_(elem.basis_functions_indexer_)
{};


template <int dim, int range, int rank>
void
ReferenceElement<dim, range, rank>::
move_to(const Index flat_index)
{
//    parent_t::move_to(flat_index);
    this->as_cartesian_grid_element_accessor().move_to(flat_index);
}

template <int dim, int range, int rank>
template <int deriv_order>
auto
ReferenceElement<dim, range, rank>::
evaluate_basis_derivatives_at_points(
    const Quadrature<dim> &points,
    const std::string &dofs_property) ->
ValueTable<
Conditional< deriv_order==0,
             Value,
             Derivative<deriv_order> > >
{
    auto elem_handler = ReferenceElementHandler<dim,range,rank>::create(this->get_space());

    ValueFlags flags;
    if (deriv_order == 0)
        flags = ValueFlags::value;
    else if (deriv_order == 1)
        flags = ValueFlags::gradient;
    else if (deriv_order == 2)
        flags = ValueFlags::hessian;
    else
    {
        Assert(false,ExcNotImplemented());
    }

    elem_handler->reset_one_element(flags,points,this->get_flat_index());
    elem_handler->template init_cache<dim>(*this);
    elem_handler->template fill_cache<dim>(*this,0);

//    Assert(false,ExcNotImplemented());

    return this->template get_values<deriv_order,dim>(0,dofs_property);
}



template <int dim, int range, int rank>
int
ReferenceElement<dim, range, rank>::
get_num_basis_comp(const int i) const
{
    return this->n_basis_direction_[i].flat_size();
}

template <int dim, int range, int rank>
auto
ReferenceElement<dim, range, rank>::
get_basis_offset() const -> OffsetTable
{
    return this->comp_offset_;
}


template <int dim, int range, int rank>
int
ReferenceElement<dim, range, rank>::
get_num_basis() const
{
    return this->n_basis_direction_.total_dimension();
}


template <int dim, int range, int rank>
void
ReferenceElement<dim, range, rank>::
print_info(LogStream &out) const
{
    parent_t::print_info(out);
    out.begin_item("Number of element basis: ");
    n_basis_direction_.print_info(out);
    out.end_item();
}

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/reference_element.inst>


