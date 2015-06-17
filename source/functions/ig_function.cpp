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

#include <igatools/functions/ig_function.h>
#include <igatools/functions/function_element.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/space_tools.h>

using std::shared_ptr;

IGA_NAMESPACE_OPEN

template<int dim,int codim,int range,int rank>
IgFunction<dim,codim,range,rank>::
IgFunction(std::shared_ptr<const Space<dim,codim,range,rank>> space,
           std::shared_ptr<const EpetraTools::Vector> coeff,
           const std::string &property)
    :
    parent_t::Function(space->get_grid()),
    space_(space),
    property_(property),
    space_elem_(space->begin()),
    space_elem_handler_(space->get_elem_handler())
{
    Assert(space_ != nullptr, ExcNullPtr());
    Assert(coeff != nullptr,ExcNullPtr());

    const auto &dof_distribution = *(space_->get_dof_distribution());
    const auto &active_dofs = dof_distribution.get_dofs_id_same_property(property);

    const auto &c = *coeff;

    const auto &epetra_map = c.Map();

    for (const auto glob_dof : active_dofs)
    {
        auto loc_id = epetra_map.LID(glob_dof);
        Assert(loc_id >= 0,
               ExcMessage("Global dof " + std::to_string(glob_dof) + " not present in the input EpetraTools::Vector."));
        coeff_[glob_dof] = c[loc_id];
    }
}


template<int dim,int codim,int range,int rank>
IgFunction<dim,codim,range,rank>::
IgFunction(std::shared_ptr<const Space<dim,codim,range,rank>> space,
           const IgCoefficients &coeff,
           const std::string &property)
    :
    parent_t::Function(space->get_grid()),
    space_(space),
    property_(property),
    space_elem_(space->begin()),
    space_elem_handler_(space->get_elem_handler())
{
    Assert(space_ != nullptr, ExcNullPtr());

#ifndef NDEBUG
    const auto &dof_distribution = *(space_->get_dof_distribution());
    const auto &active_dofs = dof_distribution.get_dofs_id_same_property(property);

    for (const auto glob_dof : active_dofs)
        coeff_[glob_dof] = coeff.at(glob_dof);
#else
    coeff_ = coeff;
#endif
}



template<int dim,int codim,int range,int rank>
IgFunction<dim,codim,range,rank>::
IgFunction(const self_t &fun)
    :
    parent_t::Function(fun.space_->get_grid()),
    space_(fun.space_),
    coeff_(fun.coeff_),
    property_(fun.property_),
    space_elem_(fun.space_->begin()),
    space_elem_handler_(fun.space_->get_elem_handler())
{}



template<int dim,int codim,int range,int rank>
auto
IgFunction<dim,codim,range,rank>::
create(std::shared_ptr<const Space<dim,codim,range,rank>> space,
       std::shared_ptr<const EpetraTools::Vector> coeff,
       const std::string &property) ->  std::shared_ptr<self_t>
{
    auto ig_func = std::make_shared<self_t>(space, coeff, property);

    Assert(ig_func != nullptr, ExcNullPtr());
#ifdef MESH_REFINEMENT
    ig_func->create_connection_for_insert_knots(ig_func);
#endif // MESH_REFINEMENT
    return ig_func;
}


template<int dim,int codim,int range,int rank>
auto
IgFunction<dim,codim,range,rank>::
create(std::shared_ptr<const Space<dim,codim,range,rank>> space,
       const IgCoefficients &coeff,
       const std::string &property) ->  std::shared_ptr<self_t>
{
    auto ig_func = std::make_shared<self_t>(space, coeff, property);

    Assert(ig_func != nullptr, ExcNullPtr());
#ifdef MESH_REFINEMENT
    ig_func->create_connection_for_insert_knots(ig_func);
#endif // MESH_REFINEMENT
    return ig_func;
}



template<int dim,int codim,int range,int rank>
void
IgFunction<dim,codim,range,rank>::
reset(const ValueFlags &flag, const eval_pts_variant &eval_pts)
{
    const std::set<int> elems_id =
        this->get_ig_space()->get_grid()->get_elements_id();

    this->reset_selected_elements(
        flag,
        eval_pts,
        SafeSTLVector<Index>(elems_id.begin(),elems_id.end()));
}



template<int dim,int codim,int range,int rank>
void
IgFunction<dim,codim,range,rank>::
reset_selected_elements(
    const ValueFlags &flag_in,
    const eval_pts_variant &eval_pts,
    const SafeSTLVector<Index> &elements_flat_id)
{
    parent_t::reset(flag_in, eval_pts);

    auto reset_dispatcher = ResetDispatcher(
                                flag_in,elements_flat_id,*space_elem_handler_,this->flags_);
    boost::apply_visitor(reset_dispatcher, eval_pts);
}



template<int dim,int codim,int range,int rank>
auto
IgFunction<dim,codim,range,rank>::
init_cache(ElementAccessor &elem, const topology_variant &k) -> void
{
    parent_t::init_cache(elem, k);

    auto init_cache_dispatcher = InitCacheDispatcher(*space_elem_handler_,*space_elem_);
    boost::apply_visitor(init_cache_dispatcher, k);
}



template<int dim,int codim,int range,int rank>
auto
IgFunction<dim,codim,range,rank>::
fill_cache(ElementAccessor &func_elem, const topology_variant &k, const int sub_elem_id) -> void
{
    parent_t::fill_cache(func_elem,k,sub_elem_id);

    space_elem_.move_to(func_elem.get_flat_index());

    const auto elem_dofs = space_elem_->get_local_to_global(property_);
    SafeSTLVector<Real> loc_coeff;
    for (auto elem_dof : elem_dofs)
        loc_coeff.push_back(coeff_[elem_dof]);


    auto fill_cache_dispatcher = FillCacheDispatcher(
        sub_elem_id,*this,*space_elem_handler_,func_elem,*space_elem_,loc_coeff,property_);

    boost::apply_visitor(fill_cache_dispatcher, k);

}



template<int dim,int codim,int range,int rank>
auto
IgFunction<dim,codim,range,rank>::
get_ig_space() const -> std::shared_ptr<const Space<dim,codim,range,rank>>
{
    return space_;
}



template<int dim,int codim,int range,int rank>
auto
IgFunction<dim,codim,range,rank>::
get_coefficients() const -> const CoeffType &
{
    return coeff_;
}



template<int dim,int codim,int range,int rank>
auto
IgFunction<dim,codim,range,rank>::
operator +=(const self_t &fun) -> self_t &
{
    Assert(space_ == fun.space_,
    ExcMessage("Functions defined on different spaces."));

    for (const auto &f_dof_value : fun.coeff_)
        coeff_[f_dof_value.first] += f_dof_value.second;

    return *this;
}


#ifdef MESH_REFINEMENT
template<int dim,int codim,int range,int rank>
void
IgFunction<dim,codim,range,rank>::
rebuild_after_insert_knots(
    const SafeSTLArray<SafeSTLVector<Real>,dim> &knots_to_insert,
    const CartesianGrid<dim> &grid_old)
{
    using std::const_pointer_cast;
    this->function_previous_refinement_ =
    		IgFunction<dim,codim,range,rank>::create(
                                              const_pointer_cast<Space<dim,codim,range,rank>>(space_->get_space_previous_refinement()),
                                              coeff_,
                                              property_);
/*
    LogStream out;
    out.begin_item("function_previous_refinement_");
    this->function_previous_refinement_->print_info(out);
    out.end_item();
//*/

    const int max_degree = space_->get_max_degree();
//    const int max_degree = 10;
    QGauss<dim> quad(max_degree+1);
    auto function_refined = space_tools::projection_l2(
                                this->function_previous_refinement_,
                                const_pointer_cast<const Space<dim,codim,range,rank>>(space_),
                                quad);
/*
    out.begin_item("function_refined");
    function_refined->print_info(out);
    out.end_item();
//*/

    this->coeff_ = std::move(function_refined->coeff_);
//    this->property_ = DofProperties::active;

}

template<int dim,int codim,int range,int rank>
void
IgFunction<dim,codim,range,rank>::
create_connection_for_insert_knots(std::shared_ptr<self_t> ig_function)
{
    Assert(ig_function != nullptr, ExcNullPtr());
    Assert(&(*ig_function) == &(*this), ExcMessage("Different objects."));

    using SlotType = typename CartesianGrid<dim>::SignalInsertKnotsSlot;

    auto func_to_connect =
        std::bind(&self_t::rebuild_after_insert_knots,
                  ig_function.get(),
                  std::placeholders::_1,
                  std::placeholders::_2);
    this->grid_wrapper_.connect_insert_knots_function(
        SlotType(func_to_connect).track_foreign(ig_function));
}
#endif // MESH_REFINEMENT

template<int dim,int codim,int range,int rank>
void
IgFunction<dim,codim,range,rank>::
print_info(LogStream &out) const
{
    using std::to_string;
    const std::string template_args= "<" + to_string(dim) + "," + to_string(codim) + ","
                                     + to_string(range) + "," + to_string(rank) + ">";

    out.begin_item("IgFunction" + template_args);

    out.begin_item("Function" + template_args);
    parent_t::print_info(out);
    out.end_item();

    out.begin_item("Reference space info:");
    space_->print_info(out);
    out.end_item();
    out << std::endl;

    out.begin_item("Coefficients (a.k.a. \"control values\"):");
    coeff_.print_info(out);
    out.end_item();

    out.end_item();
}

template<int dim,int codim,int range,int rank>
const std::string &
IgFunction<dim,codim,range,rank>::
get_property() const
{
    return property_;
}

#ifdef SERIALIZATION
template<int dim,int codim,int range,int rank>
template<class Archive>
void
IgFunction<dim,codim,range,rank>::
serialize(Archive &ar, const unsigned int version)
{
    ar &boost::serialization::make_nvp("IgFunction_base_t",
                                       boost::serialization::base_object<base_t>(*this));

    //TODO (martinelli, May 18, 2105): register PhysicalSpace, PhysicalSpaceElement and PhysSpaceElementHandler
    ar.template register_type<BSplineSpace<dim,range,rank>>();
    ar.template register_type<NURBSSpace<dim,range,rank>>();
//    ar.template register_type<PhysicalSpace<dim,range,rank,0,Transformation::h_grad>>();
    auto non_nonst_space = std::const_pointer_cast<Space<dim,codim,range,rank>>(space_);
    ar &boost::serialization::make_nvp("space_",non_nonst_space);
    space_ = non_nonst_space;
    Assert(space_ != nullptr,ExcNullPtr());

    ar &boost::serialization::make_nvp("coeff_",coeff_);

    ar &boost::serialization::make_nvp("property_",const_cast<std::string &>(property_));

    ar.template register_type<BSplineElement<dim,range,rank>>();
    ar.template register_type<NURBSElement<dim,range,rank>>();
    ar.template register_type<PhysicalSpaceElement<dim,range,rank,codim>>();
    ar &boost::serialization::make_nvp("space_elem_",space_elem_);


    ar.template register_type<BSplineElementHandler<dim,range,rank>>();
    ar.template register_type<NURBSElementHandler<dim,range,rank>>();
    ar.template register_type<PhysSpaceElementHandler<dim,range,rank,codim>>();
    ar &boost::serialization::make_nvp("space_elem_handler_",space_elem_handler_);
    Assert(space_elem_handler_ != nullptr,ExcNullPtr());
}
#endif // SERIALIZATION


IGA_NAMESPACE_CLOSE

#include <igatools/functions/ig_function.inst>


