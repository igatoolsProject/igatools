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
#include <igatools/functions/ig_function_handler.h>
#include <igatools/functions/function_element.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/physical_space_basis.h>
#include <igatools/basis_functions/space_tools.h>


using std::shared_ptr;

IGA_NAMESPACE_OPEN

template<int dim,int codim,int range,int rank>
IgFunction<dim,codim,range,rank>::
IgFunction(const SharedPtrConstnessHandler<PhysBasis> &space,
           const EpetraTools::Vector &coeff,
           const std::string &dofs_property,
           const std::string &name)
  :
  parent_t::Function(
   space.data_is_const() ?
   SharedPtrConstnessHandler<DomainType>(space.get_ptr_const_data()->get_physical_domain()) :
   SharedPtrConstnessHandler<DomainType>(space.get_ptr_data()->get_physical_domain()),
   name),
  basis_(space),
  dofs_property_(dofs_property)
{
  const auto &dof_distribution = *(basis_->get_spline_space()->get_dof_distribution());
  const auto &active_dofs = dof_distribution.get_global_dofs(dofs_property);

  const auto &epetra_map = coeff.Map();

  for (const auto glob_dof : active_dofs)
  {
    auto loc_id = epetra_map.LID(glob_dof);
    Assert(loc_id >= 0,
           ExcMessage("Global dof " + std::to_string(glob_dof) + " not present in the input EpetraTools::Vector."));
    coeff_[glob_dof] = coeff[loc_id];
  }
}


template<int dim,int codim,int range,int rank>
IgFunction<dim,codim,range,rank>::
IgFunction(const SharedPtrConstnessHandler<PhysBasis> &space,
           const IgCoefficients &coeff,
           const std::string &dofs_property,
           const std::string &name)
  :
  parent_t::Function(
   space.data_is_const() ?
   SharedPtrConstnessHandler<DomainType>(space.get_ptr_const_data()->get_physical_domain()) :
   SharedPtrConstnessHandler<DomainType>(space.get_ptr_data()->get_physical_domain()),
   name),
  basis_(space),
  dofs_property_(dofs_property)
{

#ifndef NDEBUG
  const auto &dof_distribution = *(basis_->get_spline_space()->get_dof_distribution());
  const auto &active_dofs = dof_distribution.get_global_dofs(dofs_property);

  for (const auto glob_dof : active_dofs)
    coeff_[glob_dof] = coeff.at(glob_dof);
#else
  coeff_ = coeff;
#endif
}



template<int dim,int codim,int range,int rank>
auto
IgFunction<dim,codim,range,rank>::
const_create(const std::shared_ptr<const PhysBasis> &space,
             const EpetraTools::Vector &coeff,
             const std::string &dofs_property,
             const std::string &name) ->  std::shared_ptr<const parent_t>
{
  auto ig_func = std::make_shared<self_t>(SharedPtrConstnessHandler<PhysBasis>(space),
  coeff, dofs_property,name);
  Assert(ig_func != nullptr, ExcNullPtr());

  return ig_func;
}


template<int dim,int codim,int range,int rank>
auto
IgFunction<dim,codim,range,rank>::
const_create(const std::shared_ptr<const PhysBasis> &space,
             const IgCoefficients &coeff,
             const std::string &dofs_property,
             const std::string &name) ->  std::shared_ptr<const parent_t>
{
  auto ig_func = std::make_shared<self_t>(SharedPtrConstnessHandler<PhysBasis>(space),
  coeff, dofs_property,name);
  Assert(ig_func != nullptr, ExcNullPtr());

  return ig_func;
}


template<int dim,int codim,int range,int rank>
auto
IgFunction<dim,codim,range,rank>::
create(const std::shared_ptr<PhysBasis> &space,
       const EpetraTools::Vector &coeff,
       const std::string &dofs_property,
       const std::string &name) ->  std::shared_ptr<parent_t>
{
  auto ig_func = std::make_shared<self_t>(SharedPtrConstnessHandler<PhysBasis>(space),
  coeff, dofs_property,name);

  Assert(ig_func != nullptr, ExcNullPtr());
#ifdef MESH_REFINEMENT
  ig_func->create_connection_for_insert_knots(ig_func);
#endif // MESH_REFINEMENT
  return ig_func;
}


template<int dim,int codim,int range,int rank>
auto
IgFunction<dim,codim,range,rank>::
create(const std::shared_ptr<PhysBasis> &space,
       const IgCoefficients &coeff,
       const std::string &dofs_property,
       const std::string &name) ->  std::shared_ptr<parent_t>
{
  auto ig_func = std::make_shared<self_t>(SharedPtrConstnessHandler<PhysBasis>(space),
  coeff, dofs_property,name);
  Assert(ig_func != nullptr, ExcNullPtr());

#ifdef MESH_REFINEMENT
  ig_func->create_connection_for_insert_knots(ig_func);
#endif // MESH_REFINEMENT

  return ig_func;
}



template<int dim,int codim,int range,int rank>
auto
IgFunction<dim,codim,range,rank>::
create_cache_handler() const
-> std::unique_ptr<typename parent_t::ElementHandler>
{
  return std::make_unique<IgFunctionHandler<dim,codim,range,rank>>(
    std::dynamic_pointer_cast<const self_t>(this->shared_from_this()));
}


template<int dim,int codim,int range,int rank>
auto
IgFunction<dim,codim,range,rank>::
get_basis() const -> std::shared_ptr<const PhysBasis>
{
  return basis_.get_ptr_const_data();
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
  Assert(basis_.get_ptr_const_data() == fun.basis_.get_ptr_const_data(),
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
  const Grid<dim> &grid_old)
{
  using std::const_pointer_cast;
  this->function_previous_refinement_ =
    IgFunction<dim,codim,range,rank>::const_create(
      std::dynamic_pointer_cast<const PhysBasis>(basis_->get_basis_previous_refinement()),
      coeff_,
      dofs_property_);


  const int max_degree = basis_->get_spline_space()->get_max_degree();

  const auto quad = QGauss<dim>::create(max_degree+1);
  this->coeff_ =
    space_tools::projection_l2_function<dim,codim,range,rank>(
      *(this->function_previous_refinement_),*basis_,quad);

//  this->coeff_ = std::move(function_refined->coeff_);
}

template<int dim,int codim,int range,int rank>
void
IgFunction<dim,codim,range,rank>::
create_connection_for_insert_knots(std::shared_ptr<self_t> ig_function)
{
  Assert(ig_function != nullptr, ExcNullPtr());
  Assert(&(*ig_function) == &(*this), ExcMessage("Different objects."));

  using SlotType = typename Grid<dim>::SignalInsertKnotsSlot;

  auto func_to_connect =
    std::bind(&self_t::rebuild_after_insert_knots,
              ig_function.get(),
              std::placeholders::_1,
              std::placeholders::_2);
  /*
  this->domain_.get_ptr_data()->get_grid_function()->get_grid()->connect_insert_knots(
    SlotType(func_to_connect).track_foreign(ig_function));
    //*/
  this->domain_.get_ptr_data()->
  connect_insert_knots(SlotType(func_to_connect).track_foreign(ig_function));
//  Assert(false,ExcNotImplemented());
}

#endif // MESH_REFINEMENT

template<int dim,int codim,int range,int rank>
void
IgFunction<dim,codim,range,rank>::
print_info(LogStream &out) const
{
  using std::to_string;
  const std::string func_template_args= "<" +
                                        to_string(dim) + "," + to_string(codim) + "," +
                                        to_string(range) + "," + to_string(rank) + ">";

  out.begin_item("IgFunction" + func_template_args);

  const std::string basis_template_args= "<" +
                                         to_string(dim) + "," + to_string(range) + "," +
                                         to_string(rank) + "," + to_string(codim) + ">";

  out.begin_item("PhysicalSpaceBasis" + basis_template_args);
  basis_->print_info(out);
  out.end_item();


  out.begin_item("IgCoefficients:");
  coeff_.print_info(out);
  out.end_item();

  out << "Dofs property: " << dofs_property_ << std::endl;

  out.end_item();
}



template<int dim,int codim,int range,int rank>
const std::string &
IgFunction<dim,codim,range,rank>::
get_dofs_property() const
{
  return dofs_property_;
}

#if 0
#ifdef SERIALIZATION
template<int dim,int codim,int range,int rank>
template<class Archive>
void
IgFunction<dim,codim,range,rank>::
serialize(Archive &ar, const unsigned int version)
{
  ar &boost::serialization::make_nvp("IgFunction_base_t",
                                     boost::serialization::base_object<base_t>(*this));

  ar.template register_type<BSpline<dim,range,rank>>();

#ifdef USE_NURBS
  ar.template register_type<NURBS<dim,range,rank>>();
#endif // NURBS

  ar.template register_type<PhysicalSpaceBasis<dim,range,rank,codim>>();
  auto non_nonst_space = std::const_pointer_cast<Basis<dim,codim,range,rank>>(basis_);
  ar &boost::serialization::make_nvp("basis_",non_nonst_space);
  basis_ = non_nonst_space;
  Assert(basis_ != nullptr,ExcNullPtr());

  ar &boost::serialization::make_nvp("coeff_",coeff_);

  ar &boost::serialization::make_nvp("dofs_property_",const_cast<std::string &>(dofs_property_));

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
#endif

IGA_NAMESPACE_CLOSE

#include <igatools/functions/ig_function.inst>

