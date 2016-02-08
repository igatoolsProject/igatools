//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2016  by the igatools authors (see authors.txt).
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

#include <igatools/operators/interface_operator.h>
#include <igatools/base/exceptions.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/grid_tools.h>
//#include <igatools/basis_functions/physical_space.h>


IGA_NAMESPACE_OPEN
#if 0
using namespace iga;
using std::shared_ptr;
using std::array;


template <class PhysSpace>
InterfaceOperator<PhysSpace>::
InterfaceOperator(const shared_ptr<const PhysSpace> space_0,
                  const Index &face_id_0,
                  const shared_ptr<const PhysSpace> space_1,
                  const Index &face_id_1,
                  const shared_ptr<Form> face_form_trial,
                  const shared_ptr<Form> face_form_test,
                  const Real &constant)
  :
  space_0_(space_0),
  space_1_(space_1),
  face_id_0_(face_id_0),
  face_id_1_(face_id_1),
  face_form_trial_(face_form_trial),
  face_form_test_(face_form_test),
  integration_constant_(constant)
{
  Assert(space_0 != nullptr, ExcNullPtr());
  Assert(space_1 != nullptr, ExcNullPtr());
  Assert(face_id_0_ >= 0 && face_id_0_ < dim_ * 2,
         ExcIndexRange(face_id_0_, 0, dim_ * 2));
  Assert(face_id_1_ >= 0 && face_id_1_ < dim_ * 2,
         ExcIndexRange(face_id_1_, 0, dim_ * 2));

  // TODO: to check here that the face spaces are compatible.
};



template <class PhysSpace>
void
InterfaceOperator<PhysSpace>::
evaluate(Constraints_ &constraints)
{

  vector<Index> face_dofs_0;
  const auto face_space_0 = space_0_->get_face_space(face_id_0_, face_dofs_0);

  vector<Index> face_dofs_1;
  const auto face_space_1 = space_1_->get_face_space(face_id_1_, face_dofs_1);

  //-------------------------------------------------------------------------------------------
  // Creating the quadrature on the face (using the maximum degree for all components and for all directions)
  // TODO: this should be done only with the active directions on both facesdfx

  const auto &face_0_degrees = face_space_0->get_reference_space()->get_degree();
  const auto &face_1_degrees = face_space_1->get_reference_space()->get_degree();

  shared_ptr<Quadrature<face_dim_>> face_quad;
  const auto &face_0_degrees_comp = face_0_degrees[0];

  int max_degree = 0;
  for (const auto &face_0_degrees_comp : face_0_degrees)
    for (const auto &face_0_degrees_comp_dir : face_0_degrees_comp)
      max_degree = std::max(max_degree, face_0_degrees_comp_dir);
  for (const auto &face_1_degrees_comp : face_1_degrees)
    for (const auto &face_1_degrees_comp_dir : face_1_degrees_comp)
      max_degree = std::max(max_degree, face_1_degrees_comp_dir);

  TensorSize<face_dim_> n_points(face_0_degrees_comp + 1);
  face_quad = QGauss<face_dim_>::create(n_points);

  const QGauss<dim_> elem_quad(max_degree + 1);

  auto elem_0_quad = elem_quad.collapse_to_face(face_id_0_);
  auto elem_1_quad = elem_quad.collapse_to_face(face_id_1_);
  const auto elem_0_quad_pts_unit_domain = elem_0_quad.get_points().get_flat_cartesian_product();
  const auto elem_1_quad_pts_unit_domain = elem_1_quad.get_points().get_flat_cartesian_product();

  const auto  n_qp                      = face_quad->get_num_points();
  const auto  face_w_unit_domain        = face_quad->get_weights().get_flat_tensor_product();
  const auto  face_quad_pts_unit_domain = face_quad->get_points().get_flat_cartesian_product();
  //-------------------------------------------------------------------------------------------


  //-------------------------------------------------------------------------------------------
  // Face joint grid construction -> i.e., the parametric merged mesh
  const auto grid_0 = space_0_->get_grid();
  const auto grid_1 = space_1_->get_grid();
  using FaceGridMap = typename CartesianGrid<dim_>::FaceGridMap;
  auto face_grid_to_grid_0 = std::make_shared<FaceGridMap> ();
  auto face_grid_to_grid_1 = std::make_shared<FaceGridMap> ();
  const auto face_0_grid = grid_0->get_face_grid(face_id_0_, *face_grid_to_grid_0);
  const auto face_1_grid = grid_1->get_face_grid(face_id_1_, *face_grid_to_grid_1);

  typename grid_tools::InterGridMap<face_dim_> joint_face_grid_to_face_0_grid;
  typename grid_tools::InterGridMap<face_dim_> joint_face_grid_to_face_1_grid;
  const auto face_grid = grid_tools::build_cartesian_grid_union(
                           *face_0_grid,
                           *face_1_grid,
                           joint_face_grid_to_face_0_grid,
                           joint_face_grid_to_face_1_grid);
  //-------------------------------------------------------------------------------------------


  //-------------------------------------------------------------------------------------------
  // Loop on the face joint grid, i.e., the merged grid, to lead the integration on the 0 side
  auto felem_j     = face_grid->begin();
  auto felem_j_end = face_grid->end();


  //-------------------------------------------------------------------------------------------
  for (; felem_j != felem_j_end ; ++felem_j)
  {
    //---------------------------------------------------------------------------------------
    auto face_quad_pts_ref = felem_j->transform_points_unit_to_reference(face_quad_pts_unit_domain);

    auto face_meas = static_cast<CartesianGridElement<face_dim_>&>(*felem_j).get_measure();

    // felem_0 and elem_0 -----------------begin
    const Index felem_0_id = joint_face_grid_to_face_0_grid[felem_j]->get_flat_index();

    Index elem_0_id;
    vector<Index> glob_face_0_dofs;
    PointsContainer_ curr_elem_0_quad_pts_unit;
    this->get_element_and_quad_on_face(
      face_quad_pts_ref,
      space_0_,
      face_space_0,
      face_id_0_,
      face_dofs_0,
      *face_grid_to_grid_0,
      felem_0_id,
      elem_0_id,
      glob_face_0_dofs,
      curr_elem_0_quad_pts_unit);

    auto elem_0  = space_0_->get_element(elem_0_id);
    // felem_0 and elem_0 -------------------end



    // felem_1 and elem_1 ---------------begin
    const Index felem_1_id = joint_face_grid_to_face_1_grid[felem_j]->get_flat_index();

    Index elem_1_id;
    vector<Index> glob_face_1_dofs;
    PointsContainer_ curr_elem_1_quad_pts_unit;
    this->get_element_and_quad_on_face(
      face_quad_pts_ref,
      space_1_,
      face_space_1,
      face_id_1_,
      face_dofs_1,
      *face_grid_to_grid_1,
      felem_1_id,
      elem_1_id,
      glob_face_1_dofs,
      curr_elem_1_quad_pts_unit);

    auto elem_1 = space_1_->get_element(elem_1_id);
    // felem_1 and elem_1 ---------------end

    using FuncFaceMap = Function<face_dim_, space_dim_>;
    using GradientFaceMap = typename FuncFaceMap::Gradient;
    using ValFaceMap      = typename FuncFaceMap::Value;


    //get the mapping 0 from the physical space 0
    ValueVector<GradientFaceMap> face_0_map_grad(n_qp);
    ValueVector<ValFaceMap> face_0_map_val(n_qp);
    auto face_0_map = face_space_0->get_push_forward()->get_mapping();
    face_0_map->evaluate_gradients_at_points(face_quad_pts_ref, face_0_map_grad);
    face_0_map->evaluate_at_points(face_quad_pts_ref, face_0_map_val);

    auto face_1_map = face_space_0->get_push_forward()->get_mapping();
    ValueVector<ValFaceMap> face_1_map_val(n_qp);
    face_1_map->evaluate_at_points(face_quad_pts_ref, face_1_map_val);


    //------------------------------------------------------------------------------------
    ValueVector<Real> w_0_at_points(n_qp);
    const Real eps_map_check = std::pow(10, -8);
    for (int qp = 0; qp < n_qp; ++qp)
    {
      w_0_at_points[qp] =
        face_meas * face_w_unit_domain[qp] *
        determinant<face_dim_, space_dim_>(face_0_map_grad[qp]);

      auto vv = face_0_map_val[qp] - face_1_map_val[qp];
      Assert(vv.norm() < eps_map_check,
             ExcMessage("Maps are not geometricaly conforming at the chosen faces."))
    }
    //------------------------------------------------------------------------------------

    InterfaceElementMap elem_map_0(elem_0_id, felem_j->get_flat_index(),
                                   face_id_0_, curr_elem_0_quad_pts_unit,
                                   face_quad_pts_unit_domain);
    InterfaceElementMap elem_map_1(elem_1_id, felem_j->get_flat_index(),
                                   face_id_1_, curr_elem_1_quad_pts_unit,
                                   face_quad_pts_unit_domain);


    const auto vals_trial  = face_form_trial_->evaluate(elem_map_0, elem_map_1);
    const auto vals_test  = face_form_test_->evaluate(elem_map_0, elem_map_1);

    const array<Elem_ *, 2> elems({&elem_0, &elem_1});

    for (int sp_trial = 0; sp_trial < 2; ++sp_trial)
    {
      const auto &elem_tr = *(elems[sp_trial]);
      const auto &elem_dofs_tr = elem_tr.get_local_to_global();
      const auto &vals_tr = vals_trial[sp_trial];
      const auto &n_basis_trial = vals_tr.get_num_functions();

      for (int iFun = 0; iFun < n_basis_trial; ++iFun)
      {
        const auto &dof_i = elem_dofs_tr[iFun];
        const auto &vals_tr_i = vals_tr.get_function_view(iFun);
        auto &val_i = constraints[dof_i];

        for (int sp_test = 0; sp_test < 2; ++sp_test)
        {
          const auto &elem_ts = *(elems[sp_test]);
          const auto &elem_dofs_ts = elem_ts.get_local_to_global();
          const auto &vals_ts = vals_test[sp_test];
          const auto &n_basis_test = vals_ts.get_num_functions();

          for (int jFun = 0; jFun < n_basis_test; ++jFun)
          {
            const auto &dof_j = elem_dofs_ts[jFun];
            const auto &vals_ts_j = vals_ts.get_function_view(jFun);

            Real val_ij = Real(0.0);
            auto vi = vals_tr_i.cbegin();
            auto vj = vals_ts_j.cbegin();

            for (const auto &wm : w_0_at_points)
              val_ij += scalar_product(*vi++, *vj++) * wm;

            val_i[dof_j] += val_ij * integration_constant_;
          } // jFun
        } // sp_test
      } // iFun
    } // sp_trial

  } // end loop elem
};



template <class PhysicalSpace>
void
InterfaceOperator<PhysicalSpace>::
get_element_and_quad_on_face(
  const FacePointsContainer_ &eval_pts_face_ref_domain,
  const shared_ptr<const PhysicalSpace> space,
  const shared_ptr<const PhysFaceSpace_> face_space,
  const Index &face_id,
  const vector<Index> &face_space_dofs,
  const FaceGridMap_ &face_grid_to_space_grid,
  const Index face_element_id,
  Index &element_id,
  vector<Index> &face_elem_dofs,
  PointsContainer_ &current_elem_pts_unit_domain) const
{
  Assert(space != nullptr, ExcNullPtr());
  Assert(face_space != nullptr, ExcNullPtr());

  auto face_elem = face_space->get_element(face_element_id);

  element_id = IteratorState::invalid;
  for (const auto &tmp : face_grid_to_space_grid)
  {
    if (tmp.first->get_flat_index() == face_element_id)
    {
      element_id = tmp.second->get_flat_index();
      break;
    }
  }

  Assert(element_id != IteratorState::invalid,
         ExcMessage("Element not found"));

  auto loc_face_dofs = face_elem.get_local_to_global();
  face_elem_dofs.clear();
  for (const auto &loc_dof : loc_face_dofs)
    face_elem_dofs.push_back(face_space_dofs[loc_dof]);

  const auto grid = space->get_reference_space()->get_grid();
  auto elem = space->get_element(element_id);


  auto cgea_elem         = elem.as_cartesian_grid_element_accessor();

  // Transforming point from face reference domain to space reference domain.
  const auto face_active_directions = UnitElement<dim_>::face_active_directions[face_id];
  const auto face_constant_direction = UnitElement<dim_>::face_constant_direction[face_id];
  const auto face_side = UnitElement<dim_>::face_side[face_id];
  const auto knots = grid->get_knot_coordinates(face_constant_direction);
  const auto const_knot = face_side == 0 ? knots.front() : knots.back();

  PointsContainer_ elem_quad_pts_ref(eval_pts_face_ref_domain.get_num_points());
  auto face_ref_pt_it = eval_pts_face_ref_domain.begin();
  for (auto &rpt : elem_quad_pts_ref)
  {
    const auto &face_ref_pt = *face_ref_pt_it;
    int d = 0;
    for (const auto &dir : face_active_directions)
      rpt[dir] = face_ref_pt[d++];

    rpt[face_constant_direction] = const_knot;
    ++face_ref_pt_it;
  }

  // Transforming point from space reference domain to unit domain.
  current_elem_pts_unit_domain  = cgea_elem.transform_points_reference_to_unit(elem_quad_pts_ref);
}



template <class PhysSpace>
ValueJump<PhysSpace>::
ValueJump(const shared_ptr<const PhysSpace> space_0,
          const shared_ptr<const PhysSpace> space_1)
  :
  Base_(space_0, space_1)
{
  Assert(space_0 != nullptr, ExcNullPtr());
  Assert(space_1 != nullptr, ExcNullPtr());
};


template <class PhysSpace>
auto
ValueJump<PhysSpace>::
create(const shared_ptr<const PhysSpace> space_0,
       const shared_ptr<const PhysSpace> space_1) -> BaseConstPtr_
{
  Assert(space_0 != nullptr, ExcNullPtr());
  Assert(space_1 != nullptr, ExcNullPtr());
  return BaseConstPtr_(new ValueJump(space_0, space_1));
};


template <class PhysSpace>
auto
ValueJump<PhysSpace>::
evaluate(const InterfaceElementMap_& interface_elem_0,
         const InterfaceElementMap_& interface_elem_1) ->
std::array<ValueTable<Value_>, 2>
{
  array<ValueTable<Value_>, 2> values;

  auto &elem_0 = Base_::elem_0_;
  auto &elem_1 = Base_::elem_1_;
  elem_0.move_to(interface_elem_0.get_element_id());
  elem_1.move_to(interface_elem_1.get_element_id());

  values[0] = elem_0->evaluate_basis_values_at_points(interface_elem_0.get_unit_points());

  values[1] = elem_1->evaluate_basis_values_at_points(interface_elem_1.get_unit_points());
  for (auto &v : values[1])
    v *= -1.0;

  return values;
};



template <class PhysSpace>
NormalGradientJump<PhysSpace>::
NormalGradientJump(const shared_ptr<const PhysSpace> space_0,
                   const shared_ptr<const PhysSpace> space_1)
  :
  Base_(space_0, space_1)
{
  Assert(space_0 != nullptr, ExcNullPtr());
  Assert(space_1 != nullptr, ExcNullPtr());
};


template <class PhysSpace>
auto
NormalGradientJump<PhysSpace>::
create(const shared_ptr<const PhysSpace> space_0,
       const shared_ptr<const PhysSpace> space_1) -> BaseConstPtr_
{
  Assert(space_0 != nullptr, ExcNullPtr());
  Assert(space_1 != nullptr, ExcNullPtr());
  return BaseConstPtr_(new NormalGradientJump(space_0, space_1));
};



template <class PhysSpace>
auto
NormalGradientJump<PhysSpace>::
evaluate(const InterfaceElementMap_& interface_elem_0,
         const InterfaceElementMap_& interface_elem_1) ->
std::array<ValueTable<Value_>, 2>
{
  // TODO: include assert for checking normals.

  array<ValueTable<Value_>, 2> values;

  auto &elem_0 = Base_::elem_0_;
  auto &elem_1 = Base_::elem_1_;
  elem_0.move_to(interface_elem_0.get_element_id());
  elem_1.move_to(interface_elem_1.get_element_id());
  const auto &quad_points_0 = interface_elem_0.get_unit_points();
  const auto &quad_points_1 = interface_elem_1.get_unit_points();
  const auto &face_id_0 = interface_elem_0.get_face_id();
  const auto &face_id_1 = interface_elem_1.get_face_id();

  const auto grads_0 = elem_0->evaluate_basis_gradients_at_points(quad_points_0);
  const auto normals_0 = elem_0->get_push_forward_accessor().
  evaluate_normals_at_points(quad_points_0, face_id_0);

  values[0].resize(grads_0.get_num_functions(), grads_0.get_num_points());
  for (int iFun = 0; iFun < grads_0.get_num_functions(); ++iFun)
  {
    auto gi = grads_0.get_function_view(iFun).cbegin();
    auto n = normals_0.cbegin();
    for (auto &v : values[0].get_function_view(iFun))
      v = action(*gi++, *n++);
  }

  const auto grads_1 = elem_1->evaluate_basis_gradients_at_points(quad_points_1);
  const auto normals_1 = elem_1->get_push_forward_accessor().
  evaluate_normals_at_points(quad_points_1, face_id_1);

  values[1].resize(grads_1.get_num_functions(), grads_1.get_num_points());
  for (int iFun = 0; iFun < grads_1.get_num_functions(); ++iFun)
  {
    auto gi = grads_1.get_function_view(iFun).cbegin();
    auto n = normals_1.cbegin();
    for (auto &v : values[1].get_function_view(iFun))
      v = action(*gi++, *n++);
  }

  return values;
};



template <class PhysSpace>
NormalGradientAverage<PhysSpace>::
NormalGradientAverage(const shared_ptr<const PhysSpace> space_0,
                      const shared_ptr<const PhysSpace> space_1,
                      const Real &weight_average)
  :
  Base_(space_0, space_1),
  weight_average_(weight_average)
{
  Assert(space_0 != nullptr, ExcNullPtr());
  Assert(space_1 != nullptr, ExcNullPtr());
  Assert(weight_average_ >= Real(0.0) && weight_average_ <= Real(1.0),
         ExcMessage("weight_average not in range [0, 1]."));
};



template <class PhysSpace>
auto
NormalGradientAverage<PhysSpace>::
create(const shared_ptr<const PhysSpace> space_0,
       const shared_ptr<const PhysSpace> space_1,
       const Real &weight_average) -> BaseConstPtr_
{
  Assert(space_0 != nullptr, ExcNullPtr());
  Assert(space_1 != nullptr, ExcNullPtr());
  return BaseConstPtr_(new NormalGradientAverage(space_0, space_1, weight_average));
};



template <class PhysSpace>
auto
NormalGradientAverage<PhysSpace>::
evaluate(const InterfaceElementMap_& interface_elem_0,
         const InterfaceElementMap_& interface_elem_1) ->
std::array<ValueTable<Value_>, 2>
{
  // TODO: include assert for checking normals.

  array<ValueTable<Value_>, 2> values;

  auto &elem_0 = Base_::elem_0_;
  auto &elem_1 = Base_::elem_1_;
  elem_0.move_to(interface_elem_0.get_element_id());
  elem_1.move_to(interface_elem_1.get_element_id());
  const auto &quad_points_0 = interface_elem_0.get_unit_points();
  const auto &quad_points_1 = interface_elem_1.get_unit_points();
  const auto &face_id_0 = interface_elem_0.get_face_id();
  const auto &face_id_1 = interface_elem_1.get_face_id();

  const Real l1 = weight_average_;
  const Real l2 = Real(1.0) - l1;

  const auto grads_0 = elem_0->evaluate_basis_gradients_at_points(quad_points_0);
  const auto normals_0 = elem_0->get_push_forward_accessor().
  evaluate_normals_at_points(quad_points_0, face_id_0);

  values[0].resize(grads_0.get_num_functions(), grads_0.get_num_points());
  for (int iFun = 0; iFun < grads_0.get_num_functions(); ++iFun)
  {
    auto gi = grads_0.get_function_view(iFun).cbegin();
    auto n = normals_0.cbegin();
    for (auto &v : values[0].get_function_view(iFun))
      v = l1 * action(*gi++, *n++);
  }

  const auto grads_1 = elem_1->evaluate_basis_gradients_at_points(quad_points_1);
  const auto normals_1 = elem_1->get_push_forward_accessor().
  evaluate_normals_at_points(quad_points_1, face_id_1);

  values[1].resize(grads_1.get_num_functions(), grads_1.get_num_points());
  for (int iFun = 0; iFun < grads_1.get_num_functions(); ++iFun)
  {
    auto gi = grads_1.get_function_view(iFun).cbegin();
    auto n = normals_1.cbegin();
    for (auto &v : values[1].get_function_view(iFun))
      v = l2 * action(*gi++, *n++);
  }

  return values;
};



template <class PhysSpace>
LaplacianAverage<PhysSpace>::
LaplacianAverage(const shared_ptr<const PhysSpace> space_0,
                 const shared_ptr<const PhysSpace> space_1,
                 const Real &weight_average)
  :
  Base_(space_0, space_1),
  weight_average_(weight_average)
{
  Assert(space_0 != nullptr, ExcNullPtr());
  Assert(space_1 != nullptr, ExcNullPtr());
  Assert(weight_average_ >= Real(0.0) && weight_average_ <= Real(1.0),
         ExcMessage("weight_average not in range [0, 1]."));
};



template <class PhysSpace>
auto
LaplacianAverage<PhysSpace>::
create(const shared_ptr<const PhysSpace> space_0,
       const shared_ptr<const PhysSpace> space_1,
       const Real &weight_average) -> BaseConstPtr_
{
  Assert(space_0 != nullptr, ExcNullPtr());
  Assert(space_1 != nullptr, ExcNullPtr());
  return BaseConstPtr_(new LaplacianAverage(space_0, space_1, weight_average));
};



template <class PhysSpace>
auto
LaplacianAverage<PhysSpace>::
evaluate(const InterfaceElementMap_& interface_elem_0,
         const InterfaceElementMap_& interface_elem_1) ->
std::array<ValueTable<Value_>, 2>
{
  // TODO: include assert for checking normals.

  array<ValueTable<Value_>, 2> values;

  auto &elem_0 = Base_::elem_0_;
  auto &elem_1 = Base_::elem_1_;
  elem_0.move_to(interface_elem_0.get_element_id());
  elem_1.move_to(interface_elem_1.get_element_id());
  const auto &quad_points_0 = interface_elem_0.get_unit_points();
  const auto &quad_points_1 = interface_elem_1.get_unit_points();

  const Real l1 = weight_average_;
  const Real l2 = Real(1.0) - l1;

  const auto hessians_0 = elem_0->evaluate_basis_hessians_at_points(quad_points_0);

  values[0].resize(hessians_0.get_num_functions(), hessians_0.get_num_points());

  auto h0 = hessians_0.cbegin();
  for (auto &v : values[0])
    v = l1 * trace(*h0++);

  const auto hessians_1 = elem_1->evaluate_basis_hessians_at_points(quad_points_1);

  values[1].resize(hessians_1.get_num_functions(), hessians_1.get_num_points());

  auto h1 = hessians_1.cbegin();
  for (auto &v : values[1])
    v = l2 * trace(*h1++);

  return values;
};



template <class PhysSpace>
InterfaceOperator<PhysSpace>::
InterfaceElementMap::
InterfaceElementMap(const Index &elem_id,
                    const Index &interface_elem_id,
                    const Index &face_id,
                    const PointsContainer_ &unit_points,
                    const FacePointsContainer_ &face_unit_points)
  :
  elem_id_(elem_id),
  interface_elem_id_(interface_elem_id),
  face_id_(face_id),
  unit_points_(unit_points),
  face_unit_points_(face_unit_points)
{};



template <class PhysSpace>
const Index &
InterfaceOperator<PhysSpace>::
InterfaceElementMap::
get_element_id() const
{
  return elem_id_;
};



template <class PhysSpace>
const Index &
InterfaceOperator<PhysSpace>::
InterfaceElementMap::
get_interface_element_id() const
{
  return interface_elem_id_;
};



template <class PhysSpace>
const Index &
InterfaceOperator<PhysSpace>::
InterfaceElementMap::
get_face_id() const
{
  return face_id_;
};



template <class PhysSpace>
auto
InterfaceOperator<PhysSpace>::
InterfaceElementMap::
get_unit_points() const ->
const PointsContainer_ &
{
  return unit_points_;
};



template <class PhysSpace>
auto
InterfaceOperator<PhysSpace>::
InterfaceElementMap::
get_face_unit_points() const ->
const FacePointsContainer_ &
{
  return face_unit_points_;
};



template <class PhysSpace>
InterfaceOperator<PhysSpace>::
Form::
Form(const shared_ptr<const PhysSpace> space_0,
     const shared_ptr<const PhysSpace> space_1)
  :
  space_0_(space_0),
  space_1_(space_1),
  elem_0_(space_0->begin()),
  elem_1_(space_1->begin())
{
  Assert(space_0 != nullptr, ExcNullPtr());
  Assert(space_1 != nullptr, ExcNullPtr());
};

#endif

IGA_NAMESPACE_CLOSE

//#include <igatools/operators/interface_operator.inst>
