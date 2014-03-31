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

#include <igatools/geometry/ig_mapping.h>
#include <igatools/base/exceptions.h>
#include <igatools/basis_functions/space_tools.h>
#include <igatools/basis_functions/physical_space_element_accessor.h>

using std::vector;
using std::array;
using std::shared_ptr;
using std::endl;

IGA_NAMESPACE_OPEN

namespace
{
template<class RefSpace>
StaticMultiArray<DynamicMultiArray<Real,RefSpace::dim>,RefSpace::range,1>
get_weights_from_ref_space(const RefSpace &ref_space,
                           EnableIf<RefSpace::has_weights> *hw = 0)
{
    //in the case of NURBSSpace get the weights used in the space
    return ref_space.get_weights();
}



template<class RefSpace>
StaticMultiArray<DynamicMultiArray<Real,RefSpace::dim>,RefSpace::range,1>
get_weights_from_ref_space(const RefSpace &ref_space,
                           EnableIf<!RefSpace::has_weights> *hw = 0)
{
    //in the case of BSplineSpace do nothing
    StaticMultiArray<DynamicMultiArray<Real,RefSpace::dim>,RefSpace::range,1> weights;
    return weights;
}
};



template<class RefSpace>
IgMapping<RefSpace>::
IgMapping(const std::shared_ptr<RefSpace> space,
          const std::vector<Real> &control_points)
    :
    base_t::Mapping(space->get_grid()),
    ref_space_(space),
    control_points_(control_points),
    element_(space->begin())
{
    Assert(space != nullptr, ExcNullPtr());
    Assert(space->get_num_basis() == control_points_.size(),
           ExcDimensionMismatch(space->get_num_basis(), control_points_.size()));
    Assert(RefSpace::rank == 1, ExcDimensionMismatch(RefSpace::rank,1));


    //----------------------------------
    // if RefSpace is NURBSSpace
    // save the weights in order to be used in the h-refinement algorithm
    // (the different possibilities for RefSpace are handled by specialization of the
    // function get_weights_from_ref_space() in the anonymous namespace above).
    weights_pre_refinement_ = get_weights_from_ref_space(*ref_space_);
    /*
        LogStream out ;
        for (Index comp_id = 0 ; comp_id < space_dim ; ++comp_id)
            out << "weights_pre_refinement_["<<comp_id<<"]= " << weights_pre_refinement_[comp_id] << std::endl;
    //*/
    //----------------------------------


    //----------------------------------
    // copy the knots (with repetitions) that defines the RefSpace before any refinement
    knots_with_repetitions_pre_refinement_ = space->get_knots_with_repetitions();
    //----------------------------------



    //----------------------------------
    // copy the control mesh before any refinement
    const auto &index_space = ref_space_->get_index_space();

    for (int comp_id = 0 ; comp_id < space_dim ; ++comp_id)
    {
        const auto &index_space_comp = index_space(comp_id);
        auto &ctrl_mesh_comp = ctrl_mesh_(comp_id);

        ctrl_mesh_comp.resize(index_space_comp.tensor_size());

        const Size n_dofs_comp = ref_space_->get_component_num_basis(comp_id);
//        out << "n_dofs_comp["<<comp_id<<"]= " << n_dofs_comp << endl ;

        const auto &weights_pre_refinement_comp = weights_pre_refinement_(comp_id);

        for (Index loc_id = 0 ; loc_id < n_dofs_comp ; ++loc_id)
        {
            const Index glob_id = index_space_comp(loc_id);
            if (RefSpace::has_weights)
            {
                // if NURBS, transform the control points from euclidean to projective coordinates
                const Real w = weights_pre_refinement_comp(loc_id);

                ctrl_mesh_comp(loc_id) = w * control_points_[glob_id];
            }
            else
                ctrl_mesh_comp(loc_id) = control_points_[glob_id];

        }
//        out << "ctrl_mesh_["<<comp_id<<"]= " << ctrl_mesh_[comp_id] << endl;
    }

    // create a signal and a connection for the control-mesh refinement
    this->connect_refinement_h_function(
        std::bind(
            &IgMapping<RefSpace>::refine_h_control_mesh,
            this,
            std::placeholders::_1,std::placeholders::_2));
    //----------------------------------
}



template<class RefSpace>
void
IgMapping<RefSpace>::
init_element(const ValueFlags flag,
             const Quadrature<dim> &quad)
{
    ValueFlags ref_space_flag = ValueFlags::none;

    if (contains(flag,ValueFlags::point) || contains(flag,ValueFlags::map_value))
    {
        ref_space_flag |= ValueFlags::value;
    }

    if (contains(flag,ValueFlags::face_point) || contains(flag,ValueFlags::map_face_value))
    {
        ref_space_flag |= ValueFlags::face_value;
    }

    if (contains(flag,ValueFlags::map_gradient))
    {
        ref_space_flag |= ValueFlags::gradient;
    }

    if (contains(flag,ValueFlags::map_face_gradient))
    {
        ref_space_flag |= ValueFlags::face_gradient;
    }

    if (contains(flag,ValueFlags::map_hessian))
    {
        ref_space_flag |= ValueFlags::hessian;
    }

    if (contains(flag,ValueFlags::map_face_hessian))
    {
        ref_space_flag |= ValueFlags::face_hessian;
    }

    element_->init_values(ref_space_flag, quad);
}



template<class RefSpace>
void IgMapping<RefSpace>::
set_element(const CartesianGridElementAccessor<dim> &elem)
{
    element_->reset_flat_tensor_indices(elem.get_flat_index());
    element_->fill_values();
}



template<class RefSpace>
void IgMapping<RefSpace>::
set_face_element(const Index face_id, const CartesianGridElementAccessor<dim> &elem)
{
    Assert(face_id < UnitElement<dim>::faces_per_element && face_id >= 0,
           ExcIndexRange(face_id,0,UnitElement<dim>::faces_per_element));
    element_->reset_flat_tensor_indices(elem.get_flat_index());
    element_->fill_face_values(face_id);
}



template<class RefSpace>
auto
IgMapping<RefSpace>::create(
    const std::shared_ptr<RefSpace> space,
    const std::vector<Real> &control_points) -> shared_ptr<base_t>
{
    return (shared_ptr<Mapping<dim,codim>>(
        new IgMapping<RefSpace>(space,control_points)));
}



template<class RefSpace>
auto
IgMapping<RefSpace>::
clone() const -> shared_ptr<base_t>
{
    auto map = shared_ptr<IgMapping<RefSpace>>(new IgMapping<RefSpace>(*this));
    map->element_->reset_global_cache() ;
    return map;
}



template<class RefSpace>
void
IgMapping<RefSpace>::
evaluate(vector<ValueType> &values) const
{
    const auto &local_to_global = element_->get_local_to_global();

    vector<Real> ctrl_pts_element;

    for (const auto &local_id : local_to_global)
        ctrl_pts_element.emplace_back(control_points_[local_id]);

    values = element_->evaluate_field(ctrl_pts_element);
}



template<class RefSpace>
void
IgMapping<RefSpace>::
evaluate_gradients(std::vector<GradientType> &gradients) const
{
    const auto &local_to_global = element_->get_local_to_global();

    vector<Real> ctrl_pts_element;

    for (const auto &local_id : local_to_global)
        ctrl_pts_element.emplace_back(control_points_[local_id]);

    gradients = element_->evaluate_field_gradients(ctrl_pts_element);

}


template<class RefSpace>
void
IgMapping<RefSpace>::
evaluate_hessians(std::vector<HessianType> &hessians) const
{
    const auto &local_to_global = element_->get_local_to_global();

    vector<Real> ctrl_pts_element;

    for (const auto &local_id : local_to_global)
        ctrl_pts_element.emplace_back(control_points_[local_id]);

    hessians = element_->evaluate_field_hessians(ctrl_pts_element);
}



template<class RefSpace>
void
IgMapping<RefSpace>::
evaluate_face(const Index face_id, vector<ValueType> &values) const
{
    Assert(face_id < UnitElement<dim>::faces_per_element && face_id >= 0,
           ExcIndexRange(face_id,0,UnitElement<dim>::faces_per_element));

    const auto &local_to_global = element_->get_local_to_global();

    vector<Real> ctrl_pts_element;

    for (const auto &local_id : local_to_global)
        ctrl_pts_element.emplace_back(control_points_[local_id]);

    values = element_->evaluate_field(ctrl_pts_element,FaceTopology<dim>(face_id));
}



template<class RefSpace>
void
IgMapping<RefSpace>::
evaluate_face_gradients(const Index face_id, std::vector<GradientType> &gradients) const
{
    Assert(face_id < UnitElement<dim>::faces_per_element && face_id >= 0,
           ExcIndexRange(face_id,0,UnitElement<dim>::faces_per_element));

    const auto &local_to_global = element_->get_local_to_global();

    vector<Real> ctrl_pts_element;

    for (const auto &local_id : local_to_global)
        ctrl_pts_element.emplace_back(control_points_[local_id]);

    gradients = element_->evaluate_field_gradients(ctrl_pts_element,FaceTopology<dim>(face_id));
}


template<class RefSpace>
void
IgMapping<RefSpace>::
evaluate_face_hessians(const Index face_id, std::vector<HessianType> &hessians) const
{
    Assert(face_id < UnitElement<dim>::faces_per_element && face_id >= 0,
           ExcIndexRange(face_id,0,UnitElement<dim>::faces_per_element));

    const auto &local_to_global = element_->get_local_to_global();

    vector<Real> ctrl_pts_element;

    for (const auto &local_id : local_to_global)
        ctrl_pts_element.emplace_back(control_points_[local_id]);

    hessians = element_->evaluate_field_hessians(ctrl_pts_element,FaceTopology<dim>(face_id));
}




template<class RefSpace>
void
IgMapping<RefSpace>::
set_control_points(const std::vector<Real> &control_points)
{
    Assert(control_points_.size() == control_points.size(),
           ExcDimensionMismatch(control_points_.size(), control_points.size()));

    control_points_ = control_points;
}




template<class RefSpace>
void
IgMapping<RefSpace>::
refine_h_control_mesh(
    const std::array<bool,dim> &refinement_directions,
    const typename base_t::GridType &grid_old)
{
    for (int direction_id = 0 ; direction_id < dim ; ++direction_id)
    {
        if (refinement_directions[direction_id])
        {

            // knots in the refined grid along the selected direction
            vector<Real> knots_new = ref_space_->get_grid()->get_knot_coordinates(direction_id);

            // knots in the original (unrefined) grid along the selected direction
            vector<Real> knots_old = grid_old.get_knot_coordinates(direction_id);

            vector<Real> knots_added(knots_new.size());

            // find the knots in the refined grid that are not present in the old grid
            auto it = std::set_difference(
                          knots_new.begin(),knots_new.end(),
                          knots_old.begin(),knots_old.end(),
                          knots_added.begin());

            knots_added.resize(it-knots_added.begin());


            for (int comp_id = 0 ; comp_id < space_dim ; ++comp_id)
            {

                const int p = ref_space_->get_degree()(comp_id)[direction_id];

                const auto &U = knots_with_repetitions_pre_refinement_(comp_id).get_data_direction(direction_id);
                const auto &X = knots_added;
                const auto &Ubar = ref_space_->get_knots_with_repetitions()(comp_id).get_data_direction(direction_id);

                const int m = U.size()-1;
                const int r = X.size()-1;
                const int a = space_tools::find_span(p,X[0],U);
                const int b = space_tools::find_span(p,X[r],U)+1;

                const int n = m-p-1;


                const auto &Pw = ctrl_mesh_(comp_id);
                const auto old_sizes = Pw.tensor_size();
                Assert(old_sizes[direction_id] == n+1,
                       ExcDimensionMismatch(old_sizes[direction_id],n+1));


                auto new_sizes = old_sizes;
                new_sizes[direction_id] += r+1; // r+1 new weights in the refinement direction
                Assert(new_sizes[direction_id] ==
                       ref_space_->get_component_dir_num_basis(comp_id,direction_id),
                       ExcDimensionMismatch(new_sizes[direction_id],
                                            ref_space_->get_component_dir_num_basis(comp_id,direction_id)));

                DynamicMultiArray<Real,dim> Qw(new_sizes);

                for (Index j = 0 ; j <= a-p ; ++j)
                {
                    Qw.copy_slice(direction_id,j,
                                  Pw.get_slice(direction_id,j));
                }

                for (Index j = b-1 ; j <= n ; ++j)
                {
                    Qw.copy_slice(direction_id,j+r+1,
                                  Pw.get_slice(direction_id,j));
                }

                Index i = b + p - 1;
                Index k = b + p + r;
                for (Index j = r ; j >= 0 ; --j)
                {
                    while (X[j] <= U[i] && i > a)
                    {
                        Qw.copy_slice(direction_id,k-p-1,Pw.get_slice(direction_id,i-p-1));
                        k = k-1;
                        i = i-1;
                    }
                    Qw.copy_slice(direction_id,k-p-1,
                                  Qw.get_slice(direction_id,k-p));

                    for (Index l = 1 ; l <= p ; ++l)
                    {
                        Index ind = k-p+l;

                        Real alfa = Ubar[k+l] - X[j];
                        if (fabs(alfa) == 0.0)
                        {
                            Qw.copy_slice(direction_id,ind-1,Qw.get_slice(direction_id,ind));
                        }
                        else
                        {
                            alfa = alfa / (Ubar[k+l] - U[i-p+l]);

                            Qw.copy_slice(direction_id,ind-1,
                                          alfa  * Qw.get_slice(direction_id,ind-1) +
                                          (1.0-alfa) * Qw.get_slice(direction_id,ind));
                        }
                    } // end loop l
                    k = k-1;

                } // end loop j

                ctrl_mesh_(comp_id) = Qw;
                //*/
            } // end loop comp_id
        } // end if (refinement_directions[direction_id])

    } // end loop direction_id



    //----------------------------------
    // copy the control mesh after the refinement
    control_points_.resize(ref_space_->get_num_basis());
//    control_points_ = Vector(dof_tools::get_dofs(*ref_space_));

//    const auto &index_space = ref_space_->get_index_space();

    const auto weights_after_refinement = get_weights_from_ref_space(*ref_space_);

    Index ctrl_pt_id = 0;
    for (int comp_id = 0 ; comp_id < space_dim ; ++comp_id)
    {
//        const auto &index_space_comp = index_space(comp_id);
        const auto &ctrl_mesh_comp = ctrl_mesh_(comp_id);
        const auto &weights_after_refinement_comp = weights_after_refinement(comp_id);

        const Size n_dofs_comp = ref_space_->get_component_num_basis(comp_id);
        for (Index loc_id = 0 ; loc_id < n_dofs_comp ; ++loc_id, ++ctrl_pt_id)
        {
//            const Index glob_id = index_space_comp(loc_id);
            if (RefSpace::has_weights)
            {
                // if NURBS, transform the control points from  projective to euclidean coordinates
                const Real w = weights_after_refinement_comp(loc_id);

                control_points_[ctrl_pt_id] = ctrl_mesh_comp(loc_id) / w ;
            }
            else
                control_points_[ctrl_pt_id] = ctrl_mesh_comp(loc_id);
        }
    }
    //----------------------------------
}





template<class RefSpace>
void
IgMapping<RefSpace>::
print_info(LogStream &out) const
{
    out << "Type = IgMapping" << endl;

    out.push("\t");
    out << "Reference space info:" << endl;
    ref_space_->print_info(out);
    out << endl;

    out.push("\t");
    out << "Control points info (projective coordinates):" << endl;

    out.push("\t");
    for (Index comp_id = 0 ; comp_id < space_dim ; ++comp_id)
        out << "Control mesh["<<comp_id<<"] = " << ctrl_mesh_(comp_id) << endl;
    out << endl;
    out.pop();


    out << "Control points info (euclidean coordinates): [ ";
    for (const auto &ctrl_pt : control_points_)
        out << ctrl_pt << " ";
    out << endl;
//    control_points_.print(out);

    out.pop();
}

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/ig_mapping.inst>
