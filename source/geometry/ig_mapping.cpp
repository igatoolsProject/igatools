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
#include <igatools/utils/vector_tools.h>


using std::array;
using std::shared_ptr;
using std::make_shared;
using std::endl;
using std::const_pointer_cast;

IGA_NAMESPACE_OPEN

namespace
{
//template<class RefSpace>
//typename NURBSSpace<RefSpace::dim,RefSpace::range,RefSpace::rank>::WeightsTable
//get_weights_from_ref_space(const RefSpace &ref_space,
//                           EnableIf<RefSpace::has_weights> *hw = 0)
//{
//    //in the case of NURBSSpace get the weights used in the space
//    return ref_space.get_weights();
//}
//
//
//
//template<class RefSpace>
//typename NURBSSpace<RefSpace::dim,RefSpace::range,RefSpace::rank>::WeightsTable
//get_weights_from_ref_space(const RefSpace &ref_space,
//                           EnableIf<!RefSpace::has_weights> *hw = 0)
//{
//    //in the case of BSplineSpace do nothing (it should returns all weights equal to 1.0)
//    typename NURBSSpace<RefSpace::dim,RefSpace::range,RefSpace::rank>::WeightsTable
//    weights(ref_space.get_components_map());
//
//
//    const auto &basis_tensor_size_table = ref_space.get_num_basis_table();
//
//    for (Index comp_id : weights.get_active_components_id())
//    {
//        auto &weights_component = weights[comp_id];
//
//        weights_component.resize(basis_tensor_size_table[comp_id]);
//
//        for (auto &w : weights_component)
//            w = 1.0;
//    }
//
//
//    return weights;
//}


//template<int dim, int range, int rank>
//const BSplineSpace<dim,range,rank> &
//get_bspline_space(const NURBSSpace<dim,range,rank> &nurbs_space)
//{
//    return *nurbs_space.get_spline_space();
//}

template<int dim, int range, int rank>
const BSplineSpace<dim,range,rank> &
get_bspline_space(const BSplineSpace<dim,range,rank> &bspline_space)
{
    return bspline_space;
}

};



template<class RefSpace>
IgMapping<RefSpace>::
IgMapping(const std::shared_ptr<RefSpace> space,
          const vector<Real> &control_points)
    :
    base_t::SplineMapping(space->get_grid()),
    data_(shared_ptr<IgMappingData>(new IgMappingData)),
    cache_(space->begin())
{
    Assert(space != nullptr, ExcNullPtr());
    data_->ref_space_ = space;

    data_->control_points_ = control_points;
    Assert(space->get_num_basis() == data_->control_points_.size(),
           ExcDimensionMismatch(space->get_num_basis(), data_->control_points_.size()));

    Assert(RefSpace::rank == 1, ExcDimensionMismatch(RefSpace::rank,1));

#if 0
    //----------------------------------
    // if RefSpace is NURBSSpace
    // save the weights in order to be used in the h-refinement algorithm
    // (the different possibilities for RefSpace are handled by specialization of the
    // function get_weights_from_ref_space() in the anonymous namespace above).
    data_->weights_pre_refinement_ = get_weights_from_ref_space(*(data_->ref_space_));
    //----------------------------------




    //----------------------------------
    // copy the control mesh before any refinement
    using bspline_space_t = BSplineSpace<RefSpace::dim,RefSpace::range,RefSpace::rank>;
    const bspline_space_t &bspline_space = get_bspline_space(*space);
    const auto &num_basis_table = bspline_space.get_num_basis_table();


    Index ctrl_pt_fid = 0;
    for (int comp_id  = 0 ; comp_id < space_dim ; ++comp_id)
    {
        const TensorSize<dim> &num_basis_comp = num_basis_table[comp_id];

        auto &ctrl_mesh_comp = data_->ctrl_mesh_[comp_id];

        ctrl_mesh_comp.resize(num_basis_comp);

        const Size n_dofs_comp = data_->ref_space_->get_num_basis(comp_id);

        const auto &weights_pre_refinement_comp = data_->weights_pre_refinement_[comp_id];

        for (Index loc_id = 0 ; loc_id < n_dofs_comp ; ++loc_id)
        {
            if (RefSpace::has_weights)
            {
                // If NURBS, transform the control points from euclidean to
                // projective coordinates.
                const Real &w = weights_pre_refinement_comp[loc_id];

                ctrl_mesh_comp[loc_id] = w * data_->control_points_[ctrl_pt_fid];
            }
            else
                ctrl_mesh_comp[loc_id] = data_->control_points_[ctrl_pt_fid];

            ++ctrl_pt_fid;
        }
    }

    this->connect_refinement_h_function(
        std::bind(
            &IgMapping<RefSpace>::refine_h_control_mesh,
            this,
            std::placeholders::_1,std::placeholders::_2));
#endif
}



template<class RefSpace>
IgMapping<RefSpace>::
IgMapping(const self_t &map)
    :
    base_t::SplineMapping(map),
    data_(shared_ptr<IgMappingData>(new IgMappingData(*map.data_))),
    cache_(map.cache_)
{}

template<class RefSpace>
IgMapping<RefSpace>::
IgMapping(const std::shared_ptr<IgMappingData> mapping_data)
    :
    base_t::SplineMapping(mapping_data->ref_space_->get_grid()),
    data_(mapping_data),
    cache_(mapping_data->ref_space_->begin())
{
    Assert(mapping_data != nullptr,ExcNullPtr());
}



template<class RefSpace>
auto
IgMapping<RefSpace>::
get_data() const -> shared_ptr<IgMappingData>
{
    return data_;
}



template<class RefSpace>
void
IgMapping<RefSpace>::
init_element(const ValueFlags flag,
             const Quadrature<dim> &quad)  const
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
    // TODO (pauletti, Sep 12, 2014): update to new cache
    //cache_->init_cache(ref_space_flag, quad);
}



template<class RefSpace>
void IgMapping<RefSpace>::
set_element(const GridIterator &elem) const
{
    cache_->move_to(elem.get_flat_index());
    // TODO (pauletti, Sep 12, 2014): update to new cache
    //cache_->fill_cache();
}



template<class RefSpace>
void IgMapping<RefSpace>::
set_face_element(const Index face_id, const GridIterator &elem) const
{
    Assert(face_id < UnitElement<dim>::faces_per_element && face_id >= 0,
           ExcIndexRange(face_id,0,UnitElement<dim>::faces_per_element));
    cache_->move_to(elem.get_flat_index());
    // TODO (pauletti, Sep 12, 2014): fix next line
    Assert(true, ExcMessage(" fix next line "));
    //cache_->fill_face_cache(face_id);
}



template<class RefSpace>
auto
IgMapping<RefSpace>::create(
    const std::shared_ptr<RefSpace> space,
    const vector<Real> &control_points) -> shared_ptr<Mapping<dim,codim>>
{
    return (shared_ptr<Mapping<dim,codim>>(
        new IgMapping<RefSpace>(space,control_points)));
}



template<class RefSpace>
vector<Real>
IgMapping<RefSpace>::
get_control_points_elem() const
{
    Assert(data_ != nullptr, ExcNullPtr());
    const auto &local_to_patch = cache_->get_local_to_patch();

    vector<Real> ctrl_pts_element;

    for (const auto &local_id : local_to_patch)
        ctrl_pts_element.emplace_back(data_->control_points_[local_id]);

    return ctrl_pts_element;
}

template<class RefSpace>
void
IgMapping<RefSpace>::
evaluate(ValueVector<Value> &values) const
{
    values = cache_->evaluate_field(this->get_control_points_elem());
}



template<class RefSpace>
void
IgMapping<RefSpace>::
evaluate_gradients(ValueVector<Gradient> &gradients) const
{
    gradients = cache_->evaluate_field_gradients(this->get_control_points_elem());
}


template<class RefSpace>
void
IgMapping<RefSpace>::
evaluate_hessians(ValueVector<Hessian> &hessians) const
{
    hessians = cache_->evaluate_field_hessians(this->get_control_points_elem());
}



template<class RefSpace>
void
IgMapping<RefSpace>::
evaluate_face(const Index face_id, ValueVector<Value> &values) const
{
    values = cache_->evaluate_field(this->get_control_points_elem(),FaceTopology<dim>(face_id));
}



template<class RefSpace>
void
IgMapping<RefSpace>::
evaluate_face_gradients(const Index face_id, ValueVector<Gradient> &gradients) const
{
    gradients = cache_->evaluate_field_gradients(this->get_control_points_elem(),FaceTopology<dim>(face_id));
}


template<class RefSpace>
void
IgMapping<RefSpace>::
evaluate_face_hessians(const Index face_id, ValueVector<Hessian> &hessians) const
{
    hessians = cache_->evaluate_field_hessians(this->get_control_points_elem(),FaceTopology<dim>(face_id));
}



// TODO (pauletti, Aug 7, 2014): evaluate_at_points, evaluate_gradients_at_points, etc
// should eventually be a template function template<order> eval_derivative_at_point
template<class RefSpace>
void
IgMapping<RefSpace>::
evaluate_at_points(const ValueVector<Point> &points, ValueVector<Value> &values) const
{
    Assert(points.size() > 0, ExcEmptyObject());
    Assert(values.size() == points.size(),
           ExcDimensionMismatch(values.size(),points.size()));

    auto elem_list = this->get_grid()->find_elements_of_points(points);
    auto &elem = cache_;

    for (auto p : elem_list)
    {
        elem->move_to(p.first->get_flat_index());
        ValueVector<Point> pts(p.second.size());
        for (int j=0; j<p.second.size(); ++j)
            pts[j] = points[p.second[j]];

        const auto points_unit_element =
            elem->transform_points_reference_to_unit(pts);
        const auto elem_val =
            elem->evaluate_field_values_at_points(
                get_control_points_elem(), points_unit_element);

        for (int j=0; j<p.second.size(); ++j)
            values[p.second[j]] = elem_val[j];
    }
}



template<class RefSpace>
void
IgMapping<RefSpace>::
evaluate_gradients_at_points(const ValueVector<Point> &points,
                             ValueVector<Gradient> &gradients) const
{
    Assert(points.size() > 0, ExcEmptyObject());
    Assert(gradients.size() == points.size(),
           ExcDimensionMismatch(gradients.size(),points.size()));

    auto elem_list = this->get_grid()->find_elements_of_points(points);
    auto &elem = cache_;

    for (auto p : elem_list)
    {
        elem->move_to(p.first->get_flat_index());
        ValueVector<Point> pts(p.second.size());
        for (int j=0; j<p.second.size(); ++j)
            pts[j] = points[p.second[j]];

        const auto points_unit_element =
            elem->transform_points_reference_to_unit(pts);
        const auto elem_val =
            elem->evaluate_field_gradients_at_points(
                get_control_points_elem(), points_unit_element);

        for (int j=0; j<p.second.size(); ++j)
            gradients[p.second[j]] = elem_val[j];
    }

}



template<class RefSpace>
void
IgMapping<RefSpace>::
evaluate_hessians_at_points(const ValueVector<Point> &points,
                            ValueVector<Hessian> &hessians) const
{
    Assert(points.size() > 0, ExcEmptyObject());
    Assert(hessians.size() == points.size(),
           ExcDimensionMismatch(hessians.size(),points.size()));

    auto elem_list = this->get_grid()->find_elements_of_points(points);
    auto &elem = cache_;

    for (auto p : elem_list)
    {
        elem->move_to(p.first->get_flat_index());
        ValueVector<Point> pts(p.second.size());
        for (int j=0; j<p.second.size(); ++j)
            pts[j] = points[p.second[j]];

        const auto points_unit_element =
            elem->transform_points_reference_to_unit(pts);
        const auto elem_val =
            elem->evaluate_field_hessians_at_points(
                get_control_points_elem(), points_unit_element);

        for (int j=0; j<p.second.size(); ++j)
            hessians[p.second[j]] = elem_val[j];
    }

}


template<class RefSpace>
void
IgMapping<RefSpace>::
set_control_points(const vector<Real> &control_points)
{
#if 0
    Assert(data_->control_points_.size() == control_points.size(),
           ExcDimensionMismatch(data_->control_points_.size(), control_points.size()));

    // Updating of the euclidean coordinates of the control_points.
    data_->control_points_ = control_points;


    Assert(data_->ref_space_ != nullptr, ExcNullPtr());
    //const auto weights = get_weights_from_ref_space(*(data_->ref_space_));


    // Updating of the euclidean (in case of BSpline) or projective (in case of NURBS)
    // coordinates of the control_points .

    Index ctrl_pt_fid = 0;
    for (int comp_id = 0 ; comp_id < space_dim ; ++comp_id)
    {
//        const TensorSize<dim> &num_basis_comp = num_basis_table[comp_id];

        auto &ctrl_mesh_comp = data_->ctrl_mesh_[comp_id];

        const Size n_dofs_comp = data_->ref_space_->get_num_basis(comp_id);

      // const auto &weights_after_refinement_comp = weights[comp_id];

        for (Index loc_id = 0 ; loc_id < n_dofs_comp ; ++loc_id)
        {
            if (RefSpace::has_weights)
            {
                // If NURBS, transform the control points from euclidean to
                // projective coordinates.
                const Real &w = weights_after_refinement_comp[loc_id];

                ctrl_mesh_comp[loc_id] = w * data_->control_points_[ctrl_pt_fid];
            }
            else
                ctrl_mesh_comp[loc_id] = data_->control_points_[ctrl_pt_fid];

            ++ctrl_pt_fid;

        } // end loop loc_id
    } // end loop comp_id
#endif
}

template<class RefSpace>
vector<Real>
IgMapping<RefSpace>::
get_control_points() const
{
    Assert(data_ != nullptr, ExcNullPtr());

    return data_->control_points_;
}



template<class RefSpace>
void
IgMapping<RefSpace>::
refine_h_control_mesh(
    const std::array<bool,dim> &refinement_directions,
    const typename base_t::GridType &grid_old1)
{
#if 0
    auto grid = this->get_grid();
    auto grid_old = this->get_grid()->get_grid_pre_refinement();

    auto ref_space = data_->ref_space_;

    using bspline_space_t = BSplineSpace<RefSpace::dim,RefSpace::range,RefSpace::rank>;
    const bspline_space_t &bspline_space = get_bspline_space(*ref_space);

    auto knots_with_repetitions_pre_refinement = bspline_space.get_spline_space_previous_refinement()
                                                 ->compute_knots_with_repetition(
                                                     bspline_space.get_end_behaviour());
    auto knots_with_repetitions = bspline_space.compute_knots_with_repetition(
                                      bspline_space.get_end_behaviour());

    for (int direction_id = 0 ; direction_id < dim ; ++direction_id)
    {
        if (refinement_directions[direction_id])
        {
            // knots in the refined grid along the selected direction
            vector<Real> knots_new = grid->get_knot_coordinates(direction_id);

            // knots in the original (unrefined) grid along the selected direction
            vector<Real> knots_old = grid_old->get_knot_coordinates(direction_id);

            vector<Real> knots_added(knots_new.size());

            // find the knots in the refined grid that are not present in the old grid
            auto it = std::set_difference(
                          knots_new.begin(),knots_new.end(),
                          knots_old.begin(),knots_old.end(),
                          knots_added.begin());

            knots_added.resize(it-knots_added.begin());


            for (int comp_id = 0 ; comp_id < space_dim ; ++comp_id)
            {
                const int p = ref_space->get_degree()[comp_id][direction_id];
                const auto &U = knots_with_repetitions_pre_refinement[comp_id].get_data_direction(direction_id);
                const auto &X = knots_added;
                const auto &Ubar = knots_with_repetitions[comp_id].get_data_direction(direction_id);

                const int m = U.size()-1;
                const int r = X.size()-1;
                const int a = space_tools::find_span(p,X[0],U);
                const int b = space_tools::find_span(p,X[r],U)+1;

                const int n = m-p-1;

                const auto &Pw = data_->ctrl_mesh_[comp_id];
                const auto old_sizes = Pw.tensor_size();
                Assert(old_sizes[direction_id] == n+1,
                       ExcDimensionMismatch(old_sizes[direction_id],n+1));


                auto new_sizes = old_sizes;
                new_sizes[direction_id] += r+1; // r+1 new weights in the refinement direction
                Assert(new_sizes[direction_id] ==
                       data_->ref_space_->get_num_basis(comp_id,direction_id),
                       ExcDimensionMismatch(new_sizes[direction_id],
                                            data_->ref_space_->get_num_basis(comp_id,direction_id)));

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

                data_->ctrl_mesh_[comp_id] = Qw;
                //*/
            } // end loop comp_id
        } // end if (refinement_directions[direction_id])

    } // end loop direction_id



    //----------------------------------
    // copy the control mesh after the refinement
    data_->control_points_.resize(data_->ref_space_->get_num_basis());

    const auto weights_after_refinement = get_weights_from_ref_space(*(data_->ref_space_));

    Index ctrl_pt_id = 0;
    for (int comp_id = 0 ; comp_id < space_dim ; ++comp_id)
    {
        const auto &ctrl_mesh_comp = data_->ctrl_mesh_[comp_id];
        const auto &weights_after_refinement_comp = weights_after_refinement[comp_id];

        const Size n_dofs_comp = data_->ref_space_->get_num_basis(comp_id);
        for (Index loc_id = 0 ; loc_id < n_dofs_comp ; ++loc_id, ++ctrl_pt_id)
        {
            if (RefSpace::has_weights)
            {
                // if NURBS, transform the control points from  projective to euclidean coordinates
                const Real &w = weights_after_refinement_comp[loc_id];

                data_->control_points_[ctrl_pt_id] = ctrl_mesh_comp[loc_id] / w ;
            }
            else
                data_->control_points_[ctrl_pt_id] = ctrl_mesh_comp[loc_id];
        }
    }
    //----------------------------------

//    Assert(false,ExcNotImplemented());
//    AssertThrow(false,ExcNotImplemented());
#endif
}




template<class RefSpace>
auto
IgMapping<RefSpace>::
begin() const -> ElementIterator
{
    // TODO (pauletti, Apr 23, 2014): why not use this->shared_from_this()?
    //       is this a bug?
    return ElementIterator(
               const_pointer_cast<const self_t>(make_shared<self_t>(this->get_data())),0);
}



template<class RefSpace>
auto
IgMapping<RefSpace>::
last() const -> ElementIterator
{
//    return ElementIterator(
//               const_cast<self_t &>(*(new self_t(this->get_data()))),
//               this->get_grid()->get_num_active_elems() - 1);
    return ElementIterator(
               const_pointer_cast<const self_t>(make_shared<self_t>(this->get_data())),
               this->get_grid()->get_num_active_elems() - 1);
}



template<class RefSpace>
auto
IgMapping<RefSpace>::
end() const -> ElementIterator
{
//    return ElementIterator(
//               const_cast<self_t &>(*(new self_t(this->get_data()))),
//               IteratorState::pass_the_end);
    return ElementIterator(
               const_pointer_cast<const self_t>(make_shared<self_t>(this->get_data())),
               IteratorState::pass_the_end);
}



template<class RefSpace>
void
IgMapping<RefSpace>::
print_info(LogStream &out) const
{
    out << "Type = IgMapping" << endl;

    out.push("\t");
    out << "Reference space info:" << endl;
    data_->ref_space_->print_info(out);
    out << endl;

    out.push("\t");
    out << "Control points info (projective coordinates):" << endl;

    out.push("\t");

    // TODO (pauletti, Aug 26, 2014): does not satisfy printinfo standards, correct
    for (Index comp_id = 0 ; comp_id < space_dim ; ++comp_id)
    {
        out << "Control mesh["<<comp_id<<"] = ";
        data_->ctrl_mesh_[comp_id].print_info(out);
        out << endl;
    }
    out << endl;
    out.pop();


    out << "Control points info (euclidean coordinates): [ ";
    for (const auto &ctrl_pt : data_->control_points_)
        out << ctrl_pt << " ";
    out << "]" << endl;
    out.pop();
}

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/ig_mapping.inst>
