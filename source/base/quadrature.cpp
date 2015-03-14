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

#include <igatools/base/quadrature.h>
#include <igatools/base/exceptions.h>
#include <igatools/geometry/unit_element.h>
#include <igatools/utils/multi_array_utils.h>

#include <set>

using std::array;
using std::endl;
using std::set;

IGA_NAMESPACE_OPEN

template<int dim_>
Quadrature<dim_>::
Quadrature()
{
    //  for (auto &box_direction : bounding_box_)
    //    {
    //        box_direction[0] = 0.0;
    //        box_direction[1] = 1.0;
    //    }
}



template<int dim_>
Quadrature<dim_>::
Quadrature(const BBox<dim_> &bounding_box)
    :
    is_tensor_product_(false),
    bounding_box_(bounding_box)
{}

template<int dim_>
Quadrature<dim_>::
Quadrature(const TensorSize<dim> &num_points,
           void (*quad_1d)(int, iga::vector<double> &, iga::vector<double> &))
    :
    points_1d_(num_points),
    weights_1d_(num_points),
    is_tensor_product_(true)
{

    for (int i = 0; i<dim; ++i)
    {
        iga::vector<double> pts(num_points[i]);
        iga::vector<double> w(num_points[i]);
        quad_1d(num_points[i],pts,w);
        points_1d_.copy_data_direction(i, pts);
        weights_1d_.copy_data_direction(i, w);
    }

    const auto n_pts = points_1d_.flat_size();
    for (int i = 0 ; i < n_pts ; ++i)
        map_point_id_to_coords_id_.push_back(points_1d_.flat_to_tensor(i));
}



template<int dim_>
Quadrature<dim_>::
Quadrature(const PointArray &points,
           const WeightArray &weights_1d)
    :
    points_1d_(points),
    weights_1d_(weights_1d),
    is_tensor_product_(true)
{
    const auto n_pts = points_1d_.flat_size();
    for (int i = 0 ; i < n_pts ; ++i)
        map_point_id_to_coords_id_.push_back(points_1d_.flat_to_tensor(i));
}



template<int dim_>
Quadrature<dim_>::
Quadrature(const ValueVector<Point> &pts)
    :
    is_tensor_product_(false)
{
    const int n_pts = pts.size();
    TensorSize<dim_> size(n_pts);
    WeightArray weights_1d(size);

    Assert(false, ExcMessage("put weight to 1"));
    this->reset_points_points_1d_and_weights(pts,weights_1d);
}



template<int dim_>
Quadrature<dim_>::
Quadrature(
    const PointVector &pts,
    const WeightArray &weights_1d,
    const BBox<dim_> &bounding_box)
    :
    Quadrature(bounding_box)
{
    this->reset_points_points_1d_and_weights(pts,weights_1d);
}



template<int dim_>
const BBox<dim_> &
Quadrature<dim_>::
get_bounding_box() const
{
    return bounding_box_;
}



template<int dim_>
void
Quadrature<dim_>::
reset_bounding_box(const BBox<dim_> &bounding_box)
{
    bounding_box_ = bounding_box;

#ifndef NDEBUG
    for (const auto &box_direction : bounding_box_)
        Assert(box_direction[0] <= box_direction[1],
               ExcMessage("Wrong coordinates for the bounding box."));
#endif
}



template<int dim_>
void
Quadrature<dim_>::
dilate(const Point &dilate)
{
    points_1d_.dilate(dilate);
    weights_1d_.dilate(dilate);


    for (int i = 0 ; i < dim_ ; ++i)
    {
        Assert(dilate[i] > 0., ExcMessage("Dilation factor must be positive."));
        bounding_box_[i][0] *= dilate[i];
        bounding_box_[i][1] *= dilate[i];
    }
}



template<int dim_>
void
Quadrature<dim_>::
translate(const Point &translate)
{
    points_1d_.translate(translate);

    // TODO (pauletti, Feb 27, 2015): code BBox translate
    for (int i = 0 ; i < dim_ ; ++i)
    {
        bounding_box_[i][0] += translate[i];
        bounding_box_[i][1] += translate[i];
    }
}



template<int dim_>
void
Quadrature<dim_>::
dilate_translate(const Point &dilate, const Point &translate)
{
    this->dilate(dilate);
    this->translate(translate);
}



template<int dim_>
void
Quadrature<dim_>::
reset_points_points_1d_and_weights(
    const PointVector &pts,
    const WeightArray &weights_1d)
{
    Assert(!is_tensor_product_, ExcNotImplemented());
    const int n_pts = pts.size();
    Assert(n_pts > 0 , ExcEmptyObject());

    //-----------------------------------------------------------------
    TensorSize<dim_> n_dirs;
    for (int i = 0 ; i < dim_ ; ++i)
    {
        set<Real> coords_set;
        for (const auto &pt : pts)
            coords_set.emplace(pt[i]);

        //inserting the point coordinates and removing the duplicates
        points_1d_.copy_data_direction(i,vector<Real>(coords_set.begin(),coords_set.end()));
        weights_1d_ = weights_1d;

#ifndef NDEBUG
        // check that the points coordinate are within the bounding box
        const auto box_min = bounding_box_[i][0];
        const auto box_max = bounding_box_[i][1];
        for (const auto &coord : points_1d_.get_data_direction(i))
            Assert(coord >= box_min && coord <= box_max,
                   ExcMessage("Point coordinate outside the bounding box."));
#endif


//        Assert(n_dirs[i] == weights_1d[i].size(),
//               ExcDimensionMismatch(n_dirs[i],weights_1d[i].size()));

    } // end loop i
    //-----------------------------------------------------------------



    //-----------------------------------------------------------------
    //for each point :
    //  - we retrieve the coordinate index;
    //  - we set the associated weight
    map_point_id_to_coords_id_.resize(n_pts);
    for (int j = 0 ; j < n_pts ; ++j)
    {
        auto &pt = pts[j];
        TensorIndex<dim_> coords_tensor_id;
        for (int i = 0 ; i < dim_ ; ++i)
        {
            const auto coords_begin = points_1d_.get_data_direction(i).begin();
            const auto coords_end   = points_1d_.get_data_direction(i).end();

            const auto it = std::find(coords_begin, coords_end, pt[i]);

            coords_tensor_id[i] = std::distance(coords_begin,it);
        }
        map_point_id_to_coords_id_[j]=(coords_tensor_id);
    }
    //-----------------------------------------------------------------
}



template<int dim_>
bool
Quadrature<dim_>::
is_tensor_product() const
{
    return this->is_tensor_product_;
}



template<int dim_>
const vector<Real> &
Quadrature<dim_>::
get_coords_direction(const int i) const
{
    return points_1d_.get_data_direction(i);
}



template<int dim_>
TensorIndex<dim_>
Quadrature<dim_>::
get_coords_id_from_point_id(const int point_id) const
{
    Assert(point_id >= 0 && point_id < this->get_num_points(),
           ExcIndexRange(point_id,0,this->get_num_points()));

    return map_point_id_to_coords_id_[point_id];
}



template<int dim_>
auto
Quadrature<dim_>::
get_point(const int pt_id) const -> Point
{
    const auto tensor_id = this->get_coords_id_from_point_id(pt_id);
    return points_1d_.cartesian_product(tensor_id);
}



template<int dim_>
auto
Quadrature<dim_>::
get_points() const -> ValueVector<Point>
{
    const int n_pts = this->get_num_points();
    ValueVector<Point> points(n_pts);
    for (int ipt = 0 ; ipt < n_pts ; ++ipt)
        points[ipt] = this->get_point(ipt);

    return points;
}



template<int dim_>
Real
Quadrature<dim_>::
get_weight(const int pt_id) const
{
    const auto tensor_id = this->get_coords_id_from_point_id(pt_id);
    return weights_1d_.tensor_product(tensor_id);
}



template<int dim_>
ValueVector<Real>
Quadrature<dim_>::
get_weights() const
{
    const int n_pts = this->get_num_points();
    ValueVector<Real> weights(n_pts);
    for (int ipt = 0 ; ipt < n_pts ; ++ipt)
        weights[ipt] = this->get_weight(ipt);

    return weights;
}



template<int dim_>
auto
Quadrature<dim_>::
get_weights_1d() const -> const WeightArray &
{
    return weights_1d_;
}



template<int dim_>
auto
Quadrature<dim_>::
get_points_1d() const -> const PointArray &
{
    return points_1d_;
}



template<int dim_>
int
Quadrature<dim_>::
get_num_points() const
{
    return map_point_id_to_coords_id_.size();
}



template<int dim_>
TensorSize<dim_>
Quadrature<dim_>::
get_num_coords_direction() const noexcept
{
    return points_1d_.tensor_size();
}



template<int dim_>
void
Quadrature<dim_>::
print_info(LogStream &out) const
{
    out << "Number of points:" << this->get_num_points() << endl;

    out.begin_item("Weights:");
    get_weights().print_info(out);
    out.end_item();
    out << endl;

    out.begin_item("Coordinates:");
    for (int dir = 0 ; dir < dim_ ; ++dir)
    {
        out << "Direction: " << dir << endl;
        this->get_coords_direction(dir).print_info(out);
        out << endl;
    }
    out.end_item();

    out.begin_item("Points:");
    this->get_points().print_info(out);
    out.end_item();
    out << endl;

}



template<int dim_>
template<int k>
auto
Quadrature<dim_>::
collapse_to_sub_element(const int sub_elem_id) const -> Quadrature<dim_>
{
    Assert(is_tensor_product_, ExcNotImplemented());

    auto &k_elem = UnitElement<dim_>::template get_elem<k>(sub_elem_id);

    PointArray new_coords_1d;
    WeightArray new_weights_1d;
    const int n_dir = k_elem.constant_directions.size();
    for (int j = 0 ; j < n_dir; ++j)
    {
        auto dir = k_elem.constant_directions[j];
        auto val = k_elem.constant_values[j];

        new_coords_1d.copy_data_direction(dir, vector<Real>(1, val));
        new_weights_1d.copy_data_direction(dir, vector<Real>(1, 1.));

    }

    for (auto i : k_elem.active_directions)
    {
        new_coords_1d.copy_data_direction(i, points_1d_.get_data_direction(i));
        new_weights_1d.copy_data_direction(i, weights_1d_.get_data_direction(i));
    }

    return self_t(new_coords_1d, new_weights_1d);

    return self_t();
}



template<int k, int dim>
Quadrature<dim>
extend_sub_elem_quad(const Quadrature<k> &eval_pts,
                     const int sub_elem_id)
{

    Assert(eval_pts.is_tensor_product(), ExcNotImplemented());
    using WeightArray = typename Quadrature<dim>::WeightArray;
    using PointArray = typename Quadrature<dim>::PointArray;

    auto &k_elem = UnitElement<dim>::template get_elem<k>(sub_elem_id);

    PointArray new_coords_1d;
    WeightArray new_weights_1d;


    const auto &old_weights_1d = eval_pts.get_weights_1d();
    const auto &old_points_1d = eval_pts.get_points_1d();

    const int n_dir = k_elem.constant_directions.size();
    for (int j = 0 ; j < n_dir ; ++j)
    {
        const auto dir = k_elem.constant_directions[j];
        const auto val = k_elem.constant_values[j];
        new_coords_1d.copy_data_direction(dir, vector<Real>(1, val));
        new_weights_1d.copy_data_direction(dir, vector<Real>(1, 1.));

    }
    int ind = 0;
    for (auto i : k_elem.active_directions)
    {
        new_coords_1d.copy_data_direction(i, old_points_1d.get_data_direction(ind));
        new_weights_1d.copy_data_direction(i, old_weights_1d.get_data_direction(ind));
        ++ind;
    }


    return Quadrature<dim>(new_coords_1d, new_weights_1d);
}


IGA_NAMESPACE_CLOSE

#include <igatools/base/quadrature.inst>
