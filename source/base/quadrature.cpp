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


using std::endl;
using std::set;

IGA_NAMESPACE_OPEN




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
           void (*quad_1d)(int, iga::SafeSTLVector<double> &, iga::SafeSTLVector<double> &))
  :
  points_1d_(num_points),
  weights_1d_(num_points),
  is_tensor_product_(true)
{

  for (int i = 0; i<dim; ++i)
  {
    iga::SafeSTLVector<double> pts(num_points[i]);
    iga::SafeSTLVector<double> w(num_points[i]);
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
           const WeightArray &weights_1d,
           const BBox<dim_> &bounding_box)
  :
  points_1d_(points),
  weights_1d_(weights_1d),
  is_tensor_product_(true),
  bounding_box_(bounding_box)
{

#ifndef NDEBUG
  for (int i = 0 ; i < dim_ ; ++i)
  {
    // check that the points coordinate are within the bounding box
    const auto box_min = bounding_box_[i][0];
    const auto box_max = bounding_box_[i][1];
    for (const auto &coord : points_1d_.get_data_direction(i))
      Assert(coord >= box_min && coord <= box_max,
             ExcMessage("Point coordinate outside the bounding box."));


    Assert(points_1d_.get_data_direction(i).size() == weights_1d_.get_data_direction(i).size(),
           ExcDimensionMismatch(points_1d_.get_data_direction(i).size(),weights_1d_.get_data_direction(i).size()));

  } // end loop i
#endif


  const auto n_pts = points_1d_.flat_size();
  for (int i = 0 ; i < n_pts ; ++i)
    map_point_id_to_coords_id_.push_back(points_1d_.flat_to_tensor(i));
}




template<int dim_>
Quadrature<dim_>::
Quadrature(
  const PointVector &pts,
  const BBox<dim_> &bounding_box)
  :
  is_tensor_product_(false),
  bounding_box_(bounding_box)
{
  const int n_pts = pts.size();
  Assert(n_pts > 0 , ExcEmptyObject());

  //-----------------------------------------------------------------
  SafeSTLArray<std::map<Real,set<int>>,dim_> map_coords_point_ids;
  for (int i = 0 ; i < dim_ ; ++i)
  {
    set<Real> coords_set;
    for (const auto &pt : pts)
      coords_set.emplace(pt[i]);

    //inserting the point coordinates and removing the duplicates
    points_1d_.copy_data_direction(i,SafeSTLVector<Real>(coords_set.begin(),coords_set.end()));

    weights_1d_.copy_data_direction(i,SafeSTLVector<Real>(coords_set.size(),1.0));

#ifndef NDEBUG
    // check that the points coordinate are within the bounding box
    const auto box_min = bounding_box_[i][0];
    const auto box_max = bounding_box_[i][1];
    for (const auto &coord : points_1d_.get_data_direction(i))
      Assert(coord >= box_min && coord <= box_max,
             ExcMessage("Point coordinate " + std::to_string(i) + " outside the bounding box."));
#endif
  } // end loop i
  //-----------------------------------------------------------------



  //-----------------------------------------------------------------
  //for each point, for each coordinate direction, we retrieve the coordinate index;
  map_point_id_to_coords_id_.resize(n_pts);
  for (int j = 0 ; j < n_pts ; ++j)
  {
    const auto &pt = pts[j];
    TensorIndex<dim_> coords_tensor_id;
    for (int i = 0 ; i < dim_ ; ++i)
    {
      const auto coords_begin = points_1d_.get_data_direction(i).begin();
      const auto coords_end   = points_1d_.get_data_direction(i).end();

      const auto it = std::find(coords_begin, coords_end, pt[i]);

      coords_tensor_id[i] = std::distance(coords_begin,it);
    }
    map_point_id_to_coords_id_[j] = (coords_tensor_id);
  }
  //-----------------------------------------------------------------
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
dilate(const Point &dilation_factor)
{
  points_1d_.dilate(dilation_factor);
  weights_1d_.dilate(dilation_factor);

  bounding_box_.dilate(dilation_factor);
}



template<int dim_>
void
Quadrature<dim_>::
translate(const Point &translation)
{
  points_1d_.translate(translation);

  bounding_box_.translate(translation);
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
translate_dilate(const Point &translate,const Point &dilate)
{
  this->translate(translate);
  this->dilate(dilate);
}

template<int dim_>
bool
Quadrature<dim_>::
is_tensor_product() const
{
  return this->is_tensor_product_;
}



template<int dim_>
const SafeSTLVector<Real> &
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
  for (int pt = 0 ; pt < n_pts ; ++pt)
    points[pt] = this->get_point(pt);

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
  for (int pt = 0 ; pt < n_pts ; ++pt)
    weights[pt] = this->get_weight(pt);

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

  out.begin_item("Is tensor product: " + std::string(is_tensor_product_?"TRUE":"FALSE"));
  out.end_item();

  out.begin_item("Map points id --- coordinates id:");
  map_point_id_to_coords_id_.print_info(out);
  out.end_item();

  out.begin_item("Bounding box:");
  bounding_box_.print_info(out);
  out.end_item();
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
  BBox<dim_> new_bounding_box;

  const int n_dir = k_elem.constant_directions.size();
  for (int j = 0 ; j < n_dir; ++j)
  {
    auto dir = k_elem.constant_directions[j];
    auto val = k_elem.constant_values[j];

    new_coords_1d.copy_data_direction(dir, SafeSTLVector<Real>(1, val));
    new_weights_1d.copy_data_direction(dir, SafeSTLVector<Real>(1, 1.));

    new_bounding_box[dir][0] = val;
    new_bounding_box[dir][1] = val;

  }

  for (auto i : k_elem.active_directions)
  {
    new_coords_1d.copy_data_direction(i, points_1d_.get_data_direction(i));
    new_weights_1d.copy_data_direction(i, weights_1d_.get_data_direction(i));

    new_bounding_box[i] = bounding_box_[i];
  }

  return self_t(new_coords_1d, new_weights_1d,new_bounding_box);

  return self_t();
}



template<int sub_dim, int dim>
Quadrature<dim>
extend_sub_elem_quad(const Quadrature<sub_dim> &eval_pts,const int sub_elem_id)
{
  using WeightArray = typename Quadrature<dim>::WeightArray;
  using PointArray = typename Quadrature<dim>::PointArray;

  auto &subdim_elem = UnitElement<dim>::template get_elem<sub_dim>(sub_elem_id);


  const auto &old_bounding_box = eval_pts.get_bounding_box();
  BBox<dim> new_bounding_box;

  const int n_new_dirs = subdim_elem.constant_directions.size();


  if (eval_pts.is_tensor_product())
  {
    const auto &old_points_1d = eval_pts.get_points_1d();
    PointArray new_coords_1d;

    const auto &old_weights_1d = eval_pts.get_weights_1d();
    WeightArray new_weights_1d;

    for (int j = 0 ; j < n_new_dirs ; ++j)
    {
      const auto dir = subdim_elem.constant_directions[j];
      const auto val = subdim_elem.constant_values[j];
      new_coords_1d.copy_data_direction(dir, SafeSTLVector<Real>(1, val));
      new_weights_1d.copy_data_direction(dir, SafeSTLVector<Real>(1, 1.));

      new_bounding_box[dir][0] = val;
      new_bounding_box[dir][1] = val;
    }
    int ind = 0;
    for (auto i : subdim_elem.active_directions)
    {
      new_coords_1d.copy_data_direction(i, old_points_1d.get_data_direction(ind));
      new_weights_1d.copy_data_direction(i, old_weights_1d.get_data_direction(ind));

      new_bounding_box[i] = old_bounding_box[ind];

      ++ind;
    }

    return Quadrature<dim>(new_coords_1d,new_weights_1d,new_bounding_box);
  } // if (eval_pts.is_tensor_product())
  else
  {
    const auto &old_points = eval_pts.get_points();

    const auto n_pts = old_points.size();

    ValueVector<Points<dim>> new_points(n_pts);
    for (int pt = 0 ; pt < n_pts ; ++pt)
    {
      for (int j = 0 ; j < n_new_dirs ; ++j)
      {
        const auto dir = subdim_elem.constant_directions[j];
        const auto val = subdim_elem.constant_values[j];
        new_points[pt][dir] = val;

        new_bounding_box[dir][0] = val;
        new_bounding_box[dir][1] = val;
      }
      int ind = 0;
      for (auto i : subdim_elem.active_directions)
      {
        new_points[pt][i] = old_points[pt][ind];

        new_bounding_box[i] = old_bounding_box[ind];

        ++ind;
      }

    } // end loop pt

    return Quadrature<dim>(new_points,new_bounding_box);
  } // if (!eval_pts.is_tensor_product())
}


#ifdef SERIALIZATION

template<int dim_>
template<class Archive>
void
Quadrature<dim_>::
serialize(Archive &ar, const unsigned int version)
{
  ar &boost::serialization::make_nvp("points_1d_",points_1d_);
  ar &boost::serialization::make_nvp("weights_1d_",weights_1d_);
  ar &boost::serialization::make_nvp("map_point_id_to_coords_id_",map_point_id_to_coords_id_);
  ar &boost::serialization::make_nvp("is_tensor_product_",is_tensor_product_);
  ar &boost::serialization::make_nvp("bounding_box_",bounding_box_);
};

#endif // SERIALIZATION


IGA_NAMESPACE_CLOSE

#include <igatools/base/quadrature.inst>
