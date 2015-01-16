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


template<int dim_,int sp_dim_>
EvaluationPoints<dim_,sp_dim_>::
EvaluationPoints(const ValueVector<Point> &pts)
{
    this->reset_points_coordinates(pts);
}


template<int dim_,int sp_dim_>
void
EvaluationPoints<dim_,sp_dim_>::
reset_points_coordinates(const ValueVector<Point> &pts)
{
    const int n_pts = pts.size();
    Assert(n_pts > 0 , ExcEmptyObject());

    for (int i = 0 ; i < sp_dim_ ; ++i)
    {
        set<Real> coords_set;
        for (const auto &pt : pts)
            coords_set.emplace(pt[i]);

        //inserting the point coordinates and removing the duplicates
        coordinates_[i].assign(coords_set.begin(),coords_set.end());
    } // end loop i


    //for each point we retrieve the coordinate index
    map_point_id_to_coords_id_.clear();
    for (const auto &pt : pts)
    {
        TensorIndex<sp_dim_> coords_tensor_id;
        for (int i = 0 ; i < sp_dim_ ; ++i)
        {
            const auto coords_begin = coordinates_[i].begin();
            const auto coords_end   = coordinates_[i].end();

            const auto it = std::find(coords_begin, coords_end, pt[i]);

            coords_tensor_id[i] = std::distance(coords_begin,it);
        }
        map_point_id_to_coords_id_.emplace_back(coords_tensor_id);
    }
}


template<int dim_,int sp_dim_>
TensorIndex<sp_dim_>
EvaluationPoints<dim_,sp_dim_>::
get_coords_id_from_point_id(const int point_id) const
{
    Assert(point_id >= 0 && point_id < this->get_num_points(),
           ExcIndexRange(point_id,0,this->get_num_points()));

    return map_point_id_to_coords_id_[point_id];
}

template<int dim_,int sp_dim_>
int
EvaluationPoints<dim_,sp_dim_>::
get_num_points() const
{
    return map_point_id_to_coords_id_.size();
}



template<int dim_,int sp_dim_>
TensorIndex<sp_dim_>
EvaluationPoints<dim_,sp_dim_>::
get_num_coords_direction() const noexcept
{
    TensorIndex<sp_dim_> n_coords;
    for (int i = 0 ; i < sp_dim_ ; ++i)
        n_coords[i] = coordinates_[i].size();

    return n_coords;
}


template<int dim_>
QuadratureTensorProduct<dim_>::
QuadratureTensorProduct(
    const TensorSize<dim> num_points,
    void (*compute_coords_and_weight_1d)(const int n_pts_id, vector<Real> &coords,vector<Real> &weights),
    const Real eps_scaling)
    :
    points_(num_points),
    weights_(num_points)
{
    Assert(compute_coords_and_weight_1d != nullptr,ExcNullPtr());

    Assert(eps_scaling >= Real(0.0) && eps_scaling < Real(0.5),
           ExcMessage("The scaling factor must be >= 0.0 and < 0.5"));

    array<vector<Real>,dim> coords;
    array<vector<Real>,dim> weights;
    for (int i = 0; i < dim; ++i)
    {
        const auto n_pts = num_points[i];

        coords[i].resize(n_pts);
        weights[i].resize(n_pts);
        compute_coords_and_weight_1d(n_pts, coords[i], weights[i]);

        if (eps_scaling > 0)
            for (int ip = 0; ip < n_pts; ++ip)
                coords[i][ip] = 0.5 + (coords[i][ip] / 0.5 - 1.0) * (0.5 - eps_scaling) ;

        this->points_.copy_data_direction(i,coords[i]);
        this->weights_.copy_data_direction(i,weights[i]);
    }


    const int n_pts_total = num_points.flat_size();
    ValueVector<Point> points(n_pts_total);
    const auto n_pts_w = MultiArrayUtils<dim>::compute_weight(num_points);
    for (int pt_flat_id = 0 ; pt_flat_id < n_pts_total ; ++pt_flat_id)
    {
        const auto pt_tensor_id = MultiArrayUtils<dim>::flat_to_tensor_index(pt_flat_id,n_pts_w);

        for (int i = 0 ; i < dim ; ++i)
            points[pt_flat_id][i] = coords[i][pt_tensor_id[i]];
    }
    this->reset_points_coordinates(points);
}


template<int dim_>
QuadratureTensorProduct<dim_>::
QuadratureTensorProduct(const CartesianProductArray<Real,dim> &points,
                        const TensorProductArray<dim> &weights)
    :
    EvaluationPoints<dim_,dim_>(points.get_flat_cartesian_product()),
    points_(points),
    weights_(weights)
{
    Assert(points.tensor_size() == weights.tensor_size(),
           ExcMessage("Sizes of points and weights do not match."));
}

template<int dim_>
bool
QuadratureTensorProduct<dim_>::
is_tensor_product_struct() const
{
    return true;
}


template<int dim_>
auto
QuadratureTensorProduct<dim_>::
get_points() const noexcept -> PointArray
{
    return points_;
}



template<int dim_>
auto
QuadratureTensorProduct<dim_>::
get_weights() const noexcept -> WeigthArray
{
    return weights_;
}




/*
template<int dim_>
Size
QuadratureTensorProduct<dim_>::
get_num_points() const noexcept
{
    return points_.flat_size();
}
//*/


template<int dim_>
void
QuadratureTensorProduct<dim_>::
print_info(LogStream &out) const
{
    out << "Number of points:" << this->get_num_points() << endl;

    out << "weights:" << endl;
    get_weights().print_info(out);
    out << endl;

    // TODO (pauletti, Aug 26, 2014): redundant info, remove
    out << "weights (flat tensor product):" << endl;
    get_weights().get_flat_tensor_product().print_info(out);
    out << endl;

    out << "coordinates:" << endl;
    this->get_points().print_info(out);
    out << endl;

    // TODO (pauletti, Aug 26, 2014): redundant info, remove
    out << "coordinates (flat cartesian_product):" << endl;
    get_points().get_flat_cartesian_product().print_info(out);
    out << endl;

    out << endl;
}

template<int dim_>
template<int k>
auto
QuadratureTensorProduct<dim_>::
collapse_to_sub_element(const int sub_elem_id) const -> self_t
{
    auto &k_elem = UnitElement<dim>::template get_elem<k>(sub_elem_id);

    PointArray  new_points;
    WeigthArray new_weights;

    const int n_dir = k_elem.constant_directions.size();
    for (int j=0; j<n_dir; ++j)
    {
        auto dir = k_elem.constant_directions[j];
        auto val = k_elem.constant_values[j];
        new_points.copy_data_direction(dir,vector<Real>(1, val));
        new_weights.copy_data_direction(dir,vector<Real>(1, 1.0));

    }

    for (auto i : k_elem.active_directions)
    {
        new_points.copy_data_direction(i,points_.get_data_direction(i));
        new_weights.copy_data_direction(i,weights_.get_data_direction(i));
    }

    return QuadratureTensorProduct<dim>(new_points, new_weights);
}



template<int k, int dim>
QuadratureTensorProduct<dim>
extend_sub_elem_quad(const QuadratureTensorProduct<k> &quad,
                     const int sub_elem_id)
{

    auto &k_elem = UnitElement<dim>::template get_elem<k>(sub_elem_id);

    typename QuadratureTensorProduct<dim>::PointArray  new_points;
    typename QuadratureTensorProduct<dim>::WeigthArray new_weights;

    const int n_dir = k_elem.constant_directions.size();
    for (int j=0; j<n_dir; ++j)
    {
        auto dir = k_elem.constant_directions[j];
        auto val = k_elem.constant_values[j];
        new_points.copy_data_direction(dir,vector<Real>(1, val));
        new_weights.copy_data_direction(dir,vector<Real>(1, 1.0));

    }

    int ind = 0;
    for (auto i : k_elem.active_directions)
    {
        new_points.copy_data_direction(i,quad.get_points().get_data_direction(ind));
        new_weights.copy_data_direction(i,quad.get_weights().get_data_direction(ind));
        ++ind;
    }

    return QuadratureTensorProduct<dim>(new_points, new_weights);
}

IGA_NAMESPACE_CLOSE

#include <igatools/base/quadrature.inst>
