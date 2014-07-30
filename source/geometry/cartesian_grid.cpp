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

#include <igatools/geometry/cartesian_grid.h>
#include <igatools/base/exceptions.h>
#include <igatools/base/array_utils.h>
#include <igatools/utils/vector_tools.h>
#include <igatools/utils/multi_array_utils.h>

#include <set>
#include <algorithm>

using std::endl;
using std::array;
using std::vector;
using std::set;
using std::shared_ptr;
using std::unique_ptr;
using std::make_shared;

IGA_NAMESPACE_OPEN

namespace
{

/**
 * Returns a uniform filled knots
 */
template <int dim>
CartesianProductArray<Real,dim>
filled_progression(const BBox<dim> &end_points, const TensorSize<dim> &n_knots)
{
    CartesianProductArray<Real,dim> knot_coordinates(n_knots);

    vector<Real> knots_1d;
    for (int i = 0; i < dim; ++i)
    {
        const Size n_i = n_knots(i);
        Assert(n_i > 1, ExcLowerRange(n_i,2));

        knots_1d.resize(n_i);

        const Real h=(end_points[i][1] - end_points[i][0]) /(n_i-1);

        knots_1d[0] = end_points[i][0];
        for (int j = 1; j < n_i; ++j)
            knots_1d[ j ] = knots_1d[ j-1 ] + h;

        knot_coordinates.copy_data_direction(i,knots_1d);
    }
    return knot_coordinates;
}
}



template<int dim_>
CartesianGrid<dim_>::
CartesianGrid(const Size n)
    :
    CartesianGrid(TensorSize<dim>(n), Kind::uniform)
{}



template<int dim_>
auto
CartesianGrid<dim_>::
create(const Size n) -> shared_ptr<self_t>
{
    return shared_ptr<self_t>(new self_t(n));
}



template<int dim_>
CartesianGrid<dim_>::
CartesianGrid(const TensorSize<dim> &n, const Kind kind)
    :
    CartesianGrid(filled_array<array<Real,2>,dim>(array<Real,2> {{0,1}}), n, kind)
{}



template<int dim_>
auto
CartesianGrid<dim_>::
create(const TensorSize<dim> &n) -> shared_ptr<self_t>
{
    return shared_ptr<self_t>(new self_t(n, Kind::direction_uniform));
}



template<int dim_>
CartesianGrid<dim_>::
CartesianGrid(const BBox<dim> &end_points, const Size n_knots, const Kind kind)
    :
    CartesianGrid(end_points, TensorSize<dim>(n_knots), kind)
{}



template<int dim_>
auto
CartesianGrid<dim_>::
create(const BBox<dim> &end_points, const Size n_knots) -> shared_ptr<self_t>
{
    return shared_ptr<self_t>(new self_t(end_points, n_knots, Kind::uniform));
}



template<int dim_>
CartesianGrid<dim_>::
CartesianGrid(const BBox<dim> &end_points,
              const TensorSize<dim> &n,
              const Kind kind)
    :
    CartesianGrid(filled_progression<dim>(end_points, n), kind)
{}



template<int dim_>
auto
CartesianGrid<dim_>::
create(const BBox<dim> &end_points,
       const TensorSize<dim> &n) -> shared_ptr<self_t>
{
    return shared_ptr<self_t>(new self_t(end_points, n, Kind::direction_uniform));
}



template<int dim_>
CartesianGrid<dim_>::
CartesianGrid(const KnotCoordinates &knot_coordinates,
              const Kind kind)
    :
    kind_(kind),
    boundary_id_(filled_array<int,UnitElement<dim>::faces_per_element>(0)),
    knot_coordinates_(knot_coordinates)
{
#ifndef NDEBUG
    for (int i = 0; i < dim; i++)
    {
        const auto &knots_i = knot_coordinates.get_data_direction(i);
        // checks that we have at least two knot values (i.e. one knot span) in
        // each coordinate direction

        AssertThrow(knots_i.size() > 1,
                    ExcLowerRange(knots_i.size(), 2));


        // check if the array is sorted and does not contains duplicates
        vector< Real > vec = knots_i ;
        std::sort(vec.begin(), vec.end());
        vec.erase(unique(vec.begin(), vec.end()), vec.end());
        AssertThrow(knots_i == vec,
                    ExcMessage("The knot coordinate vector is not sorted and/or contains duplicates"));

    }
#endif

    num_elem_ = knot_coordinates_.tensor_size();
    num_elem_ -= 1;
    weight_elem_id_ = MultiArrayUtils<dim>::compute_weight(num_elem_);

    influent_.resize(num_elem_, true);
    active_elems_.resize(num_elem_, true);
}



template<int dim_>
auto
CartesianGrid<dim_>::
create(const KnotCoordinates &knot_coordinates) -> shared_ptr<self_t>
{
    return shared_ptr<self_t>(new self_t(knot_coordinates, Kind::non_uniform));
}



template<int dim_>
CartesianGrid<dim_>::
CartesianGrid(const array<vector<Real>,dim> &knot_coordinates)
    :
    self_t(CartesianProductArray<Real,dim>(knot_coordinates),
           Kind::direction_uniform)
{}



template<int dim_>
auto
CartesianGrid<dim_>::
create(const array<vector<Real>,dim> &knot_coordinates) -> shared_ptr<self_t>
{
    return shared_ptr< self_t >(new self_t(knot_coordinates));
}



template<int dim_>
CartesianGrid<dim_>::
CartesianGrid(const self_t &grid)
    :
    kind_(grid.kind_),
    boundary_id_(grid.boundary_id_),
    knot_coordinates_(grid.knot_coordinates_),
    weight_elem_id_(grid.weight_elem_id_),
    num_elem_(grid.num_elem_),
    active_elems_(grid.active_elems_),
    influent_(grid.influent_)
{}



//TODO: inline this function
template<int dim_>
vector< Real > const &
CartesianGrid<dim_>::
get_knot_coordinates(const int i) const
{
    return (knot_coordinates_.get_data_direction(i));
}



//TODO: inline this function
template<int dim_>
auto
CartesianGrid<dim_>::
get_knot_coordinates() const -> CartesianProductArray<Real,dim> const &
{
    return knot_coordinates_;
}



template<int dim_>
auto
CartesianGrid<dim_>::
get_element_lengths() const -> KnotCoordinates
{
    auto const &size = get_num_elements_dim();
    KnotCoordinates length(size);
    for (int i = 0; i < dim; ++i)
    {
        const auto &knots_i = knot_coordinates_.get_data_direction(i);

        const Size size_i = size(i);

        for (int j = 0 ; j < size_i ; ++j)
        {
            length.entry(i,j) = knots_i[j+1] - knots_i[j];
        }
    }
    return length;
}



template<int dim_>
auto
CartesianGrid<dim_>::begin() const -> ElementIterator
{
    return ElementIterator(this->shared_from_this(), 0);
}



template<int dim_>
auto
CartesianGrid<dim_>::
end() const -> ElementIterator
{
    return ElementIterator(this->shared_from_this(),
                           IteratorState::pass_the_end);
}



template<int dim_>
void
CartesianGrid<dim_>::
set_boundary_id(const int face, const boundary_id id)
{
    Assert(face < UnitElement<dim>::faces_per_element,
           ExcIndexRange(face,0, UnitElement<dim>::faces_per_element));
    boundary_id_[face] = id;
}



template<int dim_>
boundary_id
CartesianGrid<dim_>::
get_boundary_id(const int face) const
{
    Assert(face < UnitElement<dim>::faces_per_element,
           ExcIndexRange(face,0, UnitElement<dim>::faces_per_element));
    return (boundary_id_[face]);
}



template<int dim_>
auto
CartesianGrid<dim_>::
get_face_normal(const int face_no) const -> Points<dim>
{
    Points<dim> normal;
    normal[UnitElement<dim>::face_to_component[face_no][0]] =
    UnitElement<dim>::face_normal_direction[face_no];

    return normal;
}



template<int dim_>
Size
CartesianGrid<dim_>::
get_num_elements() const
{
    return std::count(active_elems_.begin(), active_elems_.end(), true);
}



// TODO (pauletti, Jul 30, 2014): this should not be necesary to specialize
template<>
Size
CartesianGrid<0>::
get_num_elements() const
{
    return 1;
}



template<int dim_>
auto
CartesianGrid<dim_>::
get_num_elements_dim() const -> TensorSize<dim>
{
    return num_elem_;
}



template<int dim_>
auto
CartesianGrid<dim_>::
get_num_knots_dim() const -> TensorSize<dim>
{
    return knot_coordinates_.tensor_size();
}


template<int dim_>
auto
CartesianGrid<dim_>::
get_grid_pre_refinement() const -> shared_ptr<const self_t>
{
    return grid_pre_refinement_;
}



template <int dim_>
void
CartesianGrid<dim_>::
refine_directions(
    const array<bool,dim> &refinement_directions,
    const array<Size,dim> &n_subdivisions)
{
    // make a copy of the grid before the refinement
    grid_pre_refinement_ = make_shared<const self_t>(self_t(*this));

    for (int i = 0 ; i < dim ; ++i)
        if (refinement_directions[i])
            this->refine_knots_direction(i,n_subdivisions[i]);

    num_elem_ = knot_coordinates_.tensor_size();
    num_elem_ -= 1;
    weight_elem_id_ = MultiArrayUtils<dim>::compute_weight(num_elem_);

    // TODO (pauletti, Jul 30, 2014): this is wrong in general !!!
    influent_.resize(num_elem_, true);
    active_elems_.resize(num_elem_, true);

    // refining the objects that's are attached to the CartesianGrid
    // (i.e. that are defined using this CartesianGrid object)
    this->refine_signals_(refinement_directions,*grid_pre_refinement_);
}



template <int dim_>
void
CartesianGrid<dim_>::
refine_direction(const int direction_id, const Size n_subdivisions)
{
    Assert(direction_id >= 0 && direction_id < dim,
           ExcIndexRange(direction_id, 0, dim));

    array<bool,dim> refinement_directions = filled_array<bool,dim>(false);
    refinement_directions[direction_id] = true;

    array<Size,dim> n_subdiv;
    n_subdiv[direction_id] = n_subdivisions;

    this->refine_directions(refinement_directions,n_subdiv);
}



template <int dim_>
void
CartesianGrid<dim_>::
refine(const Size n_subdivisions)
{
    Assert(n_subdivisions >= 2, ExcLowerRange(n_subdivisions,2));

    this->refine_directions(
        filled_array<bool,dim>(true),
        filled_array<Size,dim>(n_subdivisions));
}



template <int dim_>
boost::signals2::connection
CartesianGrid<dim_>::
connect_refinement(const SignalRefineSlot &subscriber)
{
    return refine_signals_.connect(subscriber);
}


template <int dim_>
void
CartesianGrid<dim_>::
refine_knots_direction(const int direction_id,
                       const Size n_subdivisions)
{
    Assert(n_subdivisions >= 2,ExcLowerRange(n_subdivisions,2));

    Assert(direction_id >=0 && direction_id < dim,
           ExcIndexRange(direction_id,0,dim));

    const auto &knots_old = this->get_knot_coordinates(direction_id);

    const Size n_knots_old = knots_old.size();
    const Size n_knots_to_add = (n_knots_old - 1) * (n_subdivisions - 1);

    const Size n_knots_new = n_knots_old + n_knots_to_add;
    vector<Real> knots_new(n_knots_new);

    Index i_new = 0;
    for (Index i_old = 0 ; i_old < n_knots_old - 1 ; ++i_old)
    {
        const Real h = (knots_old[i_old+1] - knots_old[i_old]) / n_subdivisions ;

        for (Index j = 0 ; j < n_subdivisions ; ++j, ++i_new)
            knots_new[i_new] = knots_old[i_old] + j * h;
    }
    knots_new[n_knots_new-1] = knots_old[n_knots_old-1];

    knot_coordinates_.copy_data_direction(direction_id,knots_new);


}



template <int dim_>
void
CartesianGrid<dim_>::
print_info(LogStream &out) const
{
    out << "Number of active elements: " << get_num_elements() << endl;
    out << "Number of intervals per direction: " << num_elem_ << endl;

    out.begin_item("Knot coordinates:");
    knot_coordinates_.print_info(out);
    out.end_item();
}



template <int dim_>
auto
CartesianGrid<dim_>::
get_face_grid(const int face_id, std::map<int,int> &elem_map) const
-> shared_ptr<FaceType>
{
    Assert(dim > 0, ExcLowerRange(dim,1));

    const auto active_dirs = UnitElement<dim>::face_active_directions[face_id];
    const int const_dir = UnitElement<dim>::face_constant_direction[face_id];
    const int const_value = UnitElement<dim>::face_side[face_id];

    TensorIndex<dim> v_index;
    auto knot_coordinates_ = this->get_knot_coordinates();
    v_index[const_dir] = const_value==0 ?
    0 :
    (knot_coordinates_.tensor_size()(const_dir)-2);

    auto face_knots = knot_coordinates_.get_sub_product(active_dirs);
    auto face_grid = FaceType::create(face_knots);

    auto v_elem = this->begin();
    auto f_elem = face_grid->begin();
    auto end = face_grid->end();
    for (; f_elem != end; ++f_elem)
    {
        auto f_index = f_elem->get_tensor_index();
        for (int j=0; j<dim-1; ++j)
            v_index[active_dirs[j]] = f_index[j];
        v_elem->reset_flat_tensor_indices(v_index);
        elem_map[f_elem->get_flat_index()]=v_elem->get_flat_index();
    }

    return face_grid;
}



template <int dim_>
auto
CartesianGrid<dim_>::
get_bounding_box() const -> BBox<dim>
{
    BBox<dim> bounding_box;

    for (int i = 0 ; i < dim ; ++i)
    {
        bounding_box[i][0] = knot_coordinates_.get_data_direction(i).front();
        bounding_box[i][1] = knot_coordinates_.get_data_direction(i).back();
    }

    return bounding_box;
}


template <int dim_>
Index
CartesianGrid<dim_>::
tensor_to_flat_element_index(const TensorIndex<dim> &tensor_id) const
{
    return MultiArrayUtils<dim>::tensor_to_flat_index(tensor_id,weight_elem_id_);
}



template <int dim_>
auto
CartesianGrid<dim_>::
flat_to_tensor_element_index(const Index flat_id) const ->TensorIndex<dim>
{
    return MultiArrayUtils<dim>::flat_to_tensor_index(flat_id,weight_elem_id_);
}



template <int dim_>
Index
CartesianGrid<dim_>::
get_element_flat_id_from_point(const Points<dim> &point) const
{
#ifndef NDEBUG
    const auto bounding_box = this->get_bounding_box();
#endif

    TensorIndex<dim> elem_t_id;
    for (int i = 0 ; i < dim ; ++i)
    {
        Assert(point[i] >= bounding_box[i][0] && point[i] <= bounding_box[i][1],
               ExcMessage("Point " +
                          std::to_string(point[i]) +
                          " outside the domain [" +
                          std::to_string(bounding_box[i][0]) + "," +
                          std::to_string(bounding_box[i][1])+
                          "]"));

        const auto &knots = knot_coordinates_.get_data_direction(i);

        //find the index j in the knots for which knots[j] <= point[i]
        const auto low = std::lower_bound(knots.begin(),knots.end(),point[i]);
        const Index j = low - knots.begin();

        elem_t_id[i] = (j>0) ? j-1 : 0;
    }


    return this->tensor_to_flat_element_index(elem_t_id);
}



template <int dim>
bool
CartesianGrid<dim>::
operator==(const CartesianGrid<dim> &grid) const
{
    bool same_knots_coordinates = true;
    for (int i = 0 ; i < dim ; ++i)
    {
        const auto &knots_a =  this->knot_coordinates_.get_data_direction(i);
        const auto &knots_b =   grid.knot_coordinates_.get_data_direction(i);

        same_knots_coordinates = same_knots_coordinates && (knots_a == knots_b);
    }
    return same_knots_coordinates;
}



// TODO (pauletti, Jul 30, 2014): why is this necessary?
template <int dim_>
const int
CartesianGrid<dim_>::dim ;



template <int dim>
vector<Index>
build_map_elements_between_cartesian_grids(
    const CartesianGrid<dim> &grid_fine,
    const CartesianGrid<dim> &grid_coarse)
{
    //---------------------------------------------------------
    // checks that the grid are on the same domain
    Assert(grid_fine.get_bounding_box() == grid_coarse.get_bounding_box(),
           ExcMessage("Grids on different domains."));
    //---------------------------------------------------------



    //---------------------------------------------------------
    array<vector<int>,dim> map_interv_fid_fine_coarse;
    for (int i = 0 ; i < dim ; ++i)
    {
        const auto &coords_coarse = grid_coarse.get_knot_coordinates(i);
        const auto &coords_fine = grid_fine.get_knot_coordinates(i);

#ifndef NDEBUG
        const int n_intervals_coarse = coords_coarse.size() - 1;
#endif
        const int n_intervals_fine = coords_fine.size() - 1;

        for (int fid_fine = 0 ; fid_fine < n_intervals_fine ; ++fid_fine)
        {
            int fid_coarse = 0;
            while (!(coords_fine[fid_fine] >= coords_coarse[fid_coarse] &&
                     coords_fine[fid_fine+1] <= coords_coarse[fid_coarse+1]))
            {
                ++fid_coarse;
            }
            Assert(fid_coarse < n_intervals_coarse,
                   ExcMessage("Impossible to find an interval "
                              "on the coarse grid that fully contains the interval " +
                              std::to_string(fid_fine) + " along the direction "
                              + std::to_string(i) +
                              "of the fine grid."));

            map_interv_fid_fine_coarse[i].push_back(fid_coarse);
        }
    }

    const int n_elems_fine = grid_fine.get_num_elements();
    vector<int> map_elem_fine_to_elem_coarse(n_elems_fine);
    for (int elem_fine_fid = 0 ; elem_fine_fid < n_elems_fine ; ++elem_fine_fid)
    {
        TensorIndex<dim> elem_fine_tid =
                grid_fine.flat_to_tensor_element_index(elem_fine_fid);

        TensorIndex<dim> elem_coarse_tid;
        for (int i = 0 ; i < dim ; ++i)
            elem_coarse_tid[i] = map_interv_fid_fine_coarse[i][elem_fine_tid[i]];

        const int elem_coarse_fid =
                grid_coarse.tensor_to_flat_element_index(elem_coarse_tid);

        map_elem_fine_to_elem_coarse[elem_fine_fid] = elem_coarse_fid;
    }
    //---------------------------------------------------------


    return map_elem_fine_to_elem_coarse;
}

// TODO (pauletti, Jul 30, 2014): this should return a smart pointer not an instance
#if 0
template <int dim>
CartesianGrid<dim> build_cartesian_grid_union(
    const CartesianGrid<dim> &grid_1,
    const CartesianGrid<dim> &grid_2,
    vector<Index> &map_elem_grid_union_to_elem_grid_1,
    vector<Index> &map_elem_grid_union_to_elem_grid_2)
{
    //---------------------------------------------------------
    // checks that the grid are on the same domain
    Assert(grid_1.get_bounding_box() == grid_2.get_bounding_box(),
           ExcMessage("Grids on different domains."));
    //---------------------------------------------------------


    //---------------------------------------------------------
    // getting the coordinates from the two grids and building the grid union
    array<vector<Real>,dim> knots_union;
    for (int i = 0 ; i < dim ; ++i)
    {
        const auto &coords_grid_1 = grid_1.get_knot_coordinates(i);
        const auto &coords_grid_2 = grid_2.get_knot_coordinates(i);

        // here we remove the duplicates (if any)
        set<Real> coords_unique(coords_grid_1.begin(),coords_grid_1.end());
        std::copy(coords_grid_2.begin(), coords_grid_2.end(),
                  std::inserter(coords_unique, coords_unique.end()));

        std::copy(coords_unique.begin(),coords_unique.end(),
                  std::back_inserter(knots_union[i]));
    }
    CartesianGrid<dim> grid_union(knots_union);
    //---------------------------------------------------------


    //---------------------------------------------------------
    map_elem_grid_union_to_elem_grid_1 =
            build_map_elements_between_cartesian_grids(grid_union,grid_1);
    map_elem_grid_union_to_elem_grid_2 =
            build_map_elements_between_cartesian_grids(grid_union,grid_2);
    //---------------------------------------------------------

    return grid_union;
}
#endif

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/cartesian_grid.inst>
