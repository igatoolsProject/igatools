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

#include <igatools/geometry/grid.h>
#include <igatools/geometry/grid_tools.h>
#include <igatools/base/exceptions.h>
#include <igatools/base/array_utils.h>
#include <igatools/utils/vector_tools.h>
#include <igatools/utils/multi_array_utils.h>
#include <igatools/utils/unique_id_generator.h>
#include <igatools/utils/tensor_range.h>
#include <algorithm>

using std::endl;
using std::shared_ptr;
using std::make_shared;

IGA_NAMESPACE_OPEN

namespace
{
/**
 * Given the boundaries of a dim_-dimensional box, this function
 * computes and returns a vector of knot vectors uniformly
 * distributed with the required numbers of knots.
 */
template <int dim_>
CartesianProductArray<Real,dim_>
filled_progression(const BBox<dim_> &end_points, const TensorSize<dim_> &n_knots)
{
  CartesianProductArray<Real,dim_> knot_coordinates(n_knots);

  SafeSTLVector<Real> knots_1d;
  for (const int i : UnitElement<dim_>::active_directions)
  {
    const Size n_i = n_knots[i];
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
Grid<dim_>::
Grid(const Size n)
  :
  Grid(TensorSize<dim_>(n))
{}



template<int dim_>
auto
Grid<dim_>::
create(const Size n) -> shared_ptr<self_t>
{
  return shared_ptr<self_t>(new self_t(n));
}

template<int dim_>
auto
Grid<dim_>::
const_create(const Index n_knots) -> shared_ptr<const self_t>
{
  return create(n_knots);
}




template<int dim_>
Grid<dim_>::
Grid(const TensorSize<dim_> &n)
  :
  Grid(BBox<dim_>(), n)
{}



template<int dim_>
auto
Grid<dim_>::
create(const TensorSize<dim_> &n) -> shared_ptr<self_t>
{
  return shared_ptr<self_t>(new self_t(n));
}

template<int dim_>
auto
Grid<dim_>::
const_create(const TensorSize<dim_> &n_knots) -> shared_ptr<const self_t>
{
  return create(n_knots);
}



template<int dim_>
Grid<dim_>::
Grid(const BBox<dim_> &end_points, const Size n_knots)
  :
  Grid(end_points, TensorSize<dim_>(n_knots))
{}



template<int dim_>
auto
Grid<dim_>::
create(const BBox<dim_> &end_points, const Size n_knots) -> shared_ptr<self_t>
{
  return shared_ptr<self_t>(new self_t(end_points, n_knots));
}

template<int dim_>
auto
Grid<dim_>::
const_create(const BBox<dim_> &bbox, const Size n_knots) -> shared_ptr<const self_t>
{
  return create(bbox, n_knots);
}


template<int dim_>
Grid<dim_>::
Grid(const BBox<dim_> &end_points,
     const TensorSize<dim_> &n)
  :
  Grid(filled_progression<dim_>(end_points, n))
{}



template<int dim_>
auto
Grid<dim_>::
create(const BBox<dim_> &end_points,
       const TensorSize<dim_> &n) -> shared_ptr<self_t>
{
  return shared_ptr<self_t>(new self_t(end_points, n));
}

template<int dim_>
auto
Grid<dim_>::
const_create(const BBox<dim_> &bbox,
             const TensorSize<dim_> &n_knots) -> shared_ptr<const self_t>
{
  return create(bbox, n_knots);
}


template<int dim_>
Grid<dim_>::
Grid(const KnotCoordinates &knot_coordinates)
  :
//  TensorSizedContainer<dim_>(TensorSize<dim_>(knot_coordinates.tensor_size()-1)),
  knot_coordinates_(knot_coordinates),
  boundary_id_(0),
  object_id_(UniqueIdGenerator::get_unique_id()),
  elems_size_(TensorSize<dim_>(knot_coordinates.tensor_size()-1))
{
  elem_properties_.add_property(ElementProperties::active);
  elem_properties_[ElementProperties::active] =
    el_tensor_range<dim>(TensorIndex<dim>(), get_num_intervals());
//    elem_properties_.set_ids_property_status(ElementProperties::active,
//                                             el_tensor_range<dim>(TensorIndex<dim>(), get_num_intervals()),
//                                             true);

#ifndef NDEBUG
  for (const int i : UnitElement<dim_>::active_directions)
  {
    const auto &knots_i = knot_coordinates.get_data_direction(i);
    // checks that we have at least two knot values (i.e. one knot span) in
    // each coordinate direction
    AssertThrow(knots_i.size() > 1, ExcLowerRange(knots_i.size(), 2));

    // check if the array is sorted and does not contains duplicates
    SafeSTLVector<Real> vec = knots_i ;
    std::sort(vec.begin(), vec.end());
    vec.erase(unique(vec.begin(), vec.end()), vec.end());
    AssertThrow(knots_i == vec,
                ExcMessage("The knot coordinate vector is not sorted and/or contains duplicates"));
  }
#endif
}



template<int dim_>
auto
Grid<dim_>::
create(const KnotCoordinates &knot_coordinates) -> shared_ptr<self_t>
{
  return shared_ptr<self_t>(new self_t(knot_coordinates));
}

template<int dim_>
auto
Grid<dim_>::
const_create(const KnotCoordinates &knots) -> shared_ptr<const self_t>
{
  return create(knots);
}


template<int dim_>
Grid<dim_>::
Grid(const SafeSTLArray<SafeSTLVector<Real>,dim_> &knot_coordinates)
  :
  self_t(CartesianProductArray<Real,dim_>(knot_coordinates))
{}



template<int dim_>
auto
Grid<dim_>::
create(const SafeSTLArray<SafeSTLVector<Real>,dim_> &knot_coordinates)
-> shared_ptr<self_t>
{
  return shared_ptr<self_t>(new self_t(knot_coordinates));
}

template<int dim_>
auto
Grid<dim_>::
const_create(const SafeSTLArray<SafeSTLVector<Real>,dim_> &knots)
-> shared_ptr<const self_t>
{
  return create(knots);
}



template<int dim_>
auto
Grid<dim_>::
create(const self_t &grid) -> std::shared_ptr<self_t>
{
  return shared_ptr<self_t>(new self_t(grid));
}

template<int dim_>
auto
Grid<dim_>::
const_create(const self_t &grid) -> shared_ptr<const self_t>
{
  return create(grid);
}


template<int dim_>
Grid<dim_>::
Grid(const self_t &grid)
  :
//  TensorSizedContainer<dim_>(grid),
  knot_coordinates_(grid.knot_coordinates_),
  boundary_id_(grid.boundary_id_),
  elem_properties_(grid.elem_properties_),
  object_id_(UniqueIdGenerator::get_unique_id()),
  elems_size_(grid.elems_size_)
{}



template <int dim_>
bool
Grid<dim_>::
operator==(const Grid<dim_> &grid) const
{
  bool same_knots_coordinates = true;
  for (const auto i : UnitElement<dim_>::active_directions)
  {
    const auto &knots_a =  this->knot_coordinates_.get_data_direction(i);
    const auto &knots_b =   grid.knot_coordinates_.get_data_direction(i);

    same_knots_coordinates = same_knots_coordinates && (knots_a == knots_b);
  }
  return same_knots_coordinates;
}



template<int dim_>
Index Grid<dim_>::get_object_id() const
{
  return object_id_;
}



template<int dim_>
SafeSTLVector< Real > const &
Grid<dim_>::get_knot_coordinates(const int i) const
{
  return knot_coordinates_.get_data_direction(i);
}



template<int dim_>
auto
Grid<dim_>::
get_knot_coordinates() const -> CartesianProductArray<Real,dim_> const &
{
  return knot_coordinates_;
}



template<int dim_>
auto
Grid<dim_>::get_element_lengths() const -> KnotCoordinates
{
  auto const &size = get_num_intervals();
  KnotCoordinates length(size);

  for (auto &i : UnitElement<dim_>::active_directions)
  {
    const auto &knots_i = knot_coordinates_.get_data_direction(i);
    const auto n_elem = size[i];
    for (int j = 0 ; j < n_elem ; ++j)
      length.entry(i,j) = knots_i[j+1] - knots_i[j];
  }
  return length;
}



template<int dim_>
void
Grid<dim_>::
add_property(const PropId &property)
{
  elem_properties_.add_property(property);
}



template<int dim_>
auto
Grid<dim_>::
get_elements_with_property(const PropId &prop) const
-> const List &
{
  return elem_properties_[prop];
}



template<int dim_>
auto
Grid<dim_>::
get_elements_with_property(const PropId &prop)
-> List &
{
  return elem_properties_[prop];
}



template<int dim_>
auto
Grid<dim_>::create_cache_handler() const
-> std::unique_ptr<ElementHandler>
{
  auto handler = std::make_unique<ElementHandler>(this->shared_from_this());
  Assert(handler != nullptr,ExcNullPtr());

  return handler;
}

template<int dim_>
auto
Grid<dim_>::
create_element(const ListIt &index, const PropId &prop) const
-> std::unique_ptr<ElementAccessor>
{
  using Elem = ElementAccessor;
  auto elem = std::make_unique<Elem>(this->shared_from_this(), index, prop);
  Assert(elem != nullptr,ExcNullPtr());

  return elem;
}







template<int dim_>
auto
Grid<dim_>::
begin(const PropId &prop) const -> ElementIterator
{
  return this->cbegin(prop);
}



template<int dim_>
auto
Grid<dim_>::
end(const PropId &prop) const -> ElementIterator
{
  return this->cend(prop);
}



template<int dim_>
auto
Grid<dim_>::
cbegin(const PropId &prop) const -> ElementIterator
{
  return ElementIterator(
           this->create_element(elem_properties_[prop].begin(),prop));
}



template<int dim_>
auto
Grid<dim_>::
cend(const PropId &prop) const -> ElementIterator
{
  return ElementIterator(
           this->create_element(elem_properties_[prop].end(),prop));
}


//
//template<int dim_>
//Index
//Grid<dim_>::
//tensor_to_flat(const TensorIndex<dim_> &tensor_index) const
//{
//    return TensorSizedContainer<dim_>::tensor_to_flat(tensor_index);
//}
//






template<int dim_>
void
Grid<dim_>::
set_boundary_id(const int face, const boundary_id id)
{
  Assert(face < UnitElement<dim_>::n_faces,
         ExcIndexRange(face,0, UnitElement<dim_>::n_faces));
  boundary_id_[face] = id;
}



template<int dim_>
boundary_id
Grid<dim_>::
get_boundary_id(const int face) const
{
  Assert(face < UnitElement<dim_>::n_faces,
         ExcIndexRange(face,0, UnitElement<dim_>::n_faces));
  return (boundary_id_[face]);
}



template<int dim_>
template<int sdim>
auto
Grid<dim_>::
get_boundary_normals(const int s_id) const -> BoundaryNormal<sdim>
{
  auto all_elems = UnitElement<dim_>::all_elems;
  auto element = std::get<sdim>(all_elems)[s_id];

  BoundaryNormal<sdim> normals;
  for (int i=0; i<dim_-sdim; ++i)
  {
    auto val = 2*element.constant_values[i]-1;
    normals[i][element.constant_directions[i]] = val;
  }

  return normals;
}



template<int dim_>
Size
Grid<dim_>::
get_num_elements(const PropId &prop) const
{
  return elem_properties_[prop].size();
}



template<int dim_>
Size
Grid<dim_>::
get_num_all_elems() const
{
  return elems_size_.flat_size();
}



template<int dim_>
auto
Grid<dim_>::
get_num_intervals() const -> TensorSize<dim_>
{
  return elems_size_.tensor_size();
}



template<int dim_>
auto
Grid<dim_>::
get_num_knots_dim() const -> TensorSize<dim_>
{
  return knot_coordinates_.tensor_size();
}


#ifdef MESH_REFINEMENT

template<int dim_>
auto
Grid<dim_>::
get_grid_pre_refinement() const -> shared_ptr<const self_t>
{
  return grid_pre_refinement_;
}


template <int dim_>
void
Grid<dim_>::
refine_directions(
  const SafeSTLArray<bool,dim_> &refinement_directions,
  const SafeSTLArray<Size,dim_> &n_subdivisions)
{
  //-------------------------------------------------------------
  SafeSTLArray<SafeSTLVector<Real>,dim_> knots_to_insert;
  for (const auto dir : UnitElement<dim_>::active_directions)
  {
    if (refinement_directions[dir])
    {
      const int n_sub_interv = n_subdivisions[dir];
      Assert(n_sub_interv > 0,ExcLowerRange(n_sub_interv,1));

      const auto &knots_old = this->get_knot_coordinates(dir);

      const Size n_knots_old = knots_old.size();
      for (Index i = 0 ; i < n_knots_old - 1 ; ++i)
      {
        const Real h = (knots_old[i+1] - knots_old[i]) / n_sub_interv;

        for (Index j = 1 ; j < n_sub_interv ; ++j)
          knots_to_insert[dir].emplace_back(knots_old[i] + j * h);
      }
    }
  }
  this->insert_knots(knots_to_insert);
  //-------------------------------------------------------------
}



template <int dim_>
void
Grid<dim_>::
refine_direction(const int direction_id, const Size n_subdivisions)
{
  Assert(direction_id >= 0 && direction_id < dim_,
         ExcIndexRange(direction_id, 0, dim_));

  SafeSTLArray<bool,dim_> refinement_directions(false);
  refinement_directions[direction_id] = true;

  SafeSTLArray<Size,dim_> n_subdiv;
  n_subdiv[direction_id] = n_subdivisions;

  this->refine_directions(refinement_directions,n_subdiv);
}



template <int dim_>
void
Grid<dim_>::
refine(const Size n_subdivisions)
{
  Assert(n_subdivisions >= 2, ExcLowerRange(n_subdivisions,2));

  this->refine_directions(
    SafeSTLArray<bool,dim_>(true),
    SafeSTLArray<Size,dim_>(n_subdivisions));
}



template <int dim_>
boost::signals2::connection
Grid<dim_>::
connect_insert_knots(const SignalInsertKnotsSlot &subscriber)
{
  return insert_knots_signals_.connect(subscriber);
}


template <int dim_>
void
Grid<dim_>::
insert_knots(SafeSTLArray<SafeSTLVector<Real>,dim_> &knots_to_insert)
{
  //----------------------------------------------------------------------------------
  // make a copy of the grid before the refinement
  grid_pre_refinement_ = make_shared<self_t>(*this);
  //----------------------------------------------------------------------------------


  //----------------------------------------------------------------------------------
  // inserts the knots into the current grid --- begin
  for (const auto dir : UnitElement<dim_>::active_directions)
  {
    std::set<Real> new_coords_no_duplicates(knots_to_insert[dir].begin(),knots_to_insert[dir].end());

    const auto &old_coords = knot_coordinates_.get_data_direction(dir);
    new_coords_no_duplicates.insert(old_coords.begin(),old_coords.end());

    knot_coordinates_.copy_data_direction(
      dir,
      SafeSTLVector<Real>(new_coords_no_duplicates.begin(),
                          new_coords_no_duplicates.end()));
  }
  TensorSizedContainer<dim_>::reset_size(knot_coordinates_.tensor_size()-1);
  // inserts the knots into the current grid --- end
  //----------------------------------------------------------------------------------

//Assert(false,ExcNotImplemented());
  //----------------------------------------------------------------------------------
  // transferring the element properties from the old grid to the new grid --- begin

  elem_properties_[ElementProperties::active] =
    el_tensor_range<dim>(TensorIndex<dim>(), get_num_intervals());

  const auto fine_to_coarse_elems_id = grid_tools::build_map_elements_id_between_grids(
                                         *this,*grid_pre_refinement_);

  std::map<IndexType,std::set<IndexType>> coarse_to_fine_elems_id;
  for (const auto &fine_to_coarse_elem_id : fine_to_coarse_elems_id)
  {
    const auto &elem_id_fine   = fine_to_coarse_elem_id.first;
    const auto &elem_id_coarse = fine_to_coarse_elem_id.second;
    coarse_to_fine_elems_id[elem_id_coarse].insert(elem_id_fine);
  }

  for (auto &elems_property : elem_properties_)
  {
    const auto &property_name = elems_property.first;

    if (property_name != ElementProperties::active)
    {
      auto &elems_id_fine_with_property = elems_property.second;
      elems_id_fine_with_property.clear();

      const auto &elems_id_coarse_with_property =
        grid_pre_refinement_->get_element_property(property_name);

      for (const auto &elem_id_coarse : elems_id_coarse_with_property)
      {
        const auto &elems_id_fine = coarse_to_fine_elems_id[elem_id_coarse];
        elems_id_fine_with_property.insert(elems_id_fine.begin(),elems_id_fine.end());
      }
    }
  }

#if 0
  auto coarse_elem = grid_pre_refinement_->begin();
  for (const auto &fine_coarse_elem_id : fine_to_coarse_grid)
  {
    const auto   &fine_elem_id = fine_coarse_elem_id.first;
    const auto &coarse_elem_id = fine_coarse_elem_id.second;

    coarse_elem->move_to(coarse_elem_id);
//        const auto fine_elem_id = fine_elem->get_flat_index();

    const auto old_elem_properties = coarse_elem->get_defined_properties();

    for (const auto &property : old_elem_properties)
      this->set_element_property_status(property,fine_elem_id,true);
  }
  // transferring the element properties from the old grid to the new grid --- end
  //----------------------------------------------------------------------------------
#endif

  //----------------------------------------------------------------------------------
  // refining the objects that's are attached to the Grid
  // (i.e. that are defined using this Grid object)
  this->insert_knots_signals_(knots_to_insert,*grid_pre_refinement_);
  //----------------------------------------------------------------------------------
}



template <int dim_>
bool
Grid<dim_>::
same_knots_or_refinement_of(const Grid<dim_> &grid_to_compare_with) const
{
  bool is_refinement = true;
  for (auto dir : UnitElement<dim_>::active_directions)
  {
    const auto &knots_coarse = grid_to_compare_with.get_knot_coordinates(dir);
    const auto &knots_fine   = this->get_knot_coordinates(dir);

    //look if there is any value in knots_coarse not in knots_fine
    if (std::any_of(
          knots_coarse.begin(),
          knots_coarse.end(),
          [&knots_fine](const Real &val)
  {
    return !std::binary_search(knots_fine.begin(),knots_fine.end(),val);
    }))
    {
      is_refinement = false;
      break;
    }
  }

  return is_refinement;
}
#endif // MESH_REFINEMENT



template <int dim_>
void
Grid<dim_>::
print_info(LogStream &out) const
{
  out << "Number of intervals per direction: " << elems_size_.tensor_size() << endl;

  out.begin_item("Knot coordinates:");
  knot_coordinates_.print_info(out);
  out.end_item();


  //-------------------------------------------------------------
  if (!elem_properties_.empty())
  {
    out.begin_item("Element properties:");
    elem_properties_.print_info(out);
    out.end_item();
  }
  //-------------------------------------------------------
}



template <int dim_>
template<int sdim>
auto
Grid<dim_>::
get_sub_grid(const int s_id, SubGridMap<sdim> &elem_map) const
-> shared_ptr<Grid<sdim>>
{
  auto &s_elem = UnitElement<dim_>::template get_elem<sdim>(s_id);
  const auto active_dirs = TensorIndex<sdim>(s_elem.active_directions);
  auto sub_knots = knot_coordinates_.template get_sub_product<sdim>(active_dirs);
  auto sub_grid = Grid<sdim>::create(sub_knots);

  IndexType grid_index;
  const int n_dir = s_elem.constant_directions.size();
  for (int j = 0 ; j < n_dir ; ++j)
  {
    auto dir = s_elem.constant_directions[j];
    auto val = s_elem.constant_values[j];
    grid_index[dir] = val == 0 ? 0 : (knot_coordinates_.tensor_size()[dir]-2);
  }

  elem_map.clear();
  for (const auto &s_elem : *sub_grid)
  {
    const auto s_index = s_elem.get_index();
    for (int j = 0 ; j < sdim ; ++j)
      grid_index[active_dirs[j]] = s_index[j];

    elem_map.emplace(s_index,grid_index);
  }
  return sub_grid;
}



template <int dim_>
auto
Grid<dim_>::
get_bounding_box() const -> BBox<dim_>
{
  BBox<dim_> bounding_box;

  for (const auto i : UnitElement<dim_>::active_directions)
  {
    bounding_box[i][0] = knot_coordinates_.get_data_direction(i).front();
    bounding_box[i][1] = knot_coordinates_.get_data_direction(i).back();
  }

  return bounding_box;
}


template <int dim_>
bool
Grid<dim_>::
test_if_point_on_internal_knots_line(const Points<dim_> &point) const
{
  bool res = false;
  for (const auto i : UnitElement<dim_>::active_directions)
  {
    const auto &knots = knot_coordinates_.get_data_direction(i);
    if (std::binary_search(knots.begin()+1,knots.end()-1,point[i]))
    {
      res = true;
      break;
    }
  }
  return res;
}

template <int dim_>
auto
Grid<dim_>::
find_elements_id_of_points(const ValueVector<Points<dim_>> &points) const
-> std::map<IndexType, SafeSTLVector<int> >
{
  std::map<IndexType, SafeSTLVector<int> > res;

  const int n_points = points.size();
  for (int k = 0 ; k < n_points ; ++k)
  {
    const auto &point = points[k];

    Assert(!test_if_point_on_internal_knots_line(point),
    ExcMessage("The " + std::to_string(k) + "-th point is on an intenal knots line."));

    TensorIndex<dim_> elem_t_id;
    for (const auto i : UnitElement<dim_>::active_directions)
    {
      const auto &knots = knot_coordinates_.get_data_direction(i);

      Assert(point[i] >= knots.front() && point[i] <= knots.back(),
      ExcMessage("The point " + std::to_string(k) +
      " is not in the interval spanned by the knots along the direction " +
      std::to_string(i)));

      //find the index j in the knots for which knots[j] <= point[i]
      const auto low = std::lower_bound(knots.begin(),knots.end(),point[i]);


      const Index j = low - knots.begin();

      elem_t_id[i] = (j>0) ? j-1 : 0;
    }
    /*
            auto ans = res.emplace(
                  elem_t_id,
                           SafeSTLVector<int>(1,k));

            if (!ans.second)
                (ans.first)->second.push_back(k);
                //*/
    res[elem_t_id].push_back(k);

  }
  return res;
}


#if 0
template <int dim_>
auto
Grid<dim_>::
find_elements_id_of_points(const ValueVector<Points<dim_>> &points) const -> SafeSTLVector<IndexType>
{
//    Assert(false,ExcNotImplemented());
//    Assert(false,ExcMessage("This function is not tested at all!"));
  SafeSTLVector<IndexType> elements_id;

//    SafeSTLArray<SafeSTLVector<int>,dim> ids;

//    TensorSize<dim_> n_elems_dir;

  IndexType elem_t_id;
  for (const auto &point : points)
  {
    for (const auto dir : UnitElement<dim_>::active_directions)
    {
      const auto &knots = knot_coordinates_.get_data_direction(dir);

      const auto &p = point[dir];

      Assert(p >= knots.front() && p <= knots.back(),
      ExcMessage("The point coordinate p[" + std::to_string(dir) + "]= "
      + std::to_string(p) +
      " is not in the interval spanned by the knots along the direction " +
      std::to_string(dir)));

      //find the index j in the knots for which knots[j] <= point[dir]
      const auto low = std::lower_bound(knots.begin(),knots.end(),p);
      const Index j = low - knots.begin();

      elem_t_id[dir] = (j>0) ? j-1 : 0;
      /*
                  if (j > 0)
                  {
                      ids[dir].push_back(j-1);
                      if (p == knots[j])
                          ids[dir].push_back(j);
                  }
                  else
                  {
                      ids[dir].push_back(0);
                  }
                  n_elems_dir[dir] = ids[dir].size();
                  //*/
    }
  }

  const auto n_elems = n_elems_dir.flat_size();
  const auto w_elems_dir = MultiArrayUtils<dim_>::compute_weight(n_elems_dir);
  for (int elem = 0 ; elem < n_elems ; ++elem)
  {
    const auto t_id = MultiArrayUtils<dim_>::flat_to_tensor_index(elem,w_elems_dir);
    TensorIndex<dim_> elem_t_id;
    for (const auto dir : UnitElement<dim_>::active_directions)
      elem_t_id[dir] = ids[dir][t_id[dir]];

    elements_id.emplace_back(elem_t_id);
  }
  return elements_id;
}
#endif


template <int dim_>
SafeSTLVector<Index>
Grid<dim_>::
get_sub_elements_id(const TensorSize<dim_> &n_sub_elems, const Index elem_id) const
{
  const auto coarse_elem_tensor_id = elems_size_.flat_to_tensor(elem_id);

  const auto weight_sub_elems = MultiArrayUtils<dim_>::compute_weight(n_sub_elems);

  const TensorSize<dim_> n_elems_coarse_grid = this->get_num_intervals();
  TensorSize<dim_> n_elems_fine_grid;
  TensorIndex<dim_> first_fine_elem_tensor_id;
  for (const auto &i : UnitElement<dim_>::active_directions)
  {
    Assert(n_sub_elems[i] > 0 ,ExcLowerRange(n_sub_elems[i],1));

    n_elems_fine_grid[i] = n_elems_coarse_grid[i] * n_sub_elems[i];

    first_fine_elem_tensor_id[i] = n_sub_elems[i] * coarse_elem_tensor_id[i];
  }
  const auto weight_fine_grid = MultiArrayUtils<dim_>::compute_weight(n_elems_fine_grid);


  const int n_sub_elems_total = n_sub_elems.flat_size();
  SafeSTLVector<Index> sub_elems_id(n_sub_elems_total);
  TensorIndex<dim_> fine_elem_tensor_id;
  for (int sub_elem = 0 ; sub_elem < n_sub_elems_total ; ++sub_elem)
  {
    const auto sub_elem_tensor_id =
      MultiArrayUtils<dim_>::flat_to_tensor_index(sub_elem, weight_sub_elems);

    for (const auto &i : UnitElement<dim_>::active_directions)
      fine_elem_tensor_id[i] = first_fine_elem_tensor_id[i] + sub_elem_tensor_id[i];

    sub_elems_id[sub_elem] =
      MultiArrayUtils<dim_>::tensor_to_flat_index(fine_elem_tensor_id, weight_fine_grid);
  } // end loop sub_elem_flat_id

  return sub_elems_id;
}



template <int dim_>
void
Grid<dim_>::
set_property_status_elem(const PropId &property,
                         const IndexType &elem_id,
                         const bool property_status)
{
  Assert(dim_ > 0,ExcMessage("Setting a property for Grid<dim_> with dim_==0 has no meaning."));
  elem_properties_.set_property_status_for_id(property,elem_id,property_status);
}



#if 0

template <int dim_>
void
Grid<dim_>::
set_all_elements_property_status(const PropId &property,
                                 const bool status)
{
  for (const auto &elem : (*this))
    this->set_element_property_status(property,elem.get_flat_index(),status);
}
#endif



template <int dim_>
bool
Grid<dim_>::
element_has_property(const IndexType elem_id,
                     const PropId &prop) const
{
  const auto &elem_list = this->get_elements_with_property(prop);
  return std::binary_search(elem_list.begin(),
                            elem_list.end(),
                            elem_id);
}


#ifdef SERIALIZATION

template <int dim_>
template<class Archive>
void
Grid<dim_>::
serialize(Archive &ar)
  {
    ar &make_nvp("knot_coordinates_",knot_coordinates_);
    ar &make_nvp("boundary_id_",boundary_id_);
    ar &make_nvp("properties_elements_id_",elem_properties_);
    ar &make_nvp("object_id_",object_id_);
    ar &make_nvp("elems_size_",elems_size_);
#ifdef MESH_REFINEMENT
    ar &make_nvp("grid_pre_refinement_",grid_pre_refinement_);
#endif
  }
#endif // SERIALIZATION



IGA_NAMESPACE_CLOSE

#include <igatools/geometry/grid.inst>

