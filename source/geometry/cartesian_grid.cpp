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

#include <igatools/geometry/cartesian_grid.h>
#include <igatools/base/exceptions.h>
#include <igatools/base/array_utils.h>
#include <igatools/utils/vector_tools.h>
#include <igatools/utils/multi_array_utils.h>

#include <algorithm>
using std::endl;

using std::shared_ptr;
using std::make_shared;

IGA_NAMESPACE_OPEN

namespace
{

/**
 * Given the boundaries of a dim-dimensional box, it
 * computes and returns a vector of knot vectors uniformly
 * distrubuted with the required numbers of knots.
 */
template <int dim>
CartesianProductArray<Real,dim>
filled_progression(const BBox<dim> &end_points, const TensorSize<dim> &n_knots)
{
    CartesianProductArray<Real,dim> knot_coordinates(n_knots);

    vector<Real> knots_1d;
    for (int i = 0; i < dim; ++i)
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



template <int dim_>
constexpr int CartesianGrid<dim_>::dim;



template<int dim_>
constexpr std::array<Size, dim_> CartesianGrid<dim_>::dims;



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
    CartesianGrid(filled_array<std::array<Real,2>,dim>(std::array<Real,2> {{0,1}}), n, kind)
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
    TensorSizedContainer<dim_>(TensorSize<dim_>(knot_coordinates.tensor_size()-1)),
    kind_(kind),
    boundary_id_(filled_array<int,UnitElement<dim>::n_faces>(0)),
    knot_coordinates_(knot_coordinates)
{
    const int n_elems = this->flat_size();

    auto &active_elements = properties_elements_id_[ElementProperty::active];
    auto &marked_elements = properties_elements_id_[ElementProperty::marked];
    auto &influence_elements = properties_elements_id_[ElementProperty::influence];

    for (int elem_id = 0 ; elem_id < n_elems ; ++elem_id)
    {
        active_elements.insert(elem_id);
        marked_elements.insert(elem_id);
        influence_elements.insert(elem_id);
    }


#ifndef NDEBUG
    for (int i = 0; i < dim; i++)
    {
        const auto &knots_i = knot_coordinates.get_data_direction(i);
        // checks that we have at least two knot values (i.e. one knot span) in
        // each coordinate direction
        AssertThrow(knots_i.size() > 1, ExcLowerRange(knots_i.size(), 2));

        // check if the array is sorted and does not contains duplicates
        vector<Real> vec = knots_i ;
        std::sort(vec.begin(), vec.end());
        vec.erase(unique(vec.begin(), vec.end()), vec.end());
        AssertThrow(knots_i == vec,
                    ExcMessage("The knot coordinate vector is not sorted and/or contains duplicates"));

    }
#endif
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
CartesianGrid(const std::array<vector<Real>,dim> &knot_coordinates)
    :
    self_t(CartesianProductArray<Real,dim>(knot_coordinates),
           Kind::direction_uniform)
{}



template<int dim_>
auto
CartesianGrid<dim_>::
create(const std::array<vector<Real>,dim> &knot_coordinates) -> shared_ptr<self_t>
{
    return shared_ptr< self_t >(new self_t(knot_coordinates));
}

template<int dim_>
auto
CartesianGrid<dim_>::
create(const self_t &grid) -> std::shared_ptr<self_t>
{
    return shared_ptr<self_t>(new self_t(grid));
}



template<int dim_>
CartesianGrid<dim_>::
CartesianGrid(const self_t &grid)
    :
    TensorSizedContainer<dim_>(grid),
    kind_(grid.kind_),
    boundary_id_(grid.boundary_id_),
    knot_coordinates_(grid.knot_coordinates_),
    properties_elements_id_(grid.properties_elements_id_)
{}



template<int dim_>
vector< Real > const &
CartesianGrid<dim_>::
get_knot_coordinates(const int i) const
{
    return knot_coordinates_.get_data_direction(i);
}



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
    auto const &size = get_num_intervals();
    KnotCoordinates length(size);

    for (auto &i : dims)
    {
        const auto &knots_i = knot_coordinates_.get_data_direction(i);
        const auto n_elem = size[i];
        for (int j = 0 ; j < n_elem ; ++j)
            length.entry(i,j) = knots_i[j+1] - knots_i[j];
    }
    return length;
}


template<int dim_>
std::set<Index> &
CartesianGrid<dim_>::
get_elements_id_same_property(const ElementProperty &property)
{
    return properties_elements_id_[property];
}

template<int dim_>
const std::set<Index> &
CartesianGrid<dim_>::
get_elements_id_same_property(const ElementProperty &property) const
{
    return properties_elements_id_.at(property);
}

template<int dim_>
auto
CartesianGrid<dim_>::
begin() -> ElementIterator
{
    const auto &active_elements = get_elements_id_same_property(ElementProperty::active);
    const auto active_begin = active_elements.begin();
    const auto active_end   = active_elements.end();

    int id_first_active_elem;
    if (active_begin != active_end)
        id_first_active_elem = *active_begin;
    else
        id_first_active_elem = IteratorState::pass_the_end;

    return ElementIterator(this->shared_from_this(), id_first_active_elem);
}

template<int dim_>
auto
CartesianGrid<dim_>::
create_element(const Index flat_index) const -> std::shared_ptr<ElementAccessor>
{
    using Elem = CartesianGridElement<dim_>;
    auto elem = shared_ptr<Elem>(new Elem(this->shared_from_this(),flat_index));
    Assert(elem != nullptr,ExcNullPtr());

    return elem;
}


template<int dim_>
auto
CartesianGrid<dim_>::
end() -> ElementIterator
{
    return ElementIterator(this->create_element(IteratorState::pass_the_end));
}

template<int dim_>
auto
CartesianGrid<dim_>::
last() -> ElementIterator
{
    const auto &active_elements = get_elements_id_same_property(ElementProperty::active);
    const auto active_begin = active_elements.begin();
    auto active_end   = active_elements.end();
    int id_last_active_elem;
    if (active_begin != active_end)
        id_last_active_elem = *(--active_end);
    else
        id_last_active_elem = IteratorState::pass_the_end;

    return ElementIterator(this->shared_from_this(), id_last_active_elem);
}

template<int dim_>
auto
CartesianGrid<dim_>::
last() const -> ElementConstIterator
{
    return clast();
}

template<int dim_>
auto
CartesianGrid<dim_>::
clast() const -> ElementConstIterator
{
    const auto &active_elements = get_elements_id_same_property(ElementProperty::active);
    const auto active_begin = active_elements.begin();
    auto active_end   = active_elements.end();

    int id_last_active_elem;
    if (active_begin != active_end)
        id_last_active_elem = *(--active_end);
    else
        id_last_active_elem = IteratorState::pass_the_end;

    return ElementConstIterator(this->shared_from_this(), id_last_active_elem);
}

template<int dim_>
auto
CartesianGrid<dim_>::
begin() const -> ElementConstIterator
{
    return this->cbegin();
}



template<int dim_>
auto
CartesianGrid<dim_>::
end() const -> ElementConstIterator
{
    return this->cend();
}

template<int dim_>
auto
CartesianGrid<dim_>::
cbegin() const -> ElementConstIterator
{
    const auto &active_elements = get_elements_id_same_property(ElementProperty::active);
    const auto active_begin = active_elements.begin();
    const auto active_end   = active_elements.end();

    int id_first_active_elem;
    if (active_begin != active_end)
        id_first_active_elem = *active_begin;
    else
        id_first_active_elem = IteratorState::pass_the_end;

    return ElementConstIterator(this->shared_from_this(), id_first_active_elem);
}



template<int dim_>
auto
CartesianGrid<dim_>::
cend() const -> ElementConstIterator
{
    return ElementConstIterator(this->create_element(IteratorState::pass_the_end));
}





template<int dim_>
Index
CartesianGrid<dim_>::
tensor_to_flat(const TensorIndex<dim_> &tensor_index) const
{
    return TensorSizedContainer<dim_>::tensor_to_flat(tensor_index);
}



template<int dim_>
TensorIndex<dim_>
CartesianGrid<dim_>::
flat_to_tensor(const Index flat_index) const
{
    return TensorSizedContainer<dim_>::flat_to_tensor(flat_index);
}



template<int dim_>
void
CartesianGrid<dim_>::
set_boundary_id(const int face, const boundary_id id)
{
    Assert(face < UnitElement<dim>::n_faces,
           ExcIndexRange(face,0, UnitElement<dim>::n_faces));
    boundary_id_[face] = id;
}



template<int dim_>
boundary_id
CartesianGrid<dim_>::
get_boundary_id(const int face) const
{
    Assert(face < UnitElement<dim>::n_faces,
           ExcIndexRange(face,0, UnitElement<dim>::n_faces));
    return (boundary_id_[face]);
}



template<int dim_>
template<int sub_dim>
auto
CartesianGrid<dim_>::
get_boundary_normals(const int s_id) const -> BoundaryNormal<sub_dim>
{
    auto all_elems = UnitElement<dim>::all_elems;
    auto element = std::get<sub_dim>(all_elems)[s_id];

    BoundaryNormal<sub_dim> normals;
    for (int i=0; i<dim-sub_dim; ++i)
    {
        auto val = 2*element.constant_values[i]-1;
        normals[i][element.constant_directions[i]] = val;
    }

    return normals;
}



template<int dim_>
Size
CartesianGrid<dim_>::
get_num_active_elems() const
{
    return this->get_elements_id_same_property(ElementProperty::active).size();
}

template<int dim_>
Size
CartesianGrid<dim_>::
get_num_influence_elems() const
{
    return this->get_elements_id_same_property(ElementProperty::influence).size();
}


template<int dim_>
Size
CartesianGrid<dim_>::
get_num_all_elems() const
{
    return this->flat_size();
}



// TODO (pauletti, Jul 30, 2014): this should not be necesary to specialize
template<>
Size
CartesianGrid<0>::
get_num_active_elems() const
{
    return 1;
}



template<int dim_>
auto
CartesianGrid<dim_>::
get_num_intervals() const -> TensorSize<dim>
{
    return this->tensor_size();
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
    const std::array<bool,dim> &refinement_directions,
    const std::array<Size,dim> &n_subdivisions)
{
    // make a copy of the grid before the refinement
    grid_pre_refinement_ = make_shared<const self_t>(self_t(*this));

    TensorSize<dim_> n_sub_elems;
    for (auto i : dims)
    {
        if (refinement_directions[i])
        {
            Assert(n_subdivisions[i] > 0,ExcLowerRange(n_subdivisions[i],1));
            n_sub_elems[i] = n_subdivisions[i];
            this->refine_knots_direction(i,n_sub_elems[i]);
        }
        else
        {
            n_sub_elems[i] = 1;
        }
    }
    TensorSizedContainer<dim_>::reset_size(knot_coordinates_.tensor_size()-1);

    //-------------------------------------------------------------
    for (auto &elements_same_property : properties_elements_id_)
    {
        auto &old_elem_ids = elements_same_property.second;

        std::set<Index> new_elem_ids;
        for (const auto &old_elem_id : old_elem_ids)
        {
            const auto new_elem_ids_to_add =
                grid_pre_refinement_->get_sub_elements_id(n_sub_elems,old_elem_id);

            new_elem_ids.insert(new_elem_ids_to_add.begin(),new_elem_ids_to_add.end());
        } // end loop old_elem_id
        old_elem_ids = std::move(new_elem_ids);
    } // end loop elements_same_property
    //-------------------------------------------------------------

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

    std::array<bool,dim> refinement_directions = filled_array<bool,dim>(false);
    refinement_directions[direction_id] = true;

    std::array<Size,dim> n_subdiv;
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
        const Real h = (knots_old[i_old+1] - knots_old[i_old]) / n_subdivisions;

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
    out << "Number of intervals per direction: " << this->tensor_size() << endl;

    out.begin_item("Knot coordinates:");
    knot_coordinates_.print_info(out);
    out.end_item();


    //-------------------------------------------------------
    out.begin_item("Active elements:");
    out << "Num.: " << get_num_active_elems() << endl;
    if (get_num_active_elems() > 0)
    {
        out << "   Flat IDs:";
        const auto &active_elements = get_elements_id_same_property(ElementProperty::active);
        for (const auto &elem_id : active_elements)
            out << " " << elem_id ;
    }
    out.end_item();
    //-------------------------------------------------------


    //-------------------------------------------------------
    out.begin_item("Influence elements:");
    out << "Num.: " << get_num_influence_elems() << endl;
    if (get_num_influence_elems() > 0)
    {
        out << "   Flat IDs:";
        const auto &influence_elements = get_elements_id_same_property(ElementProperty::influence);
        for (const auto &elem_id : influence_elements)
            out << " " << elem_id ;
    }
    out.end_item();
    //-------------------------------------------------------


}



template <int dim_>
template<int k>
auto
CartesianGrid<dim_>::
get_sub_grid(const int sub_elem_id, InterGridMap<k> &elem_map) const
-> std::shared_ptr<CartesianGrid<k>>
{
    auto &k_elem = UnitElement<dim>::template get_elem<k>(sub_elem_id);
    const auto active_dirs = TensorIndex<k>(k_elem.active_directions);
    auto sub_knots = knot_coordinates_.template get_sub_product<k>(active_dirs);
    auto sub_grid = CartesianGrid<k>::create(sub_knots);

    TensorIndex<dim> grid_index;
    const int n_dir = k_elem.constant_directions.size();
    for (int j = 0 ; j < n_dir ; ++j)
    {
        auto dir = k_elem.constant_directions[j];
        auto val = k_elem.constant_values[j];
        grid_index[dir] = val == 0 ? 0 : (knot_coordinates_.tensor_size()[dir]-2);
    }

    auto v_elem = begin();
    auto s_elem = sub_grid->begin();
    auto s_end  = sub_grid->end();
    for (; s_elem != s_end; ++s_elem)
    {
        auto s_index = s_elem.get_tensor_index();
        for (int j = 0 ; j < k ; ++j)
            grid_index[active_dirs[j]] = s_index[j];

        v_elem.move_to(grid_index);
        elem_map.emplace(s_elem, v_elem);
    }

    return sub_grid;
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
auto
CartesianGrid<dim_>::
find_elements_of_points(const ValueVector<Points<dim>> &points) const
-> std::map<ElementIterator, vector<int> >
{
    std::map<ElementIterator, vector<int> > res;

    const int n_points = points.size();
    for (int k=0; k<n_points; ++k)
    {
        const auto &point = points[k];
        TensorIndex<dim> elem_t_id;
        for (int i = 0 ; i < dim ; ++i)
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

        auto ans = res.emplace(
                       ElementIterator(this->shared_from_this(), this->tensor_to_flat(elem_t_id)),
                       vector<int>(1,k));

        if (!ans.second)
            (ans.first)->second.push_back(k);
    }
    return res;
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

template <int dim_>
vector<Index>
CartesianGrid<dim_>::
get_sub_elements_id(const TensorSize<dim_> &n_sub_elems, const Index elem_id) const
{
    const auto coarse_elem_tensor_id = this->flat_to_tensor(elem_id);

    const auto weight_sub_elems = MultiArrayUtils<dim>::compute_weight(n_sub_elems);

    const TensorSize<dim> n_elems_coarse_grid = this->get_num_intervals();
    TensorSize<dim> n_elems_fine_grid;
    TensorIndex<dim> first_fine_elem_tensor_id;
    for (const auto &i : dims)
    {
        Assert(n_sub_elems[i] > 0 ,ExcLowerRange(n_sub_elems[i],1));

        n_elems_fine_grid[i] = n_elems_coarse_grid[i] * n_sub_elems[i];

        first_fine_elem_tensor_id[i] = n_sub_elems[i] * coarse_elem_tensor_id[i];
    }
    const auto weight_fine_grid = MultiArrayUtils<dim>::compute_weight(n_elems_fine_grid);


    const int n_sub_elems_total = n_sub_elems.flat_size();
    vector<Index> sub_elems_id(n_sub_elems_total);
    TensorIndex<dim> fine_elem_tensor_id;
    for (int sub_elem = 0 ; sub_elem < n_sub_elems_total ; ++sub_elem)
    {
        const auto sub_elem_tensor_id =
            MultiArrayUtils<dim>::flat_to_tensor_index(sub_elem, weight_sub_elems);

        for (const auto &i : dims)
            fine_elem_tensor_id[i] = first_fine_elem_tensor_id[i] + sub_elem_tensor_id[i];

        sub_elems_id[sub_elem] =
            MultiArrayUtils<dim>::tensor_to_flat_index(fine_elem_tensor_id, weight_fine_grid);
    } // end loop sub_elem_flat_id

    return sub_elems_id;
}


template <int dim_>
void
CartesianGrid<dim_>::
set_element_property(const ElementProperty &property,
                     const Index elem_flat_id,
                     const bool property_status)
{
    Assert(dim_ > 0,ExcMessage("Setting a property for CartesianGrid<dim> with dim==0 has no meaning."));

    auto &elems_same_property = get_elements_id_same_property(property);
    if (property_status)
    {
        elems_same_property.insert(elem_flat_id);
    }
    else
    {
        Assert(!elems_same_property.empty(),ExcEmptyObject());
        elems_same_property.erase(elem_flat_id);
    }

}

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/cartesian_grid.inst>
