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
#include <igatools/utils/vector_tools.h>
#include <algorithm>

using std::endl;
using std::array;
using std::vector;
using std::shared_ptr;

IGA_NAMESPACE_OPEN

template<int dim_>
CartesianGrid<dim_>::
CartesianGrid(const Size n)
    :
    CartesianGrid(TensorSize<dim>(n))
{
    kind_ = Kind::uniform;
}



template<int dim_>
shared_ptr< CartesianGrid<dim_> >
CartesianGrid<dim_>::
create(const Size n)
{
    return (shared_ptr< CartesianGrid<dim_> >(new CartesianGrid<dim_>(n)));
}



template<int dim_>
CartesianGrid<dim_>::
CartesianGrid(const TensorSize<dim> &n)
    :
    CartesianGrid(filled_array<array<Real,2>, dim>(array<Real,2> {{0,1}}), n)
{
    kind_ = Kind::direction_uniform ;
}



template<int dim_>
shared_ptr< CartesianGrid<dim_> >
CartesianGrid<dim_>::
create(const TensorSize<dim> &n)
{
    return (shared_ptr< CartesianGrid<dim_> >(new CartesianGrid<dim_>(n)));
}



template<int dim_>
CartesianGrid<dim_>::
CartesianGrid(const BBox<dim> &end_points,
              const TensorSize<dim> &n)
    :
    kind_ {Kind::direction_uniform},
      boundary_id_(filled_array<int,UnitElement<dim>::faces_per_element>(0)),
      knot_coordinates_(n)
{
    vector<Real> knots_1d;
    for (int i = 0; i < dim; ++i)
    {
        const Size n_i = n(i);
        Assert(n_i > 1, ExcLowerRange(n_i,2));

        knots_1d.resize(n_i);

        const Real h=(end_points[i][1] - end_points[i][0]) /(n_i-1);

        knots_1d[0] = end_points[i][0];
        for (int j = 1; j < n_i; ++j)
            knots_1d[ j ] = knots_1d[ j-1 ] + h;

        knot_coordinates_.copy_data_direction(i,knots_1d);
    }
}



template<int dim_>
shared_ptr< CartesianGrid<dim_> >
CartesianGrid<dim_>::
create(const BBox<dim> &end_points,
       const TensorSize<dim> &n)
{
    return (shared_ptr< CartesianGrid<dim_> >(new CartesianGrid<dim_>(end_points, n)));
}



template<int dim_>
CartesianGrid<dim_>::
CartesianGrid(const CartesianProductArray<Real, dim> &knot_coordinates)
    :
    kind_ {Kind::non_uniform},
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
}



template<int dim_>
shared_ptr< CartesianGrid<dim_> >
CartesianGrid<dim_>::
create(const CartesianProductArray<Real,dim> &knot_coordinates)
{
    return (shared_ptr< CartesianGrid<dim_> >(
                new CartesianGrid<dim_>(knot_coordinates)));
}

template<int dim_>
CartesianGrid<dim_>::
CartesianGrid(const CartesianGrid<dim_> &grid)
    :
    kind_ {grid.kind_},
      boundary_id_(grid.boundary_id_),
      knot_coordinates_(grid.knot_coordinates_)
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
get_element_lengths() const -> CartesianProductArray<Real,dim>
{
    auto const &size = get_num_elements_dim();
    CartesianProductArray<Real, dim> length(size);
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
    return ElementIterator(
               const_cast<CartesianGrid<dim_>&>(*this),
               0);
}



template<int dim_>
auto
CartesianGrid<dim_>::
end() const -> ElementIterator
{
    return ElementIterator(const_cast<CartesianGrid<dim_>&>(*this),
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
get_face_normal(const int face_no) const -> Point<dim>
{
    Point<dim> normal;
    normal[UnitElement<dim>::face_to_component[face_no][0]] =
    UnitElement<dim>::face_normal_direction[face_no];

    return normal;
}



template<int dim_>
Size
CartesianGrid<dim_>::
get_num_elements() const
{
    const TensorSize<dim> n_elements_dim = this->get_num_elements_dim();

    return n_elements_dim.flat_size();
}



template<int dim_>
auto
CartesianGrid<dim_>::
get_num_elements_dim() const -> TensorSize<dim>
{
    // the number of elements in each coordinate direction
    // is equal to the number of knot coordinates - 1
    TensorSize<dim> num_elements_dim = knot_coordinates_.tensor_size();
    for (int i = 0 ; i < dim ; ++i)
        num_elements_dim(i) -= 1;

    return (num_elements_dim);
}



template<int dim_>
auto
CartesianGrid<dim_>::
get_num_knots_dim() const -> TensorSize<dim>
{
    return knot_coordinates_.tensor_size();
}

template <int dim_>
void
CartesianGrid<dim_>::
refine_directions(
    const array<bool,dim> &refinement_directions,
    const array<Size,dim> &n_subdivisions)
{
    // make a copy of the grid before the refinement
    const auto grid_old = (*this);

    for (int i = 0 ; i < dim ; ++i)
        if (refinement_directions[i])
            this->refine_knots_direction(i,n_subdivisions[i]);

    // refining the objects that's are attached to the CartesianGrid
    // (i.e. that are defined using this CartesianGrid object)
    this->refine_signals_(refinement_directions,grid_old);
}

template <int dim_>
void
CartesianGrid<dim_>::
refine_direction(const int direction_id, const Size n_subdivisions)
{
    Assert(direction_id >= 0 && direction_id < dim, ExcIndexRange(direction_id,0,dim));

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
    out << "CartesianGrid<" << dim << ">" << endl;

    out.push("\t");
    out << "Knot coordinates:" << endl;


    out.push("\t");
    for (int i = 0; i < dim; i++)
    {
        std::vector <Real> vec = this->get_knot_coordinates(i);
        out << "Direction[" << i << "] = " <<  vec << endl;
    }
    out.pop();


    const int num_elements = this->get_num_elements();
    out << "Num elements: " << num_elements << endl;

    auto num_elements_dim = this->get_num_elements_dim();
    out.push("\t");
    for (int i = 0; i < dim; i++)
        out << "Direction[" << i << "] = " << num_elements_dim(i) << endl;
    out.pop();
    out.pop();
}



template <int dim_>
auto
CartesianGrid<dim_>::
get_face_grid(const int face_id, std::map<int,int> &elem_map) const -> shared_ptr<FaceType>
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
const int
CartesianGrid<dim_>::dim ;


IGA_NAMESPACE_CLOSE

#include <igatools/geometry/cartesian_grid.inst>
