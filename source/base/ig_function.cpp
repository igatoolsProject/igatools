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

#include <igatools/base/ig_function.h>
#include <igatools/base/function_element.h>

using std::shared_ptr;

IGA_NAMESPACE_OPEN

template<class Space>
IgFunction<Space>::
IgFunction(std::shared_ptr<const Space> space,
           const CoeffType &coeff)
    :
    parent_t::Function(space->get_grid()),
    space_(space),
    coeff_(coeff),
    elem_(space->begin()),
    space_filler_(space->create_elem_handler())
{
    Assert(space_ != nullptr,ExcNullPtr());
    Assert(!coeff_.empty(),ExcEmptyObject());
}



template<class Space>
IgFunction<Space>::
IgFunction(const self_t &fun)
    :
    parent_t::Function(fun.space_->get_grid()),
    space_(fun.space_),
    coeff_(fun.coeff_),
    elem_(fun.space_->begin()),
    space_filler_(fun.space_->create_elem_handler())
{
    Assert(space_ != nullptr,ExcNullPtr());
    Assert(!coeff_.empty(),ExcEmptyObject());
}



template<class Space>
auto
IgFunction<Space>::
create(std::shared_ptr<const Space> space,
       const CoeffType &coeff) ->  std::shared_ptr<self_t>
{
    return std::shared_ptr<self_t>(new self_t(space, coeff));
}



template<class Space>
auto
IgFunction<Space>::
reset(const ValueFlags &flag, const eval_pts_variant &eval_pts) -> void
{
    parent_t::reset(flag, eval_pts);
    reset_impl.flag = flag;
    reset_impl.space_handler_ = space_filler_.get();
    reset_impl.flags_ = &(this->flags_);
    boost::apply_visitor(reset_impl, eval_pts);
}



template<class Space>
auto
IgFunction<Space>::
init_cache(ElementAccessor &elem, const topology_variant &k) -> void
{
    parent_t::init_cache(elem, k);
    init_cache_impl.space_handler_ = space_filler_.get();
    init_cache_impl.space_elem = &(*elem_);
    boost::apply_visitor(init_cache_impl, k);
}



template<class Space>
auto
IgFunction<Space>::
fill_cache(ElementAccessor &elem, const topology_variant &k, const int j) -> void
{
    parent_t::fill_cache(elem,k,j);

    elem_.move_to(elem.get_flat_index());

    fill_cache_impl.space_handler_ = space_filler_.get();
    fill_cache_impl.space_elem = &(*elem_);
    fill_cache_impl.func_elem = &elem;
    fill_cache_impl.function = this;


    // TODO (pauletti, Nov 27, 2014): if code is in final state remove commented line else fix
    const auto local_ids = elem_->get_local_to_global();
    vector<Real> loc_coeff;
    for (const auto &id : local_ids)
        loc_coeff.push_back(coeff_[id]);
//    auto loc_coeff = coeff_.get_local_coefs(elem_->get_local_to_global());

    fill_cache_impl.loc_coeff = &loc_coeff;
    fill_cache_impl.j =j;

    boost::apply_visitor(fill_cache_impl, k);
}



template<class Space>
auto
IgFunction<Space>::
get_iga_space() const -> std::shared_ptr<const Space>
{
    return space_;
}



template<class Space>
auto
IgFunction<Space>::
get_coefficients() const -> const CoeffType &
{
    return coeff_;
}



template<class Space>
auto
IgFunction<Space>::
operator +=(const self_t &fun) -> self_t &
{
    const auto size = coeff_.size();
    for (int i=0; i<size; ++i)
        coeff_[i] += fun.coeff_[i];

    return *this;
}



template<class Space>
void
IgFunction<Space>::
print_info(LogStream &out) const
{
    out.begin_item("Reference space info:");
    space_->print_info(out);
    out.end_item();
    out << std::endl;

#if 0
    out << "Control points info (projective coordinates):" << endl;

    //write the projective cooridnates if the reference space is NURBS
#endif

    out.begin_item("Control points info (euclidean coordinates):");
    coeff_.print_info(out);
    out.end_item();
}

IGA_NAMESPACE_CLOSE

#include <igatools/base/ig_function.inst>
