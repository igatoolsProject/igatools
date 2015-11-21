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

#ifndef __PHYSICAL_SPACE_H_
#define __PHYSICAL_SPACE_H_

#include <igatools/base/config.h>
#include <igatools/basis_functions/basis.h>
#include <igatools/basis_functions/reference_space_basis.h>
#include <igatools/basis_functions/dof_distribution.h>
#include <igatools/geometry/domain.h>
#include <igatools/geometry/grid_iterator.h>
#include <igatools/utils/static_multi_array.h>

IGA_NAMESPACE_OPEN

class SpaceManager;

template <int,int> class PushForward;

template <int,int,int,int> class PhysicalSpaceElement;

template <int,int,int,int> class PhysSpaceElementHandler;

/**
 *
 *
 * @ingroup containers
 * @ingroup serializable
 *
 */
template <int dim_, int range_= 1, int rank_ = 1, int codim_ = 0>
class PhysicalSpaceBasis :
  public Basis<dim_,codim_,range_,rank_>
{
private:
  using base_t = Basis<dim_,codim_,range_,rank_>;
  using self_t = PhysicalSpaceBasis<dim_, range_, rank_, codim_>;

public:
  ///@{
  /**
   * See documentation in \ref Basis
   *
   * @see Basis
   */
  using PushFwd = PushForward<dim_, codim_>;

  using PhysDomain = Domain<dim_, codim_>;

  using RefBasis = ReferenceSpaceBasis<dim_,range_,rank_>;

  using GridType = Grid<dim_>;
  ///@}
  using ElementHandler = PhysSpaceElementHandler<dim_,range_,rank_,codim_>;

  static const int dim = dim_;

  static const int codim = codim_;

  static const int space_dim = dim+codim;

//  static const int range = PushFwd::template PhysRange<range_>::value;
  static const int range = range_;

  static const int rank = rank_;

  static const bool is_physical_space = true;

  static constexpr int n_components = constexpr_pow(range, rank);

  static const SafeSTLArray<int, n_components> components;


  using IndexType = typename GridType::IndexType;
  using PropertyList = PropertiesIdContainer<IndexType>;
  using List = typename PropertyList::List;
  using ListIt = typename PropertyList::List::iterator;

public:
  template <int order>
  using Derivative = typename base_t::template Derivative<order>;
  using Point = typename base_t::Point;
  using Value = typename base_t::Value;
  using Gradient = typename base_t::Gradient;
  using Div   = typename base_t::Div;

  using RefPoint = typename RefBasis::Point;


public:
  template< class T>
  using ComponentContainer = typename RefBasis::template ComponentContainer<T>;


public:


  using ElementAccessor = PhysicalSpaceElement<dim_,range_,rank_,codim_>;
  using ElementIterator = GridIterator<ElementAccessor>;

  /**
   * Default constructor. It does nothing but it is needed for the serialization
   * mechanism.
   */
  PhysicalSpaceBasis() = default;

  PhysicalSpaceBasis(const self_t &phys_space) = delete;

  virtual ~PhysicalSpaceBasis() = default;

  static std::shared_ptr<self_t>
  create(const std::shared_ptr<RefBasis> &ref_basis,
         const std::shared_ptr<PhysDomain> &phys_domain,
         const Transformation &transformation_type = Transformation::h_grad);

  static std::shared_ptr<const self_t>
  const_create(const std::shared_ptr<const RefBasis> &ref_basis,
               const std::shared_ptr<const PhysDomain> &phys_domain,
               const Transformation &transformation_type = Transformation::h_grad);

  /**
   * Create an element (defined on this grid) with a given index.
   */
  std::unique_ptr<SpaceElement<dim_,codim_,range_,rank_> >
  create_element(const ListIt &index, const PropId &property) const override final;





  template <int k>
  using SubSpace = PhysicalSpaceBasis<k, range, rank, codim + dim-k>;

  template <int sdim>
  using SubGridMap = typename RefBasis::GridType::template SubGridMap<sdim>;

//    using InterGridMap = std::map<Index,Index>;

  template <int k>
  using InterSpaceMap = typename RefBasis::template InterSpaceMap<k>;



  template<int k>
  std::shared_ptr<const SubSpace<k> >
  get_sub_space(const int s_id, InterSpaceMap<k> &dof_map,
                const std::shared_ptr<const Grid<k>> &sub_grid,
                SubGridMap<k> &elem_map,
                EnableIf<(dim_ != 0) &&(k>=0)> * = nullptr) const;







  void print_info(LogStream &out) const override final;

  std::unique_ptr<SpaceElementHandler<dim_,codim_,range_,rank_> >
  create_cache_handler() const override final;


#if 0
  /**
   * Return the maximum value of the polynomial degree, for each component, for each direction;
   */
  virtual int get_max_degree() const override final;
#endif

  std::shared_ptr<const Domain<dim_,codim_>> get_physical_domain() const;


  virtual std::shared_ptr<const Grid<dim_>> get_grid() const override final;


  std::shared_ptr<const SplineSpace<dim_,range_,rank_>>
                                                     get_spline_space() const override final;


  std::shared_ptr<const RefBasis> get_reference_basis() const;

//  std::shared_ptr<RefBasis> get_reference_basis();

  Transformation get_transformation_type() const;

private:


  PhysicalSpaceBasis(const SharedPtrConstnessHandler<RefBasis> &ref_basis,
                     const SharedPtrConstnessHandler<PhysDomain> &phys_domain,
                     const Transformation &transformation_type);



  SharedPtrConstnessHandler<RefBasis> ref_basis_;

  SharedPtrConstnessHandler<PhysDomain> phys_domain_;

  const Transformation transformation_type_ = Transformation::h_grad;


  std::shared_ptr<const self_t> phys_basis_previous_refinement_ = nullptr;


  friend ElementAccessor;


  /**
   * Returns the current object wrapped by a std::shared_ptr.
   *
   * @note Internally uses the shared_from_this() function.
   */
  std::shared_ptr<const self_t > get_this_basis() const;


#ifdef MESH_REFINEMENT

  /**
   * Rebuild the internal state of the object after an insert_knots() function is invoked.
   *
   * @pre Before invoking this function, must be invoked the function grid_->insert_knots().
   * @note This function is connected to the Grid's signal for the refinement, and
   * it is necessary in order to avoid infinite loops in the insert_knots() function calls.
   *
   * @ingroup h_refinement
   */
  void rebuild_after_insert_knots(
    const SafeSTLArray<SafeSTLVector<Real>,dim> &knots_to_insert,
    const Grid<dim> &old_grid) override final;

public:

  virtual void refine_h(const Size n_subdivisions = 2) override final;

  void create_connection_for_insert_knots(const std::shared_ptr<self_t> &space);


  std::shared_ptr<const base_t> get_basis_previous_refinement() const;

#endif // MESH_REFINEMENT

private:

#ifdef SERIALIZATION
  /**
   * @name Functions needed for the serialization
   * @see <a href="http://uscilab.github.io/cereal/">Cereal serialization library</a>
   */
  ///@{
  friend class cereal::access;

  template<class Archive>
  void
  serialize(Archive &ar)
  {
    using std::to_string;
    const std::string base_name = "Space_" +
                                  to_string(dim_) + "_" +
                                  to_string(codim_) + "_" +
                                  to_string(range_) + "_" +
                                  to_string(rank_);

    ar &make_nvp(base_name,base_class<base_t>(this));


    ar &make_nvp("ref_basis_",ref_basis_);

    ar &make_nvp("phys_domain_",phys_domain_);

    Transformation transformation_type_tmp = transformation_type_;
    ar &make_nvp("transformation_type_",transformation_type_tmp);
    const_cast<Transformation &>(transformation_type_) = transformation_type_tmp;

#ifdef MESH_REFINEMENT
    auto tmp = std::const_pointer_cast<self_t>(phys_basis_previous_refinement_);
    ar &make_nvp("phys_basis_previous_refinement_",tmp);
    phys_basis_previous_refinement_ = tmp;
#endif // MESH_REFINEMENT
  }

  ///@}
#endif // SERIALIZATION


};

IGA_NAMESPACE_CLOSE

#endif
