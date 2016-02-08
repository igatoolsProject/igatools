//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2016  by the igatools authors (see authors.txt).
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


#ifndef INTERFACE_OPERATOR_H_
#define INTERFACE_OPERATOR_H_

#include <igatools/base/config.h>
#include <igatools/base/tensor.h>
#include <memory>

IGA_NAMESPACE_OPEN

template <int dim> class CartesianGrid;
template <class T> class ValueTable;
template <class T> class ValueVector;
template <class T> class vector;

/**
 * TODO: To be documented!!
 * @brief To be documented.
 *
 * @note To be documented.
 *
 * @author P. Antolin
 * @date 12 May 2015
 */
template <class PhysSpace>
class InterfaceOperator
{
private:
  /** Type for the element accessor of the physical space. */
  using Elem_ = typename PhysSpace::ElementAccessor;

  /** Type for the element iterator of the physical space. */
  using ElemIt_ = typename PhysSpace::ElementIterator;

  /** Return type of the function of the form methods. */
  using Value_ = typename PhysSpace::Value;

  /** Dimension of the space. **/
  static const int dim_ = PhysSpace::dim;

  /** Physical dimension of the space. **/
  static const int space_dim_ = PhysSpace::space_dim;

  /** Dimension of the face. **/
  static const int face_dim_ = PhysSpace::FaceSpace::dim;

  using FaceGridMap_ = typename CartesianGrid<dim_>::FaceGridMap;

  using PhysFaceSpace_ = typename PhysSpace::FaceSpace;

  using Constraints_ = std::map<Index, std::map<Index, Real>>;

  using PointsContainer_ = ValueVector<Points<dim_>>;

  using FacePointsContainer_ = ValueVector<Points<face_dim_>>;

public:

  /**
   * TODO: To be documented!!
   * @brief To be documented.
  */
  class InterfaceElementMap
  {
  public:
    /** @name Constructors */
    ///@{
    /**
    * Default constructor.
    */
    InterfaceElementMap(const Index &elem_id,
                        const Index &interface_elem_id,
                        const Index &face_id,
                        const PointsContainer_ &unit_points,
                        const FacePointsContainer_ &face_unit_points);

    /**
    * Default constructor.
    */
    InterfaceElementMap() = delete;


    /** Copy constructor. */
    InterfaceElementMap(const InterfaceElementMap &in) = default;


    /** Move constructor. */
    InterfaceElementMap(InterfaceElementMap &&in) = default;

    /** Destructor. */
    ~InterfaceElementMap() = default;
    ///@}


    /** @name Assignment operators */
    ///@{
    /** Copy assignment operator. */
    InterfaceElementMap &
    operator=(const InterfaceElementMap &in) = default;


    /** Move assignment operator. */
    InterfaceElementMap &
    operator=(InterfaceElementMap &&in) = default;
    ///@}

    /** Return the id of the element. **/
    const Index &get_element_id() const;

    /** Return the id of the element at the interface grid. **/
    const Index &get_interface_element_id() const;

    /** Return the id of the face. **/
    const Index &get_face_id() const;

    /** Return the unit points. **/
    const PointsContainer_ &get_unit_points() const;

    /** Return the unit points at the element of the interface grid. **/
    const FacePointsContainer_ &get_face_unit_points() const;

  private:
    /** Id of the element. **/
    const Index elem_id_;

    /** Id of the element at the interface grid. **/
    const Index interface_elem_id_;

    /** Face id of th interface. **/
    const Index face_id_;

    /** Unit points at the element face. **/
    const PointsContainer_ unit_points_;

    /** Unit points at the element of the interface grid. **/
    const FacePointsContainer_ face_unit_points_;

  };


  /**
   * TODO: To be documented!!
   * @brief To be documented.
   *
   * @note This class is purely virtual.
  */
  class Form
  {
  protected:
    using InterfaceElementMap_ = typename InterfaceOperator<PhysSpace>::InterfaceElementMap;

  public:
    /** @name Constructors */
    ///@{
    /**
    * Default constructor.
    */
    Form(const std::shared_ptr<const PhysSpace> space_0,
         const std::shared_ptr<const PhysSpace> space_1);


    /** Copy constructor. */
    Form(const Form &in) = default;


    /** Move constructor. */
    Form(Form &&in) = default;

    /** Destructor. */
    ~Form() = default;
    ///@}


    /** @name Assignment operators */
    ///@{
    /** Copy assignment operator. */
    Form &
    operator=(const Form &in) = default;


    /** Move assignment operator. */
    Form &
    operator=(Form &&in) = default;
    ///@}

    /**
     * Evaluation of the form.
     * @note This method is purely virtual.
     *       It must be implemented in derived classes.
     */
    virtual std::array<ValueTable<Value_>, 2>
    evaluate(const InterfaceOperator::InterfaceElementMap &interface_elem_0,
             const InterfaceOperator::InterfaceElementMap &interface_elem_1) = 0;
  protected:

    /** Physical space 0. **/
    const std::shared_ptr<const PhysSpace> space_0_;

    /** Physical space 1. **/
    const std::shared_ptr<const PhysSpace> space_1_;

    /** Element accessor to physical space 0. **/
    ElemIt_ elem_0_;

    /** Element accessor to physical space 0. **/
    ElemIt_ elem_1_;

  };

  /** @name Constructors */
  ///@{
  /**
   * Default constructor.
   * In Debug mode, it checks if the template arguments are consistent.
   */
  InterfaceOperator(const std::shared_ptr<const PhysSpace> space_0,
                    const Index &face_id_0,
                    const std::shared_ptr<const PhysSpace> space_1,
                    const Index &face_id_1,
                    const std::shared_ptr<Form> face_form_trial,
                    const std::shared_ptr<Form> face_form_test,
                    const Real &constant);


  /** Copy constructor. */
  InterfaceOperator(const InterfaceOperator<PhysSpace> &in) = default;


  /** Move constructor. */
  InterfaceOperator(InterfaceOperator<PhysSpace> &&in) = default;

  /** Destructor. */
  ~InterfaceOperator() = default;
  ///@}


  /** @name Assignment operators */
  ///@{
  /** Copy assignment operator. */
  InterfaceOperator<PhysSpace> &
  operator=(const InterfaceOperator<PhysSpace> &in) = default;


  /** Move assignment operator. */
  InterfaceOperator<PhysSpace> &
  operator=(InterfaceOperator<PhysSpace> &&in) = default;
  ///@}


  /**
   * TODO: To be documented!!
   * \f[
        (A_e)_{ij} = \int_{\Omega_e} \sum_{s=1}^{sp\_dim}
        \phi^{e,\text{test}}_i
        \, \beta_{s}(x) \,
        \bigl( \nabla \phi^{e,\text{trial}}_j \bigr)_s \; d \Omega
        = \int_{\Omega_e}
        \phi^{e,\text{test}}_i
        \, \vec{\beta}(x) \, \cdot \,
        \nabla \phi^{e,\text{trial}}_j \; d \Omega .
     \f]
   */
  void evaluate(Constraints_ &constraints);

private:
  const std::shared_ptr<const PhysSpace> space_0_;
  const std::shared_ptr<const PhysSpace> space_1_;
  const Index face_id_0_;
  const Index face_id_1_;
  const std::shared_ptr<Form> face_form_trial_;
  const std::shared_ptr<Form> face_form_test_;
  const Real integration_constant_;

  void get_element_and_quad_on_face(
    const FacePointsContainer_ &eval_pts_ref_domain,
    const std::shared_ptr<const PhysSpace> space,
    const std::shared_ptr<const PhysFaceSpace_> face_space,
    const Index &face_id,
    const vector<Index> &face_space_dofs,
    const FaceGridMap_ &face_grid_to_space_grid,
    const Index face_element_id,
    Index &element_id,
    vector<Index> &face_elem_dofs,
    PointsContainer_ &current_elem_pts_unit_domain) const;

};



template <class PhysSpace>
class ValueJump : public InterfaceOperator<PhysSpace>::Form
{

private:
  /** Type for the element accessor of the physical space. */
  using Elem_ = typename PhysSpace::ElementAccessor;

  /** Return type of the function of the form methods. */
  using Value_ = typename PhysSpace::Value;

  /** Dimension of the space. **/
  static const int dim_ = PhysSpace::dim;

  /** Physical dimension of the space. **/
  static const int space_dim_ = PhysSpace::space_dim;

  /** Dimension of the face. **/
  static const int face_dim_ = PhysSpace::FaceSpace::dim;

  /** Base class from with this class derives. **/
  using Base_ = typename InterfaceOperator<PhysSpace>::Form;

  /** Shared pointer of a const type ot the base class. **/
  using BaseConstPtr_ = std::shared_ptr<Base_>;

  using PointsContainer_ = ValueVector<Points<dim_>>;

  using InterfaceElementMap_ = typename Base_::InterfaceElementMap_;

public:
  /** @name Constructors */
  ///@{
  /**
  * Constructor.
  */
  ValueJump(const std::shared_ptr<const PhysSpace> space_0,
            const std::shared_ptr<const PhysSpace> space_1);

  /**
  * Default constructor.
  */
  ValueJump() = delete;


  /** Copy constructor. */
  ValueJump(const ValueJump &in) = default;


  /** Move constructor. */
  ValueJump(ValueJump &&in) = default;

  /** Destructor. */
  ~ValueJump() = default;
  ///@}


  /** @name Assignment operators */
  ///@{
  /** Copy assignment operator. */
  ValueJump &
  operator=(const ValueJump &in) = default;


  /** Move assignment operator. */
  ValueJump &
  operator=(ValueJump &&in) = default;
  ///@}

  /**
    * Creates a new instantiation of the class and returns wrapped into a
    * shared pointer.
    * @note Return a shared pointer with a new object.
    */
  static BaseConstPtr_ create(const std::shared_ptr<const PhysSpace> space_0,
                              const std::shared_ptr<const PhysSpace> space_1);

  /**
    * Evaluation of the form.
    * @note This method is purely virtual.
    *       It must be implemented in derived classes.
    */
  virtual std::array<ValueTable<Value_>, 2>
  evaluate(const InterfaceElementMap_& interface_elem_0,
           const InterfaceElementMap_& interface_elem_1) override final;
};



template <class PhysSpace>
class NormalGradientJump : public InterfaceOperator<PhysSpace>::Form
{

private:
  /** Type for the element accessor of the physical space. */
  using Elem_ = typename PhysSpace::ElementAccessor;

  /** Return type of the function of the form methods. */
  using Value_ = typename PhysSpace::Value;

  /** Dimension of the space. **/
  static const int dim_ = PhysSpace::dim;

  /** Physical dimension of the space. **/
  static const int space_dim_ = PhysSpace::space_dim;

  /** Dimension of the face. **/
  static const int face_dim_ = PhysSpace::FaceSpace::dim;

  /** Base class from with this class derives. **/
  using Base_ = typename InterfaceOperator<PhysSpace>::Form;

  /** Shared pointer of a const type ot the base class. **/
  using BaseConstPtr_ = std::shared_ptr<Base_>;

  using PointsContainer_ = ValueVector<Points<dim_>>;

  using InterfaceElementMap_ = typename Base_::InterfaceElementMap_;

public:
  /** @name Constructors */
  ///@{
  /**
  * Constructor.
  */
  NormalGradientJump(const std::shared_ptr<const PhysSpace> space_0,
                     const std::shared_ptr<const PhysSpace> space_1);

  /**
  * Default constructor.
  */
  NormalGradientJump() = delete;


  /** Copy constructor. */
  NormalGradientJump(const NormalGradientJump &in) = default;


  /** Move constructor. */
  NormalGradientJump(NormalGradientJump &&in) = default;

  /** Destructor. */
  ~NormalGradientJump() = default;
  ///@}


  /** @name Assignment operators */
  ///@{
  /** Copy assignment operator. */
  NormalGradientJump &
  operator=(const NormalGradientJump &in) = default;


  /** Move assignment operator. */
  NormalGradientJump &
  operator=(NormalGradientJump &&in) = default;
  ///@}

  /**
    * Creates a new instantiation of the class and returns wrapped into a
    * shared pointer.
    * @note Return a shared pointer with a new object.
    */
  static BaseConstPtr_ create(const std::shared_ptr<const PhysSpace> space_0,
                              const std::shared_ptr<const PhysSpace> space_1);

  /**
    * Evaluation of the form.
    * @note This method is purely virtual.
    *       It must be implemented in derived classes.
    */
  virtual std::array<ValueTable<Value_>, 2>
  evaluate(const InterfaceElementMap_& interface_elem_0,
           const InterfaceElementMap_& interface_elem_1) override final;
};



template <class PhysSpace>
class NormalGradientAverage : public InterfaceOperator<PhysSpace>::Form
{

private:
  /** Type for the element accessor of the physical space. */
  using Elem_ = typename PhysSpace::ElementAccessor;

  /** Return type of the function of the form methods. */
  using Value_ = typename PhysSpace::Value;

  /** Dimension of the space. **/
  static const int dim_ = PhysSpace::dim;

  /** Physical dimension of the space. **/
  static const int space_dim_ = PhysSpace::space_dim;

  /** Dimension of the face. **/
  static const int face_dim_ = PhysSpace::FaceSpace::dim;

  /** Base class from with this class derives. **/
  using Base_ = typename InterfaceOperator<PhysSpace>::Form;

  /** Shared pointer of a const type ot the base class. **/
  using BaseConstPtr_ = std::shared_ptr<Base_>;

  using PointsContainer_ = ValueVector<Points<dim_>>;

  using InterfaceElementMap_ = typename Base_::InterfaceElementMap_;

public:
  /** @name Constructors */
  ///@{
  /**
  * Constructor.
  */
  NormalGradientAverage(const std::shared_ptr<const PhysSpace> space_0,
                        const std::shared_ptr<const PhysSpace> space_1,
                        const Real &weight = Real(0.5));

  /**
  * Default constructor.
  */
  NormalGradientAverage() = delete;


  /** Copy constructor. */
  NormalGradientAverage(const NormalGradientAverage &in) = default;


  /** Move constructor. */
  NormalGradientAverage(NormalGradientAverage &&in) = default;

  /** Destructor. */
  ~NormalGradientAverage() = default;
  ///@}


  /** @name Assignment operators */
  ///@{
  /** Copy assignment operator. */
  NormalGradientAverage &
  operator=(const NormalGradientAverage &in) = default;


  /** Move assignment operator. */
  NormalGradientAverage &
  operator=(NormalGradientAverage &&in) = default;
  ///@}

  /**
    * Creates a new instantiation of the class and returns wrapped into a
    * shared pointer.
    * @note Return a shared pointer with a new object.
    */
  static BaseConstPtr_ create(const std::shared_ptr<const PhysSpace> space_0,
                              const std::shared_ptr<const PhysSpace> space_1,
                              const Real &weight = Real(0.5));

  /**
    * Evaluation of the form.
    * @note This method is purely virtual.
    *       It must be implemented in derived classes.
    */
  virtual std::array<ValueTable<Value_>, 2>
  evaluate(const InterfaceElementMap_& interface_elem_0,
           const InterfaceElementMap_& interface_elem_1) override final;
private:
  /** Weight average */
  const Real weight_average_;
};


template <class PhysSpace>
class LaplacianAverage : public InterfaceOperator<PhysSpace>::Form
{

private:
  /** Type for the element accessor of the physical space. */
  using Elem_ = typename PhysSpace::ElementAccessor;

  /** Return type of the function of the form methods. */
  using Value_ = typename PhysSpace::Value;

  /** Dimension of the space. **/
  static const int dim_ = PhysSpace::dim;

  /** Physical dimension of the space. **/
  static const int space_dim_ = PhysSpace::space_dim;

  /** Dimension of the face. **/
  static const int face_dim_ = PhysSpace::FaceSpace::dim;

  /** Base class from with this class derives. **/
  using Base_ = typename InterfaceOperator<PhysSpace>::Form;

  /** Shared pointer of a const type ot the base class. **/
  using BaseConstPtr_ = std::shared_ptr<Base_>;

  using PointsContainer_ = ValueVector<Points<dim_>>;

  using InterfaceElementMap_ = typename Base_::InterfaceElementMap_;

public:
  /** @name Constructors */
  ///@{
  /**
  * Constructor.
  */
  LaplacianAverage(const std::shared_ptr<const PhysSpace> space_0,
                   const std::shared_ptr<const PhysSpace> space_1,
                   const Real &weight = Real(0.5));

  /**
  * Default constructor.
  */
  LaplacianAverage() = delete;


  /** Copy constructor. */
  LaplacianAverage(const LaplacianAverage &in) = default;


  /** Move constructor. */
  LaplacianAverage(LaplacianAverage &&in) = default;

  /** Destructor. */
  ~LaplacianAverage() = default;
  ///@}


  /** @name Assignment operators */
  ///@{
  /** Copy assignment operator. */
  LaplacianAverage &
  operator=(const LaplacianAverage &in) = default;


  /** Move assignment operator. */
  LaplacianAverage &
  operator=(LaplacianAverage &&in) = default;
  ///@}

  /**
    * Creates a new instantiation of the class and returns wrapped into a
    * shared pointer.
    * @note Return a shared pointer with a new object.
    */
  static BaseConstPtr_ create(const std::shared_ptr<const PhysSpace> space_0,
                              const std::shared_ptr<const PhysSpace> space_1,
                              const Real &weight = Real(0.5));

  /**
    * Evaluation of the form.
    * @note This method is purely virtual.
    *       It must be implemented in derived classes.
    */
  virtual std::array<ValueTable<Value_>, 2>
  evaluate(const InterfaceElementMap_& interface_elem_0,
           const InterfaceElementMap_& interface_elem_1) override final;

private:
  /** Weight average */
  const Real weight_average_;
};


IGA_NAMESPACE_CLOSE


#endif // #ifndef INTERFACE_OPERATOR_H_
