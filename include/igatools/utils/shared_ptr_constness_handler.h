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

#ifndef SHARED_PTR_CONSTNESS_HANDLER_H_
#define SHARED_PTR_CONSTNESS_HANDLER_H_

#include <igatools/base/config.h>
#include <igatools/base/exceptions.h>
#include <igatools/base/logstream.h>




IGA_NAMESPACE_OPEN


/**
 * @brief Class to automatically manage the constness property of an object wrapped by a shared_ptr,
 * depending on the constness property of the object passed in the constructor argument.
 *
 * This class has to constructors:
 * - one is accepting a pointer-to-const object;
 * - the other one is accepting a pointer-to-non-const object.
 *
 * Depending on which constructor is used, this class keep track of the constess property of the pointee
 * object.
 *
 * If the object used in the constructor argument is a <b>pointer-to-const</b>, then all the getters
 * returning pointer-to-non-const or non-const reference will give an error in Debug mode.
 *
 * If the object used in the constructor argument is a <b>pointer-to-non-const</b>, then all the getters
 * can be used (both the one returning the const and non-const object).
 *
 */
template<class T>
class SharedPtrConstnessHandler
{
public:
  using Ptr = std::shared_ptr<T>;
  using PtrToConst = std::shared_ptr<const T>;

  /**
   * @name Constructors and destructor
   */
  ///@{

  /**
   * Default constructor. It does nothing but it is needed for the
   * serialization mechanism.
   */
  explicit SharedPtrConstnessHandler() = default;

  /**
   * Constructs the object using a shared pointer to non-const data.
   *
   * @note In DEBUG mode is performed a check about the nullity of the input data:
   * an assertion will be raised if the input data is a nullptr.
   */
  explicit SharedPtrConstnessHandler(const std::shared_ptr<T> &data)
    :
    ptr_(data),
    ptr_to_const_(nullptr)
  {
    Assert(ptr_ != nullptr, ExcNullPtr());
  };

  /**
   * Constructs the object using a shared pointer to const data.
   */
  explicit SharedPtrConstnessHandler(const std::shared_ptr<const T> &data)
    :
    ptr_(nullptr),
    ptr_to_const_(data)
  {
    Assert(ptr_to_const_ != nullptr, ExcNullPtr());
  };

#if 0
  SharedPtrConstnessHandler(Ptr &&data)
    :
    data_is_const_(false)
  {
    ptr_ = data;
    Assert(ptr_ != nullptr, ExcNullPtr());
  };

  SharedPtrConstnessHandler(PtrToConst &&data)
    :
    data_is_const_(true)
  {
    ptr_to_const_ = data;
    Assert(ptr_to_const_ != nullptr, ExcNullPtr());
  };
#endif
  /**
   * Copy constructor;
   */
  SharedPtrConstnessHandler(const SharedPtrConstnessHandler<T> &obj) = default;

  /**
   * Move constructor;
   */
  SharedPtrConstnessHandler(SharedPtrConstnessHandler<T> &&obj) = default;


  ~SharedPtrConstnessHandler() = default;
  ///@}


  /**
   * @name Assignment operators.
   */
  ///@{
  /**
   * Copy assignment operator. Not allowed to be used.
   */
  SharedPtrConstnessHandler<T> &operator=(const SharedPtrConstnessHandler<T> &obj) = default;

  /**
   * Move assignment operator. Not allowed to be used.
   */
  SharedPtrConstnessHandler<T> &operator=(SharedPtrConstnessHandler<T> &&obj) = default;
  ///@}



  /**
   * Return a non-const reference of the data.
   */
  T &get_ref_data()
  {
    return *this->get_ptr_data();
  }

  /**
   * Return a copy of the shared pointer that is pointing to the non-const data.
   */
  std::shared_ptr<T> get_ptr_data() const
  {
    Assert(!data_is_const(),ExcMessage("Data is built as const."));
    return ptr_;
  }

  /**
   * Return the non-const reference of the shared pointer that is pointing to the non-const data.
   */
  std::shared_ptr<T> &get_ref_ptr_data()
  {
    Assert(!data_is_const(),ExcMessage("Data is built as const."));
    return ptr_;
  }


  /**
   * Return a const reference of the data.
   */
  const T &get_ref_const_data() const
  {
    return *this->get_ptr_const_data();
  }

  /**
   * Return a copy of the shared pointer that is pointing to the const data.
   */
  std::shared_ptr<const T> get_ptr_const_data() const
  {
    if (data_is_const())
      return ptr_to_const_;
    else
      return ptr_;
  }

  /**
   * Return the non-const reference of the shared pointer that is pointing to the const data.
   */
  std::shared_ptr<const T> &get_ref_ptr_const_data() const
  {
    if (data_is_const())
      return ptr_to_const_;
    else
      return ptr_;
  }

  void print_info(LogStream &out) const
  {
    if (data_is_const())
    {
      out.begin_item("Pointer to const data");
      ptr_to_const_->print_info(out);
    }
    else
    {
      out.begin_item("Pointer to non-const data");
      ptr_->print_info(out);
    }
    out.end_item();
  }

  /**
   * Dereference operator. Return a const reference of the data.
   */
  const T &operator*() const
  {
    return this->get_ref_const_data();
  }


  /**
   * Dereference operator. Return a pointer to the const data.
   */
  PtrToConst operator->() const
  {
    if (data_is_const())
      return ptr_to_const_;
    else
      return ptr_;
  }

  /**
   * Check if unique.
   *
   * Returns whether the shared_ptr object does not share ownership over its pointer with other shared_ptr objects (i.e., it is unique).
   */
  bool unique() const
  {
    if (data_is_const())
      return ptr_to_const_.unique();
    else
      return ptr_.unique();
  }

  /**
   * Returns TRUE is the pointee data is a const object.
   */
  bool data_is_const() const
  {
    return (ptr_to_const_ ? true : false);
  }


private:

  Ptr ptr_ = nullptr;

  PtrToConst ptr_to_const_ = nullptr;




#ifdef IGATOOLS_WITH_SERIALIZATION
  /**
   * @name Functions needed for boost::serialization
   * @see <a href="http://uscilab.github.io/cereal/serialization_functions.html">Cereal serialization</a>
   */
  ///@{
  friend class cereal::access;

  template<class Archive>
  void
  serialize(Archive &ar, const unsigned int version)
  {
//    ar &make_nvp("data_is_const_",data_is_const_);

    // In order to serialize the data, we need to cast it to non-const
    Ptr tmp;
    if (ptr_to_const_)
      tmp = std::const_pointer_cast<T>(ptr_to_const_);
    else
      tmp = ptr_;

    ar &make_nvp("tmp_ptr_to_data_",tmp);
    Assert(tmp != nullptr,ExcNullPtr());

    // After deserialization we need to cast the data to the correct constness
    if (ptr_to_const_)
      ptr_to_const_ = std::const_pointer_cast<const T>(tmp);
    else
      ptr_ = tmp;
  }
  ///@}
#endif // IGATOOLS_WITH_SERIALIZATION

};


IGA_NAMESPACE_CLOSE


#endif // SHARED_PTR_CONSTNESS_HANDLER_H_
