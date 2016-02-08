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



#ifndef CONCATENATED_ITERATOR_H_
#define CONCATENATED_ITERATOR_H_

#include <igatools/base/config.h>
#include <igatools/base/logstream.h>
#include <igatools/utils/container_view.h>
#include <igatools/utils/safe_stl_vector.h>


IGA_NAMESPACE_OPEN


template <class ViewType,class DerivedClass>
class ConcatenatedIteratorData
  : public std::iterator<std::forward_iterator_tag, typename ViewType::iterator::value_type>
{
public:
  using Iterator = typename ViewType::iterator;

  /**
   * Alias for specifying the value_type of the iterator used by the ViewType
   * (it's the same of the iterator specified by the alias <tt>Iterator</tt>).
   */
  using value_type = typename Iterator::value_type;


  SafeSTLVector<ViewType> get_ranges() const;

  int get_range_id() const;

  Iterator get_iterator_current() const;

  /** @name Comparison operators */
  ///@{
  /** Compare for equality.*/
  bool operator==(const ConcatenatedIteratorData<ViewType,DerivedClass> &it) const;

  /** Compare for inequality.*/
  bool operator!=(const ConcatenatedIteratorData<ViewType,DerivedClass> &it) const;


  bool operator<(const ConcatenatedIteratorData<ViewType,DerivedClass> &it) const;

  bool operator<=(const ConcatenatedIteratorData<ViewType,DerivedClass> &it) const;

  ///@}


  /** @name Advance operator */
  ///@{
  /**
   *  Prefix <tt>++</tt> operator: <tt>++i</tt>. This
   *  operator advances the iterator to
   *  the next element and returns
   *  a reference to <tt>*this</tt>.
   */
  DerivedClass &operator++();

  DerivedClass &operator+(const int n);
  ///@}


  /** Prints some information. Mostly used for debug and testing. */
  void print_info(LogStream &out) const;


  /** @name Dereferencing operators (const version) */
  ///@{
  /**
   *  Dereferencing operator, returns a
   *  const reference to the value_type.
   */
  const value_type &operator*() const;


  /**
   *  Dereferencing operator, returns a
   *  const pointer to the value_type.
   */
  const value_type *operator->() const;
  ///@}


  /**
   * Returns @p n such that <tt>a + n == (*this)</tt>, where <tt>(*this) == a + ((*this) - a)</tt>.
   */
  Size
  operator-(const ConcatenatedIteratorData<ViewType,DerivedClass> &a) const;



protected:
  ConcatenatedIteratorData();

  ConcatenatedIteratorData(const SafeSTLVector<ViewType> &views,const Index index);

  ConcatenatedIteratorData(const ConcatenatedIteratorData<ViewType,DerivedClass> &rhs) = default;

  ConcatenatedIteratorData(ConcatenatedIteratorData<ViewType,DerivedClass> &&rhs) = default;

  ConcatenatedIteratorData<ViewType,DerivedClass> &operator=(
    const ConcatenatedIteratorData<ViewType,DerivedClass> &rhs) = default;

  ConcatenatedIteratorData<ViewType,DerivedClass> &operator=(
    ConcatenatedIteratorData<ViewType,DerivedClass> &&rhs) = default;

  ~ConcatenatedIteratorData() = default;

  DerivedClass &as_derived_class();

  const DerivedClass &as_derived_class() const;


  /**
   * Vector of ranges upon which the ConcatenatedIterator is defined.
   * Each entry in the vector is a pair of objects of type <tt>Iterator</tt>
   * in the the form [begin,end), telling which is the begin of the range and which
   * is one-pass-end of the range.
   */
  SafeSTLVector<ViewType> ranges_;

  /**
   * Index used to specify which range is spanned at a given moment by the
   * iterator_current_ member variable.
   */
  int range_id_;

  /**
   * Iterator pointing to the current position.
   */
  Iterator iterator_current_;


  void get_range_id_and_entry_id_in_range(const Index id, Index &rng_id, Index &entry_id_rng) const;


  /** Returns a reference to the entry with the identifier specified by @p id. */
  const typename ViewType::reference get_entry_const_reference(const Index id) const;

  /** Returns a const reference to the entry with the identifier specified by @p id. */
  typename ViewType::reference get_entry_reference(const Index id);

private:

  bool is_comparable(const ConcatenatedIteratorData<ViewType,DerivedClass> &it) const;
};




/**
 * @brief This class represents a forward iterator made by the
 * concatenation of several View objects.
 *
 * Basically it provides the same functionality of ConcatenatedConstIterator plus
 * two methods for dereferencing the iterator to non-const value reference and value pointer.
 *
 * @author M. Martinelli
 * @date 03 June 2014
 */
template <class ViewType>
class ConcatenatedIterator
  :
  public ConcatenatedIteratorData< ViewType, ConcatenatedIterator<ViewType> >
{
public:
  using value_type = typename ViewType::iterator::value_type;


  /** @name Constructors & destructor */
  ///@{
  /**
   * Default constructor. It does nothing.
   */
  ConcatenatedIterator() = default;

  /**
   * Constructor.
   */
  ConcatenatedIterator(
    const SafeSTLVector<ViewType> &ranges,
    const Index index);

  /** Copy constructor. */
  ConcatenatedIterator(const ConcatenatedIterator<ViewType> &it) = default;

  /** Move constructor. */
  ConcatenatedIterator(ConcatenatedIterator<ViewType> &&it) = default;

  /** Destructor */
  ~ConcatenatedIterator() = default ;
  ///@}


  /** @name Assignment operators */
  ///@{
  /** Copy assignment operator. */
  ConcatenatedIterator<ViewType> &operator=(
    const ConcatenatedIterator<ViewType> &it) = default;

  /** Move assignment operator. */
  ConcatenatedIterator<ViewType> &operator=(
    ConcatenatedIterator<ViewType> &&it) = default;
  ///@}

  /** @name Dereferencing operators (non-const version) */
  ///@{
  /**
   *  Dereferencing operator, returns a
   *  reference to the value_type.
   */
  value_type &operator*();


  /**
   *  Dereferencing operator, returns a
   *  pointer to the value_type.
   */
  value_type *operator->();
  ///@}




  /** Returns a reference to the entry with the identifier specified by @p id. */
  typename ViewType::reference operator[](const Index id);

  /** Returns a const reference to the entry with the identifier specified by @p id. */
  const typename ViewType::reference operator[](const Index id) const;
};



/**
 * @brief This class represents a const iterator made by the
 * concatenation of several forward iterator.
 *
 * This class can be used to define a const iterator over elements
 * traversed by a certain number of iterators (for which the type is specified
 * by the template parameter <tt>Iterator</tt>) in order to avoid multiple loops.
 *
 * For example, if we want to iterate over the entries of the four vectors
 * @code{.cpp}
   SafeSTLVector<int> v0 = {1,2,3,4};
   SafeSTLVector<int> v1 = {5,6,7,8,9};
   SafeSTLVector<int> v2 = {10,11,12};
   SafeSTLVector<int> v3 = {13};
   @endcode
 * we can use different approaches:
 * - use one range-for-loop for each vector
 * @code{.cpp}
   for (const auto & v : v0)
       cout << v << endl;
   for (const auto & v : v1)
       cout << v << endl;
   for (const auto & v : v2)
       cout << v << endl;
   for (const auto & v : v3)
       cout << v << endl;
   @endcode
 * - define eight iterators (two for each vectors, one representing the current position and
 * the other representing the one-pass-end position) and then use four loops for the iteration
 * @code{.cpp}
   for (auto v = v0.cbegin() ; v != v0.cend() ; ++v)
       cout << *v << endl;
   for (auto v = v1.cbegin() ; v != v1.cend() ; ++v)
       cout << *v << endl;
   for (auto v = v2.cbegin() ; v != v2.cend() ; ++v)
       cout << *v << endl;
   for (auto v = v3.cbegin() ; v != v3.cend() ; ++v)
       cout << *v << endl;
   @endcode
 *
 * If we organize the vector or the iterators in another vector, we can reduce the previous 4 loops
 * to two nested loops:
 *@code
  using IteratorType = SafeSTLVector<int>::const_iterator;

  SafeSTLVector<IteratorType> begins;
  begins.push_back(v0.cbegin());
  begins.push_back(v1.cbegin());
  begins.push_back(v2.cbegin());
  begins.push_back(v3.cbegin());


  SafeSTLVector<IteratorType> begins;
  ends.push_back(v0.cend());
  ends.push_back(v1.cend());
  ends.push_back(v2.cend());
  ends.push_back(v3.cend());

  for ( i = 0 ; i < 4 ; ++i)
     for (auto v = begins[i] ; v != ends[i] ; ++v)
        cout << *v << endl;
  @endcode
 *or, more nicely using a vector of iterators pairs:
 *@code
  using IteratorType = SafeSTLVector<int>::const_iterator;
  using IteratorPair = std::pair<IteratorType,IteratorType>;

  SafeSTLVector<IteratorPair> ranges;
  ranges.push_back(IteratorPair(v0.cbegin(),v0.cend()));
  ranges.push_back(IteratorPair(v1.cbegin(),v1.cend()));
  ranges.push_back(IteratorPair(v2.cbegin(),v2.cend()));
  ranges.push_back(IteratorPair(v3.cbegin(),v3.cend()));

  for ( const auto &r : ranges)
     for (auto v = r.first ; v != r.second ; ++v)
        cout << *v << endl;
  @endcode
 *
 * The ConcatenatedConstIterator class permits to iterate over the four vectors
 * (as made in the examples above) in this way:
 *@code
  using IteratorType = SafeSTLVector<int>::const_iterator;
  using ViewType = ConstView<IteratorType>;

  SafeSTLVector<IteratorView> ranges;
  ranges.push_back(ViewType(v0.cbegin(),v0.cend()));
  ranges.push_back(ViewType(v1.cbegin(),v1.cend()));
  ranges.push_back(ViewType(v2.cbegin(),v2.cend()));
  ranges.push_back(ViewType(v3.cbegin(),v3.cend()));

  ConcatenatedConstIterator<ViewType> begin(ranges,0); // this represents the first entry
  ConcatenatedConstIterator<ViewType> end(ranges,IteratorState::pass_the_end);  // this represents the one-pass-end entry

  for (; begin != end ; ++begin)
        cout << *begin << endl;

  @endcode
  @endcode
 * avoiding the nested loop needed by the
 * previous examples.
 *
 * @author M. Martinelli
 * @date 03 June 2014
 */
template <class ViewType,class ConstViewType>
class ConcatenatedConstIterator
  : public ConcatenatedIteratorData<ConstViewType,ConcatenatedConstIterator<ViewType,ConstViewType>>
{
public:
  /** @name Constructors & destructor */
  ///@{
  /**
   * Default constructor. It does nothing.
   */
  ConcatenatedConstIterator() = default;

  /**
   * Constructor.
   */
  ConcatenatedConstIterator(
    const SafeSTLVector<ConstViewType> &ranges,
    const Index index);

  /** Builds a ConcatenatedConstIterator from a ConcatenatedIterator.*/
  ConcatenatedConstIterator(const ConcatenatedIterator<ViewType> &it);


  /** Copy constructor. */
  ConcatenatedConstIterator(const ConcatenatedConstIterator<ViewType,ConstViewType> &it) = default;

  /** Move constructor. */
  ConcatenatedConstIterator(ConcatenatedConstIterator<ViewType,ConstViewType> &&it) = default;

  /** Destructor */
  ~ConcatenatedConstIterator() = default ;
  ///@}


  /** @name Assignment operators */
  ///@{
  /** Copy assignment operator. */
  ConcatenatedConstIterator<ViewType,ConstViewType> &operator=(
    const ConcatenatedConstIterator<ViewType,ConstViewType> &it) = default;

  /** Move assignment operator. */
  ConcatenatedConstIterator<ViewType,ConstViewType> &operator=(
    ConcatenatedConstIterator<ViewType,ConstViewType> &&it) = default;
  ///@}


  /** Returns a const reference to the entry with the identifier specified by @p id. */
  const typename ConstViewType::reference operator[](const Index id) const;
};



IGA_NAMESPACE_CLOSE


#include <igatools/utils/concatenated_iterator-inline.h>

#endif // #ifndef CONCATENATED_ITERATOR_H_
