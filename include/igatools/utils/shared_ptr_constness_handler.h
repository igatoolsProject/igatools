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

#ifndef SHARED_PTR_CONSTNESS_HANDLER_H_
#define SHARED_PTR_CONSTNESS_HANDLER_H_

#include <igatools/base/config.h>
#include <igatools/base/exceptions.h>
#include <igatools/base/logstream.h>




IGA_NAMESPACE_OPEN

/*
namespace { // detail
    template <typename T>
    struct implicit_convert : boost::static_visitor<T> {
        template <typename U> T operator()(U&& u) const { return std::forward<U>(u); }
    };
}
//*/


/*
Ref getVal(std::string& name) {
    return boost::apply_visitor(implicit_convert<Ref>(), map[name]);
}
//*/

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
     * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     * mechanism.
     */
    SharedPtrConstnessHandler() = default;

    /**
     * Constructs the object using a shared pointer to non-const data.
     *
     * @note In DEBUG mode is performed a check about the nullity of the input data:
     * an assertion will be raised if the input data is a nullptr.
     */
    SharedPtrConstnessHandler(const Ptr &data)
        :
        data_is_const_(false)
    {
        Assert(data != nullptr, ExcNullPtr());

        data_.ptr_ = data;
    };

    /**
     * Constructs the object using a shared pointer to const data.
     */
    SharedPtrConstnessHandler(const PtrToConst &data)
        :
        data_is_const_(true)
    {
        Assert(data != nullptr, ExcNullPtr());
        data_.ptr_to_const_ = data;
    };

    SharedPtrConstnessHandler(Ptr &&data)
        :
        data_is_const_(false)
    {
        data_.ptr_(std::forward<Ptr>(data));
        Assert(data_.ptr_ != nullptr, ExcNullPtr());
    };

    SharedPtrConstnessHandler(PtrToConst &&data)
        :
        data_is_const_(true)
    {
        data_.ptr_to_const_(std::forward<Ptr>(data));
        Assert(data_.ptr_to_const_ != nullptr, ExcNullPtr());
    };

    /**
     * Copy constructor;
     */
    SharedPtrConstnessHandler(const SharedPtrConstnessHandler<T> &obj) = delete;

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
    SharedPtrConstnessHandler<T> &operator=(const SharedPtrConstnessHandler<T> &obj) = delete;

    /**
     * Move assignment operator. Not allowed to be used.
     */
    SharedPtrConstnessHandler<T> &operator=(SharedPtrConstnessHandler<T> &&obj) = delete;
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
    std::shared_ptr<T> get_ptr_data()
    {
        Assert(!data_is_const_,ExcMessage("Data is built as const."));
        return data_.ptr_;
    }

    /**
     * Return the non-const reference of the shared pointer that is pointing to the non-const data.
     */
    std::shared_ptr<T> &get_ref_ptr_data()
    {
        Assert(!data_is_const_,ExcMessage("Data is built as const."));
        return data_.ptr_;
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
        if (data_is_const_)
            return data_.ptr_to_const_;
        else
            return data_.ptr_;
    }

    /**
     * Return the non-const reference of the shared pointer that is pointing to the const data.
     */
    std::shared_ptr<const T> &get_ref_ptr_const_data() const
    {
        if (data_is_const_)
            return data_.ptr_to_const_;
        else
            return data_.ptr_;
    }

    void print_info(LogStream &out) const
    {
        if (data_is_const_)
        {
            out.begin_item("Pointer to const data");
            data_.ptr_to_const_->print_info(out);
        }
        else
        {
            out.begin_item("Pointer to non-const data");
            data_.ptr_->print_info(out);
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
     * Pointer operator. Return a const pointer to the const data.
     */
    const T *const operator->() const
    {
        return this->get_ptr_const_data().get();
    }

    /**
     * Check if unique.
     *
     * Returns whether the shared_ptr object does not share ownership over its pointer with other shared_ptr objects (i.e., it is unique).
     */
    bool unique() const
    {
        if (data_is_const_)
            return data_.ptr_to_const_.unique();
        else
            return data_.ptr_.unique();
    }

#if 0
    /**
     * Dereference operator. Return a non-const reference of the data.
     */
    T &operator*()
    {
        return this->get_ref_data();
    }

    /**
     * Pointer operator. Return a const pointer to the non-const data.
     */
    T *const operator->()
    {
        return this->get_ptr_data().get();
    }
#endif

private:
//    boost::variant<Ptr,PtrToConst> ptr_to_data_;


    struct DataConstnessVariants
    {
        Ptr ptr_;
        PtrToConst ptr_to_const_;
    } data_;


    bool data_is_const_ = false;



#ifdef SERIALIZATION
    /**
     * @name Functions needed for boost::serialization
     * @see <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     */
    ///@{
    friend class boost::serialization::access;

    template<class Archive>
    void
    serialize(Archive &ar, const unsigned int version)
    {
        ar &boost::serialization::make_nvp("data_is_const_",data_is_const_);


        // In order to serialize the data, we need to cast it to non-const
        Ptr tmp;
        if (data_is_const_)
            tmp = std::const_pointer_cast<T>(data_.ptr_to_const_);
        else
            tmp = data_.ptr_;

        ar &boost::serialization::make_nvp("tmp_ptr_to_data_",tmp);
        Assert(tmp != nullptr,ExcNullPtr());

        // After deserialization we need to cast the data to the correct constness
        if (data_is_const_)
            data_.ptr_to_const_ = std::const_pointer_cast<const T>(tmp);
        else
            data_.ptr_ = tmp;

    }
    ///@}
#endif // SERIALIZATION

};


IGA_NAMESPACE_CLOSE


#endif // SHARED_PTR_CONSTNESS_HANDLER_H_
