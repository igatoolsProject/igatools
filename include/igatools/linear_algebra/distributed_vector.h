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


#ifndef TRILINOS_VECTOR_H_
#define TRILINOS_VECTOR_H_

#include <igatools/base/config.h>
#include <igatools/base/logstream.h>
#include <igatools/utils/vector.h>

#include <igatools/linear_algebra/dense_vector.h>


#ifdef USE_TRILINOS
#include <igatools/linear_algebra/trilinos_tools.h>
#include <Tpetra_Vector.hpp>
#include <Epetra_Vector.h>
#endif

#ifdef USE_PETSC
#include <petscvec.h>
#endif


#include <memory>


IGA_NAMESPACE_OPEN

template < LAPack la_pack>
class Vector;




#ifdef USE_TRILINOS



/**
 *
 * @ingroup linear_algebra
 */
template<TrilinosImpl trilinos_impl>
class VectorTrilinos
{
public:
    /** Alias for the types of objects proper to the chosen TrilinosImpl (i.e. Tpetra or Epetra). */
    using Types = TrilinosTypes<trilinos_impl>;

    /** Type alias for the type of vector wrapped by the VectorTrilinos class. */
    using WrappedVector = typename Types::Vector;

    /** Type alias for the type of pointer to the matrix wrapped by the VectorTrilinos class. */
    using WrappedVectorPtr = typename Types::VectorPtr;

    /** Type alias for the pointer to the communicator. */
    using CommPtr = typename Types::CommPtr;

    /** Type alias for the pointer to the map. */
    using MapPtr = typename Types::MapPtr;

    /** Type alias for the TrilinosTool class that provides useful function for creating maps, graph, etc. */
    using Tools = TrilinosTools<trilinos_impl>;


    /**@name Constructor and destructor */
    ///@{
    /** Default constructor */
    VectorTrilinos() = delete;

    /**
     * Construct a distributed matrix with the dof distribution described by @p map.
     */
    VectorTrilinos(const MapPtr map);

    /**
     * Copy constructor. Not allowed to be used.
     */
    VectorTrilinos(const VectorTrilinos<trilinos_impl> &vector) = delete;

    /**
     * Move constructor.
     */
    VectorTrilinos(VectorTrilinos<trilinos_impl> &&vector) = default;

    /** Destructor */
    ~VectorTrilinos() = default;
    ///@}


    /** @name Assignment operators */
    ///@{
    /**
     * Copy assignment operator. Not allowed to be used.
     */
    VectorTrilinos<trilinos_impl> &operator=(const VectorTrilinos<trilinos_impl> &matrix) = delete;

    /**
     * Move assignment operator. Not allowed to be used.
     */
    VectorTrilinos<trilinos_impl> &operator=(VectorTrilinos<trilinos_impl> &&matrix) = delete;
    ///@}



    /** @name Methods for retrieving the Trilinos objects wrapped by this class */
    ///@{
    /**
     * Return the Trilinos RCP (Reference-Counted-Pointer, i.e. a smart pointer) wrapping the
     * concrete Trilinos distributed vector. Const version.
     */
    Teuchos::RCP<const WrappedVector> get_trilinos_vector() const ;

    /**
     * Return the Trilinos RCP (Reference-Counted-Pointer, i.e. a smart pointer) wrapping the
     * concrete Trilinos distributed vector. Non-const version.
     */
    Teuchos::RCP<WrappedVector> get_trilinos_vector() ;
    ///@}


protected:

    /**
     * Trilinos RCP (Reference-Counted-Pointer, i.e. a smart pointer) wrapping the
     * concrete Trilinos distributed matrix.
     */
    WrappedVectorPtr vector_ ;
};






/**
 * Numerical distributed Vector.
 * It's a wrapper to a Trilinos-Tpetra distributed vector.
 *
 * \author M. Martinelli 2013, 2014
 * \author pauletti 2013
 *
 * @ingroup linear_algebra
 */
template <>
class Vector<LAPack::trilinos_tpetra> : public VectorTrilinos<TrilinosImpl::tpetra>
{
public:
    using self_t = Vector<LAPack::trilinos_tpetra>;

    using Types = TrilinosTypes<TrilinosImpl::tpetra>;

    using CommPtr = typename Types::CommPtr;
    using MapPtr = typename Types::MapPtr;

    using LO = Index;
    using GO = Index;

    /** @name Constructor and destructor */
    ///@{
    /**
     * Default constructor not allowed.
     */
    Vector() = delete;

    /**
     * Construct a vector that gets a consecutive indexing
     * for @p n degrees of freedom.
     * All entries are set to zero.
     */
    Vector(const Index n,CommPtr comm = Teuchos::createSerialComm<int>());

    /**
     * Construct a vector that gets a non consecutive indexing
     * for @p dof_ids degrees of freedom. An example could be the usage of
     * dof numbering provided from some external library.
     * All entries are set to zero.
     */
    Vector(const vector<Index> &dof_ids,CommPtr comm = Teuchos::createSerialComm<int>());

    Vector(const MapPtr map);

    /**
     * Copy constructor. Performs a shallow copy of the object (i.e.)
     * only the smart pointers are copied, not the pointed objects.
     */
    Vector(const self_t &v) = default;

    /** Move constructor. */
    Vector(self_t &&v) = default;


    /** Destructor */
    ~Vector() = default;
    ///@}



    /** @name Assignment operators */
    ///@{

    /** Copy assignment operator */
    self_t &operator=(const self_t &v) = default;

    /** Move assignment operator */
    self_t &operator=(self_t &&v) = default;

    ///@}


    /**
     * @name Function for creating Vector objects
     * @note These methods implement the "create idiom" and return
     * a Vector object wrapped by a std::shared_ptr
     */
    ///@{

    /**
     * Create a vector that gets an consecutive indexing
     * for @p n degrees of freedom.
     * Initializing all entries to zero.
     */
    static std::shared_ptr<self_t> create(const Index n);

    /**
     * Create a vector that gets an non consecutive indexing
     * for @p dof_ids degrees of freedom. An example could be the usage of
     * dof numbering provided from some external library, IRIT as an example.
     * Initializing all entries to zero.
     */
    static std::shared_ptr<self_t> create(const vector<Index> &dof_ids);
    ///@}



    /** @name Methods for getting and/or modifying the vector entries */
    ///@{

    /**
     * Add the value @p input to the vector entry @p global_id.
     * @note @p global_id is the global index of the entry.
     */
    void add_entry(const Index global_id, const Real input);


    /**
     * \brief This function add the @p local_vector values to the global vector.
     * The local-to-global ids are passed as input argument.
     * \param local_to_global The vector containing the global ids associated to the local vector entries
     * \param local_vector The local vector that must be added to the global vector.
     * \note The size of the local vector must be equal to the dimension of the
     * local_to_global vector, otherwise an exception will be raised.
     * \author M. Martinelli
     * \date 29 Jan 2013
     */
    void add_block(
        const vector< Index > &local_to_global,
        const DenseVector &local_vector);


    /** Add vec onto current vector */
    self_t &operator+=(const self_t &vec);

    /**
     * Update vector values with scaled values of A, this = scalar_this*this + scalar_A * A.
     */
    self_t &update(const Real scalar_A, const self_t &A, const Real scalar_this);


    /** norm2 of the vector */
    Real norm2() const;

    /** Sets all the entries to zero. */
    void clear();

    /**
     * Computes and returns the dot product between *this and A.
     */
    Real dot(const self_t &A) const;

    /**
     * Returns the const reference to the entry of the vector identified by the @p global_id
     */
    const Real &operator()(const Index global_id) const;

    /**
     * Returns the reference to the entry of the vector identified by the @p global_id
     */
    Real &operator()(const Index global_id);

    /** Returns the number of entries in the vector. */
    Index size() const;

    ///@}



    /**
     * Returns the local coefficients of the distributed vector,
     * from the vector of local-to-global indices.
     */
    vector<Real>
    get_local_coefs(const vector<Index> &local_to_global_ids) const;


    // TODO (pauletti, Nov 27, 2014): Document this
    vector<Real> get_as_vector() const;



    /**
     * Print the content of the vector, mostly for debug purposes.
     * @param out
     */
    void print_info(LogStream &out) const;
};



/**
 * Numerical distributed Vector.
 * It's a wrapper to a Trilinos-Epetra distributed vector.
 *
 * \author M. Martinelli 2014
 *
 * @ingroup linear_algebra
 */
template <>
class Vector<LAPack::trilinos_epetra> : public VectorTrilinos<TrilinosImpl::epetra>
{
public:
    using self_t = Vector<LAPack::trilinos_epetra>;

    using Types = TrilinosTypes<TrilinosImpl::epetra>;

    using CommPtr = typename Types::CommPtr;
    using MapPtr = typename Types::MapPtr;
    using Map = typename Types::Map;


    /** @name Constructor and destructor */
    ///@{
    /**
     * Default constructor not allowed.
     */
    Vector() = delete;

    /**
     * Construct a vector that gets a consecutive indexing
     * for @p n degrees of freedom.
     * All entries are set to zero.
     */
    Vector(const Index n,CommPtr comm = Teuchos::rcp(new Epetra_SerialComm()));

    /**
     * Construct a vector that gets a non consecutive indexing
     * for @p dof_ids degrees of freedom. An example could be the usage of
     * dof numbering provided from some external library.
     * All entries are set to zero.
     */
    Vector(const vector<Index> &dof_ids,CommPtr comm = Teuchos::rcp(new Epetra_SerialComm()));

    Vector(const MapPtr map);

    /**
     * Copy constructor. Performs a shallow copy of the object (i.e.)
     * only the smart pointers are copied, not the pointed objects.
     */
    Vector(const self_t &v) = default;

    /** Move constructor. */
    Vector(self_t &&v) = default;


    /** Destructor */
    ~Vector() = default;
    ///@}



    /** @name Assignment operators */
    ///@{

    /** Copy assignment operator */
    self_t &operator=(const self_t &v) = default;

    /** Move assignment operator */
    self_t &operator=(self_t &&v) = default;

    ///@}


    /**
     * @name Function for creating Vector objects
     * @note These methods implement the "create idiom" and return
     * a Vector object wrapped by a std::shared_ptr
     */
    ///@{

    /**
     * Create a vector that gets an consecutive indexing
     * for @p n degrees of freedom.
     * Initializing all entries to zero.
     */
    static std::shared_ptr<self_t> create(const Index n);

    /**
     * Create a vector that gets an non consecutive indexing
     * for @p dof_ids degrees of freedom. An example could be the usage of
     * dof numbering provided from some external library, IRIT as an example.
     * Initializing all entries to zero.
     */
    static std::shared_ptr<self_t> create(const vector<Index> &dof_ids);
    ///@}



    /** @name Methods for getting and/or modifying the vector entries */
    ///@{

    /**
     * Add the value @p input to the vector entry @p global_id.
     * @note @p global_id is the global index of the entry.
     */
    void add_entry(const Index global_id, const Real input);


    /**
     * \brief This function add the @p local_vector values to the global vector.
     * The local-to-global ids are passed as input argument.
     * \param local_to_global The vector containing the global ids associated to the local vector entries
     * \param local_vector The local vector that must be added to the global vector.
     * \note The size of the local vector must be equal to the dimension of the
     * local_to_global vector, otherwise an exception will be raised.
     * \author M. Martinelli
     * \date 29 Jan 2013
     */
    void add_block(
        const vector< Index > &local_to_global,
        const DenseVector &local_vector);

    /** Add vec onto current vector */
    self_t &operator+=(const self_t &vec);


    /**
     * Update vector values with scaled values of A, this = scalar_this*this + scalar_A * A.
     */
    self_t &update(const Real scalar_A, const self_t &A, const Real scalar_this);

    /**
     * Returns the const reference to the entry of the vector identified by the @p global_id
     */
    const Real &operator()(const Index global_id) const;

    /**
     * Returns the reference to the entry of the vector identified by the @p global_id
     */
    Real &operator()(const Index global_id);

    /** Returns the number of entries in the vector. */
    Index size() const;
    ///@}

    /** norm2 of the vector */
    Real norm2() const;

    /**
     * Computes and returns the dot product between *this and A.
     */
    Real dot(const self_t &A) const;

    /** Sets all the entries to zero. */
    void clear();

    /**
     * Returns the local coefficients of the distributed vector,
     * from the vector of local-to-global indices.
     */
    vector<Real>
    get_local_coefs(const vector<Index> &local_to_global_ids) const;


    /**
     * Print the content of the vector, mostly for debug purposes.
     * @param out
     */
    void print_info(LogStream &out) const;

    /** Return the vector as the type use for coefficient of an IgFunction */
    vector<Real> get_as_vector() const;

};




#endif //#ifdef USE_TRILINOS



#ifdef USE_PETSC

/**
 * Numerical distributed Vector.
 * It's a wrapper to a PETSc distributed vector.
 *
 * \author M. Martinelli 2013, 2014
 *
 */
template <>
class Vector<LAPack::petsc>
{
public:
    using self_t = Vector<LAPack::petsc>;


    /** @name Constructor and destructor */
    ///@{
    /**
     * Default constructor not allowed.
     */
    Vector() = delete;

    /**
     * Construct a vector that gets an consecutive indexing
     * for @p dof_ids degrees of freedom.
     * All entries are set to zero.
     */
    Vector(const Index size);

    /**
     * Construct a vector that gets an non consecutive indexing
     * for @p dof_ids degrees of freedom. An example could be the usage of
     * dof numbering provided from some external library.
     * All entries are set to zero.
     */
    Vector(const vector< Index > &dof_ids);

    /**
     * Copy constructor. Performs a shallow copy of the object (i.e.)
     * only the smart pointers are copied, not the pointed objects.
     */
    Vector(const self_t &v) = default;

    /** Move constructor. */
    Vector(self_t &&v) = default;


    /** Destructor */
    ~Vector() = default;
    ///@}



    /** @name Assignment operators */
    ///@{

    /** Copy assignment operator */
    self_t &operator=(const self_t &v) = default;

    /** Move assignment operator */
    self_t &operator=(self_t &&v) = default;

    ///@}


    /**
     * @name Function for creating Vector objects
     * @note These methods implement the "create idiom" and return
     * a Vector object wrapped by a std::shared_ptr
     */
    ///@{

    /**
     * Create a vector that gets an consecutive indexing
     * for @p dof_ids degrees of freedom.
     * Initializing all entries to zero.
     */
    static std::shared_ptr<self_t> create(const Index size);

    /**
     * Create a vector that gets an non consecutive indexing
     * for @p dof_ids degrees of freedom. An example could be the usage of
     * dof numbering provided from some external library, IRIT as an example.
     * Initializing all entries to zero.
     */
    static std::shared_ptr<self_t> create(const vector<Index> &dof_ids);


    /**
     * Create a distributed vector with its index distribution
     * specified by the Space @p space.
     */
    template<class Space>
    static inline std::shared_ptr<self_t> create(const Space &space);
    ///@}



    /** @name Methods for getting and/or modifying the vector entries */
    ///@{

    /**
     * Add the value @p input to the vector entry (i).
     * @note @p i is the global index of the entry (i).
     */
    void add_entry(const Index i, const Real input);


    /**
     * \brief This function add the @p local_vector values to the global vector.
     * The local-to-global ids are passed as input argument.
     * \param local_to_global The vector containing the global ids associated to the local vector entries
     * \param local_vector The local vector that must be added to the global vector.
     * \note The size of the local vector must be equal to the dimension of the
     * local_to_global vector, otherwise an exception will be raised.
     * \author M. Martinelli
     * \date 29 Jan 2013
     */
    void add_block(
        const vector< Index > &local_to_global,
        const DenseVector &local_vector);



    /**
     * Returns the const reference to the i-th entry of the vector.
     * @note @p i is the global index of the i-th entry.
     */
    const Real &operator()(const Index global_id) const;

    /**
     * Returns the reference to the i-th entry of the vector.
     * @note @p i is the global index of the i-th entry.
     */
    Real &operator()(const Index global_id);

    /** Returns the number of entries in the vector. */
    Index size() const;

    ///@}



    /**
     * Returns the local coefficients of the distributed vector,
     * from the vector of local-to-global indices.
     */
    vector<Real>
    get_local_coefs(const vector<Index> &local_to_global_ids) const;


    /**
     * Print the content of the vector, mostly for debug purposes.
     * @param out
     */
    void print_info(LogStream &out) const;



public:

    /** Type of Petsc object wrapped by this class. */
//    using WrappedVectorType = Tpetra::MultiVector<Real,Index,Index>;



    /** @name Methods for retrieving the PETSc objects wrapped by this class */
    ///@{
    /**
     * todo: add documentation
     */
    Vec get_petsc_vector() const;

    /**
     * todo: add documentation
     */
    Vec get_petsc_vector();
    ///@}

private:
    MPI_Comm comm_;
    Vec vector_;

    Real real_tmp_;
};


#endif //#ifdef USE_PETSC



IGA_NAMESPACE_CLOSE

#endif /* TRILINOS_VECTOR_H_ */
