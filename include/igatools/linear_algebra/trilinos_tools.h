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

#ifndef TRILINOS_TOOLS_H_
#define TRILINOS_TOOLS_H_

#include <igatools/base/config.h>

#include <igatools/basis_functions/space_manager.h>
#include <igatools/linear_algebra/trilinos_types.h>


#ifdef USE_TRILINOS



IGA_NAMESPACE_OPEN


/**
 * @brief This class provides some methods for easy construction of Trilinos objects.
 */
template <TrilinosImpl trilinos_impl>
class TrilinosTools
{
public:
    using Types = TrilinosTypes<trilinos_impl>;

    using CommPtr = typename Types::CommPtr;
    using MapPtr = typename Types::MapPtr;
    using Graph = typename Types::Graph;
    using GraphPtr = typename Types::GraphPtr;
    using Matrix = typename Types::Matrix;
    using MatrixPtr = typename Types::MatrixPtr;
    using Vector = typename Types::Vector;
    using VectorPtr = typename Types::VectorPtr;
    using OpPtr = typename Types::OpPtr;

    static MapPtr build_row_map(const SpaceManager &space_manager, const CommPtr comm);

    static MapPtr build_col_map(const SpaceManager &space_manager, const CommPtr comm);

    static GraphPtr build_graph(const SpaceManager &space_manager,const MapPtr row_map,const MapPtr col_map);

    static MatrixPtr build_matrix(const GraphPtr graph);

    static VectorPtr build_vector(const MapPtr map);

    static OpPtr compute_preconditioner_ILU(const MatrixPtr matrix, const int fill_level = 0);


private:

    /**
     *
     * The ILUKPreconditionerFactory class encapsulates creation of an Ifpack or Ifpack2
     * preconditioner (returned as an OpPtr) from a MatrixPtr.
     * Change the
     * preconditioner parameters by changing the function that creates a
     * ParameterList for Ifpack/Ifpack2.
     * See the Ifpack/Ifpack2 documentation for a
     * list and description of the parameters that it accepts.
     *
    */
    class ILUKPreconditionerFactory
    {
    public:
        // The constructor doesn't do anything, since this factory doesn't
        // keep any state.
        ILUKPreconditionerFactory() {};

        /**
         * Return a ParameterList for asking Ifpack or Ifpack2 to create an RILUK
         * incomplete factorization preconditioner with a certain @p fill_level, drop
         * tolerance 0.0, and absolute threshold 0.0
         */
        Teuchos::RCP<Teuchos::ParameterList>
        get_parameter_list_for_preconditioner_creation(const int fill_level) const;


        // Compute and return a preconditioner.
        OpPtr
        create(const Teuchos::RCP<const Matrix> &A,
               const Teuchos::RCP<Teuchos::ParameterList> &plist) const;

    };
};



IGA_NAMESPACE_CLOSE

#endif // #ifdef USE_TRILINOS

#endif // #ifndef TRILINOS_TOOLS_H_
