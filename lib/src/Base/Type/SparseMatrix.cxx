//                                               -*- C++ -*-
/**
 *  @brief SparseMatrix implements Eigen sparse matrix class interface
 *
 *  Copyright 2005-2019 Airbus-EDF-IMACS-Phimeca
 *
 *  This library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include "openturns/SparseMatrix.hxx"

BEGIN_NAMESPACE_OPENTURNS

CLASSNAMEINIT(SparseMatrix)

/* Default constructor */
SparseMatrix::SparseMatrix()
  : Eigen::SparseMatrix<Scalar>()
{
  // Nothing to do
}


/* Constructor with dimensions */
SparseMatrix::SparseMatrix( const UnsignedInteger nbRows,
                            const UnsignedInteger nbCols)
  : Eigen::SparseMatrix<Scalar>(nbRows, nbCols)
{
  // Nothing to do
}

/* Virtual constructor */
SparseMatrix * SparseMatrix::clone() const
{
  return new SparseMatrix(*this);
}

/** Get the dimensions of the matrix */
/** Number of rows */
UnsignedInteger SparseMatrix::getNbRows() const
{
  return rows();
}

/** Number of columns */
UnsignedInteger SparseMatrix::getNbColumns() const
{
  return cols();
}


/** Sparse <> dense conversions */
Matrix SparseMatrix::asDenseMatrix()
{
  Matrix output(rows(), cols());
  
  for (int k=0; k<outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(*this,k); it; ++it)
    {
      UnsignedInteger i(it.row());
      UnsignedInteger j(it.col());
      output(i,j) = it.value();
    }

  return output;
}
    

END_NAMESPACE_OPENTURNS
