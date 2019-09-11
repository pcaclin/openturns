//                                               -*- C++ -*-
/**
 *  @brief The class CovarianceMatrix implements blank free samples
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
#include "openturns/CovarianceMatrix.hxx"

BEGIN_NAMESPACE_OPENTURNS

CLASSNAMEINIT(CovarianceMatrix)

/* Default constructor */
CovarianceMatrix::CovarianceMatrix()
  : SymmetricMatrix(0)
{
  // Nothing to do
}

/* Constructor with implementation */
CovarianceMatrix::CovarianceMatrix(const Implementation & i)
  : SymmetricMatrix(i)
{
  // Nothing to do
}

/* Constructor with implementation */
CovarianceMatrix::CovarianceMatrix(const MatrixImplementation & i)
  : SymmetricMatrix(i)
{
  // Nothing to do
}

/* Constructor with size (dim, which is the same for nbRows_ and nbColumns_ )*/
/* CovarianceMatrix is initialized to identity matrix */
CovarianceMatrix::CovarianceMatrix(const UnsignedInteger dim)
  : SymmetricMatrix(MatrixImplementation::identityMatrix(dim))
{
  // Nothing to do
}

/* Constructor from external collection */
/* If the dimensions of the matrix and of the collection */
/* do not match, either the collection is truncated */
/* or the rest of the matrix is filled with zeros */
CovarianceMatrix::CovarianceMatrix(const UnsignedInteger dim,
                                   const Collection<Scalar> &elementsValues)
  : SymmetricMatrix(dim, elementsValues)
{
  // Nothing to do
}

/** Create identity matrix as CovarianceMatrix */
CovarianceMatrix CovarianceMatrix::identityMatrix(const UnsignedInteger dimension)
{
  return CovarianceMatrix(MatrixImplementation::identityMatrix(dimension));
}


/* String converter */
String CovarianceMatrix::__repr__() const
{
  checkSymmetry();
  return OSS() << "class=" << getClassName()
         << " dimension=" << this->getDimension()
         << " implementation=" << getImplementation()->__repr__();
}

/* CovarianceMatrix transpose */
CovarianceMatrix CovarianceMatrix::transpose () const
{
  return *this;
}

/* CovarianceMatrix addition (must have the same dimensions) */
CovarianceMatrix CovarianceMatrix::operator + (const CovarianceMatrix & m) const
{
  return Implementation((*getImplementation() + * (m.getImplementation()) ).clone());
}

/* Check if the matrix is SPD */
Bool CovarianceMatrix::isPositiveDefinite() const
{
  return getImplementation()->isPositiveDefinite();
}

/* Build the Cholesky factorization of the matrix */
TriangularMatrix CovarianceMatrix::computeCholesky(const Bool keepIntact)
{
  return Implementation(getImplementation()->computeCholesky(keepIntact).clone());
}


/* Resolution of a linear system */
Point CovarianceMatrix::solveLinearSystem(const Point & b,
    const Bool keepIntact)
{
  return getImplementation()->solveLinearSystemCov(b, keepIntact);
}

Matrix CovarianceMatrix::solveLinearSystem(const Matrix & b,
    const Bool keepIntact)
{
  return Implementation(getImplementation()->solveLinearSystemCov(*b.getImplementation(), keepIntact).clone());
}

END_NAMESPACE_OPENTURNS
