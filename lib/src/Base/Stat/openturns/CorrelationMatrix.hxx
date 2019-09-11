//                                               -*- C++ -*-
/**
 *  @brief The class CorrelationMatrix implements correlation matrices
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
#ifndef OPENTURNS_CORRELATIONMATRIX_HXX
#define OPENTURNS_CORRELATIONMATRIX_HXX

#include "openturns/OTprivate.hxx"
#include "openturns/CovarianceMatrix.hxx"

BEGIN_NAMESPACE_OPENTURNS



/**
 * @class CorrelationMatrix
 */

class OT_API CorrelationMatrix
  : public CovarianceMatrix
{
  CLASSNAME

public:
  typedef Collection<Scalar>       ScalarCollection;


  /** Default constructor */
  CorrelationMatrix();

  /** Constructor with implementation */
  CorrelationMatrix(const Implementation & i);

  /** Constructor with implementation */
  CorrelationMatrix(const MatrixImplementation & i);

  /** Constructor with size */
  explicit CorrelationMatrix(const UnsignedInteger dim);

  /** Constructor from external collection */
  CorrelationMatrix(const UnsignedInteger dim,
                    const ScalarCollection & elementsValues);

  /** String converter */
  virtual String __repr__() const;

  /** CorrelationMatrix transpose */
  CorrelationMatrix transpose () const;
  
  /** Create identity matrix as CorrelationMatrix */
  static CorrelationMatrix identityMatrix(const UnsignedInteger dimension);

  /** Multiplications */ //TODO
#ifdef _MSC_VER   // VS2010 does not like 'using' being called after overloads
  using CovarianceMatrix::operator *;
#endif

  Matrix operator * (const Matrix & m) const;
  SquareMatrix operator * (const SquareMatrix & m) const;
  SymmetricMatrix operator * (const SymmetricMatrix & m) const;
  CovarianceMatrix operator * (const CovarianceMatrix & m) const;
  CorrelationMatrix operator * (const CorrelationMatrix & m) const;
#ifndef _MSC_VER
  using CovarianceMatrix::operator *;
#endif
  
  
protected:


private:

}; /* class CorrelationMatrix */

END_NAMESPACE_OPENTURNS

#endif /* OPENTURNS_COVARIANCEMATRIX_HXX */
