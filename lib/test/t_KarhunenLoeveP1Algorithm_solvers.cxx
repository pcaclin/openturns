//                                               -*- C++ -*-
/**
 *  @brief The test file of class KarhunenLoeveP1Algorithm
 *
 *  Copyright 2005-2019 Airbus-EDF-IMACS-ONERA-Phimeca
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
#include "openturns/OT.hxx"
#include "openturns/OTtestcode.hxx"
#include <cmath>

using namespace OT;
using namespace OT::Test;

KarhunenLoeveResult testConfiguration(const String & covarianceMatrixStorage,
                                      const String & eigenvaluesSolver,
                                      const String & dim)
{
  // Set KarhunenLoeveP1Algorithm parameters
  ResourceMap::SetAsString("KarhunenLoeveP1Algorithm-CovarianceMatrixStorage", covarianceMatrixStorage);
  ResourceMap::SetAsString("KarhunenLoeveP1Algorithm-EigenvaluesSolver", eigenvaluesSolver);
  
  // Prepare and run KarhunenLoeveP1Algorithm
  Mesh mesh(IntervalMesher(Indices(1, 9)).build(Interval(-1.0, 1.0)));
  KarhunenLoeveP1Algorithm algo;
  if (dim == "2D")
  {
    CorrelationMatrix R(2);
    R(0, 1) = 0.5;
    Point scale(1, 1.0);
    Point amplitude(2);
    amplitude[0] = 1.0;
    amplitude[1] = 2.0;
    ExponentialModel cov2D(scale, amplitude, R);
    algo = KarhunenLoeveP1Algorithm(mesh, cov2D);
  }
  else if (dim == "1D")
  {
    AbsoluteExponential cov1D(Point(1, 1.0));
    algo = KarhunenLoeveP1Algorithm(mesh, cov1D);
  }
  else
    throw InternalException(HERE) << "invalid dimension: " << dim;
  
  // Solve problem
  algo.setThreshold(1e-3);
  algo.run();
  
  return algo.getResult();
}
  
Bool testPointVsReference(const Point & testPoint,
                          const Point & refPoint,
                          const Scalar tol = 1e-6)
{
  Point resizedRefPoint(testPoint.getSize());
  
  if (testPoint.getSize() > refPoint.getSize())
    return false;
  else
    std::copy(refPoint.begin(), refPoint.begin()+testPoint.getSize(), resizedRefPoint.begin());
  
  return (testPoint - resizedRefPoint).norm() < tol;
}
  
int main(int, char *[])
{
  TESTPREAMBLE;
  OStream fullprint(std::cout);
  setRandomGenerator();
  
  Bool success = true;
  
  Description eigenValuesSolvers(3);
  eigenValuesSolvers[0] = "LAPACK";
  eigenValuesSolvers[1] = "ARPACK";
  eigenValuesSolvers[2] = "DUMMY_EV_SOLVER";
  
  Description covarianceMatrixStorageFormats(3);
  covarianceMatrixStorageFormats[0] = "LAPACK";
  covarianceMatrixStorageFormats[1] = "HMAT";
  covarianceMatrixStorageFormats[2] = "DUMMY_STORAGE_FORMAT";
  
  Description dimension(2);
  dimension[0] = "1D";
  dimension[1] = "2D";
  
  KarhunenLoeveResult resultReference1D(testConfiguration("LAPACK","LAPACK","1D"));
  KarhunenLoeveResult resultReference2D(testConfiguration("LAPACK","LAPACK","2D"));
  
  for (UnsignedInteger i=0; i<covarianceMatrixStorageFormats.getSize(); ++i)
    for (UnsignedInteger j=0; j<eigenValuesSolvers.getSize(); ++j)
      for (UnsignedInteger k=0; k<dimension.getSize(); ++k)
      {
        fullprint << "# TESTING CONFIGURATION: covarianceMatrixStorage=" << covarianceMatrixStorageFormats[i] <<
        ", eigenValuesSolver=" << eigenValuesSolvers[j] << ", dimension=" << dimension[k] << std::endl;
        try
        {        
          // Test configuration
          KarhunenLoeveResult result;
          result = testConfiguration(covarianceMatrixStorageFormats[i], eigenValuesSolvers[j], dimension[k]);
          
          // Computation successful
          fullprint << "  - Computation successfully completed, eigenValues=" << result.getEigenValues().__str__() << std::endl;
          
          // Retrieving reference eigenvalues
          Point eigenValuesReference;
          if (dimension[k] == "1D")
            eigenValuesReference = resultReference1D.getEigenValues();
          else
            eigenValuesReference = resultReference2D.getEigenValues();
          
          // Comparing computed eigenvalues to reference
          if (testPointVsReference(result.getEigenValues(), eigenValuesReference))
            fullprint << "  => Configuration passed";
          else
            fullprint << "  => Configuration failed";
          
          fullprint << std::endl << std::endl;
        }
        catch (InternalException & ex)
        {
          fullprint << "  - " << ex << std::endl;
          fullprint << "  => Configuration failed" << std::endl << std::endl;
        }  
      }


}

