//                                               -*- C++ -*-
/**
 *  @brief The test file of class KarhunenLoeveP1Algorithm
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
#include "openturns/OT.hxx"
#include "openturns/OTtestcode.hxx"
#include <cmath>

using namespace OT;
using namespace OT::Test;

Scalar solveProblem(const UnsignedInteger size,
                    const String & format,
                    const String & solver,
                    const UnsignedInteger nbModes)
{ 
  ResourceMap::SetAsString("KarhunenLoeveP1Algorithm-CovarianceMatrixStorage", format);
  ResourceMap::SetAsString("KarhunenLoeveP1Algorithm-EigenvaluesSolver", solver);
  
  std::clock_t c_start = std::clock();
  Mesh mesh(IntervalMesher(Indices(1, size-1)).build(Interval(-1.0, 1.0)));
    
  AbsoluteExponential cov1D(Point(1, 1.0));
  KarhunenLoeveP1Algorithm algo(mesh, cov1D, 1e-3);

  algo.setNbModes(nbModes);
  algo.run();
  
  KarhunenLoeveResult result(algo.getResult());
  Point lambda(result.getEigenValues());
  ProcessSample KLModes(result.getModesAsProcessSample());
  GaussianProcess process(cov1D, KLModes.getMesh());
  Sample coefficients(result.project(process.getSample(10)));
  Basis KLFunctions(result.getModes());

  std::clock_t c_stop = std::clock();
  Scalar time_elapsed_ms = 1000.0 * (c_stop-c_start) / CLOCKS_PER_SEC;
  
  return time_elapsed_ms;
}

////////////////////////////////////

int main(int, char *[])
{
  TESTPREAMBLE;
  OStream fullprint(std::cout);
  setRandomGenerator();
  
  Description solvers(3);
  solvers[0] = "LAPACK";
  solvers[1] = "ARPACK";
  solvers[2] = "SPECTRA";

  Description formats(2);
  formats[0] = "LAPACK";
  formats[1] = "HMAT";
 
  Indices sizes(7);
  sizes[0] =    100;
  sizes[1] =    200;
  sizes[2] =    500;
  sizes[3] =   1000;
  sizes[4] =   2000;
  sizes[5] =  10000;
  sizes[6] = 100000;
  
  UnsignedInteger nev = 50;
  
  std::cout << "Size,Cov. mat. format,EV solver,Run time" << std::endl; 
  
  Scalar runtime = 0.0;
  
  for (UnsignedInteger i=0; i<sizes.getSize(); ++i)
    for (UnsignedInteger j=0; j<solvers.getSize(); ++j)
      for (UnsignedInteger k=0; k<formats.getSize(); ++k)
        {
          std::cout << sizes[i]    << ","
                    << formats[k] << ","
                    << solvers[j]   << ",";

          if (solvers[j] == "LAPACK" && (formats[k] == "HMAT" || sizes[i] > 1000))
            std::cout << "N/A" << std::endl << std::flush;
          else
          {
            runtime = solveProblem(sizes[i], formats[k], solvers[j], nev);
            std::cout << runtime << std::endl << std::flush;
          }
        }
        

  return ExitCode::Success;
}

