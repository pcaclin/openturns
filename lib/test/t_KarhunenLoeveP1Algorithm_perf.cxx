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
                    const Scalar threshold,
                    const Scalar nevRatio,
                    const Scalar ncvRatio)
{  
  std::clock_t c_start = std::clock();
  Mesh mesh(IntervalMesher(Indices(1, size-1)).build(Interval(-1.0, 1.0)));
    
  AbsoluteExponential cov1D(Point(1, 1.0));
  KarhunenLoeveP1Algorithm algo(mesh, cov1D, threshold);

  algo.runWithParameters(nevRatio, ncvRatio);
  
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

  Indices sizes(5);
  sizes[0] =   100;
  sizes[1] =   200;
  sizes[2] =   500;
  sizes[3] =  1000;
  sizes[4] =  2000;
  
  Point nevRatios(10);
  nevRatios[0] = 0.01;
  nevRatios[1] = 0.02;
  nevRatios[2] = 0.05;
  nevRatios[3] = 0.10;
  nevRatios[4] = 0.15;
  nevRatios[5] = 0.20;
  nevRatios[6] = 0.25;
  nevRatios[7] = 0.30;
  nevRatios[8] = 0.50;
  nevRatios[9] = 1.00;
  
  Point ncvRatios(7);
  ncvRatios[0] =  1.0;
  ncvRatios[1] =  1.2;
  ncvRatios[2] =  1.5;
  ncvRatios[3] =  2.0;
  ncvRatios[4] =  3.0;
  ncvRatios[5] =  5.0;
  ncvRatios[6] = 10.0;
  
  Point thresholds(4);
  thresholds[0] = 1e-1;
  thresholds[1] = 1e-2;
  thresholds[2] = 1e-3;
  thresholds[3] = 1e-4;
  
  std::cout << "Size,NEV ratio,NCV ratio,Threshold,nev,ncv,Selected EV,Cumulated variance,Computed variance,Selected variance,Comp. var. ratio,Sel. var. ratio,Run time" << std::endl; 
  
  Scalar runtime = 0.0;
  
  for (UnsignedInteger i=0; i<sizes.getSize(); ++i)
    for (UnsignedInteger j=0; j<nevRatios.getSize(); ++j)
      for (UnsignedInteger k=0; k<ncvRatios.getSize(); ++k)
        for (UnsignedInteger l=0; l<thresholds.getSize(); ++l)
        {
          std::cout << sizes[i]    << "," 
                    << nevRatios[j]   << ","
                    << ncvRatios[k]   << ","
                    << thresholds[l] << ",";
                    
          runtime = solveProblem(sizes[i], thresholds[l], nevRatios[j], ncvRatios[k]);
          
          std::cout << runtime << std::endl << std::flush;
        }
        

  return ExitCode::Success;
}

