//                                               -*- C++ -*-
/**
 *  @brief This class implements the computation of the Karhunen-Loeve
 *         basis and eigenvalues of a given covariance model based on
 *         P1 Lagrange approximation.
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
#include "openturns/KarhunenLoeveP1Algorithm.hxx"
#include "openturns/SquareComplexMatrix.hxx"
#include "openturns/P1LagrangeEvaluation.hxx"
#include "openturns/PiecewiseLinearEvaluation.hxx"
#include "openturns/PersistentObjectFactory.hxx"
#include "openturns/SpecFunc.hxx"
#include "openturns/SparseMatrix.hxx"
#include <algorithm>
#include <arpack.hpp>


BEGIN_NAMESPACE_OPENTURNS

typedef Collection<Complex> ComplexCollection;
typedef Eigen::Triplet<Scalar>  Triplet;
typedef std::vector<Triplet>    TripletVector;


// Nothing

/** Defining class KLGenMatProd **/
class KLGenMatProd
{
  public:
    
  KLGenMatProd( const CovarianceMatrix & C,
                const SparseMatrix & G)
  : C_(C)
  , G_(G)
  , rows_(C_.getNbRows())
  , cols_(C_.getNbColumns())
  {    
    // Nothing to do
  }
    
  int rows()
  {
    return rows_;
  }
  
  int cols()
  {
    return cols_;
  }
  
  void av(const Scalar * x_in, Scalar * y_out)
  {
    // Convert double array to Eigen::VectorXd
    Eigen::VectorXd v(rows_);
    std::copy(x_in, x_in+rows_, v.data());
    
    // Compute product with sparse matrix G       
    Eigen::VectorXd w(G_*v);

    // Compute product with dense matrix C
    Point OTw(rows_);
    std::copy(w.data(), w.data()+rows_, OTw.begin());
    
    Point OTv(C_*OTw);
        
    // Output double array
    std::copy(OTv.begin(), OTv.end(), y_out);
  }
    
  private:
    CovarianceMatrix C_;
    SparseMatrix G_;
    int rows_;
    int cols_;
};


/**
 * @class KarhunenLoeveP1Algorithm
 */

CLASSNAMEINIT(KarhunenLoeveP1Algorithm)

static const Factory<KarhunenLoeveP1Algorithm> Algorithm_KarhunenLoeveP1Algorithm;

/* Constructor without parameters */
KarhunenLoeveP1Algorithm::KarhunenLoeveP1Algorithm()
  : KarhunenLoeveAlgorithmImplementation()
  , mesh_()
{
  // Nothing to do
}

/* Constructor with parameters */
KarhunenLoeveP1Algorithm::KarhunenLoeveP1Algorithm(const Mesh & mesh,
    const CovarianceModel & covariance,
    const Scalar threshold)
  : KarhunenLoeveAlgorithmImplementation(covariance, threshold)
  , mesh_(mesh)
{
  // Nothing to do
}

/* Virtual constructor */
KarhunenLoeveP1Algorithm * KarhunenLoeveP1Algorithm::clone() const
{
  return new KarhunenLoeveP1Algorithm( *this );
}

/* Compute P1 gram as a vector of Eigen::Triplet<Scalar> */
TripletVector computeP1GramAsTriplets(const Mesh & mesh)
{
  // If no simplex, the P1 gram matrix is null
  if (mesh.getSimplicesNumber() == 0)
    return std::vector<Eigen::Triplet<Scalar>>();
  
  const UnsignedInteger simplexSize = mesh.getVertices().getDimension() + 1;
  
  SquareMatrix elementaryGram(simplexSize, Point(simplexSize * simplexSize, 1.0 / SpecFunc::Gamma(simplexSize + 2.0)));
  for (UnsignedInteger i = 0; i < simplexSize; ++i) 
    elementaryGram(i, i) *= 2.0;
  
  const UnsignedInteger simplicesSize = mesh.getSimplicesNumber();
  const Point simplexVolume(mesh.computeSimplicesVolume());
  
  typedef Eigen::Triplet<Scalar> Triplet;
  std::vector<Triplet> triplets;

  for (UnsignedInteger i = 0; i < simplicesSize; ++i)
  {
    const Indices simplex(mesh.getSimplex(i));
    const Scalar delta = simplexVolume[i];
        
    if (delta != 0)
      for (UnsignedInteger j = 0; j < simplexSize; ++j)
        for (UnsignedInteger k = 0; k < simplexSize; ++k)
          triplets.push_back(Triplet(simplex[j], simplex[k], delta * elementaryGram(j,k)));
  } // Loop over simplices
  
  return triplets;
}


/* Call to arpack routines to compute EV */
void computeEVWithArpack(const int nev,
                         const int ncv,
                         KLGenMatProd & op,
                         Point & eigenvalues,
                         Matrix & eigenvectors)
{    
  int const N = op.rows();
  int const ldv = N;
  int const ldz = N + 1;
  int const lworkl = 3 * (ncv * ncv) + 6 * ncv;

  Scalar const tol = 1e-3;    // TODO: add ResourceMap item?
  Scalar const sigmar(0.0);
  Scalar const sigmai(0.0);
  
  bool const rvec = true;

  std::vector<Scalar> resid(N);
  std::vector<Scalar> V(ncv * N);
  std::vector<Scalar> workd(3 * N);
  std::vector<Scalar> workl(lworkl);
  std::vector<Scalar> dr(nev + 1);
  std::vector<Scalar> di(nev + 1);
  std::vector<Scalar> z((N + 1) * (nev + 1));
  std::vector<Scalar> rwork(ncv);

  std::array<int, 11> iparam{};
  iparam[0] = 1;
  iparam[2] = 10 * N;
  iparam[3] = 1;
  iparam[4] = 0;  // number of ev found by arpack.
  iparam[6] = 1;

  std::array<int, 14> ipntr{};

  int info = 0, ido = 0;

  while (ido != 99) 
  {
    arpack::naupd(ido,
                  arpack::bmat::identity,
                  N,
                  arpack::which::largest_magnitude,
                  nev,
                  tol,
                  resid.data(),
                  ncv,
                  V.data(),
                  ldv,
                  iparam.data(),
                  ipntr.data(),
                  workd.data(),
                  workl.data(),
                  lworkl,
                  info);
    
    op.av(&(workd[ipntr[0] - 1]),
          &(workd[ipntr[1] - 1]));
  }
  
  // check number of ev found by arpack.
  if (    iparam[4] < nev     // arpack may succeed to compute more EV than expected
      ||  info != 0)
    throw InternalException(HERE) << "Error inside ARPACK routines: iparam[4]=" << iparam[4] << ", nev=" << nev << ", ncv=" << ncv << ", info=" << info;

  std::vector<Scalar> workev(3*ncv);
  std::vector<int> select(ncv);

  arpack::neupd(rvec,
                arpack::howmny::ritz_vectors,
                select.data(),
                dr.data(),
                di.data(),
                z.data(),
                ldz,
                sigmar,
                sigmai,
                workev.data(),
                arpack::bmat::identity,
                N,
                arpack::which::largest_magnitude,
                nev,
                tol,
                resid.data(),
                ncv,
                V.data(),
                ldv, 
                iparam.data(),
                ipntr.data(),
                workd.data(),
                workl.data(),
                lworkl,
                info);
  
  // Post-process eigenvalues and eigenvectors
  std::copy(dr.data(),dr.data()+nev,eigenvalues.begin());
  
  Collection<Scalar> eigenvectorsData((N+1)*(nev+1));
  
  for (int i=0; i<nev; ++i)
    std::copy(z.data()+i*(N+1), z.data()+i*(N+1)+N, eigenvectorsData.begin()+i*N);
  
  eigenvectors = Matrix(N, nev,eigenvectorsData);
}

/* Here we discretize the following Fredholm problem:
   \int_{\Omega}C(s,t)\phi_n(s)ds=\lambda_n\phi_n(t)
   using a P1 approximation of C and \phi_n:
   C(s,t)=\sum_{i,j}C(s_i,s_j)\theta_i(s)\theta_j(t)
   \phi_n(t)=\sum_k\alpha_k^n\theta_k(t)
   where s_i,s_j are vertices of the mesh, C(s_i,s_j)\in\cS^+_d(\R), \alpha_k^n\in\R^d
   leading to:
   \forall t\in\Omega, \sum_{i,j}C(s_i,s_j)\theta_j(t)\int_{\Omega}\theta_i(s)\sum_k\alpha_k^n\theta_k(s)=\lambda_n\sum_l\alpha_k^l\theta_l(t)
   For each value of n we get d equations in \alpha^n. We write these equations for t=s_1,...,s_N the N vertices of the mesh:
   \sum_{i,j}C(s_i,s_j)\theta_j(s_m)\int_{\Omega}\theta_i(s)\sum_k\alpha_k^n\theta_k(s)=\lambda_n\sum_l\alpha_l^n\theta_l(s_m)
   ie:
   \sum_i C(s_i,s_m)\int_{\Omega}\theta_i(s)\sum_k\alpha_k^n\theta_k(s)=\lambda_n\alpha_m^n\theta_m(s_m)
   In a block-matrix form we get:
   [C(s_1,s_1) ... C(s_1,s_N)][K_11 ... K_1N][\alpha_1]             [\alpha_1]
   [    ...            ...   ][ ...      ...][   ...  ] = \lambda_n [   ...  ]
   [C(s_N,s_1) ... C(s_N,s_N)][K_N1 ... K_NN][\alpha_N]             [\alpha_N]
   Where:
   K_ij = \int_{\Omega}\theta_i(s)\theta_j(s)ds I
   with I the dxd identity matrix
*/

void KarhunenLoeveP1Algorithm::run()
{
  runWithParameters(1.0, 1.0);
}

void KarhunenLoeveP1Algorithm::runWithParameters( const Scalar nevRatio,
                                                  const Scalar ncvRatio,
                                                  const Bool writeCsv)
{ 
  // Compute the gram of the mesh
  LOGINFO("Build the Gram matrix");
  TripletVector gram(computeP1GramAsTriplets(mesh_));
  
  const UnsignedInteger numVertices = mesh_.getVerticesNumber();    
  const Scalar epsilon = ResourceMap::GetAsScalar("KarhunenLoeveP1Algorithm-RegularizationFactor");
  
  if (epsilon > 0.0)
    for (UnsignedInteger i = 0; i < mesh_.getVerticesNumber(); ++i)
      gram.push_back(Triplet(i,i,epsilon));
  
  // Extend the Gram matrix of the mesh
  const UnsignedInteger dimension = covariance_.getOutputDimension();
  const UnsignedInteger augmentedDimension = dimension * numVertices;
  
  const UnsignedInteger nev = std::min(augmentedDimension*nevRatio + 1, augmentedDimension*1.0-2);
  const UnsignedInteger ncv = std::min(std::max(nev + 2, static_cast<UnsignedInteger>(nev*ncvRatio)),
                                       augmentedDimension);
  
  TripletVector tripletList;
  
  if (dimension == 1)
    tripletList = gram;
  else
  {
    for (UnsignedInteger p = 0; p < gram.size(); ++p)
    {
      Triplet g(gram[p]);
      for (UnsignedInteger k = 0; k < dimension; ++k)
        tripletList.push_back(Triplet(g.row()*dimension+k, g.col()*dimension+k, g.value()));
    }
  }
  
  SparseMatrix GSparse(augmentedDimension, augmentedDimension);
  GSparse.setFromTriplets(tripletList.begin(), tripletList.end());
  
  // Discretize the covariance model
  LOGINFO("Discretize the covariance model");
  CovarianceMatrix C(covariance_.discretize(mesh_.getVertices()));

  LOGINFO("Solve the eigenvalue problem");
  Matrix eigenVectors(augmentedDimension,nev);
  Point eigenValues(nev);
  
  KLGenMatProd op(C,GSparse);
  
  computeEVWithArpack(nev,
                      ncv,
                      op,
                      eigenValues,
                      eigenVectors);
        
  // Computing computed variance (i.e sum of computed eigenvalues)
  LOGINFO("Post-process the eigenvalue problem");  
  Scalar computedVariance = 0.0;
  for (UnsignedInteger k=0; k<eigenValues.getSize(); ++k)
    computedVariance += eigenValues[k];
      
  // Computing cumulated variance (i.e. sum of all eigenvalues)
  SquareMatrix G(augmentedDimension);
  for (UnsignedInteger k=0; k<tripletList.size(); ++k)
    G(tripletList[k].row(), tripletList[k].col()) += tripletList[k].value();

  Scalar cumulatedVariance = (C*G).computeTrace();
   
  LOGDEBUG(OSS(false) << "eigenValues=" << eigenValues);
  
  // Applying cut-off on spectrum
  LOGINFO("Extract the relevant eigenpairs");
  
  UnsignedInteger K = 0;
  
  Scalar selectedVariance = 0.0;
  while ((K < eigenValues.getSize()) && (eigenValues[K] >= threshold_ * cumulatedVariance)) 
  {
    selectedVariance += eigenValues[K];
    ++K;
  }
  
  LOGINFO(OSS() << "Selected " << K << " eigenvalues");
  
  // Reduce and rescale the eigenvectors
  MatrixImplementation projection(K, augmentedDimension);
  Point selectedEV(K);
  Collection<Function> modes(0);
  ProcessSample modesAsProcessSample(mesh_, 0, dimension);
  const UnsignedInteger meshDimension = mesh_.getDimension();
  SampleImplementation values(numVertices, dimension);
  Pointer<PiecewiseLinearEvaluation> evaluation1D;
  Pointer<P1LagrangeEvaluation> evaluationXD;
  if (meshDimension == 1)
    evaluation1D = new PiecewiseLinearEvaluation(mesh_.getVertices().getImplementation()->getData(), values);
  else
    evaluationXD = new P1LagrangeEvaluation(Field(mesh_, dimension));
  Point a(augmentedDimension);
  for (UnsignedInteger k = 0; k < K; ++k)
  {
    selectedEV[k] = eigenValues[k];
    for (UnsignedInteger i = 0; i < augmentedDimension; ++i)
      a[i] = eigenVectors(i, k);
    
    Eigen::VectorXd aSparse(augmentedDimension);
    std::copy(a.begin(), a.end(), aSparse.data());
    
    const Eigen::VectorXd GaSparse(GSparse * aSparse);
    Point Ga(augmentedDimension);
    std::copy(GaSparse.data(), GaSparse.data()+GaSparse.size(), Ga.begin());
    
    const Scalar norm = std::sqrt(a.dot(Ga));
    const Scalar factor = a[0] < 0.0 ? -1.0 / norm : 1.0 / norm;
    // Store the eigen modes in two forms
    values.setData(a * factor);
    modesAsProcessSample.add(values);
    if (meshDimension == 1)
    {
      evaluation1D->setValues(values);
      modes.add(*evaluation1D);
    }
    else
    {
      evaluationXD->setValues(values);
      modes.add(*evaluationXD);
    }

    // Build the relevant row of the projection matrix
    const Point b(Ga * (factor / sqrt(selectedEV[k])));
    for(UnsignedInteger i = 0; i < augmentedDimension; ++i)
      projection(k, i) = b[i];
  }
  
  if (writeCsv)
    std::cout << nev                                  << ","
              << ncv                                  << ","
              << K                                    << ","
              << cumulatedVariance                    << ","
              << computedVariance                     << ","
              << selectedVariance                     << ","
              << computedVariance / cumulatedVariance << ","
              << selectedVariance / cumulatedVariance << "," << std::flush;
     
  result_ = KarhunenLoeveResultImplementation(covariance_, threshold_, selectedEV, modes, modesAsProcessSample, projection);
}

/* Mesh accessor */
Mesh KarhunenLoeveP1Algorithm::getMesh() const
{
  return mesh_;
}

/* String converter */
String KarhunenLoeveP1Algorithm::__repr__() const
{
  OSS oss(true);
  oss << "class=" << KarhunenLoeveP1Algorithm::GetClassName()
      << ", mesh=" << mesh_;
  return oss;
}

/* String converter */
String KarhunenLoeveP1Algorithm::__str__(const String & ) const
{
  OSS oss(false);
  oss << "class=" << KarhunenLoeveP1Algorithm::GetClassName()
      << ", mesh=" << mesh_;
  return oss;
}

/* Method save() stores the object through the StorageManager */
void KarhunenLoeveP1Algorithm::save(Advocate & adv) const
{
  KarhunenLoeveAlgorithmImplementation::save(adv);
  adv.saveAttribute( "mesh_", mesh_ );
}

/* Method load() reloads the object from the StorageManager */
void KarhunenLoeveP1Algorithm::load(Advocate & adv)
{
  KarhunenLoeveAlgorithmImplementation::load(adv);
  adv.loadAttribute( "mesh_", mesh_ );
}

END_NAMESPACE_OPENTURNS
