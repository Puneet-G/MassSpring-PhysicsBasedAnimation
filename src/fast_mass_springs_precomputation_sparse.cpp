#include "fast_mass_springs_precomputation_sparse.h"
#include "signed_incidence_matrix_sparse.h"
#include <vector>

bool fast_mass_springs_precomputation_sparse(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & E,
  const double k,
  const Eigen::VectorXd & m,
  const Eigen::VectorXi & b,
  const double delta_t,
  Eigen::VectorXd & r,
  Eigen::SparseMatrix<double>  & M,
  Eigen::SparseMatrix<double>  & A,
  Eigen::SparseMatrix<double>  & C,
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > & prefactorization)
{
  /////////////////////////////////////////////////////////////////////////////
  // Replace with your code
  std::vector<Eigen::Triplet<double> > ijv;
  const int n = V.rows();
  for(int i = 0;i<n;i++) ijv.emplace_back(i,i,1);
  Eigen::SparseMatrix<double> Q(n,n);
  Q.setFromTriplets(ijv.begin(),ijv.end());
  //Compute A
  signed_incidence_matrix_sparse(V.rows(), E, A);

  //Compute r (the edge lengths), the difference betweent the vertex positions of an edge
  r.resize(E.rows());
  //r = (A * V);
  for (int i = 0; i < E.rows(); i++) {
	  r(i) = (V.row(E(i, 0)) - V.row(E(i, 1))).norm();
  }

  //Build the M matrix, square Matrix with the masses across the diagonal
  M.resize(V.rows(), V.rows());
  std::vector<Eigen::Triplet<double> > ijvMass;
  for (int i = 0; i < V.rows(); i++) ijvMass.emplace_back(i, i, m(i));
  M.setFromTriplets(ijvMass.begin(), ijvMass.end());

  //Build the C matrix, of selected pinned masses
  std::vector<Eigen::Triplet<double> > ijvPinned;
  C.resize(b.rows(), V.rows());
  for (int j = 0; j < b.rows(); j++) {
	  ijvPinned.emplace_back(j, b(j), 1.0);
  }
  C.setFromTriplets(ijvPinned.begin(), ijvPinned.end());

  //Build Q Matrix
  Q = (k * A.transpose() * A) + ((1.0 / pow(delta_t, 2)) * M) + (1e10 * C.transpose() * C);
  /////////////////////////////////////////////////////////////////////////////
  prefactorization.compute(Q);
  return prefactorization.info() != Eigen::NumericalIssue;
}
