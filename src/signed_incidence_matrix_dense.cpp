#include "signed_incidence_matrix_dense.h"

// Construct the sparse incidence matrix for a given edge network.
//
// Inputs: 
//   n  number of vertices (#V)
//   E  #E by 2 list of edge indices into rows of V
// Outputs:
//   A  #E by n signed incidence matrix

void signed_incidence_matrix_dense(
  const int n,
  const Eigen::MatrixXi & E,
  Eigen::MatrixXd & A)
{
  //////////////////////////////////////////////////////////////////////////////
  // Replace with your code
  A = Eigen::MatrixXd::Zero(E.rows(),n);
  for (int i = 0; i < E.rows(); i++) {
	  A(i, E(i, 0)) = 1.0;
	  A(i, E(i, 1)) = -1.0;
  }
  //////////////////////////////////////////////////////////////////////////////
}
