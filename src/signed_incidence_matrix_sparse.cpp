#include "signed_incidence_matrix_sparse.h"
#include <vector>

void signed_incidence_matrix_sparse(
	const int n,
	const Eigen::MatrixXi& E,
	Eigen::SparseMatrix<double>& A)
{
	//////////////////////////////////////////////////////////////////////////////
	// Replace with your code
	std::vector<Eigen::Triplet<double> > ijv;
	for (int i = 0; i < E.rows(); i++) {
		ijv.emplace_back(i, E(i, 0), 1.0);
		ijv.emplace_back(i, E(i, 1), -1.0);
	}
	//ijv.emplace_back(0, 0, 1234.5678);
	A.resize(E.rows(), n);
	A.setFromTriplets(ijv.begin(), ijv.end());
	//////////////////////////////////////////////////////////////////////////////
}
