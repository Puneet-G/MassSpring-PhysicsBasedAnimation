#include "fast_mass_springs_precomputation_dense.h"
#include "signed_incidence_matrix_dense.h"
#include <Eigen/Dense>

// Precompute matrices and factorizations necessary for the "Fast Simulation of
// Mass-Spring Systems" method.
//
// Inputs: 
//   V  #V by 3 list of vertex positions
//   E  #E by 2 list of edge indices into rows of V
//   k  spring stiffness
//   m  #V list of masses 
//   b  #b list of indices of fixed vertices as indices into rows of V
//   delta_t  time step in seconds
// Outputs:
//   r  #E list of edge lengths
//   M  #V by #V mass matrix
//   A  #E by #V signed incidence matrix
//   C  #b by #V selection matrix
//   prefactorization  LLT prefactorization of energy's quadratic matrix

bool fast_mass_springs_precomputation_dense(
	const Eigen::MatrixXd& V,
	const Eigen::MatrixXi& E,
	const double k,
	const Eigen::VectorXd& m,
	const Eigen::VectorXi& b,
	const double delta_t,
	Eigen::VectorXd& r,
	Eigen::MatrixXd& M,
	Eigen::MatrixXd& A,
	Eigen::MatrixXd& C,
	Eigen::LLT<Eigen::MatrixXd>& prefactorization)
{
	/////////////////////////////////////////////////////////////////////////////
	// Replace with your code
	//Compute A
	signed_incidence_matrix_dense(V.rows(), E, A);

	Eigen::MatrixXd Q = Eigen::MatrixXd::Identity(V.rows(), V.rows());

	//Compute r (the edge lengths), the difference betweent the vertex positions of an edge
	r.resize(E.rows());
	//r = (A * V);
	for (int i = 0; i < E.rows(); i++) {
		r(i) = (V.row(E(i, 0)) - V.row(E(i, 1))).norm();
	}

	//Build the M matrix, square Matrix with the masses across the diagonal
	M.resize(V.rows(), V.rows());
	M = m.asDiagonal();

	//Build the C matrix, of selected pinned masses
	C = Eigen::MatrixXd::Zero(b.rows(), V.rows());
	for (int j = 0; j < b.rows(); j++) {
		C(j, b(j)) = 1.0;
	}

	//Build Q Matrix
	Q = (k * A.transpose() * A) + ((1.0 / pow(delta_t, 2)) * M) + (1e10 * C.transpose()*C);
	/////////////////////////////////////////////////////////////////////////////
	prefactorization.compute(Q);
	return prefactorization.info() != Eigen::NumericalIssue;
}
