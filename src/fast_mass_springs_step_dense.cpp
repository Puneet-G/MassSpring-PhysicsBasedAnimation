#include "fast_mass_springs_step_dense.h"
#include <igl/matlab_format.h>
// Conduct a single step of the "Fast Simulation of Mass-Spring Systems" method.
//
// Inputs: 
//   V  #V by 3 list of **rest** vertex positions
//   E  #E by 2 list of edge indices into rows of V
//   k  spring stiffness
//   b  #b list of indices of fixed vertices as indices into rows of V
//   delta_t  time step in seconds
//   fext  #V by 3 list of external forces
//   r  #E list of edge lengths
//   M  #V by #V mass matrix
//   A  #E by #V signed incidence matrix
//   C  #b by #V selection matrix
//   prefactorization  LLT prefactorization of energy's quadratic matrix
//   Uprev  #V by 3 list of previous vertex positions (at time t-delta_t)
//   Ucur  #V by 3 list of current vertex positions (at time t)
// Outputs:
//   Unext #V by 3 list of previous vertex positions (at time t+delta_t)
void fast_mass_springs_step_dense(
	const Eigen::MatrixXd& V,
	const Eigen::MatrixXi& E,
	const double k,
	const Eigen::VectorXi& b,
	const double delta_t,
	const Eigen::MatrixXd& fext,
	const Eigen::VectorXd& r,
	const Eigen::MatrixXd& M,
	const Eigen::MatrixXd& A,
	const Eigen::MatrixXd& C,
	const Eigen::LLT<Eigen::MatrixXd>& prefactorization,
	const Eigen::MatrixXd& Uprev,
	const Eigen::MatrixXd& Ucur,
	Eigen::MatrixXd& Unext)
{
	//////////////////////////////////////////////////////////////////////////////
	// Replace with your code
	//Build Matrix of Pinned rest Positions
	Eigen::MatrixXd l = Ucur;
	for (int iter = 0; iter < 50; iter++)
	{	
		//Step 1 Find D
		Eigen::MatrixXd d(E.rows(), 3);
		for (int i = 0;  i < E.rows(); i++) {
			d.row(i) = (l.row(E(i, 0)) - l.row(E(i, 1))).normalized() * r(i);
		}
		//Step 2:
		// Build y
		Eigen::MatrixXd y = Eigen::MatrixXd::Zero(V.rows(), 3);
		y = ((1.0 / pow(delta_t, 2)) * M) * (2.0 * Ucur - Uprev) + fext + 1e10 * C.transpose() * C * V;
		
		//Build B
		Eigen::MatrixXd B = Eigen::MatrixXd::Zero(V.rows(), 3);
		B = k * A.transpose() * d + y;
		
		Unext = prefactorization.solve(B);
		l = Unext;
	}
	//////////////////////////////////////////////////////////////////////////////
}
