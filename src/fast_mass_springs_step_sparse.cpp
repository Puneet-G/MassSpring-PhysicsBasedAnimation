#include "fast_mass_springs_step_sparse.h"
#include <igl/matlab_format.h>

void fast_mass_springs_step_sparse(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & E,
  const double k,
  const Eigen::VectorXi & b,
  const double delta_t,
  const Eigen::MatrixXd & fext,
  const Eigen::VectorXd & r,
  const Eigen::SparseMatrix<double>  & M,
  const Eigen::SparseMatrix<double>  & A,
  const Eigen::SparseMatrix<double>  & C,
  const Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > & prefactorization,
  const Eigen::MatrixXd & Uprev,
  const Eigen::MatrixXd & Ucur,
  Eigen::MatrixXd & Unext)
{
  //////////////////////////////////////////////////////////////////////////////
  // Replace with your code
	//Build Matrix of Pinned rest Positions
	//std::cout << "Input\n" << "A\n" << A << "\n\nM\n" << M << "\n\nV\n" << V << "\n\nC\n" << C;
	//Eigen::MatrixXd pinV(b.rows(), 3);
	//for (int i = 0; i < b.rows(); i++) {
	//	pinV.row(i) = V.row(b(i));
	//}
	Eigen::MatrixXd l = Ucur;
	for (int iter = 0; iter < 50; iter++)
	{	
		//Step 1 Find D
		Eigen::MatrixXd d(E.rows(), 3);
		for (int i = 0;  i < E.rows(); i++) {
			d.row(i) = (l.row(E(i, 0)) - l.row(E(i, 1))).normalized() * r(i);
		}
		//std::cout << "\n\nl - Pt Iterated\n" << l << "\n\nPt-step\n" << Uprev << "\n\nd\n" << d;
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
