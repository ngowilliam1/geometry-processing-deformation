#include "../include/arap_single_iteration.h"
#include <igl/polar_svd3x3.h>
#include <igl/min_quad_with_fixed.h>

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
  // REPLACE WITH YOUR CODE
  Eigen::MatrixXd CTranspose = (U.transpose() * K).transpose();
  Eigen::MatrixXd R(3*data.n, 3);

  for (int i = 0; i < data.n; i++){
    Eigen::Matrix3d M_k, R_k;
	  M_k = CTranspose.block(3*i, 0, 3, 3);
    igl::polar_svd3x3(M_k, R_k);
    R.block(3*i, 0, 3, 3) = R_k;
  }

  Eigen::MatrixXd Beq;
  igl::min_quad_with_fixed_solve(data, K*R, bc, Beq, U);
}
