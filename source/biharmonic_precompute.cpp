#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>

void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{
  // REPLACE WITH YOUR CODE

  // L is a #V by #V cotangent matrix, each row i corresponding to V(i,:)
  Eigen::SparseMatrix<double> L(V.rows(),V.rows());
  igl::cotmatrix(V, F, L);

  // We now need M inverse, #V by #V mass matrix
  Eigen::SparseMatrix<double> M, MInverse;
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
  igl::invert_diag(M, MInverse);

  // Q is L.transpose * M.inverse * L
  Eigen::SparseMatrix<double> Q = L.transpose() * MInverse * L;

  // Aeq is m by n list of linear equality constraint coefficients, but we don't use it
  Eigen::SparseMatrix<double> Aeq;
  igl::min_quad_with_fixed_precompute(Q, b, Aeq, false, data);

}
