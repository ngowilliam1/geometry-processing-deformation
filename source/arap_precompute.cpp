#include "arap_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
  // REPLACE WITH YOUR CODE

  // Calculating data is similar to biharmonic_precompute.cpp 
  // L is a #V by #V cotangent matrix, each row i corresponding to V(i,:)
  Eigen::SparseMatrix<double> L(V.rows(),V.rows());
  igl::cotmatrix(V, F, L);
  Eigen::SparseMatrix<double> Aeq;
  igl::min_quad_with_fixed_precompute(L, b, Aeq, false, data);

  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  // 54 = 3*3*3*2
  tripletList.reserve(F.rows() * 54);

  // K represents the bilinear form that combines unknown vertex positions and unknown rotations
  for (int f = 0; f< F.rows(); f++){
    for (int edge =0; edge<3; edge++){
      int i = F(f, edge % 3);
      int j = F(f, (edge+1) % 3);
      Eigen::Vector3d vij = V.row(i) - V.row(j);
      Eigen::Vector3d eij = L.coeff(i, j) * vij;
      
      for (int ki = 0; ki<3; ki++){ 
        int k = F(f, ki);
        for (int beta = 0; beta < 3; beta++){
          tripletList.push_back(T(i, 3*k+beta, eij[beta]));
          tripletList.push_back(T(j, 3*k+beta, -eij[beta]));
        }
      }
    }
  }
  K.resize(V.rows(), 3*V.rows());
  K.setFromTriplets(tripletList.begin(), tripletList.end());
  K = K / 6.0;

}
