#ifndef L_Diag
#define L_Diag
#include "header.h"
using namespace Eigen;


auto Diagonal_Matrix_lambda(Matrix<std::complex<double>, Dynamic, Dynamic> MatrixS, int nm) {

Eigen::SelfAdjointEigenSolver < Eigen::Matrix < std::complex<double>  , Dynamic, Dynamic >> orthogonal(MatrixS);
Eigen::Matrix<std::complex<double>, Dynamic, Dynamic> MatrixU;
Eigen::Matrix<std::complex<double>, Dynamic, Dynamic> MatrixS_1_2 ;
Eigen::Matrix<std::complex<double>, Dynamic, Dynamic> MatrixS1_2;

MatrixS_1_2.resize(nm, nm);
MatrixS1_2.resize(nm, nm);

for (int j = 0; j!=nm; ++j){
    for (int i = 0; i!=nm; ++i){
        MatrixS_1_2(i,j) = std::complex<double> (0,0);
        MatrixS1_2(i,j) = std::complex<double> (0,0);

    }
}
MatrixU = orthogonal.eigenvectors();



auto lambdas = orthogonal.eigenvalues();

//std::cout << lambdas << std::endl;
    //getchar();


for (int i = 0; i < nm; i++) {

    MatrixS_1_2 += MatrixU.col(i) * MatrixU.col(i).transpose() / sqrt(lambdas(i));
    MatrixS1_2 += MatrixU.col(i) * MatrixU.col(i).transpose() * lambdas(i);

}
    //std::cout << MatrixS1_2 << std::endl << std::endl<< MatrixS_1_2;
    //getchar();
std::array<Eigen::Matrix<std::complex<double>, Dynamic, Dynamic>, 2> ans = {MatrixS_1_2, MatrixS1_2};

return  ans;
}

#endif