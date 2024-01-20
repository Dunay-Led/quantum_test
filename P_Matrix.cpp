#ifndef P_Matrix
#define P_Matrix
#include "header.hpp"
using namespace Eigen;

extern int nm;
extern int Number_e;


Matrix<std::complex<double>, Dynamic, Dynamic> P_Matrix_form_C_count(Matrix<std::complex<double>, Dynamic, Dynamic> MatrixC) {

    
    Matrix<std::complex<double>, Dynamic, Dynamic> MatrixP;
    MatrixP.resize(nm, nm);




    for (int j = 0; j < nm; ++j) {
        for (int m = 0; m < nm; ++m) {
            std::complex<double> P(0, 0);

            for (int a = 0; a < Number_e / 2; ++a) {

                P +=   MatrixC(j, a) * conj(MatrixC(m, a));

            }
            MatrixP(j, m) = P; //

        }
    }
   

   
    return MatrixP;

    

}

#endif