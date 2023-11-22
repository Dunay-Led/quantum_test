//#include "header.h"
#ifndef Schmidt_orth
#define Schmidt_orth
#include <iostream>
#include "Eigen/Dense"
#include <complex>
#include <math.h>
#define Mat_init Matrix<std::complex<double>, Dynamic, Dynamic>

using namespace Eigen;
using namespace std;

// Takes Matrix as a first arg and then orthogonalising it

auto Schmidt_orthonormal(Mat_init MatrixX){

    int Matrix_size = MatrixX.col(1).size();

    //GRAMM-SCHMIDT CYCLE
    for (int i = 0; i != Matrix_size; ++i){
        Matrix<complex<double>, Dynamic, 1> extra_v;
        extra_v.resize(Matrix_size, 1);

        extra_v.setZero();


        for (int j = 0; j!=i; ++j){
            extra_v += (MatrixX.col(i).transpose()*MatrixX.col(j))*MatrixX.col(j)/(MatrixX.col(j).transpose()*MatrixX.col(j)); // Main Schmidt formula


        }
        MatrixX.col(i) = (MatrixX.col(i)-extra_v);


    }


    for (int i = 0; i != Matrix_size ; ++i){

            complex<double> len (sqrt((MatrixX.col(i).transpose()*MatrixX.col(i))(0, 0).real()), 0);
            MatrixX.col(i) = MatrixX.col(i)/len;

    }

    return MatrixX; // return orthogonal matrix
}

#endif