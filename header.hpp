
#ifndef HARTRY_FOCKC_HEADER_H
#define HARTRY_FOCKC_HEADER_H

#define Mat_init Matrix<std::complex<double>, Dynamic, Dynamic>
#define Mat_init_r Matrix<double, Dynamic, Dynamic>

#define RHF_Mem_RI_Compressed
#include <string>
#include <tuple>
#include <map>
#include <chrono>
#include "Eigen/Dense"
#include <iostream>
#include <vector>
#include <fstream>
#include <libint2.hpp>
#include <typeinfo>
#include <cmath>
#include <complex>
using namespace libint2;
using namespace Eigen;
using namespace std;
using namespace std::chrono;
using libint2::Shell;
using libint2::BasisSet;
using libint2::Operator;


Mat_init P_Matrix_form_C_count(Matrix<std::complex<double>, Dynamic, Dynamic> MatrixC);
std::array<Eigen::Matrix<std::complex<double>, Dynamic, Dynamic>, 2> Diagonal_Matrix_lambda(Matrix<std::complex<double>, Dynamic, Dynamic> MatrixS, int nm);
Mat_init Extended_Huckel_Guess();
complex<double> J_contribution_calc(Mat_init MatrixP, Mat_init MatrixC, Mat_init V_R_1, Mat_init_r integrals, int alpha, int beta, vector<int> cumulant, vector<int> compressed_num, vector<double> compressed_data, Mat_init i_a_B);
Mat_init Fock_Matrix_calc(Mat_init MatrixHCore, Mat_init MatrixC, Mat_init_r integrals, Mat_init V_RI_1, vector<int> cumulant, vector<int> compressed_num, vector<double> compressed_data);
Mat_init Schmidt_orthonormal(Mat_init MatrixX);
int SCF_procedure(Mat_init MatrixC, Mat_init MatrixHCore, Mat_init MatrixX, Mat_init_r integrals, Mat_init V_RI_1, vector<int> cumulant, vector<int> compressed_num, vector<double> compressed_data);
int SOSCF_calc(Mat_init MatrixC, Mat_init MatrixHCore, Mat_init MatrixX1, Mat_init_r integrals, Mat_init V_RI_1, vector<int> cumulant, vector<int> compressed_num, vector<double> compressed_data, int Number_e);


#endif //HARTRY_FOCKC_HEADER_H
