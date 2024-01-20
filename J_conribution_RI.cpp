#ifndef J_contr
#define J_contr
#include <libint2.hpp>
#include "header.hpp"
#define Mat_init Matrix<std::complex<double>, Dynamic, Dynamic>
#define Mat_init_r Matrix<double, Dynamic, Dynamic>
using namespace Eigen;
using namespace std;
using libint2::Shell;
using libint2::BasisSet;
using libint2::Operator;
using namespace libint2;


extern int nm_aux;
extern int nm;
extern string file_path;
extern string basis;
extern int Number_e;

complex<double> J_contribution_calc(Mat_init MatrixP, Mat_init MatrixC, Mat_init V_R_1, Mat_init_r integrals, int alpha, int beta, vector<int> cumulant, vector<int> compressed_num, vector<double> compressed_data, Mat_init i_a_B){

int n_occ = Number_e/2;

Mat_init D;
Mat_init DV_1;
Mat_init Sqrt_Lambda;
Mat_init V_R_1_2; 
Mat_init M;
Mat_init M_first;

M.resize(n_occ*nm_aux, nm);
M_first.resize(n_occ*nm, nm_aux);
D.resize(1, nm_aux);
Sqrt_Lambda.resize(nm_aux, nm_aux);
V_R_1_2.resize(nm_aux, nm_aux);
DV_1.resize(nm_aux, 1);

D.setZero();
DV_1.setZero();
Sqrt_Lambda.setZero();
V_R_1_2.setZero();
M.setZero();
M_first.setZero();

complex<double> a_J_b (0, 0);


// Taking Sqrt(V^-1)


SelfAdjointEigenSolver<Eigen::Matrix<complex<double>, Dynamic, Dynamic>> eigensolver(V_R_1);
    for (int i = 0; i!=nm_aux; ++i){
        Sqrt_Lambda(i, i) = sqrt(eigensolver.eigenvalues()(i, 0));
    }
    V_R_1_2 = eigensolver.eigenvectors() * Sqrt_Lambda * eigensolver.eigenvectors().transpose();




for (int f_aux = 0; f_aux != nm*nm; ++f_aux){
    for (int first = cumulant[f_aux]; first != cumulant[f_aux+1]; ++first){
        
        
        
        D(0, compressed_num[first]) += MatrixP(f_aux/nm, f_aux%nm) * compressed_data[first];
        
    }
}

DV_1 = V_R_1 * D.transpose();


int x_num = alpha*nm+beta;

for (int first = cumulant[x_num]; first != cumulant[x_num+1]; ++first){

    a_J_b += compressed_data[first] * DV_1(compressed_num[first], 0); // Speed can be increased

}


M_first = i_a_B * V_R_1_2;


complex<double> a_K_b (0, 0); 
// << alpha << beta << endl;
for (int i = 0; i != n_occ; ++i){
    for (int A = 0; A != nm_aux; ++A){

        //// << a_K_b << endl;
        a_K_b += M_first(i * nm + alpha, A) * M_first(i * nm + beta, A);
        
        //

    }
}



return 2.0 * a_J_b - a_K_b ;

}



#endif 