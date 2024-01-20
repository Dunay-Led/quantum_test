#ifndef SCF_c
#define SCF_c
#include "header.hpp"
using namespace Eigen;
using namespace std;

extern int nm;

int SCF_procedure(Mat_init MatrixC, Mat_init MatrixHCore, Mat_init MatrixX, Mat_init_r integrals, Mat_init V_RI_1, vector<int> cumulant, vector<int> compressed_num, vector<double> compressed_data){

libint2::initialize();

    Matrix<std::complex<double>, Dynamic, Dynamic> MatrixP;
    Matrix<complex<double>, Dynamic, Dynamic> MatrixF;
    Matrix<complex<double>, Dynamic, Dynamic> _MatrixF_;
    
    Matrix<double, Dynamic, 1> e;
    double deltaE_tot;
    MatrixP.resize(nm, nm);
    MatrixF.resize(nm, nm);
    _MatrixF_.resize(nm, nm);
    e.resize(nm, 1);
    complex<double> Efirst (-10, -10);
    int number_iter = 1;

    


do {
    
        MatrixP = P_Matrix_form_C_count(MatrixC);

        MatrixF = Fock_Matrix_calc(MatrixHCore, MatrixC, integrals, V_RI_1, cumulant, compressed_num, compressed_data);



        _MatrixF_ = MatrixX*MatrixF*MatrixX;
        SelfAdjointEigenSolver<Eigen::Matrix<complex<double>, Dynamic, Dynamic>> eigensolver(_MatrixF_);
        MatrixC = MatrixX*eigensolver.eigenvectors();
        e = 27.2113961318065*eigensolver.eigenvalues();
        //e = eigensolver.eigenvalues();


    complex<double> E0;
    for (int d=0; d < nm; ++d){
        for (int y = 0; y < nm; ++y){
            E0 += MatrixP(y, d)*(MatrixHCore(d, y)+MatrixF(d, y));
        }
    }

    deltaE_tot = abs(Efirst-E0);

    cout << "dE\n" << deltaE_tot << endl;

    Efirst = E0;

    ++number_iter;

    if (number_iter > 31){

        cout <<  "Iterations haven't converged" << endl;

        break;

    }

    E0 = 0;
    for (int d=0; d < nm; ++d){
        for (int y = 0; y < nm; ++y){
            E0 += MatrixP(y, d)*(MatrixHCore(d, y)+MatrixF(d, y));
        }
    }
    //cout << "F\n" << MatrixF << endl;
    //cout << "F\n" << MatrixF << endl;
    getchar();
    printf("Total Energy is      %9.9f eV\n", E0.real());
    cout << "number iter\n" << number_iter << endl;

    }   while (1 == 1);


    // Total Energy count

    complex<double> E0;
    for (int d=0; d < nm; ++d){
        for (int y = 0; y < nm; ++y){
            E0 += MatrixP(y, d)*(MatrixHCore(d, y)+MatrixF(d, y));
        }
    }

    printf("Total Energy is      %9.9f eV\n", E0.real());
    std::cout <<  endl << e << endl << E0 << endl << "Iterations number:  " << number_iter << endl ;
    cout << MatrixC;
    //std::cout <<  "F : " << endl<< MatrixF << endl << "F MO : " << MatrixC.inverse()*MatrixF*MatrixC << endl;

    libint2::finalize();
    return 1;
}
#endif