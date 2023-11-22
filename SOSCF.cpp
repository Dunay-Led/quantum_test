#ifndef SOSCF
#define SOSCF
#include <libint2.hpp>
#include "header.h"
#include <math.h>
#define Mat_init Matrix<std::complex<double>, Dynamic, Dynamic>
using namespace Eigen;
using namespace std;
using libint2::Shell;
using libint2::BasisSet;
using libint2::Operator;
using namespace libint2;
using namespace std::chrono;

extern int nm;

int SOSCF_calc(Mat_init MatrixC, Mat_init MatrixHCore, Mat_init MatrixX1, double *integrals, Mat_init V_RI_1){

    libint2::initialize();

    // Initializing data
     

    Matrix<complex<double>, Dynamic, Dynamic> MatrixP;
    Matrix<complex<double>, Dynamic, Dynamic> MatrixA;
    Matrix<complex<double>, Dynamic, Dynamic> MatrixG;
    Matrix<complex<double>, Dynamic, Dynamic> MatrixI;
    Matrix<complex<double>, Dynamic, Dynamic> MatrixI_exp;
    Matrix<complex<double>, Dynamic, Dynamic> MatrixC_New;
    Matrix<complex<double>, Dynamic, Dynamic> MatrixF;
    Matrix<complex<double>, Dynamic, Dynamic> MatrixF_MO;
    Matrix<complex<double>, Dynamic, Dynamic> Matrix_F_;
    Matrix<complex<double>, Dynamic, Dynamic> MatrixF_MO_New;
    Matrix<complex<double>, Dynamic, Dynamic> _MatrixF_;
    Matrix<complex<double>, Dynamic, Dynamic> MatrixGradient;
    Matrix<complex<double>, Dynamic, Dynamic> MatrixGradient_New;
    Matrix<complex<double>, Dynamic, Dynamic> MatrixF_New;
    Matrix<complex<double>, Dynamic, Dynamic> MatrixX;
    Matrix<complex<double>, Dynamic, Dynamic> MatrixXT;
    Matrix<complex<double>, Dynamic, Dynamic> MatrixB;
    Matrix<complex<double>, Dynamic, Dynamic> vector_s;
    Matrix<complex<double>, Dynamic, Dynamic> vector_y;
    Matrix<complex<double>, Dynamic, Dynamic> MatrixU;
    Matrix<complex<double>, Dynamic, Dynamic> a;
    Matrix<complex<double>, Dynamic, Dynamic> q;
    Matrix<complex<double>, Dynamic, Dynamic> z;
    Matrix<double, Dynamic, 1> e;
    e.resize(nm,1);
    MatrixP.resize(nm, nm);
    MatrixC_New.resize(nm, nm);
    MatrixG.resize(nm, nm);
    MatrixF.resize(nm, nm);
    MatrixF_MO.resize(nm, nm);
    MatrixU.resize(nm, nm);
    _MatrixF_.resize(nm, nm);
    MatrixF_New.resize(nm, nm);
    MatrixF_MO_New.resize(nm, nm);
    MatrixI_exp.resize(nm, nm);
    MatrixA.resize(nm, nm);
    
    complex<double> Efirst (-10, -10);
    double Occupied_orb_num = Number_e/2;
    double Virtuall_orb_num = nm-Occupied_orb_num;
    
    complex<double> k;
    complex<double> four (4,0);
    complex<double> one (1,0);
    complex<double> deltaE0;
    complex<double> E0 = 0;
    int number_iter = 0;
    int step_in_mem = 10;        // Saved steps in l-bfgs
    double conv_crit = 0.00000001;  // Convergence criterion
    double deltaE_tot;


    int NO_NV = Occupied_orb_num*Virtuall_orb_num;      

    MatrixB.resize(NO_NV, NO_NV);
    MatrixB.setZero();  
    MatrixI.resize(NO_NV, NO_NV); 
    MatrixGradient.resize(NO_NV, 1); 
    vector_y.resize(NO_NV, step_in_mem);
    vector_s.resize(NO_NV, step_in_mem);
    q.resize(NO_NV, 1);
    z.resize(NO_NV, 1);
    a.resize(1, step_in_mem);
    
    MatrixGradient_New.resize(NO_NV, 1);       
    vector_s.setZero();
    vector_y.setZero();
    MatrixI.setIdentity();
    MatrixI_exp.setIdentity();

    // Matrix of step is saved as folows: Row-major form, ( (o1v1, o1v2 ... o1v(Virtuall_orb_num) ... (o(Oc_num)v1, ... o(Occ_num)v(Virtuall_orb_num)) )

    MatrixF = Fock_Matrix_calc(MatrixHCore, MatrixC, integrals, V_RI_1);
    MatrixF_MO = MatrixC.transpose()*MatrixF*MatrixC;
   

    // MatrixB initiall guess
   
    for (int i = 0; i!= Occupied_orb_num; ++i){
                for (int a = 0; a!= Virtuall_orb_num; ++a){
                    MatrixB(i*Virtuall_orb_num+a, i*Virtuall_orb_num+a) = four*MatrixF_MO(Occupied_orb_num+a, Occupied_orb_num+a) - four*MatrixF_MO(i, i);
                }
            }

 

    // Calculating Gradient Matrix (Row-major form)

    for (int i = 0; i!= Occupied_orb_num; ++i){
            for (int a = 0; a!= Virtuall_orb_num; ++a){
                MatrixGradient(i*Virtuall_orb_num+a) = four * MatrixF_MO(i, a+Occupied_orb_num);
            }
        }
        
        
    z = MatrixB.inverse()*MatrixGradient;


    cout << "Hello" << endl;
    do {

    auto start = steady_clock::now(); // Timer start


    
    double min_f = 1000;

    


    MatrixF_MO = MatrixC.transpose()*MatrixF*MatrixC;       //Fock Matrix in MO basis


    MatrixXT.resize(NO_NV, 1);

    MatrixXT = - z;


    
    double help  = 0;

    MatrixXT.resize(Virtuall_orb_num, Occupied_orb_num);


    // Stating, that step cannot be more than 0.1
    
    double mins = 0.1;  // Max step = 0.1

    for (int i = 0; i != Virtuall_orb_num; ++i){
        for (int j = 0; j != Occupied_orb_num; ++j){
            
            if (abs(MatrixXT(i, j).real()) > help){
                help = abs(MatrixXT(i, j).real());
            }

        }
    }

    double help_m = mins/help;
    complex<double> help_1 (help_m, 0);


    if (help > mins){
        MatrixXT = help_1*MatrixXT;
    }






    MatrixX = MatrixXT.transpose();


    // Matrix A Generation

    MatrixA.setZero();

 
    MatrixA.block(0, Occupied_orb_num, Occupied_orb_num, Virtuall_orb_num) = MatrixX;
    MatrixA.block(Occupied_orb_num, 0, Virtuall_orb_num, Occupied_orb_num) = -MatrixX.transpose();
   

    MatrixU = Schmidt_orthonormal(MatrixI_exp + MatrixA);

    


    MatrixP = P_Matrix_form_C_count(MatrixC);







    // Total Energy count

    deltaE0 = -E0;

    E0 = 0;

    for (int d=0; d < nm; ++d){
        for (int y = 0; y < nm; ++y){
            E0 += MatrixP(y, d)*(MatrixHCore(d, y)+MatrixF(d, y));
        }
    }

    deltaE0 += E0;





    // New Matricies obtaining
    

    MatrixC_New = MatrixC*MatrixU.transpose();      // New Coefficient matrix creation


    MatrixF_New = Fock_Matrix_calc(MatrixHCore, MatrixC_New, integrals, V_RI_1);

    MatrixF_MO_New = MatrixC_New.transpose()*MatrixF_New*MatrixC_New;

     for (int i = 0; i!= Occupied_orb_num; ++i){
            for (int a = 0; a!= Virtuall_orb_num; ++a){
                MatrixGradient_New(i*Virtuall_orb_num+a) = four * MatrixF_MO_New(i, a+Occupied_orb_num);
            }
        }
    

    for (int i = 0; i!= Occupied_orb_num; ++i){
                for (int a = 0; a!= Virtuall_orb_num; ++a){
                    MatrixB(i*Virtuall_orb_num+a, i*Virtuall_orb_num+a) = four*MatrixF_MO_New(Occupied_orb_num+a, Occupied_orb_num+a) - four*MatrixF_MO_New(i, i);
                }
            }






    MatrixXT.resize(NO_NV, 1);


    // Updating saved l-bfgs information
   

    vector_s.block(0, 0, NO_NV, step_in_mem-1) = vector_s.block(0, 1, NO_NV, step_in_mem-1);
    vector_s.col(step_in_mem-1) = MatrixXT;
    vector_y.block(0, 0, NO_NV, step_in_mem-1) = vector_y.block(0, 1, NO_NV, step_in_mem-1);
    vector_y.col(step_in_mem-1) = MatrixGradient_New - MatrixGradient;


    // Updating Hessian


    q = MatrixGradient_New;


    int m;
    if (number_iter > step_in_mem-1){
        m = step_in_mem;
    }

    else {
        m = number_iter+1;
    }


    for (int j = 0; j!= m; ++j){
        
        auto p = one/(vector_y.col(step_in_mem-1-j).transpose()*vector_s.col(step_in_mem-1-j))(0, 0);
        
        a(0, step_in_mem-1 - j) = p*(vector_s.col(step_in_mem-1-j).transpose()*q)(0, 0);
        
        
        q = q - a(0,step_in_mem-1-j)*vector_y.col(step_in_mem-1-j);

        

    }

    
    z = MatrixB.inverse()*q;

    
    for (int j = step_in_mem-m; j!= step_in_mem; ++j){
        auto p = one/(vector_y.col(j).transpose()*vector_s.col(j))(0, 0);
        auto b = p*(vector_y.col(j).transpose()*z)(0, 0);
        z = z + vector_s.col(j)*(a(0, j)-b);
        
    }



    // Measuring memory consuption


    struct rusage usage;

    getrusage(RUSAGE_SELF, &usage);

    cout << "Memory used: " << usage.ru_maxrss << " kilobytes" << endl;






    // Updating data for next cycle

    MatrixF = MatrixF_New;
    MatrixC = MatrixC_New;
    MatrixGradient = MatrixGradient_New;


  
    
    number_iter += 1;




    auto time = steady_clock::now()-start;  // Timer End


    printf("Total Energy is      %9.9f eV\n", E0.real());
    printf("Time to step      %12d ms\n", duration_cast<milliseconds>(time).count());
    getchar();



    
    }   while (abs(deltaE0.real()) > conv_crit);

    
    Matrix_F_ = Fock_Matrix_calc(MatrixHCore, MatrixC, integrals, V_RI_1);
    SelfAdjointEigenSolver<Eigen::Matrix<complex<double>, Dynamic, Dynamic>> eigensolver(MatrixX1 * Matrix_F_ * MatrixX1);
    e = 27.2113961318065*eigensolver.eigenvalues();


    printf("Total Energy is      %9.9f eV\n", E0.real());
    printf("Iterations number:    %d\n", number_iter);
    cout << "Orbitals Energy:\n" << e << endl;
    //cout << "Coefficient matrix C:\n" << MatrixC << endl;

    return 1;
}

#endif