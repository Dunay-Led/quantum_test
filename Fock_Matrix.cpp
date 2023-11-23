#ifndef FockC
#define FockC
#include <libint2.hpp>
#include "header.h"
#define Mat_init Matrix<std::complex<double>, Dynamic, Dynamic>
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

auto Fock_Matrix_calc(Mat_init MatrixHCore, Mat_init MatrixC, double *integrals, Mat_init V_RI_1){

    
    Matrix<complex<double>, Dynamic, Dynamic> MatrixF;
    Matrix<complex<double>, Dynamic, Dynamic> MatrixG;
    Matrix<complex<double>, Dynamic, Dynamic> MatrixP;
    
    MatrixP.resize(nm, nm);
    MatrixG.resize(nm, nm);
    MatrixF.resize(nm, nm);
    complex<double> a_b_y_d;    // (ab|yd)
    
    
    

    MatrixG.setZero();
    

    MatrixP = P_Matrix_form_C_count(MatrixC);

    #ifdef RHF_Mem_Free

    string ExploringMolecule = file_path;

    ifstream input_file(ExploringMolecule);

    vector<Atom> atoms = read_dotxyz(input_file);

    // Basis is here

    BasisSet obs(basis, atoms);

    auto shell2bf = obs.shell2bf();

    Engine eri_engine(Operator::coulomb,  // will compute overlap ints
                      obs.max_nprim(),    // max # of primitives in shells this engine will accept
                      obs.max_l()         // max angular momentum of shells this engine will accept
    );



        const auto& buf = eri_engine.results();

        // loop over shell pairs of the Fock matrix, {s1,s2}
        // Fock matrix is symmetric, but skipping it here for simplicity (see compute_2body_fock)
        for(auto s1=0; s1!=obs.size(); ++s1) {

            auto bf1_first = shell2bf[s1]; // first basis function in this shell
            auto n1 = obs[s1].size();

            for(auto s2=0; s2!=obs.size(); ++s2) {

                auto bf2_first = shell2bf[s2];
                auto n2 = obs[s2].size();

                // loop over shell pairs of the density matrix, {s3,s4}
                // again symmetry is not used for simplicity
                for(auto s3=0; s3!=obs.size(); ++s3) {

                    auto bf3_first = shell2bf[s3];
                    auto n3 = obs[s3].size();

                    for(auto s4=0; s4!=obs.size(); ++s4) {

                        auto bf4_first = shell2bf[s4];
                        auto n4 = obs[s4].size();

                        // Coulomb contribution to the Fock matrix is from {s1,s2,s3,s4} integrals
                        eri_engine.compute(obs[s1], obs[s2], obs[s3], obs[s4]);
                        const auto* buf_1234 = buf[0];
                        if (buf_1234 == nullptr)
                            continue; // if all integrals screened out, skip to next quartet

                        for(auto f1=0, f1234=0; f1!=n1; ++f1) {
                            const auto bf1 = f1 + bf1_first;
                            for(auto f2=0; f2!=n2; ++f2) {
                                const auto bf2 = f2 + bf2_first;
                                for(auto f3=0; f3!=n3; ++f3) {
                                    const auto bf3 = f3 + bf3_first;
                                    for(auto f4=0; f4!=n4; ++f4, ++f1234) {
                                        const auto bf4 = f4 + bf4_first;


                                        MatrixG(bf1,bf2) += MatrixP(bf3,bf4) * 2.0 * buf_1234[f1234];
                                        MatrixG(bf1,bf4) -= MatrixP(bf3,bf2) * buf_1234[f1234];

                                       
                                    }
                                }
                            }
                        }



                    }
                }
            }
        }

    #endif

    

    #ifdef RHF_Mem

    for (int phi_first = 0; phi_first != nm; ++phi_first){
        for (int psi_first = 0; psi_first != nm; ++psi_first){
            for (int phi_second = 0; phi_second != nm; ++phi_second){
                for (int psi_second = 0; psi_second != nm; ++psi_second){

                    
                    MatrixG(phi_first, psi_first) += 2.0*MatrixP(psi_second, phi_second)*integrals[nm*nm*nm*phi_first + nm*nm*psi_first + nm*phi_second + psi_second];
                    MatrixG(phi_first, psi_second) -= MatrixP(psi_first, phi_second)*integrals[nm*nm*nm*phi_first + nm*nm*psi_first + nm*phi_second + psi_second];


                    //cout << MatrixG << endl;
                    //printf("%d_%d_%d_%d = %6.6f", phi_first, psi_first, phi_second, psi_second, integrals[nm*nm*nm*phi_first + nm*nm*psi_first + nm*phi_second + psi_second]);
                    //getchar();

                }
            }
        }
    }



    #endif

    #ifdef RHF_Mem_RI


    for (int phi_first = 0; phi_first != nm; ++phi_first){
        for (int psi_first = 0; psi_first != nm; ++psi_first){
            for (int phi_second = 0; phi_second != nm; ++phi_second){
                for (int psi_second = 0; psi_second != nm; ++psi_second){

                a_b_y_d = 0;

                for (int B = 0; B != nm_aux; ++B){
                    for (int C = 0; C != nm_aux; ++C){
                
                        a_b_y_d += integrals[nm*nm_aux*phi_first+nm_aux*psi_first+B] * V_RI_1(B, C) * integrals[nm*nm_aux*phi_second+nm_aux*psi_second+C];
        
                    
                    }
                }

                //printf("%d_%d_%d_%d = %6.6f", phi_first, psi_first, phi_second, psi_second, a_b_y_d);
                //getchar();

                MatrixG(phi_first, psi_first) += 2.0*MatrixP(psi_second, phi_second)*a_b_y_d;
                MatrixG(phi_first, psi_second) -= MatrixP(psi_first, phi_second)*a_b_y_d;

                     ////cout << MatrixG << endl;
                     //getchar();


                }
            }
        }
    }
   

    #endif
    

    #ifdef RHF_Mem_RI_Compressed


    int start_point;
    cout << 1;
    int non_zero_num;
    
    int col_num;
     
    double element_data;
    cout << 1;


    for (int phi_first = 0; phi_first != nm; ++phi_first){
        for (int psi_first = 0; psi_first != nm; ++psi_first){
            for (int phi_second = 0; phi_second != nm; ++phi_second){
                for (int psi_second = 0; psi_second != nm; ++psi_second){

                a_b_y_d = 0;
cout << 1;
                for (int i = 0; i!= nm_aux; ++i){
                    abB[i] = 0;
                    ydC[i] = 0;
                }
                
           cout << 1;
                start_point = cumulanta[psi_first+phi_first*nm];
                
                
                non_zero_num = cumulanta[psi_first+phi_first*nm+1] - cumulanta[psi_first+phi_first*nm];



                for (int i = 0; i!= non_zero_num; ++i){

                
                    
                    col_num = comp_num[start_point + i];
                
                
                    
                    element_data = comp_data[start_point + i];
                    abB[col_num] = element_data;

                }
cout << 2;
                start_point = cumulanta[phi_second*nm + psi_second];
                cout << 2;
                

                non_zero_num = cumulanta[phi_second + psi_second*nm+1] - cumulanta[phi_second + psi_second*nm];
                cout << 2;

                for (int i = 0; i!= non_zero_num; ++i){
                    

                    col_num = comp_num[start_point + i];
                    
                    element_data = comp_data[start_point + i];
                    
                    abB[col_num] = element_data;

                }
           

                for (int B = 0; B != nm_aux; ++B){
                    for (int C = 0; C != nm_aux; ++C){
                
                        a_b_y_d += abB[B] * V_RI_1(B, C) * ydC[C];
        
                    
                    }
                }


                MatrixG(phi_first, psi_first) += 2.0*MatrixP(psi_second, phi_second)*a_b_y_d;
                MatrixG(phi_first, psi_second) -= MatrixP(psi_first, phi_second)*a_b_y_d;

                
                }
            }
        }
    }
   

    #endif
    

    MatrixF = MatrixHCore + MatrixG;


    //cout << MatrixF << endl;
    //getchar();

   

    
    return MatrixF;
}

#endif