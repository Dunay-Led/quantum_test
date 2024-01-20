#ifndef FockC
#define FockC
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

Mat_init Fock_Matrix_calc(Mat_init MatrixHCore, Mat_init MatrixC, Mat_init_r integrals, Mat_init V_RI_1, vector<int> cumulant, vector<int> compressed_num, vector<double> compressed_data){

    
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

                    
                    MatrixG(phi_first, psi_first) += 2.0*MatrixP(psi_second, phi_second)*integrals(0, nm*nm*nm*phi_first + nm*nm*psi_first + nm*phi_second + psi_second);
                    MatrixG(phi_first, psi_second) -= MatrixP(psi_first, phi_second)*integrals(0, nm*nm*nm*phi_first + nm*nm*psi_first + nm*phi_second + psi_second);


                    
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
                
                        a_b_y_d += integrals(0, nm*nm_aux*phi_first+nm_aux*psi_first+B) * V_RI_1(B, C) * integrals(0, nm*nm_aux*phi_second+nm_aux*psi_second+C);
        
                    
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
    

    #ifdef RHF_Mem_RI_Compressed_mb


    int start_point;
    
    int non_zero_num;
    
    int col_num;
     
    double element_data;
    


    for (int phi_first = 0; phi_first != nm; ++phi_first){
        for (int psi_first = 0; psi_first != nm; ++psi_first){
            for (int phi_second = 0; phi_second != nm; ++phi_second){
                for (int psi_second = 0; psi_second != nm; ++psi_second){

                a_b_y_d = 0;

                for (int i = 0; i!= nm_aux; ++i){
                    abB[i] = 0;
                    ydC[i] = 0;
                }
                
           
                start_point = cumulanta[psi_first+phi_first*nm];
                
                
                non_zero_num = cumulanta[psi_first+phi_first*nm+1] - cumulanta[psi_first+phi_first*nm];



                for (int i = 0; i!= non_zero_num; ++i){

                
                    
                    col_num = comp_num[start_point + i];
                
                
                    
                    element_data = comp_data[start_point + i];
                    abB[col_num] = element_data;

                }

                start_point = cumulanta[phi_second*nm + psi_second];
                
                

                non_zero_num = cumulanta[phi_second + psi_second*nm+1] - cumulanta[phi_second + psi_second*nm];
                

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
 
    #ifdef RHF_Mem_RI_Compressed

    int n_occ = Number_e/2;
    Mat_init i_a_B;
    i_a_B.resize(nm*n_occ, nm_aux);
    i_a_B.setZero();

    for (int f_aux = 0; f_aux != nm*nm; ++f_aux){

    int a = f_aux % nm;

    for (int first = cumulant[f_aux]; first != cumulant[f_aux+1]; ++first){

        for (int i = 0; i != n_occ; ++i){


       
            i_a_B(i * nm+a, compressed_num[first]) += MatrixC(f_aux / nm, i) * compressed_data[first];


            }
        }
    }
 
    for (int alpha = 0; alpha != nm; ++alpha){
        for (int beta = 0; beta != nm; ++beta){

 
            MatrixG(alpha, beta) = J_contribution_calc(MatrixP, MatrixC,  V_RI_1, integrals, alpha, beta, cumulant, compressed_num, compressed_data, i_a_B);


        }
    }
    
















    #endif
    

    MatrixF = MatrixHCore + MatrixG;



    
    return MatrixF;
}

#endif