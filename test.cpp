#include "header.h"

int nm_aux = 0; // Number of basis functions of auxiliar basis
int nm = 0;
int Number_e;
string file_path;
string basis;


int main(int argc, char* argv[]) {

    


    // Get vital information about molecule

    ifstream data_file;
    
    data_file.open("/Users/vyasma/CLionProjects/Hartry_FockC/RHF_Mem_Free/RHF_Data.txt");

    string data_str;
    string aux_basis;

    getline(data_file, data_str);
    getline(data_file, file_path);
    getline(data_file, basis);

    Number_e = stoi(data_str);





    libint2::initialize();  // safe to use libint now .. do `libint2::initialize(true)` to produce diagnostic messages


    string ExploringMolecule = file_path; // see XYZ format description at http://en.wikipedia.org/wiki/XYZ_file_format



    ifstream input_file(ExploringMolecule);

    vector<Atom> atoms = read_dotxyz(input_file);

    // Basis is here

    BasisSet obs(basis, atoms);
    


    for (int i = 0; i!=obs.size(); ++i){
        nm += obs[i].size();
    }


    Matrix<complex<double>, Dynamic, Dynamic> MatrixS;
    Matrix<complex<double>, Dynamic, Dynamic> MatrixT;
    Matrix<complex<double>, Dynamic, Dynamic> MatrixV;
    Matrix<complex<double>, Dynamic, Dynamic> MatrixC;
    Matrix<complex<double>, Dynamic, Dynamic> MatrixF;
    Matrix<complex<double>, Dynamic, Dynamic> MatrixU;
    Matrix<complex<double>, Dynamic, Dynamic> MatrixSD;
    Matrix<complex<double>, Dynamic, Dynamic> Huckel;
    Matrix<complex<double>, Dynamic, Dynamic> Huckel_d;
    Matrix<complex<double>, Dynamic, Dynamic> V_RI_1;  // 2-electron repulsion matrix inverted

    Huckel.resize(nm, nm);
    Huckel_d.resize(nm, nm);
    MatrixS.resize(nm, nm);
    MatrixT.resize(nm, nm);
    MatrixV.resize(nm, nm);
    MatrixC.resize(nm, nm);
    MatrixF.resize(nm, nm);
    MatrixU.resize(nm, nm);
    MatrixSD.resize(nm, nm);
    V_RI_1.resize(nm_aux, nm_aux);

    #ifdef RHF_Mem
    double *integrals {new double[nm*nm*nm*nm]};
    #endif

    #ifdef RHF_Mem_RI

    getline(data_file, aux_basis);

    BasisSet obs_aux(aux_basis, atoms);

    auto shell2bf_aux = obs_aux.shell2bf();

    for (int i = 0; i!=obs_aux.size(); ++i){
        nm_aux += obs_aux[i].size();
    }

    Matrix<complex<double>, Dynamic, Dynamic> two_aux_eri;  // 2-electron repulsion matrix

    two_aux_eri.resize(nm_aux, nm_aux);
    

    double *integrals {new double[nm*nm*nm_aux]};

    #endif





    #ifdef RHF_Mem_RI_Compressed

    getline(data_file, aux_basis);

    BasisSet obs_aux(aux_basis, atoms);

    auto shell2bf_aux = obs_aux.shell2bf();

    for (int i = 0; i!=obs_aux.size(); ++i){
        nm_aux += obs_aux[i].size();
    }

    Matrix<complex<double>, Dynamic, Dynamic> two_aux_eri;  // 2-electron repulsion matrix

    two_aux_eri.resize(nm_aux, nm_aux);
    

    double *integrals {new double[nm_aux]};
    

    vector<double> compressed_data {};
    vector<int> compressed_num {};
    vector<int> cumulant ;
    //cumulant.push_back(0);
    
    

    double *abB {new double[nm_aux]};
    double *ydC {new double[nm_aux]};

    #endif





    #ifdef RHF_Mem_Free
    double integrals[1];
    #endif
    


    Engine s_engine(Operator::overlap,  // will compute overlap ints
                    obs.max_nprim(),    // max # of primitives in shells this engine will accept;
                    obs.max_l()         // max angular momentum of shells this engine will accept
    );


    Engine t_engine(Operator::kinetic,  // will compute overlap ints
                    obs.max_nprim(),    // max # of primitives in shells this engine will accept
                    obs.max_l()         // max angular momentum of shells this engine will accept
    );
    Engine v_engine(Operator::nuclear,  // will compute overlap ints
                    obs.max_nprim(),    // max # of primitives in shells this engine will accept
                    obs.max_l()         // max angular momentum of shells this engine will accept
    );
    Engine eri_engine(Operator::coulomb,  // will compute overlap ints
                    obs.max_nprim(),    // max # of primitives in shells this engine will accept
                    obs.max_l()         // max angular momentum of shells this engine will accept
    );

    v_engine.set_params(make_point_charges(atoms));
    // can repeat the libint2::initialize() ... finalize() cycle as many times as
    // necessary
    auto shell2bf = obs.shell2bf(); // maps shell index to basis function index
    // shell2bf[0] = index of the first basis function in shell 0
    // shell2bf[1] = index of the first basis function in shell 1
    // ...

    const auto& buf_vec = s_engine.results(); // will point to computed shell sets
    // const auto& is very important!


    for(auto s1=0; s1!=obs.size(); ++s1) {
        for(auto s2=0; s2!=obs.size(); ++s2) {


            s_engine.compute(obs[s1], obs[s2]);

            auto ints_shellset = buf_vec[0];  // location of the computed integrals
            if (ints_shellset == nullptr)
                continue;  // nullptr returned if the entire shell-set was screened out
            auto bf1 = shell2bf[s1];  // first basis function in first shell
            auto n1 = obs[s1].size(); // number of basis functions in first shell
            auto bf2 = shell2bf[s2];  // first basis function in second shell
            auto n2 = obs[s2].size(); // number of basis functions in second shell


            for(auto f1=0; f1!=n1; ++f1)
                for(auto f2=0; f2!=n2; ++f2)
                    MatrixS(bf1+f1, bf2+f2) = ints_shellset[f1*n2+f2];

        }
    }

    const auto& Tbuf_vec = t_engine.results();
    for(auto s1=0; s1!=obs.size(); ++s1) {
        for(auto s2=0; s2!=obs.size(); ++s2) {


            t_engine.compute(obs[s1], obs[s2]);

            auto ints_shellset = Tbuf_vec[0];  // location of the computed integrals
            if (ints_shellset == nullptr)
                continue;  // nullptr returned if the entire shell-set was screened out

            auto bf1 = shell2bf[s1];  // first basis function in first shell
            auto n1 = obs[s1].size(); // number of basis functions in first shell
            auto bf2 = shell2bf[s2];  // first basis function in second shell
            auto n2 = obs[s2].size(); // number of basis functions in second shell

            // integrals are packed into ints_shellset in row-major (C) form
            // this iterates over integrals in this order
            for(auto f1=0; f1!=n1; ++f1)
                for(auto f2=0; f2!=n2; ++f2)

                    MatrixT(bf1+f1, bf2+f2) = ints_shellset[f1*n2+f2];
        }
    }
    
    const auto& Vbuf_vec = v_engine.results();
    for(auto s1=0; s1!=obs.size(); ++s1) {
        for(auto s2=0; s2!=obs.size(); ++s2) {


            v_engine.compute(obs[s1], obs[s2]);

            auto ints_shellset = Vbuf_vec[0];  // location of the computed integrals
            if (ints_shellset == nullptr)
                continue;  // nullptr returned if the entire shell-set was screened out

            auto bf1 = shell2bf[s1];  // first basis function in first shell
            auto n1 = obs[s1].size(); // number of basis functions in first shell
            auto bf2 = shell2bf[s2];  // first basis function in second shell
            auto n2 = obs[s2].size(); // number of basis functions in second shell

            // integrals are packed into ints_shellset in row-major (C) form
            // this iterates over integrals in this order
            for(auto f1=0; f1!=n1; ++f1)
                for(auto f2=0; f2!=n2; ++f2)

                    MatrixV(bf1+f1, bf2+f2) = ints_shellset[f1*n2+f2];
        }
    }


    #ifdef RHF_Mem

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

                                    integrals[nm*nm*nm*(bf1_first+f1)+nm*nm*(bf2_first+f2)+nm*(bf3_first+f3)+(bf4_first+f4)] = buf_1234[f1234];


                //printf("%d_%d_%d_%d = %6.6f", bf1, bf2, bf3, bf4, buf_1234[f1234]);
                //getchar();
                                
                                }
                            }
                        }
                    }



                }
            }
        }
    }

    #endif

    #ifdef RHF_Mem_RI
    double  c = 0;
    double c1 = 0;

    Engine three_eri(Operator::coulomb,obs_aux.max_nprim(),    // max # of primitives in shells this engine will accept
                    obs_aux.max_l());
    three_eri.set(BraKet::xs_xx);

    const auto& buf = three_eri.results();

    
    for(auto s1=0; s1!=obs.size(); ++s1) {

        auto bf1_first = shell2bf[s1]; // first basis function in this shell
        auto n1 = obs[s1].size();

        for(auto s2=0; s2!=obs.size(); ++s2) {

            auto bf2_first = shell2bf[s2];
            auto n2 = obs[s2].size();

            // loop over shell pairs of the density matrix, {s3,s4}
            // again symmetry is not used for simplicity
            for(auto s_aux=0; s_aux!=obs_aux.size(); ++s_aux) {

                auto bf3_first = shell2bf_aux[s_aux];
                auto n3 = obs_aux[s_aux].size();

                    // Coulomb contribution to the Fock matrix is from {x(a),x(b),aux(A)} integrals
                    three_eri.compute(obs_aux[s_aux], obs[s1], obs[s2]);
                    const auto* buf_1234 = buf[0];
                    if (buf_1234 == nullptr)
                        continue; // if all integrals screened out, skip to next quartet

                    for(auto f3=0, f1234=0; f3!=n3; ++f3) {
                        const auto bf3 = f3 + bf3_first;
                        for(auto f1=0; f1!=n1; ++f1) {
                            const auto bf1 = f1 + bf1_first;
                            for(auto f2=0; f2!=n2; ++f2, ++f1234) {
                                const auto bf2 = f2 + bf2_first;
                            
                                
                                    integrals[nm*nm_aux*(bf1_first+f1)+nm_aux*(bf2_first+f2)+(bf3_first+f3)] = buf_1234[f1234];
                                    //printf("%d_%d_%d = %6.6f", bf1, bf2, bf3, buf_1234[f1234]);
                                    //getchar();
                                    if (abs(buf_1234[f1234]) < 0.001){

                                        c += 1;

                                    }
                                    c1 += 1;
                                }
                            }
                        }
                    }



                }
            }

        cout << c/c1;
        getchar();
    

    Engine two_eri(Operator::coulomb,obs_aux.max_nprim(),    // max # of primitives in shells this engine will accept
                    obs_aux.max_l());
    two_eri.set(BraKet::xs_xs);

    const auto& buf_aux = two_eri.results();
    for(auto s1=0; s1!=obs_aux.size(); ++s1) {
        for(auto s2=0; s2!=obs_aux.size(); ++s2) {


            two_eri.compute(obs_aux[s1], obs_aux[s2]);

            auto ints_shellset = buf_aux[0];  // location of the computed integrals
            if (ints_shellset == nullptr)
                continue;  // nullptr returned if the entire shell-set was screened out

            auto bf1 = shell2bf_aux[s1];  // first basis function in first shell
            auto n1 = obs_aux[s1].size(); // number of basis functions in first shell
            auto bf2 = shell2bf_aux[s2];  // first basis function in second shell
            auto n2 = obs_aux[s2].size(); // number of basis functions in second shell

            // integrals are packed into ints_shellset in row-major (C) form
            // this iterates over integrals in this order
            for(auto f1=0; f1!=n1; ++f1){
                for(auto f2=0; f2!=n2; ++f2){

                    two_aux_eri(bf1+f1, bf2+f2) = ints_shellset[f1*n2+f2];
                }
            }
    
        }
    }
    
    V_RI_1 = two_aux_eri.inverse(); // inverting V - RI matri
c = 0;
c1 = 0;
for(auto f1=0; f1!=V_RI_1.col(0).size(); ++f1){
            for(auto f2=0; f2!=V_RI_1.col(0).size(); ++f2){

                if (abs(V_RI_1(f1, f2).real())<0.0001){
                    c += 1;
                }
    c1+=1;
                }
}
cout << c/c1;


    #endif



#ifdef RHF_Mem_RI_Compressed

    double  c = 0;
    double c1 = 0;

    Engine three_eri(Operator::coulomb,obs_aux.max_nprim(),    // max # of primitives in shells this engine will accept
                    obs_aux.max_l());
    three_eri.set(BraKet::xs_xx);

    const auto& buf = three_eri.results();

    int non_zero = 0;

   
        
    for(auto s1=0; s1!=obs.size(); ++s1) {

        auto bf1_first = shell2bf[s1]; // first basis function in this shell
        auto n1 = obs[s1].size();

        for(auto f1=0; f1!=n1; ++f1) {
                    const auto bf1 = f1 + bf1_first;


            for(auto s2=0; s2!=obs.size(); ++s2) {

                auto bf2_first = shell2bf[s2];
                auto n2 = obs[s2].size();


            
                    for(auto f2=0; f2!=n2; ++f2) {
                        const auto bf2 = f2 + bf2_first;



                 for(auto s_aux=0; s_aux!=obs_aux.size(); ++s_aux) {

                    auto bf3_first = shell2bf_aux[s_aux];
                    auto n3 = obs_aux[s_aux].size();


                    three_eri.compute(obs_aux[s_aux], obs[s1], obs[s2]);

                    const auto* buf_1234 = buf[0];
                    int f1234 = 0;
                    if (buf_1234 == nullptr)
                        continue; // if all integrals screened out, skip to next quartet

                    

                            for(auto f3=0; f3!=n3; ++f3, ++f1234) {
                                const auto bf3 = f3 + bf3_first;
                        
                                
                                integrals[bf3] = buf_1234[f1234*n1*n2+f1*n2+f2]; 
                                cout << bf3 << endl;

                                if (abs(buf_1234[f1234*n1*n2+f1*n2+f2]) > 0.000001){

                                    c += 1;

                                }
                                c1 += 1;
                                
                                }
                                
                                for (int i = 0 ; i!= nm_aux; ++i){
                                    cout << integrals[i];
                                }
                                
                            }
                            
                           

                        for (int i = 0; i != nm_aux; ++i){

                    if (abs(integrals[i]) > 0.000001){
                        
                        ++non_zero;

                        compressed_data.push_back(integrals[i]);
                        compressed_num.push_back(i);
                        

                        
                    }
                   

                    

                }

                cumulant.push_back(non_zero);

                        }
                    }
                   

            cout <<"End << endl" << endl;
                            getchar();
                            

                
                

                }

                


            }
        

            print(cumulant);
    cout << endl;
    print(compressed_num);
    getchar();

        cout << c/c1;
        getchar();
        
    cout << endl <<  c/c1 << endl;

    Engine two_eri(Operator::coulomb,obs_aux.max_nprim(),    // max # of primitives in shells this engine will accept
                    obs_aux.max_l());
    two_eri.set(BraKet::xs_xs);

    const auto& buf_aux = two_eri.results();
    for(auto s1=0; s1!=obs_aux.size(); ++s1) {
        for(auto s2=0; s2!=obs_aux.size(); ++s2) {


            two_eri.compute(obs_aux[s1], obs_aux[s2]);

            auto ints_shellset = buf_aux[0];  // location of the computed integrals
            if (ints_shellset == nullptr)
                continue;  // nullptr returned if the entire shell-set was screened out

            auto bf1 = shell2bf_aux[s1];  // first basis function in first shell
            auto n1 = obs_aux[s1].size(); // number of basis functions in first shell
            auto bf2 = shell2bf_aux[s2];  // first basis function in second shell
            auto n2 = obs_aux[s2].size(); // number of basis functions in second shell

            // integrals are packed into ints_shellset in row-major (C) form
            // this iterates over integrals in this order
            for(auto f1=0; f1!=n1; ++f1){
                for(auto f2=0; f2!=n2; ++f2){

                    two_aux_eri(bf1+f1, bf2+f2) = ints_shellset[f1*n2+f2];
                }
            }
    
        }
    }
    
    V_RI_1 = two_aux_eri.inverse(); // inverting V - RI matri
//cout << V_RI_1;
c = 0;
c1 = 0;
for(auto f1=0; f1!=V_RI_1.col(0).size(); ++f1){
            for(auto f2=0; f2!=V_RI_1.col(0).size(); ++f2){

                if (abs(V_RI_1(f1, f2).real())<0.0001){
                    c += 1;
                }
    c1+=1;
                }
}
cout << c/c1;


    #endif
    







    auto MatrixX = Diagonal_Matrix_lambda(MatrixS, nm)[0]; // Get S^-1/2 Matrix


    // First Guess of F

    MatrixF = MatrixV + MatrixT;

    Matrix<complex<double>, Dynamic, Dynamic> MatrixHCore = MatrixF;

    Matrix<complex<double>, Dynamic, Dynamic> _MatrixHCore_ = MatrixX*MatrixF*MatrixX;
    _MatrixHCore_.resize(nm, nm);

    SelfAdjointEigenSolver<Eigen::Matrix<complex<double>, Dynamic, Dynamic>> eigensolver(_MatrixHCore_);
    MatrixC = MatrixX*eigensolver.eigenvectors();


    MatrixC = Extended_Huckel_Guess();


    Huckel = Fock_Matrix_calc(MatrixHCore, MatrixC, integrals, V_RI_1);
    Huckel_d = MatrixX * Huckel * MatrixX;
    SelfAdjointEigenSolver<Eigen::Matrix<complex<double>, Dynamic, Dynamic>> eigensolver1(Huckel_d);
    MatrixC = MatrixX*eigensolver1.eigenvectors();

    
    SOSCF_calc(MatrixC, MatrixHCore, MatrixX, integrals, V_RI_1); // Main SCF procedure

    


    libint2::finalize();

    return 0;
}