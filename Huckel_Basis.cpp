#ifndef HuckelABC
#define HuckelABC
#include "header.hpp"
#include "Hamilt_guess.cpp"
#include <tuple>
#include <map>
#include <chrono>
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <fstream>
#include <libint2.hpp>
#include <typeinfo>
#include <cmath>
#include <complex>
#define Mat_init Matrix<std::complex<double>, Dynamic, Dynamic>
using namespace Eigen;
using namespace std;
using libint2::Shell;
using libint2::BasisSet;
using libint2::Operator;
using namespace libint2;


Matrix<complex<double>, Dynamic, Dynamic> MatrixS;
Matrix<complex<double>, Dynamic, Dynamic> MatrixS_min;
Matrix<complex<double>, Dynamic, Dynamic> MatrixS_trans;
Matrix<complex<double>, Dynamic, Dynamic> MatrixC;
Matrix<complex<double>, Dynamic, Dynamic> Hamiltonian;
Matrix<complex<double>, Dynamic, Dynamic> MatrixX;
Matrix<complex<double>, Dynamic, Dynamic> MatrixHCore;
Matrix<complex<double>, Dynamic, Dynamic> MatrixC_final;

int mini_nm = 0;
extern int nm;
extern string basis;
extern string file_path;
extern int Number_e;



Mat_init Extended_Huckel_Guess(){ 
    libint2::initialize();

    string ExploringMolecule = file_path;

    ifstream input_file(ExploringMolecule);

    vector<Atom> atoms = read_dotxyz(input_file);


    // Basis is here

    BasisSet obs(basis, atoms);
    BasisSet Huckel_shell("Huckel", atoms);

    

    for (int i = 0; i!=Huckel_shell.size(); ++i){
        mini_nm += Huckel_shell[i].size();

    }
  
    

    MatrixS.resize(nm, nm);
    MatrixC_final.resize(nm, nm);
    MatrixX.resize(mini_nm, mini_nm);
    MatrixC.resize(mini_nm, mini_nm);
    MatrixS_min.resize(mini_nm, mini_nm);
    MatrixHCore.resize(mini_nm, mini_nm);
    Hamiltonian.resize(mini_nm, mini_nm);
    MatrixS_trans.resize(mini_nm, nm);

    MatrixC_final.setZero();


    Engine s_engine(Operator::overlap,  // will compute overlap ints
                    15,    // max # of primitives in shells this engine will accept;
                    5         // max angular momentum of shells this engine will accept
    );

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

    auto shell2bf_1 = Huckel_shell.shell2bf(); // maps shell index to basis function index
    // shell2bf[0] = index of the first basis function in shell 0
    // shell2bf[1] = index of the first basis function in shell 1
    // ...

    const auto& buf_vec_1 = s_engine.results(); // will point to computed shell sets
    // const auto& is very important!

    for(auto s1=0; s1!=Huckel_shell.size(); ++s1) {
        for(auto s2=0; s2!=Huckel_shell.size(); ++s2) {


            s_engine.compute(Huckel_shell[s1], Huckel_shell[s2]);

            auto ints_shellset = buf_vec_1[0];  // location of the computed integrals
            if (ints_shellset == nullptr)
                continue;  // nullptr returned if the entire shell-set was screened out
            auto bf1 = shell2bf_1[s1];  // first basis function in first shell
            auto n1 = Huckel_shell[s1].size(); // number of basis functions in first shell
            auto bf2 = shell2bf_1[s2];  // first basis function in second shell
            auto n2 = Huckel_shell[s2].size(); // number of basis functions in second shell


            for(auto f1=0; f1!=n1; ++f1)
                for(auto f2=0; f2!=n2; ++f2)
                    MatrixS_min(bf1+f1, bf2+f2) = ints_shellset[f1*n2+f2];

        }
    }


    const auto& buf_vec_2 = s_engine.results(); // will point to computed shell sets
    // const auto& is very important!

    for(auto s1=0; s1!=Huckel_shell.size(); ++s1) {
            for(auto s2=0; s2!=obs.size(); ++s2) {


                s_engine.compute(Huckel_shell[s1], obs[s2]);

                auto ints_shellset = buf_vec_2[0];  // location of the computed integrals
                if (ints_shellset == nullptr)
                    continue;  // nullptr returned if the entire shell-set was screened out
                auto bf1 = shell2bf_1[s1];  // first basis function in first shell
                auto n1 = Huckel_shell[s1].size(); // number of basis functions in first shell
                auto bf2 = shell2bf[s2];  // first basis function in second shell
                auto n2 = obs[s2].size(); // number of basis functions in second shell


                for(auto f1=0; f1!=n1; ++f1)
                    for(auto f2=0; f2!=n2; ++f2)
                        MatrixS_trans(bf1+f1, bf2+f2) = ints_shellset[f1*n2+f2];

            }
        }



    
    Hamiltonian = Huckel_Hamiltonian_Guess(MatrixS_min, mini_nm, file_path);       // Guessed Hamiltonian
    





    //  Obtaining Coefficient Matrix from guessed Hamiltonian in mini_basis

    MatrixX = Diagonal_Matrix_lambda(MatrixS_min, mini_nm)[0]; // Get S^-1/2 Matrix
    
    Matrix<complex<double>, Dynamic, Dynamic> MatrixHCore = MatrixX * Hamiltonian * MatrixX;

    SelfAdjointEigenSolver<Eigen::Matrix<complex<double>, Dynamic, Dynamic>> eigensolver(MatrixHCore);
    

    MatrixC = MatrixX * eigensolver.eigenvectors();



    // Obtaining Coefficient Matrix from guessed Hamiltonian in normal basis


    
    MatrixC_final.block(0, 0, nm, mini_nm) = (MatrixS.inverse() * MatrixS_trans.transpose() * MatrixC).rowwise().reverse();

  
    

    

    return MatrixC_final;
}

#endif