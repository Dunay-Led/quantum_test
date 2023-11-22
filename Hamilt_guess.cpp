#include <iostream>
#include <fstream>
#include <string>
#include <complex>
#include "Eigen/Dense"
using namespace std;
using namespace Eigen;




// Function returns guessed Hamilonian


auto Huckel_Hamiltonian_Guess(Matrix<complex<double>, Dynamic, Dynamic> MatrixS1, int nm1, string file_path){

Matrix<complex<double>, Dynamic, Dynamic> Hamiltonian;
Hamiltonian.resize(nm1, nm1);
ifstream atom_file(file_path);
string atom_name;
string atom_total;
complex<double> max_ham;
int count = 0;

// Ionisation potentials data

double data_C[5] = {11.3, 0.706, 0.433, 0.433, 0.433};
double data_N[5] = {15.6, 0.945, 0.568, 0.568, 0.568};
double data_O[5] = {20.7, 1.244, 0.632, 0.632, 0.632};
double data_F[5] = {26.4, 1.573, 0.730, 0.730, 0.730};
double data_H[1] = {0.500};



getline(atom_file, atom_total);     //  Nuber of atoms
getline(atom_file, atom_name);      //  Empty string




    //  Diagonal elements of Hamiltonian

for (int i = 0; i!=stoi(atom_total); ++i){

    getline(atom_file, atom_name);  // Get information about atom type


    if (atom_name[0] == 'C'){

    for (int j =0; j!=5; ++j){

        Hamiltonian(count, count) = data_C[j];
        ++count;
        }

    }

    if (atom_name[0] == 'F'){

    for (int j =0; j!=5; ++j){

        Hamiltonian(count, count) = data_F[j];
        ++count;
        }

    }


    if (atom_name[0] == 'N'){

    for (int j =0; j!=5; ++j){

        Hamiltonian(count, count) = data_N[j];
        ++count;
        }

    }

    if (atom_name[0] == 'O'){

    for (int j =0; j!=5; ++j){

        Hamiltonian(count, count) = data_O[j];
        ++count;
        }

    }

    


    if (atom_name[0] == 'H'){

    for (int j =0; j!=1; ++j){

        Hamiltonian(count, count) = data_H[j];
        ++count;

        }

    }

}






    // All non-diagonal Hamiltonian elements

for (int i =0; i!=nm1; ++i){
    for (int j =0; j!=nm1; ++j){

        if (i == j){
            continue;
        }

        else{

        if (Hamiltonian(j, j).real() > Hamiltonian(i, i).real()){

            max_ham = Hamiltonian(j, j);

        }

        else{

            max_ham = Hamiltonian(i, i);
        }
         

        if (max_ham.real()  > 2 ){


            Hamiltonian(j, i)  = 0.5 * 1.75 * 0.05 * MatrixS1(j, i)*(Hamiltonian(i, i)+Hamiltonian(j, j));
            
        }

        else{

        Hamiltonian(j, i)  = 0.5 * 1.75 * MatrixS1(j, i)*(Hamiltonian(i, i)+Hamiltonian(j, j));

        

        }
     
        }
          

    }

    
}




return Hamiltonian;

}