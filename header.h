//
// Created by Vyasma on 09.10.2023.
//

#ifndef HARTRY_FOCKC_HEADER_H
#define HARTRY_FOCKC_HEADER_H

#define RHF_Mem_RI
#include <string>
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
#include "Schmidt_orthogonal.cpp"
#include "L_diag.cpp"
#include "P_Matrix.cpp"
#include "Fock_Matrix.cpp"
#include "SCF_Simple_iter.cpp"
#include "SOSCF.cpp"
#include "Huckel_Basis.cpp"
using namespace libint2;
using namespace Eigen;
using namespace std;
using namespace std::chrono;
using libint2::Shell;
using libint2::BasisSet;
using libint2::Operator;




#endif //HARTRY_FOCKC_HEADER_H
