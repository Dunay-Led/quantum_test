#include "Eigen/Dense"
#include <iostream>
using namespace Eigen;
using namespace std;
#define Mat_init Matrix<double, Dynamic, Dynamic>

void print(vector<int> mas){
    for (int i=0; i!=(sizeof(mas)/sizeof(mas[0])); ++i){
        cout << mas[i] << " ";
        
    }
}

int main(){

    vector<double> compressed_data {};
    vector<int> compressed_num {};
    vector<int> cumulant {};
    Mat_init Sparse;
    Sparse.resize(2, 2);
    Sparse.setZero();
    Sparse << 1, 2, 3, 0;
    int non_zero = 0;
    for (int i=0; i!=Sparse.col(0).size(); ++i){
        
        for (int j=0; j!=Sparse.row(0).size(); ++j){
            if (abs(Sparse(i, j)) > 0.000001){
                non_zero += 1;
                compressed_num.push_back(j);
                compressed_data.push_back(Sparse(i, j));
            }
            
        }

        cumulant.push_back(non_zero);

    }
    print(cumulant);
    cout << endl;
    print(compressed_num);




    return 1;
}