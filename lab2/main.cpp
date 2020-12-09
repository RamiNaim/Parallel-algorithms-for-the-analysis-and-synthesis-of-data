#include "include/utils.h"
#include <ctime>
#include <cstdlib>
#include <mpi.h>

#define N 300

using namespace std;

vector<vector<double>> serialMatrixMult(vector<vector<double>>& A,  vector<vector<double>>& B){
    vector<vector<double>> C;
    vector<double> row;
    double elem;

    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){
            elem = 0;
            for (int k=0;k<N;k++){
                elem += A.at(i).at(k) * B.at(k).at(j);
            }
            row.push_back(elem);
        }
        C.push_back(row);
        row.clear();
    }

    return C;
}

int main(int argc, char* argv[]) {
    srand (time(0));

    vector<vector<double>> A, B;

    readMatrix("./input/A.txt", &A);
    readMatrix("./input/B.txt", &B);

    auto start = MPI_Wtime();
    auto C = serialMatrixMult(A, B);
    auto end = MPI_Wtime();

    cout << "Elapsed time: " << end - start << endl;

    std::ofstream c_file("./output/C_serial.txt");

    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){
            c_file << std::to_string(C[i][j]) << " ";
        }
        c_file << std::endl;
    }

    return 0;
}
