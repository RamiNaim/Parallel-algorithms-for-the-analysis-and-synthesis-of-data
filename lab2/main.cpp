#include "include/utils.h"
#include <ctime>
#include <cstdlib>
#include <chrono>
#include <mpi.h>

#define N 200

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
    srand (static_cast <unsigned> (time(0)));

    vector<vector<double>> A, B;

    generateRandomMatrix("input/A.txt", N);
    generateRandomMatrix("input/B.txt", N);

    readMatrix("input/A.txt", &A);
    readMatrix("input/B.txt", &B);

    // std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    auto start = MPI_Wtime();
    auto C = serialMatrixMult(A, B);
    auto end = MPI_Wtime();
    // std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    // std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
    // std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count() << "[ns]" << std::endl;

    cout << "Elapsed time: " << end - start << endl;

    /*
    for (auto& row:C){
        for (auto& elem:row ){
            cout << elem << " ";
        }
        cout << "\n";
    }
    */

    return 0;
}
