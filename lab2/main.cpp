#include <iostream>
#include "mpi.h"
#include "include/utils.h"
#include <ctime>
#include <cstdlib>

#define N 10

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

    /*
    if (!is_file_exist("input/A.txt")){
        generateRandomMatrix("input/A.txt", N);
    }

    if (!is_file_exist("input/B.txt")){
        generateRandomMatrix("input/B.txt", N);
    }
    */

    readMatrix("input/A.txt", &A);
    readMatrix("input/B.txt", &B);

    auto C = serialMatrixMult(A, B);

    for (auto& row:C){
        for (auto& elem:row ){
            cout << elem << " ";
        }
        cout << "\n";
    }

    /*
    int ProcNum, ProcRank, RecvRank;

    MPI_Status Status;
    MPI_Init(&argc, &argv); //init parallel block
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    if ( ProcRank == 0 ) { // Process with Rank=0
        printf ("\n Hello from process %3d", ProcRank);
        std::cout << "Master Thread with rank = " << ProcRank << " with total number of threads of " << ProcNum << std::endl;

        for ( int i=1; i < ProcNum; i++ ) {
            MPI_Recv(&RecvRank, 1, MPI_INT, MPI_ANY_SOURCE,

                     MPI_ANY_TAG, MPI_COMM_WORLD, &Status);

            printf("\n Hello from process %3d", RecvRank);
        }
    }
    else // All other processes
        MPI_Send(&ProcRank,1,MPI_INT,0,0,MPI_COMM_WORLD);

    MPI_Finalize(); // terminate parallel block
    */

    return 0;
}
