#include <cstdio>
#include <mpi.h>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>

#define N 300 //dimensions of matrices
#define MASTER_TO_SLAVE_TAG 1 //tag for messages sent from master to slaves
#define SLAVE_TO_MASTER_TAG 4 //tag for messages sent from slaves to master
void makeAB(); //makes the [A] and [B] matrixes
void printArray(); //print the content of output matrix [C];
void writeMatrices();

int rank; //process rank
int size; //number of processes
int i, j, k; //helper variables
double mat_a[N][N]; //declare input [A]
double mat_b[N][N]; //declare input [B]
double mat_c[N][N]; //declare output [C]
double start_time; //hold start time
double end_time; // hold end time
int low_bound; //low bound of the number of rows of [A] allocated to a slave
int upper_bound; //upper bound of the number of rows of [A] allocated to a slave
int portion; //portion of the number of rows of [A] allocated to a slave
MPI_Status status; // store status of a MPI_Recv
MPI_Request request; //capture request of a MPI_Isend

int main(int argc, char *argv[])
{

    MPI_Init(&argc, &argv); //initialize MPI operations
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //get the rank
    MPI_Comm_size(MPI_COMM_WORLD, &size); //get number of processes

    /* master initializes work*/
    if (rank == 0) {
        makeAB();
        start_time = MPI_Wtime();
        for (i = 1; i < size; i++) {//for each slave other than the master
            portion = (N / (size - 1)); // calculate portion without master
            low_bound = (i - 1) * portion;
            if (((i + 1) == size) && ((N % (size - 1)) != 0)) {//if rows of [A] cannot be equally divided among slaves
                upper_bound = N; //last slave gets all the remaining rows
            } else {
                upper_bound = low_bound + portion; //rows of [A] are equally divisable among slaves
            }
            //send the low bound first without blocking, to the intended slave
            MPI_Isend(&low_bound, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD, &request);
            //next send the upper bound without blocking, to the intended slave
            MPI_Isend(&upper_bound, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG + 1, MPI_COMM_WORLD, &request);
            //finally send the allocated row portion of [A] without blocking, to the intended slave
            MPI_Isend(&mat_a[low_bound][0], (upper_bound - low_bound) * N, MPI_DOUBLE, i, MASTER_TO_SLAVE_TAG + 2, MPI_COMM_WORLD, &request);
        }
    }
    //broadcast [B] to all the slaves
    MPI_Bcast(&mat_b, N*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /* work done by slaves*/
    if (rank > 0) {
        //receive low bound from the master
        MPI_Recv(&low_bound, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD, &status);
        //next receive upper bound from the master
        MPI_Recv(&upper_bound, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG + 1, MPI_COMM_WORLD, &status);
        //finally receive row portion of [A] to be processed from the master
        MPI_Recv(&mat_a[low_bound][0], (upper_bound - low_bound) * N, MPI_DOUBLE, 0, MASTER_TO_SLAVE_TAG + 2, MPI_COMM_WORLD, &status);
        for (i = low_bound; i < upper_bound; i++) {//iterate through a given set of rows of [A]
            for (j = 0; j < N; j++) {//iterate through columns of [B]
                for (k = 0; k < N; k++) {//iterate through rows of [B]
                    mat_c[i][j] += (mat_a[i][k] * mat_b[k][j]);
                }
            }
        }
        //send back the low bound first without blocking, to the master
        MPI_Isend(&low_bound, 1, MPI_INT, 0, SLAVE_TO_MASTER_TAG, MPI_COMM_WORLD, &request);
        //send the upper bound next without blocking, to the master
        MPI_Isend(&upper_bound, 1, MPI_INT, 0, SLAVE_TO_MASTER_TAG + 1, MPI_COMM_WORLD, &request);
        //finally send the processed portion of data without blocking, to the master
        MPI_Isend(&mat_c[low_bound][0], (upper_bound - low_bound) * N, MPI_DOUBLE, 0, SLAVE_TO_MASTER_TAG + 2, MPI_COMM_WORLD, &request);
    }

    /* master gathers processed work*/
    if (rank == 0) {
        for (i = 1; i < size; i++) {// untill all slaves have handed back the processed data
            //receive low bound from a slave
            MPI_Recv(&low_bound, 1, MPI_INT, i, SLAVE_TO_MASTER_TAG, MPI_COMM_WORLD, &status);
            //receive upper bound from a slave
            MPI_Recv(&upper_bound, 1, MPI_INT, i, SLAVE_TO_MASTER_TAG + 1, MPI_COMM_WORLD, &status);
            //receive processed data from a slave
            MPI_Recv(&mat_c[low_bound][0], (upper_bound - low_bound) * N, MPI_DOUBLE, i, SLAVE_TO_MASTER_TAG + 2, MPI_COMM_WORLD, &status);
        }
        end_time = MPI_Wtime();
        printf("\nRunning Time = %f\n\n", end_time - start_time);
        // printArray();
        writeMatrices();
    }
    MPI_Finalize(); //finalize MPI operations
    return 0;
}

void makeAB()
{
    srand(time(0));

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            mat_a[i][j] = 1.0*rand() / RAND_MAX;
        }
    }
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            mat_b[i][j] = 1.0*rand() / RAND_MAX;
        }
    }
}

void printArray()
{
    for (i = 0; i < N; i++) {
        printf("\n");
        for (j = 0; j < N; j++)
            printf("%8.2f  ", mat_a[i][j]);
    }
    printf("\n\n\n");
    for (i = 0; i < N; i++) {
        printf("\n");
        for (j = 0; j < N; j++)
            printf("%8.2f  ", mat_b[i][j]);
    }
    printf("\n\n\n");
    for (i = 0; i < N; i++) {
        printf("\n");
        for (j = 0; j < N; j++)
            printf("%8.2f  ", mat_c[i][j]);
    }
    printf("\n\n");
}


void writeMatrices()
{
    std::ofstream a_file("./input/A.txt");
    std::ofstream b_file("./input/B.txt");
    std::ofstream c_file("./output/C_mpi1.txt");

    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){
            a_file << std::to_string(mat_a[i][j]) << " ";
            b_file << std::to_string(mat_b[i][j]) << " ";
            c_file << std::to_string(mat_c[i][j]) << " ";
        }
        a_file << std::endl;
        b_file << std::endl;
        c_file << std::endl;
    }

}