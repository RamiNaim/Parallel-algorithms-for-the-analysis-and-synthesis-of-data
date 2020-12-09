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
int i, j, k, block; //helper variables
double mat_a[N][N]; //declare input [A]
double mat_b[N][N]; //declare input [B]
double mat_c[N][N]; //declare output [C]
double column_b[N];
double column_c[N];
int column_idx;
double start_time; //hold start time
double end_time; // hold end time
int low_bound; //low bound of the number of rows of [A] allocated to a slave
int upper_bound; //upper bound of the number of rows of [A] allocated to a slave
int portion; //portion of the number of rows of [A] allocated to a slave
MPI_Status status; // store status of a MPI_Recv
MPI_Request request; //capture request of a MPI_Isend

#pragma clang diagnostic push
#pragma ide diagnostic ignored "EndlessLoop"
int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv); //initialize MPI operations
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //get the rank
    MPI_Comm_size(MPI_COMM_WORLD, &size); //get number of processes

    MPI_Datatype columntype;
    MPI_Type_vector(N, 1, N, MPI_DOUBLE, &columntype);
    MPI_Type_commit(&columntype);

    /* master initializes work*/
    if (rank == 0) {
        makeAB();
    }

    //broadcast [B] to all the slaves
    MPI_Bcast(&mat_a, N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        start_time = MPI_Wtime();
        for (block = 0; block < N;) {
            for (i = 1; i < size; i++) {//for each slave other than the master
                column_idx = block + i - 1;
                //finally send the allocated row portion of [A] without blocking, to the intended slave
                MPI_Isend(&mat_b[0][column_idx], 1, columntype, i, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD, &request);
            }
            for (i = 1; i < size; i++) {//for each slave other than the master
                //finally send the allocated row portion of [A] without blocking, to the intended slave
                MPI_Recv(&column_c, N, MPI_DOUBLE, i, SLAVE_TO_MASTER_TAG, MPI_COMM_WORLD, &status);
                column_idx = block + i - 1;
                for (j = 0; j < N; j++) {
                    mat_c[j][column_idx] = column_c[j];
                }
            }
            block += size - 1;
        }
    }

    /* work done by slaves*/
    if (rank > 0) {
        for (k=0;k<N;)
        {
            //receive low bound from the master
            MPI_Recv(&column_b, N, MPI_DOUBLE, 0, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD, &status);

            // MPI_Recv(&column_idx, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG + 1, MPI_COMM_WORLD, &status);
            for (i = 0; i < N; i++) {//iterate through a given set of rows of [A]
                column_c[i] = 0;
                for (j = 0; j < N; j++) {//iterate through columns of [B]
                    column_c[i] += (mat_a[i][j] * column_b[j]);
                }
            }

            //send back the low bound first without blocking, to the master
            MPI_Send(&column_c, N, MPI_DOUBLE, 0, SLAVE_TO_MASTER_TAG, MPI_COMM_WORLD);
            k += size - 1;
        }
    }

    /* master gathers processed work*/
    if (rank == 0) {
        end_time = MPI_Wtime();
        printf("\nRunning Time = %f\n\n", end_time - start_time);
        writeMatrices();
    }
    MPI_Finalize(); //finalize MPI operations
    return 0;
}
#pragma clang diagnostic pop

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
    std::ofstream c_file("./output/C_mpi2.txt");

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

