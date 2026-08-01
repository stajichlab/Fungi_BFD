#include <mpi.h>
#include <stdio.h>
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    char name[MPI_MAX_PROCESSOR_NAME]; int len;
    MPI_Get_processor_name(name, &len);
    printf("MPI_HELLO rank=%d of %d on host=%s\n", rank, size, name);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
