#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <sys/stat.h>
#include <limits.h>
#include <assert.h>
#include <time.h>
#include <mpi.h>

struct data
{
    long int intervals;
    int min, max;
    float* endpoints;
    unsigned long long int* occurences;
    unsigned long long int* local_occurences;
    unsigned long long int size;
    unsigned long long int local_size;
    int* buffer;
    int* local_buffer;

    // Pointers for MPI
    long int* interval;
} node;

void construct(data node)
{
    intervals = malloc(sizeof(long int));
}

void initialize(data node, char* s, char* filename, int my_rank, MPI_Comm comm, int comm_sz);
{
    if(my_rank == 0)
    {
        char* temp;
        node.intervals = strtol(s, &temp, 10);
        node.interval = &intervals;

    }
}
int main(int argc, char* argv[])
{
    int my_rank;
    int comm_sz;
    MPI_Comm comm;

    //Call constructor to initialize values to arrays.
    construct(data node);

    //Initialize mpi
    e(MPI_Init(&argc, &argv));
    comm = MPI_COMM_WORLD;
    e(MPI_Comm_size(comm, &comm_sz));
    e(MPI_Comm_rank(comm, &my_rank));

    initialize(node, argv[1], argv[2], my_rank, comm, comm_sz);
    return 0;
}
