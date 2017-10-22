#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <sys/stat.h>
#include <limits.h>
#include <assert.h>
#include <time.h>
#include <mpi.h>

struct Data
{
    long int intervals;
    int min;
    int max;
    float* endpoints;
    unsigned long long int* occurences;
    unsigned long long int* local_occurences;
    unsigned long long int size;
    unsigned long long int local_size;
    int* buffer;
    int* local_buffer;

} node;

void construct(struct Data node)
{
    node.endpoints = malloc( node.intervals * sizeof(float));
    node.local_occurences = malloc( node.intervals * sizeof(unsigned long long int));
    node.occurences = malloc( node.intervals * sizeof(unsigned long long int));
    node.buffer = malloc (node.size * sizeof(int));
    node.local_buffer = malloc (node.size * sizeof(int));
}

void destroy(struct Data node)
{
    free(node.endpoints);
    free(node.local_occurences);
    free(node.occurences);
    free(node.buffer);
    free(node.local_buffer);
}

void endpoints(struct Data node)
{
    float length = (node.max - node.min) / (float) node.intervals;
    float temp = node.min;
    for(size_t i = 0; i < node.intervals; i++)
    {
        node.endpoints[i] = temp;
        temp += length;
    }
}

size_t determine_index(int temp, float* endpoints, long int intervals)
{
    assert(endpoints != NULL);
    size_t index;
    for( index =0; index < intervals-1; index++)
    {
        if( temp <= endpoints[index])
        {
            break;
        }
    }
    return index;
}

void count_occurences(struct Data node)
{
    assert (node.buffer != NULL);
    assert (node.endpoints != NULL);

    #pragma omp parallel for
    for(int i=0; i< node.local_size; i++)
    {
        size_t index = determine_index(node.local_buffer[i], node.endpoints, node.intervals);
        node.occurences[index]++;
    }
}

void e(int error)
{
    if(error != MPI_SUCCESS)
    {
        fprintf(stderr, "Error starting MPI program! \n");
        MPI_Abort(MPI_COMM_WORLD, error);
        MPI_Finalize();
        exit(1);
    }
}

void initialize(struct Data* node, char* s, char* filename, int my_rank, MPI_Comm comm, int comm_sz)
{
    if(my_rank == 0)
    {
        char* temp;
        node->intervals = strtol(s, &temp, 10);

        FILE *fp;
        node->min = INT_MAX;
        node->max = INT_MIN;
        struct stat file_stat;
        unsigned long long int amount;

        fp = fopen(filename, "r");
        if(fp == NULL)
        {
            printf("\nFile doesn't exist");
            exit(0);
        }
        int result = stat(filename, &file_stat);
        if(result == -1)
        {
            printf("\n File invalid.");
            exit(0);
        }
        node->size = file_stat.st_size;
        node->size /= sizeof(int);
        node->local_size = node->size / comm_sz;
        node->buffer = malloc ( node->size *sizeof(int));
        if(node->buffer)
        {
            amount = fread(node->buffer, sizeof(int), node->size, fp);
            if(amount == 0)
            {
                printf("\nCouldn't read the file.");
                exit(0);
            }
        }
        else
        {
            printf("\nMalloc didn't succed.");
        }
        for(unsigned long long int i = 0; i < node->size; i++)
		{
			if((node->buffer[i]) < node->min)
			{
				node->min = node->buffer[i];
			}
			if((node->buffer[i]) > node->max)
			{
				node->max = node->buffer[i];
			}
		}
    }
    long int* interval = &node->intervals;
    int* mins = &node->min;
    int* maxs = &node->max;
    unsigned long long int* sizes = &node->size;
    unsigned long long int* local_sizes = &node->local_size;
    e(MPI_Bcast(interval, 1, MPI_LONG_INT, 0, comm));
	e(MPI_Bcast(mins, 1, MPI_INT, 0, comm));
	e(MPI_Bcast(maxs, 1, MPI_INT, 0, comm));
	e(MPI_Bcast(sizes, 1, MPI_LONG_LONG_INT, 0, comm));
	e(MPI_Bcast(local_sizes, 1, MPI_LONG_LONG_INT, 0, comm));
    if( my_rank == 0)
    {
        free(node->buffer);
    }
}

void display_histogram(struct Data node)
{
    assert (node.occurences != NULL);
    assert (node.endpoints != NULL);

    float length = node.endpoints[1] - node.endpoints[0];
    for(size_t i = 0; i < node.intervals ; ++i)
    {
        printf("%f - %f", node.endpoints[i], node.endpoints[i] + length );
        printf("            %lld\n", node.occurences[i] );
    }
}
int main(int argc, char* argv[])
{
    int my_rank;
    int comm_sz;
    MPI_Comm comm;

    //Initialize mpi
    e(MPI_Init(&argc, &argv));
    comm = MPI_COMM_WORLD;
    e(MPI_Comm_size(comm, &comm_sz));
    e(MPI_Comm_rank(comm, &my_rank));

    initialize(&node, argv[1], argv[2], my_rank, comm, comm_sz);
    construct(node);

    endpoints(node);
    if( my_rank == 0)
    {
        FILE *fp;
        unsigned long long int amount;

        fp = fopen(argv[2], "r");
        node.buffer = malloc ( node.size * sizeof(int));
        if(node.buffer)
        {
            amount = fread(node.buffer, sizeof(int), node.size, fp);
            if(amount == 0)
            {
                printf("\nCouldn't read the file.");
                exit(0);
            }
        }
        else
        {
            printf("\nMalloc didn't succed.");
        }
    }
    e(MPI_Scatter(node.buffer, node.local_size, MPI_INT, node.local_buffer, node.local_size, MPI_INT, 0, comm));

    count_occurences(node);
    e(MPI_Reduce(node.local_occurences, node.occurences, node.intervals, MPI_LONG_INT, MPI_SUM, 0, comm));
    if(my_rank == 0)
    {
        display_histogram(node);
    }

    destroy(node);
    MPI_Finalize();
    return 0;
}
