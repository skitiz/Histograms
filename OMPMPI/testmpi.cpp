/*
        Kshitij Bantupalli
        Parallel and Distributed Systems
        Assignment 4

        The program is written in C++.

        To compile : mpiCC -g -filename -fopenmp -o mpi testmpi.cpp
        To run : mpiexec -n <n> ./mpi <intervals> <filename> <threads>


        Build status : Working

*/

#include <sstream>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <assert.h>
#include <mpi.h>
#ifdef _OPENMP
    #include <omp.h>
#endif
#include <sys/stat.h>


using namespace std;


//          Define the struct for MPI       //
typedef struct
{
    //      MPI Variables       //
    int my_rank;
    int comm_sz;

    //      Root Variables      //
    int intervals;
    char* filename;
    int min;
    int max;
    long int size;
    long int local_size;
    int* buffer;
    int* local_buffer;
    float* endpoints;
    int* occurences;
    int* local_occurences;

} node;

//          Broadcasting data to all nodes of MPI       //
void build_mpi_data_type(int* data_1, int* data_2, int* data_3, long int* data_4, long int* data_5)
{
    MPI_Bcast(data_1, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(data_2, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(data_3, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(data_4, 1, MPI_LONG_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(data_5, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
}


//          Converting user arguments for char* to usable values    //
int user_arguments(char *argv) {
    char *endptr;
    int intervalSize = strtol(argv, &endptr, 10);
    if (!*argv || *endptr)
        cerr << "Invalid number " << argv << '\n';
    return intervalSize;
}


//          Root node will display histogram        //
void display_histogram (void *ptr)
{
    node *data = (node *) ptr;
    assert(data->occurences != NULL);
    assert(data->endpoints != NULL);

    int i;

    float length = data->endpoints[1] - data->endpoints[0];

    for(i = 0; i < data->intervals; i++ )
    {
        cout<< data->endpoints[i] << " - " << data->endpoints[i] + length;
        cout << "    " << data->occurences[i] << "\n";
    }
}


//              Determine the index of the elements             //
size_t determine_index(int temp, float* endpoints, long int intervals)
{
    assert(endpoints != NULL);
    size_t index;
    for( index = 0; index < intervals -1 ; index++)
    {
        if(temp <= endpoints[index])
        {
            break;
        }
    }
    return index;
}


//     Counting occurances for histogram. Currently sorts everything into one interval.     //
void count_occurences(void *ptr, int numThreads)
{
    node *data = (node *) ptr;
    assert(data->buffer != NULL);
    assert(data->endpoints != NULL);
    //Edits

    //data->occurences = malloc (data->intervals, sizeof(int));
    if(data->occurences == NULL)
    {
        cout<<"\nMemory allocation failed....Exiting." ;
        MPI_Finalize();
        exit(0);
    }
    #pragma omp parallel for
    for(long int i = 0; i < data->local_size; i++)
    {
        size_t index = determine_index (data->local_buffer[i], data->endpoints, data->intervals);
        data->local_occurences[index]++;
    }
}


//          Reading the file.       //
void read_file(void *ptr)
{
    node *data = (node *) ptr;
    if(data->my_rank == 0)
    {
        FILE* fp;
        struct stat file_stat;
		unsigned long long int amount;

		fp = fopen(data->filename, "r");
		if(fp == NULL)
		{
			printf("\nFile doesn't exist.");
			exit(0);
		}
		int result = stat(data->filename, &file_stat);
		if(result == -1)
		{
			printf("\nFile invalid.");
			exit(0);
		}
		data->size = file_stat.st_size;
		data->size /= sizeof(int);
		if(data->buffer)
		{
			amount = fread(data->buffer, sizeof(int), data->size, fp);
			if(amount == 0)
			{
				printf("\nCouldn't read.");
				exit(0);
			}
		}
		else
		{
			printf("\nValue of malloc didn't succed.");
		}

        for(unsigned long long int i = 0; i < data->size; i++)
        {
            if(data->buffer[i] < data->min)
            {
                data->min = data->buffer[i];
            }
            if(data->buffer[i] > data->max)
            {
                data->max = data->buffer[i];
            }
        }
            //fp->fileLength = fp.tellg();
            //fp->iterations = (int) ceil((double) data->fileLength / data->bufferSize);

            //if (data->fileLength < data->bufferSize)
            //{
            //    data->minMessage = (data->fileLength / data->unit) / data->comm_sz;
            //}
    }
}
/*
    MPI_Bcast(data->intervals, 1, MPI_LONG_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(data->min, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(data->max, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(data->size, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(data->local_size, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
*/


//             Determine intervals from # of intervals and max, min         //
void determine_intervals(void *ptr)
{
    node *data = (node *) ptr;
    //float length = (data->max-data->min) / (float) data->intervals;
    //float temp = data->min;

    size_t i = 0;
    for(i =0; i < data->intervals; i++)
    {
        data->endpoints[i] = data->intervals * (i + 1);
        data->local_occurences[i] = 0;
    }
}



int main(int argc, char* argv[])
{
    clock_t begin = clock();

    node data;
    data.size = 1000000;
    data.local_size = 1000000;
    data.buffer = NULL;
    data.local_buffer = NULL;
    data.endpoints = NULL;
    data.occurences = NULL;
    data.local_occurences = NULL;

    //File for reading
    data.filename = argv[2];

    //Get intervals for Histogram
    data.intervals = user_arguments(argv[1]);
    assert(data.intervals > 0);

    //Get threads for OMP
    int numThreads = user_arguments(argv[3]);
    assert(numThreads > 0);

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &data.comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &data.my_rank);


    data.local_buffer = (int *) malloc(sizeof(int) * data.local_size);
    data.endpoints = (float *) malloc( data.intervals * sizeof(float));
    data.buffer = (int *) malloc(data.size * sizeof(int));
    data.local_buffer = (int *) malloc(data.local_size * sizeof(int));
    data.local_occurences = (int *) malloc( data.intervals * sizeof(int));
    data.occurences = (int *) malloc (data.intervals * sizeof(int));

    read_file(&data);

    build_mpi_data_type(&data.intervals, &data.min, &data.max, &data.size, &data.local_size);
    MPI_Scatter(data.buffer, data.local_size, MPI_INT, data.local_buffer, data.local_size, MPI_INT, 0, MPI_COMM_WORLD);

    determine_intervals(&data);

    count_occurences(&data, numThreads);
    MPI_Reduce(data.local_occurences, data.occurences, data.intervals, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if(data.my_rank == 0)
    {
        display_histogram(&data);
    }

    free(data.buffer);
    free(data.local_buffer);
    free(data.endpoints);
    free(data.occurences);

    MPI_Finalize();

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    cout<< "\nTime spent = " << time_spent;

    cout<<"\n";

    return 0;
}
