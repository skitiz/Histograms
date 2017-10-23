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
    unsigned long long int size;
    unsigned long long int local_size;
    int* buffer;
    int* local_buffer;
    float* endpoints;
    int* occurences;
    int* local_occurences;

} node;

void display_histogram (void *ptr)
{
    node *data = (node *) ptr;
    assert(data->occurences != NULL);
    assert(data->endpoints != NULL);

    int width = 40;
    int max = 0;
    int row_width;
    int i;
    int j;

    for(i = 0; i < data->intervals; i++)
    {
        if(data->occurences[i] > max)
            max = data->occurences[i];
    }

    for(i = 0; i < data->intervals; i++ )
    {
        cout<<" |" << data->endpoints[i];
        row_width = data->occurences[i] /  max * width;
        for(j=0; j< row_width; j++)
        {
            cout<< "#";
        }
        cout << "   \n" << data->occurences[i];
    }
}

size_t determine_index(int temp, float* endpoints, long int intervals)
{
    assert(endpoints != NULL);
    size_t index;
    for( index = 0; index < intervals -1 ; index++)
    {
        if(temp <= endpoints[index]) break;
    }
    return index;
}

void count_occurences(void *ptr)
{
    node *data = (node *) ptr;
    assert(data->buffer != NULL);
    assert(data->endpoints != NULL);

    //data->occurences = malloc (data->intervals, sizeof(int));
    if(data->occurences == NULL)
    {
        cout<<"\nMemory allocation failed....Exiting." ;
        MPI_Finalize();
        exit(0);
    }
    #pragma omp parallel for
    for(unsigned long long int i = 0; i < data->local_size; i++)
    {
        size_t index = determine_index (data->local_buffer[i], data->endpoints, data->intervals);
        node->local_occurences[index]++;
    }
}


void read_file(void *ptr)
{
    node *data = (node *) ptr;
    if(data->my_rank == 0)
    {
        ifstream fp;
        struct stat file_stat;
        vector <int> v;
        unsigned long long int amount;

        fp.open(data->filename, ios::binary);
        if(fp.is_open())
        {
            int result = stat(data->filename, &file_stat);
            if(result == -1)
            {
                cout << "\nFile Invalid.";
                MPI_Finalize();
                exit(0);
            }
            data->size = file_stat.st_size;
            data->size /= sizeof(int);
            data->local_size = data->size/ data->comm_sz;
            //data->buffer = malloc( data->size * sizeof(int));
            if(data->buffer)
            {
                fp.read(&v[0], data->size);
                amount = fp.gcount();
                if(amount == 0)
                {
                    cout << "\nCouldn't read the file.";
                    MPI_Finalize();
                    exit(0);
                }
            }
            else
            {
                cout << "\nValue of malloc didn't succed.";
                MPI_Finalize();
                exit(0);
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
            fp->fileLength = fp.tellg();
            fp->iterations = (int) ceil((double) data->fileLength / data->bufferSize);

            if (data->fileLength < data->bufferSize)
            {
                data->minMessage = (data->fileLength / data->unit) / data->comm_sz;
            }
        }
        else
        {
            cout << "\nCannot open file.";
        }
    }
    MPI_Bcast(data->intervals, 1, MPI_LONG_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(data->min, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(data->max, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(data->size, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(data->local_size, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
}

void determine_intervals(void *ptr)
{
    node *data = (node *) ptr;
    float length = (data->max-data->min) / (float) data->intervals;
    float temp = data->min;

    size_t i = 0;
    for(i =0; i < data->intervals; i++)
    {
        data->endpoints[i] = length;
        length += length;
        data->local_occurences[i] = 0;
    }
}

int main(int argc, char* argv[])
{
    clock_t begin = clock();

    node data;
    data.size = 10000;
    data.local_size = 10000;
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
    int numThreads = user_arguments(argv[2]);
    assert(numThreads > 0);

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &data.comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &data.my_rank);


    data.local_buffer = malloc(local_size * sizeof(int));
    data.endpoints = malloc( data.intervals * sizeof(float));
    data.buffer = malloc(data.size = * sizeof(int));
    data.local_buffer = malloc(local_size * sizeof(int));
    data.local_occurences = malloc( data.intervals * sizeof(int));
    data.occurences = malloc (data.intervals * sizeof(int));

    omp_set_dynamic(0);
    read_file(&data);

    MPI_Scatter(data->buffer, data->local_size, MPI_INT, data->local_buffer, data->local_size, MPI_INT, 0, MPI_COMM_WORLD));

    determine_intervals(&data);

    count_occurences(&data);
    MPI_Reduce(data->local_occurences, data->occurences, data->intervals, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD));

    if(node.my_rank == 0)
    {
        display_histogram(&data);
    }

    free(data->buffer);
    free(data->local_buffer);
    free(data->endpoints);
    free(data->occurences);

    MPI_Finalize();

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    cout<< "\nTime spent = " << time_spent;

    return 0;
}
