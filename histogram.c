#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
//#include "histogram.h"
#include <string.h>
#include <limits.h>
#include <sys/stat.h>
#include <util.h>
#include <assert.h>

//Functions

long int get_intervals(const char *s);

int* getdatacount(const char* s);

void getmaxmin(char* filename, unsigned long long int* size, int* min, int* max);

void broadcastdata(unsigned int int bin_count, int min, int max, int comm_sz, int rank, MPI_Comm comm);

long int* setbins(int min, int max, long int bin_count);

int* findbins(int buffer[], long int bin_max[], long int bin_count, int min);

int whichbin(int buffer, long int bin_max[], long int bin_count);

void print_histogram(long int bin_max, int* bin_cts, long int bin_count, int* loc_bin_cts);


int main(int argc, char** argv[])
{
	long int bin_count;
	long int* bin_max;
	unsigned long data_count;
	unsigned long local_data_count;
	int rank;
	int* buffer;
	int* local_data;
	float min, max;
	int* bin_cts, loc_bin_cts;
	int comm_sz;
	const char* filename;
	MPI_Comm comm;

	//Initialize MPI
	MPI_Init(&argc, &argv);
	comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm, &comm_sz);
	MPI_Comm_rank(comm, &rank);


	filename = argv[1];
	bin_count = get_intervals(argv[2]);
	buffer = getdatacount(filename, &data_count);
	local_data_count = data_count/comm_sz;
	data_count = comm_sz * local_data_count;
	getmaxmin(buffer, data_count, &min, &max);
	broadcastdata(bin_count, min, max, data_count, local_data,count, comm_sz, rank, comm);
	

	bin_max = malloc(bin_count*sizeof(long int));
	bin_cts = malloc(bin_count*sizeof(int));
	loc_bin_cts = malloc(bin_count*sizeof(int));
	local_data = malloc(local_data_count*sizeof(int));


	MPI_Scatter(buffer,local_data_count,MPI_UNSIGNED_LONG,local_data,local_data_count,MPI_UNSIGNED_LONG, 0, comm)
	
	free(buffer);

	bin_max = setbins(min, max, bin_count,loc_bin_cts);
	bin_cts = findbins(local_data, bin_max, bin_count, loc_bin_cts);
	MPI_Reduce(loc_bin_cts, bin_cts, bin_count, MPI_LONG_INT, MPI_SUM, 0, comm);
	if(rank==0)
	print_histogram(bin_max,bin_cts,bin_count,min);

	free(local_data);
	free(bin_max);
	free(bin_cts);
	free(loc_bin_cts);
	free(filename);
	MPI_Finalize();
	return 0;
}

long int get_intervals(const char* s)
{
	if(s == NULL)
	{
		printf("Improper arguments.")
	}
	char* temp;
	long int num = strtol(s, &temp, 10);
	return num;
}

void getmaxmin(int buffer[], unsigned int size, int* min, int* max)
{
	*max = INT_MAX;
	*min = INT_MIN;
	for(unsigned long long int i = 0; i< size; ++i)
	{
		if(buffer[i] < *min)
		{
			*min = buffer[i];
		}
		if(buffer[i] > max)
		{
			*min = buffer[i];
		}
	}
}

int* getdatacount(const char* s, unsigned long* size)
{
	FILE* fp = fopen(filename, "r");
	struct stat file_stat;
	stat(s, &file_stat);
	*size = file_stat.st_size;
	*size = *size/sizeof(int);
	int* buffer = malloc(*size * sizeof(int));
	fread(buffer, sizeof(int), *size, fp);
	return buffer;
}

void broadcastdata(long int bin_count, int min, int max, unsigned long local_data_count, int comm_sz, int rank, MPI_Comm comm)
{
	MPI_Bcast(bin_count, 1, MPI_LONG_INT, 0, comm);
	MPI_Bcast(min, 1, MPI_INT, 0, comm);
	MPI_Bcast(max, 1, MPI_INT, 0, comm);
	MPI(data_count, 1, MPI_UNSIGNED_LONG, 0, comm);
	MPI(loc_data_count, 1, MPI_UNSIGNED_LONG, 0, comm);
}

long int* setbins(int min, int max, long int bin_count, int loc_bin_cts[])
{
	long int range = max - min;
	long int interval =  range/bin_count;
	float* endpoints;
	for(int i = 0; i<bin_count, i++)
	{
		endpoints[i] = interval * (long int)(i+1);
		loc_bin_cts[i] = 0;
	}
	return endpoints;
}

int* findbins(int local_data[], long int bin_max[], long int bin_count, int min, int loc_bin_cts[])
{
	int* bin_cts; 
	int bin;
	for(long int i =0; i<data_count; i++)
	{
		bin = whichbin(local_data[i], bin_max, bin_count, min);
		loc_bin_cts[i]++;
	}
	return bin_cts;
}

int whichbin(int data, long int bin_max[], long int bin_count)
{
	long int i;
	for(i = 0; i<bin_count; i++)
	{
		if(data <=bin_max[i]) break;
	}
	return i;
}

void print_histogram(long int bin_max[], int bin_cts[], long int bin_count, int min)
{
	int width=40;
	int max = 0;
	int row_width;
	long int i;
	int j;

	for(i=0; i<bin_count; i++)
	{
		if(bin_cts[i] > max)
		{
			max = bin_cts[i];
		}
	}
	for(i=0; i<bin_count; i++)
	{
		printf("%ld |", bin_max[i]);
		row_width = (float) bin_cts[i] / (float) max * (float) width;
		for(j=0; j<row_width; j++)
		{
			printf("#");
		}
		printf("   %d\n", bin_cts[i]);
	}
}
}
