#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "histogram.h"

void get_input( 
	int* bin_count;
	int rank;
	int comm_sz;
	MPI_Comm comm);

void setbins(	
);

long int get_intervals(const char* s);

int main(int argc, char** argv[])
{
	long int bin_count;
	float* bin_max;
	unsigned long data_count;
	unsigned long local_data_count;
	int rank;
	int* data;
	float* data;
	float* local_data
	float min;
	float max;
	int* bin_cts;
	int comm_sz;
	const char* filename;
	MPI_Comm comm;

	//Initialize MPI
	MPI_Init(&argc, &argv);
	comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm, &comm_sz);
	MPI_Comm_rank(comm, &rank);

	filename = argv[1];
	bin_count = getintervals(argv[2]);
	data_count = getdatacount(argv[1]);
	local_data_count = data_count/comm_sz;
	data_count = comm_sz * local_data_count;
	getmaxmin(argv[1], &data_count, &min, &max);
	broadcastdata(bin_count, min, max, data_count, comm_sz, rank, comm);
	bin_max = malloc(bin_count*sizeof(float));
	data = malloc(data_count*sizeof(int));


	FILE* fp = fopen(filename, "r");
	struct stat file_stat;
	stat(filename, &file_stat);
	unsigned long long int *size = file_stat.st_size;
	int* buffer = malloc(*size * sizeof(int));
	fread(buffer, sizeof(int), *size, fp);


	MPI_Scatter(buffer,local_data_count,MPI_UNSIGNED_LONG,local_data,local_data_count,MPI_UNSIGNED_LONG, 0, comm)
	bin_max = setbins(min, max, bin_count);
	bin_cts = findbins(buffer, bin_max, bin_count);
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

void getmaxmin(char* filename, unsigned long long int* size, int* min, int* max)
{
	*max = INT_MAX;
	*min = INT_MIN;
	FILE* fp = fopen(filename,"r");
	struct stat file_stat;
	*size = file_stat.st_size;
	*size = *size/sizeof(int);
	int* buffer = malloc(*size *sizeof(int));
	fread(buffer, sizeof(int), *size, f);
	for(unsigned long long int i = 0; i<*size; ++i)
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

unsigned long data_count(char* s)
{
	unsigned long* size;
	struct stat file_stat;
	stat(s, &file_stat);
	*size = file_stat.st_size;
	*size = *size/sizeof(int);
	return size;
}

void broadcastdata(unsigned int int bin_count, int min, int max, int comm_sz, int rank, MPI_Comm comm)
{
	MPI_Bcast(bin_count, 1, MPI_LONG_INT, 0, comm);
	MPI_Bcast(min, 1, MPI_INT, 0, comm);
	MPI_Bcast(max, 1, MPI_INT, 0, comm);
	MPI(data_count, 1, MPI_UNSIGNED_LONG, comm);
	//MPI(loc_data_count, 1, MPI_UNSIGNED_LONG, comm);
}

float* setbins(int min, int max, long int bin_count)
{
	long int range = max - min;
	long int interval =  range/bin_count;
	float* endpoints;
	for(int i = 0; i<bin_count, i++)
	{
		endpoints[i] = interval * (long int)(i+1);
		bin_cts[i] = 0;
	}
	return endpoints;
}

int* findbins(int buffer[], float bin_max[], long int bin_count, int min)
{
	int* bin_cts; 
	int bin;
	for(long int i =0; i<data_count; i++)
	{
		bin = whichbin(buffer[i], bin_max, bin_count, min);
		bin_cts[i]++;
	}
	return bin_cts;
}

int whichbin(int buffer, float bin_max[], long int bin_count)
{
	long int i;
	for(i = 0; i<bin_count; i++)
	{
		if(buffer <=bin_max[i]) break;
	}
	return i;
}
