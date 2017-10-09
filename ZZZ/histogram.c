/* Possible solutions
	1. Change data types to have consistency
	2. Write function for filename as input
	3. Problem with pointers
*/

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include <limits.h>
#include <sys/stat.h>
#include <assert.h>

//Functions

/*long int get_intervals(const char *s);*/
int get_intervals(const char *s);

/*int* getdatacount(const char* s, unsigned long* data_count);*/
void getdatacount(const char* s, int* data_count);

/*void getmaxmin(int* buffer, unsigned long size, int* min, int* max);*/
void getmaxmin(int* buffer, int size, float* min, float* max);

/*void broadcastdata(long int* bin_count, int* min, int* max, unsigned long* data_count, unsigned long* local_data_count, int comm_sz, int rank, MPI_Comm comm);*/
void broadcastdata(int* bin_count, float* min, float* max, int* local_data_count, int* data_count, int comm_sz, int rank, MPI_Comm comm);

/*long int* setbins(int min, int max, long int bin_count, int loc_bin_cts[]);*/
void setbins(float bin_max[], float min, float max, int bin_count, int loc_bin_cts[], int rank, MPI_Comm comm);

/*int* findbins(int buffer[], long int bin_max[], long int bin_count, int loc_bin_cts[], unsigned long local_data_count, int min);*/
void findbins(int bin_counts[], int local_data[], float bin_max[], int bin_count, int loc_bin_cts[], int local_data_count, float min, MPI_Comm comm);

int whichbin(float data, float bin_max[], int bin_count, float min);


/*void print_histogram(long int* bin_max, int* bin_cts, long int bin_count, int* loc_bin_cts);*/
void print_histogram(float bin_max[], int bin_cts[], int bin_count, float min);


void e (int error);


int main(int argc, char* argv[])
{
	int bin_count;
	float min;
	float max;
	float* bin_max;
	int* bin_cts;
	int* loc_bin_cts;
	int data_count;
	int local_data_count;
	int* data;
	int* local_data;
	int rank;
	int comm_sz;
	char* filename;
	MPI_Comm comm;

	//Initialize MPI
	e(MPI_Init(&argc, &argv));
	comm = MPI_COMM_WORLD;
	e(MPI_Comm_size(comm, &comm_sz));
	e(MPI_Comm_rank(comm, &rank));
	filename =  malloc(sizeof(char) * 128);

/* Getting all inputs */
	if(rank == 0)
	{
		printf("What is the data file? :")
		scanf("%126s", filename);
		printf("\nWhat are the intervals? : ");
		scanf("%d", bin_count);
		//bin_count = get_intervals();
		getdatacount(filename, &data_count);
		local_data_count = data_count/comm_sz;
		data_count = comm_sz * local_data_count;
		//How much data you have //
		data = malloc(data_count * sizeof(int));
		//											//
		getmaxmin(data, data_count, &min, &max);
		broadcastdata(&bin_count, &min, &max, &data_count, &local_data_count, comm_sz, rank, comm);
	}
/* Getting all inputs */


	bin_max = malloc(bin_count*sizeof(long int));
	bin_cts = malloc(bin_count*sizeof(int));
	loc_bin_cts = malloc(bin_count*sizeof(int));
	local_data = malloc(local_data_count*sizeof(int));

/* Generating data seperatly */
	e(MPI_Scatter(data,local_data_count,MPI_FLOAT,local_data,local_data_count,MPI_FLOAT, 0, comm));
	if(rank  == 0)
		free(data);
/*													*/

/* Setting bins */
	setbins(bin_max, min, max, bin_count, loc_bin_cts, rank, comm);
/*							*/

	findbins(bin_cts, local_data, bin_max, bin_count, loc_bin_cts, local_data_count, min, comm);
	e(MPI_Reduce(loc_bin_cts, bin_cts, bin_count, MPI_INT, MPI_SUM, 0, comm));
	if(rank==0)
	print_histogram(bin_max,bin_cts,bin_count,min);

	free(local_data);
	free(bin_max);
	free(bin_cts);
	free(loc_bin_cts);
	//free(filename);
	MPI_Finalize();
	return 0;
}

/*long */int get_intervals(const char* s)
{
	if(s == NULL)
	{
		printf("Improper arguments.");
	}
	char* temp;
	long int num = strtol(s, &temp, 10);
	return num;
}

/*void getmaxmin(int buffer[], unsigned long size, int* min, int* max)*/
void getmaxmin(int buffer[], int size, float* min, float* max)
{
	*max = INT_MAX;
	*min = INT_MIN;
	for(unsigned long int i = 0; i< size; ++i)
	{
		if(buffer[i] < *min)
		{
			*min = buffer[i];
		}
		if(buffer[i] > *max)
		{
			*min = buffer[i];
		}
	}
}

void getdatacount(const char* s, int* data_count)
{
	FILE* fp = fopen(s, "r");
	struct stat file_stat;
	stat(s, &file_stat);
	*data_count = file_stat.st_size;
	*data_count = *data_count/sizeof(int);
	/*int* buffer = malloc(*data_count * sizeof(int));
	fread(buffer, sizeof(int), *size, fp);
	free(buffer);
	*/
}

/*void broadcastdata(long int* bin_count, int* min, int* max, unsigned long* local_data_count, unsigned long* data_count, int comm_sz, int rank, MPI_Comm comm)*/
void broadcastdata(int* bin_count, float* min, float* max, int* local_data_count, int* data_count, int comm_sz, int rank, MPI_Comm comm)
{
	e(MPI_Bcast(bin_count, 1, MPI_INT, 0, comm));
	e(MPI_Bcast(min, 1, MPI_FLOAT, 0, comm));
	e(MPI_Bcast(max, 1, MPI_FLOAT, 0, comm));
	e(MPI_Bcast(data_count, 1, MPI_INT, 0, comm));
	e(MPI_Bcast(local_data_count, 1, MPI_INT, 0, comm));
}

/*long int* setbins(int min, int max, long int bin_count, int loc_bin_cts[])*/
void setbins(float bin_max[], float min, float max, int bin_count, int loc_bin_cts[], int rank, MPI_Comm comm)
{
	float range = max - min;
	float interval =  range/bin_count;
	for(int i = 0; i<bin_count; i++)
	{
		bin_max[i] = interval * (float)(i+1);
		loc_bin_cts[i] = 0;
	}
}

/*int* findbins(int local_data[], long int bin_max[], long int bin_count, int loc_bin_cts[], unsigned long local_data_count, int min)*/
void findbins(int bin_counts[], int local_data[], float bin_max[], int bin_count, int loc_bin_cts[], int local_data_count, float min, MPI_Comm comm)
{
	int i=0;
	int bin;

	for(i =0; i<local_data_count; i++)
	{
		bin = whichbin(local_data[i], bin_max, bin_count, min);
		loc_bin_cts[i]++;
	}
}

/*int whichbin(int data, long int bin_max[], long int bin_count, int min)*/
int whichbin(float data, float bin_max[], int bin_count, float min)
{
	long int i;
	for(i = 0; i<bin_count-1; i++)
	{
		if(data <=bin_max[i]) break;
	}
	return i;
}

/*void print_histogram(long int bin_max[], int bin_cts[], long int bin_count, float min)*/
void print_histogram(float bin_max[], int bin_cts[], int bin_count, float min)

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
		printf("%10.3f |", bin_max[i]);
		row_width = (float) bin_cts[i] / (float) max * (float) width;
		for(j=0; j<row_width; j++)
		{
			printf("#");
		}
		printf("   %d\n", bin_cts[i]);
	}
}

void e(int error)
{
	if(error != MPI_SUCCESS)
	{
		fprintf(stderr, "Error starting MPI program. \n");
		MPI_Abort(MPI_COMM_WORLD, error);
		MPI_Finalize();
		exit(1);
	}
}
