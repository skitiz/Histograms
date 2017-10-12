#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>
#include <sys/stat.h>
#include <limits.h>
#include <assert.h>

long int get_intervals(char* s);
int* max_min(char* filename, unsigned long long int* size, int* min, int* max);
float* determine_intervals(int min, int max, long int intervals);
unsigned long long int* occurences(int* buffer, unsigned long long int size, long int intervals, float* endpoints);
size_t determine_index(int temp, float* endpoints, long int intervals);

int main(int argc, char* argv[])
{
		long int intervals;
		int max, min;
		unsigned long long int size;
		int* buffer = NULL;
		float* endpoints = NULL;
		unsigned long long int* occurences = NULL;

		int thread_count = strtol(argv[1], NULL, 10);

		//File name
		char* s = argv[2];
		if(s == NULL)
		{
			printf("Invalid file name.");
			exit(0);
		}

		//Get the number of intervals
		intervals = get_intervals(argv[3]);
		//printf("%ld", intervals);

		//Determine max-min
		buffer = max_min(s, &size, &max, &min);

		//Determine intervals
		endpoints = determine_intervals(min, max, intervals);

		//Determine occurances of each number
		occurences = count_occurrences(buffer, size, intervals, endpoints);

		printf("\nEnd of program.");
		return 0;
}

unsigned long long int* count_occurrences(int* buffer, unsigned long long int size, long int intervals, float* endpoints)
{
	// Error checks
	assert (buffer != NULL);
	assert (endpoints != NULL);
	//

	unsigned long long int* occurences = calloc (intervals, sizeof(unsigned long long int));
	if( occurences == NULL)
	{
		printf("Memory allocation failed. Exiting...");
		exit(0);
	}
	#pragma omp parallel for
	for(unsigned long int i = 0; i< size; i++)
	{
		size_t index = determine_index (buffer[i], endpoints, intervals);
		num_occurrences[index]++;
	}
	return num_occurrences;
}

size_t determine_index(int temp, float* endpoints, long int intervals)
{
	assert (endpoints != NULL);
	size_t index = 0;
	while (index < intervals && endpoints[index] > temp)
		++index;
	return index;
}

long int get_intervals(char* s)
{
		char* temp;
		long int num = strtol(s, &temp, 10);
		return num;
}

int* max_min(char* filename, unsigned long long int* size, int *max, int *min)
{
		FILE *fp;
		*max = INT_MIN;
		*min = INT_MAX;
		struct stat file_stat;
		unsigned long long int amount;
		int* buffer = NULL;

		fp = fopen(filename, "r");
		if(fp == NULL)
		{
			printf("\nFile doesn't exist.");
			exit(0);
		}
		int result = stat(filename, &file_stat);
		if(result == -1)
		{
			printf("\nFile invalid.");
			exit(0);
		}
		*size = file_stat.st_size;
		*size /= sizeof(int);
		buffer = malloc(*size *sizeof(int));
		if(buffer)
		{
			amount = fread(buffer, sizeof(int), *size, fp);
			if(amount == 0)
			{
				printf("\nCouldn't read.");
				exit(0);
			}
		}
		else
			printf("\nValue of malloc didn't succed.");
			for(unsigned long long int i = 0; i < *size; i++)
			{
				if((buffer[i]) < *min)
				{
					*min = buffer[i];
				}
				if((buffer[i]) > *max)
				{
					*max = buffer[i];
				}
			}
		return buffer;
}

float* determine_intervals(int min, int max, long int intervals)
{
	float* endpoints = malloc(intervals * sizeof(float));
	float length = (max-min) / (float) intervals;
	float temp = min;
	//#pragma omp parallel for
	for(size_t i = 0; i < intervals; i++)
	{
		endpoints[i] = temp;
		temp += length;
	}
	return endpoints;
}
