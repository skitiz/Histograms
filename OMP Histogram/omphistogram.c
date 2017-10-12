#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>
#include <sys/stat.h>
#include <limits.h>

long int get_intervals(char* s);
int* max_min(char* filename, unsigned long long int* size, int* min, int* max);
float* determine_intervals(int min, int max, long int intervals);

int main(int argc, char* argv[])
{
		long int intervals;
		int max, min;
		unsigned long long int size;
		int* buffer = NULL;
		float* endpoints = NULL;

		int thread_count = strtol(argv[1], NULL, 10);

		//File name
		char* s = argv[2];

		//Get the number of intervals
		intervals = get_intervals(argv[3]);
		printf("%ld", intervals);

		//Determine max-min
		buffer = max_min(s, &size, &max, &min);

		//Determine intervals
		endpoints = determine_intervals(min, max, intervals);

		for(long int i=0; i< intervals; i++)
		{
			printf("%6f , ", endpoints[i]);
		}

		printf("\nEnd of program.");
		return 0;
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
			printf("\nFile doesn't exist.")
		}
		int result = stat(filename, &file_stat);
		if(result == -1)
		{
			printf("\nFile invalid.")
		}
		*size = file_stat.st_size;
		*size /= sizeof(int);
		buffer = malloc(*size *sizeof(int));
		if(buffer)
		{
			amount = fread(buffer, sizeof(int), *size, fp);
			if(amount == NULL)
			{
				printf("\nCouldn't read.")
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
