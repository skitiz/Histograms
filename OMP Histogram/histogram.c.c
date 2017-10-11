#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>
#include <sys/stat.h>
#include <limits.h>

long int get_intervals(char* s);
int* max_min(char* filename, unsigned long long int* size, int* min, int* max);

int main(int argc, char* argv[]) 
{
	long int intervals;
	int max, min;
	unsigned long long int size;

	int thread_count = strtol(argv[1], NULL, 10);

	//File name
	char* s = argv[2];

	//Get the number of intervals
	intervals = get_intervals(argv[3]);
	printf("%ld", intervals);

	//Determine max-min
	int* buffer = max_min(s, &size, &max, &min);

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

	fp = fopen(filename, "r");
	stat(filename, &file_stat);
	*size = file_stat.st_size;
	*size /= sizeof(int);
	int* buffer = malloc(*size *sizeof(int));
	amount = fread(buffer, sizeof(int), *size, fp);
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
