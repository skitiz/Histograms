#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <time.h>

//Gonna leave error handling for now.

size_t data_count (const char* error_message);

void create_file (const char* file_name, size_t num);

int main(int argc, char** argv)
{
	if(argc < 3)
	{
		//Error handling stuff.
	}
	size_t num = data_count (argv[2]);
	create_file(argv[1], num);
}

size_t data_count (const char* arg)
{
	char* temp;
	size_t num = strtol (arg, &temp, 10);
	//Error handling for temp
	return num;
}

void create_file (const char* file_name, size_t num)
{
	FILE* ptr = fopen(file_name, "wb");
	if(ptr == NULL)
	{
		printf("Unable to create file.");
	}
	srand(time(NULL));
	for(size_t i=0; i<num; i++)
	{
	    int value = rand() %100;
	    fwrite(&value, sizeof(value), 1, ptr);
	}
	fclose(ptr);
}