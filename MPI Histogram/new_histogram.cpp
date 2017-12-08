#include <iostream.h>
#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include <assert.h>
#include "histfunc.h"

int get_args(char* s) {
    int temp;
    temp = strtol(s, NULL, 10);
    return temp;
}

read_file(char *filename, int &min, int &max, long int &size, long int &local_size, int comm_sz,
            int my_rank, int bin_counts);

void main(int argc, char *argv[]) {
    //
    char* filename;
    int min, max, comm_sz, my_rank, bin_counts;
    long int size, local_size;
    my_rank == 0;

    if(argc<3){
        cout<<"\tMissing arguments. Please reenter arguments.";
    }

    //Get user inputs.
    filename = argv[1];
    bin_counts = get_args(argv[2]);

    read_file(filename, &min, &max, &size, &local_size, comm_sz, my_rank, bin_counts);
    
}
