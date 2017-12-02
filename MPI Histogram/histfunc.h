//
// Created by Kshitij Bantupalli on 11/27/2017.
//

#ifndef HISTOGRAMS_HISTFUNC_H
#define HISTOGRAMS_HISTFUNC_H

#endif //HISTOGRAMS_HISTFUNC_H

#include <stdio.h>
#include <stdlib.h>


typedef struct data
{
    long int size, local_size;
    int comm_sz, my_rank, bin_counts;
    long fsize;
};

read_file(char *filename, int *min, int *max, long int *size, long int *local_size,
int comm_sz, int my_rank, int bin_counts) {
    int i = 0;
    if(my_rank == 0){
        FILE *f = fopen(filename, "r");
        fseek(f, 0, SEEK_END);
        *fsize = ftell(f);
        fseek(f, 0, SEEK_SET);
        char* buffer = malloc(fsize + 1);
        size = buffer/sizeof(int);
        local_size = *size/comm_sz;
        fread(buffer, fsize, 1, f);
        while( i < *size) {
            if(buffer[i] < *min){
                *min = buffer[i];
            }
            else if(buffer[i] >= *max) {
                *max = buffer[i]
            }
            i++;
        }
    }
}