void read_file(char* filename, int* min, int* max, long int* size, long int* local_size, int comm_sz, int my_rank, int bin_counts) {
  FILE* fp;
  *max = INT_MIN;
  *min = INT_MAX;
  struct stat file_stat;
  unsigned long long int amount;
  int* buffer = NULL;

  fp = fopen(filename, "r");
  if(fp == NULL) {
    printf("\n File doesn't exist.")
  }

  int result = stat(filename, &file_stat);
  if(result == -1) {
    printf("\nFile invalid.")
    exit(0);
  }
  *size = file_stat.st_size;
  *size /= sizeof(int);
  *local_size = *size / comm_sz;
  buffer = malloc(*size *sizeof(int));
}
