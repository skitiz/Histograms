Histograms
==========
Calculating the run time of various ways of generating a histogram using different parallel and distributed computing algorithms. Using the serial version as a base for comparision.


_Note : I haven't uploaded the serial version yet. And there are a few seg faults in MPI for now. I'll redo the whole MPI thing._

1. MPI : Using a distributed memory system to distribute data over nodes and simultaneously reduce and generate a histogram.
2. OMP : Using a shared memory system to parallelize blocks of code to optimize performance.
3. OMP + MPI : Using parallelized `for` loops over distributed memory to generate a histogram.

# How do I run these?
Clone the repo and generate a data file for the program to read.

`gcc -g -Wall -o create_file create_file.c -std=c99`    (I'm using for loops so I need c++ implementation.)

./create_file 'name of data file' 'number of elements'


To compile the exec files :

`gcc -g -Wall -o -mpi -fopenmp execfile mpi_file.c`   (Pick and choose the optional depending on the program you're running.)

All the programs take similar command line arguments. 
+ Serial version: ./nameofexec 'data file to read from' '# of intervals of histogram'
+ MPI version: ./nameofexec '# of nodes of MPI' 'data_file_to_read_from' '# of intervals'
+ OMP version: ./nameofexec '# of threads' 'data_file_to_read_from' '# of intervals'

# Error checks
I've added error checks in the programs to ensure that the minimum # of command line arguments are provided by the user.

~Also painstakingly made sure that there are no segfaults.~

Please ensure you generate a data file __BEFORE__ you run the exec file.

# Updates
I've tried to add verbose comments about functions wherever I can, if you feel I'm lacking submit a pull request and I'll take a look.
