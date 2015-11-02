//calculates the estimated arithmetic time
// multiplies a large array with a constant

#include<stdio.h>
#include<stdlib.h>
#include "mpi.h"

void main(int argc, char** argv)
{
  int* huge_array;
  int constant = 10;

  double start_time,end_time;
  MPI_Init(NULL,NULL);

 start_time = MPI_Wtime();
  huge_array=malloc(sizeof(int)*1024);
  for(int i =0;i<1024;i++)
    {
      huge_array[i] = huge_array[i]*constant;
    }
 end_time = MPI_Wtime();
 printf("time taken %f\n",end_time-start_time);

}
