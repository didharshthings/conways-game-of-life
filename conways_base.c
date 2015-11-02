// Conway's Game of Life
// Main Executable Program
//
// CSCI 4576/5576 High Performance Scientific Computing
// Michael Oberg, modified from code supplied by Dr. Matthew Woitaszek

// siddharth singh, added logic for the assignment
// to run: mpirun -n 9 ./conway_base <file_name> <distribution> <iteration> <count_bugs> <count_at> <write_from> <write_to>
// mpirun -n 9 ./conway_base life.pgm 1 100 1 10 (0 for slice, 1 for checkerboard)


// System includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"

//check for np=4,9,25,36

// Include global variables. Only this file needs the #define
#define __MAIN
#include "globals.h"
#undef __MAIN

// User includes
#include "pprintf.h"
#include "pgm.h"
#include <math.h>
#include "papi.h"



enum direction //used for checkerboard distribution
{
  UP,DOWN,LEFT,RIGHT,
  UPLEFT,UPRIGHT,DOWNLEFT,DOWNRIGHT
};


void calc_neighbours(MPI_Comm comm_cart, int rank, int dimx, int dimy, int source_ranks[8], int dest_ranks[8])
{

  int my_coords[2];
  int n_coords[2];
  int temp_ranks;

  MPI_Cart_coords(comm_cart, rank, 2 ,my_coords);

  //Up
  MPI_Cart_shift(comm_cart,1,1,&source_ranks[UP],&dest_ranks[DOWN]);
  //Down
  MPI_Cart_shift(comm_cart,1,-1,&source_ranks[DOWN],&dest_ranks[UP]);
  //Right
  MPI_Cart_shift(comm_cart,0,-1,&source_ranks[LEFT],&dest_ranks[RIGHT]);
  //Left
  MPI_Cart_shift(comm_cart,0,1,&source_ranks[RIGHT], &dest_ranks[LEFT]);
  //upeft
  n_coords[1] = (my_coords[1] + 1) % dimx;
  n_coords[0] = (my_coords[0] - 1 + dimy) % dimy;
  MPI_Cart_rank(comm_cart, n_coords, &dest_ranks[UPLEFT]);
  source_ranks[DOWNRIGHT] = dest_ranks[UPLEFT];
  // upright
  n_coords[1] = (my_coords[1] + 1) % dimx;
  n_coords[0] = (my_coords[0] + 1) % dimy;
  MPI_Cart_rank(comm_cart, n_coords, &dest_ranks[UPRIGHT]);
  source_ranks[DOWNLEFT] = dest_ranks[UPRIGHT];
  // downleft
  n_coords[1] = (my_coords[1] -1 + dimx) % dimx;
  n_coords[0] = (my_coords[0] - 1 + dimy) % dimy;
  MPI_Cart_rank(comm_cart, n_coords, &dest_ranks[DOWNLEFT]);
  source_ranks[UPRIGHT] = dest_ranks[DOWNLEFT];
  // downright
  n_coords[1] = (my_coords[1] -1 + dimx) % dimx;
  n_coords[0] = (my_coords[0] + 1) % dimy;
  MPI_Cart_rank(comm_cart, n_coords, &dest_ranks[DOWNRIGHT]);
  source_ranks[UPLEFT] = dest_ranks[DOWNRIGHT];

  //int my_coords[2];
  MPI_Cart_coords(comm_cart, rank, 2 ,my_coords);

  //boundary cases
  /*
  0,0 1,0  2,0
  0,0| 0 | 3 | 6 |
  0,1| 1 | 4 | 7 |
  0,2| 2 | 5 | 8 |
  */

  //print final neighbours

  //if (rank == 1)
  //{
  pprintf("x - [%i], y-[%i] \n", my_coords[0], my_coords[1]);

  if(my_coords[0] == 0 && my_coords[1] == 0)
  {
    //pprintf("dont send left \n");
      dest_ranks[RIGHT] = MPI_PROC_NULL;
      dest_ranks[DOWNLEFT] = MPI_PROC_NULL;
      dest_ranks[UPLEFT] = MPI_PROC_NULL;
      source_ranks[RIGHT] = MPI_PROC_NULL;
      source_ranks[UPRIGHT] = MPI_PROC_NULL;
      source_ranks[DOWNRIGHT] = MPI_PROC_NULL;

      //pprintf("dont send up\n");
      dest_ranks[UP] = MPI_PROC_NULL;
      dest_ranks[DOWNRIGHT] = MPI_PROC_NULL;
      dest_ranks[DOWNLEFT] = MPI_PROC_NULL;
      source_ranks[UP] = MPI_PROC_NULL;
      source_ranks[UPLEFT] = MPI_PROC_NULL;
      source_ranks[UPRIGHT] = MPI_PROC_NULL;
    }


  if(my_coords[0] == 0 && my_coords[1] == nrows-1)
  {
    //pprintf("dont send left \n");
    dest_ranks[RIGHT] = MPI_PROC_NULL;
    dest_ranks[DOWNLEFT] = MPI_PROC_NULL;
    dest_ranks[UPLEFT] = MPI_PROC_NULL;
    source_ranks[RIGHT] = MPI_PROC_NULL;
    source_ranks[UPRIGHT] = MPI_PROC_NULL;
    source_ranks[DOWNRIGHT] = MPI_PROC_NULL;


    dest_ranks[DOWN] = MPI_PROC_NULL;
    dest_ranks[UPLEFT] = MPI_PROC_NULL;
    dest_ranks[UPRIGHT] = MPI_PROC_NULL;
    source_ranks[DOWN] = MPI_PROC_NULL;
    source_ranks[DOWNRIGHT] = MPI_PROC_NULL;
    source_ranks[DOWNLEFT] = MPI_PROC_NULL;

    //pprintf("dont send down\n");

  }


  if(my_coords[0] == nrows-1  && my_coords[1] == 0)
  {
    dest_ranks[LEFT] = MPI_PROC_NULL;
    dest_ranks[UPRIGHT] = MPI_PROC_NULL;
    dest_ranks[DOWNRIGHT] = MPI_PROC_NULL;
    source_ranks[LEFT] = MPI_PROC_NULL;
    source_ranks[DOWNLEFT] = MPI_PROC_NULL;
    source_ranks[UPLEFT] = MPI_PROC_NULL;
    //pprintf("dont send right\n");

    //pprintf("dont send up\n");
    dest_ranks[UP] = MPI_PROC_NULL;
    dest_ranks[DOWNRIGHT] = MPI_PROC_NULL;
    dest_ranks[DOWNLEFT] = MPI_PROC_NULL;
    source_ranks[UP] = MPI_PROC_NULL;
    source_ranks[UPLEFT] = MPI_PROC_NULL;
    source_ranks[UPRIGHT] = MPI_PROC_NULL;

  }

  if(my_coords[0] == nrows-1 && my_coords[1] == nrows-1)
  {

    dest_ranks[LEFT] = MPI_PROC_NULL;
    dest_ranks[DOWNRIGHT] = MPI_PROC_NULL;
    dest_ranks[UPRIGHT] = MPI_PROC_NULL;
    source_ranks[LEFT] = MPI_PROC_NULL;
    source_ranks[DOWNLEFT] = MPI_PROC_NULL;
    source_ranks[UPLEFT] = MPI_PROC_NULL;
    //pprintf("dont send right\n");


    dest_ranks[DOWN] = MPI_PROC_NULL;
    dest_ranks[UPLEFT] = MPI_PROC_NULL;
    dest_ranks[UPRIGHT] = MPI_PROC_NULL;
    source_ranks[DOWN] = MPI_PROC_NULL;
    source_ranks[DOWNRIGHT] = MPI_PROC_NULL;
    source_ranks[DOWNLEFT] = MPI_PROC_NULL;
    //pprintf("dont send down\n");

  }


  //columns
  if(my_coords[0] == 0)
  {
    //pprintf("dont send left \n");
    dest_ranks[RIGHT] = MPI_PROC_NULL;
    dest_ranks[DOWNLEFT] = MPI_PROC_NULL;
    dest_ranks[UPLEFT] = MPI_PROC_NULL;
    source_ranks[RIGHT] = MPI_PROC_NULL;
    source_ranks[UPRIGHT] = MPI_PROC_NULL;
    source_ranks[DOWNRIGHT] = MPI_PROC_NULL;
  }

  if(my_coords[0] == nrows-1)
  {
    dest_ranks[LEFT] = MPI_PROC_NULL;
    dest_ranks[DOWNRIGHT] = MPI_PROC_NULL;
    dest_ranks[UPRIGHT] = MPI_PROC_NULL;
    source_ranks[LEFT] = MPI_PROC_NULL;
    source_ranks[DOWNLEFT] = MPI_PROC_NULL;
    source_ranks[UPLEFT] = MPI_PROC_NULL;
    //pprintf("dont send right\n");

  }

  if(my_coords[1] == 0)
  {
    //pprintf("dont send up\n");
    dest_ranks[UP] = MPI_PROC_NULL;
    dest_ranks[DOWNRIGHT] = MPI_PROC_NULL;
    dest_ranks[DOWNLEFT] = MPI_PROC_NULL;
    source_ranks[UP] = MPI_PROC_NULL;
    source_ranks[UPLEFT] = MPI_PROC_NULL;
    source_ranks[UPRIGHT] = MPI_PROC_NULL;

  }

  if(my_coords[1] == nrows-1)
  {
      dest_ranks[DOWN] = MPI_PROC_NULL;
      dest_ranks[UPLEFT] = MPI_PROC_NULL;
      dest_ranks[UPRIGHT] = MPI_PROC_NULL;
      source_ranks[DOWN] = MPI_PROC_NULL;
      source_ranks[DOWNRIGHT] = MPI_PROC_NULL;
      source_ranks[DOWNLEFT] = MPI_PROC_NULL;

      //pprintf("dont send down\n");
    }

/*
    //print final neighbours
    if (rank==4){
    pprintf("rank %i sends upright rank %i - recieves from rank %i \n",rank,dest_ranks[DOWNRIGHT],source_ranks[UPLEFT]);
    pprintf("rank %i sends upleft rank %i - recieves from rank %i \n",rank,dest_ranks[DOWNLEFT],source_ranks[UPRIGHT]);
    pprintf("rank %i sends downright rank %i - recieves from rank %i \n",rank,dest_ranks[UPRIGHT],source_ranks[DOWNLEFT]);
    pprintf("rank %i sends downleft rank %i - recieves from rank %i \n",rank,dest_ranks[UPLEFT],source_ranks[DOWNRIGHT]);

    pprintf("rank %i sends down rank %i -recieves from rank %i \n",rank,dest_ranks[DOWN],source_ranks[DOWN]);
    pprintf("rank %i sends up rank %i - recieves from rank %i \n",rank,dest_ranks[UP],source_ranks[UP]);
    pprintf("rank %i sends left rank %i - recieves from rank %i \n",rank,dest_ranks[RIGHT],source_ranks[RIGHT]);
    pprintf("rank %i sends right rank %i - recieves from rank %i\n",rank,dest_ranks[LEFT],source_ranks[LEFT]);
  }
*/
  }//calculate neigbours

void update(int iteration) {

  int sum;
  int *pointer_old = (iteration%2==0)?field_a:field_b;
  int *pointer_new = (iteration%2==0)?field_b:field_a;

  for (int y = 0; y < local_height; y++)
  {
    for (int x = 0; x < local_width; x++)
    {
      int yb = y+1;
      int xb = x+1;

      sum  = 0;
      sum = pointer_old[(yb-1)*field_width+xb+1] + pointer_old[(yb-1)*field_width+xb] + pointer_old[(yb-1)*field_width+xb-1] + pointer_old[(yb)*field_width+xb+1]  +pointer_old[(yb)*field_width+xb-1] +pointer_old[(yb+1)*field_width+xb+1] + pointer_old[(yb+1)*field_width+xb] +pointer_old[(yb+1)*field_width+xb-1];

      //copy previous value
      pointer_new[yb*field_width+xb] = pointer_old[yb*field_width+xb];
      if(pointer_old[yb*field_width+xb]) //alive cell
      {
        if(sum == 2 || sum == 3)
        {
          pointer_new[yb*field_width+xb] = 1;
        }
        if(sum < 2 || sum > 3)
        {
          pointer_new[yb*field_width+xb] = 0;
        }
      }
      else //dead cell
      {
        if(sum == 3)
        {
          pointer_new[yb*field_width+xb] = 1;
        }
      }

    }
  }
}
void write_to_file (char *in_file, int iteration, int offset)
{
  char filename[1000];
  offset = 0;
  int a;
  MPI_Aint extent;
  MPI_Datatype filetype, contig;
  MPI_Offset disp;
  MPI_File fh;

  sprintf(filename, "%d_test_", iteration);
  strcat(filename, in_file);
  int header_offset;
  //write header
  if (rank == 0)
  {
    FILE *fp = fopen( filename, "w+" );
    header_offset = fprintf( fp, "%s\n%i %i\n%i\n", header, width, height, depth );
    fclose(fp);
    offset = header_offset;
  }

  MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDWR | MPI_MODE_APPEND, MPI_INFO_NULL, &fh);

  int *pointer = (iteration%2==0) ? (field_a) : (field_b);

  char *write_buffer=(char *)malloc( local_width * local_height * sizeof(char));

  //fill temp matrix
  for (int y = 1; y < local_height + 1; y++ )
  {
    for(int x = 1; x < local_width + 1; x++ )
    {
      if (pointer[ y * field_width + x])
      {
        a = 0;
      }
      else
      {
        a = 0xFF;
      }
      write_buffer[x-1+(y-1)*local_width] = (char)a;
    }
  }

  disp = offset;
  disp += my_col * width * local_height + (my_row) * local_width;

  MPI_Type_contiguous(local_width, MPI_CHAR, &contig);
  extent = width * sizeof(char);
  MPI_Type_create_resized(contig, 0, extent, &filetype);
  MPI_Type_commit(&filetype);

  // writes to the file
  MPI_File_set_view(fh, disp, MPI_CHAR, filetype, "native", MPI_INFO_NULL);
  MPI_File_write_all(fh, write_buffer, local_width*local_height, MPI_CHAR, MPI_STATUS_IGNORE);
  MPI_File_close(&fh);

  /*
  int gsizes[2],distribs[2],psizes[2],dargs[2];
  gsizes[0] = local_height;
  gsizes[1] = local_width;
  distribs[0] = MPI_DISTRIBUTE_BLOCK;
  distribs[1] = MPI_DISTRIBUTE_BLOCK;
  dargs[0] = MPI_DISTRIBUTE_DFLT_DARG;
  dargs[1] = MPI_DISTRIBUTE_DFLT_DARG;
  psizes[0] = sqrt(np);
  psizes[1] = sqrt(np);

  MPI_Type_create_darray(np,rank,2,gsizes,distribs,dargs,psizes,MPI_ORDER_C,MPI_INT,&filetype);
  MPI_Type_commit(&filetype);
  MPI_File_set_view(fh, 0, MPI_CHAR, filetype, "native", MPI_INFO_NULL);
  int local_array_size = local_width*local_height;
  MPI_File_write_all(fh, write_buffer, local_array_size, MPI_CHAR, MPI_STATUS_IGNORE);
  MPI_File_close(&fh);
  */

}


void sync_checkerboard(int iteration, int dims[2], int source_rank[8], int dest_rank[8],MPI_Comm comm_cart)
{
  MPI_Datatype coltype;
  MPI_Type_vector(local_width, 1, field_height, MPI_INT, &coltype);
  MPI_Type_commit(&coltype);

  int *pointer = (iteration%2==0)?field_a:field_b;

  /*
    |GUL|GUP|.......|GUR|
    |GLF| a |...| b |GRT|
    .....................
    .....................
    .....................
    ....| c |...| d |....
    |GDL|GDW|.......|GDR|

    *(array + i*c + j)
  */

  int* gul = pointer; //size 1 [0][0]
  int* gup = pointer + 1; //size local_width [0][1]
  int* gur = pointer + field_height - 1;//size 1 [0][ncolumns-1]
  int* glf = pointer + field_width; //col_type [1][0]
  int* grt = pointer + field_width + (local_height + 1);//col_type [1][ncolumns-1]
  int* gdl = pointer + field_width*(local_height + 1); //size 1 [nrows-1][0]
  int* gdw = pointer + field_width*(local_height + 1) + 1;//local_width [nrows-1][1]
  int* gdr = pointer + field_width*(local_height + 1) + (local_height + 1); //size  [nrows-1][ncolumns-1]

  int* a = pointer + field_width + 1; //[1][1]
  int* b = pointer + field_width + local_width; //[1][jmax]
  int* c = pointer + field_width*(local_height) + 1; //[imax][1]
  int* d = pointer + field_width*(local_height) + local_width;  //[imax][jmax]

  MPI_Request reqs[8];

  MPI_Isend(a, local_width, MPI_INT, dest_rank[UP], 0, comm_cart,reqs);
  MPI_Isend(c, local_width, MPI_INT, dest_rank[DOWN], 0, comm_cart,reqs+1);
  MPI_Isend(b, 1, coltype, dest_rank[LEFT],0, comm_cart,reqs+2);
  MPI_Isend(a, 1, coltype, dest_rank[RIGHT],0, comm_cart,reqs+3);

  MPI_Isend(b, 1, MPI_INT, dest_rank[DOWNRIGHT], 0, comm_cart,reqs+4);
  MPI_Isend(c, 1, MPI_INT, dest_rank[UPLEFT], 0, comm_cart,reqs+5);
  MPI_Isend(d, 1, MPI_INT, dest_rank[UPRIGHT],0, comm_cart,reqs+6);
  MPI_Isend(a, 1, MPI_INT, dest_rank[DOWNLEFT],0, comm_cart,reqs+7);

  MPI_Waitall(8,reqs,MPI_STATUS_IGNORE);
  //pprintf("wait done for sending rank \n");

  MPI_Irecv(gul,1,MPI_INT, source_rank[UPRIGHT],0,comm_cart,reqs);
  MPI_Irecv(gur,1,MPI_INT, source_rank[UPLEFT],0,comm_cart,reqs+1);
  MPI_Irecv(gdl,1,MPI_INT, source_rank[DOWNRIGHT],0,comm_cart,reqs+2);
  MPI_Irecv(gdr,1,MPI_INT, source_rank[DOWNLEFT],0,comm_cart,reqs+3);

  MPI_Irecv(gdw,local_width,MPI_INT, source_rank[DOWN],0,comm_cart,reqs+4);
  MPI_Irecv(gup,local_width,MPI_INT, source_rank[UP],0,comm_cart,reqs+5);
  MPI_Irecv(grt,1,coltype, source_rank[LEFT],0,comm_cart,reqs+6);
  MPI_Irecv(glf,1,coltype, source_rank[RIGHT],0,comm_cart,reqs+7);

  MPI_Waitall(8,reqs,MPI_STATUS_IGNORE);
  //pprintf("wait done for receiving \n");

/*if(rank == 4){
  pprintf("%i sends up to %i\n",rank, dest_rank[UP]);
  pprintf("%i sends down to %i\n",rank, dest_rank[DOWN]);
  pprintf("%i sends right to %i\n",rank, dest_rank[LEFT]);
  pprintf("%i sends left to %i\n",rank, dest_rank[RIGHT]);
  pprintf("%i sends upright to %i\n",rank, dest_rank[DOWNRIGHT]);
  pprintf("%i sends downleft to %i\n",rank, dest_rank[UPLEFT]);
  pprintf("%i sends downright to %i\n",rank, dest_rank[UPRIGHT]);
  pprintf("%i sends upleft to %i\n",rank, dest_rank[DOWNLEFT]);

  pprintf("%i recieves up from %i\n",rank, source_rank[UP]);
  pprintf("%i recieves down from %i\n",rank, source_rank[DOWN]);
  pprintf("%i recieves left from %i\n",rank, source_rank[RIGHT]);
  pprintf("%i recieves right from %i\n",rank, source_rank[LEFT]);
  pprintf("%i recieves downleft from %i\n",rank, source_rank[DOWNRIGHT]);
  pprintf("%i recieves upright from %i\n",rank, source_rank[UPLEFT]);
  pprintf("%i recieves upleft from %i\n",rank, source_rank[UPRIGHT]);
  pprintf("%i recieves downright from %i\n",rank, source_rank[DOWNLEFT]);
}
*/

}//sync


void sync(int i)
{

  MPI_Datatype row;
  MPI_Type_vector(field_width, 1, 1, MPI_INT, &row);
  MPI_Type_commit(&row);

  //pprintf("updating ghost rows\n");
  int *pointer = (i%2==0)?field_a:field_b;


  int *ghup = pointer ;
  int *up = pointer + field_width ;
  int *dw = pointer + field_width * local_height ;
  int *ghdw = pointer + field_width * (local_height + 1);

  /*
    |GU|.............
    |UP|.............
    .................
    .................
    .................
    |DW|.............
    |GD|.............

  */

    //calculate neighbors

  int dest_up = rank - ncols;
  int dest_down = rank + ncols;


  // to prevent deadlocks
  if ((rank)%2 == 0)
    {
      if (rank != nrows -1)
	{
	  MPI_Send(dw , 1, row, dest_down , 0, MPI_COMM_WORLD);
	  MPI_Recv(ghdw , 1, row, dest_down  , 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

      if (rank != 0 )
	{
	  MPI_Recv(ghup , 1, row, dest_up , 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  MPI_Send(up, 1, row, dest_up , 0, MPI_COMM_WORLD);
	}
    }
  else
    {

      MPI_Recv(ghup , 1, row, dest_up , 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(up , 1, row, dest_up   , 0, MPI_COMM_WORLD);
      if (rank != nrows -1 )
	{
	  MPI_Send(dw , 1, row, dest_down , 0, MPI_COMM_WORLD);
	  MPI_Recv(ghdw , 1, row, dest_down , 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

    }
}

void measure(int iteration, int count_at)
{
  int total;
  int i;
  int *pointer = (iteration%2==0)?field_a:field_b;

  if (count_at == 0) //for performance
  {
   i = 0;
    for( int y=1; y<local_height+1; y++ )
    {
      for( int x=1; x<local_width+1; x++ )
      {
        if( pointer[ y * field_width + x ] )
        {
          i++;
        }
      }
    }
  }
  else
  {
    if (iteration%count_at == 0 )
    {
      i = 0;
      for( int y=1; y<local_height+1; y++ )
      {
        for( int x=1; x<local_width+1; x++ )
        {
          if( pointer[ y * field_width + x ] )
          {
            i++;
          }
        }
      }
      MPI_Allreduce( &i, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
      if( rank==0 ) pprintf( "%i buggies after iteration %i \n", total,iteration );
    }
  }
}

int main(int argc, char* argv[])
{
  char input_file[50];
  int distribution; // 0 -slice, 1- checkerboard
  int iterations ;
  int count_bugs;
  int count_at;

  // Initialize MPI
  MPI_Init(&argc, &argv);

  // Get the communicator and process information
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  MPI_Comm comm_cart;
  int dims[2];
  int periods[2] = {0,0};
  int source_ranks[8], dest_ranks[8];
  int offset;
  int write_from;
  int write_to;

  // parse args
  if (argc == 8) // will fix the parsing later
  {
    strcpy(input_file,  argv[1]);
    distribution = atoi (argv[2]);
    iterations = atoi(argv[3]);
    count_bugs = atoi(argv[4]);
    count_at = atoi(argv[5]);
    write_from = atoi(argv[6]);
    write_to = atoi(argv[7]);
  }
  else
  {
    if (rank == 0)
    {
      printf("mpirun -n 9 ./conway_base <file_name> <distribution> <iteration> <count_bugs> <count_at> <write_from> <write_to> \n");
      printf("distribution => 0 for slice, 1 for checkerboard\n");
      printf("count bugs => 0 for not counting bugs, 1 for counting bugs\n");
      printf("count at => iteration at which each count is to be displayed\n");
      printf("write_from => iteration from where it is to be written in the file\n");
      printf("write_to => iteration till where it is to be written\n");
       MPI_Finalize();
      exit(0);
    }
    else
    {
      MPI_Finalize();
      exit(0);
    }
  }

  // Print rank and hostname
  MPI_Get_processor_name(my_name, &my_name_len);
  printf("Rank %i is running on %s\n", rank, my_name );

  // Initialize the pretty printer
  init_pprintf( rank );
  pp_set_banner( "main" );
  if( rank==0 )
    pprintf( "Welcome to Conway's Game of Life!\n" );

    if (!distribution)
    {
      if (rank == 0) pprintf("\n \n SLICE distribution\n\n");
      ncols = 1;
      nrows = np;

      my_col = 0;
      my_row = rank;
      comm_cart = MPI_COMM_WORLD;
    }

  else
    {
      if (rank == 0) pprintf("\n \n CHECKERBOARD distribution\n \n");

      nrows = sqrt(np);
      ncols = sqrt(np);


      my_col = (rank/ncols );
      my_row = (rank%ncols); ///required for reading from pgm file

      dims[0] = nrows;
      dims[1] = ncols;

      MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,1,&comm_cart);
      calc_neighbours(comm_cart,rank, dims[0], dims[1],source_ranks,dest_ranks);
    }

  if (np == 1)
    {
      if (rank == 0) pprintf("single process\n");
      ncols = 1;
      nrows = 1;
      my_col = 0;
      my_row = 0;
      comm_cart = MPI_COMM_WORLD;
    }

  // Read the PGM file. The readpgm() routine reads the PGM file and, based
  // on the previously set nrows, ncols, my_row, and my_col variables, loads
  // just the local part of the field onto the current processor. The
  // variables local_width, local_height, field_width, field_height, as well
  // as the fields (field_a, field_b) are allocated and filled.
  if(!readpgm(input_file))
    {
      if( rank==0 )
      pprintf( "An error occured while reading the pgm file\n" );
      MPI_Finalize();
      return 1;
    }

    double start_measure,end_measure,start_comm,end_comm,start_update,end_update;
    MPI_Barrier(comm_cart);
    start_measure = MPI_Wtime();
    for (int i = 0; i <= iterations; i++)
    {
      // MEASURE
     // start_measure = MPI_Wtime();
      if (count_bugs)    measure(i,count_at);
      //end_measure = MPI_Wtime();

      //start_comm = MPI_Wtime();
      if(!distribution) //slice
      {
        sync(i);
      }
      else  //checkerboard
      {
        sync_checkerboard(i,dims,source_ranks,dest_ranks,comm_cart);
      }
      //end_comm = MPI_Wtime();

      //WRITE FILE
      if (i <= write_to && i >= write_from) write_to_file(input_file, i, offset);
      //UPDATE STATE
      //start_update = MPI_Wtime();
      update(i);
      //end_update = MPI_Wtime();
    }
  double local_comm,local_update,local_measure;
  double global_comm, global_update, global_measure;
//  local_comm = end_comm - start_comm;
//  local_update = end_update - start_update;
  MPI_Barrier(comm_cart);
  end_measure = MPI_Wtime();
  local_measure = end_measure - start_measure;
  //MPI_Reduce(&local_comm,&global_comm,1,MPI_DOUBLE,MPI_SUM,0,comm_cart);
  //MPI_Reduce(&local_update,&global_update,1,MPI_DOUBLE,MPI_SUM,0,comm_cart);
  MPI_Reduce(&local_measure,&global_measure,1,MPI_DOUBLE,MPI_SUM,0,comm_cart);

if (rank == 0)
  {
//    pprintf("time taken for update - %f\n", global_update);
  //  pprintf("time taken for communication - %f\n", global_comm);
   pprintf("time taken for measure - %f\n", global_measure);
    //pprintf("total time = %f\n",global_update + global_comm + global_measure);
  }
  // Free the fields
  if( field_a != NULL ) free( field_a );
  if( field_b != NULL ) free( field_b );

  // Finalize MPI and terminate
  if( rank==0 )
    pprintf( "Terminating normally\n" );

  PAPI_shutdown();
  MPI_Finalize();
  return 0;
}
