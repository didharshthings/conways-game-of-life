// Conway's Game of Life
// Main Executable Program
//
// CSCI 4576/5576 High Performance Scientific Computing
// Michael Oberg, modified from code supplied by Dr. Matthew Woitaszek

// siddharth singh, added logic for the assignment
// to run: mpirun -n 9 ./conway_base <file_name> <distribution> <iteration> <count_bugs> <count_at> 
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


enum direction //used for checkerboard distribution
  {
    UP,DOWN,LEFT,RIGHT,
    UPLEFT,UPRIGHT,DOWNLEFT,DOWNRIGHT
  };


void calc_neighbours(MPI_Comm comm_cart, int rank, int dimx, int dimy, int source_ranks[8], int dest_ranks[8])
{

  
  int my_coords[2];
  int neighbor_coords[2];
  int temp_ranks;
  
  MPI_Cart_coords(comm_cart, rank, 2 ,my_coords);
  
  //Up
  MPI_Cart_shift(comm_cart,1,1,&source_ranks[UP],&dest_ranks[DOWN]);
  //pprintf("rank %i sends up rank %i\n",source_ranks[UP],dest_ranks[UP]);
  
  //Down
  MPI_Cart_shift(comm_cart,1,-1,&source_ranks[DOWN],&dest_ranks[UP]);
  //pprintf("rank %i sends down rank %i\n",source_ranks[DOWN],dest_ranks[DOWN]);
  
  //Right
  MPI_Cart_shift(comm_cart,0,-1,&source_ranks[LEFT],&dest_ranks[RIGHT]);
  //pprintf("rank %i sends right rank %i\n",source_ranks[LEFT],dest_ranks[LEFT]);
  
  //Left
  MPI_Cart_shift(comm_cart,0,1,&source_ranks[RIGHT], &dest_ranks[LEFT]);
  //pprintf("rank %i sends left rank%i \n",source_ranks[RIGHT],dest_ranks[RIGHT]);
  
  
  
  //Up and Left 
  neighbor_coords[1] = (my_coords[1] + 1) % dimx;        
  neighbor_coords[0] = (my_coords[0] - 1 + dimy) % dimy;   
  MPI_Cart_rank(comm_cart, neighbor_coords, &dest_ranks[UPLEFT]);
  source_ranks[DOWNRIGHT] = dest_ranks[UPLEFT];
  // Up and Right 
  neighbor_coords[1] = (my_coords[1] + 1) % dimx;        
  neighbor_coords[0] = (my_coords[0] + 1) % dimy;        
  MPI_Cart_rank(comm_cart, neighbor_coords, &dest_ranks[UPRIGHT]);
  source_ranks[DOWNLEFT] = dest_ranks[UPRIGHT];
  // down left 
  neighbor_coords[1] = (my_coords[1] -1 + dimx) % dimx;    
  neighbor_coords[0] = (my_coords[0] - 1 + dimy) % dimy;   
  MPI_Cart_rank(comm_cart, neighbor_coords, &dest_ranks[DOWNLEFT]);
  source_ranks[UPRIGHT] = dest_ranks[DOWNLEFT];
  // down right 
  neighbor_coords[1] = (my_coords[1] -1 + dimx) % dimx;   
  neighbor_coords[0] = (my_coords[0] + 1) % dimy;        
  MPI_Cart_rank(comm_cart, neighbor_coords, &dest_ranks[DOWNRIGHT]);
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
  

  //print final neighbours
  /*
  pprintf("rank %i sends upright rank %i - recieves from rank %i \n",rank,dest_ranks[DOWNRIGHT],source_ranks[UPLEFT]);
  pprintf("rank %i sends upleft rank %i - recieves from rank %i \n",rank,dest_ranks[DOWNLEFT],source_ranks[UPRIGHT]);
  pprintf("rank %i sends downright rank %i - recieves from rank %i \n",rank,dest_ranks[UPRIGHT],source_ranks[DOWNLEFT]);
  pprintf("rank %i sends downleft rank %i - recieves from rank %i \n",rank,dest_ranks[UPLEFT],source_ranks[DOWNRIGHT]);
  
  pprintf("rank %i sends down rank %i -recieves from rank %i \n",rank,dest_ranks[DOWN],source_ranks[DOWN]);
  pprintf("rank %i sends up rank %i - recieves from rank %i \n",rank,dest_ranks[UP],source_ranks[UP]);
  pprintf("rank %i sends left rank %i - recieves from rank %i \n",rank,dest_ranks[RIGHT],source_ranks[LEFT]);
  pprintf("rank %i sends right rank %i - recieves from rank %i\n",rank,dest_ranks[LEFT],source_ranks[RIGHT]);
  */
//}


}//calculate neigbours

void filewrite (char *in_file, int iteration, int offset) {
  
  char filename[1000];
  sprintf(filename, "%d_", iteration); 
  strcat(filename, in_file);
  offset = 0;
  
  MPI_File fh;                                          
  MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
  int char_written;
  if (rank == 0) 
    {   
      FILE *fp = fopen( filename, "w+" );
      char_written = fprintf( fp, "%s\n%i %i\n%i\n", header, width, height, depth );
      fclose(fp);
      offset = char_written;
    }
  
  int *field_pointer = (iteration%2==0) ? (field_a) : (field_b);
  
  char *temp=(char *)malloc( local_width * local_height * sizeof(char));
  
  int x, y, b, ll;
  for (y = 0; y < field_height; y++ ) {             
    for(x = 0; x < field_width; x++ ) {
      if ((x != 0) && (y != 0) && (x != field_width - 1) && (y != field_height - 1)) {
        ll = (y * field_width + x);                     
        b = field_pointer[ll];
        b = (b==0)?0xFF:0;                              
        temp[x-1+(y-1)*local_width] = (char)b;
      }
    }
  }
  
  MPI_Aint extent;                                      
  MPI_Datatype etype, filetype, contig;                 
  MPI_Offset disp = offset;                             

  // this needs to be added to the displacement so each processor starts reading from
  // the right point of the file
  disp += (rank / ncols) * width * local_height + (rank % ncols) * local_width;

  etype = MPI_CHAR;                                     
  MPI_Type_contiguous(local_width, etype, &contig);     
  extent = width * sizeof(char);                        
  MPI_Type_create_resized(contig, 0, extent, &filetype); 
  MPI_Type_commit(&filetype);                           
  // writes the file
  MPI_File_set_view(fh, disp, etype, filetype, "native", MPI_INFO_NULL);
  MPI_File_write_all(fh, temp, local_width*local_height, MPI_CHAR, MPI_STATUS_IGNORE);
  MPI_File_close(&fh);

  
}


void sync_checkerboard(int iteration, int dims[2], int source_rank[8], int dest_rank[8],MPI_Comm comm_cart)
{
  
  int *pointer = (iteration%2==0)?field_a:field_b;        
  
  /*
    |UL|UP|..... |UR|
    |LF|........ |RT|
    .................
    .................
    .................
    .................
    |DL|DW|......|DR|
    
  */

  
  
  MPI_Datatype coltype;
  MPI_Type_vector(local_width, 1, field_height, MPI_INT, &coltype);
  MPI_Type_commit(&coltype);

  MPI_Datatype rowtype;
  MPI_Type_vector(local_width, 1, 1, MPI_INT, &rowtype);
  MPI_Type_commit(&rowtype);
  
  int* ul = pointer; //size 1
  int* up = pointer + 1; //size local_width
  int* lf = pointer + field_width; //col_type
  int* dl = pointer + field_width*(local_height + 1); //size 1
  int* dw = pointer + field_width*(local_height + 1) + 1;//local_width
  int* dr = pointer + field_width*field_height; //size 1
  int* ur = pointer + local_width + 1;//size 1
  int* rt = pointer + field_width + local_width + 1;//col_type
  
  

  MPI_Request reqs[8];
  
  MPI_Isend(up, 1, rowtype, dest_rank[UP], 0, comm_cart,reqs);  
  MPI_Isend(dw, 1, rowtype, dest_rank[DOWN], 0, comm_cart,reqs+1);
  MPI_Isend(rt, 1, coltype, dest_rank[LEFT],0, comm_cart,reqs+2);
  MPI_Isend(lf, 1, coltype, dest_rank[RIGHT],0, comm_cart,reqs+3);
  
  MPI_Isend(ur, 1, MPI_INT, dest_rank[DOWNRIGHT], 0, comm_cart,reqs+4);
  MPI_Isend(dl, 1, MPI_INT, dest_rank[UPLEFT], 0, comm_cart,reqs+5);
  MPI_Isend(dr, 1, MPI_INT, dest_rank[UPRIGHT],0, comm_cart,reqs+6);
  MPI_Isend(ul, 1, MPI_INT, dest_rank[DOWNLEFT],0, comm_cart,reqs+7);
  
  
  MPI_Waitall(8,reqs,MPI_STATUS_IGNORE);
  //pprintf("wait done for sending rank \n");
  
  
  MPI_Irecv(ul,1,MPI_INT, source_rank[UPRIGHT],0,comm_cart,reqs);
  MPI_Irecv(ur,1,MPI_INT, source_rank[UPLEFT],0,comm_cart,reqs+1);
  MPI_Irecv(dl,1,MPI_INT, source_rank[DOWNRIGHT],0,comm_cart,reqs+2);
  MPI_Irecv(dr,1,MPI_INT, source_rank[DOWNLEFT],0,comm_cart,reqs+3);

  MPI_Irecv(dw,1,rowtype, source_rank[DOWN],0,comm_cart,reqs+4);
  MPI_Irecv(up,1,rowtype, source_rank[UP],0,comm_cart,reqs+5);
  MPI_Irecv(rt,1,coltype, source_rank[RIGHT],0,comm_cart,reqs+6);
  MPI_Irecv(lf,1,coltype, source_rank[LEFT],0,comm_cart,reqs+7);
  
  MPI_Waitall(8,reqs,MPI_STATUS_IGNORE);
  //pprintf("wait done for receiving \n");
  
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
  
  if (iteration%count_at == 0 )
    {
      
      i = 0;
      for( int y=1; y<local_height+1; y++ )
	{
	  for( int x=1; x<local_width+1; x++ )
	    {
	      if( field_a[ y * field_width + x ] )
		{
		  i++;
		}
	    }
	}
      
      MPI_Allreduce( &i, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
      if( rank==0 )
	pprintf( "%i buggies after iteration %i \n", total,iteration );
    }
}

int main(int argc, char* argv[])
{
  
  int sum;
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
  
  
  for (int i = 1; i <= iterations; i++)
    
    {
      if(!distribution) //slice
    	{
	       sync(i);
  	   }
      
      else  //checkerboard
	     {
	       sync_checkerboard(i,dims,source_ranks,dest_ranks,comm_cart);
      }
      // MEASURE
      
      if (count_bugs)    measure(i,count_at);
      
      //WRITE FILE
      filewrite(input_file, i, offset);
      //UPDATE STATE
      
      int *pointer_old = (i%2==0)?field_a:field_b;
      int *pointer_new = (i%2==0)?field_b:field_a;        
      
      for (int y = 0; y < local_height; y++)
	{
	  for (int x = 0; x < local_width; x++)
	    {
	      int yb = y+1;
	      int xb = x+1;
	      sum  = 0;
	      sum = pointer_old[(yb-1)*field_width+xb+1] + pointer_old[(yb-1)*field_width+xb] + pointer_old[(yb-1)*field_width+xb-1]
                +pointer_old[(yb)*field_width+xb+1]  +pointer_old[(yb)*field_width+xb-1] +pointer_old[(yb+1)*field_width+xb+1]
                + pointer_old[(yb+1)*field_width+xb] +pointer_old[(yb+1)*field_width+xb-1];
	      
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
  
  
  // Free the fields
  if( field_a != NULL ) free( field_a );
  if( field_b != NULL ) free( field_b );
  
  // Finalize MPI and terminate
  if( rank==0 )
    pprintf( "Terminating normally\n" );
  MPI_Finalize();
  return 0;
}
