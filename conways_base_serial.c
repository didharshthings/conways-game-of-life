// Conway's Game of Life - serial
// Main Executable Program
//
// CSCI 4576/5576 High Performance Scientific Computing
// Michael Oberg, modified from code supplied by Dr. Matthew Woitaszek

// siddharth singh, added logic for the assignment
// to run: mpirun -n 9 ./conway_base <file_name> <distribution> <iteration> <count_bugs> <count_at> <write_from> <write_to>
// mpirun -n 9 ./conway_base life.pgm 0 100 1 10 0 0


// System includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// User includes
#include <math.h>

int nrows;          // Number of rows in our partitioning
int ncols;          // Number of columns in our partitioning
int my_row;         // My row number
int my_col;         // My column number

// Local logical game size
int local_width;    // Width and height of game on this processor
int local_height;

// Local physical field size
int field_width;        // Width and height of field on this processor
int field_height;       // (should be local_width+2, local_height+2)
int *field_a = NULL;      // The local data fields
int *field_b = NULL;
char header[10];
int height;
int width;
int depth;

typedef enum { false, true } bool; // Provide C++ style 'bool' type in C

bool readpgm( char *filename )
{
  // Open the file
  printf( "Opening file %s\n", filename );
  FILE *fp = fopen( filename, "r" );
  if( !fp )
  {
    printf( "Error: The file '%s' could not be opened.\n", filename );
    return false;
  }

  // Read the PGM header, which looks like this:
  //  |P5        magic version number
  //  |900 900       width height
  //  |255         depth
  //char header[10];
  //int width, height, depth;
  int rv = fscanf( fp, "%6s\n%i %i\n%i\n", header, &width, &height, &depth );
  if( rv != 4 )
  {
    printf( "Error: The file '%s' did not have a valid PGM header\n", filename );
    return false;
  }
  printf( "%s: %s %i %i %i\n", filename, header, width, height, depth );
  // Make sure the header is valid
  if( strcmp( header, "P5") )
  {
    printf( "Error: PGM file is not a valid P5 pixmap.\n" );
    return false;
  }
  if( depth != 255 )
  {
    printf( "Error: PGM file has depth=%i, require depth=255 \n", depth );
    return false;
  }

  // Make sure that the width and height are divisible by the number of
  // processors in x and y directions

  if( width % ncols )
  {
    printf( "Error: %i pixel width cannot be divided into %i cols\n", width, ncols );
    return false;
  }
  if( height % nrows )
  {
    printf( "Error: %i pixel height cannot be divided into %i rows\n",height, nrows );
    return false;
  }

  // Divide the total image among the local processors
  local_width = width / ncols;
  local_height = height / nrows;

  // Find out where my starting range it
  int start_x = local_width * my_col;
  int start_y = local_height * my_row;

  printf( "Hosting data for x:%03i-%03i y:%03i-%03i\n",start_x, start_x + local_width,
      start_y, start_y + local_height );

  // Create the array!
  field_width = local_width + 2;
  field_height = local_height + 2;
  field_a = (int *)malloc( field_width * field_height * sizeof(int));
  field_b = (int *)malloc( field_width * field_height * sizeof(int));

  // Read the data from the file. Save the local data to the local array.
  int b, ll, lx, ly;
  for( int y=0; y<height; y++ )
  {
    for( int x=0; x<width; x++ )
    {
      // Read the next character
      b = fgetc( fp );
      if( b == EOF )
      {
        printf( "Error: Encountered EoF at [%i,%i]\n", y,x );
        return false;
      }
      // From the PGM, black cells (b=0) are bugs, all other
      // cells are background
      if( b==0 )
      {
        b=1;
      }
      else
      {
        b=0;
      }
      // If the character is local, then save it!
      if( x >= start_x && x < start_x + local_width &&
        y >= start_y && y < start_y + local_height )
        {
          // Calculate the local pixels (+1 for ghost row,col)
          lx = x - start_x + 1;
          ly = y - start_y + 1;
          ll = (ly * field_width + lx );
          field_a[ ll ] = b;
          field_b[ ll ] = b;
        } // save local point
      } // for x
    } // for y
    fclose( fp );
    return true;
  }
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

void measure(int iteration, int count_at)
{
  int total;
  int i;
  int *pointer = (iteration%2==0)?field_a:field_b;

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
    printf( "%i buggies after iteration %i \n", total,iteration );
  }
}

int main(int argc, char* argv[])
{
  char input_file[50];
  int distribution; // 0 -slice, 1- checkerboard
  int iterations ;
  int count_bugs;
  int count_at;

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
    printf("mpirun -n 9 ./conway_base <file_name> <distribution> <iteration> <count_bugs> <count_at> <write_from> <write_to> \n");
    printf("distribution => 0 for slice, 1 for checkerboard\n");
    printf("count bugs => 0 for not counting bugs, 1 for counting bugs\n");
    printf("count at => iteration at which each count is to be displayed\n");
    printf("write_from => iteration from where it is to be written in the file\n");
    printf("write_to => iteration till where it is to be written\n");
    exit(0);
  }

  ncols = 1;
  nrows = 1;
  my_col = 0;
  my_row = 0;

  // Read the PGM file. The readpgm() routine reads the PGM file and, based
  // on the previously set nrows, ncols, my_row, and my_col variables, loads
  // just the local part of the field onto the current processor. The
  // variables local_width, local_height, field_width, field_height, as well
  // as the fields (field_a, field_b) are allocated and filled.
  if(!readpgm(input_file))
    {
      printf( "An error occured while reading the pgm file\n" );
      return 1;
    }

    for (int i = 0; i <= iterations; i++)
    {
      // MEASURE
      if (count_bugs)    measure(i,count_at);
      //UPDATE STATE
      update(i);
    }

    // Free the fields
  if( field_a != NULL ) free( field_a );
  if( field_b != NULL ) free( field_b );
  return 0;
}
