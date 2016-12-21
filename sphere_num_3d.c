#include <stdio.h>
#include <mpi.h>

int x = 7;
int y = 7;
int z = 7;
int r = 3;

int main(int argc, char *argv[])
{
  int i, j, k;
  int var;
  int center_x = x - r;
  int center_y = y - r;
  int center_z = z - r;
  //int check =

  for(i = 0; i < z; i++)
  {
    for(j = 0; j < y; j++)
    {
      for(k = 0; k < x; k++)
      {
        if((((i-3)*(i-3)) + ((j-3)*(j-3)) + ((k-3)*(k-3))) <= r*r)
          var = 1;
        else
          var = 0;

        printf("%f ", (double)((k-3)*(k-3)));
      }
      printf("\n");
    }
    printf("\n");
  }
  printf("\n");


  return 0;
}
