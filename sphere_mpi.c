#include <stdio.h>
#include <mpi.h>
#include "common.h"

int x = 7;
int y = 7;
int z = 7;

void decompose(int n, int dim, int coord, int* start, int* end)
{
    int length, rest;

    length = n/dim;
    rest = n%dim;
    *start = coord * length + (coord < rest ? coord : rest);
    *end = *start + length - (coord < rest ? 0 : 1);

    /*if((*end >= n) || (coord == dim-1))
        *end = n-1;*/
}

int main(int argc, char *argv[])
{
  int i, j, k, in, jn, kn, nx, ny, nz;
  int var, size, rank;
  int left, right, down, up, z_up, z_down;
  int x0, x1, y0, y1, z0, z1, start_x, start_y, start_z;
  int result = 0;
  int inside;
  int coords[3];
  int periods[3];
  int dims[3];
  int nd;
  int procs_x, procs_y, procs_z;
  int reorder = 0;
  //double ***Unew, ***Uold, ***test_arr;

  MPI_Comm comm3d;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int r = 5;
  int n = r*2+1;
  int d = r*2+1;

  int center_x = r+1;
  int center_y = r+1;
  int center_z = r+1;

  procs_x = 2;
  procs_y = 2;
  procs_z = 2;

  dims[0] = procs_x;
  dims[1] = procs_y;
  dims[2] = procs_z;

  if(procs_x*procs_y*procs_z != size)
  {
    if(rank == 0)
      fprintf(stderr, "Product of grid block dimensions must match the number of processes\n");
  }

  MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, reorder, &comm3d);
  MPI_Cart_get(comm3d, 3, dims, periods, coords);
  MPI_Cart_shift(comm3d, 0, 1, &left, &right);
  MPI_Cart_shift(comm3d, 1, 1, &up, &down);
  MPI_Cart_shift(comm3d, 2, 1, &z_up, &z_down);

  decompose(n, dims[0], coords[0], &x0, &x1);
  decompose(n, dims[1], coords[1], &y0, &y1);
  decompose(n, dims[2], coords[2], &z0, &z1);

  start_x = 1;
  start_y = 1;
  start_z = 1;

  if(z1 == d-1)
    z1++;
  if(y1 == d-1)
    y1++;
  if(x1 == d-1)
    x1++;


  if(x0 == 0)
    x0 = 1;
  if(y0 == 0)
    y0 = 1;
  if(z0 == 0)
    z0 = 1;

  if(left >= 0)
  {
    x0--;
    start_x = 0;
  }
  if(right >= 0)
    x1++;
  if(up >= 0)
  {
    y0--;
    start_y = 0;
  }
  if(down >= 0)
  {
    y1++;
  }
  if(z_up >= 0)
  {
    z0--;
    start_z = 0;
  }
  if(z_down >= 0)
    z1++;

  nx = x1 - x0 +1;
  ny = y1 - y0 +1;
  nz = z1 - z0 +1;



  double ***Unew, ***Uold, ***test_arr;
  double ***tensor_x, ***tensor_y, ***tensor_z;
  Unew = dallocate_3d(nx+2, ny+2, nz+2);
  Uold = dallocate_3d(nx+2, ny+2, nz+2);
  test_arr = dallocate_3d(n+2, n+2, n+2);
  tensor_x = dallocate_3d(nx+2, ny+2, nz+2);
  tensor_y = dallocate_3d(nx+2, ny+2, nz+2);
  tensor_z = dallocate_3d(nx+2, ny+2, nz+2);
  dinit_3d(Unew, nx+2, ny+2, nz+2);
  dinit_3d(Uold, nx+2, ny+2, nz+2);
  dinit_3d(test_arr, n+2, n+2, n+2);
  dinit_3d(tensor_x, nx+2, ny+2, nz+2);
  dinit_3d(tensor_y, nx+2, ny+2, nz+2);
  dinit_3d(tensor_z, nx+2, ny+2, nz+2);

  in = start_z;
  for(i = z0; i <= z1; i++)
  {
    jn = start_y;
    for(j = y0; j <= y1; j++)
    {
      kn = start_x;
      for(k = x0; k <= x1; k++)
      {
        inside = (((i-center_x)*(i-center_x)) + ((j-center_y)*(j-center_y)) + ((k-center_z)*(k-center_z)));

        if(inside <= r*r)
          Uold[in][jn][kn] = 1.0;
        else
          Uold[in][jn][kn] = 0.0;

        //Uold[in][jn][kn] = test_arr[i][j][k];

        kn++;
      }

      jn++;
    }

    in++;
  }

  if(left >= 0)
  {
    x0++;
    start_x = 0;
  }
  if(right >= 0)
    x1--;
  if(up >= 0)
  {
    y0++;
    start_y = 0;
  }
  if(down >= 0)
  {
    y1--;
  }
  if(z_up >= 0)
  {
    z0++;
    start_z = 0;
  }
  if(z_down >= 0)
    z1--;

  nx = x1 - x0 +1;
  ny = y1 - y0 +1;
  nz = z1 - z0 +1;

  for(i = 1; i <= nz; i++)
  {
    for(j = 1; j <= ny; j++)
    {
      for(k = 1; k <= nx; k++)
      {
        inside = (((i-center_x)*(i-center_x)) + ((j-center_y)*(j-center_y)) + ((k-center_z)*(k-center_z)));

        //if(inside <= r*r)
        //{
          if((Uold[i-1][j][k] != 0 && Uold[i+1][j][k] != 0) && (Uold[i][j-1][k] != 0 && Uold[i][j+1][k] != 0) && (Uold[i][j][k-1] != 0 && Uold[i][j][k+1] != 0))
          {
            tensor_x[i][j][k] = 23;
            tensor_y[i][j][k] = 23;
            tensor_z[i][j][k] = 23;
          }
          else
          {
            tensor_x[i][j][k] = 0.0;
            tensor_y[i][j][k] = 0.0;
            tensor_z[i][j][k] = 0.0;
          }
        /*}
        else
        {
          tensor_x[i][j][k] = 0.0;
          tensor_y[i][j][k] = 0.0;
          tensor_z[i][j][k] = 0.0;
        }*/
      }
    }
  }

  if(rank == 6)
  {
    for(i = 1; i <= nz; i++)
    {
      for(j = 1; j <= ny; j++)
      {
        for(k = 1; k <= nx; k++)
        {
          printf("%0.1f \t", tensor_x[i][j][k]);
        }
        printf("\n");
      }
      printf("\n");
    }
    printf("\n");

    /*for(i = 0; i <= d; i++)
    {
      for(j = 0; j <= d; j++)
      {
        for(k = 0; k <= d; k++)
        {
          printf("%0.1f ", test_arr[i][j][k]);
        }
        printf("\n");
      }
      printf("\n");
    }
    printf("\n");*/

    printf("rank: %d, l,r,u,d,zup,zdown(%d, %d, %d, %d, %d, %d) x0:%d x1:%d \n", rank, left, right, up, down, z_up, z_down, y0, y1);
    printf("rank: %d \t x0,x1: (%d,%d) \t y0,y1: (%d,%d) \t z0,z1:(%d,%d) nx,ny,nz:(%d,%d,%d)\n", rank, x0, x1, y0, y1, z0, z1, nx, ny, nz);
  }

  //printf("rank: %d \t x0,x1: (%d,%d) \t y0,y1: (%d,%d) \t z0,z1:(%d,%d) \n", rank, x0, x1, y0, y1, z0, z1);


  //printf("rank: %d, %d %d %d \n", rank, in, jn, kn);
  //int check =

  /*for(i = 0; i <= d; i++)
  {
    for(j = 0; j <= d; j++)
    {
      for(k = 0; k <= d; k++)
      {
        if((((i-3)*(i-3)) + ((j-3)*(j-3)) + ((k-3)*(k-3))) <= r*r)
        {
          Uold[i][j][k] = 1.0;
        }
        else
          Uold[i][j][k] = 0.0;

        if(k >= 0 && j == 0 && i == 0)
          result = var * stress_x;

        if(k == 0 && j >= 0 && i == 0)
          result = var * stress_y;

        if(k == 0 && j == 0 && i >= 0)
          result = var * stress_z;

        printf("%0.1f ", Uold[i][j][k]);
      }
      printf("\n");
    }
    printf("\n");
  }
  printf("\n");

  for(i = 1; i <= d; i++)
  {
    for(j = 1; j <= d; j++)
    {
      for(k = 1; k <= d; k++)
      {
        if((((i-3)*(i-3)) + ((j-3)*(j-3)) + ((k-3)*(k-3))) <= r*r)
        {
          Unew[i][j][k] = stress_x/2*d*d*(Uold[i+1][j][k] + Uold[i-1][j][k]) +
            stress_y/2*d*d*(Uold[i][j+1][k] + Uold[i][j-1][k]) +
            stress_z/2*d*d*(Uold[i][j][k+1] + Uold[i][j][k-1]);
        }
      }
    }
  }*/

  MPI_Finalize();

  return 0;
}
