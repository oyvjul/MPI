#include <stdio.h>
#include <mpi.h>
#include "common.h"
#include "math.h"

int x = 7;
int y = 7;
int z = 7;

void calculatel2Norm(double*** E, int nx, int ny, int nz, int r, int iters, int rank, int d)
{
  int i, j, k;

  float mx = -1;
  float l2norm = 0;
  float l2norm_all = 0;

  for (i = 1; i <= nz; i++)
  {
      for (j = 1; j <= ny; j++)
      {
          for (k = 1; k <= nx; k++)
          {
            l2norm += E[i][j][k]*E[i][j][k];

            /*if (E[i][j][k] > mx)
            mx = E[i][j][k];*/
            }
        }
  }

  MPI_Allreduce(&l2norm, &l2norm_all, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  if(rank == 0)
    printf("Before sqrt, %20.12e \n", l2norm_all);

  l2norm_all = sqrt(l2norm_all);
  //l2norm_all /= (float) ((d)*(d)*(d));
  //l2norm_all /= (float) ((nx)*(ny)*(nz));

  if(rank == 0)
  {
    printf("radius: %d, iteration %d \n", r, iters);
    printf("max: %20.12e, l2norm: %20.12e \n", mx, l2norm_all);
  }
}

void show_grid_left_right1(double ***E, int x, int y, int z)
{
  int i, j, k;

  for(i = 0; i <= z+1; i++)
  {
    for(j = 0; j <= y+1; j++)
    {
      for(k = 0; k <= x+1; k+=x+1)
      {
        printf("%0.1f\t", E[i][j][k]);
      }
      printf("\n");
    }
    printf("\n");
  }
  printf("\n");
}

void show_grid_left_right(double ***E, int x, int y, int z)
{
  int i, j, k;

  for(i = 0; i <= z+1; i++)
  {
    for(j = 0; j <= y+1; j++)
    {
      for(k = 0; k <= x+1; k+=x+1)
      {
        E[i][j][k] = 0.0;
      }
    }
  }
}

void show_grid_up_down(double ***E, int x, int y, int z)
{
  int i, j, k;

  for(i = 0; i <= z+1; i++)
  {
    for(j = 0; j <= y+1; j+=y+1)
    {
      for(k = 0; k <= x+1; k++)
      {
        E[i][j][k] = 0.0;
      }
    }
  }
}

void show_grid_zup_zdown(double ***E, int x, int y, int z)
{
  int i, j, k;

  for(i = 0; i <= z+1; i+=z+1)
  {
    for(j = 0; j <= y+1; j++)
    {
      for(k = 0; k <= x+1; k++)
      {
        E[i][j][k] = 0.0;
      }
    }
  }
}

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
  MPI_Request req[12];

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  nd = 0;
  int r = 6;
  int n = r*2+1;
  int d = r*2+1;

  //MPI_Allreduce(&r, &nd, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  double size_x = d;
  double size_y = d;
  double size_z = d;

  double dx = 1.0/size_x;
  double dy = 1.0/size_y;
  double dz = 1.0/size_z;

  int center_x = r+1;
  int center_y = r+1;
  int center_z = r+1;

  procs_x = 4;
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



  double ***Unew, ***Uold, ***test_arr, ***temp;
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

  double *sbuf_left = (double*)calloc(1,ny*nz*sizeof(double));
  double *sbuf_right = (double*)calloc(1,ny*nz*sizeof(double));
  double *sbuf_down = (double*)calloc(1,nx*nz*sizeof(double));
  double *sbuf_up = (double*)calloc(1,nx*nz*sizeof(double));
  double *sbuf_zdown = (double*)calloc(1,nx*ny*sizeof(double));
  double *sbuf_zup = (double*)calloc(1,nx*ny*sizeof(double));

  double *rbuf_left = (double*)calloc(1,ny*nz*sizeof(double));
  double *rbuf_right = (double*)calloc(1,ny*nz*sizeof(double));
  double *rbuf_down = (double*)calloc(1,nx*nz*sizeof(double));
  double *rbuf_up = (double*)calloc(1,nx*nz*sizeof(double));
  double *rbuf_zdown = (double*)calloc(1,nx*ny*sizeof(double));
  double *rbuf_zup = (double*)calloc(1,nx*ny*sizeof(double));

  double count = 1.0;
  int test = 1;
  for(i = 1; i <= d; i++)
  {
    for(j = 1; j <= d; j++)
    {
      for(k = 1; k <= d; k++)
      {
        test_arr[i][j][k] = count;
        count++;
      }
    }
  }
  double sum = 0;
  double sum_res = 0;

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
        {
          Uold[in][jn][kn] = test_arr[i][j][k];
          //Unew[in][jn][kn] = test_arr[i][j][k];
          //Uold[in][jn][kn] = 1;

          count++;
        }
        else
        {
          Uold[in][jn][kn] = 0.0;
          //Unew[in][jn][kn] = 0.0;
        }

        //Uold[in][jn][kn] = test_arr[i][j][k];
        //Unew[in][jn][kn] = test_arr[i][j][k];
        //sum += Uold[in][jn][kn];
        //count++;
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
            tensor_x[i][j][k] = 59.10;
            tensor_y[i][j][k] = 2;
            tensor_z[i][j][k] = 2;
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

        //Uold[in][jn][kn] = test_arr[i][j][k];
      }
    }
  }



  //show_grid_left_right(Uold, nx, ny, nz);
  //show_grid_up_down(Uold, nx, ny, nz);
  //show_grid_zup_zdown(Uold, nx, ny, nz);
  double l2_new = 0;
  double all_sum = 0;

  if(rank == 0)
  {
    for(i = 0; i <= nz+1; i++)
    {
      for(j = 0; j <= ny+1; j++)
      {
        for(k = 0; k <= nx+1; k++)
        {
          printf("%0.1f \t", Unew[i][j][k]);
        }
        printf("\n");
      }
      printf("\n");
    }
    printf("\n");

    printf("rank: %d, l,r,u,d,zup,zdown(%d, %d, %d, %d, %d, %d) x0:%d x1:%d \n", rank, left, right, up, down, z_up, z_down, y0, y1);
    printf("rank: %d \t x0,x1: (%d,%d) \t y0,y1: (%d,%d) \t z0,z1:(%d,%d) nx,ny,nz:(%d,%d,%d)\n", rank, x0, x1, y0, y1, z0, z1, nx, ny, nz);
  }

  int max_time = 12;
  int time_iter = 0;

  while(time_iter < max_time)
  {


    for(i = 1; i <= nz; i++)
    {
      for(j = 1; j <= ny; j++)
      {
        for(k = 1; k <= nx; k+=nx)
        {
          Unew[i][j][k] = ((tensor_x[i][j][k]/2*dx*dx)*(Uold[i+1][j][k] - Uold[i-1][j][k])) +
            ((tensor_y[i][j][k]/2*dy*dy)*(Uold[i][j+1][k] - Uold[i][j-1][k])) +
            ((tensor_z[i][j][k]/2*dz*dz)*(Uold[i][j][k+1] - Uold[i][j][k-1]));

        }
      }
    }

    if(left >= 0)
    for(i = 0; i < nz; i++)
      for(j = 0; j < ny; j++)
        sbuf_left[i*ny+j] = Unew[i+1][j+1][1];

    MPI_Isend(sbuf_left, nz*ny, MPI_DOUBLE, left, TAG1, comm3d, &req[0]);
    MPI_Irecv(rbuf_right, nz*ny, MPI_DOUBLE, right, TAG1, comm3d, &req[1]);

    if(right >= 0)
    for(i = 0; i < nz; i++)
      for(j = 0; j < ny; j++)
        sbuf_right[i*ny+j] = Unew[i+1][j+1][nx];

    MPI_Isend(sbuf_right, nz*ny, MPI_DOUBLE, right, TAG2, comm3d, &req[2]);
    MPI_Irecv(rbuf_left, nz*ny, MPI_DOUBLE, left, TAG2, comm3d, &req[3]);

    for(i = 1; i <= nz; i++)
    {
      for(j = 1; j <= ny; j+=ny)
      {
        for(k = 1; k <= nx; k++)
        {
          Unew[i][j][k] = ((tensor_x[i][j][k]/2*dx*dx)*(Uold[i+1][j][k] - Uold[i-1][j][k])) +
            ((tensor_y[i][j][k]/2*dy*dy)*(Uold[i][j+1][k] - Uold[i][j-1][k])) +
            ((tensor_z[i][j][k]/2*dz*dz)*(Uold[i][j][k+1] - Uold[i][j][k-1]));
        }
      }
    }

    if(down >= 0)
    for(i = 0; i < nz; i++)
      for(j = 0; j < nx; j++)
        sbuf_down[i*nx+j] = Unew[i+1][ny][j+1];

    MPI_Isend(sbuf_down, nx*nz, MPI_DOUBLE, down, TAG3, comm3d, &req[4]);
    MPI_Irecv(rbuf_up, nx*nz, MPI_DOUBLE, up, TAG3, comm3d, &req[5]);

    if(up >= 0)
    for(i = 0; i < nz; i++)
      for(j = 0; j < nx; j++)
        sbuf_up[i*nx+j] = Unew[i+1][1][j+1];

    MPI_Isend(sbuf_up, nx*nz, MPI_DOUBLE, up, TAG4, comm3d, &req[6]);
    MPI_Irecv(rbuf_down, nx*nz, MPI_DOUBLE, down, TAG4, comm3d, &req[7]);

    for(i = 1; i <= nz; i+=nz)
    {
      for(j = 1; j <= ny; j++)
      {
        for(k = 1; k <= nx; k++)
        {
          Unew[i][j][k] = ((tensor_x[i][j][k]/2*dx*dx)*(Uold[i+1][j][k] - Uold[i-1][j][k])) +
            ((tensor_y[i][j][k]/2*dy*dy)*(Uold[i][j+1][k] - Uold[i][j-1][k])) +
            ((tensor_z[i][j][k]/2*dz*dz)*(Uold[i][j][k+1] - Uold[i][j][k-1]));
        }
      }
    }

    if(z_down >= 0)
    for(i = 0; i < ny; i++)
      for(j = 0; j < nx; j++)
        sbuf_zdown[i*nx+j] = Unew[nz][i+1][j+1];

    MPI_Isend(sbuf_zdown, nx*ny, MPI_DOUBLE, z_down, TAG5, comm3d, &req[8]);
    MPI_Irecv(rbuf_zup, nx*ny, MPI_DOUBLE, z_up, TAG5, comm3d, &req[9]);

    if(z_up >= 0)
    for(i = 0; i < ny; i++)
      for(j = 0; j < nx; j++)
        sbuf_zup[i*nx+j] = Unew[1][i+1][j+1];

    MPI_Isend(sbuf_zup, nx*ny, MPI_DOUBLE, z_up, TAG6, comm3d, &req[10]);
    MPI_Irecv(rbuf_zdown, nx*ny, MPI_DOUBLE, z_down, TAG6, comm3d, &req[11]);

    //Main computation
    for (i = 2; i < nz; i++)
    {
      for (j= 2; j < ny; j++)
      {
        for (k = 2; k < nx; k++)
        {
          Unew[i][j][k] = ((tensor_x[i][j][k]/2*dx*dx)*(Uold[i+1][j][k] - Uold[i-1][j][k])) +
            ((tensor_y[i][j][k]/2*dy*dy)*(Uold[i][j+1][k] - Uold[i][j-1][k])) +
            ((tensor_z[i][j][k]/2*dz*dz)*(Uold[i][j][k+1] - Uold[i][j][k-1]));
        }
      }
    }

    MPI_Waitall(12, req, MPI_STATUSES_IGNORE);

    if(right >= 0)
    for(i = 0; i < nz; i++)
      for(j = 0; j < ny; j++)
        Unew[i+1][j+1][nx+1] = rbuf_right[i*ny+j];

    if(left >= 0)
    for(i = 0; i < nz; i++)
      for(j = 0; j < ny; j++)
        Unew[i+1][j+1][0] = rbuf_left[i*ny+j];

    if(up >= 0)
    for(i = 0; i < nz; i++)
      for(j = 0; j < nx; j++)
        Unew[i+1][0][j+1] = rbuf_up[i*nx+j];

    if(down >= 0)
    for(i = 0; i < nz; i++)
      for(j = 0; j < nx; j++)
        Unew[i+1][ny+1][j+1] = rbuf_down[i*nx+j];

    if(z_up >= 0)
    for(i = 0; i < ny; i++)
      for(j = 0; j < nx; j++)
        Unew[0][i+1][j+1] = rbuf_zup[i*nx+j];

    if(z_down >= 0)
    for(i = 0; i < ny; i++)
      for(j = 0; j < nx; j++)
        Unew[nz+1][i+1][j+1] = rbuf_zdown[i*nx+j];

    temp = Uold;
    Uold = Unew;
    Unew = temp;

    time_iter++;
  }

  for(i = 1; i <= nz; i++)
  {
    for(j = 1; j <= ny; j++)
    {
      for(k = 1; k <= nx; k++)
      {
        l2_new += Uold[i][j][k]*Uold[i][j][k];
      }
    }
  }

  MPI_Allreduce(&l2_new, &all_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if(rank == 0)
  {
    printf("NEW SUM: %f \n", all_sum);
    all_sum = sqrt(all_sum);
    printf("after sqrt: %0.9f \n", all_sum);
  }

  calculatel2Norm(Uold, nx, ny, nz, r, max_time, rank, d);

  MPI_Finalize();

  return 0;
}
