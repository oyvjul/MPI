#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

#define chunk 64
#define TAG1 1
#define TAG2 2
#define TAG3 3
#define TAG4 4
#define TAG5 5
#define TAG6 6
const double kMicro = 1.0e-6;

double ***allocate_3d(int x, int y, int z)
{
    int i, j;
    double *storage = malloc(x * y * z * sizeof(*storage));
    double *alloc = storage;
    double ***matrix;
    matrix = malloc(z * sizeof(double**));

    for (i = 0; i < z; i++)
    {
        matrix[i] = malloc(y * sizeof(**matrix));

        for (j = 0; j < y; j++)
        {
            matrix[i][j] = alloc;
            alloc += x;
        }
    }

    return matrix;
}

void free3D(double*** E)
{
    free(E[0][0]);
    free(E[0]);
    free(E);
}

void init(double*** E, int rank, int x, int y, int z)
{
    int i, j, k;
    double count = 0.0;

    for(i = 0; i < z; i++)
    {
      for(j = 0; j < y; j++)
      {
        for(k = 0; k < x; k++)
        {
          //if(i == 0 || i == z-1|| j == 0 || j == y-1|| k == 0 || k == x-1)
            E[i][j][k] = 2.0;
          /*else
          {
          count++;
          E[i][j][k] = rank;
        }*/
      }
    }
  }
}

void set_values(double*** E, int rank, int x, int y, int z)
{
    int i, j, k, start_x, start_y, start_z;
    double count = 0.0;

    for(i = 1; i < z-1; i++)
        for(j = 1; j < y-1; j++)
            for(k = 1; k < x-1; k++)
            {
              //count++;
              //E[i][j][k] = count;
              E[i][j][k] = 2.0;
            }
}

void decompose(int n, int dim, int coord, int* start, int* end)
{
    int length, rest;

    length = n/dim;
    rest = n%dim;
    *start = coord * length + (coord < rest ? coord : rest);
    *end = *start + length - (coord < rest ? 0 : 1);

    if((*end >= n) || (coord == dim-1))
        *end = n-1;
}

void show_grid_all(double ***E, int x, int y, int z)
{
  int i, j, k;

  for(i = 0; i < z; i++)
  {
    for(j = 0; j < y; j++)
    {
      for(k = 0; k < x; k++)
      {
        printf("%0.1f\t", E[i][j][k]);
      }
      printf("\n");
    }
    printf("\n");
  }
  printf("\n");
}

void show_grid_interior(double ***E, int x, int y, int z)
{
  int i, j, k;

  for(i = 1; i < z-1; i++)
  {
    for(j = 1; j < y-1; j++)
    {
      for(k = 1; k < x-1; k++)
      {
        printf("%0.1f\t", E[i][j][k]);
      }
      printf("\n");
    }
    printf("\n");
  }
  printf("\n");
}

int main(int argc, char *argv[])
{

  int size, rank;
  int coords[3];
  int periods[3];
  int dims[3];
  int reorder = 0;
  int x = 6;
  int y = 6;
  int z = 6;
  int nx, ny, nz;
  int i, j, k;
  int num_iters, t;
  double result;
  int sum_iters;
  double test = 0;
  double c0=0.5;
  double c1=-0.25;
  int left, right, down, up, z_down, z_up;
  int x0, y0, z0, x1, y1, z1;
  int procs_x, procs_y, procs_z;
  double ***Unew, ***Uold, ***global;
  double *sendbuf;
  MPI_Comm comm3d;
  MPI_Status status;
  MPI_Request sendreq[6], recvreq[6], req[12];

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  double t1, t2;
  t1 = MPI_Wtime();

  procs_x = 2;
  procs_y = 2;
  procs_z = 2;

  dims[0] = procs_x;
  dims[1] = procs_y;
  dims[2] = procs_z;

  periods[0] = 0;
  periods[1] = 0;
  periods[2] = 0;

  if(procs_x*procs_y*procs_z != size)
  {
    if(rank == 0)
      fprintf(stderr, "Product of grid block dimensions must match the number of processes\n");
  }

  MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, reorder, &comm3d);
  MPI_Cart_get(comm3d, 3, dims, periods, coords);

  MPI_Cart_shift(comm3d, 0, 1, &left, &right);
  MPI_Cart_shift(comm3d, 1, 1, &down, &up);
  MPI_Cart_shift(comm3d, 2, 1, &z_up, &z_down);

  decompose(x, dims[0], coords[0], &x0, &x1);
  decompose(y, dims[1], coords[1], &y0, &y1);
  decompose(z, dims[2], coords[2], &z0, &z1);

  nx = x1 - x0 + 1;
  ny = y1 - y0 + 1;
  nz = z1 - z0 + 1;
  //printf("nx: %f", c1);

  Uold = allocate_3d(nx+2, ny+2, nz+2);
  init(Uold, rank, nx+2, ny+2, nz+2);
  set_values(Uold, rank, nx+2, ny+2, nz+2);

  Unew = allocate_3d(nx+2, ny+2, nz+2);
  init(Unew, rank, nx+2, ny+2, nz+2);
  //set_values(Unew, rank, nx+2, ny+2, nz+2);
  //set_values(Unew, rank, nx+2, ny+2, nz+2);

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

  /*if(left >= 0)
  for(i = 0; i < nz; i++)
    for(j = 0; j < ny; j++)
      sbuf_left[i*ny+j] = Unew[i+1][j+1][1];

  if(right >= 0)
  for(i = 0; i < nz; i++)
    for(j = 0; j < ny; j++)
      sbuf_right[i*ny+j] = Unew[i+1][j+1][nx];

  if(down >= 0)
  for(i = 0; i < nz; i++)
    for(j = 0; j < nx; j++)
      sbuf_down[i*nx+j] = Unew[i+1][ny][j+1];

  if(up >= 0)
  for(i = 0; i < nz; i++)
    for(j = 0; j < nx; j++)
      sbuf_up[i*nx+j] = Unew[i+1][1][j+1];

  if(z_down >= 0)
  for(i = 0; i < ny; i++)
    for(j = 0; j < nx; j++)
      sbuf_zdown[i*nx+j] = Unew[ny][i+1][j+1];

  if(z_up >= 0)
  for(i = 0; i < ny; i++)
    for(j = 0; j < nx; j++)
      sbuf_zup[i*nx+j] = Unew[1][i+1][j+1];

  MPI_Isend(sbuf_left, nz*ny, MPI_DOUBLE, left, TAG1, comm3d, &req[0]);
  MPI_Isend(sbuf_right, nz*ny, MPI_DOUBLE, right, TAG2, comm3d, &req[1]);
  MPI_Isend(sbuf_down, nx*nz, MPI_DOUBLE, down, TAG3, comm3d, &req[2]);
  MPI_Isend(sbuf_up, nx*nz, MPI_DOUBLE, up, TAG4, comm3d, &req[3]);
  MPI_Isend(sbuf_zdown, nx*ny, MPI_DOUBLE, z_down, TAG5, comm3d, &req[4]);
  MPI_Isend(sbuf_zup, nx*ny, MPI_DOUBLE, z_up, TAG6, comm3d, &req[5]);

  MPI_Irecv(rbuf_right, nz*ny, MPI_DOUBLE, right, TAG1, comm3d, &req[6]);
  MPI_Irecv(rbuf_left, nz*ny, MPI_DOUBLE, left, TAG2, comm3d, &req[7]);
  MPI_Irecv(rbuf_up, nx*nz, MPI_DOUBLE, up, TAG3, comm3d, &req[8]);
  MPI_Irecv(rbuf_down, nx*nz, MPI_DOUBLE, down, TAG4, comm3d, &req[9]);
  MPI_Irecv(rbuf_zup, nx*ny, MPI_DOUBLE, z_up, TAG5, comm3d, &req[10]);
  MPI_Irecv(rbuf_zdown, nx*ny, MPI_DOUBLE, z_down, TAG6, comm3d, &req[11]);

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
      Unew[nz+1][i+1][j+1] = rbuf_zdown[i*nx+j];*/

  t = 20;
  num_iters = 0;
  result = 0.0;
  sum_iters = 0;
  double reduce_res = 0;
  int iters = 0;

  while(num_iters < t)
  {
    num_iters++;

    for (i = 1; i < nz+1; i++)
    {
      for (j= 1; j < ny+1; j++)
      {
        for (k = 1; k < nx+1; k++)
        {
          Unew[i][j][k] = c0* Uold[i][j][k]  + c1 * (Uold[i][j][k-1] + Uold[i][j][k+1] +
              Uold[i][j-1][k] + Uold[i][j+1][k] +
              Uold[i-1][j][k] + Uold[i+1][j][k]);

          //result += Unew[i][j][k];
          iters++;
          //if(rank == 0)
            //printf("rank 0 + %d \n", num_iters);
        }
      }
    }

    if(left >= 0)
    for(i = 0; i < nz; i++)
      for(j = 0; j < ny; j++)
        sbuf_left[i*ny+j] = Unew[i+1][j+1][1];

    if(right >= 0)
    for(i = 0; i < nz; i++)
      for(j = 0; j < ny; j++)
        sbuf_right[i*ny+j] = Unew[i+1][j+1][nx];

    if(down >= 0)
    for(i = 0; i < nz; i++)
      for(j = 0; j < nx; j++)
        sbuf_down[i*nx+j] = Unew[i+1][ny][j+1];

    if(up >= 0)
    for(i = 0; i < nz; i++)
      for(j = 0; j < nx; j++)
        sbuf_up[i*nx+j] = Unew[i+1][1][j+1];

    if(z_down >= 0)
    for(i = 0; i < ny; i++)
      for(j = 0; j < nx; j++)
        sbuf_zdown[i*nx+j] = Unew[nz][i+1][j+1];

    if(z_up >= 0)
    for(i = 0; i < ny; i++)
      for(j = 0; j < nx; j++)
        sbuf_zup[i*nx+j] = Unew[1][i+1][j+1];

    MPI_Isend(sbuf_left, nz*ny, MPI_DOUBLE, left, TAG1, comm3d, &req[0]);
    MPI_Isend(sbuf_right, nz*ny, MPI_DOUBLE, right, TAG2, comm3d, &req[1]);
    MPI_Isend(sbuf_down, nx*nz, MPI_DOUBLE, down, TAG3, comm3d, &req[2]);
    MPI_Isend(sbuf_up, nx*nz, MPI_DOUBLE, up, TAG4, comm3d, &req[3]);
    MPI_Isend(sbuf_zdown, nx*ny, MPI_DOUBLE, z_down, TAG5, comm3d, &req[4]);
    MPI_Isend(sbuf_zup, nx*ny, MPI_DOUBLE, z_up, TAG6, comm3d, &req[5]);

    MPI_Irecv(rbuf_right, nz*ny, MPI_DOUBLE, right, TAG1, comm3d, &req[6]);
    MPI_Irecv(rbuf_left, nz*ny, MPI_DOUBLE, left, TAG2, comm3d, &req[7]);
    MPI_Irecv(rbuf_up, nx*nz, MPI_DOUBLE, up, TAG3, comm3d, &req[8]);
    MPI_Irecv(rbuf_down, nx*nz, MPI_DOUBLE, down, TAG4, comm3d, &req[9]);
    MPI_Irecv(rbuf_zup, nx*ny, MPI_DOUBLE, z_up, TAG5, comm3d, &req[10]);
    MPI_Irecv(rbuf_zdown, nx*ny, MPI_DOUBLE, z_down, TAG6, comm3d, &req[11]);

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

    for(i = 1; i < nz+1; i++)
    {
      for(j = 1; j < ny+1; j++)
      {
        for(k = 1; k < nx+1; k++)
        {
          rk = Uold[i][j][k] - Unew[i][j][k];
        }
      }
    }

    test += result;

    double*** tmp;
    tmp = Uold;
    Uold = Unew;
    Unew = tmp;
  }

  MPI_Reduce(&test, &reduce_res, 1, MPI_DOUBLE, MPI_SUM, 0, comm3d);
  //MPI_Reduce(&iters, &sum_iters, 1, MPI_INT, MPI_SUM, 0, comm3d);



  if(rank == 0)
  {
    show_grid_interior(Unew, nx+2, ny+2, nz+2);
    printf("num_iters: %d, result per proc: %0.1f \n", sum_iters, reduce_res);
    //printf("rank: %d, coord: %d, (%d, %d) (%d, %d) (%d, %d)\n", rank, coords[0], x0, x1, y0, y1, z0, z1);
    printf("rank %d, l,r,d,u,zu,zd= (%d, %d, %d, %d, %d, %d)\n", rank, left, right, down, up, z_up, z_down);
    printf("x:(%d, %d), y:(%d, %d), z:(%d, %d), nx:%d, ny:%d, nz%d\n", x0, x1, y0, y1, z0, z1, nx, ny, nz);
  }

  t2 = MPI_Wtime();
  printf( "Elapsed time is: %f \n", t2 - t1);

  free(sbuf_left);
  free(sbuf_right);
  free(sbuf_down);
  free(sbuf_up);
  free(sbuf_zup);
  free(sbuf_zdown);
  free(rbuf_left);
  free(rbuf_right);
  free(rbuf_down);
  free(rbuf_up);
  free(rbuf_zup);
  free(rbuf_zdown);
  free3D(Uold);
  MPI_Finalize();

  return 0;
}
