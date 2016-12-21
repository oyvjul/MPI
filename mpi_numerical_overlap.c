#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>

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

void init(double ***E, int rank, int x, int y, int z)
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

void set_values(double ***E, int rank, int x, int y, int z)
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

void init_values(double ***E, int left, int right, int up, int down, int z_up,
  int z_down, int x, int y, int z, int temp_edge, int temp_interior)
{
  int i, j, k;

  if(right < 0)
  {
    for(i = 0; i < z; i++)
      for(j = 0; j < y; j++)
        E[i+1][j+1][x+1] = temp_edge;
  }
  else
  {
    for(i = 0; i < z; i++)
      for(j = 0; j < y; j++)
        E[i+1][j+1][x+1] = 0;
  }

  if(left < 0)
  {
  for(i = 0; i < z; i++)
    for(j = 0; j < y; j++)
      E[i+1][j+1][0] = temp_edge;
  }
  else
  {
    for(i = 0; i < z; i++)
      for(j = 0; j < y; j++)
        E[i+1][j+1][0] = 0;
  }

  if(up < 0)
  for(i = 0; i < z; i++)
    for(j = 0; j < x; j++)
      E[i+1][0][j+1] = temp_edge;

  if(down < 0)
  for(i = 0; i < z; i++)
    for(j = 0; j < x; j++)
      E[i+1][y+1][j+1] = temp_edge;

  if(z_up < 0)
  for(i = 0; i < y; i++)
    for(j = 0; j < x; j++)
      E[0][i+1][j+1] = temp_edge;

  if(z_down < 0)
  for(i = 0; i < y; i++)
    for(j = 0; j < x; j++)
      E[z+1][i+1][j+1] = temp_edge;


  for(i = 1; i < z+1; i++)
  {
    for(j = 1; j < y+1; j++)
    {
      for(k = 1; k < x+1; k++)
      {
        E[i][j][k] = temp_interior;
      }
    }
  }
}

void decompose(int n, int dim, int coord, int* start, int* end)
{
    int length, rest;

    length = n/dim;
    rest = n%dim;
    //printf("rest: %d, coord: %d \n", rest, coord);
    *start = coord * length + (coord < rest ? coord : rest);
    *end = *start + length - (coord < rest ? 0 : 1);

    if((*end >= n) || (coord == dim-1))
        *end = n-1;
}

double min(double a, double b)
{
    return (a<b) ? a : b;
}

void show_grid_all(double ***E, int x, int y, int z)
{
  int i, j, k;

  for(i = 0; i < z+2; i++)
  {
    for(j = 0; j < y+2; j++)
    {
      for(k = 0; k < x+2; k++)
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

  for(i = 1; i < z+1; i++)
  {
    for(j = 1; j < y+1; j++)
    {
      for(k = 1; k < x+1; k++)
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

  for(i = 1; i < z+1; i++)
  {
    for(j = 1; j < y+1; j++)
    {
      for(k = 1; k < x+1; k+=x-1)
      {
        printf("%0.1f\t", E[i][j][k]);
      }
      printf("\n");
    }
    printf("\n");
  }
  printf("\n");
}

void show_grid_up_down(double ***E, int x, int y, int z)
{
  int i, j, k;

  for(i = 1; i < z+1; i++)
  {
    for(j = 1; j < y+1; j+=2)
    {
      for(k = 1; k < x+1; k++)
      {
        printf("%0.1f\t", E[i][j][k]);
      }
      printf("\n");
    }
    printf("\n");
  }
  printf("\n");
}

void show_grid_zup_zdown(double ***E, int x, int y, int z)
{
  int i, j, k;

  for(i = 1; i < z+1; i+=2)
  {
    for(j = 1; j < y+1; j++)
    {
      for(k = 1; k < x+1; k++)
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
  int x = 32;
  int y = 32;
  int z = 32;
  int nx, ny, nz;
  int i, j, k;
  int num_iters, t;
  double result;
  int sum_iters;
  double test = 0;
  double c0=0.5;
  double c1=-0.25;
  int step;
  double dt, dt2, hx, hy, hz, min1, min2;
  double diagx, diagy, diagz, weightx, weighty, weightz, rk;

  double dt1 = 1.0e-1;

  /* temp1_init: temperature init on edges */
  double temp1_init = 10.0;

  /* temp2_init: temperature init inside */
  double temp2_init = -10.0;

  int convergence = 0;
  double k0 = 1;

  double epsilon = 1.00e-1;

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

  /* SET NUMBER OF PROCESSES IN X, Y & Z DIMENSION. MUST MATCH -NP *NUMBER_OF_PROCS* */
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

  hx = 1.0/(double)(x+2);
  hy = 1.0/(double)(y+2);
  hz = 1.0/(double)(z+2);
  min1 = min(hx,hy);
  min2 = min(min1,hz);
  dt2  = 0.125*min2*min2*min2/k0;

  if(dt1 >= dt2)
  {
      dt = dt2;
      if(rank == 0)
        printf("Time step too large,Taking convergence criterion\n");

  }
  else
    dt=dt1;

  //if(rank == 0)
    //printf("min1: %f, min2 %f, dt1: %f, dt2: %f, hx: %d", min1, min2, dt1, dt2, hx);

  Uold = allocate_3d(nx+2, ny+2, nz+2);
  init(Uold, rank, nx+2, ny+2, nz+2);
  //set_values(Uold, rank, nx+2, ny+2, nz+2);
  init_values(Uold, left, right, up, down, z_up, z_down, nx, ny, nz, temp1_init, temp2_init);

  Unew = allocate_3d(nx+2, ny+2, nz+2);
  init(Unew, rank, nx+2, ny+2, nz+2);
  //set_values(Unew, rank, nx+2, ny+2, nz+2);
  //init_values(Unew, left, right, up, down, z_up, z_down, nx, ny, nz, temp1_init, temp2_init);

  //sender
  double *sbuf_left = (double*)calloc(1,ny*nz*sizeof(double));
  double *sbuf_right = (double*)calloc(1,ny*nz*sizeof(double));
  double *sbuf_down = (double*)calloc(1,nx*nz*sizeof(double));
  double *sbuf_up = (double*)calloc(1,nx*nz*sizeof(double));
  double *sbuf_zdown = (double*)calloc(1,nx*ny*sizeof(double));
  double *sbuf_zup = (double*)calloc(1,nx*ny*sizeof(double));

  //receiver
  double *rbuf_left = (double*)calloc(1,ny*nz*sizeof(double));
  double *rbuf_right = (double*)calloc(1,ny*nz*sizeof(double));
  double *rbuf_down = (double*)calloc(1,nx*nz*sizeof(double));
  double *rbuf_up = (double*)calloc(1,nx*nz*sizeof(double));
  double *rbuf_zdown = (double*)calloc(1,nx*ny*sizeof(double));
  double *rbuf_zup = (double*)calloc(1,nx*ny*sizeof(double));

  t = 0.0;
  num_iters = 0;
  result = 0.0;
  sum_iters = 0;
  double reduce_res = 0.0;
  int iters = 0;
  step = 0;
  int max_step = 30;
  rk = 0.0;

  while(result < epsilon || step < max_step)
  {
    num_iters++;
    step = step + 1;
    t = t + dt ;
    /*if(rank == 0)
      printf("%f \n", dt);

    if(rank == 0)
      printf("%f %f %f %f %f \n", k0, hz, hx, hy, dt);*/

      diagx = - 2.0 + hx*hx/(3*k0*dt);
      weightx = k0* dt/(hx*hx);
      diagy = - 2.0 + hy*hy/(3*k0*dt);
      weighty = k0* dt/(hy*hy);
      diagz = - 2.0 + hz*hz/(3*k0*dt);
      weightz = k0* dt/(hz*hz);

    for(i = 1; i < nz+1; i++)
    {
      for(j = 1; j < ny+1; j++)
      {
        for(k = 1; k < nx+1; k+=nx-1)
        {
          Unew[i][j][k] = weightz *( Uold[i-1][j][k] + Uold[i+1][j][k] +
            Uold[i][j][k]*diagx) + weighty *( Uold[i][j-1][k] + Uold[i][j+1][k] +
            Uold[i][j][k]*diagy) + weightx *( Uold[i][j][k-1] + Uold[i][j][k+1] +
            Uold[i][j][k]*diagz);
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

    for(i = 1; i < nz+1; i++)
    {
      for(j = 1; j < ny+1; j+=ny-1)
      {
        for(k = 1; k < nx+1; k++)
        {
          Unew[i][j][k] = weightz *( Uold[i-1][j][k] + Uold[i+1][j][k] +
            Uold[i][j][k]*diagx) + weighty *( Uold[i][j-1][k] + Uold[i][j+1][k] +
            Uold[i][j][k]*diagy) + weightx *( Uold[i][j][k-1] + Uold[i][j][k+1] +
            Uold[i][j][k]*diagz);
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

    for(i = 1; i < nz+1; i+=nz-1)
    {
      for(j = 1; j < ny+1; j++)
      {
        for(k = 1; k < nx+1; k++)
        {
          Unew[i][j][k] = weightz *( Uold[i-1][j][k] + Uold[i+1][j][k] +
            Uold[i][j][k]*diagx) + weighty *( Uold[i][j-1][k] + Uold[i][j+1][k] +
            Uold[i][j][k]*diagy) + weightx *( Uold[i][j][k-1] + Uold[i][j][k+1] +
            Uold[i][j][k]*diagz);
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
          Unew[i][j][k] = weightz *( Uold[i-1][j][k] + Uold[i+1][j][k] +
            Uold[i][j][k]*diagx) + weighty *( Uold[i][j-1][k] + Uold[i][j+1][k] +
            Uold[i][j][k]*diagy) + weightx *( Uold[i][j][k-1] + Uold[i][j][k+1] +
            Uold[i][j][k]*diagz);
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

    rk = 0.0;

    for(i = 1; i < nz+1; i++)
    {
      for(j = 1; j < ny+1; j++)
      {
        for(k = 1; k < nx+1; k++)
        {
          rk = Uold[i][j][k] - Unew[i][j][k];
          result = result + rk * rk;
        }
      }
    }

    MPI_Allreduce(&result, &reduce_res, 1, MPI_DOUBLE, MPI_SUM, comm3d);

    reduce_res= sqrt(reduce_res);


    double*** tmp;
    tmp = Uold;
    Uold = Unew;
    Unew = tmp;

    if ((result<epsilon) || (step>max_step)) break;
  }

  //MPI_Reduce(&test, &reduce_res, 1, MPI_DOUBLE, MPI_SUM, 0, comm3d);

  t2 = MPI_Wtime();

  if(rank == 0)
  {
    //show_grid_all(Unew, nx, ny, nz);
    printf("\n");
    printf("  Time step = %3.18f\n", dt);
    printf("number of iters in while: %d \n", step);
    printf("num_iters: %d, result per proc: %0.1f \n", iters, reduce_res);
    //printf("rank: %d, coord: %d, (%d, %d) (%d, %d) (%d, %d)\n", rank, coords[0], x0, x1, y0, y1, z0, z1);
    printf("rank %d, l,r,d,u,zu,zd= (%d, %d, %d, %d, %d, %d)\n", rank, left, right, down, up, z_up, z_down);
    printf("rank :%d x:(%d, %d), y:(%d, %d), z:(%d, %d), nx:%d, ny:%d, nz%d\n",rank, x0, x1, y0, y1, z0, z1, nx, ny, nz);
  }

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
