#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include <mpi.h>
#include <unistd.h>
#include <omp.h>

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
    matrix = malloc(x * sizeof(double**));

    for (i = 0; i < x; i++)
    {
        matrix[i] = malloc(y * sizeof(**matrix));

        for (j = 0; j < y; j++)
        {
            matrix[i][j] = alloc;
            alloc += z;
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

void init(double*** E, int N, int M, int K)
{
    int i,j,k;

    for(k=0 ; k < K ; k++)
        for(i=0 ; i < M ; i++)
            for(j=0 ; j < N ; j++)
            {
                E[k][i][j] = 1.0;
            }
}

void init2(double*** E, int N, int M, int K)
{
    int i,j,k;
    double count = 0.0;

    for(k=0 ; k < K ; k++)
        for(i=0 ; i < M ; i++)
            for(j=0 ; j < N ; j++)
            {
                E[k][i][j] = 0;
                if(i==0 || j == 0 || k==0)
                  E[k][i][j] = 0.0;
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

void create_datatype_array(int n, int m, int k, int nx, int ny, int nz,
  int start_x, int start_y, int start_z, MPI_Datatype *matrix1)
{
  int starts[3], size[3], sub_size[3];
  MPI_Datatype matrix;

  starts[0] = start_x;
  starts[1] = start_y;
  starts[2] = start_z;

  size[0] = k;
  size[1] = m;
  size[2] = n;

  sub_size[0] = nz;
  sub_size[1] = ny;
  sub_size[2] = nx;

  MPI_Type_create_subarray(3,  size, sub_size, starts, MPI_ORDER_C, MPI_DOUBLE, &matrix);
  MPI_Type_commit(&matrix);

  *matrix1 = matrix;
}

int main (int argc, char* argv[])
{

    int n = 256;
    int m = 256;
    int k = 256;
    int i, j;
    int x, y, z;
    int nx, ny, nz;
    double c0=0.5;
    double c1=-0.25;
    int x0, x1, y0, y1, z0, z1;
    int left, right, up, down, z_down, z_up;
    int size, rank;
    int dims[3];
    int periods[3];
    int coords[3];
    int sizes[3];
    int sub_sizes[3];
    int subsize1[3];
    int subsize2[3];
    int starts[3];
    int sub_starts[3];
    double*** global_array;
    double*** split_old;
    double*** split_new;
    double*** new_global;
    int displs_x;
    int displs_y;
    int displs_z;

    MPI_Comm comm3d;
    MPI_Status status;
    MPI_Datatype xz_slice, yz_slice, xy_slice, xyz_slice, nxnynz_slice;
    MPI_Request sendreq[6], recvreq[6];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double t1, t2;
    t1 = MPI_Wtime();

    global_array = allocate_3d(n, m, k);
    init(global_array, n, m, k);

    // set number of procs in x, y & z domain, in total p*p*p processes
    int x_domain = 2;
    int y_domain = 2;
    int z_domain = 2;

    periods[0] = 0;
    periods[1] = 0;
    periods[2] = 0;

    dims[0] = x_domain;
    dims[1] = y_domain;
    dims[2] = z_domain;

    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &comm3d);
    MPI_Cart_get(comm3d, 3, dims, periods, coords);

    //get
    decompose(n, dims[0], coords[0], &x0, &x1);
    decompose(m, dims[1], coords[1], &y0, &y1);
    decompose(k, dims[2], coords[2], &z0, &z1);

    nx = x1 - x0 + 1;
    ny = y1 - y0 + 1;
    nz = z1 - z0 + 1;

    split_old = allocate_3d(nx+2, ny+2, nz+2);
    init2(split_old, nx+2, ny+2, nz+2);
    split_new = allocate_3d(nx+2, ny+2, nz+2);
    init2(split_new, nx+2, ny+2, nz+2);

    //create datatype for gathering each local array back to global array
    create_datatype_array(n, m, k, nx, ny, nz, x0, y0, z0, &xyz_slice);
    create_datatype_array(nx+2, ny+2, nz+2, nx, ny, nz, 1, 1, 1, &nxnynz_slice);

    /* Left/West and right/Est neigbors */
    MPI_Cart_shift(comm3d, 0, 1, &left, &right);

    /* Bottom/South and Upper/North neigbors */
    MPI_Cart_shift(comm3d, 1, 1, &down, &up);

    /* Zdown/South and Zup/North neigbors */
    MPI_Cart_shift(comm3d, 2, 1, &z_down , &z_up);

    //split global array onto p*p*p processes
    int zk = 1;
    int ym;
    int xn;
    int count = 0;
    for(z = z0; z <= z1; z++)
    {
        ym = 1;
        for(y = y0; y <= y1; y++)
        {
             xn = 1;
            for(x = x0; x <= x1; x++)
            {
                split_old[zk][ym][xn] = global_array[z][y][x];
                xn++;
            }
            ym++;
        }
        zk++;
    }

    nx = nx+2;
    ny = ny+2;
    nz = nz+2;

    //create datatype for communication
    create_datatype_array(nx, ny, nz, 1, ny-2, nz-2, 0, 0, 0, &yz_slice);
    create_datatype_array(nx, ny, nz, nx-2, 1, nz-2, 0, 0, 0, &xz_slice);
    create_datatype_array(nx, ny, nz, nx-2, nz-2, 1, 0, 0, 0, &xy_slice);

    int T=20;
    int nIters = 0;
    int t=0;

    while( t < T )
    {
      t++;

      #pragma omp parallel for shared(split_old, split_new, nz, ny, nx, c0, c1) private(z, y, x)
      for (z = 1; z < nz-1; z++)
      {
        for (y= 1; y < ny-1; y++)
        {
          for (x = 1; x < nx-1; x++)
          {
            split_new[z][y][x] = c0* split_old[z][y][x]  + c1 * (split_old[z][y][x-1] + split_old[z][y][x+1] +
              split_old[z][y-1][x] + split_old[z][y+1][x] +
              split_old[z-1][y][x] + split_old[z+1][y][x]);

          }
        }
      }

      MPI_Isend(&split_old[1][ny-2][1], 1, xz_slice, up, TAG1, comm3d, &sendreq[0]);
      MPI_Irecv(&split_old[1][0][1], 1, xz_slice, down, TAG1, comm3d, &recvreq[0]);
      // //
      //DOWN - UP, send XZ-plane
      MPI_Isend(&split_old[1][1][1], 1, xz_slice, down, TAG2, comm3d, &sendreq[1]);
      MPI_Irecv(&split_old[1][ny-1][1], 1, xz_slice, up, TAG2, comm3d, &recvreq[1]);
      // //
      //RIGHT - LEFT, send YZ-plane
      MPI_Isend(&split_old[1][1][nx-2], 1, yz_slice, right, TAG3, comm3d, &sendreq[2]);
      MPI_Irecv(&split_old[1][1][0], 1, yz_slice, left, TAG3, comm3d, &recvreq[2]);
      // //
      //LEFT - RIGHT, send YZ-plane
      MPI_Isend(&split_old[1][1][1], 1, yz_slice, left, TAG4, comm3d, &sendreq[3]);
      MPI_Irecv(&split_old[1][1][nx-1], 1, yz_slice, right, TAG4, comm3d, &recvreq[3]);
      // //
      //z-up - z-down, send XY-plane
      MPI_Isend(&split_old[nz-2][1][1], 1, xy_slice, z_up, TAG5, comm3d, &sendreq[4]);
      MPI_Irecv(&split_old[0][1][1], 1, xy_slice, z_down, TAG5, comm3d, &recvreq[4]);
      // //
      //z-down - z-up, send XY-plane
      MPI_Isend(&split_old[1][1][1], 1, xy_slice, z_down, TAG6, comm3d, &sendreq[5]);
      MPI_Irecv(&split_old[nz-1][1][1], 1, xy_slice, z_up, TAG6, comm3d, &recvreq[5]);

      MPI_Wait(&recvreq[0], MPI_STATUS_IGNORE);
      MPI_Wait(&recvreq[1], MPI_STATUS_IGNORE);
      MPI_Wait(&recvreq[2], MPI_STATUS_IGNORE);
      MPI_Wait(&recvreq[3], MPI_STATUS_IGNORE);
      MPI_Wait(&recvreq[4], MPI_STATUS_IGNORE);
      MPI_Wait(&recvreq[5], MPI_STATUS_IGNORE);

      double*** tmp;
      tmp = split_old;
      split_old = split_new;
      split_new = tmp;
      nIters = t;
    }

    nx = nx-2;
    ny = ny-2;
    nz = nz-2;

    //collect local arrays back to global array
    if(rank == 0)
    {
      new_global = allocate_3d(n, m, k);
      init2(new_global, n, m, k);
      int kz = 1;
      for(z = 0; z < nz; z++)
      {
          j = 1;
          for(y = 0; y < ny; y++)
          {
              i = 1;
              for(x = 0; x < nx; x++)
              {
                new_global[z][y][x] = split_old[kz][j][i];
                i++;
              }
              j++;
          }
          kz++;
      }

      for(i = 1; i < size; i++)
      {
        MPI_Recv(&displs_x, 1, MPI_INT, i, 111, MPI_COMM_WORLD, &status);
        MPI_Recv(&displs_y, 1, MPI_INT, i, 222, MPI_COMM_WORLD, &status);
        MPI_Recv(&displs_z, 1, MPI_INT, i, 333, MPI_COMM_WORLD, &status);
        MPI_Recv(&new_global[displs_z][displs_y][displs_x], 1, xyz_slice, i, 444, MPI_COMM_WORLD, &status);
      }
    }
    else if(rank > 0)
    {
      MPI_Send(&x0, 1, MPI_INT, 0, 111, MPI_COMM_WORLD);
      MPI_Send(&y0, 1, MPI_INT, 0, 222, MPI_COMM_WORLD);
      MPI_Send(&z0, 1, MPI_INT, 0, 333, MPI_COMM_WORLD);
      MPI_Send(&split_old[0][0][0], 1, nxnynz_slice, 0, 444, MPI_COMM_WORLD);
    }

    t2 = MPI_Wtime();
    printf( "Elapsed time is: %f, iters: %d \n", t2 - t1 , nIters);

    MPI_Type_free(&xyz_slice);
    MPI_Type_free(&nxnynz_slice);
    MPI_Type_free(&xz_slice);
    MPI_Type_free(&yz_slice);
    MPI_Type_free(&xy_slice);
    if(rank == 0)
      free3D(new_global);
    free3D(split_old);
    free3D(split_new);
    free3D(global_array);

    MPI_Finalize();

    return 0;
}
