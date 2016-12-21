/*
 by Didem Unat
 3D 7-point jacobi
 Written to be used as an input program to mint translator

 See the alloc2D function, which allocates contiguous memory space to
 the array.
 */

//#include "common.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include "mpi.h"
#include <unistd.h>

#define chunk 64
#define TAG1 1
#define TAG2 2
#define TAG3 3
#define TAG4 4
#define TAG5 5
#define TAG6 6
const double kMicro = 1.0e-6;

double ***alloc3D(int n, int m,int k)
{
    double ***m_buffer=NULL;

    int nx=n, ny=m, nk = k;

    m_buffer = (double***)malloc(sizeof(double**)* nk);
    assert(m_buffer);

    double **m_tempzy = (double**)malloc(sizeof(double*)* nk * ny);
    double *m_tempzyx = (double*)malloc(sizeof(double)* nx * ny * nk );

    for ( int z = 0 ; z < nk ; z++, m_tempzy += ny )
    {
        m_buffer[z] = m_tempzy;

        for ( int y = 0 ; y < ny ; y++, m_tempzyx += nx )
        {
            m_buffer[z][y] = m_tempzyx;
        }
    }

    return m_buffer;
}

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

double getTime()
{
    struct timeval TV;

    const int RC = gettimeofday(&TV, NULL);
    if(RC == -1)
    {
        printf("ERROR: Bad call to gettimeofday\n");
        return(-1);
    }

    return( ((double)TV.tv_sec) + kMicro * ((double)TV.tv_usec) );

}  // end getTime()

//allocate 3D array
double ***alloc3D_(int n, int m,int k){

    double ***E=NULL;

    int nx=n, ny=m, nk = k;

    E = (double***)malloc(sizeof(double**)* nk);
    assert(E);

    E[0] = (double**)malloc(sizeof(double*)* nk * ny);
    E[0][0] = (double*)malloc(sizeof(double)*nx * ny * nk );

    int jj,kk;

    for(kk=0 ; kk < nk ; kk++){

        if(kk > 0)
        {
            E[kk] = E[kk-1] + ny ;
            E[kk][0] = E[kk-1][0] + ny*nx ;
        }

        for(jj=1; jj< ny; jj++) {
            E[kk][jj] = E[kk][jj-1] + nx ;
        }
    }
    return(E);
}

void free3D(double*** E)
{
    //int k=0;
    /*  for(k=0 ; k < m ; k++)
     {
     free(E[k]);
     }*/
    free(E[0][0]);
    free(E[0]);
    free(E);

}
void init(double*** E, int N, int M, int K)
{
    int i,j,k;
    double count = 0.0;

    for(k=0 ; k < K ; k++)
        for(i=0 ; i < M ; i++)
            for(j=0 ; j < N ; j++)
            {
                count++;
                E[j][i][k] = count;

                //if(i==0 || i == M-1 || j == 0 || j == N-1 || k==0 || k == K-1 )
                //E[k][i][j] = 0.0;
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
                E[j][i][k] = 0;
                if(i==0 || j == 0 ||  k==0)
                  E[j][i][k] = 0.0;
            }
}


//calculate l2norm for comparison
void calculatel2Norm(double*** E, int N, int M, int K, int nIters)
{
    int i, j, k  =0;

    float mx = -1;
    float l2norm = 0;

    for (k=1; k<= K ; k++){
        for (j=1; j<= M; j++){
            for (i=1; i<= N; i++) {
                l2norm += E[k][j][i]*E[k][j][i];

                if (E[k][j][i] > mx)
                mx = E[k][j][i];
            }
        }
    }
    l2norm /= (float) ((N)*(M)*(K));
    l2norm = sqrt(l2norm);
    printf(":N %d M %d K %d , iteration %d\n", N, M, K , nIters);
    printf(":max: %20.12e, l2norm: %20.12e\n",mx,l2norm);
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

    //if(*start == 0)
        //*start = 1;

}

int main (int argc, char* argv[])
{

    int n = 256;
    int m = 256;
    int k = 256;
    int idx;
    int i, j, d;
    int x, y, z;
    int nx, ny, nz;
    double c0=0.5;
    double c1=-0.25;
    int x0, x1, y0, y1, z0, z1;
    int left, right, up, down, z_down, z_up;
    int total_x, total_y, total_z, x_cell, y_cell, z_cell;
    int size, rank;      /* holds initial grid block specification */
    int dims[3];
    int periods[3];
    int coords[3];
    int sizes[3];
    int sub_sizes[3];
    int subsize1[3];
    int subsize2[3];
    int subsize3[3];
    int *xs, *ys, *zs, *xe, *ye, *ze;
    int starts[3];
    int sub_starts[3];
    double*** global_array;
    double*** Uold;
    double*** split_old;
    double*** split_new;
    double*** new_global;
    double* sendbuf = NULL;
    int* recvcount = NULL;
    int* displs = NULL;
    int displs_x;
    int displs_y;
    int displs_z;

    MPI_Comm comm3d;
    MPI_Status status;
    MPI_Datatype xz_slice, yz_slice, xy_slice, xyz_slice, nxnynz_slice;
    MPI_Request sendreq[6], recvreq[6];
    //MPI_Request *sendreq, *recvreq;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double t1, t2;
    //t1 = MPI_Wtime();

    global_array = allocate_3d(n, m, k);
    init(global_array, n, m, k);

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

    sizes[0] = n;
    sizes[1] = m;
    sizes[2] = k;

    subsize1[0] = nx;
    subsize1[1] = ny;
    subsize1[2] = nz;

    // HER ER PROBLEMET
    starts[0] = x0;
    starts[1] = y0;
    starts[2] = z0;

    MPI_Type_create_subarray(3,  sizes, subsize1, starts, MPI_ORDER_C, MPI_DOUBLE, &xyz_slice);
    MPI_Type_commit(&xyz_slice);

    sub_sizes[0] = nx+2;
    sub_sizes[1] = ny+2;
    sub_sizes[2] = nz+2;

    subsize2[0] = nx;
    subsize2[1] = ny;
    subsize2[2] = nz;

    // HER ER PROBLEMET
    sub_starts[0] = 1;
    sub_starts[1] = 1;
    sub_starts[2] = 1;

    MPI_Type_create_subarray(3,  sub_sizes, subsize2, sub_starts, MPI_ORDER_C, MPI_DOUBLE, &nxnynz_slice);
    MPI_Type_commit(&nxnynz_slice);
    //MPI_Comm_rank(comm3d, &rank);

    /* Left/West and right/Est neigbors */
    MPI_Cart_shift(comm3d, 0, 1, &left, &right);

    /* Bottom/South and Upper/North neigbors */
    MPI_Cart_shift(comm3d, 1, 1, &down, &up);

    /* Zdown/South and Zup/North neigbors */
    MPI_Cart_shift(comm3d, 2, 1, &z_down , &z_up);

    //printf("rank: %d coords: (%d %d %d)   (%d %d %d)   (%d %d %d)  (%d %d %d)\n",
           //rank, coords[0], coords[1], coords[2], x0, y0, z0, x1, y1, z1, nx, ny, nz);

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
                //printf("%0.1f ", global_array[z][y][x]);
                split_old[xn][ym][zk] = global_array[x][y][z];
                //split[xn][ym][zk] = rank+1;
                //printf("%0.1f ", split[zk][ym][xn]);
                xn++;
            }
            ym++;
        }
        zk++;
    }

    /*if(rank == 0)
    {
      new = (double***) realloc(split, (n*m*k)*sizeof(double**));
    }*/

    /*recvcount = (int*)malloc(size*sizeof(int));
    displs = (int*)malloc(size*sizeof(int));

    for(i = 0; i < size; i++)
    {
      recvcount[i] = 1;
      displs[i] = 0;
    }*/

    //if(rank < 7)
    //MPI_Gatherv(&split[0][0][0], nx*ny*nz, MPI_DOUBLE, &new[x0][y0][z0], recvcount, &displs[rank],  xyz_slice, 0, MPI_COMM_WORLD);

    //MPI_Gather(&split[0][0][0], nx*ny*nz, MPI_DOUBLE, &new[x0][y0][z0], 0, xyz_slice, 0, MPI_COMM_WORLD);
    nx = nx+2;
    ny = ny+2;
    nz = nz+2;

    MPI_Type_vector(nx-2, nx, (ny*nz), MPI_DOUBLE, &xz_slice);
    MPI_Type_commit(&xz_slice);

    MPI_Type_vector(ny-2, ny, nz, MPI_DOUBLE, &yz_slice);
    MPI_Type_commit(&yz_slice);

    MPI_Type_vector(((nz-1)*(nz-1))-1, 1, nz,  MPI_DOUBLE, &xy_slice);
    MPI_Type_commit(&xy_slice);

    //UP - DOWN, send XZ-plane



    //MPI_Isend(&split[0][0][0], 1, xz_slice, down, TAG2, comm3d, &sendreq[1]);
    //MPI_Irecv(&split[0][ny-1][0], 1, xz_slice, up, TAG2, comm3d, &recvreq[1]);

    //RIKTIG DOWN - UP, send XZ-plane
    //MPI_Sendrecv(&split[1][1][1], 1, xz_slice, down,   TAG1,
                 //&split[1][ny-1][1],    1, xz_slice, up, TAG1, comm3d, MPI_STATUS_IGNORE);

    //RIKTIG UP - DOWN, send XZ-plane
    //MPI_Sendrecv(&split[0][ny-1][0], 1, xz_slice, up,   TAG1,
                 //&split[0][0][0],    1, xz_slice, down, TAG1, comm3d, MPI_STATUS_IGNORE);

    //RIKTIG LEFT - RIGHT, send yz-plane
    //MPI_Sendrecv(&split[1][1][0], 1, yz_slice, left,   TAG2,
                 //&split[nx-1][1][0],    1, yz_slice, right, TAG2, comm3d, MPI_STATUS_IGNORE);
    //RIKTIG RIGHT - LEFT, send YZ-plane
    //MPI_Sendrecv(&split[nx-1][0][0], 1, yz_slice, right,   TAG2,
                 //&split[0][0][0],    1, yz_slice, left, TAG2, comm3d, MPI_STATUS_IGNORE);

    //RIKTIG Z-down - Z-UP, send xy-plane
    //MPI_Sendrecv(&split[1][1][1], 1, xy_slice, z_down,   TAG3,
                 //&split[1][1][nz-1],    1, xy_slice, z_up, TAG3, comm3d, MPI_STATUS_IGNORE);
                 //RIKTIG Z-down - Z-UP, send xy-plane
    //RIKTIG Z-up - Z-down, send XY-plane
    //MPI_Sendrecv(&split[0][0][nz-1], 1, xy_slice, z_up,   TAG3,
                 //&split[0][0][0],    1, xy_slice, z_down, TAG3, comm3d, MPI_STATUS_IGNORE);

    int T=20;

    int nIters = 0;

    int t=0;

    t1 = MPI_Wtime();
    while( t < T )
    {
        t++;

        MPI_Isend(&split_old[1][ny-2][1], 1, xz_slice, up, TAG1, comm3d, &sendreq[0]);
        MPI_Irecv(&split_old[1][0][1], 1, xz_slice, down, TAG1, comm3d, &recvreq[0]);

        //DOWN - UP, send XZ-plane
        MPI_Isend(&split_old[1][1][1], 1, xz_slice, down, TAG2, comm3d, &sendreq[1]);
        MPI_Irecv(&split_old[1][ny-1][1], 1, xz_slice, up, TAG2, comm3d, &recvreq[1]);

        //RIGHT - LEFT, send YZ-plane
        MPI_Isend(&split_old[nx-2][1][0], 1, yz_slice, right, TAG3, comm3d, &sendreq[2]);
        MPI_Irecv(&split_old[0][1][0], 1, yz_slice, left, TAG3, comm3d, &recvreq[2]);

        //LEFT - RIGHT, send YZ-plane
        MPI_Isend(&split_old[1][1][0], 1, yz_slice, left, TAG4, comm3d, &sendreq[3]);
        MPI_Irecv(&split_old[nx-1][1][0], 1, yz_slice, right, TAG4, comm3d, &recvreq[3]);

        //z-up - z-down, send XY-plane
        MPI_Isend(&split_old[1][1][nz-2], 1, xy_slice, z_up, TAG5, comm3d, &sendreq[4]);
        MPI_Irecv(&split_old[1][1][0], 1, xy_slice, z_down, TAG5, comm3d, &recvreq[4]);

        //z-down - z-up, send XY-plane
        MPI_Isend(&split_old[1][1][1], 1, xy_slice, z_down, TAG6, comm3d, &sendreq[5]);
        MPI_Irecv(&split_old[1][1][nz-1], 1, xy_slice, z_up, TAG6, comm3d, &recvreq[5]);

        for (z = 1; z < nz-1; z++)
        {
            for (y= 1; y < ny-1; y++)
            {
                for (x = 1; x < nx-1; x++)
                {
                    split_new[x][y][z] = c0* split_old[x][y][z]  + c1 * (split_old[x-1][y][z] + split_old[x+1][y][z] +
                        split_old[x][y-1][z] + split_old[x][y+1][z] +
                        split_old[x-1][y][z] + split_old[x+1][y][z]);
                }
            }
        }

        double*** tmp;
        tmp = split_old;
        split_old = split_new;
        split_new = tmp;
        nIters = t;
    }

    MPI_Wait(&recvreq[0], MPI_STATUS_IGNORE);
    MPI_Wait(&recvreq[1], MPI_STATUS_IGNORE);
    MPI_Wait(&recvreq[2], MPI_STATUS_IGNORE);
    MPI_Wait(&recvreq[3], MPI_STATUS_IGNORE);
    MPI_Wait(&recvreq[4], MPI_STATUS_IGNORE);
    MPI_Wait(&recvreq[5], MPI_STATUS_IGNORE);

    nx = nx-2;
    ny = ny-2;
    nz = nz-2;

    if(rank == 0)
    {
      new_global = allocate_3d(n, m, k);
      init2(new_global, n, m, k);
      int count = 0;
      int kz = 1;
      for(z = 0; z < nz; z++)
      {
          j = 1;
          for(y = 0; y < ny; y++)
          {
              i = 1;
              for(x = 0; x < nx; x++)
              {
                new_global[x][y][z] = split_old[i][j][kz];
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
        MPI_Recv(&new_global[displs_x][displs_y][displs_z], 1, xyz_slice, i, 444, MPI_COMM_WORLD, &status);
      }

      /*for (z=0; z < k; z++)
      {
          for (y=0; y < m; y++)
          {
              for (x=0; x < n; x++)
              {
                  printf("%0.1f\t", new[x][y][z]);
              }
              printf("\n");
          }
          printf("\n");
      }
      printf("\n  ");*/
    }
    else if(rank > 0)
    {
      MPI_Send(&x0, 1, MPI_INT, 0, 111, MPI_COMM_WORLD);
      MPI_Send(&y0, 1, MPI_INT, 0, 222, MPI_COMM_WORLD);
      MPI_Send(&z0, 1, MPI_INT, 0, 333, MPI_COMM_WORLD);
      MPI_Send(&split_old[0][0][0], 1, nxnynz_slice, 0, 444, MPI_COMM_WORLD);
    }

    t2 = MPI_Wtime();
    printf( "Elapsed time is %f\n", t2 - t1 );
    //calculatel2Norm(Uold, n, m, k, T);
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
