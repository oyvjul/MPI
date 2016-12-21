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
    matrix = malloc(y * sizeof(*matrix));
    
    for (i = 0; i < y; i++) 
    {
        matrix[i] = malloc(x * sizeof(**matrix));
        
        for (j = 0; j < x; j++) 
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
    
    for(k=0 ; k < K ; k++)
        for(i=0 ; i < M ; i++)
            for(j=0 ; j < N ; j++)
            {
                E[k][i][j]=1.0;
        
                if(i==0 || i == M-1 || j == 0 || j == N-1 || k==0 || k == K-1 )
                E[k][i][j]=0.0;
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

void decompose(int n, int dim, int rank, int* start, int* end)
{
    int length, rest;
    
    length = n/dim;
    rest = n%dim;
    *start = rank * length + (rank < rest ? rank : rest);
    *end = *start + length - (rank < rest ? 0 : 1);
    
    if((*end >= n) || (rank == dim-1))
        *end = n-2;
    
    if(*start == 0)
        *start = 1;
        
}

void processToMap(int me, int *xs, int *xe, int *ys, int *ye, int *zs, int *ze, int xcell, 
                  int ycell, int zcell, int x_domains, int y_domains, int z_domains)
{
  /* Index variables */
  int i, j, k, l, m, p, v;

  /* Computation of ys and ye with topoly process */
  for(i=1;i<=y_domains;i++) {
    ys[(i-1)*z_domains]=y_domains*(ycell+2)-ycell*i-2*(i-1);
    ye[(i-1)*z_domains]=ys[(i-1)*z_domains]+ycell-1;

    for(l=1;l<=z_domains-1;l++) {
      ys[(i-1)*z_domains+l]=ys[(i-1)*z_domains];
      ye[(i-1)*z_domains+l]=ys[(i-1)*z_domains+l]+ycell-1;
    }
  }

  /* Prolongation along x_domain */
  for(m=1;m<=y_domains;m++) {
    ys[(m-1)*z_domains]=y_domains*(ycell+2)-ycell*m-2*(m-1);
    ye[(m-1)*z_domains]=ys[(m-1)*z_domains]+ycell-1;

    for(i=1;i<=x_domains-1;i++) {
      ys[i*(y_domains*z_domains)+(m-1)*z_domains]=ys[(m-1)*z_domains];
      ye[i*(y_domains*z_domains)+(m-1)*z_domains]=ys[i*(y_domains*z_domains)+(m-1)*z_domains]+ycell-1;

      for(l=1;l<=z_domains-1;l++) {
        ys[i*(y_domains*z_domains)+(m-1)*z_domains+l]=ys[i*(y_domains*z_domains)+(m-1)*z_domains];
        ye[i*(y_domains*z_domains)+(m-1)*z_domains+l]=ys[i*(y_domains*z_domains)+(m-1)*z_domains+l]+ycell-1;
      }
    }
  }

  /* Computation of xs and xe with topoly process */
  for(i=0;i<=(z_domains*y_domains)-1;i++) {
    xs[i]=2;
    xe[i]=xs[i]+xcell-1;
  }

  for(j=1;j<=x_domains-1;j++) {
    for(k=0;k<=(z_domains*y_domains-1);k++) {    
      xs[j*(z_domains*y_domains)+k]=xs[(j-1)*(z_domains*y_domains)]+xcell+2;
      xe[j*(z_domains*y_domains)+k]=xs[j*(z_domains*y_domains)]+xcell-1;
    }
  }

  /* Computation of zs and ze with topoly process */
  for(k=0;k<=y_domains-1;k++) {
    v=k*z_domains;
    zs[v]=2;
    ze[v]=2+zcell-1;

    for(p=1;p<=x_domains-1;p++) {
      zs[v+p*(y_domains*z_domains)]=zs[v];
      ze[v+p*(y_domains*z_domains)]=ze[v];
    }
  }

  /* Prolongation along z_domain */
  for(m=1;m<=z_domains-1;m++) {
    for(i=0;i<=y_domains-1;i++) {
      l=m+i*z_domains;
      zs[l]=zs[l-1]+zcell+2;
      ze[l]=zs[l]+zcell-1;

      for(v=1;v<=x_domains-1;v++) {
        zs[l+v*(y_domains*z_domains)]=zs[l];
        ze[l+v*(y_domains*z_domains)]=zs[l+v*(y_domains*z_domains)]+zcell-1;
      }
    }
  }
}

int main (int argc, char* argv[])
{
    
    int n = 8;
    int m = 8;
    int k = 8;
    int x, y, z;
    double c0=0.5;
    double c1=-0.25;
    int x0, x1, y0, y1, z0, z1;
    int left, right, up, down, z_down, z_up;
    int total_x, total_y, total_z, x_cell, y_cell, z_cell;
    int size, rank;
    int dims[3];
    int periods[3];
    int coords[3];
    int sizes[3];
    int subsize1[3];
    int subsize2[3];
    int subsize3[3];
    int *xs, *ys, *zs, *xe, *ye, *ze;
    int starts[3];
    double*** Unew; 
    double*** Uold;
    
    MPI_Comm comm3d;
    MPI_Datatype matrix_xz, matrix_yz, matrix_xy;
    MPI_Status status;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    Unew= allocate_3d(n+2, m+2, k+2);
    Uold= allocate_3d(n+2, m+2, k+2);
    
    init(Unew, n+2, m+2, k+2);
    init(Uold, n+2, m+2, k+2);
    
    int x_domain = 2;
    int y_domain = 2;
    int z_domain = 2; 
    
    total_x = n+2*x_domain+2;
    total_y = m+2*y_domain+2;
    total_z = k+2*z_domain+2;  
    
    periods[0] = 0;
    periods[1] = 0;
    periods[2] = 0;  

    dims[0] = x_domain;
    dims[1] = y_domain;
    dims[2] = z_domain;  

    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &comm3d);
    MPI_Cart_get(comm3d, 3, dims, periods, coords);
    //MPI_Comm_rank(comm3d, &rank);
    
    /* Left/West and right/Est neigbors */
    MPI_Cart_shift(comm3d, 0, 1, &left, &right);

    /* Bottom/South and Upper/North neigbors */
    MPI_Cart_shift(comm3d, 1, 1, &down, &up);

    /* Zdown/South and Zup/North neigbors */
    MPI_Cart_shift(comm3d, 2, 1,&z_down ,&z_up);
    
    decompose(n, dims[0], coords[0], &x0, &x1);
    decompose(m, dims[1], coords[1], &y0, &y1);
    decompose(k, dims[2], coords[2], &z0, &z1);
    
    xs = malloc(size*sizeof(int));
    xe = malloc(size*sizeof(int));
    ys = malloc(size*sizeof(int));
    ye = malloc(size*sizeof(int));
    zs = malloc(size*sizeof(int));
    ze = malloc(size*sizeof(int));  
    
    x_cell = n/x_domain;
    y_cell = m/y_domain;
    z_cell = k/z_domain;  
    
    processToMap(rank, xs, xe, ys, ye, zs, ze, x_cell, y_cell, z_cell, x_domain, y_domain, z_domain);
    
    sizes[0] = total_x;
    sizes[1] = total_y;
    sizes[2] = total_z;

    
    starts[0]=0;
    starts[1]=0;
    starts[2]=0;
    
    /* Create matrix data type to communicate on vertical Oxz plan */
    subsize1[0]=x_cell;
    subsize1[1]=1;
    subsize1[2]=z_cell;

    MPI_Type_create_subarray(3, sizes, subsize1, starts, MPI_ORDER_C, MPI_DOUBLE, &matrix_xz);
    MPI_Type_commit(&matrix_xz);

    /* Create matrix data type to communicate on vertical Oyz plan */
    subsize2[0]=1;
    subsize2[1]=y_cell;
    subsize2[2]=z_cell;

    MPI_Type_create_subarray(3, sizes, subsize2, starts, MPI_ORDER_C, MPI_DOUBLE, &matrix_yz);
    MPI_Type_commit(&matrix_yz);

    /* Create matrix data type to communicate on vertical Oxy plan */
    subsize3[0]=x_cell;
    subsize3[1]=y_cell;
    subsize3[2]=1;

    MPI_Type_create_subarray(3, sizes, subsize3, starts, MPI_ORDER_C, MPI_DOUBLE, &matrix_xy); 
    MPI_Type_commit(&matrix_xy);

    
    printf("rank: %d    starts: (%d, %d, %d)    ends: (%d, %d, %d)\n", rank, x0, y0, z0, x1, y1, z1);
    MPI_Barrier(comm3d);
    usleep(1000);
    //printf("\nrank: %d    starts: (%d, %d, %d)    ends: (%d, %d, %d)\n", rank, xs[rank], ys[rank], zs[rank], xe[rank], ye[rank], ze[rank]);
    
    int T=20;
    
    int nIters = 0;
    
    int t=0;
        
    while( t < T )
    {
        t++;
        
        for (z = z0; z <= z1; z++)
        {
            for (y= y0; y <= y1; y++)
            {
                for (x = x0; x <= x1; x++) 
                {
                    Unew[z][y][x] = c0* Uold[z][y][x]  + c1 * (Uold[z][y][x-1] + Uold[z][y][x+1] +
                        Uold[z][y-1][x] + Uold[z][y+1][x] +
                        Uold[z-1][y][x] + Uold[z+1][y][x]); 
                }
            }
        }
        
        //MPI_Sendrecv(&Unew[x0][y0][z0], 1, matrix_xz , up, TAG1, 
                     //&Unew[x0][y1][z0], 1, matrix_xz, down , TAG1, comm3d, &status);
        
        //MPI_Sendrecv(&Unew[x0+1][y1+1][z0+1], 1, matrix_xz, down, TAG1, 
                     //&Unew[x0+1][y0+1][z0+1], 1, matrix_xz, up , TAG1, comm3d, &status);

        //MPI_Sendrecv(&Unew[x1+1][y0+1][z0+1], 1, matrix_yz, right, TAG2, 
                     //&Unew[x0+1][y0+1][z0+1], 1, matrix_yz, left, TAG2, comm3d, &status);

        //MPI_Sendrecv(&Unew[x0+1][y0+1][z0+1], 1, matrix_yz, left, TAG2, 
                     //&Unew[x1+1][y0+1][z0+1], 1, matrix_yz, right, TAG2, comm3d, &status);

        //MPI_Sendrecv(&Unew[x0+1][y0+1][z1+1], 1, matrix_xy, z_up, TAG3, 
                     //&Unew[x0+1][y0+1][z0+1], 1, matrix_xy, z_down, TAG3, comm3d, &status);

        //MPI_Sendrecv(&Unew[x0+1][y0+1][z0+1], 1, matrix_xy, z_down, TAG3, 
                     //&Unew[x0+1][y0+1][z1+1], 1, matrix_xy, z_up, TAG3, comm3d, &status);

        double*** tmp;
        tmp = Uold; 
        Uold = Unew; 
        Unew = tmp;
        nIters = t; 
    }

    int i;
    int j;
    
    /*for(z=0 ; z < k ; z++)
    {
        for(y=0 ; y < m ; y++)
        {
            for(x=0 ; x < n ; x++)
            {
                printf("%f ", Unew[z][y][x]);
            }
            printf("\n");
        }
    }*/
    //calculatel2Norm(Uold, n, m, k, T);
    
    free3D(Uold);
    free3D(Unew);
    
    MPI_Finalize();
    
    return 0;
}
