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
#include <omp.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include <time.h>

#define chunk 64
const double kMicro = 1.0e-6;

double ***alloc3D(int n, int m,int k)
{
    double ***m_buffer=NULL;

    int nx=n, ny=m, nk = k;

    m_buffer = (double***)malloc(sizeof(double**)* nk);
    assert(m_buffer);

    double** m_tempzy = (double**)malloc(sizeof(double*)* nk * ny);
    double *m_tempzyx = (double*)malloc(sizeof(double)* nx * ny * nk );

    for ( int z = 0 ; z < nk ; z++, m_tempzy += ny ) {
        m_buffer[z] = m_tempzy;
        for ( int y = 0 ; y < ny ; y++, m_tempzyx += nx ) {
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
                E[k][i][j]=2.0;

                //if(i==0 || i == M-1 || j == 0 || j == N-1 || k==0 || k == K-1 )
                //E[k][i][j]=0.0;
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

int main (int argc, char* argv[])
{

    int n = 6;
    int m = 6;
    int k = 6;
    double result = 0.0;
    int sum_iters = 0;
    int x, y, z;
    clock_t tt;
    tt = clock();

    double c0=0.5;
    double c1=-0.25;

    double*** Unew; double*** Uold;

    Unew= allocate_3d(n+2, m+2, k+2);
    Uold= allocate_3d(n+2, m+2, k+2);

    init(Unew, n+2, m+2, k+2);
    init(Uold, n+2, m+2, k+2);

    int T=20;

    int nIters = 0;
    double test = 0;
    //printf("test\n");

    int t=0;

    while( t < T )
    {
        t++;

        //7-point stencil
        for (z=1; z<= k; z++)
        {
            for (y=1; y<= m; y++)
            {
                for (x=1; x<= n; x++)
                {
                    Unew[z][y][x] = c0* Uold[z][y][x]  + c1 * (Uold[z][y][x-1] + Uold[z][y][x+1] +
                        Uold[z][y-1][x] + Uold[z][y+1][x] +
                        Uold[z-1][y][x] + Uold[z+1][y][x]);

                    result += Unew[z][y][x];
                    sum_iters++;
                }
            }
        }

        test += result;

        double*** tmp;
        tmp = Uold;
        Uold = Unew;
        Unew = tmp;
        nIters = t;
    }

    /*for (z=1; z<= k; z++)
    {
        for (y=1; y<= m; y++)
        {
            for (x=1; x<= n; x++)
            {
              printf("%f", Uold[z][y][x]);
            }
            printf("\n");
        }
        printf("\n");
    }*/

    int i;
    int j;
    printf("num iters: %d result: %0.1f \n", sum_iters, test);
    //calculatel2Norm(Uold, n, m, k, T);
    tt = clock() - tt;
    double time_taken = ((double)tt)/CLOCKS_PER_SEC; // in seconds

    printf("fun() took %f seconds to execute \n", time_taken);
    free3D(Uold);
    free3D(Unew);

    return 0;
}
