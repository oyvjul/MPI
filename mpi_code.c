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

int main(int argc, char *argv[])
{

  int size, rank;
  int coords[3];
  int periods[3];
  int x, y, z, i, j, k;
  double ***Unew, ***Uold, ***global;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);




  MPI_Finalize();
  return 0;
}
