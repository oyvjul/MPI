#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

typedef struct
{
    int* elements;
    int* neighbours;
    double* Gx;
    double* Gy;
    double* Gz;

    double* nodes;
    double* centroid;
    double* area;
    double* Nx;
    double* Ny;
    double* Nz;
    double* Sx;
    double* Sy;
    double* Sz;

    double* tensor;
    double* volume;
    double totalVolume;
    int numtet;
    int numnodes;
} meshdata;

typedef struct
{
  double x, y, z;

} cube;

double ***dallocate_3d(int x, int y, int z)
{
  int i, j;
  double *storage = (double*)malloc(x * y * z * sizeof(*storage));
  double *alloc = storage;
  double ***matrix;
  matrix = (double***)malloc(z * sizeof(double**));

  for (i = 0; i < z; i++)
  {
    matrix[i] = (double**)malloc(y * sizeof(**matrix));

    for (j = 0; j < y; j++)
    {
      matrix[i][j] = alloc;
      alloc += x;
    }
  }

  return matrix;
}

void dinit_3d(double*** matrix, int x, int y, int z)
{
  int i, j, k;

  for(i = 0; i < z; i++)
  {
    for(j = 0; j < y; j++)
    {
      for(k = 0; k < x; k++)
      {
        matrix[i][j][k] = 0.0;
      }
    }
  }
}

void readmesh(char* infile, meshdata *M)
{
    char fname[200];
    char suffix[20];
    FILE *fpdata;
    int dump;
    double in;

    int err;

    strcpy(fname,infile);
    strcpy(suffix,".node");
    strcat(fname,suffix);
    fpdata = fopen(fname, "r");
    if(fpdata==NULL)
    {
        printf("\nFailure to open input file.\n");
        exit(0);
    }
    err = fscanf(fpdata, "%d", &M->numnodes);
    M->nodes = (double*)calloc(3*M->numnodes,sizeof(double));
    for(int j=0;j<3;++j)
        err = fscanf(fpdata, "%d", &dump);
    for (int i = 0; i < M->numnodes; i++)
    {
        err = fscanf(fpdata, "%d", &dump);
        for(int j=0;j<3;++j)
        {
            err = fscanf(fpdata, "%lf", &in);
            M->nodes[i*3+j]=in;///10+0.5;  //Scaling the mesh can be done here
        }
    }
    fclose(fpdata);

    strcpy(fname,infile);
    strcpy(suffix,".ele");
    strcat(fname,suffix);
    fpdata = fopen(fname, "r");
    if(fpdata==NULL)
    {
        printf("\nFailure to open input file.\n");
        exit(0);
    }
    err = fscanf(fpdata, "%d", &M->numtet);
    M->elements = (int*)calloc(4*M->numtet,sizeof(int));
    for(int j=0;j<2;++j)
        err = fscanf(fpdata, "%d", &dump);
    for (int i = 0; i < M->numtet; i++)
    {
        err = fscanf(fpdata, "%d", &dump);
        for(int j=0;j<4;++j)
            err = fscanf(fpdata, "%d", &M->elements[i*4+j]);
    }
    fclose(fpdata);

    strcpy(fname,infile);
    strcpy(suffix,".neigh");
    strcat(fname,suffix);
    fpdata = fopen(fname, "r");
    if(fpdata==NULL)
    {
        printf("\nFailure to open input file.\n");
        exit(0);
    }
    err = fscanf(fpdata, "%d", &dump);
    assert(dump==M->numtet);
    M->neighbours = (int*)calloc(4*M->numtet,sizeof(int));
    for(int j=0;j<1;++j)
        err = fscanf(fpdata, "%d", &dump);
    for (int i = 0; i < M->numtet; i++)
    {
        err = fscanf(fpdata, "%d", &dump);
        for(int j=0;j<4;++j)
            err = fscanf(fpdata, "%d", &M->neighbours[i*4+j]);
    }
    fclose(fpdata);
}

/*void computeVolumes(meshdata* M)
{

    //#pragma omp parallel for reduction(+ : M->totalVolume)
    //for(int i=0; i<M->numtet; i++)
    for(int i=0; i<20; i++)
    {
        //compute indices to access nodes array
        int A = M->elements[i*4]*3;
        int B = M->elements[i*4+1]*3;
        int C = M->elements[i*4+2]*3;
        int D = M->elements[i*4+3]*3;

    }
    //printf("Volume %d %lf   tot %lf \n",INSPECT,M->volume[INSPECT],M->totalVolume);
}*/

void compute_minmax(double *x_max, double *x_min, double *y_max, double *y_min, double *z_max, double *z_min, meshdata *m)
{
  int i;
  *x_max = 0;
  *x_min = 0;
  *y_max = 0;
  *y_min = 0;
  *z_max = 0;
  *z_min = 0;

  for(i = 0; i < m->numnodes; i++)
  {
    if(m->nodes[i*3] < *x_min)
      *x_min = m->nodes[i*3];

    if(m->nodes[i*3] > *x_max)
      *x_max = m->nodes[i*3];

    if(m->nodes[i*3+1] < *y_min)
      *y_min = m->nodes[i*3+1];

    if(m->nodes[i*3+1] > *y_max)
      *y_max = m->nodes[i*3+1];

    if(m->nodes[i*3+2] < *z_min)
      *z_min = m->nodes[i*3+2];

    if(m->nodes[i*3+2] > *z_max)
      *z_max = m->nodes[i*3+2];
  }
}

void calculate_centroid(meshdata *m)
{
  int i,j;
  m->centroid = (double*)calloc(3*m->numtet, sizeof(double));

  for(i = 0; i < m->numtet; i++)
  {
    int A = m->elements[i*4];
    int B = m->elements[i*4+1];
    int C = m->elements[i*4+2];
    int D = m->elements[i*4+3];

    for(j = 0; j < 3; j++)
    {
      m->centroid[i*3+j] = ((m->nodes[A*3+j]*m->nodes[B*3+j]*m->nodes[C*3+j]*m->nodes[D*3+j])*0.25);
    }

    //printf("%f \t %f \n ",m->nodes[i*3], m->nodes[i*3+1]);
  }
}

void determinant(double *nodes, int a, int b, int c, int d)
{
  const double constant = 1;
  double det_0, det_1, det_2, det_3, det_4;

  //det_0 =

}

int main(int argc, char *argv[])
{
  int i, j, k;
  meshdata m;
  double x_max;
  double x_min;
  double y_max;
  double y_min;
  double z_max;
  double z_min;
  double ***u_old, ***u_new;
  double det_0, det_1, det_2, det_3, det_4;
  const double constant = 1.0;
  int x = 10;
  int y = 10;
  int z = 10;
  double x_step, y_step, z_step;

  cube *grid, point;
  grid = calloc(x*y*z, sizeof(*grid));

  double a[3] = {1.1, 1.2, 1.3};
  double b[3] = {2.1, 2.2, 2.3};
  double c[3] = {3.1, 3.2, 3.3};
  double d[3] = {3.1, 3.2, 3.3};

  /*det_a = (a[0]*b[1]*c[2]*constant) + (a[0]*b[2]*constant*d[1]) + (a[0]*constant*c[1]*d[3])
        + (a[1]*b[0]*constant*d[3]) + (a[1]*b[2]*c[0]*constant) + (a[1]*constant*c[2]*d[0])
        + (a[2]*b[0]*c[1]*constant) + (a[2]*b[1]*constant*d[0]) + (a[2]*constant*c[0]*d[1])
        + (constant*b[0]*c[2]*d[1]) + (constant*b[1]*c[0]*d[2]) + (constant*b[2]*c[1]*d[0])
        //- (a[0]*b[1]*constant*d[2]) - (a[0]*b[2]*constant*d[1]) - (a[0]*constant*c[1]*d[3])
        - (a[1]*b[0]*constant*d[3]) - (a[1]*b[2]*c[0]*constant) - (a[1]*constant*c[2]*d[0])
        - (a[2]*b[0]*c[1]*constant) - (a[2]*b[1]*constant*d[0]) - (a[2]*constant*c[0]*d[1])
        - (constant*b[0]*c[2]*d[1]) - (constant*b[1]*c[0]*d[2]) - (constant*b[2]*c[1]*d[0]);*/



  readmesh("mesh_new/3Dheart.1", &m);
  compute_minmax(&x_max, &x_min, &y_max, &y_min, &z_max, &z_min, &m);

  x_step = (x_max - x_min)/(double)x;
  y_step = (y_max - y_min)/(double)y;
  z_step = (z_max - z_min)/(double)z;

  printf("x_max: %f \t x_min: %f \t y_max: %f \t y_min: %f \t z_max: %f \t z_min: %f\n", x_max, x_min, y_max, y_min, z_max, z_min);
  //printf("%f \n", step);
  for(i = 0; i < z; i++)
  {
    for(j = 0; j < y; j++)
    {
      for(k = 0; k < x; k++)
      {
        point.x = x_min + x_step*i;
        point.y = y_min + y_step*j;
        point.z = z_min + z_step*k;

        grid[i*y*z+j*y+k] = point;
      }
    }
  }

  printf("%f \n", grid[0].z);


  /*for(i = 0; i < 5; i++)
  {
    for(j = 0; j < 5; j++)
    {
      for(k = 1; k < 5; k++)
      {
        printf("%f ", arr[i][j][k]);
      }
      printf("\n");
    }
    printf("\n");
  }
  printf("\n");*/

  //calculate_centroid(&m);

}
