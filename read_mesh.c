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

double determinant(double *a, double *b, double *c, double *d)
{
  const double constant = 1.0;
  double determinant;

  determinant = (a[0]*b[1]*c[2]*constant) + (a[0]*b[2]*constant*d[1]) + (a[0]*constant*c[1]*d[2])
        + (a[1]*b[0]*constant*d[2]) + (a[1]*b[2]*c[0]*constant) + (a[1]*constant*c[2]*d[0])
        + (a[2]*b[0]*c[1]*constant) + (a[2]*b[1]*constant*d[0]) + (a[2]*constant*c[0]*d[1])
        + (constant*b[0]*c[2]*d[1]) + (constant*b[1]*c[0]*d[2]) + (constant*b[2]*c[1]*d[0])
        - (a[0]*b[1]*constant*d[2]) - (a[0]*b[2]*c[1]*constant) - (a[0]*constant*c[2]*d[1])
        - (a[1]*b[0]*c[2]*constant) - (a[1]*b[2]*constant*d[0]) - (a[1]*constant*c[0]*d[2])
        - (a[2]*b[0]*constant*d[1]) - (a[2]*b[1]*c[0]*constant) - (a[2]*constant*c[1]*d[0])
        - (constant*b[0]*c[1]*d[2]) - (constant*b[1]*c[2]*d[0]) - (constant*b[2]*c[0]*d[1]);

  return determinant;
}

int inside(meshdata *m, double point_x, double point_y, double point_z)
{
  int i, j;
  double a_vector[3], b_vector[3], c_vector[3], d_vector[3];
  int a, b, c, d;
  double point_vector[3];
  double det_0, det_1, det_2, det_3, det_4;

  point_vector[0] = point_x;
  point_vector[1] = point_y;
  point_vector[2] = point_z;

  for(i = 0; i < m->numtet; i++)
  {
    a = m->elements[i*4];
    b = m->elements[i*4+1];
    c = m->elements[i*4+2];
    d = m->elements[i*4+3];

    for(j = 0; j < 3; j++)
    {
      a_vector[j] = m->nodes[a*3+j];
      b_vector[j] = m->nodes[b*3+j];
      c_vector[j] = m->nodes[c*3+j];
      d_vector[j] = m->nodes[d*3+j];
    }

    det_0 = determinant(a_vector, b_vector, c_vector, d_vector);
    det_1 = determinant(point_vector, b_vector, c_vector, d_vector);
    det_2 = determinant(a_vector, point_vector, c_vector, d_vector);
    det_3 = determinant(a_vector, b_vector, point_vector, d_vector);
    det_4 = determinant(a_vector, b_vector, c_vector, point_vector);

    if((det_0 == (det_1+det_2+det_3+det_4)) || (det_1==0 && det_2==0 && det_3==0 && det_4==0))
      return 1;
    else
      return 0;
  }

  return 0;
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
  int x = 100;
  int y = 100;
  int z = 100;
  int is_inside;
  double x_step, y_step, z_step;
  double det_a;
  double *grid_x, *grid_y, *grid_z;
  grid_x = (double*)calloc(x+2, sizeof(double));
  grid_y = (double*)calloc(y+2, sizeof(double));
  grid_z = (double*)calloc(z+2, sizeof(double));
  u_new = dallocate_3d(x+2, y+2, z+2);
  u_old = dallocate_3d(x+2, y+2, z+2);
  dinit_3d(u_new, x+2, y+2, z+2);
  dinit_3d(u_old, x+2, y+2, z+2);

  cube *grid, point;
  //grid = calloc(x*y*z, sizeof(*grid));

  double a[3] = {3, 3, 3};
  double b[3] = {0, 1, 0};
  double c[3] = {0, 0, 1};
  double d[3] = {0, 0, 0};

  /*det_a = (a[0]*b[1]*c[2]*constant) + (a[0]*b[2]*constant*d[1]) + (a[0]*constant*c[1]*d[2])
        + (a[1]*b[0]*constant*d[2]) + (a[1]*b[2]*c[0]*constant) + (a[1]*constant*c[2]*d[0])
        + (a[2]*b[0]*c[1]*constant) + (a[2]*b[1]*constant*d[0]) + (a[2]*constant*c[0]*d[1])
        + (constant*b[0]*c[2]*d[1]) + (constant*b[1]*c[0]*d[2]) + (constant*b[2]*c[1]*d[0])
        - (a[0]*b[1]*constant*d[2]) - (a[0]*b[2]*c[1]*constant) - (a[0]*constant*c[2]*d[1])
        - (a[1]*b[0]*c[2]*constant) - (a[1]*b[2]*constant*d[0]) - (a[1]*constant*c[0]*d[2])
        - (a[2]*b[0]*constant*d[1]) - (a[2]*b[1]*c[0]*constant) - (a[2]*constant*c[1]*d[0])
        - (constant*b[0]*c[1]*d[2]) - (constant*b[1]*c[2]*d[0]) - (constant*b[2]*c[0]*d[1]);*/



  readmesh("mesh_new/3Dheart.1", &m);
  compute_minmax(&x_max, &x_min, &y_max, &y_min, &z_max, &z_min, &m);

  x_step = (x_max - x_min)/(double)x;
  y_step = (y_max - y_min)/(double)y;
  z_step = (z_max - z_min)/(double)z;

  printf("x_max: %f \t x_min: %f \t y_max: %f \t y_min: %f \t z_max: %f \t z_min: %f\n", x_max, x_min, y_max, y_min, z_max, z_min);
  printf("%0.30f \n", det_a);

  //for(i = 0; i <= x; i++)
    //printf("%f \n", x_min + x_step*i);


  /*for(i = 0; i <= z; i++)
  {
    for(j = 0; j <= y; j++)
    {
      for(k = 0; k <= x; k++)
      {
        point.x = x_min + x_step*k;
        point.y = y_min + y_step*j;
        point.z = z_min + z_step*i;

        grid[i*x*y+j*x+k] = point;
        /*grid_x[i] = x_min + x_step*i;
        grid_y[i] = y_min + y_step*i;
        grid_z[i] = z_min + z_step*i;*/
      /*}
    }
  }*/

  for(i = 0; i <= z; i++)
  {
    grid_x[i] = x_min + x_step*i;
    grid_y[i] = y_min + y_step*i;
    grid_z[i] = z_min + z_step*i;
  }

  int count_inside = 0;
  int count_outside = 0;

  for(i = 0; i <= z; i++)
  {
    for(j = 0; j <= y; j++)
    {
      for(k = 0; k <= x; k++)
      {
        //is_inside = inside(&m, grid[i*x*y+j*x+k].x, grid[i*x*y+j*x+k].y, grid[i*x*y+j*x+k].z);
        is_inside = inside(&m, grid_x[k], grid_y[j], grid_z[i]);

        if(is_inside == 1)
        {
          u_old[i][j][k] = 1;
          count_inside++;
        }
        else
        {
          u_old[i][j][k] = 0;
          count_outside++;
        }
      }
    }
  }

  printf("total inside points: %d, total outside points: %d , num iters: %d\n", count_inside, count_outside, x*y*z);

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
