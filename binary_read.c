#include <stdio.h>
#include <stdlib.h>

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

void allocate_matrix_part_from_read_matrix_method(double*** matrix, int* num_rows, int* num_cols)
{
    int i;
    /* storage allocation of the matrix */
    *matrix = (double**)malloc((*num_rows)*sizeof(double*));
    (*matrix)[0] = (double*)malloc((*num_rows)*(*num_cols)*sizeof(double));
    for (i=1; i<(*num_rows); i++)
        (*matrix)[i] = (*matrix)[i-1]+(*num_cols);
        /* read in the entire matrix */
}

void read_matrix_binaryformat(char* filename, double ****matrix, int *x, int *y, int *z)
{
    int i;
    FILE* fp = fopen (filename,"rb");
    fread (x, sizeof(int), 1, fp);
    fread (y, sizeof(int), 1, fp);
    fread (z, sizeof(int), 1, fp);
    /* storage allocation of the matrix */
        /* read in the entire matrix */
    fread ((*matrix)[0][0], sizeof(double), (*x)*(*y)*(*z), fp);
    fclose (fp);
}

void write_matrix_binaryformat(char* filename, double ***matrix, int x, int y, int z)
{
    FILE *fp = fopen (filename,"wb");
    fwrite (&x, sizeof(int), 1, fp);
    fwrite (&y, sizeof(int), 1, fp);
    fwrite (&z, sizeof(int), 1, fp);
    fwrite (matrix[0][0], sizeof(double), x*y*z, fp);
    fclose (fp);
}

void read_matrix(char* filename, double ***m1, int x, int y, int z)
{
    int i, j, k;
    FILE* fp = fopen (filename,"r");

    for(i = 1; i <= z+1; i++)
    {
      for(j = 1; j <= y+1; j++)
      {
        for(k = 1; k <= x+1; k++)
        {
          //m1[i][j][k] = x_min + x_step*(k-1);

          fscanf(fp, "%lf", &m1[i][j][k]);
        }
        fscanf(fp, "\n");
      }
      fscanf(fp, "\n");
    }

    fclose (fp);
}

void write_matrix(char* filename, double ***m1, int x, int y, int z)
{
  int i, j, k;
  FILE *fp = fopen (filename,"w");
  double x_max = 5.0;
  double x_min = 1.0;

  double x_step = (x_max - x_min)/(double)x;

  for(i = 1; i <= z+1; i++)
  {
    for(j = 1; j <= y+1; j++)
    {
      for(k = 1; k <= x+1; k++)
      {
        m1[i][j][k] = x_min + x_step*(k-1);

        fprintf(fp, "%0.20f ", m1[i][j][k]);
      }
      fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
  }

  fclose (fp);
}

int main(int argc, char *argv[])
{
  int i, j, k, ii, jj, kk;
  double count = 1;
  int x = 5;
  int y = 5;
  int z = 5;

  double ***m1, ***m2;

  FILE *file, *file_read;
  double **m, **m_new;
  double *u =(double*)malloc(10*sizeof(double));
  double *u_new =(double*)malloc(10*sizeof(double));

  allocate_matrix_part_from_read_matrix_method(&m, &x, &y);
  m1 = dallocate_3d(x+3, y+3, z+3);
  m2 = dallocate_3d(x+3, y+3, z+3);

  double x_max = 5.0;
  double x_min = 1.0;

  double x_step = (x_max - x_min)/(double)x;

  /*for(i = 1; i <= z+1; i++)
  {
    for(j = 1; j <= y+1; j++)
    {
      for(k = 1; k <= x+1; k++)
      {
        m1[i][j][k] = x_min + x_step*(k-1);
      }
    }
  }*/

  //write_matrix_binaryformat("test.bin", m1, x+3, y+3, z+3);
  write_matrix("heart.txt", m1, x, y, z);

  read_matrix("heart.txt", m2, x, y, z);

  //read_matrix_binaryformat("test.bin", &m2, &x, &y, &z);

  for(i = 1; i <= z+1; i++)
  {
    for(j = 1; j <= y+1; j++)
    {
      for(k = 1; k <= x+1; k++)
      {
        printf("%f ", m1[i][j][k]);
      }
      printf("\n");
    }
    printf("\n");
  }


  return 0;
}
