#include <stdio.h>
#include <mpi.h>
#include "common.h"
#include "math.h"
#include "time.h"

void calculatel2Norm(double*** E, int d, int r, int iters)
{
  int i, j, k;

  float mx = -1;
  double l2norm = 0;

  for (i = 1; i <= d; i++)
  {
      for (j = 1; j <= d; j++)
      {
          for (k = 1; k <= d; k++)
          {
            l2norm += E[i][j][k]*E[i][j][k];

            if (E[i][j][k] > mx)
            mx = E[i][j][k];
            }
        }
  }
  printf("Before sqrt: %f \n", l2norm);
  //l2norm /= (float) ((d)*(d)*(d));
  l2norm = sqrt(l2norm);
  printf("radius: %d, iteration %d \n", r, iters);
  printf("max: %20.12e, l2norm: %0.10f \n", mx, l2norm);
}

int main(int argc, char *argv[])
{
  int i, j, k;
  int inside, inside_x, inside_xx;
  double ***temp, ***Unew, ***Uold;
  double ***tensor_x, ***tensor_y, ***tensor_z;

  clock_t begin = clock();

  int r = 100;
  int d = r*2+1;

  int center_x = r+1;
  int center_y = r+1;
  int center_z = r+1;

  /*double kxx = 0.04;
  double kyy = 0.01;
  double kzz = 0.02;*/

  //for r = 29
  double kxx = 0.000099;
  double kyy = 0.000082;
  double kzz = 0.00009;

  srand((unsigned)time(NULL));
  //int rand_x = rand()%1+1;

  double size_x = d;
  double size_y = d;
  double size_z = d;

  double dx = 1.0/size_x;
  double dy = 1.0/size_y;
  double dz = 1.0/size_z;

  double stress_x = kxx/(2*dx*dx);
  double stress_y = kyy/(2*dy*dy);
  double stress_z = kzz/(2*dz*dz);

  Unew = dallocate_3d(d+2, d+2, d+2);
  Uold = dallocate_3d(d+2, d+2, d+2);
  tensor_x = dallocate_3d(d+2, d+2, d+2);
  tensor_y = dallocate_3d(d+2, d+2, d+2);
  tensor_z = dallocate_3d(d+2, d+2, d+2);
  dinit_3d(Unew, d+2, d+2, d+2);
  dinit_3d(Uold, d+2, d+2, d+2);
  dinit_3d(tensor_x, d+2, d+2, d+2);
  dinit_3d(tensor_y, d+2, d+2, d+2);
  dinit_3d(tensor_z, d+2, d+2, d+2);
  double count = 1.0;

  //fill interior points, 1 if inside, else 0
  int testing = 1;
  for(i = 1; i <= d; i++)
  {
    for(j = 1; j <= d; j++)
    {
      for(k = 1; k <= d; k++)
      {
        inside = (((i-center_x)*(i-center_x)) + ((j-center_y)*(j-center_y)) + ((k-center_z)*(k-center_z)));

        if(inside <= r*r)
        {
          //if(count == 3)
            //count = 1.0;

          Uold[i][j][k] = 50.0;
          /*if(testing % 2 == 1)
            Uold[i][j][k] = 2;
          else
            Uold[i][j][k] = 3;

          testing++;*/
          //count++;
        }
        else
        {
          Uold[i][j][k] = 0.0;
        }

        count++;
      }
    }
  }

  for(i = 1; i <= d; i++)
  {
    for(j = 1; j <= d; j++)
    {
      for(k = 1; k <= d; k++)
      {
        inside = (((i-center_x)*(i-center_x)) + ((j-center_y)*(j-center_y)) + ((k-center_z)*(k-center_z)));

        //if(inside <= r*r)
        //{
          if((Uold[i-1][j][k] != 0 && Uold[i+1][j][k] != 0) && (Uold[i][j-1][k] != 0 && Uold[i][j+1][k] != 0) && (Uold[i][j][k-1] != 0 && Uold[i][j][k+1] != 0))
          {
            tensor_x[i][j][k] = 159.0;
            tensor_y[i][j][k] = 160;
            tensor_z[i][j][k] = 160;
          }
          else
          {
            tensor_x[i][j][k] = 0.0;
            tensor_y[i][j][k] = 0.0;
            tensor_z[i][j][k] = 0.0;
          }
        /*}
        else
        {
          tensor_x[i][j][k] = 0.0;
          tensor_y[i][j][k] = 0.0;
          tensor_z[i][j][k] = 0.0;
        }*/
      }
    }
  }

  int max_time = 50;
  int time_iter = 0;

  while(time_iter < max_time)
  {
    //interior points from 1 to d+1
    for(i = 1; i <= d; i++)
    {
      for(j = 1; j <= d; j++)
      {
        for(k = 1; k <= d; k++)
        {
          /*Unew[i][j][k] = stress_x*(Uold[i+1][j][k] - Uold[i-1][j][k]) +
            stress_y*(Uold[i][j+1][k] - Uold[i][j-1][k]) +
            stress_z*(Uold[i][j][k+1] - Uold[i][j][k-1]);*/

          Unew[i][j][k] = ((tensor_x[i][j][k]/2*dx*dx)*(Uold[i+1][j][k] - Uold[i-1][j][k])) +
            ((tensor_y[i][j][k]/2*dy*dy)*(Uold[i][j+1][k] - Uold[i][j-1][k])) +
            ((tensor_z[i][j][k]/2*dz*dz)*(Uold[i][j][k+1] - Uold[i][j][k-1]));
          //if(Unew[i][j][k] >= 1.0)
            //printf("%f \n", Unew[i][j][k]);

            //printf(" %f \n", tensor_x[i][j][k]);
        }
      }
    }

    temp = Uold;
    Uold = Unew;
    Unew = temp;

    time_iter++;
  }

  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  printf("Elapsed time: %f \n", time_spent);

  /*int xx = 1;
  for(i = 1; i <= d; i++)
  {
    for(j = 1; j <= d; j++)
    {
      for(k = 1; k <= d; k++)
      {
        //printf("");
        printf("%0.1f \t", Uold[i][j][k]);
      }
      printf("\n");
    }
    printf("\n");
  }
  printf("\n");*/

  //srand((unsigned) time(&rand_time));

  //printf("%f \n", rand_x);


  calculatel2Norm(Uold, d, r, max_time);

  return 0;
}
