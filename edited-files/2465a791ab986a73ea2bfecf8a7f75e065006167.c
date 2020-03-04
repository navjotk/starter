#define _POSIX_C_SOURCE 200809L
#include "stdlib.h"
#include "math.h"
#include "sys/time.h"
#include "omp.h"

struct dataobj
{
  void *restrict data;
  int * size;
  int * npsize;
  int * dsize;
  int * hsize;
  int * hofs;
  int * oofs;
} ;

struct profiler
{
  double section0;
} ;


int padfunc(struct dataobj *restrict phi_vec, const int x_M, const int y_M, const int z_M, const int abc_x_l_ltkn, const int abc_x_r_rtkn, const int abc_y_l_ltkn, const int abc_y_r_rtkn, const int abc_z_l_ltkn, const int abc_z_r_rtkn, struct profiler * timers, const int x_m, const int y_m, const int z_m)
{
  float (*restrict phi)[phi_vec->size[1]][phi_vec->size[2]] __attribute__ ((aligned (64))) = (float (*)[phi_vec->size[1]][phi_vec->size[2]]) phi_vec->data;
  #pragma omp target enter data map(to: phi[0:phi_vec->size[0]][0:phi_vec->size[1]][0:phi_vec->size[2]])
  struct timeval start_section0, end_section0;
  gettimeofday(&start_section0, NULL);
  /* Begin section0 */
  for (int abc_x_l = x_m; abc_x_l <= abc_x_l_ltkn + x_m - 1; abc_x_l += 1)
  {
    #pragma omp target teams distribute parallel for collapse(2)
    for (int y = y_m; y <= y_M; y += 1)
    {
      for (int z = z_m; z <= z_M; z += 1)
      {
        phi[abc_x_l + 2][y + 2][z + 2] = phi[12][y + 2][z + 2];
      }
    }
  }
  for (int abc_x_r = -abc_x_r_rtkn + x_M + 1; abc_x_r <= x_M; abc_x_r += 1)
  {
    #pragma omp target teams distribute parallel for collapse(2)
    for (int y = y_m; y <= y_M; y += 1)
    {
      for (int z = z_m; z <= z_M; z += 1)
      {
        phi[abc_x_r + 2][y + 2][z + 2] = phi[x_M - 8][y + 2][z + 2];
      }
    }
  }
  #pragma omp target teams distribute parallel for collapse(1)
  for (int x = x_m; x <= x_M; x += 1)
  {
    for (int abc_y_l = y_m; abc_y_l <= abc_y_l_ltkn + y_m - 1; abc_y_l += 1)
    {
      for (int z = z_m; z <= z_M; z += 1)
      {
        phi[x + 2][abc_y_l + 2][z + 2] = phi[x + 2][12][z + 2];
      }
    }
    for (int abc_y_r = -abc_y_r_rtkn + y_M + 1; abc_y_r <= y_M; abc_y_r += 1)
    {
      for (int z = z_m; z <= z_M; z += 1)
      {
        phi[x + 2][abc_y_r + 2][z + 2] = phi[x + 2][y_M - 8][z + 2];
      }
    }
    for (int y = y_m; y <= y_M; y += 1)
    {
      for (int abc_z_l = z_m; abc_z_l <= abc_z_l_ltkn + z_m - 1; abc_z_l += 1)
      {
        phi[x + 2][y + 2][abc_z_l + 2] = phi[x + 2][y + 2][12];
      }
      for (int abc_z_r = -abc_z_r_rtkn + z_M + 1; abc_z_r <= z_M; abc_z_r += 1)
      {
        phi[x + 2][y + 2][abc_z_r + 2] = phi[x + 2][y + 2][z_M - 8];
      }
    }
  }
  /* End section0 */
  gettimeofday(&end_section0, NULL);
  timers->section0 += (double)(end_section0.tv_sec-start_section0.tv_sec)+(double)(end_section0.tv_usec-start_section0.tv_usec)/1000000;
  #pragma omp target update from(phi[0:phi_vec->size[0]][0:phi_vec->size[1]][0:phi_vec->size[2]])
  #pragma omp target exit data map(release: phi[0:phi_vec->size[0]][0:phi_vec->size[1]][0:phi_vec->size[2]])
  return 0;
}
