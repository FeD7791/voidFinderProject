#include <stdio.h>
#include <stdlib.h>
//#include "user.h"

struct tracers {
  float Pos[3];
  float Vel[3];
  float Cen[3];
  float Delta;
  float Volume;
  float size;
}; 


int main(){
    return 0;
}

void c_binary_writter(
    double *arr_x,
    double *arr_y,
    double *arr_z,
    double *arr_vx,
    double *arr_vy,
    double *arr_vz,
    double *arr_m , int size
){
    
    FILE *fout = fopen("tracers_zobov.raw","wb");
    fwrite(&size,sizeof(int),1,fout);
    fwrite(arr_x,sizeof(double),size,fout);
    fwrite(arr_y,sizeof(double),size,fout);
    fwrite(arr_z,sizeof(double),size,fout);
    fclose(fout);
    FILE *fout2 = fopen("tracers_zobov.txt","w");
    fprintf(fout2, "%d\n", size);
    int i;
    for (i = 0; i < size; i++) {
    fprintf(fout2, "%.10f %.10f %.10f\n", arr_x[i], arr_y[i], arr_z[i]);
  }
    fclose(fout2);

}
