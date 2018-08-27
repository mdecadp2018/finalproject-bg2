#include <stdio.h>

#define problem_order 2

double x_prime(double x1[], double tt, int i);
double integrate(double x[], double tt, int i);
double t, step;

// 資料檔案存檔用
FILE *output;
// 利用 pipe 呼叫 gnuplot 繪圖
FILE *pipe;

int main(){
  double initial_time, final_time;
  int num_of_data, j;
  double x[problem_order];
  output=fopen("osc.dat", "w");
  // initial setup
  initial_time = 0;
  final_time = 30;
  num_of_data = 500;
  step = (final_time - initial_time)/num_of_data;
  t = initial_time;
  // set the initial conditions
  x[0] = 1;
  x[1] = 0;
  // x[2] = 0.1;
  for(t = initial_time; t < final_time; t+=step){
    printf("time = %-5.5f, displacement = %-5.5f, velocity = %-5.5f\n", t, x[0], x[1]);
    for(j =0; j < problem_order; j++){
      x[j] = integrate(x, t, j);
    }
    fprintf(output, "%f\t%f\n", t, x[0]);
  }
  fclose(output);

  pipe = popen("gnuplot -persist","w");
  //fprintf(pipe,"set term png enhanced font \"v:/fireflysung.ttf\" 18 \n");
  fprintf(pipe,"set term png enhanced font \"y:/wqy-microhei.ttc\" 18 \n");
  //fprintf(pipe,"set yrange [68:70]\n");
  fprintf(pipe,"set output \"c_rk_4.png\"\n");
  fprintf(pipe, "plot \"osc.dat\" title \"C 以 Runge-Kutta 解微分方程式\" with lines\n");
  fprintf(pipe,"quit\n");
  pclose(pipe);
return 0;
}

// x is an array which represents the big matrix of the state variables
double x_prime(double x[], double tt, int i){
    double k, b, m, xp[problem_order];
    k=1;
    b=0;
    m=1;
    xp[0]=x[1];
    xp[1]=-x[0]-0.2*x[1];
    // xp[2]=x[0];
    return xp[i];
}


// Integrate the differential equation using the Runge-Kutta Method

// x is an array of the matrix which represents the state variables

double integrate(double x[], double tt, int i){

// the purpose of this function is to calculate the value of the $x for each step
// upon time series, calculate corresponding $x array
// there we introduct the array k1,k2,k3,k4
//k1[0],k2[0],k3[0],k4[0] is corresponding to x[0]
//k1[1],k2[1],k3[1],k4[1] is corresponding to x[1]
// initialize $order
// order=0;
// step=0.1;
double k1[problem_order], k2[problem_order], k3[problem_order], k4[problem_order];
double x1[problem_order], x2[problem_order], x3[problem_order], x4[problem_order];
double t2, t3, t4;
double h = step/2;
int k, j;
for(j = 0; j < problem_order; j++){
  k1[j] = x_prime(x, tt, j);
}

for(k = 0; k < problem_order; k++){
  x2[k] = x[k] + (h/2)*k1[k];
}

t2 = tt + (h/2);

for(j = 0; j < problem_order; j++){
  k2[j] = x_prime(x2, t2, j);
}

for(k = 0; k < problem_order; k++){
  x3[k] = x[k] + (h/2)*k2[k];
}

t3 = tt + (h/2);

for(j = 0; j < problem_order; j++){
  k3[j] = x_prime(x3, t3, j);
}
for(k = 0; k < problem_order; k++){
  x4[k] = x[k] + (h)*k3[k];
}
t4 = tt + (h);
for(j = 0; j < problem_order; j++){
  k4[j] = x_prime(x4, t4, j);
}

for(k = 0;k < problem_order; k++){
  x[k] = x[k] + (h/6)*(k1[k] + (2.0*k2[k]) + (2.0*k3[k]) + k4[k]);
}
return x[i];

}
