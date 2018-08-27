#include <stdio.h>

main(){
    
  // 利用 pipe 呼叫 gnuplot 繪圖
  FILE *pipe;
  pipe = popen("gnuplot -persist","w");
  //fprintf(pipe,"set term png enhanced font \"v:/fireflysung.ttf\" 18 \n");
  fprintf(pipe,"set term png enhanced font \"y:/wqy-microhei.ttc\" 12 \n");
  //fprintf(pipe,"set yrange [68:70]\n");
  fprintf(pipe,"set output \"ex1_3d.png\"\n");
  //fprintf(pipe, "plot \"osc.dat\" title \"Runge-Kutta 解微分方程式\" with lines\n");
  //fprintf(pipe, "plot [0:3.14159*2] sin(x) with linespoints\n");
  fprintf(pipe,"set hidden3d\n");
  fprintf(pipe,"set isosamples 30\n");
  fprintf(pipe, "splot [-2:2][-2:2] exp(-(x**2 + y**2))*cos(x/4)\n");
  fprintf(pipe,"quit\n");
  pclose(pipe);
}