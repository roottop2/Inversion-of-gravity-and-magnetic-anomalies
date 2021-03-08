#pragma once
#include<stdio.h>
#include<math.h>
//基础变量
#define PI 3.1415926
#define miu0 4*PI*1e-7
#define I 90*PI/180	//磁倾角
#define A 0			//磁偏角
#define M cos(I)*sin(A)
#define N sin(I)
#define L cos(I)*cos(A)
#define k1 2*M*N
#define k2 2*N*L
#define k3 2*M*L
#define k4 L*L
#define k5 M*M
#define k6 -N*N

/*0<x<2000,0<y<2000,0<z<1000*/
/*每层分为21*21个体元*/
/*则dx=100,dy=100*/
/*输入对角点的坐标(x_0,y_0,z_0),(x_1,y_1,z_1)*/
double Forward(double x_0, double y_0, double z_0, double x_1, double y_1, double z_1, double x, double y, double M_1, int times,double *G,int part)
{
	double x1[2], y1[2], z1[2], R;
	extern num;

	int divide_xy,divide_z;
	if (part != num)
		divide_xy = 5,divide_z = 2;
	else
		divide_xy = 21, divide_z = 10;

	double dx = (x_1 - x_0) / divide_xy, dy = (y_1 - y_0) / divide_xy, dz = (z_1 - z_0) / divide_z;
	double T = 0, mid = 0;
	int i, j, k, l, m, n, point = 0;	
	for (i = 0; i < divide_xy; i++){
		x1[0] = i * dx + x_0 - x, x1[1] = (i + 1) * dx + x_0 - x;
		for (j = 0; j < divide_xy; j++){
			y1[0] = j * dy + y_0 - y, y1[1] = (j + 1) * dy + y_0 - y;
			for (k = 0; k < divide_z; k++)	{
				z1[0] = k * dz + z_0, z1[1] = (k + 1) * dz + z_0;

				//一个剖分单元对测点的影响值
				for (l = 2; l > 0; l--)
					for (m = 2; m > 0; m--)
						for (n = 2; n > 0; n--){
								R = sqrt(x1[l - 1] * x1[l - 1] + y1[m - 1] * y1[m - 1] + z1[n - 1] * z1[n - 1]);
								if (l * m * n == 8 || l * m * n == 2){
;									T += miu0 / (4 * PI) * M_1 * (k1 * log(R + x1[l - 1]) + k2 * log(R + y1[m - 1]) + k3 * log(R + z1[n - 1]) + k4 * atan(x1[l - 1] * y1[m - 1] / (pow(x1[l - 1], 2) + R * z1[n - 1] + pow(z1[n - 1], 2))) + k5 * atan(x1[l - 1] * y1[m - 1] / (pow(y1[m - 1], 2) + R * z1[n - 1] + pow(z1[n - 1], 2))) + k6 * atan(x1[l - 1] * y1[m - 1] / (R * z1[n - 1]))) * 1e9;
									mid += miu0 / (4 * PI) * M_1 * (k1 * log(R + x1[l - 1]) + k2 * log(R + y1[m - 1]) + k3 * log(R + z1[n - 1]) + k4 * atan(x1[l - 1] * y1[m - 1] / (pow(x1[l - 1], 2) + R * z1[n - 1] + pow(z1[n - 1], 2))) + k5 * atan(x1[l - 1] * y1[m - 1] / (pow(y1[m - 1], 2) + R * z1[n - 1] + pow(z1[n - 1], 2))) + k6 * atan(x1[l - 1] * y1[m - 1] / (R * z1[n - 1]))) * 1e9;
								}
								else {
									T -= miu0 / (4 * PI) * M_1 * (k1 * log(R + x1[l - 1]) + k2 * log(R + y1[m - 1]) + k3 * log(R + z1[n - 1]) + k4 * atan(x1[l - 1] * y1[m - 1] / (pow(x1[l - 1], 2) + R * z1[n - 1] + pow(z1[n - 1], 2))) + k5 * atan(x1[l - 1] * y1[m - 1] / (pow(y1[m - 1], 2) + R * z1[n - 1] + pow(z1[n - 1], 2))) + k6 * atan(x1[l - 1] * y1[m - 1] / (R * z1[n - 1]))) * 1e9;
									mid -= miu0 / (4 * PI) * M_1 * (k1 * log(R + x1[l - 1]) + k2 * log(R + y1[m - 1]) + k3 * log(R + z1[n - 1]) + k4 * atan(x1[l - 1] * y1[m - 1] / (pow(x1[l - 1], 2) + R * z1[n - 1] + pow(z1[n - 1], 2))) + k5 * atan(x1[l - 1] * y1[m - 1] / (pow(y1[m - 1], 2) + R * z1[n - 1] + pow(z1[n - 1], 2))) + k6 * atan(x1[l - 1] * y1[m - 1] / (R * z1[n - 1]))) * 1e9;
								}
						}
				if (part == num){
					*(G + times * 4410 + k + j * divide_z + i * divide_xy * divide_z) = mid;
					//double* G = (double*)malloc(sizeof(double) * 441*4410);
					mid = 0;
				}
			}			
		}
	}
	return T;
}