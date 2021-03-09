#pragma once
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
//剖分网格个数 20*20*10
//长度DX=DY=DZ=100;
//
double Matrix_multiplication(double* matrix_1, double* matrix_2, double* output_matrix, int line1, int column1，int line2, int column2)
{
	int i, j, k;
	double sum = 0;
        for (i = 0; i < line; i++)
		for (j = 0; j < column2; j++)
			for (k = 0; k < column1; k++)
			{
                                sum += *(matrix_1 + i * column1 + k) * *(matrix_2 + k * column2 + j);
	                 	*(output_matrix + i * column2 + j) = sum;
		                sum = 0;
		        }
	return 0;
}
//A正演矩阵
//反演初值m 0     
//迭代次数k<N
//正则化参数0.6<miu<0.9
//聚焦因子e = 10e-10
//m = Wm逆*We逆*mw
//Aw = Wd*A*Wm逆*We逆
//dw = Wd*d
//mw = We*Wm*m
//Q = We逆的平方*Aw的转置*Aw+miu0*e*e*I
//q = We逆的平方*Aw的转置*dw 
//φ在mw0处的导数为f0 = Q*mw-q
//沿mw0的初始搜索方向选择目标函数φ的最速下降方向d0 = -f0
//对应搜索步长t = (d0转置*f0)/(d0转置*Q*d0)
//背景模型mb = a  , db = Amb  ,d = dobs+db

//正演磁化强度mc = 1A/m  背景磁化强度1A/m   约束M_min= 0 M_max = 2;

//最大迭代次数
#define  Maximum_number_of_iterations 100

//聚焦因子
#define e 10e-10
//参数 G 核矩阵  ALL_T总磁异常  约束条件M_min= 0 M_max = 2
double *Inversion(double* G, double* ALL_T, double M_max, double M_min)
{
	//迭代次数
	int Iterations = 0;
	//循环变量
	int i, j;
	//正则化参数
	double miu = 0.6;//(0.6~0.9)

	double sum = 0;

	//数据加权矩阵Wd;聚焦权矩阵We;模型加权矩阵Wd;背景模型mb(磁化强度，此处为1A/m),d为观测数据;
	//Gm = d
	//核矩阵G 441行 4410列
	double* Wd, * We, * Wm, * Aw, t0, beta;
	Aw = (double*)malloc(sizeof(double) * 441*4410);
	Wd = (double*)malloc(sizeof(double) * 441);	
	We = (double*)malloc(sizeof(double) * 4410);
	Wm = (double*)malloc(sizeof(double) * 4410);
	//double Wd[441], We[4410], Wm[4410],mb,d;
	//We = diag(1/sqrt(m*m+e*e))
	//Wd = diag(sqrt(A*A转置))
	//Wm = diag(sqrt(A转置*A))

	double* G_T = (double*)malloc(sizeof(double) * 4410 * 441);

	double* Aw_T = (double*)malloc(sizeof(double) * 4410 * 441);

	double* Q = (double*)malloc(sizeof(double) * 4410 * 4410);

	double* q = (double*)malloc(sizeof(double) * 4410);

	double* m = (double*)malloc(sizeof(double) * 4410);

	double* dw = (double*)malloc(sizeof(double) * 441);

	double* mid = (double*)malloc(sizeof(double) * 4410);

	double* f0 = (double*)malloc(sizeof(double) * 4410);

	double* f1 = (double*)malloc(sizeof(double) * 4410);

	double* d0 = (double*)malloc(sizeof(double) * 4410);

	memset(m,1,sizeof(double)*4410);
	memset(dw,0,sizeof(double)*441);
	memset(mid,0,sizeof(double)*4410);
	memset(d0,0,sizeof(double)*4410);
	memset(Q,0,sizeof(double)*4410);
	memset(f1,0,sizeof(double)*4410);

	for (i = 0; i < 4410; i++)
		for (j = 0; j < 441; j++)
			*(G_T + 441 * i + j) = *(G + 4410 * j + i);

	//We Wm Wd计算
	
	for (i = 0; i < 441; i++)
	{
		for (j = 0; j < 4410; j++)
			sum += *(G + i * 4410 + j) * *(G_T + j * 441 + i);
		*(Wd + i) = sqrt(sum);
		sum = 0;
	}
	for (i = 0; i < 4410; i++)
	{
		for (j = 0; j < 441; j++)
			sum += *(G_T + i * 441 + j) * *(G + j * 4410 + i);
		*(Wm + i) = sqrt(sum);
		sum = 0;
	}

	//m0w
	
	
	//Aw与Aw_T
	
	
	//计算Q矩阵
	
	//计算dw  dw = Wd * d;   441 441   441 1  =441 1
	for (i = 0; i < 441; i++)
		*(dw + i) = *(ALL_T + i) * *(Wd + i);
	//计算q
	

	//计算f0 = Q*mw - q


	//计算t0
	do {
		//更新We 与 mwk;
		for (i = 0; i < 4410; i++){
			We[i] = 1.0 / sqrt(m[i] * m[i] + e * e);
			m[i] = We[i] * Wm[i] * m[i];
		}
		//更新Aw
		for (i = 0; i < 441; i++)
			for (j = 0; j < 4410; j++)
				*(Aw + i * 4410 + j) = *(G + i * 4410 + j) * *(Wd + i) * 1.0 / *(Wm + j) * 1.0 / *(We + j);
		for (i = 0; i < 4410; i++)
			for (j = 0; j < 441; j++)
				*(Aw_T + i * 441 + j) = *(Aw + j * 4410 + i);
		//更新Q
		Matrix_multiplication(Aw_T, Aw, Q, 4410, 441);
		for (i = 0; i < 4410; i++)
			for (j = 0; j < 4410; j++) {
				*(Q + i * 4410 + j) *= 1.0 / *(We + i) * 1.0 / *(We + i);
				while (i == j)
					*(Q + i * 4410 + j) += miu * e * e;
			}
		//更新q
		for (i = 0; i < 4410; i++) {
			for (j = 0; j < 441; j++)
				*(mid + i) += *(Aw_T + 441 * i + j) * *(dw + j);
			*(q + i) = 1.0 / *(We + i) * 1.0 / *(We + i) * *(mid + i);
		}
		//fk
		for (i = 0; i < 4410; i++)
			for (j = 0; j < 4410; j++)
				*(f1 + i) += *(Q + i * 4410 + j) * *(m + j) - *(q + i);
		double molecule = 0, denominator = 0;
		for (i = 0; i < 4410; i++) {
			if(Iterations == 0)
				*(d0 + i) = -*(f1 + i);
			else
				for (j = 0; j < 4410; j++)
				{
					molecule += *(f1 + j) * *(f1 + j);
					denominator += *(f0 + j) * *(f0 + j);
					beta = molecule / denominator;
					*(d0 + i) = -*(f1 + i) + beta * *(d0 + i);
				}
			molecule = 0, denominator = 0;
			molecule += *(d0 + i) * *(f0 + i);
			for (j = 0; j < 4410; j++)
				*(mid + i) += *(d0 + j) * *(Q + i * 4410 + j);
			denominator += *(mid + i) * *(d0 + i);
		}
		t0 = molecule / denominator;

		for (i = 0; i < 4410; i++)
			*(f0 + i) = *(f1 + i);

		Iterations++;
		for (i = 0; i < 4410; i++){
			*(m + i) = *(m + i) + t0 * *(d0 + i);
			*(m + i) = 1.0 / (*(Wm + i) * *(We + i)) * *(m + i);
			if (*(m + i) > M_max)
				*(m + i) = M_max;
			else if (*(m + i) < M_min)
				*(m + i) = M_min;
		}

	} while (Iterations < Maximum_number_of_iterations);
	for (i = 0; i < 4410; i++)
		*(m + i) -= 1;
	return m;
}


//转置与非转置矩阵相乘
