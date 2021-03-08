#include<stdio.h>
#include<math.h>
#include<memory.h>
#include"Forward.h"
#include"Inversion.h"

//double Forward(double x_0, double y_0, double z_0, double x_1, double y_1, double z_1, double x, double y, double M_1,double times)
//					第一个角点三个坐标				第二个角点三个坐标				所求测点的坐标		磁化强度     第几个剖分单元

#define dx 100	//x间距
#define dy 100	//y间距
#define M1 1	//A/m

#pragma warning(disable:4996)

//M个模型，N个测点

int num = 6;
//第几次读取角点坐标，最后一行角点坐标为计算背景场B_T的数据

int main()
{
	FILE* fp1, * fp2, * fp3, * fp4, * fp5;
	double x_0[2], y_0[2], z_0[2];
	double R_T[441] = { 0 };//d
	//R_T模型磁异常

	int line = 21, point = 21;	//测线和测点数
	int times = 0, part = 0;	
	//times 表示体元对第times个测点的异常影响，
	int i, j, x, y;				

	//double G[21 * 21][5 * 5 * 2 * 6];
	double* G = (double*)malloc(sizeof(double) * 441*4410);
	//核矩阵按照先正演第一层的x方向的体元
	//正演出的核矩阵;利用核矩阵和背景场的大小计算出背景异常;
	
	//背景磁异常
	double* B_T = (double*)malloc(sizeof(double) * 441);	
	//将动态数组G和B_T初始化为0;

	//背景+正演磁异常
	double* ALL_T = (double*)malloc(sizeof(double) * 441);
	//反演出的m
	double* m = (double*)malloc(sizeof(double) * 4410);

	memset(G, 0, sizeof(double) * 441 * 300);
	memset(B_T, 0, sizeof(double) * 441);
	memset(ALL_T, 0, sizeof(double) * 441);

	fp1 = fopen("C:\\Users\\wc\\Desktop\\台阶正演参数.txt", "rb+");
	fp2 = fopen("C:\\Users\\wc\\Desktop\\台阶正演.txt", "wb+");
	fp3 = fopen("C:\\Users\\wc\\Desktop\\背景磁异常.txt", "wb+");
	fp4 = fopen("C:\\Users\\wc\\Desktop\\背景+模型磁异常.txt", "wb+");
	fp5 = fopen("C:\\Users\\wc\\Desktop\\反演.txt", "wb+");
	fprintf(fp2, "x\ty\tT\n");	
	fprintf(fp3, "x\ty\tT\n");	
	fprintf(fp4, "x\ty\tT\n");	
	do {
		int ret = fscanf_s(fp1, "%lf%lf%lf%lf%lf%lf", &x_0[0], &y_0[0], &z_0[0], &x_0[1], &y_0[1], &z_0[1]);		
		if (ret > 0){	
			times = 0;						
			for (i = 0; i < 21; i++){
				x = dx * i;
				for (j = 0; j < 21; j++){
					y = dy * j;																	//440
					if (part != num){
						R_T[times] += Forward(x_0[0], y_0[0], z_0[0], x_0[1], y_0[1], z_0[1], x, y, M1, times, G, part);
						times++;
					}
					else {
						B_T[times] += Forward(x_0[0], y_0[0], z_0[0], x_0[1], y_0[1], z_0[1], x, y, M1, times, G, part);
						fprintf(fp2, "%d\t%d\t%lf\n", x, y, R_T[times]);
						fprintf(fp3, "%d\t%d\t%lf\n", x, y, B_T[times]);
						ALL_T[times] = R_T[times] + B_T[times];
						fprintf(fp4, "%d\t%d\t%lf\n", x, y, ALL_T[times++]);
					}
				}
			}			
			part++;
		}
		if (feof(fp1)) break;
	}while(1);
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	fclose(fp4);
	m = Inversion(G, ALL_T, 2, 1);
	fprintf(fp5, "T\n");
	for (i = 0; i < 4410; i++)
	{
		fprintf(fp5, "%lf\n", *(m + i));
	}

	fclose(fp5);
}