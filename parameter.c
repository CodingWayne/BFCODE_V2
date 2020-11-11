#include <stdio.h>
#include <stdlib.h>
#include "uclib.h"
#include "parameter.h"

void Load_bfkin(struct BFKIN* bfkin)
{
	int i;
	FILE* fp;
	fp = fopen("bf.kin", "r+");
	fscanf(fp, "%d%d%d%d%d%d", &bfkin->np, &bfkin->itestar, &bfkin->ci, &bfkin->wake_offon, &bfkin->wake_ci, &bfkin->wake_onite);  //從bf->kin讀入; np:所使用的核心數; itestar:開始計算ci前的先行計算步數; ci:每幾步迭代一次; wake_offon:是否開啟wake輸出; wake_ci:每幾步輸出wake一次; wake_onite:開始輸出wake的步數; 
	fscanf(fp, "%d%d%lf%d%d", &bfkin->self_offon, &bfkin->firstJ_offon, &bfkin->firstJ, &bfkin->bodyforce_offon, &bfkin->P_3);		//self_offon:是否開啟self;	  firstJ_offon:是否開啟自定義第一次J值;  firstJ:自定義的第一次J值;   bodyforce_offon:是否開啟bodyforce;   
	fscanf(fp, "%d%d%d%lf%d%d%d", &bfkin->Unsteady_offon, &bfkin->unsteady_c, &bfkin->UN_ci, &bfkin->theta_one_blade, &bfkin->tnb, &bfkin->disk_offen, &bfkin->UNB);	     // Blade_offon:是否開啟葉片模式 ;  Unsteady_offon:是否開啟Unsteady bodyforce ; unsteady_c: ubf開始進行步數 ; theta_one_blade: 投影螺槳面積  ;  tnb: 每一次迭代螺槳葉面轉動圈數
	fscanf(fp, "%lf%lf%lf%d%d%d", &bfkin->R, &bfkin->x, &bfkin->xp, &bfkin->plane_r_numb, &bfkin->plane_x_numb, &bfkin->plane_c_numb); //R:圓盤半徑; x:圓盤厚度; xp:在x軸上的位置; plane_r_numb:圓盤r方向之網格分割數; plane_x_numb:圓盤x方向之網格分割數; plane_c_numb:圓盤周方向之網格分割數
	fscanf(fp, "%lf%lf%lf%lf%lf", &bfkin->Vs, &bfkin->rho, &bfkin->wr, &bfkin->ar, &bfkin->sfc);					//Vs:船速; rho:密度; WR:興波阻力; AR:波浪附加阻力; SFC:摩擦阻力修正量;
	fscanf(fp, "%d%lf%lf%lf%d", &bfkin->bow_dir, &bfkin->p_xo, &bfkin->p_yo, &bfkin->p_zo, &bfkin->st);	//圓盤中心座標 ; 船艏方向
	bfkin->bfkin_hub = (double*)malloc(bfkin->plane_x_numb * sizeof(double));
	for (i = 0; i != bfkin->plane_x_numb; ++i)fscanf(fp, "%lf", &bfkin->bfkin_hub[i]);  //x方向劃分之hub位置
	bfkin->st = bfkin->st * 1000;
	fclose(fp);
}

void Output_bfkin(struct BFKIN* bfkin, double* reitestar, double* pitch)
{
	int i;
	FILE* fp3;
	fp3 = fopen("Check_bfkin.dat", "w");
	fprintf(fp3, "np=%d\n", bfkin->np);
	fprintf(fp3, "itestart=%d\n", bfkin->itestar);
	fprintf(fp3, "ci=%d\n", bfkin->ci);
	fprintf(fp3, "wake_offon=%d\n", bfkin->wake_offon);
	fprintf(fp3, "wake_ci=%d\n", bfkin->wake_ci);
	fprintf(fp3, "wake_onite=%d\n", bfkin->wake_onite);
	fprintf(fp3, "firstJ_offon=%d\n", bfkin->firstJ_offon);
	fprintf(fp3, "firstJ=%lf\n", bfkin->firstJ);
	fprintf(fp3, "R=%lf\n", bfkin->R);
	fprintf(fp3, "x=%lf\n", bfkin->x);
	fprintf(fp3, "xp=%lf\n", bfkin->xp);
	fprintf(fp3, "P_3=%d\n", bfkin->P_3);
	fprintf(fp3, "Unsteady_offon=%d\n", bfkin->Unsteady_offon);
	fprintf(fp3, "disk_offen=%d\n", bfkin->disk_offen);
	fprintf(fp3, "plane_r_numb=%d\n", bfkin->plane_r_numb);
	fprintf(fp3, "plane_x_numb=%d\n", bfkin->plane_x_numb);
	fprintf(fp3, "plane_c_numb=%d\n", bfkin->plane_c_numb);
	fprintf(fp3, "Vs=%lf\n", bfkin->Vs);
	fprintf(fp3, "rho=%lf\n", bfkin->rho);
	fprintf(fp3, "Wave Resistance=%lf\n", bfkin->wr);
	fprintf(fp3, "Added Resistance=%lf\n", bfkin->ar);
	fprintf(fp3, "SFC=%lf\n", bfkin->sfc);
	fprintf(fp3, "reitestar=%d\n", *reitestar);
	fprintf(fp3, "bfkin_hub={");
	for (i = 0; i != bfkin->plane_x_numb; ++i) {
		fprintf(fp3, " %lf", bfkin->bfkin_hub[i]);
	}
	fprintf(fp3, "}\n");
	fprintf(fp3, "Pitch=%f\n", *pitch);
	fprintf(fp3, "bow_dir=%d\n", bfkin->bow_dir);
	fprintf(fp3, "P_X = %lf\tP_Y = %lf\tP_Z = %lf\n", bfkin->p_xo, bfkin->p_yo, bfkin->p_zo);
	fclose(fp3);
}
