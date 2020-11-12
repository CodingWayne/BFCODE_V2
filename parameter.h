#ifndef test_H
#define test_H

#include <stdio.h>
#include <stdlib.h>
#include "uclib.h"

struct BFKIN
{
	int np;
	int itestar;
	int ci;
	int wake_offon;
	int wake_ci;
	int wake_onite;

	int self_offon;
	int firstJ_offon;
	double firstJ;
	int bodyforce_offon;
	int P_3;

	int Unsteady_offon;
	int unsteady_c;
	int UN_ci;
	double theta_one_blade;
	int tnb;
	int disk_offen;
	int UNB;

	double R;
	double x;
	double xp;
	int plane_r_numb;
	int plane_x_numb;
	int plane_c_numb;

	double Vs;
	double rho;
	double wr;
	double ar;
	double sfc;

	int bow_dir;
	double p_xo;
	double p_yo;
	double p_zo;
	int st;

	double* bfkin_hub;

};
double reitestar, pitch;

int ite;
char str[50];
int i, j;
/*-------------------------網格變數---------------------------
 *meshload為判斷變數，抓完網格後=1
 *minCellid為最小網格編號
 *threadboy為判斷變數，等於1時表此時用單一process開啟對應之mesh檔，此時建立newsize矩陣及threadnumber矩陣，分別開啟各thread%d_%d.txt，輸入pid到threadnumber,size到newsize裡
 *allsize=newsize矩陣各元素加總，但是不等於圓盤網格數(會有重複項)
 *kk為占=暫存變數
 *numb,non僅為記錄用
 */
int  meshload,minCellid,np, threadboy, * newsize, * threadnumber, all_size ;
double** mesh;
int numb, kk;
double non;
//-------------
void Load_bfkin(struct BFKIN* bfkin);
void Output_bfkin(struct BFKIN* bfkin, double* reitestar, double* pitch);


#endif